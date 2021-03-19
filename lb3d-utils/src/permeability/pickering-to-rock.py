#!/usr/bin/env python
# -*- coding: utf-8 -*-

# These 2d and 3d functions should work for generic multidimensional arrays.
def map_2d(fun, d2):
    return map(lambda d1: map(fun, d1), d2)

def map_3d(fun, d3):
    return map(lambda d2: map_2d(fun, d2), d3)

# This only works for numpy arrays (which are used by h5py anyway).
def map_nd(fun, arr):
    if ( len(arr.shape) > 1 ):
        return map(lambda ld: map_nd(fun, ld), arr)
    return map(fun, arr)

# Parser
import argparse
parser = argparse.ArgumentParser(description="Convert an HDF5 colour file (plus optional HDF5 cluster-index file) to an HDF5 rock file. Sites with colour >= min and colour <= max (if set) are set to 1, fluid sites are set to 0. The cluster-index file is used to calculate some effective quantities.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-c", "--colourfile", dest="colourfile", help="colour file", type=str, required=True)
parser.add_argument("-o", "--outfile", dest="outfile", help="output file", type=str, required=True)
parser.add_argument("-k", "--clusterfile", dest="clusterfile", help="cluster file", type=str)
parser.add_argument("--colourds", help="colour dataset name", type=str, default="OutArray")
parser.add_argument("--outds", help="output dataset name", type=str, default="OutArray")
parser.add_argument("--clusterds", help="cluster dataset name", type=str, default="OutArray")
parser.add_argument("--gzip", help="gzip compression level", type=int, default=9)
parser.add_argument("--min", help="minimum colour value to turn into rock", type=float, default=0.0)
parser.add_argument("--max", help="maximum colour value to turn into rock", type=float)
options = parser.parse_args()

# Check for fatal errors.
if ( ( options.colourfile == options.outfile ) or ( options.clusterfile == options.outfile ) ):
    print("Please specify different files for input and output.")
    exit(-1)

print("Colour: from dataset '%s' in file '%s'." % ( options.colourds, options.colourfile ) )
if ( options.clusterfile ):
    print("Clusters: from dataset '%s' in file '%s'." % ( options.clusterds, options.clusterfile ) )
print("Output: to dataset '%s' in file '%s'." % ( options.outds, options.outfile ) )
print("Compression level: %d." % ( options.gzip ) )

import h5py
from numpy import *
import numpy as np
from scipy import weave
import time

# Open files.
try:
    f = h5py.File(options.colourfile,"r")
    colarr = np.array(f[options.colourds])
except IOError:
    print("Failed to open input file '%s'." % ( options.colourfile ) )
    exit(-1)

if ( options.clusterfile is not None ):
    try:
        kf = h5py.File(options.clusterfile,"r")
        hkarr = np.array(kf[options.clusterds])
    except IOError:
        print("Failed to open input file '%s'." % ( options.clusterfile ) )
        exit(-1)

try:
    f2 = h5py.File(options.outfile,"w")
except IOError:
    print("Failed to open output file '%s'." % ( options.outfile ) )
    exit(-1)

# Sizes
nx = colarr.shape[0]
ny = colarr.shape[1]
nz = colarr.shape[2]
ns = nx * ny * nz

def colour_to_rock(a, cmin, cmax):
    # C code to speed up the loop by a lot...
    code = """
    const int rows = Na[0];
    const int cols = Na[1];
    const int depth = Na[2];
    for (int i=0; i<rows; i++) {
      for (int j=0; j<cols; j++) {
        for (int k=0; k<depth; k++) {
          int val = (i*cols + j)*depth + k;
            a[val] = ( ( ( a[val] >= cmin )
    """
    # If cmax is not set we don't need and can't use this extra condition.
    if ( cmax is not None ):
        code = code + "&& ( a[val] <= cmax )"
    # Close the loop.
    code = code + """
          ) ? 1 : 0);
        }
      }
    }
    """
    weave.inline(code, ['a', 'cmin', 'cmax'])
    return a

def surface_of_rock(a):
    # C code to speed up the loop by a lot...
    # Look for neighbour sites that are fluid sites.
    code = """
    const int rows = Na[0];
    const int cols = Na[1];
    const int depth = Na[2];
    const int ns = rows*cols*depth;
    int surf = 0;
    for (int i=0; i<rows; i++) {
      for (int j=0; j<cols; j++) {
        for (int k=0; k<depth; k++) {
          int val = (i*cols + j)*depth + k;
          if ( a[val] == 1 ) {
            surf += ( a[(i*cols + j)*depth + ( (k + 1 + depth) % depth )] == 0 );
            surf += ( a[(i*cols + j)*depth + ( (k - 1 + depth) % depth )] == 0 );
            surf += ( a[(i*cols + ( (j + 1 + cols) % cols ) )*depth + k] == 0 );
            surf += ( a[(i*cols + ( (j - 1 + cols) % cols ) )*depth + k] == 0 );
            surf += ( a[(( (i + 1 + rows) % rows)*cols + j)*depth + k] == 0 );
            surf += ( a[(( (i - 1 + rows) % rows)*cols + j)*depth + k] == 0 );
          }
        }
      }
    }
    return_val = surf;
    """
    return weave.inline(code, ['a'])

def calc_hk_fs(a, b, pids):
    # C code to speed up the loop by a lot...
    code = """
    const int rows = Na[0];
    const int cols = Na[1];
    const int depth = Na[2];
    const int npid = Npids[0];
    int fs = 0;
    for (int i=0; i<rows; i++) {
      for (int j=0; j<cols; j++) {
        for (int k=0; k<depth; k++) {
          int val = (i*cols + j)*depth + k;
          for (int p=0; p<npid; p++) {
            if ( a[val] == 0 && b[val] == pids[p] ) {
              fs++;
            }
          }
        }
      }
    }
    return_val = fs;
    """
    return weave.inline(code, ['a', 'b', 'pids'])

def surface_of_rock_hk(a, b, pids):
    # C code to speed up the loop by a lot...
    # Look at neighbour sites and check conditions for both colour array (a) and HK array (b).
    code = """
    const int rows = Na[0];
    const int cols = Na[1];
    const int depth = Na[2];
    const int ns = rows*cols*depth;
    const int npid = Npids[0];
    int surf = 0;
    for (int i=0; i<rows; i++) {
      for (int j=0; j<cols; j++) {
        for (int k=0; k<depth; k++) {
          int val = (i*cols + j)*depth + k;
          if ( a[val] == 1 ) {
            for(int p=0; p<npid; p++) {
              surf += ( ( a[(i*cols + j)*depth + ( (k + 1 + depth) % depth )] == 0 ) &&
                        ( b[(i*cols + j)*depth + ( (k + 1 + depth) % depth )] == pids[p] ) );
              surf += ( ( a[(i*cols + j)*depth + ( (k - 1 + depth) % depth )] == 0 ) &&
                        ( b[(i*cols + j)*depth + ( (k - 1 + depth) % depth )] == pids[p] ) );
              surf += ( ( a[(i*cols + ( (j + 1 + cols) % cols ) )*depth + k] == 0 ) &&
                        ( b[(i*cols + ( (j + 1 + cols) % cols ) )*depth + k] == pids[p] ) );
              surf += ( ( a[(i*cols + ( (j - 1 + cols) % cols ) )*depth + k] == 0 ) &&
                        ( b[(i*cols + ( (j - 1 + cols) % cols ) )*depth + k] == pids[p] ) );
              surf += ( ( a[(( (i + 1 + rows) % rows)*cols + j)*depth + k] == 0 ) && 
                        ( b[(( (i + 1 + rows) % rows)*cols + j)*depth + k] == pids[p] ) );
              surf += ( ( a[(( (i - 1 + rows) % rows)*cols + j)*depth + k] == 0 ) &&
                        ( b[(( (i - 1 + rows) % rows)*cols + j)*depth + k] == pids[p] ) );
            }
          }
        }
      }
    }
    return_val = surf;
    """
    return weave.inline(code, ['a', 'b', 'pids'])

# Convert colour to rock using the upper and lower bounds specified on the command line.
print("Mapping colour... via C weave")
t_start = time.time()
cmapped = colarr.copy() # Don't touch the original array
cmapped = colour_to_rock(cmapped, options.min, options.max)
cmapped = cmapped.astype(int)
rs = float(sum(cmapped))
fs = ns - rs
porosity = fs / ns
t_end = time.time()
print("Took %s seconds" % str(t_end-t_start) )

print("Calculating surface area... via C weave")
t_start = time.time()
surface = surface_of_rock(cmapped)
t_end = time.time()
stv = surface / fs
print("Took %s seconds" % str(t_end-t_start) )

# Simple approximation, which fails for bijels, or any system with a bijel component.
# All values larger than min are now considered to be rock, which correctly excludes
# Pickering droplets, but also erroneously excludes the oil-filled bijel domains.
print("Remapping without upper bound to correct porosity for Pickering emulsions... via C weave")
t_start = time.time()
pmapped = colarr.copy() # Don't touch the original array
pmapped = colour_to_rock(pmapped, options.min, options.max)
pmapped = pmapped.astype(int)
rs_eff = float(sum(pmapped))
fs_eff = ns - rs_eff
porosity_eff = fs_eff / ns
t_end = time.time()
print("Took %s seconds" % str(t_end-t_start))

print("Calculating surface area... via C weave")
t_start = time.time()
surface_eff = surface_of_rock(pmapped)
stv_eff = surface_eff / fs_eff
t_end = time.time()
print("Took %s seconds" % str(t_end-t_start) )

if ( options.clusterfile ):
    # Make a mask having True values for sites that are part of a percolating cluster,
    # by checking which cluster IDs are present at every slice in z-direction.
    # This is not guaranteed to be 100% correct (cluster could theoretically not line up for example:
    #
    # 00100000
    # 00100000
    # 00000100
    # 00000100
    # 00000100
    # 00100000
    # 00100000
    # 00000100
    #
    # For now it seems to be good enough, however, this situation should normally not occur in
    # the bijels / Pickering emulsions of interest (which are also 3D of course).
    print("Finding percolating clusters using HK cluster file.")
    hkslabs = hkarr[:,:]
    percids = list(set(hkslabs[0].flatten())) # Set to remove duplicates.
    for slab in range(1, len(hkslabs)):
        for i in percids:
            if i not in hkslabs[slab].flatten():
                percids.remove(i)

    percids = np.array(percids)
    print("Creating percolation mask from %s..." % str(percids) )
    print("Calculating number of fluid sites... via C weave")
    fs_hk = float(calc_hk_fs(cmapped, hkarr, percids))
    porosity_hk = fs_hk / ns
    print("Calculating surface area... via C weave")
    surface_hk = surface_of_rock_hk(cmapped, hkarr, percids)
    stv_hk = surface_hk / fs_hk

# HDF5 helper functions.
def h5_add_attributes(parent, attr):
    if attr and type(attr) == type({}):
        # attr is a dictionary of attributes
        for k, v in attr.items():
            if v is not None:
                parent.attrs[k] = v

def h5_make_group(parent, name, attr):
    obj = parent.create_group(name)
    h5_add_attributes(obj, attr)
    return obj

# Write output.
print("Writing file...")
if ( options.gzip > 0 ):
    dset = f2.create_dataset(options.outds, data=cmapped, compression="gzip", compression_opts=options.gzip, chunks=True, shuffle=True)
else:
    dset = f2.create_dataset(options.outds, data=cmapped)

# Reproduce the source and input attributes from the colour file.
for a in f[options.colourds].attrs:
    dset.attrs['source'] = a
    dset.attrs['input'] = f[options.colourds].attrs[a]

# Add extra attributes.
dset.attrs['porosity'] = porosity
dset.attrs['porosity_eff'] = porosity_eff
dset.attrs['porosity_hk'] = porosity_hk

dset.attrs['surface'] = surface
dset.attrs['surface_eff'] = surface_eff
dset.attrs['surface_hk'] = surface_hk

dset.attrs['stv'] = stv
dset.attrs['stv_eff'] = stv_eff
dset.attrs['stv_hk'] = stv_hk

# Save the options used to call the script as well.
h5_make_group(f2, "options", vars(options))

f.close()
f2.close()

