#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Parser
import argparse
parser = argparse.ArgumentParser(description="Merge two HDF5 files with binary operation.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-i", "--infiles", dest="infiles", help="input files", nargs="+", type=str, required=True)
parser.add_argument("-o", "--outfile", dest="outfile", help="output file", type=str, required=True)
parser.add_argument("--inds", help="input dataset name", type=str, default="OutArray")
parser.add_argument("--outds", help="output dataset name", type=str, default="OutArray")
parser.add_argument("--gzip", help="gzip compression level", type=int, default=9)
options = parser.parse_args()

# Check for fatal errors.
if ( options.outfile  in options.infiles ):
    print("Please specify different files for input and output.")
    exit(-1)

if ( len(options.infiles) != 2) :
    print("Please specify exactly two input files.")
    exit(-1)

print("Input 1: from dataset '%s' in file '%s'." % ( options.inds, options.infiles[0] ) )
print("Input 2: from dataset '%s' in file '%s'." % ( options.inds, options.infiles[1] ) )
print("Output: to dataset '%s' in file '%s'." % ( options.outds, options.outfile ) )
print("Compression level: %d." % ( options.gzip ) )

import h5py
from numpy import *
import numpy as np
from scipy import weave
import time

# Open files.
try:
    fin1 = h5py.File(options.infiles[0],"r")
    finarr1 = np.array(fin1[options.inds])
except IOError:
    print("Failed to open input file '%s'." % ( options.infiles[0] ) )
    exit(-1)

try:
    fin2 = h5py.File(options.infiles[1],"r")
    finarr2 = np.array(fin2[options.inds])
except IOError:
    print("Failed to open input file '%s'." % ( options.infiles[1] ) )
    exit(-1)

if ( finarr1.shape != finarr2.shape ):
    print("Error: dataset shape mismatch.")
    exit(-1)

try:
    fout = h5py.File(options.outfile,"w")
except IOError:
    print("Failed to open output file '%s'." % ( options.outfile ) )
    exit(-1)

# Sizes
nx = finarr1.shape[0]
ny = finarr1.shape[1]
nz = finarr1.shape[2]
ns = nx * ny * nz

def merge_data(a, b):
    # C code to speed up the loop by a lot...
    code = """
    const int rows = Na[0];
    const int cols = Na[1];
    const int depth = Na[2];
    for (int i=0; i<rows; i++) {
      for (int j=0; j<cols; j++) {
        for (int k=0; k<depth; k++) {
          int val = (i*cols + j)*depth + k;
          if ( a[val] != 0 || b[val] != 0 ) {
            a[val] = 1;
          }
        }
      }
    }
    return_val = 0;
    """
    weave.inline(code, ['a', 'b'])
    return a

print("Merging files... via C weave")
t_start = time.time()
cmapped = finarr1.copy() # Don't touch the original array
cmapped = merge_data(cmapped, finarr2)
cmapped = cmapped.astype(int)
t_end = time.time()
print("Took %s seconds" % str(t_end-t_start) )

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
    dset = fout.create_dataset(options.outds, data=cmapped, compression="gzip", compression_opts=options.gzip, chunks=True, shuffle=True)
else:
    dset = fout.create_dataset(options.outds, data=cmapped)

# Reproduce the input file from the rock file.
dset.attrs['input'] = fin1[options.inds].attrs['input']

fin1.close()
fin2.close()
fout.close()

