#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Parser
import argparse
parser = argparse.ArgumentParser(description="Convert a rock file by switching ones and zeroes.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-i", "--infile", dest="infile", help="input file", type=str, required=True)
parser.add_argument("-o", "--outfile", dest="outfile", help="output file", type=str, required=True)
parser.add_argument("--inds", help="input dataset name", type=str, default="OutArray")
parser.add_argument("--outds", help="output dataset name", type=str, default="OutArray")
parser.add_argument("--gzip", help="gzip compression level", type=int, default=9)
options = parser.parse_args()

# Check for fatal errors.
if ( options.infile == options.outfile ):
    print("Please specify different files for input and output.")
    exit(-1)

print("Input: from dataset '%s' in file '%s'." % ( options.inds, options.infile ) )
print("Output: to dataset '%s' in file '%s'." % ( options.outds, options.outfile ) )
print("Compression level: %d." % ( options.gzip ) )

import h5py
from numpy import *
import numpy as np
from scipy import weave
import time

# Open files.
try:
    f = h5py.File(options.infile,"r")
    colarr = np.array(f[options.inds])
except IOError:
    print("Failed to open input file '%s'." % ( options.infile ) )
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

def invert_rock(a):
    # C code to speed up the loop by a lot...
    code = """
    const int rows = Na[0];
    const int cols = Na[1];
    const int depth = Na[2];
    for (int i=0; i<rows; i++) {
      for (int j=0; j<cols; j++) {
        for (int k=0; k<depth; k++) {
          int val = (i*cols + j)*depth + k;
            a[val] = 1 - a[val];
        }
      }
    }
    """
    weave.inline(code, ['a'])
    return a

# Convert colour to rock using the upper and lower bounds specified on the command line.
print("Mapping colour... via C weave")
t_start = time.time()
cmapped = colarr.copy() # Don't touch the original array
cmapped = invert_rock(cmapped)
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
    dset = f2.create_dataset(options.outds, data=cmapped, compression="gzip", compression_opts=options.gzip, chunks=True, shuffle=True)
else:
    dset = f2.create_dataset(options.outds, data=cmapped)

# Reproduce the input file from the rock file.
dset.attrs['input'] = f[options.inds].attrs['input']

f.close()
f2.close()

