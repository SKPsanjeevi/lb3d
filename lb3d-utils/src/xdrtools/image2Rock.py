#!/usr/bin/env python
# -*- coding: utf-8 -*-

#####################################
# Convert a black and white picture to a h5 rock file.
# Written by Dennis Hessling
#####################################

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
import re
parser = argparse.ArgumentParser(description="Convert a black and white image to a rock file. Height of the image and the respective colours may be given. The size of the rock file is given by the size of the image file, the remaining direction is defined through the remainingsize parameter.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-c", "--colourfile", dest="colourfile", help="Imagefile to use for color", type=str, required=True)
parser.add_argument("-o", "--outfile", dest="outfile", help="Output file", type=str, required=True)
parser.add_argument("--colourthreshhold", help="Threshhold value of black colour in image, needed for eventual grayscale conversion or compression artifacts", type=float, default=5)
parser.add_argument("--gzip", help="gzip compression level", type=int, default=9)
parser.add_argument("--colourds", help="colour dataset name", type=str, default="OutArray")
parser.add_argument("--outds", help="output dataset name", type=str, default="OutArray")
parser.add_argument("--bcolour", help="Colour of black lattice sites", type=float, default=0)
parser.add_argument("--wcolour", help="Colour of white lattice sites", type=float, default=6)
parser.add_argument("--bgcolour", help="Colour of the empty void, 0 is a fluid site", type=int, default=0)
parser.add_argument("--remainingsize", help="Distance to grow normal to the picture", type=int, default=64)
parser.add_argument("--colourthick", help="Thicknes of the layer from the picutre", type=int, default=2)
parser.add_argument("--colourfloat", help="Height the layer from the picture is floating", type=int, default=0)
options = parser.parse_args()

print("Colour: from dataset '%s' in file '%s'." % ( options.colourds, options.colourfile ) )
print("Output: to dataset '%s' in file '%s'." % ( options.outds, options.outfile ) )
print("Compression level: %d." % ( options.gzip ) )

import h5py
from numpy import *
import numpy as np
import time
from scipy import misc

# Open files.
try:
    # f = h5py.File(options.colourfile,"r")
    # colarr = np.array(f[options.colourds])
    readimage=misc.imread(options.colourfile)
    readimage.shape
except IOError:
    print("Failed to open input file '%s'." % ( options.colourfile ) )
    exit(-1)
 
try:
    f2 = h5py.File(options.outfile,"w")
except IOError:
    print("Failed to open output file '%s'." % ( options.outfile ) ) 
    exit(-1)

# Flip the image to have it right as the bottom
readimage = np.fliplr(readimage)

# Sizes
nx = readimage.shape[0]
ny = readimage.shape[1]
nz = options.remainingsize

print("Mapping colour...")
t_start = time.time()
cmapped = np.ones((nx,ny,nz))*options.bgcolour

for inx in range (0,nx):
   for iny in range (0,ny):
      for inz in range(options.colourfloat,options.colourfloat+options.colourthick-1):
         if readimage.ndim<=3:
            # This is a rgb image, threat the first channel only
            if readimage[inx,iny,0]<=options.colourthreshhold:
               cmapped[inx,iny,inz]=options.wcolour
            else:
               cmapped[inx,iny,inz]=options.bcolour
         else:
            # This is a bw image, all is fine
            if readimage[inx,iny]<=options.colourthreshhold:
               cmapped[inx,iny,inz]=options.wcolour
            else:
               cmapped[inx,iny,inz]=options.bcolour

cmapped=swapaxes(cmapped,1,2)
nz = cmapped.shape[0]
ny = cmapped.shape[1]
nx = cmapped.shape[2]
print("system size after",nx,ny,nz)

t_end = time.time()
print("Took %s seconds" % str(t_end-t_start) )
print(cmapped)

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

#!!!WARNING
#this one is dangerous, if this is a fluid file
dset.attrs['shear_sum'] = 0

# Save the options used to call the script as well.
h5_make_group(f2, "options", vars(options))

f2.close()
