#!/usr/bin/env python
# -*- coding: utf-8 -*-

###########################
# Prints out the position data of the interface along the  direction. The position data is rescaled as described in Eddi Winkels Snoeijer 2013 10.1103/PhysRevLett.111.144502
# The interface is simply the colorfield value between xmin and xmax; the default has proven to be usefull for every system so far.
###########################

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
parser = argparse.ArgumentParser(description="Convert an HDF5 colour file (plus optional HDF5 cluster-index file) to an HDF5 rock file. Sites with colour >= min and colour <= max (if set) are set to rvalue, particle sites are set to pvalue. The cluster-index file is used to calculate some effective quantities.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-c", "--colourfile", dest="colourfile", help="colour file", type=str, required=True)
parser.add_argument("--colourds", help="colour dataset name", type=str, default="OutArray")
parser.add_argument("--clusterds", help="cluster dataset name", type=str, default="OutArray")
parser.add_argument("--gzip", help="gzip compression level", type=int, default=9)
parser.add_argument("--min", help="minimum colour value to turn into rock", type=float, default=-0.01)
parser.add_argument("--max", help="maximum colour value to turn into rock", type=float, default=0.01)
parser.add_argument("--rnot", help="Initial droplet radius", type=float, default=1.00)
options = parser.parse_args()

print("Colour: from dataset '%s' in file '%s'." % ( options.colourds, options.colourfile ) )
print("Compression level: %d." % ( options.gzip ) )

import h5py
from numpy import *
import numpy as np
import time

# Open files.
try:
    f = h5py.File(options.colourfile,"r")
    colarr = np.array(f[options.colourds])
except IOError:
    print("Failed to open input file '%s'." % ( options.colourfile ) )
    exit(-1)


# Sizes
nz = colarr.shape[0]
ny = colarr.shape[1]
nx = colarr.shape[2]
ns = nx * ny * nz

def colour_to_rock(a, cmin, cmax):
    # back to python code, since the inline C just won't work, for unknown reasons.
    # find h_0
    for hnoty in range (0,ny):
 #       print("Thing ",hnoty,"is",a[floor(nx/2)][hnoty][floor(nz/2)])
        if a[floor(nz/2)][hnoty][floor(nx/2)]>1:
            break
    #print("the largest thingy is hnot=",hnoty,"==", a[floor(nx/2)][hnoty][floor(nz/2)],"hnot-1=",hnoty-1,"==",a[floor(nx/2)][hnoty-1][floor(nz/2)] )
    hnot=ny-(hnoty+((1-a[floor(nz/2)][hnoty-1][floor(nx/2)])/a[floor(nz/2)][hnoty][floor(nx/2)]-a[floor(nz/2)][hnoty-1][floor(nx/2)]))
    print("result hnot",hnot)
    f = open('%s.csv' % options.colourfile,'w')
    f.write('# hnot %s \t rnot %s \n'% (hnot,options.rnot))
    for fx in range (0,nx):
        for fy in range (2,ny-1):
            if a[floor(nz/2)][fy][fx]>cmin:
                if a[floor(nz/2)][fy][fx]<cmax:
                    temp=fx,"\t",ny-fy,"\n"
                    tempx=float(fx-nx/2)*options.rnot/pow(hnot,2)
                    print('Pos %s\t%s'% (tempx,(ny-fy)/hnot))
                    f.write('%s\t\t%s\n'% (tempx,(ny-fy)/hnot))
    print("system size",nx,ny,nz)
    f.close()
                   
#Go go go
t_start = time.time()
cmapped = colarr.copy() # Don't touch the original array
cmapped = colour_to_rock(cmapped, options.min, options.max)
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

