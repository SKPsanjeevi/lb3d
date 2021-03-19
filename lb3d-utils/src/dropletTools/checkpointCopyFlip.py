#!/usr/bin/env python
# -*- coding: utf-8 -*-

#####################################
# Takes a checkpoint file and copies it next to itself, mirrored in x, stitched in x, around x=0
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
parser = argparse.ArgumentParser(description="Convert an HDF5 colour file (plus optional HDF5 cluster-index file) to an HDF5 rock file. Sites with colour >= min and colour <= max (if set) are set to rvalue, particle sites are set to pvalue. The cluster-index file is used to calculate some effective quantities.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-c", "--colourfile", dest="colourfile", help="colour file", type=str, required=True)
parser.add_argument("-o", "--outfile", dest="outfile", help="output file", type=str, required=True)
parser.add_argument("--gzip", help="gzip compression level", type=int, default=9)
parser.add_argument("--colourds", help="colour dataset name", type=str, default="OutArray")
parser.add_argument("--outds", help="output dataset name", type=str, default="OutArray")
parser.add_argument("--xshift", help="Shift the droplet in x", type=int)
options = parser.parse_args()

print("Colour: from dataset '%s' in file '%s'." % ( options.colourds, options.colourfile ) )
print("Output: to dataset '%s' in file '%s'." % ( options.outds, options.outfile ) )
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
 
try:
    f2 = h5py.File(options.outfile,"w")
except IOError:
    print("Failed to open output file '%s'." % ( options.outfile ) ) 
    exit(-1)


if "colour" in options.colourfile:
    # check how much to fit
    # Sizes
    nz = colarr.shape[0]
    ny = colarr.shape[1]
    nx = colarr.shape[2]
    bp=0
    print("system size",nx,ny,nz)
    for hnotx in range ((nx-1)/2,0,-1):
	bp=bp+1
        print colarr[floor(nz/2)][ny-2][hnotx]
        if colarr[floor(nz/2)][ny-2][hnotx]<0:
            print 'now'
            print colarr[floor(nz/2)][2][hnotx]
            print hnotx
            print ((nx-1)/2)-bp
            #break
else:
    # Sizes
    nz = colarr.shape[0]
    ny = colarr.shape[1]
    nx = colarr.shape[2]
    nq = colarr.shape[3]
    ns = nx * ny * nz
    
    def colour_to_rock(a, cmin, cmax):
        # back to python code, since the inline C just won't work, for unknown reasons.
        # find h_0
        for hnoty in range (ny,0):
     #       print("Thing ",hnoty,"is",a[floor(nx/2)][hnoty][floor(nz/2)])
            if a[floor(nz/2)][hnoty][floor(nx/2)]>0:
                break
        #print("the largest thingy is hnot=",hnoty,"==", a[floor(nx/2)][hnoty][floor(nz/2)],"hnot-1=",hnoty-1,"==",a[floor(nx/2)][hnoty-1][floor(nz/2)] )
        hnot=ny-(hnoty+((1-a[floor(nz/2)][hnoty-1][floor(nx/2)])/a[floor(nz/2)][hnoty][floor(nx/2)]-a[floor(nz/2)][hnoty-1][floor(nx/2)]))
        print("result hnot",hnot,"timestep",re.sub(r".*t([0-9]*)-.*",r"\1",options.colourfile))
        f = open('%s.csv' % options.outfile,'a')
        f.write('%s\t%s\n' % (re.sub(r".*t([0-9]*)-.*",r"\1",options.colourfile),hnot))
        print("system size",nx,ny,nz,nq)
        f.close()
    
    # Convert colour to rock using the upper and lower bounds specified on the command line.
    print("Mapping colour... via C weave")
    t_start = time.time()
    cmapped = colarr.copy() # Don't touch the original array
    print("system size before",nx,ny,nz,nq)
    if options.xshift:
        cmapped=np.roll(cmapped,options.xshift,axis=2)
    cmappedflip = cmapped[:,:,::-1,:]
    cmapped=np.append(cmapped,cmappedflip,axis=2)
    nz = cmapped.shape[0]
    ny = cmapped.shape[1]
    nx = cmapped.shape[2]
    nq = cmapped.shape[3]
    print("system size after",nx,ny,nz,nq)
    
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
    #this one is dangerous
    dset.attrs['shear_sum'] = 0
    
    # Save the options used to call the script as well.
    h5_make_group(f2, "options", vars(options))
    
    f.close()
    f2.close()
    
