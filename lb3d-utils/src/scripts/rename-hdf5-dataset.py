#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Parser
from optparse import OptionParser, OptionGroup
parser = OptionParser(description='HDF5 converter helper')
parser.add_option("-i", "--infile", dest="infile", help="Input file", type="string", default="")
parser.add_option("-o", "--outfile", dest="outfile", help="Output file", type="string", default="")
parser.add_option("--inds", help="Input dataset name", type="string", default="OutArray")
parser.add_option("--outds", help="Output dataset name", type="string", default="OutArray")
parser.add_option("--gzip", help="GZIP compression level", type="int", default=0)
parser.add_option("--xmin", help="Minimum x-coordinate", type="int", default=0)
parser.add_option("--xmax", help="Maximum x-coordinate", type="int", default=-1)
parser.add_option("--ymin", help="Minimum y-coordinate", type="int", default=0)
parser.add_option("--ymax", help="Maximum y-coordinate", type="int", default=-1)
parser.add_option("--zmin", help="Minimum z-coordinate", type="int", default=0)
parser.add_option("--zmax", help="Maximum z-coordinate", type="int", default=-1)
(options, args) = parser.parse_args()

if ( not ( options.infile and options.outfile ) ):
    print "Please specify both input and output filenames."
    exit(-1)

if ( options.infile == options.outfile ):
    print "Please specify different files for input and output."
    exit(-1)

print "Input:  dataset '" + options.inds  +"' in file '" + options.infile  + "'."
print "Output: dataset '" + options.outds +"' in file '" + options.outfile + "'."
print "Compression level: " + str(options.gzip)

import h5py

# Read files
try:
    f = h5py.File(options.infile,'r')
except IOError:
    print "Failed to open input file '" + options.infile + "'."
    exit(-1)

# Set ranges
xsize = len(f[options.inds][0][0])
ysize = len(f[options.inds][0])
zsize = len(f[options.inds])

if ( options.xmax < 0 ):
    options.xmax = xsize
else:
    options.xmax = min(options.xmax, xsize)

if ( options.ymax < 0 ):
    options.ymax = ysize
else:
    options.ymax = min(options.ymax, ysize)

if ( options.zmax < 0 ):
    options.zmax = zsize
else:
    options.zmax = min(options.zmax, zsize)

options.xmin = max(options.xmin, 0)
options.ymin = max(options.ymin, 0)
options.zmin = max(options.zmin, 0)

# Write output
try:
    f2 = h5py.File(options.outfile,'w')
except IOError:
    print "Failed to open output file '" + options.outfile + "'."
    exit(-1)

if ( options.gzip > 0 ):
    try:
        dset = f2.create_dataset(options.outds, data=f[options.inds][options.zmin:options.zmax,options.ymin:options.ymax,options.xmin:options.xmax], compression='gzip', compression_opts=options.gzip, chunks=True, shuffle=True)
    except KeyError:
        print "Failed to open dataset."
        exit(-1)
else:
    try:
        dset = f2.create_dataset(options.outds, data=f[options.inds][options.zmin:options.zmax,options.ymin:options.ymax,options.xmin:options.xmax])
    except KeyError:
        print "Failed to open dataset."
        exit(-1)

f.close()
f2.close()

