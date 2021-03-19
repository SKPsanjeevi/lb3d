#!/usr/bin/env python
# -*- coding: utf-8 -*-

import h5py
import os
import re
import subprocess
import tempfile

# Parser.
import argparse
parser = argparse.ArgumentParser(description="Wrapper for the structure factor calculation. This script processes a single HDF5 file, taking care of calling the binary with the required arguments and providing a temporary input file based on the HDF5 metadata of the input file. It also assumes that only a single .siz file will be created (i.e. only one timestep will be considered), and will prepend a header to that file.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-i", "--infile", help="input file (HDF5), relative path only", type=str, required=True)
parser.add_argument("--inds", help="dataset name", type=str, default="OutArray")
parser.add_argument("-o", "--outputpath", help="output path, relative path only", type=str, default="./sf")
parser.add_argument("-b", "--binary", help="binary executable", type=str, default="./struct_hdf")
header_grp = parser.add_mutually_exclusive_group()
# Order matters! Switch default values by swapping these two lines.
header_grp.add_argument('--no-header', action='store_false', dest='header', help='do not add header information')
header_grp.add_argument('--header', action='store_true', dest='header', help='add header information')
options = parser.parse_args()

# infile and outputpath should be relative.
if ( options.outputpath[0] == "/" or options.infile[0] == "/" ):
    print("ERROR: --infile and --outputpath should be relative - this is a restriction of the binary.")
    exit(-1)

# Extract information from filename.
(path, fn) = os.path.split(options.infile)
if ( path == "" ):
    path = "./"
m = re.match(r"(?P<ftype>.*)_(?P<gr_out_file>.*)_t(?P<ts>[0-9]*).h5", fn)
if m:
    ftype = m.group('ftype')
    ts = m.group('ts')
    gr_out_file = m.group('gr_out_file')
else:
    print("ERROR: Failed to parse filename '%s'. Please note, the unique ID should be stripped from the filename - this is a restriction of the binary." % ( options.infile ) )
    exit(-1)
print("Found input file of type '%s' at timestep '%d' in folder '%s'." % ( ftype, int(ts), path ) )

# If necessary, create output directory.
if ( options.outputpath == "" ):
    options.outputpath = "./"

if not os.path.exists(options.outputpath):
    print("Creating output path '%s'." % options.outputpath)
    os.makedirs(options.outputpath)

# If necessary, add header information to .siz file.
if ( options.header ):
    sizfn = "%s/%s_%s.siz" % ( options.outputpath, ftype, gr_out_file )
    if os.path.isfile(sizfn):
        has_header = False
        for l in open(sizfn):
            if ( l[0:2] == "#?" ):
                has_header = True
                break
        if not has_header:
            print("Adding header information to '%s'." % sizfn)
            sizf = open(sizfn, 'w')
            sizf.write("#? ts k1 errk1 k2 errk2 R1 errR1 R2 errR2 Rx errRx Ry errRy Rz errRz min max avg\n")
            sizf.close()
    else:
        print("Writing header information to '%s'." % sizfn)
        sizf = open(sizfn, 'w')
        sizf.write("#? ts k1 errk1 k2 errk2 R1 errR1 R2 errR2 Rx errRx Ry errRy Rz errRz min max avg\n")
        sizf.close()

# Open HDF5 file.
try:
    f = h5py.File(options.infile,"r")
except IOError:
    print("ERROR: Failed to open HDF5 file '%s'." % ( options.infile ) )
    exit(-1)

# Write temporary input-file, overriding n_sci_start, folder, and gr_out_file parameters.
tmpf, tmpfn = tempfile.mkstemp()
try:
    for a in f[options.inds].attrs:
        for l in f[options.inds].attrs[a]:
            # Override parameters.
            if "n_sci_start" in l: l = "n_sci_start = %d" % int(ts)
            if "folder" in l:      l = "folder = '%s'" % path
            if "gr_out_file" in l: l = "gr_out_file = '%s'" % gr_out_file # To fix broken input files.
            os.write(tmpf, l+"\n")
    os.close(tmpf)
    print("Calling binary '%s'" % options.binary)
    subprocess.call([options.binary, "-file", ftype, "-nt", ts, "-meas", ts, "-inpfile", tmpfn, "-data", os.getcwd(), "-res", options.outputpath])
finally:
    os.remove(tmpfn)

