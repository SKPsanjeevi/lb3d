#!/usr/bin/env python

import sys
import xdrlib

if len(sys.argv) != 3:
 print "usage: ",sys.argv[0],"<n> <file>"
 print "interprets the <n>th 32bits of <file> as an int and dump it to"
 print "standard output"
 exit -1

nint = int(sys.argv[1])
filename = sys.argv[2]

data = open(filename,"r").read()
u = xdrlib.Unpacker(data)

for i in range(nint):
 dummy = u.unpack_uint()
 if i == nint - 1:
  print dummy
