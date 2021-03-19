#!/usr/bin/python

from optparse import OptionParser
import re

rockcolour = 5.0
domaincolour = 0.0
scale = 1.0

usage = "Usage: %prog [options] file_in file_out"
parser = OptionParser(usage=usage)
parser.add_option("-f", type="float", dest="fluidcolour", default=domaincolour, help="colour of fluid sites (default 0.0).")
parser.add_option("-r", type="float", dest="rockcolour", default=rockcolour, help="colour of rock sites (default 5.0).")
parser.add_option("-s", type="float", dest="scale", default=scale, help="upscale the system with this factor (default 1.0).")
(options, args) = parser.parse_args()

fname_in = args[0]
fname_out = args[1]

scale = options.scale

f = open(fname_in,"r")
tf = open(fname_out,"w")
for ln, line in enumerate(f):
    line = re.sub(r"\n","",line)
    line = re.sub(r"\r","",line)
    line = re.sub("\t"," ",line)
    match = re.match(r"^#\?",line)
    #print line
    if match is None:
        if (ln == 0):
            data = line.split(' ')[1:]
            # First line contains system size and such
            np = int(data[0])
            t = int(data[1])
            x0 = float(data[2])
            y0 = float(data[3])
            z0 = float(data[4])
            x1 = float(data[5])
            y1 = float(data[6])
            z1 = float(data[7])
            # Rescale to get size 
            nx = int((x1 - x0)*scale)
            ny = int((y1 - y0)*scale)
            nz = int((z1 - z0)*scale)
            s = "%d\n%d\n%d\n\n%f\n\n" % (nx, ny, nz, options.fluidcolour)
            tf.write(s)
            print "System size: ", nx, ny, nz

        else:
            data = line.split(' ')[:]
            # All subsequent lines contain sphere data
            xf = float(data[0])
            yf = float(data[1])
            zf = float(data[2])
            rf = float(data[6])
            x = (xf - x0)*scale
            y = (yf - y0)*scale
            z = (zf - z0)*scale
            r = rf*scale
            s = "sphere\n%f %f %f\n%f\n%f\n\n" % (x, y, z, r, options.rockcolour)
            tf.write(s)

tf.write("end\n\n")

f.close()
tf.close()
