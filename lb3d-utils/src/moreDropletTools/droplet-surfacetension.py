#!/usr/bin/python

import ConfigParser
import glob
from optparse import OptionParser
import os
import math
import re

fpath = "../results/"
fname_config = "droplet-surfacetension.cfg"
fname_raw = "droplet-surfacetension.txt"
t_min = 0
t_max = -1

parser = OptionParser()
parser.add_option("--fpath", type="string", dest="fpath", default=fpath)
parser.add_option("--fname_config", type="string", dest="fname_config",default=fname_config)
parser.add_option("--fname_raw", type="string", dest="fname_raw", default=fname_raw)
parser.add_option("--t_min", type="int", dest="t_min", default=t_min)
parser.add_option("--t_max", type="int", dest="t_max", default=t_max)
(options, args) = parser.parse_args()

if (options.t_max < options.t_min and (options.t_max > 0 )):
    print "ERROR: t_max < t_min"
    exit()

# Concatenate full paths
fname_config = options.fpath + options.fname_config
fname_raw = options.fpath + options.fname_raw

def read_file_columns(fname_raw):
    f = open(fname_raw,"r")
    for line in f:
        line = re.sub(r"\n","",line)
        match = re.match(r"^#\?",line)
        # print line
        if match is not None:
            line = re.sub(r"[ \n]+"," ",line)
            headers = line.split(' ')[1:]
            f.close()
            return headers
    f.close()

def calculate_stdev(data):
    n = 0
    Sum = 0
    Sum_sqr = 0
 
    for x in data:
        n = n + 1
        Sum = Sum + x
        Sum_sqr = Sum_sqr + x*x
 
    mean = Sum/n

    if ( n == 1):
        variance = 0
    else:
        variance = (Sum_sqr - Sum*mean)/(n - 1)

    if (variance < 0):
        variance = 0

    return mean, math.sqrt(variance)


# Set up empty path for the new surfacetension file

if not ( os.path.isdir(fpath) ):
    os.mkdir(fpath)

if ( os.path.isfile(fname_raw) ):
    os.remove(fname_raw)

# Create raw surfacetension file

header = True

for fname_od in sorted(glob.glob('od*.h5')):
    # print fname_od

    fname_wd = re.sub("^od","wd",fname_od)
    # print fname_wd
    exists_wd = os.path.isfile(fname_wd)

    fname_sur = re.sub("^od","sur",fname_od)
    # print fname_sur
    exists_sur = os.path.isfile(fname_sur)

    if ( exists_wd ):
        time_match = re.match(r".*_t(\d{8})-.*",fname_od)
        time = time_match.group(1)

        process = "droplet-surfacetension -o " + fname_od + " -w " + fname_wd + " -t " + time

        if ( exists_sur ):
            process += " -s " + fname_sur
        if ( header ):
            process += " -h "
            header = False

        process += " | tee -a " + fname_raw
        
        # print process
        os.system(process)

    else:
        print "WARNING: No matching wd file for od file " + fname_od + " , skipping"

# Get positions of the needed data columns

headers = read_file_columns(fname_raw)
pos_t = headers.index("t")
pos_R = headers.index("R")
pos_DPc = headers.index("DPc")
pos_sigma = headers.index("sigma")
pos_rho = headers.index("rho_m")

# Fill data arrays

R = []
sigma = []
DPc = []
rho = []

f = open(fname_raw,"r")
for line in f:
    line = re.sub(r"\n","",line)
    match = re.match(r"^#",line)
    if match is None:
        line = re.sub(r"[ \n]+"," ",line)
        data = line.split(' ')
        t = int(data[pos_t])
        if ( (t >= options.t_min and t <= options.t_max ) or ( t>= options.t_min and options.t_max < 0 ) or ( options.t_min < 0 and t <= options.t_max ) ):
            R.append(float(data[pos_R]))
            sigma.append(float(data[pos_sigma]))
            DPc.append(float(data[pos_DPc]))
            rho.append(float(data[pos_rho]))

f.close()

# Get statistics

R_mean, R_stdev = calculate_stdev(R)
# print R_mean, R_stdev

sigma_mean, sigma_stdev = calculate_stdev(sigma)
# print sigma_mean, sigma_stdev

DPc_mean, DPc_stdev = calculate_stdev(DPc)
# print DPc_mean, DPc_stdev

rho_mean, rho_stdev = calculate_stdev(rho)
# print rho_mean, rho_stdev

# Write processed data to file

config = ConfigParser.RawConfigParser()
sname = "Errors"
config.add_section(sname)
config.set(sname, "R", R_stdev)
config.set(sname, "sigma", sigma_stdev)
config.set(sname, "DPc", DPc_stdev)
config.set(sname, "rho", rho_stdev)
sname = "Averages"
config.add_section(sname)
config.set(sname, "R", R_mean)
config.set(sname, "sigma", sigma_mean)
config.set(sname, "DPc", DPc_mean)
config.set(sname, "rho", rho_mean)
sname = "Ranges"
config.add_section(sname)
config.set(sname, "t_min", options.t_min)
config.set(sname, "t_max", options.t_max)

configfile = open(fname_config, 'wb')
config.write(configfile)


