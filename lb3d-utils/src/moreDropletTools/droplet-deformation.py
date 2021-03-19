#!/usr/bin/python

import ConfigParser
import glob
from optparse import OptionParser
import os
import math
import re

fpath = "../results/"
fname_config = "droplet-deformation.cfg"
fname_raw = "droplet-deformation.txt"
t_min = 0
t_max = -1

parser = OptionParser()
parser.add_option("--cutoff", type="string", dest="cutoff",default="0.0")
parser.add_option("--parameters", type="string", dest="fname_parameters",default="droplet-surfacetension.cfg")
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
fname_parameters = options.fpath + options.fname_parameters
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

# Write processed data to file


config = ConfigParser.RawConfigParser()
config.read(fname_parameters)

sname = 'Errors'
R_stdev = config.get(sname, "R")
sigma_stdev = config.get(sname, "sigma")
rho_stdev = config.get(sname, "rho")

sname = "Averages"
R = config.get(sname, "R")
sigma = config.get(sname, "sigma")
rho = config.get(sname, "rho")

# Create raw surfacetension file

header = True

for fname_od in sorted(glob.glob('od*.h5')):
    # print fname_od

    fname_vel = re.sub("^od","velz",fname_od)
    # print fname_sur
    exists_vel = os.path.isfile(fname_vel)

    time_match = re.match(r".*_t(\d{8})-.*",fname_od)
    time = time_match.group(1)

    process  = "droplet-deformation -o " + fname_od + " -t " + time + " -c " + options.cutoff
    process += " -radius " + R + " -Dradius " + R_stdev
    process += " -rho " + rho + " -Drho " + rho_stdev
    process += " -sigma " + sigma + " -Dsigma " + sigma_stdev

    if ( exists_vel ):
        process += " -v " + fname_vel
    if ( header ):
        process += " -h "
        header = False

    process += " | tee -a " + fname_raw
        
    # print process
    os.system(process)

# Get positions of the needed data columns

headers = read_file_columns(fname_raw)
pos_t = headers.index("t")
pos_D = headers.index("D")
pos_theta = headers.index("theta")
pos_Ca = headers.index("Ca")
pos_Ca_gdavg = headers.index("Ca_gdavg")
pos_Ca_gdmax = headers.index("Ca_gdmax")
pos_Re = headers.index("Re")
pos_chi = headers.index("chi")

# Fill data arrays

D = []
theta = []
Ca = []
Ca_gdavg = []
Ca_gdmax = []
Re = []
chi = []

f = open(fname_raw,"r")
for line in f:
    line = re.sub(r"\n","",line)
    match = re.match(r"^#",line)
    if match is None:
        line = re.sub(r"[ \n]+"," ",line)
        data = line.split(' ')
        t = int(data[pos_t])
        if ( (t >= options.t_min and t <= options.t_max ) or ( t>= options.t_min and options.t_max < 0 ) or ( options.t_min < 0 and t <= options.t_max ) ):
            D.append(float(data[pos_D]))
            theta.append(float(data[pos_theta]))
            Ca.append(float(data[pos_Ca]))
            Ca_gdavg.append(float(data[pos_Ca_gdavg]))
            Ca_gdmax.append(float(data[pos_Ca_gdmax]))
            Re.append(float(data[pos_Re]))
            chi.append(float(data[pos_chi]))

f.close()

Re_gdavg = []
Re_gdmax = []
for r in Re:
    Re_gdavg.append(r * ( Ca_gdavg[0] / Ca[0] ))
    Re_gdmax.append(r * ( Ca_gdmax[0] / Ca[0] ))

# Get statistics

D_mean, D_stdev = calculate_stdev(D)
# print D_mean, D_stdev

theta_mean, theta_stdev = calculate_stdev(theta)
# print D_mean, D_stdev

Ca_mean, Ca_stdev = calculate_stdev(Ca)
# print Ca_mean, Ca_stdev

Ca_gdavg_mean, Ca_gdavg_stdev = calculate_stdev(Ca_gdavg)
# print Ca_gdavg_mean, Ca_gdavg_stdev

Ca_gdmax_mean, Ca_gdmax_stdev = calculate_stdev(Ca_gdmax)
# print Ca_gdmax_mean, Ca_gdmax_stdev

Re_mean, Re_stdev = calculate_stdev(Re)
# print Ca_mean, Ca_stdev

Re_gdavg_mean, Re_gdavg_stdev = calculate_stdev(Re_gdavg)
# print Ca_gdavg_mean, Ca_gdavg_stdev

Re_gdmax_mean, Re_gdmax_stdev = calculate_stdev(Re_gdmax)
# print Ca_gdmax_mean, Ca_gdmax_stdev

chi_mean, chi_stdev = calculate_stdev(chi)
# print Ca_chi_mean, Ca_chi_stdev

config = ConfigParser.RawConfigParser()
sname = "Errors"
config.add_section(sname)
config.set(sname, "D", D_stdev)
config.set(sname, "theta", theta_stdev)
config.set(sname, "Ca", Ca_stdev)
config.set(sname, "Ca_gdavg", Ca_gdavg_stdev)
config.set(sname, "Ca_gdmax", Ca_gdmax_stdev)
config.set(sname, "Re", Re_stdev)
config.set(sname, "Re_gdavg", Re_gdavg_stdev)
config.set(sname, "Re_gdmax", Re_gdmax_stdev)
config.set(sname, "chi", chi_stdev)
sname = "Averages"
config.add_section(sname)
config.set(sname, "D", D_mean)
config.set(sname, "theta", theta_mean)
config.set(sname, "Ca", Ca_mean)
config.set(sname, "Ca_gdavg", Ca_gdavg_mean)
config.set(sname, "Ca_gdmax", Ca_gdmax_mean)
config.set(sname, "Re", Re_mean)
config.set(sname, "Re_gdavg", Re_gdavg_mean)
config.set(sname, "Re_gdmax", Re_gdmax_mean)
config.set(sname, "chi", chi_mean)
sname = "Ranges"
config.add_section(sname)
config.set(sname, "t_min", options.t_min)
config.set(sname, "t_max", options.t_max)

configfile = open(fname_config, 'wb')
config.write(configfile)

