#!/usr/bin/env python
import os
import math
import sys

if len(sys.argv)!=2:
 print 'usage: '+sys.argv[0]+' <dx>'
 print
 print 'Calculates conversion factors from lattice to physical units for'
 print 'various quantities based on the length conversion factor dx given'
 print 'as argument, choosing the relaxation time tau=1, and matching the'
 print 'resulting lattice viscosity and mass density with the respective'
 print 'values for blood plasma as given by the constants in this script.'
 os._exit(1)

# physical kinematic viscosity and mass density of blood plasma
nu_plasma = 1.09e-6
rho_plasma = 1.03e3

dx = float(sys.argv[1])
dt = dx * dx / (6.0 * nu_plasma)
dm = dx * dx * dx * rho_plasma

print "dx = ", dx, " m"
print "dt = ", dt, " s"
print "dm = ", dm, " kg"
print
print "dE = ", dm*dx*dx/(dt*dt), " J"
print
print "dv = ", dx/dt, " m/s"
print
print "dQ = ", dx*dx*dx/dt, " m^3/s"
print
print "dgammadot = ", 1.0/dt, " 1/s"
print
print "dsigma = dp = ", dm/(dx*dt*dt), " N/m^2"
print
print "dF = ", dm*dx/(dt*dt), " N"
print "df = ", dm/(dx*dx*dt*dt), " N/m^3"
