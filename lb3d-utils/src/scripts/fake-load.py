#!/usr/bin/env python

import sys

from multiprocessing import Pool

def f(x):
    # Put any cpu (only) consuming operation here. I have given 1 below -
    while True:
        x * x

# decide how many cpus you need to load with.
no_of_cpu_to_be_consumed = int(sys.argv[1])

p = Pool(processes=no_of_cpu_to_be_consumed)
p.map(f, range(no_of_cpu_to_be_consumed))

