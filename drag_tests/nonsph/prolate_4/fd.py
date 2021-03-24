import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import glob as gb
import math

rho = 1.
U = 0.09
R_eq = 23.811
A = math.pi*R_eq**2.

# Wildcard search for filenames
wildcard = './Re1000ang30/Production/md-cfg_out_p*.asc'
for f in gb.glob(wildcard):
	print(f)
	df = pd.read_csv(f, usecols=[0,3], names=['timestep', 'Fd'],
	 header=None, delim_whitespace=True)
df['Cd'] = df['Fd']/(0.5*rho*A*U**2.)

# Compute drag for last N*dt steps
N = 2000  # dt = 10; So last 20000 timesteps
print(df['Cd'].tail(2000).mean(axis=0))

df.loc[1000:].plot(x='timestep', y='Cd')
plt.show()
