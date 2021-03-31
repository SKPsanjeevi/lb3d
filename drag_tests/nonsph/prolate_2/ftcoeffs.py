import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import glob as gb
import math

Re = 1000
rho = 1.
U = 0.09
R_eq = 25.198
A = math.pi*R_eq**2.
ang = [0, 10, 30, 45, 60, 80, 90]
Cd = np.zeros(len(ang))

for i in range(len(ang)):
	# Wildcard search for filenames
	wildcard = './Re'+str(Re)+'ang'+str(ang[i])+'/Production/md-cfg_out_p*.asc'
	for f in gb.glob(wildcard):
		print(f)
		df = pd.read_csv(f, usecols=[0,3], names=['timestep', 'Fd'],
		 header=None, delim_whitespace=True)
	df['Cd'] = df['Fd']/(0.5*rho*A*U**2.)	
	# Compute drag for last N*dt steps
	N = 1000  # dt = 10; So last 10000 timesteps
	Cd[i] = df['Cd'].tail(1000).mean(axis=0)
	print(Cd[i])
	df.loc[1000:].plot(x='timestep', y='Cd')
	plt.show()

df2 = pd.DataFrame({'Angle':ang, 'Cd':Cd})
fname = 'Re'+str(Re)+'_coeffs.txt'
df2.to_csv(fname, index=False, sep='\t', float_format='%.3f')
