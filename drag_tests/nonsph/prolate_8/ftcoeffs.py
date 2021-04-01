import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import glob as gb
import math

Re = 1000
rho = 1.
U = 0.09
R_eq = 30. #23.811
A = math.pi*R_eq**2.
ang = [0, 10, 30, 45, 60, 80, 90]
Cd = np.zeros(len(ang))
Cl = np.zeros(len(ang))
Ct = np.zeros(len(ang))

for i in range(len(ang)):
	# Wildcard search for filenames
	wildcard = './Re'+str(Re)+'ang'+str(ang[i])+'/Production/md-cfg_out_p*.asc'
	for f in gb.glob(wildcard):
		print(f)
		df = pd.read_csv(f, usecols=[0, 2, 3, 4], names=['timestep', 'Fl', 'Fd', 'Tp'],
		 header=None, delim_whitespace=True)
	df['Cd'] = df['Fd']/(0.5*rho*A*U**2.)
	df['Cl'] = df['Fl']/(0.5*rho*A*U**2.)
	df['Ct'] = df['Tp']/(0.5*rho*A*R_eq*U**2.)
	# Compute drag for last N*dt steps
	N = 2000  # dt = 10; So last 20000 timesteps
	Cd[i] = df['Cd'].tail(N).mean(axis=0)
	Cl[i] = df['Cl'].tail(N).mean(axis=0)
	Ct[i] = df['Ct'].tail(N).mean(axis=0)
	np.set_printoptions(precision=3)
	print(Cd[i])
	df.loc[1000:].plot(x='timestep', y='Cd')
	plt.show()

df2 = pd.DataFrame({'Angle':ang, 'Cd':Cd, 'Cl':Cl, 'Ct':Ct})
fname = 'Re'+str(Re)+'_coeffs.txt'
df2.to_csv(fname, index=False, sep='\t', float_format='%.3f')
