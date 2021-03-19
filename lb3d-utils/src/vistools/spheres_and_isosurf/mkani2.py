#!/usr/bin/python
#Beware: This is badly written.

import os
import sys
import spheres_and_isosurf2


h5tovtk="/usr/local/bin/h5tovtk"
#sai="/data/data2/fabian/tools/bin/spheres_and_isosurf.py"
env=os.environ
env['LD_LIBRARY_PATH']='/usr/local64/stow/VTK-5.0.3/lib:/usr/local64/stow/VTK-5.0.3/lib64:'+env['LD_LIBRARY_PATH']
env['PYTHONPATH']='/usr/local64/stow/VTK-5.0.3/lib64/python2.5/site-packages'

fps='5'
contactangle=0
timestamps=[]
for arg in sys.argv[1:]:
	print arg
	if arg.startswith('-fps='):
		fps=arg.split('=')[1]
	else:
		if arg.startswith('-car='):
			contactangle=1
			import contact_angle
			r=int(arg.split('=')[1])
		else:
			timestamps.append(arg)
files=os.listdir(".")
try:
	paramfile=open("mkani_params","r")
	params=[param for param in paramfile]
	paramfile.close()
except IOError:
	params=[]
if len(timestamps):
	chars_to_remove=len(timestamps[0])+4
	os.mkdir(timestamps[0])
	os.chdir(timestamps[0])
	for file in files:
		for timestamp in timestamps:
			if file.endswith(timestamp+".h5"):
				newname=file[:-chars_to_remove]+".h5"
				try:
					os.symlink("../"+file, newname)
				except OSError:
					print "Couldn't create symlink for "+file+", perhaps because we are making a single movie out of several simulations?"
				else:
					if os.spawnl(os.P_WAIT,h5tovtk,h5tovtk,'--',newname)!=0:
						print "error converting "+file+" to vtk"
					os.remove(newname)
			if file.endswith(timestamp+".asc"):
				newname=file[:-(chars_to_remove+1)]+".asc"
				try:
					os.symlink("../"+file, newname)
				except OSError:
					print "Couldn't create symlink for "+file+", perhaps because we are making a single movie out of several simulations?"




if contactangle:
	cafile=open('contact_angle.dat','a')

first=1
for file in sorted(os.listdir(".")):
	if file.endswith(".vtk") and file.startswith("colour"):
		if not contactangle:		
			last=file
			fid="_".join(file.split('_')[1:])[:-4]
			if first:
				first=0
				sai=spheres_and_isosurf2.spheres_and_isosurf(fid,params)
			else:
				sai.update(fid)
				
			montageparams=['/usr/bin/montage','rendered_left_'+fid+'.png','rendered_right_'+fid+'.png','-tile','2x1','-geometry','+0+0','-title',fid,'rendered_both_'+fid+'.jpg']
			os.spawnv(os.P_WAIT,'/usr/bin/montage',montageparams)

		else:
			mdfile='md-cfg'+file[len('colour'):-len('.vtk')]+'.asc'
			angle=contact_angle.get_contact_angle(file,mdfile,r)
			t=file.split('_')[-1].split('.')[0][1:]
			cafile.write(t+' '+str(angle)+'\n')
			print t,angle
		#file_params=[sai,'-id='+id,'-nointeractive','-write']
		#os.spawnve(os.P_WAIT,sai,file_params+params,env)
if not contactangle:
	mencoderparams=['/usr/local/bin/mencoder', 'mf://rendered_left_*.png', '-mf', 'fps='+fps, '-o', 'rendered_left.avi', '-ovc', 'lavc', '-lavcopts', 'vcodec=mpeg4:mbd=1:vbitrate=6000']
	print os.spawnv(os.P_WAIT,'/usr/local/bin/mencoder',mencoderparams)
	mencoderparams[1]='mf://rendered_right_*.png'
	mencoderparams[5]='rendered_right.avi'
	print os.spawnv(os.P_WAIT,'/usr/local/bin/mencoder',mencoderparams)
	mencoderparams[1]='mf://rendered_both_*.jpg'
	mencoderparams[5]='rendered_both.avi'
	print os.spawnv(os.P_WAIT,'/usr/local/bin/mencoder',mencoderparams)

