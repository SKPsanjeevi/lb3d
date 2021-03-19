#!/usr/bin/python

# Based on mkani2.py, but should work with omnirender.py
# Contact angle and timestamps stuff are not yet implemented

import os
import sys
import omnirender
import fnmatch

h5tovtk="/usr/local/bin/h5tovtk"

mencoder="/usr/local/bin/mencoder"
#montage="/usr/bin/montage"
#mogrify="/usr/bin/mogrify"

# Set default values
fps = '5'
zstep = '0'
#contactangle = 0
#timestamps = []

if ( len(sys.argv) < 2 ):
	print "Usage: ./omniani.py <filestring> [-fps=n] [-zstep=n]"
	sys.exit()

fncontains = sys.argv[1]

print "OMNIANI initializing parameters"
for arg in sys.argv[2:]:
	if arg.startswith('-fps'):
		fps = arg.split('=')[1]
		print "OMNIANI set FPS =", fps
	elif arg.startswith('-zstep'):
		zstep = arg.split('=')[1]
		print "OMNIANI set zstep =", zstep
	elif arg.startswith('-zmin'):
		zmin = arg.split('=')[1]
		print "OMNIANI set zmin =", zmin
	elif arg.startswith('-zmax'):
		zmax = arg.split('=')[1]
		print "OMNIANI set zmax =", zmax
#	elif arg.startswith('-car'):
#		contactangle=1
#		import contact_angle
#		r = int(arg.split('=')[1])

#	else:
#		timestamps.append(arg)

try:
	paramfile = open("omniani.params","r")
	params = [param for param in paramfile]
	paramfile.close()
except IOError:
	print "OMNIANI no param file found"
	params=[]

params.append("-writeimage")

# files = os.listdir(".")
#if len(timestamps):
#	chars_to_remove=len(timestamps[0])+4
#	os.mkdir(timestamps[0])
#	os.chdir(timestamps[0])
#	for file in files:
#		for timestamp in timestamps:
#			if file.endswith(timestamp+".h5"):
#				newname=file[:-chars_to_remove]+".h5"
#				try:
#					os.symlink("../"+file, newname)
#				except OSError:
#					print "Couldn't create symlink for "+file+", perhaps because we are making a single movie out of several simulations?"
#				else:
#					if os.spawnl(os.P_WAIT,h5tovtk,h5tovtk,'--',newname)!=0:
#						print "error converting "+file+" to vtk"
#					os.remove(newname)
#			if file.endswith(timestamp+".asc"):
#				newname=file[:-(chars_to_remove+1)]+".asc"
#				try:
#					os.symlink("../"+file, newname)
#				except OSError:
#					print "Couldn't create symlink for "+file+", perhaps because we are making a single movie out of several simulations?"

#if contactangle:
#	cafile = open('contact_angle.dat','a')

if (int(zstep) == 0):
	first = 1
	for file in sorted(os.listdir('.')):
		if fnmatch.fnmatch(file,'*'+fncontains+'*'):
			if (file.endswith(".vtk") and ( file.startswith("colour") or file.startswith("rock") ) or file.endswith(".asc") ):
				fid="_".join(file.split('_')[1:])[:-4]
				pngfilename = 'rendered-'+fid+'.png'
				if ( not os.path.isfile(pngfilename) ):
		#			if not contactangle:		
					last = file

					print "OMNIANI calling OmniRender for fid =", fid
					if first:
						first = 0
						sai=omnirender.omnirender(fid,params)
						# sai.update(fid)
					else:
						sai.update(fid)
						
		#			montageparams=[montage,'rendered_left_'+fid+'.png','rendered_right_'+fid+'.png','-tile','2x1','-geometry','+0+0','-title',fid,'rendered_both_'+fid+'.png']
		#			print "OMNIANI calling montage"
		#			os.spawnv(os.P_WAIT,montage,montageparams)
			
		#			else:
		#				mdfile='md-cfg'+file[len('colour'):-len('.vtk')]+'.asc'
		#				angle=contact_angle.get_contact_angle(file,mdfile,r)
		#				t=file.split('_')[-1].split('.')[0][1:]
		#				cafile.write(t+' '+str(angle)+'\n')
		#				print t,angle
		#			#file_params=[sai,'-id='+id,'-nointeractive','-write']
		#			#os.spawnve(os.P_WAIT,sai,file_params+params,env)
				else:
					print "Skipping "+pngfilename+", file already exists."
	
else:
	# Rendering slices along the z-axis, starting with [0:zstep)
	# Need to specify xmax and ymax in the params file, or OmniRender will not work
	params.append('-zmin='+str(max(0,zmin)))
	params.append('-zmax=0')
	first = 1
	for file in sorted(os.listdir('.')):
		if fnmatch.fnmatch(file,'*'+fncontains+'*'):
			if (file.endswith(".vtk") and ( file.startswith("colour") or file.startswith("rock") ) ):
				last = file
				fid="_".join(file.split('_')[1:])[:-4]
				for i in range(int(zstep),int(zmax)+int(zstep),int(zstep)):
					# Replace the previous '-zmax=n' parameter to be passed to OmniRender with a new value
					print "OMNIANI calling OmniRender with slicing for fid =", fid
					findstr = '-zmax='+str(i-int(zstep))
					zmaxindex = params.index(findstr)
					params[zmaxindex]='-zmax='+str(i)
					# Call OmniRender with the new parameter list
					sai=omnirender.omnirender(fid,params)
					# sai.update(fid)
					# Since OmniRender always saves to files 'rendered_left_'+fid+'.png', rename the files now to include the zmax parameter
					#os.rename('rendered_left_'+fid+'.png','rendered_left_'+str(i).zfill(len(zmax))+'_'+fid+'.png')
					#os.rename('rendered_right_'+fid+'.png','rendered_right_'+str(i).zfill(len(zmax))+'_'+fid+'.png')
					os.rename('rendered_'+fid+'.png','rendered_'+str(i).zfill(len(zmax))+'_'+fid+'.png')
					#montageparams=[montage,'rendered_left_'+str(i)+'_'+fid+'.png','rendered_right_'+str(i)+'_'+fid+'.png','-tile','2x1','-geometry','+0+0','-title',fid,'rendered_both_'+str(i).zfill(len(zmax))+'_'+fid+'.png']
					#print "OMNIANI calling montage for slices"
					#os.spawnv(os.P_WAIT,montage,montageparams)

# This part just creates a movie of all the images in the current directory, so if you make multiple movies, move the images of the first somewhere else!
#if not contactangle:
#mencoderparams=[mencoder, 'mf://rendered_left_*.png', '-mf', 'fps='+fps,'-o', 'rendered_left.avi', '-ovc', 'lavc', '-lavcopts','vcodec=mpeg4:mbd=1:vbitrate=6000']
id=fid.split('-')[1]
name=fid.split('_')[0]
mencoderparams=[mencoder, 'mf://rendered-*.png', '-mf', 'fps='+fps,'-o', 'rendered-'+name+'-'+id+'.avi', '-ovc', 'lavc', '-lavcopts','vcodec=mpeg4:mbd=1:vbitrate=6000']
os.spawnv(os.P_WAIT,mencoder,mencoderparams)
#mencoderparams[1]='mf://rendered_right_*.png'
#mencoderparams[5]='rendered_right.avi'
#os.spawnv(os.P_WAIT,mencoder,mencoderparams)
# For some reason, mencoder doesn't like the PNG files from montage, so we convert them to jpg (ew!)
#mogrifyparams=[mogrify, '-format','jpg','rendered_both_*.png']
#os.spawnv(os.P_WAIT,mogrify,mogrifyparams)
#mencoderparams[1]='mf://rendered_both_*.jpg'
#mencoderparams[5]='rendered_both.avi'
#os.spawnv(os.P_WAIT,mencoder,mencoderparams)

