#!/usr/bin/python
#This skript tries to compile LBE with different flags. The flags are read form compile_flags.dat
import os
import sys
try:
	import subprocess
#JUMP uses a very ancient version of Python where the subprocess module is not available, the following is a very ugly hack around this
except ImportError:
	class subprocess:
	   class PIPE:
		def __init__(self):
			pass
	   class Popen:
		def __init__(self,params,stdout,stderr):
			self.params=params
		def wait(self):
			return os.spawnv(os.P_WAIT, self.params[0], self.params)
		def communicate(self):
			(stdin, stdout, stderr)=os.popen3(' '.join(self.params))
			return [stdout.read(),stderr.read()]
if len(sys.argv)!=2:
	print 'Usage:'
	print sys.argv[0]+' PLATFORM'
	print 'where PLATFORM is the same as in LBECONFIG.sh'
	os._exit(1)

platform=sys.argv[1]
cwd=os.getcwd()
if not 'LBECONFIG.sh' in os.listdir(cwd):
	print 'Please call this skript from the directory containing LBECONFIG.sh'
	os._exit(1)
lbeconfig=os.path.join(cwd,'LBECONFIG.sh')

supported_platform=False

#Now come the parameters needed so LBECONFIG.sh is called in the right way on different platforms

if platform.startswith('LINUX64'):
	flags_start='' # prepend this to to the list of flags
	flags_prepend='-D' #prepend this to every single flag
	flags_join=False # join flags together to a single argument
	execute_before_clean='' # execute this shell command before running LBECONFIG.sh CLEAN
	execute_before_config='' # execute this shell command before running LBECONFIG.sh CONFIG
	execute_before_make='' # execute this shell command before running LBECONFIG.sh MAKE
	supported_platform=True

if platform.startswith('SP2MPI'):
	flags_start='-WF' # prepend this to to the list of flags
	flags_prepend=',-D' #prepend this to every single flag
	flags_join=True # join flags together to a single argument
	execute_before_clean='' # execute this shell command before running LBECONFIG.sh CLEAN
	execute_before_config='cp code/xdrf/Makefile code/xdrf/Makefile_working' # execute this shell command before running LBECONFIG.sh CONFIG
	execute_before_make='mv code/xdrf/Makefile_working code/xdrf/Makefile' # execute this shell command before running LBECONFIG.sh MAKE
	supported_platform=True


if not supported_platform:
	raise exception('This script was not adapted to the current platform. Please add the right flags values for flags_start, flags_prepend, flags_join and execute_before_make.')

#Open a list of whitespace-separated flags to test compilation for
flagsfile=open(os.path.join(os.path.dirname(sys.argv[0]),'compile_flags.dat'),'r')
	
	
for flagsline in flagsfile:
	if not flagsline.startswith('#'):
		flags=[flags_prepend+i for i in flagsline.split()]
		# Clean up
		if len(execute_before_clean)>0:
			os.system(execute_before_clean)
		if not subprocess.Popen([lbeconfig,'CLEAN'],stdout=subprocess.PIPE,stderr=subprocess.PIPE).wait()==0:
			raise Exception('CLEAN failed')
		# Put flags together in the way this platforms expects them
		if len(execute_before_config)>0:
			os.system(execute_before_config)
		lbeconfig_params=[lbeconfig,'CONFIG',platform]
		if (len(flags_start)>0)&(not flags_join):
			lbeconfig_params+=flags_start
		if flags_join:
			flags=flags_start+''.join(flags)
			lbeconfig_params.append(flags)
		else:
			lbeconfig_params+=flags
		# ... and call LBECONFIG.sh CONFIG
		if not subprocess.Popen(lbeconfig_params,stdout=subprocess.PIPE,stderr=subprocess.PIPE).wait()==0:
			raise Exception('CONFIG failed for'+flagsline)
		if len(execute_before_make)>0:
			os.system(execute_before_make)
		print '\nCalling LBECONFIG.sh MAKE for the following flags:'
		print flagsline
		# ... MAKE
		lbeconfigoutput=['','']
		for dir in ['code/xdrf','code']:
			os.chdir(os.path.join(cwd,dir))
			makeoutput=subprocess.Popen(['make'],stdout=subprocess.PIPE,stderr=subprocess.PIPE).communicate()
			os.chdir(cwd)
			lbeconfigoutput=[lbeconfigoutput[i]+makeoutput[i] for i in range(len(makeoutput))]
		print lbeconfigoutput[1]
		# Try to find out if compilation worked, sadly the return value is meaningless.
		#This is certainly not the best way to do it.
		if 'rror' in lbeconfigoutput[1]:
			print 'Compilation failed! Flags were:'
			print flagsline
			os._exit(1)
		#Warnings may be called different in other locales, but the detection is not that important
		if 'arning' in lbeconfigoutput[1]:
			print 'Compilation successfull, warnings above!'
		else:
			print 'Compilation successfull.'

print 'All tests compiled successfully'
print 'Remember du run LBECONFIG.sh as usual to get your normal working lbe executable'
