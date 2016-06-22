import os, subprocess, sys
bold = '\033[1m'
normal = '\033[0m'

#############Settings##############
examples_dir = '/afs/crc.nd.edu/user/e/emarinri/Git/Cassandra/Examples'
cassandra = '/afs/crc.nd.edu/user/e/emarinri/Git/Cassandra/Src/cassandra_gfortran.exe'
fraggen = '/afs/crc.nd.edu/user/e/emarinri/Git/Cassandra/Scripts/Frag_Library_Setup/library_setup.py'
libgenfile = open('/afs/crc.nd.edu/user/e/emarinri/Git/Cassandra/Examples/libgenout.dat','a')
cassoutfile = open('/afs/crc.nd.edu/user/e/emarinri/Git/Cassandra/Examples/sim.out','a')
###################################

#for root, dirs, files in os.walk(examples_dir):
#	for thisfile in files:
#		if thisfile.endswith('.out.chk') or \
#                   thisfile.endswith('.log') or \
#                   thisfile.endswith('.BAK') or \
#                   '.box' in thisfile or \
#                   '.out' in thisfile:
#			subprocess.call(['rm',root+'/'+thisfile])
#
#	try:
#		
#		if 'species' in root:
#			subprocess.call(['rm','-r',root])
#	except:
#		pass
#print bold + 'Done clean up'

pdbfiles = []
inputfiles = []
execute = 0
for root, dirs, files in os.walk(examples_dir):
	for thisfile in files:
		if thisfile.endswith('pdb') and \
                   'GCMC' not in root and \
                   'molecule.pdb' not in thisfile:
			pdbfiles.append(thisfile)
			execute = execute + 1
		if thisfile.endswith('inp') and 'species' not in thisfile and \
                   'frag' not in thisfile and \
                   'GCMC' not in root:
			inputfiles.append(thisfile)
			execute = execute + 1
		if execute == 2:
			execute = 0
			os.chdir(root)
			print bold + inputfiles[-1] + ' ' + pdbfiles[-1] + normal
			libgenfile.write('\n\n' + inputfiles[-1] + \
                           ' ' + pdbfiles[-1] + '\n\n')
			thisprocess = subprocess.Popen \
                           ([fraggen, cassandra, inputfiles[-1], pdbfiles[-1]],\
                            stdout = subprocess.PIPE)
			stdoutvalue = thisprocess.communicate()[0]
			libgenfile.write(stdoutvalue)
			os.chdir(examples_dir)
libgenfile.close()
for root, dirs, files in os.walk(examples_dir):
	for thisfile in files:
		if thisfile.endswith('inp') and \
                'GCMC' not in root and \
                'species' not in thisfile and \
                'frag' not in thisfile:
			os.chdir(root)
			print  cassandra, root+'/'+thisfile
			cassoutfile.write('\n\n' + cassandra + ' ' + \
                           root+'/'+thisfile + '\n\n')
			thisprocess = subprocess.Popen([cassandra, root+\
                           '/'+thisfile], stdout = subprocess.PIPE)
			stdoutvalue = thisprocess.communicate()[0]
			cassoutfile.write(stdoutvalue)
			os.chdir(examples_dir)
cassoutfile.close()
