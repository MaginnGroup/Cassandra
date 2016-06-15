import os, subprocess, sys
bold = '\033[1m'
normal = '\033[0m'


#######SETTINGS#######

#Where are the examples located?
examples_dir = '/afs/crc.nd.edu/user/e/emarinri/Git/Cassandra_V1.1/Examples/'
#Which is Cassandra executable will I use?
cassandra =    '/afs/crc.nd.edu/user/e/emarinri/Git/Cassandra_V1.1/Examples/cassandra_intel_openMP.exe'
#Where is library_setup.py script?
fraggen =      '/afs/crc.nd.edu/user/e/emarinri/Git/Cassandra_V1.1/Examples/library_setup.py'
#Output file for frag generation outputs
libgenfile = open(examples_dir + 'libgen.out','a')
#Output file for cassandra simulation outputs
cassoutfile = open(examples_dir + 'sim.out','a')

#######SETTINGS#######


#Remove old files
for root, dirs, files in os.walk(examples_dir):
	for thisfile in files:
		if thisfile.endswith('.chk') or \
                   thisfile.endswith('.log') or \
                   thisfile.endswith('.BAK') or \
                   '.box' in thisfile or \
                   '.out' in thisfile:
			subprocess.call(['rm',root+'/'+thisfile])

	try:
		
		if 'species' in root:
			subprocess.call(['rm','-r',root])
	except:
		pass

print bold + 'Done clean up'

pdbfiles = []
inputfiles = []
execute = 0
for root, dirs, files in os.walk(examples_dir):
	for thisfile in files:
                #Append each PDB files for each species
		if thisfile.endswith('pdb') and \
                   'GCMC' not in root and \
                   'molecule.pdb' not in thisfile:
			pdbfiles.append(thisfile)
			execute = execute + 1

                #Append input files for simulation. Exclude input files 'species' not in thisfile and for fragment/mcf generation.
		if thisfile.endswith('inp') and \
                   'frag' not in thisfile and \
                   'GCMC' not in root:
			inputfiles.append(thisfile)
			execute = execute + 1

                #We should have all necessary information to generate frag libs
		if execute == 2: 
			execute = 0
			os.chdir(root)
			print bold + inputfiles[-1] + ' ' + pdbfiles[-1] + normal
			libgenfile.write('\n\n' + inputfiles[-1] + ' ' + pdbfiles[-1] + '\n\n')
			thisprocess = subprocess.Popen([fraggen, cassandra, inputfiles[-1], pdbfiles[-1]], stdout = subprocess.PIPE)
			stdoutvalue = thisprocess.communicate()[0]
			libgenfile.write(stdoutvalue)
			os.chdir(examples_dir)
libgenfile.close()
#for root, dirs, files in os.walk(examples_dir):
#	for thisfile in files:
#		if thisfile.endswith('inp') and \
#                'GCMC' not in root and \
#                'species' not in thisfile and \
#                'frag' not in thisfile:
#			os.chdir(root)
#			print  cassandra, root+'/'+thisfile
#			cassoutfile.write('\n\n' + cassandra + ' ' + root+'/'+thisfile + '\n\n')
#			thisprocess = subprocess.Popen([cassandra, root+'/'+thisfile], stdout = subprocess.PIPE)
#			stdoutvalue = thisprocess.communicate()[0]
#			cassoutfile.write(stdoutvalue)
#			os.chdir(examples_dir)
#cassoutfile.close()
