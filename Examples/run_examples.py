import os, subprocess, sys
import argparse
bold = '\033[1m'
normal = '\033[0m'

#*******************************************************************************
# ARGUMENT PARSE
#*******************************************************************************

parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=
"""DESCRIPTION:
Runs all examples using the Cassandra executable specified, including path. 
Also cleans all previous results, and the fragment libraries of all non-GCMC 
examples.

Note: It is neccessary to speficy an *EXACT* path. The use of ../ to specify 
a relative path will prohibit the script from running properly. 

EXAMPLES:
To run the examples using a Cassandra executable elsewhere:

	> python run_examples.py /home/applications/cassandra.exe 
	> python run_examples.py /home/applications/cassandra_gfortran.exe

""")
parser.add_argument('cassandra_exe', 
                help="Cassandra executable, including path. To call an executable in the same"+
                "Cassandra package's Src folder, utilize ../Src/cassandra.exe as the executable path.")

args = parser.parse_args()

#############Settings##############
examples_dir = os.getcwd()
cassandra = args.cassandra_exe
fraggen = examples_dir[0:-9]+'/Scripts/Frag_Library_Setup/library_setup.py'
libgenfile = open('libgenout.dat','a')
cassoutfile = open('sim.out','a')
###################################

for root, dirs, files in os.walk(examples_dir):
	for thisfile in files:
		if (thisfile.endswith('.out.chk') or \
                   thisfile.endswith('.log') or \
                   thisfile.endswith('.BAK') or \
                   '.box' in thisfile or \
                   '.out' in thisfile):
			subprocess.call(['rm',root+'/'+thisfile])

	try:
		if 'species' in root and \
                   'GCMC' not in root:
			subprocess.call(['rm','-r',root])
	except:
		pass
print bold + 'Done clean up'

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
                           ([fraggen, cassandra, inputfiles[-1], pdbfiles[-1],'-n 100'],\
                            stdout = subprocess.PIPE)
			stdoutvalue = thisprocess.communicate()[0]
			libgenfile.write(stdoutvalue)
			os.chdir(examples_dir)
libgenfile.close()
for root, dirs, files in os.walk(examples_dir):
	for thisfile in files:
		if thisfile.endswith('inp') and \
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
