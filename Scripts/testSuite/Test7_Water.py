#!/usr/bin/env python

#*******************************************************************************
# SCRIPT:  Test7_Water.py
# VERSION: 2.0
# FEATURES: Compute the electrostatic energy between two point charges 
#           at a range of distances
#*******************************************************************************

#*******************************************************************************
# IMPORT MODULES
#*******************************************************************************
import subprocess as sp #This module lets us run cassandra from python
import numpy as np #this module s the package for scientific computing in python
import os, sys

#*******************************************************************************
# ARGUMENT PARSE
#*******************************************************************************

#*******************************************************************************
# VARIABLE DEFINITIONS
#*******************************************************************************
# NIST answers
kBoltz = 1.3806488 * 6.02214129 / 1000 # kJ / mol K
nistAnswer = ((9.95387E+04,-8.23715E+02,-5.58889E+05+2.80999E+06,6.27009E+03,-2.84469E+06,-4.88604E+05),
              (1.93712E+05,-3.29486E+03,-1.19295E+06+5.61998E+06,6.03495E+03,-5.68938E+06,-1.06590E+06),
              (3.54344E+05,-7.41343E+03,-1.96297E+06+8.42998E+06,5.24461E+03,-8.53407E+06,-1.71488E+06),
              (4.48593E+05,-1.37286E+04,-3.57226E+06+1.41483E+07,7.58785E+03,-1.42235E+07,-3.20501E+06))
# index = [check][property]

# Simulation parameters
nSpecies = 1
nAtoms = (3,) # index = [species]
atomParms = ((('O','O',15.999,-0.84760,'LJ',78.19743111,3.16555789),
              ('H','H',1.008,0.42380,'LJ',0.,0.),
              ('H','H',1.008,0.42380,'LJ',0.,0.)),) # index = [species][atom][parm]
bondParms = (((1,2,'fixed',1.),
              (1,3,'fixed',1.)),) # index = [species][bond][parm]
angleParms = (((1,2,3,'fixed',109.47),),) # index = [species][angle][parm]
prpList = ("energy_intervdw", "energy_lrc", "energy_interq", "energy_reciprocal", 
           "energy_self", "energy_total") # index = [property]
vdwStyle = 'lj cut_tail 10.'

# Check parameters
# params for each simulation
numChecks = 4 # number of simulations to run
nMols = (100, 200, 300, 750) # index = [check]
box = (20., 20., 20., 30.)   # index = [check]
chargeStyle = ('coul ewald 10. 0.000393669','coul ewald 10. 0.000393669',
               'coul ewald 10. 0.000393669','coul ewald 10. 0.0306708') # index = [check]
cassAnswer     = [None] * numChecks #this list will hold cassandra's answers
passCheck      = [None] * numChecks
errorTol = 1e-5

# Formatting variables
bold = '\033[1m' #Will make text bold
normal = '\033[0m' #Will make the next text normal(ie. unbold)

#*******************************************************************************
# FILE MANAGEMENT
#*******************************************************************************
cassDir = "/Users/rmullen2/dev/cassandra/Src/"
cassExe = "cassandra.dev"
cassRun = "test7.out"
inpName = "test7.inp"
xyzName = ("inputFiles/test7.water1.xyz", "inputFiles/test7.water2.xyz",
           "inputFiles/test7.water3.xyz", "inputFiles/test7.water4.xyz") # index = [check]
mcfName = ("test7.water.mcf",) # index = [species]


#*******************************************************************************
# FUNCTION DEFINITIONS
#*******************************************************************************
# This functon writes over specific lines of the input file created above using a function, called replace_line that is created below. 
# This function takes three inputs: the name of the file where you would like to replace a line, the line number you would like to replace, and the text you want to replace the old text with. 
def replace_line(file_name, line_num, text):
	lines = open(file_name, 'r').readlines()
	lines[line_num] = text
	out = open(file_name, 'w')
	out.writelines(lines)
	out.close() # Closes the file so that the program doesn't explode. 

# Error function
def erf(x):
	# constants
	a1 =  0.254829592
	a2 = -0.284496736
	a3 =  1.421413741
	a4 = -1.453152027
	a5 =  1.061405429
	p  =  0.3275911

	# Save the sign of x
	sign = 1
	if x < 0:
		sign = -1
	x = abs(x)

	# A & S 7.1.26
	t = 1.0/(1.0 + p*x)
	y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*np.exp(-x*x)

	return sign*y
#*******************************************************************************
# MAIN PROGRAM BEGINS HERE
#*******************************************************************************
#
# Step 1) Write input files
#
# Loop through the checks
# 
# Step 2) Run Cassandra to get its answer
# Step 3) Compare answers
#

#This prints the starting line.
print "\n\n"+bold+"Test 7: Energy of water configs" + normal 

# Step 1) Write input files
# 1.1) Write MCF for cation, anion
for s in range(nSpecies):
	mcf = open(mcfName[s],"w") #Creates mcf file and allows for edits
	nAtoms = len(atomParms[s])
	mcf.write("# Atom_Info\n%d\n" % (nAtoms))
	for i in range(nAtoms):
		mcf.write("%d %s %s %.1f %.4f %s %.8f %.8f\n" % ((i+1,)+atomParms[s][i]))
	nBonds = len(bondParms[s])
	mcf.write("\n# Bond_Info\n%d\n" % (nBonds))
	for i in range(nBonds):
		mcf.write("%d %d %d %s %.8f\n" % ((i+1,)+bondParms[s][i]))
	nAngles = len(angleParms[s])
	mcf.write("\n# Angle_Info\n%d\n" % (nAngles))
	for i in range(nAngles):
		mcf.write("%d %d %d %d %s %.8f\n" % ((i+1,)+angleParms[s][i]))
	mcf.write("\n# Dihedral_Info\n0\n")
	mcf.write("\n# Improper_Info\n0\n")
	mcf.write("\n# Fragment_Info\n0\n")
	mcf.write("\n# Fragment_Connectivity\n0\n")
	mcf.write("\n# Intra_Scaling\n0. 0. 0. 0.\n0. 0. 0. 0.\n")
	mcf.write("\nEND\n")
	mcf.close()

# Loop through checks
print "%-25s %18s %18s %18s %8s" % ("Property","Cassandra","NIST","Relative_Err","Pass")
for i in range(numChecks):

	# Step 2) Run Cassandra to get its answer
	# 2.1) Write inp file

	# Combine file names and nmols
	mcfStr = mcfName[0] + " " + str(nMols[i]) + '\n'

	# Write inp file
	inp = open(inpName,"w")
	inp.write("# Run_Name\n%s\n\n" % (cassRun))
	inp.write("# Sim_Type\nnvt\n\n")
	inp.write("# Nbr_Species\n%d\n\n" % (nSpecies))
	inp.write("# VDW_Style\n%s\n\n" % (vdwStyle))
	inp.write("# Charge_Style\n%s\n\n" % (chargeStyle[i]))
	inp.write("# Seed_Info\n1 2\n\n")
	inp.write("# Rcutoff_Low\n0.1\n\n")
	inp.write("# Molecule_Files\n%s\n" % (mcfStr))
	inp.write("# Box_Info\n1\ncubic\n%.1f\n\n" % (box[i]))
	inp.write("# Temperature_Info\n1.0\n\n")
	inp.write("# Move_Probability_Info\n\n")
	inp.write("# Prob_Translation\n1.0\n0.0 0.0\n\n")
	inp.write("# Done_Probability_Info\n\n")
	inp.write("# Start_Type\nread_config %d %s\n\n" % (nMols[i],xyzName[i]))
	inp.write("# Run_Type\nEquilibration 100\n\n")
	inp.write("# Simulation_Length_Info\nunits steps\nprop_freq 1\ncoord_freq 1\nrun 0\n\n")
	inp.write("# Property_Info 1\n")
	for prp in prpList:
		inp.write("%s\n" % (prp))
	inp.write("\n")
	inp.write("END\n")
	inp.close()

	# 3.3) Run Cassandra Jobs
	proc = sp.Popen([cassDir + cassExe + " " + inpName], stdout=sp.PIPE, shell=True)
	(out, err) = proc.communicate()

	if err is not None:
		print("Error.Abort.")

	# 3.3) Read logfile
	prpFile = open(cassRun + ".prp", "r")
	# search line by line in log for the words "Total system energy"
	nPrp = len(prpList)
	prpFile.readline() # header1
	prpFile.readline() # header2
	prpFile.readline() # header3
	cassAnswer[i] = [float(prp) for prp in prpFile.readline().split()[1:]]
	prpFile.close()

	# Step 4) Compare answers
	for j in range(nPrp):
		nist = nistAnswer[i][j]*kBoltz
		if (nist == 0.):
			passCheck = cassAnswer[i][j] == 0.
			print "%-25s %18.6g %18.6g %18s %8s" % (cassStr[j],
						cassAnswer[i][j],nist,'',passCheck)
		else:
			errorRel = abs(cassAnswer[i][j] - nist)/nist
			passCheck = abs(errorRel) <= errorTol
			print "%-25s %18.6g %18.6g %18.6g %8s" % (cassStr[j],
						cassAnswer[i][j],nist,errorRel,passCheck)
	print ""

# Clean up scratch files
os.system('rm ' + inpName)
os.system('rm ' + ' '.join(mcfName))
os.system('rm ' + cassRun + '*')
