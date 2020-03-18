#!/usr/bin/env python

#*******************************************************************************
# SCRIPT:  Test7_Water.py
# VERSION: 2.0
# FEATURES: Water box, electrostatic energy
#*******************************************************************************

#*******************************************************************************
# IMPORT MODULES
#*******************************************************************************
from __future__ import print_function
import subprocess as sp #This module lets us run cassandra from python
import numpy as np #this module s the package for scientific computing in python
import os, sys
import argparse

#*******************************************************************************
# ARGUMENT PARSE
#*******************************************************************************

parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=
"""DESCRIPTION:
Runs the test suite using the Cassandra executable specified, including path.

EXAMPLES:
To run the test suite using a Cassandra executable inside of Cassandra/Src/:

	> python testSuite.py ../../Src/cassandra.exe
	> python testSuite.py ../../Src/cassandra_gfortran.exe

To run a test using a Cassandra executable elsewhere:

	> python Test#_Description.py /home/applications/cassandra.exe 
	> python Test#_Description.py /home/applications/cassandra_gfortran.exe

""")
parser.add_argument('cassandra_exe', 
                help="Cassandra executable, including path. To call an executable in the same"+
                "Cassandra package's Src folder, utilize ../../Src/cassandra.exe as the executable path.")

args = parser.parse_args()

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

test_no = 7
test_desc = "Water Energy"

# Simulation parameters
nSpecies = 1
nAtoms = (3,) # index = [species]
atomParms = ((('O','O',15.999,-0.84760,'LJ',78.19743111,3.16555789),
              ('H','H',1.008,0.42380,'LJ',0.,0.),
              ('H','H',1.008,0.42380,'LJ',0.,0.)),) # index = [species][atom][parm]
bondParms = (((1,2,'fixed',1.),
              (1,3,'fixed',1.)),) # index = [species][bond][parm]
angleParms = (((2,1,3,'fixed',109.47),),) # index = [species][angle][parm]
cassStr = (("Inter molecule vdw","Long range correction","Inter molecule q","Reciprocal ewald","Self ewald","Total system energy"),) * 4 # one tuple for each check

cassPrint = (("Inter molecule vdw [kJ/mol-Ext]","Long range correction [kJ/mol-Ext]","Inter molecule q [kJ/mol-Ext]","Reciprocal ewald [kJ/mol-Ext]","Self ewald [kJ/mol-Ext]","Total system energy [kJ/mol-Ext]"),) * 4 # one tuple for each check

prpList = ("energy_intervdw", "energy_lrc", "energy_interq", "energy_recip", 
           "energy_self", "energy_total") # index = [property]
vdwStyle = 'lj cut_tail 10.'

# Check parameters
# params for each simulation
numChecks = 4 # number of simulations to run
nMols = (100, 200, 300, 750) # index = [check]
box = (20., 20., 20., 30.)   # index = [check]
title = ('100 molecules','200 molecules', '300 molecules','750 molecules',)
chargeStyle = ('coul ewald 10. 0.000393669','coul ewald 10. 0.000393669',
               'coul ewald 10. 0.000393669','coul ewald 10. 0.0306708') # index = [check]
#chargeStyle = ('coul ewald 10. 1e-5','coul ewald 10. 1e-5',
 #              'coul ewald 10. 1e-5','coul ewald 10. 3.05e-2') # index = [check]
cassAnswer     = [None] * numChecks #this list will hold cassandra's answers
passCheck      = [None] * numChecks
errorTol = (1e-5,5e-5,1e-5,8e-2, 1e-5, 5e-4)

# Formatting variables
bold = '\033[1m' #Will make text bold
normal = '\033[0m' #Will make the next text normal(ie. unbold)

#*******************************************************************************
# FILE MANAGEMENT
#*******************************************************************************
testSuiteFolder = os.getcwd()
MainDir 	= testSuiteFolder[0:len(testSuiteFolder)-len('Scripts/testSuite')]
resourceDir = MainDir + "Scripts/testSuite/Resources"
cassExe     = args.cassandra_exe
cassRun 	= "test7.out"
inpName 	= "test7.inp"
xyzName = (resourceDir+"/inputFiles/test7.water1.xyz", resourceDir+"/inputFiles/test7.water2.xyz",
           resourceDir+"/inputFiles/test7.water3.xyz", resourceDir+"/inputFiles/test7.water4.xyz") # index = [check]
xyzSimpleName = ("test7.water1.xyz", "test7.water2.xyz",
           "test7.water3.xyz", "test7.water4.xyz") # index = [check]
mcfName 	= ("test7.water.mcf",) # one string per species

#*******************************************************************************
# FUNCTION DEFINITIONS
#*******************************************************************************

from testSuiteFunctions import replace_line

from testSuiteFunctions import erf
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
print("\n\n"+bold+"Test 7: Energy of Water Configurations" + normal) 


FailCount = 0; 

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
	mcf.write("\n# Intra_Scaling\n0. 0. 0. 1.0\n0. 0. 0. 1.0\n")
	mcf.write("\nEND\n")
	mcf.close()

# Loop through checks
print("%-30s %-35s %18s %18s %18s %8s" % ("Title", "Property","Cassandra","NIST","Relative_Err","Pass"))
for i in range(numChecks):

	# Step 2) Run Cassandra to get its answer
	# Pull in resource xyz file
	os.system('cp '+xyzName[i]+' ' +xyzSimpleName[i])

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
	inp.write("# Rcutoff_Low\n0.850\n\n")
	inp.write("# Molecule_Files\n%s\n" % (mcfStr))
	inp.write("# Box_Info\n1\ncubic\n%.1f\n\n" % (box[i]))
	inp.write("# Temperature_Info\n300.0\n\n")
	inp.write("# Move_Probability_Info\n\n")
	inp.write("# Prob_Translation\n1.0\n0.0 0.0\n\n")
	inp.write("# Done_Probability_Info\n\n")
	inp.write("# Start_Type\nread_config %d %s\n\n" % (nMols[i],xyzSimpleName[i]))
	inp.write("# Run_Type\nEquilibration 100\n\n")
	inp.write("# Simulation_Length_Info\nunits steps\nprop_freq 1\ncoord_freq 1\nrun 0\n\n")
	inp.write("END\n")
	inp.close()

	# 3.3) Run Cassandra Jobs
	proc = sp.Popen([cassExe + " " + inpName], stdout=sp.PIPE, shell=True)
	(out, err) = proc.communicate()

	if err is not None:
		print("Error.Abort.")

	os.system('rm '+xyzSimpleName[i])

	# 3.3) Read logfile
	log = open(cassRun + ".log", "r")
	# search line by line in log for the words "Total system energy"
	nPrp = len(cassStr[i])
	cassAnswer[i] = [None] * nPrp
	for line in log:
		for j in range(nPrp):
			if (cassStr[i][j] in line):
				cassAnswer[i][j] = float(line.split()[-1])

	# Step 4) Compare answers
	for j in range(nPrp):
		nist = nistAnswer[i][j]*kBoltz
		if (nist == 0.):
			passCheck = cassAnswer[i][j] == 0.
			if passCheck == 0:
				FailCount = FailCount+1;
				failureOutString = MainDir+ 'Scripts/testSuite/failureLog/test' + str(test_no) + '_check' + str(i+1)
				os.system('mkdir -p ' + failureOutString)
				os.system('cp ' + inpName + ' ' + failureOutString )
				os.system('cp ' + xyzName[i] + ' ' + failureOutString )
				os.system('cp ' + ' '.join(mcfName) + ' ' + failureOutString )
				os.system('cp ' + cassRun + '*' + ' ' + failureOutString )
			if (j == 0):
				print("%-30s %-35s %18.6g %18.6g %18s %8s" % (title[i],cassPrint[i][j],
						cassAnswer[i][j],nist,'',passCheck))
			else: 
				print("%-30s %-36s %17.6g %18.6g %18s %8s" % ('',cassPrint[i][j],
						cassAnswer[i][j],nist,'',passCheck))

		else:
			errorRel = abs(cassAnswer[i][j] - nist)/nist
			passCheck = abs(errorRel) <= errorTol[j]
			if passCheck == 0:
				FailCount = FailCount+1;
				failureOutString = MainDir+ 'Scripts/testSuite/failureLog/test' + str(test_no) + '_check' + str(i+1)
				os.system('mkdir -p ' + failureOutString)
				os.system('cp ' + inpName + ' ' + failureOutString )
				os.system('cp ' + xyzName[i] + ' ' + failureOutString )
				os.system('cp ' + ' '.join(mcfName) + ' ' + failureOutString )
				os.system('cp ' + cassRun + '*' + ' ' + failureOutString )
			if (j == 0):
				print("%-30s %-35s %18.6g %18.6g %18.6g %8s" % (title[i],cassPrint[i][j],
						cassAnswer[i][j],nist,errorRel,passCheck))
			else: 
				print("%-30s %-36s %17.6g %18.6g %18.6g %8s" % ('',cassPrint[i][j],
						cassAnswer[i][j],nist,errorRel,passCheck))

if (FailCount != 0):
	PassState = "False"
else:
	PassState = "True"

LastTest = open(MainDir+ 'Scripts/testSuite/testOutput/LastTest.txt',"w")
LastTest.write(PassState)
LastTest.close()


# Clean up scratch files
os.system('rm ' + inpName)
os.system('rm ' + ' '.join(mcfName))
os.system('rm ' + cassRun + '*')
