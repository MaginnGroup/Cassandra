#!/usr/bin/env python 
#*******************************************************************************
# SCRIPT:  Test6_NIST_LennardJonesEnergy.py
# VERSION: 2.0
# FEATURES: Compute the LJ energy of a sample of molecules and compare to NIST:
#            (1) with a cut_tail of 3*sigma
#			 (2) with a cut_tail of 4*sigma
#*******************************************************************************

#*******************************************************************************
# IMPORT MODULES
#*******************************************************************************
from __future__ import print_function
import subprocess as sp #This module lets us run cassandra from python
import numpy as np #this module's the package for scientific computing in python
import random #This allows us to run random numbers
import os, sys
import re
from itertools import *
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
# Test description
test_no = 6
test_desc = "LJ NIST Energy"

# Simulation parameters
# Parms that will be the same for every check run with this script
nSpecies = 1
nAtoms = (1,) # one integer per species
atomParms = [None] * nSpecies # one list per species
s = 0
atomParms[s] = [None] * nAtoms[s] # one dictionary per atom in species
a = 0
atomParms[s][a] = { 'name':'LJ', 'type':'LJ', 'element':'LJ', 'mass':'1.0' }
atomParms[s][a]['charge'] = '0.' # store as a string
atomParms[s][a]['vdw'] = ('LJ', '120.272', '3.540')
bondList = [None]*nSpecies; bondParms = {} # no bonds
angleList = [None]*nSpecies; angleParms = {} # no angles
dihedList = [None]*nSpecies; dihedParms = {} # no diheds
improperList = [None]*nSpecies; improperParms = {} # no impropers
box = (35.4,35.4,28.32,28.32,35.4,35.4,28.32,28.32,)
vdwStyle = ('lj cut_tail 10.620','lj cut_tail 10.620','lj cut_tail 10.620','lj cut_tail 10.620',
		'lj cut_tail 14.16','lj cut_tail 14.16','lj cut_tail 14.16','lj cut_tail 14.16',)
chargeStyle = 'none'

# Check parameters
# Parms that will change
nChecks = 8 # number of simulations to run
title = ("Cut 3*sigma, 800 molecules", "Cut 3*sigma, 400 molecules", "Cut 3*sigma, 200 molecules","Cut 3*sigma, 30 molecules",
		"Cut 4*sigma, 800 molecules", "Cut 4*sigma, 400 molecules", "Cut 4*sigma, 200 molecules","Cut 4*sigma, 30 molecules")
analyticAnswer = [None] * nChecks # list to hold analytic answers
cassAnswer     = [None] * nChecks # list to hold cassandra's answers
passCheck      = [None] * nChecks # list to hold if cassandra passed each check
cassStr = (("Total system energy","Inter molecule vdw","Long range correction"),) * nChecks # one tuple for each check
cassPrint = (("Total system energy [kJ/mol-Ext]","Inter molecule vdw [kJ/mol-Ext]","Long range correction [kJ/mol-Ext]"),) * nChecks # one tuple for each check
errorTol = .01

# Number of moelcules changing
nMolsByCheck= [None] * nChecks
nMolsByCheck[0] = (800,)
nMolsByCheck[1] = (400,)
nMolsByCheck[2] = (200,)
nMolsByCheck[3] = (30,)
nMolsByCheck[4] =nMolsByCheck[0]
nMolsByCheck[5] =nMolsByCheck[1]
nMolsByCheck[6] =nMolsByCheck[2]
nMolsByCheck[7] =nMolsByCheck[3]

# Formatting variables
bold = '\033[1m' #Will make text bold
normal = '\033[0m' #Will make the next text normal(ie. unbold)

#*******************************************************************************
# FILE MANAGEMENT
#*******************************************************************************
testSuiteFolder = os.getcwd()
MainDir 	= testSuiteFolder[0:len(testSuiteFolder)-len('Scripts/testSuite')]
resourceDir = MainDir + "Scripts/testSuite/Resources/"
cassExe     = args.cassandra_exe
cassRun 	= "test6.out"
inpName 	= "test6.inp"
mcfName 	= ("test6.mcf",) # one string per species

# Coordinate Files
xyzfileByCheck= [None] * nChecks
xyzfileByCheck[0] = resourceDir + "testConfigurations/nist1.xyz"
xyzfileByCheck[1] = resourceDir + "testConfigurations/nist3.xyz"
xyzfileByCheck[2] = resourceDir + "testConfigurations/nist2.xyz" 
xyzfileByCheck[3] = resourceDir + "testConfigurations/nist4.xyz"
xyzfileByCheck[4] = xyzfileByCheck[0]
xyzfileByCheck[5] = xyzfileByCheck[1]
xyzfileByCheck[6] = xyzfileByCheck[2]
xyzfileByCheck[7] = xyzfileByCheck[3]
xyzSimplefileByCheck= [None] * nChecks
xyzSimplefileByCheck[0] = "nist1.xyz"
xyzSimplefileByCheck[1] = "nist3.xyz"
xyzSimplefileByCheck[2] = "nist2.xyz" 
xyzSimplefileByCheck[3] = "nist4.xyz"
xyzSimplefileByCheck[4] = xyzSimplefileByCheck[0]
xyzSimplefileByCheck[5] = xyzSimplefileByCheck[1]
xyzSimplefileByCheck[6] = xyzSimplefileByCheck[2]
xyzSimplefileByCheck[7] = xyzSimplefileByCheck[3]

#*******************************************************************************
# FUNCTION DEFINITIONS
#*******************************************************************************

#*******************************************************************************
# MAIN PROGRAM BEGINS HERE
#*******************************************************************************
# Step 1) Calculate the correct answer
# Step 2) Write input files
# Step 3) Run Cassandra to get its answer
# Step 4) Compare answers
#

#This prints the starting line.
print("\n\n"+bold+"Test " + str(test_no) +": " + test_desc + normal)

FailCount = 0; 

analyticAnswer[0] = (-4548.72 , -4351.5 , -198.49 , )
analyticAnswer[1] = (-1194.59, -1144.96, -49.622)
analyticAnswer[2] = (-710.96, -690.0, -24.230, )
analyticAnswer[3] = (-17.46, -16.92, -.54517, )
analyticAnswer[4] = (-4549.95 , -4467.5 , -83.769 , )
analyticAnswer[5] = (-1194.62, -1175.4, -20.942)
analyticAnswer[6] = (-711.56, -704.60 , -10.210, )
analyticAnswer[7] = (-17.41 , -17.060 , -.23008 , )

# Loop through checks
print("%-30s %-35s %18s %18s %18s %8s" % ("Title", "Property","Cassandra","NIST","Relative_Err","Pass"))
for i in range(nChecks):

	# Step 2) Write input files
	# Pull in resource xyz file
	os.system('cp '+xyzfileByCheck[i]+' ' +xyzSimplefileByCheck[i])
	
	# 2.1) Write MCF
	for s in range(nSpecies):
		mcf = open(mcfName[s],"w") #Creates mcf file and allows for edits
		mcf.write("# Atom_Info\n%d\n" % (nAtoms[s]))
		for j,a in enumerate(atomParms[s]):
			mcf.write("%d %s %s %s %s" % (j+1,a['name'],a['element'],a['mass'],a['charge']))
			for parm in a['vdw']:
				mcf.write(" %s" % (parm))
			mcf.write("\n")
		mcf.write("\n# Bond_Info\n0\n")
		mcf.write("\n# Angle_Info\n0\n")
		mcf.write("\n# Dihedral_Info\n0\n")
		mcf.write("\n# Improper_Info\n0\n")
		mcf.write("\n# Fragment_Info\n0\n")
		mcf.write("\n# Fragment_Connectivity\n0\n")
		mcf.write("\n# Intra_Scaling\n0. 0. 0. 0.\n0. 0. 0. 0.\n")
		mcf.write("\nEND\n")
		mcf.close()

	# 2.2) Write inp file
	# Combine file names and nmols
	mcfStr = ''
	for s in range(nSpecies):
		mcfStr = mcfStr + mcfName[s] + " " + str(nMolsByCheck[i][s]) + '\n'

	# Write inp file
	inp = open(inpName,"w")
	inp.write("# Run_Name\n%s\n\n" % (cassRun))
	inp.write("# Sim_Type\nnvt\n\n")
	inp.write("# Nbr_Species\n%d\n\n" % (nSpecies))
	inp.write("# VDW_Style\n%s\n\n" % (vdwStyle[i]))
	inp.write("# Charge_Style\n%s\n\n" % (chargeStyle))
	inp.write("# Intra_Scaling\n0.0 0.0 0.0 1.0\n\n")
	inp.write("# Mixing_Rule\nLB\n\n")
	inp.write("# Seed_Info\n1 2\n\n")
	inp.write("# Rcutoff_Low\n0.0\n\n")
	inp.write("# Molecule_Files\n%s\n" % (mcfStr))
	inp.write("# Box_Info\n1\ncubic\n%.2f\n\n" % (box[i]))
	inp.write("# Temperature_Info\n300.0\n\n")
	inp.write("# Move_Probability_Info\n\n")
	inp.write("# Prob_Translation\n1.0\n0.79 0.5\n\n")
	inp.write("# Done_Probability_Info\n\n")
	inp.write("# Start_Type\nread_config %s %s\n\n" % (' '.join(str(x) for x in nMolsByCheck[i]),xyzSimplefileByCheck[i]))
	inp.write("# Run_Type\nEquilibration 100\n\n")
	inp.write("# Average_Infor\n1\n\n")
	inp.write("# Simulation_Length_Info\nunits steps\nprop_freq 1\ncoord_freq 1\nrun 0\n\n")
	inp.write("END\n")
	inp.close()	

	proc = sp.Popen([cassExe + " " + inpName], stdout=sp.PIPE, shell=True)
	(out, err) = proc.communicate()

	if err is not None:
		print("Error.Abort.")

	os.system('rm '+xyzSimplefileByCheck[i])

	# 3.3) Read logfile
	log = open(cassRun + ".log", "r")
	# search line by line in log for each of the three properties
	nPrp = len(cassStr[i])
	cassAnswer[i] = [None] * nPrp
	for line in log:
		for j in range(nPrp):
			if (cassStr[i][j] in line):
				cassAnswer[i][j] = float(line.split()[-1])

	for j in range(nPrp):
		if (analyticAnswer[i][j] == 0.):
			passCheck = cassAnswer[i][j] == 0.
			if passCheck == 0:
				FailCount = FailCount+1;
				failureOutString = MainDir+ 'Scripts/testSuite/failureLog/test' + str(test_no) + '_check' + str(i+1)
				os.system('mkdir -p ' + failureOutString)
				os.system('cp ' + inpName + ' ' + failureOutString )
				os.system('cp ' + xyzfileByCheck[i] + ' ' + failureOutString )
				os.system('cp ' + ' '.join(mcfName) + ' ' + failureOutString )
				os.system('cp ' + cassRun + '*' + ' ' + failureOutString )
			if (j == 0):
				print("%-30s %-35s %18.6g %18.6g %18s %8s" % (title[i],cassPrint[i][j],
						cassAnswer[i][j],analyticAnswer[i][j],'',passCheck))
			else: 
				print("%-30s %-36s %17.6g %18.6g %18s %8s" % ('',cassPrint[i][j],
						cassAnswer[i][j],analyticAnswer[i][j],'',passCheck))

		else:
			errorRel = abs(cassAnswer[i][j] - analyticAnswer[i][j])/analyticAnswer[i][j]
			passCheck = abs(errorRel) <= errorTol
			if passCheck == 0:
				FailCount = FailCount+1;
				failureOutString = MainDir+ 'Scripts/testSuite/failureLog/test' + str(test_no) + '_check' + str(i+1)
				os.system('mkdir -p ' + failureOutString)
				os.system('cp ' + inpName + ' ' + failureOutString )
				os.system('cp ' + xyzfileByCheck[i] + ' ' + failureOutString )
				os.system('cp ' + ' '.join(mcfName) + ' ' + failureOutString )
				os.system('cp ' + cassRun + '*' + ' ' + failureOutString )
			if (j == 0):
				print("%-30s %-35s %18.6g %18.6g %18.6g %8s" % (title[i],cassPrint[i][j],
						cassAnswer[i][j],analyticAnswer[i][j],errorRel,passCheck))
			else: 
				print("%-30s %-36s %17.6g %18.6g %18.6g %8s" % ('',cassPrint[i][j],
						cassAnswer[i][j],analyticAnswer[i][j],errorRel,passCheck))

if (FailCount != 0):
	PassState = "False"
else:
	PassState = "True"

LastTest = open(MainDir+ 'Scripts/testSuite/testOutput/LastTest.txt',"w")
LastTest.write(PassState)
LastTest.close()

#*******************************************************************************
# CLEAN UP SCRATCH FILES
#*******************************************************************************
os.system('rm ' + inpName)
os.system('rm ' + ' '.join(mcfName))
os.system('rm ' + cassRun + '*')

