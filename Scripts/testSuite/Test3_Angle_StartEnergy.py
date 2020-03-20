#!/usr/bin/env python

#*******************************************************************************
# SCRIPT:  Test3_Angle_StartEnergy.py
# VERSION: 2.0
# FEATURES: Compute the angle energy using harmonic potential at
#            (1) 110 degrees
#            (2) 114 degrees
#            (3) 118 degrees
#*******************************************************************************

#*******************************************************************************
# IMPORT MODULES
#*******************************************************************************
from __future__ import print_function
import subprocess as sp #This module lets us run cassandra from python
import numpy as np #this module s the package for scientific computing in python
import random #This allows us to run random numbers
import os, sys
import re
import argparse

# if user has matplotlib, utilize it at end
try:
	import matplotlib.pyplot as pyplot
except ImportError:
	plotting=0
else:
	plotting=1

#*******************************************************************************
# IMPORT FUNCTIONS
#*******************************************************************************

from testSuiteFunctions import xyzFromAngleRandom
from testSuiteFunctions import xyzFromDistanceRandom
from testSuiteFunctions import AngleEnergyHarmonic

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
test_no = 3
test_desc = "Angle energy - Harmonic potential"


# Simulation parameters
# Parms that will be the same for every check run with this script

nSpecies = 1
nMols = (1,) # one integer per species
nAtoms = (3,) # one integer per species
atomParms = [None] * nSpecies # one list per species
s = 0
atomParms[s] = [None] * nAtoms[s] # one dictionary per atom in species
a = 0
atomParms[s][0] = { 'name':'CH3_1', 'element':'C3', 'mass':'15.034' }
atomParms[s][0]['charge'] = '0.' # store as a string
atomParms[s][0]['vdw'] = ('LJ', '98.000', '3.750')
atomParms[s][1] = { 'name':'CH2_2', 'element':'C2', 'mass':'14.027' }
atomParms[s][1]['charge'] = '0.' # store as a string
atomParms[s][1]['vdw'] = ('LJ', '46.000', '3.950')
atomParms[s][2] = { 'name':'CH3_3', 'element':'C3', 'mass':'15.034' }
atomParms[s][2]['charge'] = '0.' # store as a string
atomParms[s][2]['vdw'] = ('LJ', '98.000', '3.750')

nBonds = (2,)
bondParms = [None] * nSpecies 
bondParms[s] = [None] * nBonds[s]
bondParms[s][0] = { 'index':'1', 'molec1':'1', 'molec2':'2', 'char':'fixed', 'len':'1.540' }
bondParms[s][1] = { 'index':'2', 'molec1':'2', 'molec2':'3', 'char':'fixed', 'len':'1.540' }

# Angle values change with check
nAngles = (1,)
nChecks = 3 # number of simulations to run

angleParmsByCheck = [None] * nChecks
for c in range(nChecks):
	angleParmsByCheck[c] = [None] * nSpecies
	for s in range(nSpecies):
		angleParmsByCheck[c][s] = [None] * nMols[s]
		for m in range(nMols[s]):
			angleParmsByCheck[c][s][m] = [None] * nAngles[s]

angleParmsByCheck[0][0][0][0] = { 'index':'1', 'molec1':'1', 'molec2':'2', 'molec3':'3', 'type':'harmonic', 'force': '31250', 'nomangle': '114.00'}
angleParmsByCheck[1][0][0][0] = { 'index':'1', 'molec1':'1', 'molec2':'2', 'molec3':'3', 'type':'harmonic', 'force': '31250', 'nomangle': '114.00'}
angleParmsByCheck[2][0][0][0] = { 'index':'1', 'molec1':'1', 'molec2':'2', 'molec3':'3', 'type':'harmonic', 'force': '31250', 'nomangle': '114.00'}

dihedList = [None]*nSpecies; dihedParms = {} # no diheds
improperList = [None]*nSpecies; improperParms = {} # no impropers

box = 100.0
vdwStyle = 'lj cut_tail 14.0'
chargeStyle = 'none'

# Check parameters
# Parms that will change
nChecks = 3 # number of simulations to run
anglesToCheck = [None] * nChecks
anglesToCheck = (110.0,114.0,118.0)
title = (str(anglesToCheck[0]) + " Degrees", str(anglesToCheck[1]) + " Degrees", str(anglesToCheck[2]) + " Degrees")
analyticAnswer = [None] * nChecks # list to hold analytic answers
cassAnswer     = [None] * nChecks # list to hold cassandra's answers
passCheck      = [None] * nChecks # list to hold if cassandra passed each check
cassStr = (("Bond angle energy",),) * nChecks # one tuple for each check
cassPrint = (("Bond angle energy [kJ/mol-Ext]",),) * nChecks # one tuple for each check
errorTol = 1e-5

# the atomic coordinates are stored in set of nested lists
# each atom's coords can be accessed using the check, species, molecule and atom indices
# e.g. atomCoordsByCheck[c][s][m][a] = (x,y,z)
atomCoordsByCheck = [None] * nChecks
for c in range(nChecks):
	atomCoordsByCheck[c] = [None] * nSpecies
	for s in range(nSpecies):
		atomCoordsByCheck[c][s] = [None] * nMols[s]
		for m in range(nMols[s]):
			atomCoordsByCheck[c][s][m] = [None] * nAtoms[s]

#check 1
atomCoordsByCheck[0][0][0][0] = (0 , 0 , 0)
atomCoordsByCheck[0][0][0][1] = xyzFromDistanceRandom(atomCoordsByCheck[0][0][0][0], 1.540)
atomCoordsByCheck[0][0][0][2] = xyzFromAngleRandom(atomCoordsByCheck[0][0][0][0],atomCoordsByCheck[0][0][0][1],1.540,anglesToCheck[0])
#check 2
atomCoordsByCheck[1][0][0][0] = (0 , 0 , 0)
atomCoordsByCheck[1][0][0][1] = xyzFromDistanceRandom(atomCoordsByCheck[1][0][0][0], 1.540)
atomCoordsByCheck[1][0][0][2] = xyzFromAngleRandom(atomCoordsByCheck[1][0][0][0],atomCoordsByCheck[1][0][0][1],1.540,anglesToCheck[1])
#check 3
atomCoordsByCheck[2][0][0][0] = (0 , 0 , 0)
atomCoordsByCheck[2][0][0][1] = xyzFromDistanceRandom(atomCoordsByCheck[2][0][0][0], 1.540)
atomCoordsByCheck[2][0][0][2] = xyzFromAngleRandom(atomCoordsByCheck[2][0][0][0],atomCoordsByCheck[2][0][0][1],1.540,anglesToCheck[2])
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
cassRun 	= "test3.out"
inpName 	= "test3.inp"
xyzName 	= "test3.inp.xyz"
mcfName 	= ("test3.mcf",) # one string per species

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

# Step 1) Calculate the analytical answer
analyticAnswer[0] = (AngleEnergyHarmonic(31250,114.0,anglesToCheck[0]),)
analyticAnswer[1] = (AngleEnergyHarmonic(31250,114.0,anglesToCheck[1]),)
analyticAnswer[2] = (AngleEnergyHarmonic(31250,114.0,anglesToCheck[2]),)

print("%-30s %-35s %18s %18s %18s %8s" % ("Title", "Property","Cassandra","Analytic","Relative_Err","Pass"))
for i in range(nChecks):
	#variables that change from one check to the next
	atomCoords = atomCoordsByCheck[i] # atom coords

	# Step 2) Write input files
	# 2.1) Write MCF
	for s in range(nSpecies):
		mcf = open(mcfName[s],"w") #Creates mcf file and allows for edits
		mcf.write("# Atom_Info\n%d\n" % (nAtoms[s]))
		for j,a in enumerate(atomParms[s]):
			mcf.write("%d %s %s %s %s" % (j+1,a['name'],a['element'],a['mass'],a['charge']))
			for parm in a['vdw']:
				mcf.write(" %s" % (parm))
			mcf.write("\n")
		mcf.write("\n# Bond_Info\n%d\n" % (nBonds))
		for j,a in enumerate(bondParms[s]):
			mcf.write("%s %s %s %s %s " % (a['index'],a['molec1'],a['molec2'],a['char'],a['len']))
			mcf.write("\n")
		mcf.write("\n# Angle_Info\n%d\n" % (nAngles)) 
		for m in range(nMols[s]):
			for a in range(nAngles[s]):
				mcf.write("%s %s %s %s %s %s %s " % (angleParmsByCheck[i][s][m][a]['index'], angleParmsByCheck[i][s][m][a]['molec1'],angleParmsByCheck[i][s][m][a]['molec2'],
						angleParmsByCheck[i][s][m][a]['molec3'],angleParmsByCheck[i][s][m][a]['type'],angleParmsByCheck[i][s][m][a]['force'],angleParmsByCheck[i][s][m][a]['nomangle']))
		mcf.write("\n")
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
		mcfStr = mcfStr + mcfName[s] + " " + str(nMols[s]) + '\n'

	# Write inp file
	inp = open(inpName,"w")
	inp.write("# Run_Name\n%s\n\n" % (cassRun))
	inp.write("# Sim_Type\nnvt\n\n")
	inp.write("# Nbr_Species\n%d\n\n" % (nSpecies))
	inp.write("# VDW_Style\n%s\n\n" % (vdwStyle))
	inp.write("# Charge_Style\n%s\n\n" % (chargeStyle))
	inp.write("# Mixing_Rule\nLB\n\n")
	inp.write("# Seed_Info\n1 2\n\n")
	inp.write("# Rcutoff_Low\n1.0\n\n")
	inp.write("# Molecule_Files\n%s\n" % (mcfStr))
	inp.write("# Box_Info\n1\ncubic\n%.1f\n\n" % (box))
	inp.write("# Temperature_Info\n300.0\n\n")
	inp.write("# Move_Probability_Info\n\n")
	inp.write("# Prob_Translation\n1.0\n1.00\n\n")
	inp.write("# Done_Probability_Info\n\n")
	inp.write("# Start_Type\nread_config %s %s\n\n" % (' '.join(str(x) for x in nMols),xyzName))
	inp.write("# Run_Type\nEquilibration 100\n\n")
	inp.write("# Simulation_Length_Info\nunits steps\nprop_freq 1\ncoord_freq 1\nrun 0\n\n")
	inp.write("END")
	inp.close()

	# 2.3) Write xyz
	xyz = open(xyzName,"w")
	xyz.write("%d\n\n" % (np.dot(nMols,nAtoms))) # This is the number of atoms in the simulation 
	for s in range(nSpecies):
		for m in range(nMols[s]):
			for a in range(nAtoms[s]):
				(x,y,z) = atomCoords[s][m][a]
				xyz.write("%s %.6f %.6f %.6f\n" % (atomParms[s][a]['name'], x, y, z))
	xyz.close()


	# Step 3) Run Cassandra to get its answer
	proc = sp.Popen([cassExe + " " + inpName], stdout=sp.PIPE, shell=True)
	(out, err) = proc.communicate()

	if err is not None:
		print("Error.Abort.")


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
		if (analyticAnswer[i][j] == 0.):
			passCheck = cassAnswer[i][j] == 0.
			if passCheck == 0:
				FailCount = FailCount+1;
				failureOutString = MainDir+ 'Scripts/testSuite/failureLog/test' + str(test_no) + '_check' + str(i+1)
				os.system('mkdir -p ' + failureOutString)
				os.system('cp ' + inpName + ' ' + failureOutString )
				os.system('cp ' + xyzName + ' ' + failureOutString )
				os.system('cp ' + ' '.join(mcfName) + ' ' + failureOutString )
				os.system('cp ' + cassRun + '*' + ' ' + failureOutString )
			print("%-30s %-35s %18.6g %18.6g %18s %8s" % (title[i],cassPrint[i][j],
						cassAnswer[i][j],analyticAnswer[i][j],'',passCheck))
		else:
			errorRel = abs(cassAnswer[i][j] - analyticAnswer[i][j])/analyticAnswer[i][j]
			passCheck = abs(errorRel) <= errorTol
			if passCheck == 0:
				FailCount = FailCount+1;
				failureOutString = MainDir+ 'Scripts/testSuite/failureLog/test' + str(test_no) + '_check' + str(i+1)
				os.system('mkdir -p ' + failureOutString)
				os.system('cp ' + inpName + ' ' + failureOutString )
				os.system('cp ' + xyzName + ' ' + failureOutString )
				os.system('cp ' + ' '.join(mcfName) + ' ' + failureOutString )
				os.system('cp ' + cassRun + '*' + ' ' + failureOutString )
			print("%-30s %-35s %18.6g %18.6g %18.6g %8s" % (title[i],cassPrint[i][j],
						cassAnswer[i][j],analyticAnswer[i][j],errorRel,passCheck))


if (FailCount != 0):
	PassState = "False"
else:
	PassState = "True"

LastTest = open(MainDir+ 'Scripts/testSuite/testOutput/LastTest.txt',"w")
LastTest.write(PassState)
LastTest.close()

#*******************************************************************************
# DATA OUTPUT
#*******************************************************************************
if plotting==1:
	xPlot = [anglesToCheck[0], anglesToCheck[1], anglesToCheck[2]]
	yPlot = [cassAnswer[0][0], cassAnswer[1][0], cassAnswer[2][0]]

	pyplot.plot(xPlot,yPlot)
	pyplot.savefig(MainDir+'Scripts/testSuite/testOutput/test3_angle.png')

#*******************************************************************************
# CLEAN UP SCRATCH FILES
#*******************************************************************************
os.system('rm ' + inpName)
os.system('rm ' + xyzName)
os.system('rm ' + ' '.join(mcfName))
os.system('rm ' + cassRun + '*')
