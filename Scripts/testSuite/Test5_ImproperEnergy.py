#!/usr/bin/env python 
#*******************************************************************************
# SCRIPT:  Test5_ImproperEnergy.py
# VERSION: 2.0
# FEATURES: Compute the improper angle energy of a benzene molecule 
#            (1) At a list of ~1000 xyz positions
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
test_no = 5
test_desc = "Improper angle start energy"

# Simulation parameters
# Parms that will be the same for every check run with this script
nSpecies = 1
nMols = (1,) # one integer per species
nAtoms = (12,) # one integer per species
atomParms = [None] * nSpecies # one list per species
s = 0
atomParms[s] = [None] * nAtoms[s] # one dictionary per atom in species
a = 0
for a in range(nAtoms[s]):
	if a <=5:
		atomParms[s][a] = { 'name':'C%d_s1' % (a+1), 'element':'C', 'mass':'12.0101' , 'addParam':'ring'}
		atomParms[s][a]['charge'] = '-0.115000' # store as a string
		atomParms[s][a]['vdw'] = ('LJ', '35.226', '3.550')
	else:
		atomParms[s][a] = { 'name':'H%d_s1' % (a-5), 'element':'H', 'mass':'1.008' , 'addParam':''}
		atomParms[s][a]['charge'] = '0.115000' 
		atomParms[s][a]['vdw'] = ('LJ', '15.097', '2.420')	

nBonds = (12,)
bondParms = [None] * nSpecies 
bondParms[s] = [None] * nBonds[s]
for b in range(nBonds[s]):
	if b <= 4:
		bondParms[s][b] = { 'index':'%d' % (b+1), 'molec1':'%d' % (b+1), 'molec2':'%d' % (b+2), 'char':'fixed', 'len':'1.394'}
	elif b ==5:
		bondParms[s][b] = { 'index':'%d' % (b+1), 'molec1':'1', 'molec2':'6', 'char':'fixed', 'len':'1.394'}
	else:
		bondParms[s][b] = { 'index':'%d' % (b+1), 'molec1':'%d' % (b-5), 'molec2':'%d' % (b+1), 'char':'fixed', 'len':'1.080'}

nAngles = (18,)
angleParms = [None] * nSpecies 
angleParms[s] = [None] * nAngles[s]
for a in range(nAngles[s]):
	if a <= 3:
		angleParms[s][a] = { 'index':'%d' % (a+1), 'molec1':'%d' % (a+1), 'molec2':'%d' % (a+2),'molec3':'%d' % (a+3), 'type':'harmonic', 'force': '40257.8', 'nomangle': '120.00'}
	elif a == 4:
		angleParms[s][a] = { 'index':'%d' % (a+1), 'molec1':'%d' % (5), 'molec2':'%d' % (6),'molec3':'%d' % (1), 'type':'harmonic', 'force': '40257.8', 'nomangle': '120.00'}
	elif a == 5:
		angleParms[s][a] = { 'index':'%d' % (a+1), 'molec1':'%d' % (6), 'molec2':'%d' % (1),'molec3':'%d' % (2), 'type':'harmonic', 'force': '40257.8', 'nomangle': '120.00'}
	elif a <= 10:
		angleParms[s][a] = { 'index':'%d' % (a+1), 'molec1':'%d' % (a-5), 'molec2':'%d' % (a-4),'molec3':'%d' % (a+2), 'type':'harmonic', 'force': '30193.4', 'nomangle': '120.00'}
	elif a == 11:
		angleParms[s][a] = { 'index':'%d' % (a+1), 'molec1':'%d' % (6), 'molec2':'%d' % (1),'molec3':'%d' % (7), 'type':'harmonic', 'force': '30193.4', 'nomangle': '120.00'}
	elif a <= 16:
		angleParms[s][a] = { 'index':'%d' % (a+1), 'molec1':'%d' % (a-10), 'molec2':'%d' % (a-11),'molec3':'%d' % (a-5), 'type':'harmonic', 'force': '30193.4', 'nomangle': '120.00'}
	else:	
		angleParms[s][a] = { 'index':'%d' % (a+1), 'molec1':'%d' % (1), 'molec2':'%d' % (6),'molec3':'%d' % (12), 'type':'harmonic', 'force': '30193.4', 'nomangle': '120.00'}

# NOTE: To structure the nested list in the same way and install consistency, there are four parameters after type for each.
# For OPLS: parm1 = a0, parm2 = a1, parm3 = a2, parm4 = a3
# For harmonic: parm1 = KPhi, parm2 = Phi0, parm3/parm4 unused - set to 0
# For CHARMM: parm1 = a0, parm2 = a1, parm3 = delta, parm4 unused - set to 0
nDiheds = (24,)
dihedParms = [None] * nSpecies 
dihedParms[s] = [None] * nDiheds[s]
for d in range(nDiheds[s]):
	if d <= 2:
		dihedParms[s][d] = { 'index':'%d' % (d+1), 'molec1':'%d' % (d+1), 'molec2':'%d' % (d+2),'molec3':'%d' % (d+3),  'molec4':'%d' % (d+4),
			 'type':'CHARMM', 'parm1': '12.9700774517', 'parm2': '2.0', 'parm3': '180.0', 'parm4': '0'}
	elif d == 3:
		dihedParms[s][d] = { 'index':'%d' % (d+1), 'molec1':'%d' % (4), 'molec2':'%d' % (5), 'molec3':'%d' % (6),  'molec4':'%d' % (1),
			 'type':'CHARMM', 'parm1': '12.9700774517', 'parm2': '2.0', 'parm3': '180.0', 'parm4': '0'}
	elif d == 4:
		dihedParms[s][d] = { 'index':'%d' % (d+1), 'molec1':'%d' % (5), 'molec2':'%d' % (6),'molec3':'%d' % (1),  'molec4':'%d' % (2),
			 'type':'CHARMM', 'parm1': '12.9700774517', 'parm2': '2.0', 'parm3': '180.0', 'parm4': '0'}
	elif d == 5:
		dihedParms[s][d] = { 'index':'%d' % (d+1), 'molec1':'%d' % (6), 'molec2':'%d' % (1),'molec3':'%d' % (2),  'molec4':'%d' % (d+3),
			 'type':'CHARMM', 'parm1': '12.9700774517', 'parm2': '2.0', 'parm3': '180.0', 'parm4': '0'}
	elif d <= 9:
		dihedParms[s][d] = { 'index':'%d' % (d+1), 'molec1':'%d' % (d+1), 'molec2':'%d' % (d-5),'molec3':'%d' % (d-4),  'molec4':'%d' % (d-3),
			 'type':'CHARMM', 'parm1': '17.5731271258', 'parm2': '2.0', 'parm3': '180.0', 'parm4': '0'}
	elif d == 10:
		dihedParms[s][d] = { 'index':'%d' % (d+1), 'molec1':'%d' % (11), 'molec2':'%d' % (5),'molec3':'%d' % (6),  'molec4':'%d' % (1),
			 'type':'CHARMM', 'parm1': '17.5731271258', 'parm2': '2.0', 'parm3': '180.0', 'parm4': '0'}
	elif d == 11:
		dihedParms[s][d] = { 'index':'%d' % (d+1), 'molec1':'%d' % (12), 'molec2':'%d' % (6),'molec3':'%d' % (1),  'molec4':'%d' % (2),
			 'type':'CHARMM', 'parm1': '17.5731271258', 'parm2': '2.0', 'parm3': '180.0', 'parm4': '0'}
	elif d <= 15:
		dihedParms[s][d] = { 'index':'%d' % (d+1), 'molec1':'%d' % (d-11), 'molec2':'%d' % (d-10),'molec3':'%d' % (d-9),  'molec4':'%d' % (d-3),
			 'type':'CHARMM', 'parm1': '17.5731271258', 'parm2': '2.0', 'parm3': '180.0', 'parm4': '0'}
	elif d == 16:
		dihedParms[s][d] = { 'index':'%d' % (d+1), 'molec1':'%d' % (5), 'molec2':'%d' % (6),'molec3':'%d' % (1),  'molec4':'%d' % (7),
			 'type':'CHARMM', 'parm1': '17.5731271258', 'parm2': '2.0', 'parm3': '180.0', 'parm4': '0'}
	elif d == 17:
		dihedParms[s][d] = { 'index':'%d' % (d+1), 'molec1':'%d' % (6), 'molec2':'%d' % (1),'molec3':'%d' % (2),  'molec4':'%d' % (8),
			 'type':'CHARMM', 'parm1': '17.5731271258', 'parm2': '2.0', 'parm3': '180.0', 'parm4': '0'}
	elif d <= 22:
		dihedParms[s][d] = { 'index':'%d' % (d+1), 'molec1':'%d' % (d-11), 'molec2':'%d' % (d-17),'molec3':'%d' % (d-16),  'molec4':'%d' % (d-10),
			 'type':'CHARMM', 'parm1': '10.0420779356', 'parm2': '2.0', 'parm3': '180.0', 'parm4': '0'}
	else:
		dihedParms[s][d] = { 'index':'%d' % (d+1), 'molec1':'%d' % (12), 'molec2':'%d' % (6),'molec3':'%d' % (1),  'molec4':'%d' % (7),
			 'type':'CHARMM', 'parm1': '10.0420779356', 'parm2': '2.0', 'parm3': '180.0', 'parm4': '0'}

nImpropers = (6,)
improperParms = [None] * nSpecies 
improperParms[s] = [None] * nImpropers[s]
for i in range(nImpropers[s]):
	if i <= 3:
		improperParms[s][i] = { 'index':'%d' % (i+1), 'molec1':'%d' % (i+1), 'molec2':'%d' % (i+3),'molec3':'%d' % (i+2),  'molec4':'%d' % (i+8),
			 'type':'cvff', 'parm1': '17.5731271258', 'parm2': '-1.0', 'parm3': '2.0'}
	elif i == 4:
		improperParms[s][i] = { 'index':'%d' % (i+1), 'molec1':'%d' % (2), 'molec2':'%d' % (6),'molec3':'%d' % (1),  'molec4':'%d' % (7),
			 'type':'cvff', 'parm1': '17.5731271258', 'parm2': '-1.0', 'parm3': '2.0'}
	else:
		improperParms[s][i] = { 'index':'%d' % (i+1), 'molec1':'%d' % (1), 'molec2':'%d' % (5),'molec3':'%d' % (6),  'molec4':'%d' % (12),
			 'type':'cvff', 'parm1': '17.5731271258', 'parm2': '-1.0', 'parm3': '2.0'}

box = 50.0
vdwStyle = 'lj cut_tail 14.0'
chargeStyle = 'coul ewald 14.0 0.0000001'

# Check parameters
nChecks = 25 # number of simulations to run
title = ['None']*nChecks
for c in range(nChecks):
	title[c] = ['Test %d' % (c+1)]

analyticAnswer = [None] * nChecks # list to hold analytic answers
cassAnswer     = [None] * nChecks # list to hold cassandra's answers
passCheck      = [None] * nChecks # list to hold if cassandra passed each check
errorRel       = [None] * nChecks # list to hold err0r each check
cassSum 	   = 0
cassStr = (("Improper angle energy",),) * nChecks # one tuple for each check
cassPrint = (("Improper angle energy [kJ/mol-Ext]",),) * nChecks # one tuple for each check
errorTol = 1e-4 # Error Percent Tolerated
FailOut="n"


# Coordinates will change with each check, but will be loaded into file using a list of xyz coordinates inside
# the main loop - doing so here would require more calculations and end up hindering program and using more memory.

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
cassRun 	= "test5.out"
inpName 	= "test5.inp"
xyzName 	= "test5.inp.xyz"
mcfName 	= ("test5.mcf",) # one string per species
outputName  = MainDir + "Scripts/testSuite/testOutput/test5Results.txt"

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

FailCount = 0; # Counting number of failed tests
PassCount = 0; # Counting number of passed tests
if (FailOut=="y"):
	os.system('mkdir -p ' + MainDir + 'Scripts/testSuite/failureLog/test5/')
outText= open(outputName,"w")

k=14; # Initialize the count for xyz file input in loop

# Analytical results
givenResults = open(resourceDir + "testResults/test5_givenResults.txt", "r")
givenNumbers = [float(i) for i in givenResults]

# Distance File:
with open(resourceDir+'testConfigurations/distance.txt') as distancelist:
    distances = distancelist.readlines()

# Loop through checks
print("%-30s %-35s %18s %18s %18s %8s" % ("Passes","Failures","Pass Percent","Average Energy","Average Error","Pass"))

for i in range(nChecks):

	# Step 2) Write input files
	# 2.1) Write - only required to write on first check - mcf does not change

	for s in range(nSpecies):
		mcf = open(mcfName[s],"w") #Creates mcf file and allows for edits
		mcf.write("# Atom_Info\n%d\n" % (nAtoms[s]))
		for j,a in enumerate(atomParms[s]):	
			mcf.write("%d %s %s %s %s" % (j+1,a['name'],a['element'],a['mass'],a['charge']))
			for parm in a['vdw']:
				mcf.write(" %s" % (parm))
			mcf.write("\n")
		mcf.write("\n# Bond_Info\n%d\n" % (nBonds[s]))
		for j,a in enumerate(bondParms[s]):
			mcf.write("%s %s %s %s %s " % (a['index'],a['molec1'],a['molec2'],a['char'],str(distances[12*i+j])))
			#mcf.write("\n")
		mcf.write("\n# Angle_Info\n%d\n" % (nAngles[s]))
		for j,a in enumerate(angleParms[s]):
			mcf.write("%s %s %s %s %s %s %s " % (a['index'], a['molec1'],a['molec2'],a['molec3'],a['type'],a['force'],a['nomangle']))
			mcf.write("\n")
		mcf.write("\n# Dihedral_Info\n%d\n" % (nDiheds[s]))
		for j,a in enumerate(dihedParms[s]):
			mcf.write("%s %s %s %s %s %s %s %s %s " % (a['index'], a['molec1'],a['molec2'],a['molec3'],a['molec4'],a['type'],a['parm1'],a['parm2'],a['parm3']))
			mcf.write("\n")
		mcf.write("\n# Improper_Info\n%d\n" % (nImpropers[s]))
		for j,a in enumerate(improperParms[s]):
			mcf.write("%s %s %s %s %s %s %s %s %s " % (a['index'], a['molec1'],a['molec2'],a['molec3'],a['molec4'],a['type'],a['parm1'],a['parm2'],a['parm3']))
			mcf.write("\n")
		mcf.write("\n# Fragment_Info\n0\n")
		mcf.write("\n# Fragment_Connectivity\n0\n")
		mcf.write("\n# Intra_Scaling\n0. 0. 0.0000 1.\n0. 0. 0.0000 1.\n")
		mcf.write("\nEND\n")
		mcf.close()

	# 2.2) Write inp file - only required to write on first check - does not change
	if(i == 0):
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
		inp.write("# Intra_Scaling\n0.0 0.0 0.0 1.0\n\n")
		inp.write("# Seed_Info\n1 2\n\n")
		inp.write("# Rcutoff_Low\n1.0\n\n")
		inp.write("# Molecule_Files\n%s\n" % (mcfStr))
		inp.write("# Box_Info\n1\ncubic\n%.1f\n\n" % (box))
		inp.write("# Temperature_Info\n510.0\n\n")
		inp.write("# Move_Probability_Info\n\n")
		inp.write("# Prob_Translation\n0.40\n0.50\n14.0\n\n")
		inp.write("# Done_Probability_Info\n\n")
		inp.write("# Start_Type\nread_config %s %s\n\n" % (' '.join(str(x) for x in nMols),xyzName))
		inp.write("# Run_Type\nEquilibration 100\n\n")
		inp.write("# Simulation_Length_Info\nunits steps\nprop_freq 1\ncoord_freq 1\nrun 0\n\n")
		inp.write("# Property_Info 1\nEnergy_Total\nPressure\nVolume\nNmols\nDensity\n\n")
		inp.write("# Fragment_Files\n\n")
		inp.write("# CBMC_Info\nkappa_ins 12\nkappa_rot 0\nkappa_dih 12\nrcut_cbmc 6.5 6.5\n\n")
		inp.write("# Pair_Energy\nTRUE\n\n")
		inp.write("END")
		inp.close()	

	# 2.3) Write xyz
	xyz = open(xyzName,"w")
	with open(resourceDir + "testConfigurations/trajectoryImproper.xyz", "r") as points:
		data = iter(points)
		xyz.write(str(nAtoms[0]) + '\n\n')
		xyz.write(''.join(list(islice(data, k-12, k))))
		k = k + 14
	xyz.close()

	# Step 3) Run Cassandra to get its answer
	FNULL = open(os.devnull, 'w')
	proc = sp.Popen([cassExe + " " + inpName], stdout=sp.PIPE, stderr=FNULL, shell=True)
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

	# Write number to file
	outText.write(str(cassAnswer[i][j])+"\n")
	
	cassSum += cassAnswer[i][j]


	# Step 4) Compare answers
	for j in range(nPrp):
		if (givenNumbers[i] == 0.):
			passCheck = cassAnswer[i][j] == 0.
			if passCheck == 0:
				FailCount = FailCount+1;
				if (FailOut=="y"):
					failureOutString = MainDir+ 'Scripts/testSuite/failureLog/test5/test' + str(test_no) + '_check' + str(i+1)
					os.system('mkdir -p ' + failureOutString)
					os.system('cp ' + inpName + ' ' + failureOutString )
					os.system('cp ' + xyzName + ' ' + failureOutString )
					os.system('cp ' + ' '.join(mcfName) + ' ' + failureOutString )
					os.system('cp ' + cassRun + '*' + ' ' + failureOutString )
			else:
				PassCount = PassCount + 1

		else:
			errorRel[i] = abs(cassAnswer[i][j] - givenNumbers[i]/100)/givenNumbers[i]/100
			passCheck = abs(errorRel[i]) <= errorTol
			if passCheck == 0:
				FailCount = FailCount+1
				if (FailOut=="y"):
					failureOutString = MainDir+ 'Scripts/testSuite/failureLog/test5/test' + str(test_no) + '_check' + str(i+1)
					os.system('mkdir -p ' + failureOutString)
					os.system('cp ' + inpName + ' ' + failureOutString )
					os.system('cp ' + xyzName + ' ' + failureOutString )
					os.system('cp ' + ' '.join(mcfName) + ' ' + failureOutString )
					os.system('cp ' + cassRun + '*' + ' ' + failureOutString )			
			else:
				PassCount = PassCount + 1

if (FailCount != 0):
	PassState = "False"
else:
	PassState = "True"

LastTest = open(MainDir+ 'Scripts/testSuite/testOutput/LastTest.txt',"w")
LastTest.write(PassState)
LastTest.close()


averageError = round((float(sum(errorRel)) / nChecks),10)
averageEnergy = round(((cassSum) / nChecks),4)

print("%-30s %-35s %18s %18s %18s %8s" % (str(PassCount),str(FailCount),'%'+str(PassCount*100/(PassCount+FailCount)),str(averageEnergy),str(averageError),PassState))
#*******************************************************************************
# DATA OUTPUT
#*******************************************************************************

#*******************************************************************************
# CLEAN UP SCRATCH FILES
#*******************************************************************************
os.system('rm ' + inpName)
os.system('rm ' + xyzName)
os.system('rm ' + ' '.join(mcfName))
os.system('rm ' + cassRun + '*')








