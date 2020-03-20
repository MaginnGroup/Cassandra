#!/usr/bin/env python

#*******************************************************************************
# SCRIPT:  Test4_Dihedral_StartEnergy.py
# VERSION: 2.0
# FEATURES: Compute the dihedral angle energy using 
#            (1) OPLS at (a) 0, (b) 60, (c) 120, (d) 180
#            (2) harmonic at (a) 0, (b) 60, (c) 120, (d) 180
#            (3) CHARMM at (a) 0, (b) 60, (c) 120, (d) 180
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

from testSuiteFunctions import DihedEnergyOLPS
from testSuiteFunctions import DihedEnergyCHARMM
from testSuiteFunctions import DihedEnergyHarmonic
from testSuiteFunctions import xyzFromDistanceRandom
from testSuiteFunctions import xyzFromAngleRandom
from testSuiteFunctions import xyzFromDihedral

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
test_no = 4
test_desc = "Dihedral angle start energy"

# Simulation parameters
# Parms that will be the same for every check run with this script
nSpecies = 1
nMols = (1,) # one integer per species
nAtoms = (4,) # one integer per species
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
atomParms[s][2] = { 'name':'CH2_3', 'element':'C2', 'mass':'14.027' }
atomParms[s][2]['charge'] = '0.' # store as a string
atomParms[s][2]['vdw'] = ('LJ', '46.000', '3.950')
atomParms[s][3] = { 'name':'CH3_4', 'element':'C3', 'mass':'15.034' }
atomParms[s][3]['charge'] = '0.' # store as a string
atomParms[s][3]['vdw'] = ('LJ', '98.000', '3.750')

nBonds = (3,)
bondParms = [None] * nSpecies 
bondParms[s] = [None] * nBonds[s]
bondParms[s][0] = { 'index':'1', 'molec1':'1', 'molec2':'2', 'char':'fixed', 'len':'1.540' }
bondParms[s][1] = { 'index':'2', 'molec1':'2', 'molec2':'3', 'char':'fixed', 'len':'1.540' }
bondParms[s][2] = { 'index':'3', 'molec1':'3', 'molec2':'4', 'char':'fixed', 'len':'1.540' }

nAngles = (2,) # Angle parameters below with check
nDiheds = (1,) # Dihed parameters below with check
improperList = [None]*nSpecies; improperParms = {} # no impropers

box = 100.0
vdwStyle = 'lj cut_tail 14.0'
chargeStyle = 'none'

# Check parameters
nChecks = 12 # number of simulations to run
dihedToCheck = [None] * 4
dihedToCheck = (0.0,60.0,120.0,180.0)
typeToCheck = [None] * 3
typeToCheck = ('OPLS ', 'Harmonic ', 'CHARMM ')
title = (typeToCheck[0] + str(dihedToCheck[0])+ " degrees", typeToCheck[0] + str(dihedToCheck[1])+ " degrees", typeToCheck[0] + str(dihedToCheck[2])+ " degrees",typeToCheck[0] + str(dihedToCheck[3])+ " degrees",
	typeToCheck[1] + str(dihedToCheck[0])+ " degrees",typeToCheck[1] + str(dihedToCheck[1])+ " degrees",typeToCheck[1] + str(dihedToCheck[2])+ " degrees",typeToCheck[1] + str(dihedToCheck[3])+ " degrees",
	typeToCheck[2] + str(dihedToCheck[0])+ " degrees",typeToCheck[2] + str(dihedToCheck[1])+ " degrees",typeToCheck[2] + str(dihedToCheck[2])+ " degrees",typeToCheck[2] + str(dihedToCheck[3])+ " degrees")
analyticAnswer = [None] * nChecks # list to hold analytic answers
cassAnswer     = [None] * nChecks # list to hold cassandra's answers
passCheck      = [None] * nChecks # list to hold if cassandra passed each check
cassStr = (("Dihedral angle energy",),) * nChecks # one tuple for each check
cassPrint = (("Dihedral angle energy [kJ/mol-Ext]",),) * nChecks # one tuple for each check
errorTol = 5e-4

# Initialize
# each paramater can be accessed using the check, species, molecule and type of parameter index
# e.g. atomCoordsByCheck[c][s][m][a] = (x,y,z)
angleParmsByCheck = [None] * nChecks
dihedParmsByCheck = [None] * nChecks
atomCoordsByCheck = [None] * nChecks
for c in range(nChecks):
	angleParmsByCheck[c] = [None] * nSpecies
	dihedParmsByCheck[c] = [None] * nSpecies
	atomCoordsByCheck[c] = [None] * nSpecies
	for s in range(nSpecies):
		angleParmsByCheck[c][s] = [None] * nMols[s]
		dihedParmsByCheck[c][s] = [None] * nMols[s]
		atomCoordsByCheck[c][s] = [None] * nMols[s]
		for m in range(nMols[s]):
			angleParmsByCheck[c][s][m] = [None] * nAngles[s]
			dihedParmsByCheck[c][s][m] = [None] * nDiheds[s]
			atomCoordsByCheck[c][s][m] = [None] * nAtoms[s]

# Angle values (two per check, same all check)
angleParmsByCheck[0][0][0][0] = { 'index':'1', 'molec1':'1', 'molec2':'2', 'molec3':'3', 'type':'harmonic', 'force': '31250.0', 'nomangle': '114.00'}
angleParmsByCheck[0][0][0][1] = { 'index':'2', 'molec1':'2', 'molec2':'3', 'molec3':'4', 'type':'harmonic', 'force': '31250.0', 'nomangle': '114.00'}
angleParmsByCheck[1] = angleParmsByCheck[0]
angleParmsByCheck[2] = angleParmsByCheck[0]
angleParmsByCheck[3] = angleParmsByCheck[0]
angleParmsByCheck[4] = angleParmsByCheck[0]
angleParmsByCheck[5] = angleParmsByCheck[0]
angleParmsByCheck[6] = angleParmsByCheck[0]
angleParmsByCheck[7] = angleParmsByCheck[0]
angleParmsByCheck[8] = angleParmsByCheck[0]
angleParmsByCheck[9] = angleParmsByCheck[0]
angleParmsByCheck[10] = angleParmsByCheck[0]
angleParmsByCheck[11] = angleParmsByCheck[0]

# Dihed values change with check (change every three simulations)
# NOTE: To structure the nested list in the same way, there must be four parameters after type for each.
# For OPLS: parm1 = a0, parm2 = a1, parm3 = a2, parm4 = a3
# For harmonic: parm1 = KPhi, parm2 = Phi0, parm3/parm4 unused - set to 0
# For CHARMM: parm1 = 1st set of a0, a1, delta, parm2 = second set, parm3 = third set
nDiheds = (1,) 
dihedParmsByCheck[0][0][0][0] = { 'index':'1', 'molec1':'1', 'molec2':'2', 'molec3':'3', 'molec4':'4', 'type':'OPLS', 'parm1': '0.000', 'parm2': '2.952', 'parm3': '-0.567', 'parm4': '6.579'}
dihedParmsByCheck[1][0][0][0]= dihedParmsByCheck[0][0][0][0]
dihedParmsByCheck[2][0][0][0]= dihedParmsByCheck[0][0][0][0]
dihedParmsByCheck[3][0][0][0]= dihedParmsByCheck[0][0][0][0]
dihedParmsByCheck[4][0][0][0] = { 'index':'1', 'molec1':'1', 'molec2':'2', 'molec3':'3', 'molec4':'4', 'type':'harmonic', 'parm1': '232.4', 'parm2': '180.0', 'parm3': '0.000', 'parm4': '0.000'}
dihedParmsByCheck[5][0][0][0]= dihedParmsByCheck[4][0][0][0]
dihedParmsByCheck[6][0][0][0]= dihedParmsByCheck[4][0][0][0]
dihedParmsByCheck[7][0][0][0]= dihedParmsByCheck[4][0][0][0]
dihedParmsByCheck[8][0][0][0] = { 'index':'1', 'molec1':'1', 'molec2':'2', 'molec3':'3', 'molec4':'4', 'type':'CHARMM', 'parm1': '2.952 1.0 0.0', 'parm2': '-0.567 2.0 180.0', 'parm3': '6.579 3.0 0.0'}
dihedParmsByCheck[9][0][0][0]= dihedParmsByCheck[8][0][0][0]
dihedParmsByCheck[10][0][0][0]= dihedParmsByCheck[8][0][0][0]
dihedParmsByCheck[11][0][0][0]= dihedParmsByCheck[8][0][0][0]

# Coordinate values change with check (change every time, cycle after three)
#check 1
atomCoordsByCheck[0][0][0][0] = (0,0,0)
atomCoordsByCheck[0][0][0][1] = xyzFromDistanceRandom(atomCoordsByCheck[0][0][0][0],1.540)
atomCoordsByCheck[0][0][0][2] = xyzFromAngleRandom(atomCoordsByCheck[0][0][0][0],atomCoordsByCheck[0][0][0][1],1.540,114.0)
atomCoordsByCheck[0][0][0][3] = xyzFromDihedral(atomCoordsByCheck[0][0][0][0],atomCoordsByCheck[0][0][0][1],atomCoordsByCheck[0][0][0][2],1.540,114.0,dihedToCheck[0])
#check 2
atomCoordsByCheck[1][0][0][0] = (0,0,0)
atomCoordsByCheck[1][0][0][1] = xyzFromDistanceRandom(atomCoordsByCheck[1][0][0][0],1.540)
atomCoordsByCheck[1][0][0][2] = xyzFromAngleRandom(atomCoordsByCheck[1][0][0][0],atomCoordsByCheck[1][0][0][1],1.540,114.0)
atomCoordsByCheck[1][0][0][3] = xyzFromDihedral(atomCoordsByCheck[1][0][0][0],atomCoordsByCheck[1][0][0][1],atomCoordsByCheck[1][0][0][2],1.540,114.0,dihedToCheck[1])
#check 3
atomCoordsByCheck[2][0][0][0] = (0,0,0)
atomCoordsByCheck[2][0][0][1] = xyzFromDistanceRandom(atomCoordsByCheck[2][0][0][0],1.540)
atomCoordsByCheck[2][0][0][2] = xyzFromAngleRandom(atomCoordsByCheck[2][0][0][0],atomCoordsByCheck[2][0][0][1],1.540,114.0)
atomCoordsByCheck[2][0][0][3] = xyzFromDihedral(atomCoordsByCheck[2][0][0][0],atomCoordsByCheck[2][0][0][1],atomCoordsByCheck[2][0][0][2],1.540,114.0,dihedToCheck[2])
# check 4 
atomCoordsByCheck[3][0][0][0] = (0,0,0)
atomCoordsByCheck[3][0][0][1] = xyzFromDistanceRandom(atomCoordsByCheck[3][0][0][0],1.540)
atomCoordsByCheck[3][0][0][2] = xyzFromAngleRandom(atomCoordsByCheck[3][0][0][0],atomCoordsByCheck[3][0][0][1],1.540,114.0)
atomCoordsByCheck[3][0][0][3] = xyzFromDihedral(atomCoordsByCheck[3][0][0][0],atomCoordsByCheck[3][0][0][1],atomCoordsByCheck[3][0][0][2],1.540,114.0,dihedToCheck[3])
# check 5 - 12
atomCoordsByCheck[4] = atomCoordsByCheck[0]
atomCoordsByCheck[5] = atomCoordsByCheck[1]
atomCoordsByCheck[6] = atomCoordsByCheck[2]
atomCoordsByCheck[7] = atomCoordsByCheck[3]

atomCoordsByCheck[8] = atomCoordsByCheck[0]
atomCoordsByCheck[9] = atomCoordsByCheck[1]
atomCoordsByCheck[10] = atomCoordsByCheck[2]
atomCoordsByCheck[11] = atomCoordsByCheck[3]
 

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
cassRun 	= "test4.out"
inpName 	= "test4.inp"
xyzName 	= "test4.inp.xyz"
mcfName 	= ("test4.mcf",) # one string per species

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

# Step 1) Input the analytical answers
analyticAnswer[0] = (DihedEnergyOLPS(0.0, 2.952, -0.567,  6.579, dihedToCheck[0]),)
analyticAnswer[1] = (DihedEnergyOLPS(0.0, 2.952, -0.567,  6.579, dihedToCheck[1]),)
analyticAnswer[2] = (DihedEnergyOLPS(0.0, 2.952, -0.567,  6.579, dihedToCheck[2]),)
analyticAnswer[3] = (DihedEnergyOLPS(0.0, 2.952, -0.567,  6.579, dihedToCheck[3]),)
analyticAnswer[4] = (DihedEnergyHarmonic(232.4, 180.0, dihedToCheck[0]),)
analyticAnswer[5] = (DihedEnergyHarmonic(232.4, 180.0, dihedToCheck[1]),)
analyticAnswer[6] = (DihedEnergyHarmonic(232.4, 180.0, dihedToCheck[2]),)
analyticAnswer[7] = (DihedEnergyHarmonic(232.4, 180.0, dihedToCheck[3]),)
analyticAnswer[8] = (DihedEnergyCHARMM(2.952, 1.0, 0.0, dihedToCheck[0])+DihedEnergyCHARMM(-0.567, 2.0, 180.0, dihedToCheck[0])+DihedEnergyCHARMM(6.579, 3.0, 0.0, dihedToCheck[0]),)
analyticAnswer[9] = (DihedEnergyCHARMM(2.952, 1.0, 0.0, dihedToCheck[1])+DihedEnergyCHARMM(-0.567, 2.0, 180.0, dihedToCheck[1])+DihedEnergyCHARMM(6.579, 3.0, 0.0, dihedToCheck[1]),)
analyticAnswer[10] = (DihedEnergyCHARMM(2.952, 1.0, 0.0, dihedToCheck[2])+DihedEnergyCHARMM(-0.567, 2.0, 180.0, dihedToCheck[2])+DihedEnergyCHARMM(6.579, 3.0, 0.0, dihedToCheck[2]),)
analyticAnswer[11] = (DihedEnergyCHARMM(2.952, 1.0, 0.0, dihedToCheck[3])+DihedEnergyCHARMM(-0.567, 2.0, 180.0, dihedToCheck[3])+DihedEnergyCHARMM(6.579, 3.0, 0.0, dihedToCheck[3]),)

# Loop through checks
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
		mcf.write("\n# Bond_Info\n%d\n" % (nBonds[s]))
		for j,a in enumerate(bondParms[s]):
			mcf.write("%s %s %s %s %s " % (a['index'],a['molec1'],a['molec2'],a['char'],a['len']))
			mcf.write("\n")
		mcf.write("\n# Angle_Info\n%d\n" % (nAngles[s]))
		for m in range(nMols[s]):
			for a in range(nAngles[s]):
				mcf.write("%s %s %s %s %s %s %s " % (angleParmsByCheck[i][s][m][a]['index'], angleParmsByCheck[i][s][m][a]['molec1'],angleParmsByCheck[i][s][m][a]['molec2'],
						angleParmsByCheck[i][s][m][a]['molec3'],angleParmsByCheck[i][s][m][a]['type'],angleParmsByCheck[i][s][m][a]['force'],angleParmsByCheck[i][s][m][a]['nomangle']))
				mcf.write("\n")
		if (i<8):
			mcf.write("\n# Dihedral_Info\n%d\n" % (nDiheds[s]))
		else:
			mcf.write("\n# Dihedral_Info\n%d\n" % (3))
		for m in range(nMols[s]):
			for a in range(nDiheds[s]):
				if (dihedParmsByCheck[i][s][m][a]['type'] == 'OPLS'):
					mcf.write("%s %s %s %s %s %s %s %s %s %s " % (dihedParmsByCheck[i][s][m][a]['index'], dihedParmsByCheck[i][s][m][a]['molec1'],dihedParmsByCheck[i][s][m][a]['molec2'],
						dihedParmsByCheck[i][s][m][a]['molec3'],dihedParmsByCheck[i][s][m][a]['molec4'],dihedParmsByCheck[i][s][m][a]['type'],dihedParmsByCheck[i][s][m][a]['parm1'],
						dihedParmsByCheck[i][s][m][a]['parm2'],dihedParmsByCheck[i][s][m][a]['parm3'],dihedParmsByCheck[i][s][m][a]['parm4']))
				elif (dihedParmsByCheck[i][s][m][a]['type'] == 'harmonic'):
					mcf.write("%s %s %s %s %s %s %s %s " % (dihedParmsByCheck[i][s][m][a]['index'], dihedParmsByCheck[i][s][m][a]['molec1'],dihedParmsByCheck[i][s][m][a]['molec2'],
						dihedParmsByCheck[i][s][m][a]['molec3'],dihedParmsByCheck[i][s][m][a]['molec4'],dihedParmsByCheck[i][s][m][a]['type'],dihedParmsByCheck[i][s][m][a]['parm1'],
						dihedParmsByCheck[i][s][m][a]['parm2']))
				elif (dihedParmsByCheck[i][s][m][a]['type'] == 'CHARMM'):
					mcf.write("%s %s %s %s %s %s %s \n" % (1, dihedParmsByCheck[i][s][m][a]['molec1'],dihedParmsByCheck[i][s][m][a]['molec2'],
						dihedParmsByCheck[i][s][m][a]['molec3'],dihedParmsByCheck[i][s][m][a]['molec4'],dihedParmsByCheck[i][s][m][a]['type'],dihedParmsByCheck[i][s][m][a]['parm1']))
					mcf.write("%s %s %s %s %s %s %s \n" % (2, dihedParmsByCheck[i][s][m][a]['molec1'],dihedParmsByCheck[i][s][m][a]['molec2'],
						dihedParmsByCheck[i][s][m][a]['molec3'],dihedParmsByCheck[i][s][m][a]['molec4'],dihedParmsByCheck[i][s][m][a]['type'],dihedParmsByCheck[i][s][m][a]['parm2']))
					mcf.write("%s %s %s %s %s %s %s" % (3, dihedParmsByCheck[i][s][m][a]['molec1'],dihedParmsByCheck[i][s][m][a]['molec2'],
						dihedParmsByCheck[i][s][m][a]['molec3'],dihedParmsByCheck[i][s][m][a]['molec4'],dihedParmsByCheck[i][s][m][a]['type'],dihedParmsByCheck[i][s][m][a]['parm3']))
			mcf.write("\n")
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
	inp.write("# Intra_Scaling\n0.0 0.0 0.0 1.0\n\n")
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
	inp.write("# Property_Info 1\nEnergy_Total\nPressure\nVolume\nNmols\nDensity\n!---------------\n\n")
	inp.write("# Fragment_Files\n\n")
	inp.write("# CBMC_Info\nkappa_ins 12\nkappa_rot 0\nkappa_dih 10\nrcut_cbmc 6.5\n\n")
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
			print("%-30s %-36s %17.6g %18.6g %18s %8s" % (title[i],cassPrint[i][j],
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
			print("%-30s %-36s %17.6g %18.6g %18.6g %8s" % (title[i],cassPrint[i][j],
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

if plotting:
	xPlot = [0.0, 60.0, 120.0, 180.0]
	yPlot1 = [cassAnswer[0][0], cassAnswer[1][0], cassAnswer[2][0],cassAnswer[3][0]]
	yPlot2 = [cassAnswer[4][0], cassAnswer[5][0], cassAnswer[6][0],cassAnswer[7][0]]
	yPlot3 = [cassAnswer[8][0], cassAnswer[9][0], cassAnswer[10][0],cassAnswer[11][0]]

	pyplot.figure(1)
	pyplot.plot(xPlot,yPlot1)
	pyplot.savefig(MainDir+'Scripts/testSuite/testOutput/test4_OLPS.png')

	pyplot.figure(2)
	pyplot.plot(xPlot,yPlot2)
	pyplot.savefig(MainDir+'Scripts/testSuite/testOutput/test4_harmonic.png')

	pyplot.figure(3)
	pyplot.plot(xPlot,yPlot3)
	pyplot.savefig(MainDir+'Scripts/testSuite/testOutput/test4_CHARMM.png')

#*******************************************************************************
# CLEAN UP SCRATCH FILES
#*******************************************************************************
os.system('rm ' + inpName)
os.system('rm ' + xyzName)
os.system('rm ' + ' '.join(mcfName))
os.system('rm ' + cassRun + '*')
