#!/usr/bin/env python

#This is Test 2 in a series of tests in order to check the functionality of Cassandra. 
#Test 2: This test tests 4 different configurations of two molecules in a box of size 100 Angstroms using the Mie potential equation. The first configuration is of two atoms a distance of sigma apart. The second is with each of the atoms a distance of 0.5*sigma away from the edge of opposite walls. The third test is a test where the two atoms are a distance of 2^(1/6)*sigma apart. The last test is of two atoms that are a distance of 5 angstroms apart which is out of bounds and thus the energy should be zero. 
#Note: The test uses a search to extract the energies from the log file created when running Cassandra based on the name of the energy. Furthermore, the extracted energies where than compared to the energies calculated by hand using the Mie potential equation. The exponents 14 and 6 were used.

#*******************************************************************************
# SCRIPT:  Test2_Mie_StartEnergy.py
# VERSION: 2.0
# FEATURES: Compute the Mie start energy energy between two atoms
#            (1) at a distance of 1.0 sigma (energy = 0)
#            (2) at a distance of 1.0 sigma, across period-boundary conditions
#            (3) at a distance of 2^(1/6) sigma
#			 (4) at a distance of 5 Angstrom, out of bounds
#*******************************************************************************

#*******************************************************************************
# IMPORT MODULES
#*******************************************************************************
from __future__ import print_function
import subprocess as sp #This module lets us run cassandra from python
import numpy as np #This module s the package for scientific computing in python
import random #This allows us to run random numbers
import os, sys 
import inspect
import re
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
test_no = 2
test_desc = "MIE Start Energy"

# Simulation parameters
# Parms that will be the same for every check run with this script
nSpecies = 1
nMols = (2,) # one integer per species
nAtoms = (1,) # one integer per species
atomParms = [None] * nSpecies # one list per species
s = 0
atomParms[s] = [None] * nAtoms[s] # one dictionary per atom in species
a = 0
atomParms[s][a] = { 'name':'Mie', 'type':'Mie', 'element':'Mie', 'mass':'1.0' }
atomParms[s][a]['charge'] = '0.' # store as a string
atomParms[s][a]['vdw'] = ('Mie', '120.272', '1.000', '14.000', '6.000')
bondList = [None]*nSpecies; bondParms = {} # no bonds
angleList = [None]*nSpecies; angleParms = {} # no angles
dihedList = [None]*nSpecies; dihedParms = {} # no diheds
improperList = [None]*nSpecies; improperParms = {} # no impropers
box = 100
vdwStyle = 'Mie cut_tail 7.0'
chargeStyle = 'none'

# Check parameters
# Parms that will change
nChecks = 4 # number of simulations to run
title = ("Distance 1.0 sigma", "Distance 1.0 sigma, PBC","Distance 2^(1/6) sigma", "Distance 5 A, Out of Bounds")
analyticAnswer = [None] * nChecks # list to hold analytic answers
cassAnswer     = [None] * nChecks # list to hold cassandra's answers
passCheck      = [None] * nChecks # list to hold if cassandra passed each check
cassStr = (("Total system energy",),) * nChecks # one tuple for each check
cassPrint = (("Total system energy [kJ/mol-Ext]",),) * nChecks # one tuple for each check
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
atomCoordsByCheck[0][0][0][0] = (0.,0.,0.)
atomCoordsByCheck[0][0][1][0] = (1.,0.,0.)
#check 2
atomCoordsByCheck[1][0][0][0] = (-49.5,0.,0.)
atomCoordsByCheck[1][0][1][0] = (49.5,0.,0.)
#check 3
atomCoordsByCheck[2][0][0][0] = (0.,0.,0.)
atomCoordsByCheck[2][0][1][0] = (1.12246205,0.,0.)
#check 4
atomCoordsByCheck[3][0][0][0] = (0.,0.,0.)
atomCoordsByCheck[3][0][1][0] = (5.,0.,0.)

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
cassRun 	= "test2.out"
inpName 	= "test2.inp"
xyzName 	= "test2.inp.xyz"
mcfName 	= ("test2.mcf",) # one string per species

#*******************************************************************************
# FUNCTION DEFINITIONS
#*******************************************************************************

from testSuiteFunctions import Mie

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


# Step 1) Calculate the analytical answer
analyticAnswer[0] = (Mie(1.,1.,1.,14.,6.),)
analyticAnswer[1] = (Mie(1.,1.,1.,14.,6.),)
analyticAnswer[2] = ( round(Mie(2.**(1./6.),1.,1.,14.,6.),3 ) ,)
analyticAnswer[3] = (0,)

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
		mcfStr = mcfStr + mcfName[s] + " " + str(nMols[s]) + '\n'

	# Write inp file
	inp = open(inpName,"w")
	inp.write("# Run_Name\n%s\n\n" % (cassRun))
	inp.write("# Sim_Type\nnvt\n\n")
	inp.write("# Nbr_Species\n%d\n\n" % (nSpecies))
	inp.write("# VDW_Style\n%s\n\n" % (vdwStyle))
	inp.write("# Charge_Style\n%s\n\n" % (chargeStyle))
	inp.write("# Seed_Info\n1 2\n\n")
	inp.write("# Rcutoff_Low\n1.0\n\n")
	inp.write("# Molecule_Files\n%s\n" % (mcfStr))
	inp.write("# Box_Info\n1\ncubic\n%.1f\n\n" % (box))
	inp.write("# Temperature_Info\n300.0\n\n")
	inp.write("# Move_Probability_Info\n\n")
	inp.write("# Prob_Translation\n1.0\n0.0 0.0\n\n")
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
# CLEAN UP SCRATCH FILES
#*******************************************************************************
os.system('rm ' + inpName)
os.system('rm ' + xyzName)
os.system('rm ' + ' '.join(mcfName))
os.system('rm ' + cassRun + '*')
