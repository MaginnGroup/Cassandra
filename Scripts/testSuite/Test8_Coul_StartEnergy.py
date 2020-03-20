#!/usr/bin/env python

#*******************************************************************************
# SCRIPT:  Test8_Coul_StartEnergy.py
# VERSION: 2.0
# FEATURES: Compute the electrostatic energy between two point charges 
#           at a range of distances
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
# Test description
test_no = 8
test_desc = "Coulombic energy of an ion pair"

# Physical constants
# to calculate the analytic answer
elementary_charge = 1.602176565 # 10^-19, C
epsilon_0 = 8.854187817 # 10^-12, C^2 / N m^2
avogadro = 6.02214129 # 10^23, /mol
charge_factor = elementary_charge**2 * avogadro/(4*np.pi*epsilon_0) * 10000 # kJ A / mol

# Simulation parameters
nSpecies = 2
nMols = (1, 1)
nAtoms = (1, 1)
atomName = (('C',), ('A',))
atomCharge = ((1.0,), (-1.0,))
box = 39

# Check parameters
# params for each simulation
numChecks = 4 # number of simulations to run
distList = (1.0,11.0,1.0,11.)
vdwStyle = ('lj minimum_image', 'lj minimum_image', 'none', 'none')
chargeStyle = ('coul minimum_image', 'coul minimum_image', 'coul ewald 10. 1e-5', 'coul ewald 10. 1e-5')
title = ("1.0 dist [LJ/coul min]", "11.0 dist [LJ/coul min]","1.0 dist [none/coul ewald]","11.0 dist [none/coul ewald]")
cassStr = (("Total system energy",),
           ("Total system energy",),
           ("Total system energy", "Inter molecule q", "Reciprocal ewald", "Self ewald"),
           ("Total system energy", "Inter molecule q", "Reciprocal ewald", "Self ewald"))
cassPrint = (("Total system energy [kJ/mol-Ext]",),
           ("Total system energy [kJ/mol-Ext]",),
           ("Total system energy [kJ/mol-Ext]", "Inter molecule q [kJ/mol-Ext]", "Reciprocal ewald [kJ/mol-Ext]", "Self ewald [kJ/mol-Ext]"),
           ("Total system energy [kJ/mol-Ext]", "Inter molecule q [kJ/mol-Ext]", "Reciprocal ewald [kJ/mol-Ext]", "Self ewald [kJ/mol-Ext]"))
analyticAnswer = [None] * numChecks #this list will hold the analytic answers
cassAnswer     = [None] * numChecks #this list will hold cassandra's answers
passCheck      = [None] * numChecks
errorTol = 1e-5

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
cassRun = "test8.out"
inpName = "test8.inp"
xyzName = "test8.inp.xyz"
mcfName = ("test8.cation.mcf", 'test8.anion.mcf')


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
# Step 2) Calculate the correct answer
# Step 3) Run Cassandra to get its answer
# Step 4) Compare answers
#

#This prints the starting line.
print("\n\n"+bold+"Test " + str(test_no) +": " + test_desc + normal)

FailCount = 0; 

# Step 1) Write input files
# 1.1) Write MCF for cation, anion
for s in range(nSpecies):
	mcf = open(mcfName[s],"w") #Creates mcf file and allows for edits
	mcf.write("# Atom_Info\n%d\n" % (nAtoms[s]))
	for a in range(nAtoms[s]):
		mcf.write("%d    %s    %s   1.0 %.1f   NONE\n" % (a+1,atomName[s][a],atomName[s][a],atomCharge[s][a]))
	mcf.write("\n# Bond_Info\n0\n")
	mcf.write("\n# Angle_Info\n0\n")
	mcf.write("\n# Dihedral_Info\n0\n")
	mcf.write("\n# Improper_Info\n0\n")
	mcf.write("\n# Fragment_Info\n0\n")
	mcf.write("\n# Fragment_Connectivity\n0\n")
	mcf.write("\n# Intra_Scaling\n0. 0. 0. 0.\n0. 0. 0. 0.\n")
	mcf.write("\nEND\n")
	mcf.close()

# Loop through checks
print("%-30s %-35s %18s %18s %18s %8s" % ("Title", "Property","Cassandra","Analytic","Relative_Err","Pass"))
for i in range(numChecks):
	#variables that change from one check to the next
	d = distList[i] # atomic separation

	# Step 2) Calculate the correct answer
	if (chargeStyle[i] == 'coul minimum_image'):
		# Coulomb's law
		analyticAnswer[i] = [atomCharge[0][0]*atomCharge[1][0] / d * charge_factor]
	elif ('ewald' in chargeStyle[i]):
		# Ewald summation
		ewald_tol = float(chargeStyle[i].split()[-1])
		rcut = float(chargeStyle[i].split()[2])
		alpha = np.sqrt(-np.log(ewald_tol)) / rcut
		k2_cut = (np.log(ewald_tol) / rcut / np.pi)**2
		nmax = int( -2.0 * np.log(ewald_tol) / rcut * box / (2 * np.pi)) + 1

		analyticAnswer[i] = [0.] * 4
		# Real
		if (d < rcut):
			analyticAnswer[i][1] = atomCharge[0][0] * atomCharge[1][0] * charge_factor * (1 - erf(alpha * d))/d
		# Reciprocal
		recip = 0.
		for nx in range(-nmax,nmax):
			kx = float(nx) / box
			real2 = (atomCharge[0][0] + atomCharge[1][0] * np.cos(2 * np.pi * kx * d))**2
			im2 = (atomCharge[1][0] * np.sin(2 * np.pi * kx * d))**2
			for ny in range(-nmax,nmax):
				ky = float(ny) / box
				for nz in range(-nmax,nmax):
					kz = float(nz) / box
					k2 = kx**2 + ky**2 + kz**2
					if (k2 == 0. or k2 > k2_cut):
						continue
					prefactor = 1.0/float(k2) * np.exp(-(np.pi/alpha)**2 * float(k2))
					recip = recip + prefactor * (real2 + im2)
		analyticAnswer[i][2] = recip * charge_factor / (2 * np.pi * box**3)
		# Self
		analyticAnswer[i][3] = - alpha / np.sqrt(np.pi) * (atomCharge[0][0]**2 + atomCharge[1][0]**2) * charge_factor
		# Total
		analyticAnswer[i][0] = sum(analyticAnswer[i][1:4])

	# Step 3) Run Cassandra to get its answer
	# 3.1) Write inp file

	# Combine file names and nmols
	mcfStr = ''
	for s in range(nSpecies):
		mcfStr = mcfStr + mcfName[s] + " " + str(nMols[s]) + '\n'

	# Write inp file
	inp = open(inpName,"w")
	inp.write("# Run_Name\n%s\n\n" % (cassRun))
	inp.write("# Sim_Type\nnvt\n\n")
	inp.write("# Nbr_Species\n%d\n\n" % (nSpecies))
	inp.write("# VDW_Style\n%s\n\n" % (vdwStyle[i]))
	inp.write("# Charge_Style\n%s\n\n" % (chargeStyle[i]))
	inp.write("# Seed_Info\n1 2\n\n")
	inp.write("# Rcutoff_Low\n0.1\n\n")
	inp.write("# Molecule_Files\n%s\n" % (mcfStr))
	inp.write("# Box_Info\n1\ncubic\n%.1f\n\n" % (box))
	inp.write("# Temperature_Info\n1.0\n\n")
	inp.write("# Move_Probability_Info\n\n")
	inp.write("# Prob_Translation\n1.0\n0.0 0.0\n\n")
	inp.write("# Done_Probability_Info\n\n")
	inp.write("# Start_Type\nread_config %s %s\n\n" % (' '.join(str(x) for x in nMols),xyzName))
	inp.write("# Run_Type\nEquilibration 100\n\n")
	inp.write("# Simulation_Length_Info\nunits steps\nprop_freq 1\ncoord_freq 1\nrun 0\n\n")
	inp.write("END\n")
	inp.close()

	# 3.1) Write xyz
	xyz = open(xyzName,"w")
	xyz.write("2\n\n") # This is the number of atoms in the simulation 
	xyz.write("C   0.0  0.0   0.0\n") #Location of atom 1
	xyz.write("A   %.1f  0.0  0.0\n" % (d)) #Location of atom 2
	xyz.close()

	# 3.2) Run Cassandra Jobs
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
			if (j == 0):
				print("%-30s %-35s %18.6g %18.6g %18s %8s" % (title[i],cassStr[i][j],
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
				os.system('cp ' + xyzName + ' ' + failureOutString )
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

# Clean up scratch files
os.system('rm ' + inpName)
os.system('rm ' + xyzName)
os.system('rm ' + ' '.join(mcfName))
os.system('rm ' + cassRun + '*')
