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
import subprocess as sp #This module lets us run cassandra from python
import numpy as np #this module s the package for scientific computing in python
import os
import re

#*******************************************************************************
# ARGUMENT PARSE
#*******************************************************************************

#*******************************************************************************
# VARIABLE DEFINITIONS
#*******************************************************************************
# Physical constants
elementary_charge = 1.60217662 # 10^-19, C
epsilon_0 = 8.85418782 # 10^-12, C^2 / N m^2
avogadro = 6.0221409 # 10^23, /mol
charge_factor = elementary_charge**2 * avogadro/(4*pi*epsilon_0) * 10000 # kJ A / mol

# Atomic parameters
charge = [1.0, 1.0]

# Check parameters
numChecks = 2
distList = [1.0,10.]
analyticList = [] #this list will hold the analytic answers
cassTestList = [] #this list will hold cassandra's answers

# Formatting variables
bold = '\033[1m' #Will make text bold
normal = '\033[0m' #Will make the next text normal(ie. unbold)

#*******************************************************************************
# FILE MANAGEMENT
#*******************************************************************************
cassDir = "/Users/rmullen2/dev/cassandra/"
cassExe = "cassandra.dev"
cassRun = "test8.out"
inpName = "inp"
xyzName = "inp.xyz"
mcfName = "inp.mcf"


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
print "\n\n"+bold+"Test 8: Coulombic energy" + normal 

# Step 1) Write input files
# 1.1) Write MCF with two point charges
mcf = open(mcfName,"w") #Creates mcf file and allows for edits
mcf.write("# Atom_Info\n2\n")
mcf.write("1    A    A   %.3f 1.0   NONE\n\n" % (charge[0]))
mcf.write("2    A    A   %.3f 1.0   NONE\n\n" % (charge[1]))
mcf.write("# Bond_Info\n0\n\n")
mcf.write("# Angle_Info\n0\n\n")
mcf.write("# Dihedral_Info\n0\n\n")
mcf.write("# Improper_Info\n0\n\n")
mcf.write("# Fragment_Info\n0\n\n")
mcf.write("# Fragment_Connectivity\n0\n\n")
mcf.write("# Intra_Scaling\n0. 0. 0. 0.\n0. 0. 0. 0.\n\n\n")
mcf.write("END")
mcf.close()

# 1.2) Write inp file
inp = open(inpName,"w")
inp.write("# Run_Name\ntest8.out\n\n")
inp.write("# Sim_Type\nnvt\n\n")
inp.write("# Nbr_Species\n1\n\n")
inp.write("# VDW_Style\nnone\n\n")
inp.write("# Charge_Style\nminimum_image\n\n")
inp.write("# Seed_Info\n1 2\n\n")
inp.write("# Rcutoff_Low\n0.1\n\n")
inp.write("# Molecule_Files\n%s 2\n\n" % (mcfName))
inp.write("# Box_Info\n1\ncubic\n100.0\n\n")
inp.write("# Temperature_Info\n0.0\n\n")
inp.write("# Move_Probability_Info\n\n")
inp.write("# Prob_Translation\n1.0\n1.00\n\n")
inp.write("# Done_Probability_Info\n\n")
inp.write("# Start_Type\nread_config 2 %s\n\n" % (xyzName))
inp.write("# Run_Type\nEquilibration 100\n\n")
inp.write("# Simulation_Length_Info\nunits steps\nprop_freq 1\ncoord_freq 1\nrun 0\n\n")
inp.write("END")
inp.close()

# Loop through checks
print "%8s %12s %12s %12s\n" % ("Distance","Casssandra","Analytic","Relative_Err")
for i in range(numChecks):
	#variables that change from one check to the next
	d = distList[i] # atomic separation

	# Step 2) Calculate the correct answer
  analyticList[i] = charge[0]*charge[1] / d * charge_factor

	# Step 3) Run Cassandra to get its answer
	# 3.1) Write xyz
	xyz = open(xyzName,"w")
	xyz.write("2\n\n") # This is the number of atoms in the simulation 
	xyz.write("A   0.0  0.0   0.0\n") #Location of atom 1
	xyz.write("A   %.1f  0.0   0.0\n" % (d)) #Location of atom 2
	xyz.close()

	# 3.2) Run Cassandra Jobs
	proc = sp.Popen([cassDir + cassExe + " " + inpName], stdout=sp.PIPE, shell=True)
	(out, err) = proc.communicate()

	if err is not None:
		print("Error.Abort. ")

	# 3.3) Read logfile
	log = open(cassRun + ".log", "r")
	# search line by line in log for the words "Total system energy"
	for line in log:
		if re.match("(.*)" + cassStr + "(.*)",line):
			lineStr = line
	    break
 
	# extract the numbers from each line of script extracted from the files above. 
	for t in lineStr.split():
		try:
			cassTestList[1] = float(t)
		except ValueError:
			pass

	# Step 4) Compare answers
	errorRel = abs(cassTestList[i] - analyticList[i])/analyticList[i]
	passCheck[i] = errorRel <= errorTol
	if passCheck[i]:
		print "%-8.1f %12.6e %12.6e %12.6e\n" % (d,cassTestList[i],analyticList[i],errorRel)

if all(passCheck):
		print bold + "Test 8: Coulombic Energy" + normal
	else:
		print bold + "Test Fails" + normal
