#!/usr/bin/env python

# This is Test 5 in a series of tests for a testSuite in order to check updates made to the Cassandra program. 
# Test 1 is an energy test which tests the starting energies for a molecule using the Lennard-Jones equation. This test does not perform any moves in the Cassandra program, it is simply an initial testing of energy. 
# Test 1 - Checks the starting energy for 4 different scenarios.
# Import Modules
import subprocess as sp #This module lets us run cassandra from python
import numpy as np #this module s the package for scientific computing in python
import random #This allows us to run random numbers
import os #idk what this one does yet
import re
import matplotlib.pyplot as pyplot

#We will now create an input file (.inp file)
bold = '\033[1m' #Will make text bold
normal = '\033[0m' #Will make the next text normal(ie. unbold)

#This prints the starting line.
print "\n\n"+bold+"Test 5: Improper Angle Energy Commencing ... " + normal 
#Creates(opens if already existent) the file file.inp which will be our input file while running cassandra
# The input files are all angle. because it deals with the angle test.
input_inp = open("improper.inp","w")
# Read in mcf and xyz files to use in inp generation
input_mcf = open("improper.mcf","w") #Creates mcf file and allows for edits
input_xyz = open("improper.xyz","w") #Creates an input xyz file (used for read_config)

#Write info for the MCF file.
input_mcf.write("# Atom_Info\n4\n1   CH3   C3   15.034   0.0   LJ   98.000   3.750\n2   CH   C   13.019   0.0   LJ   10.000   4.680\n3   CH3   C3   15.034   0.0   LJ   98.000   3.750\n4   CH3   C3   15.034   0.0   LJ   98.000   3.750\n\n")
input_mcf.write("# Bond_Info\n3\n1   1   2   fixed   1.540\n2   2   3   fixed   1.540\n3   2   4   fixed   1.540\n")
input_mcf.write("# Angle_Info\n3\n1   1   2   3   harmonic   31250.0   112.0\n2   1   2   4   harmonic   31250.0   112.0\n3   3   2   4   harmonic   31250.0   112.0\n\n")
input_mcf.write("# Dihedral_Info\n0\n\n")
input_mcf.write("# Improper_Info\n1\n1   1   2   3   harmonic   31250.0   112.0\n")
input_mcf.write("# Fragment_Info\n0\n\n")
input_mcf.write("# Fragment_Connectivity\n0\n\n")
input_mcf.write("# Intra_Scaling\n0. 0. 0. 0.\n0. 0. 0. 0.\n\n\n")
input_mcf.write("END")
input_mcf.close()

#Write info for the XYZ file. 
input_xyz.write("4\n\n") # This is the number of atoms in the simulation 
input_xyz.write("C3      0.116    0.116    0.000 \n") #Location of atom 1
input_xyz.write("C       0.687    0.832    1.239 \n") #Location of atom 2
input_xyz.write("C3      0.630   -1.336    0.000 \n") #Location of atom 3
input_xyz.write("C3     -1.423    0.086   -0.053 \n") #Location of atom 4
input_xyz.close()

#Write info into the file - this will create each section for the .inp file
# This input file is populated with numbers for an LJ simulation with argon
input_inp.write("# Run_Name\ntest1improper.out\n!---------------\n\n")
input_inp.write("# Sim_Type\nNVT_MC\n!---------------\n\n")
input_inp.write("# Nbr_Species\n1\n!---------------\n\n")
input_inp.write("# VDW_Style\nLJ cut_tail 14.0\n!---------------\n\n")
input_inp.write("# Charge_Style\nNONE\n!---------------\n\n")
input_inp.write("# Intra_Scaling\n0.0 0.0 0.0 1.0\n!---------------vdw, coul lines for each species\n\n")
input_inp.write("# Mixing_Rule\nLB\n!---------------\n\n")
input_inp.write("# Seed_Info\n21498 489625\n!---------------\n\n")
input_inp.write("# Rcutoff_Low\n1.0\n!---------------\n\n")
input_inp.write("# Molecule_Files\nimproper.mcf 1\n!----------------\n\n")
input_inp.write("# Box_Info\n1\nCUBIC\n100.0\n!---------------\n\n")
input_inp.write("# Temperature_Info\n300.0\n!---------------\n\n")
input_inp.write("# Move_Probability_Info\n\n")
input_inp.write("# Prob_Translation\n1.0\n1.00\n\n")
input_inp.write("# Done_Probability_Info\n!----------------\n\n")
input_inp.write("# Start_Type\nread_config 1 improper.xyz\n!---------------\n\n")
input_inp.write("# Run_Type\nEquilibration 100\n!---------------\n\n")
input_inp.write("# Average_Info\n1\n!---------------(0 = yes, 1 = no)\n\n")
input_inp.write("# Simulation_Length_Info\nUnits Steps\nProp_Freq 1\nCoord_Freq 1\nRun 0\n!---------------\n\n")
input_inp.write("# Property_Info 1\nEnergy_Total\n!---------------\n\n")
input_inp.write("# Fragment_Files\n!---------------one line per fragment\n\n")
input_inp.write("# CBMC_Info\nkappa_ins 12\nkappa_rot 0\nkappa_dih 10\nrcut_cbmc 6.5\n!---------------\n\n")
input_inp.write("END")
input_inp.close()

# Read in input files
inp = open("improper.inp").read() #This command reads in the input file
mcf = open("improper.mcf").read() #Reads in the mcf file
xyz = open("improper.xyz").read() #Reads in the orginal xyz file (used for read_config)

# The following line will print the input file when the hashtag is removed from in front of it.
#print str(inp)

# Run Cassandra Jobs
# Allows Cassandra run through Python
# This first subtest runs Cassandra with two molecules in a box of length 100, with the two molecules at a distance of sigma apart. Sigma in this case is 1.0 and epsilon is 120.0. The energy output here should match the Energy output of subtest 2, this is how we will check the accuracy of this test later on in this script.
print "Running"
proc = sp.Popen(["/afs/crc.nd.edu/x86_64_linux/c/cassandra/src/Cassandra_V1.2/Src/cassandra.exe " + "improper.inp"], stdout=sp.PIPE, shell=True)
(out, err) = proc.communicate()
print "done"

if err is not None:
	print("Error.Abort. ")

# Next, we will write over specific lines of the input file created above using a function, called replace_line that is created below. 
# This fucntion takes three inputs: the name of the file where you would like to replace a line, the line number you would like to replace, and the text you want to replace the old text with. 
def replace_line(file_name, line_num, text):
	lines = open(file_name, 'r').readlines()
	lines[line_num] = text
	out = open(file_name, 'w')
	out.writelines(lines)
	out.close() # Closes the file so that the program doesn't explode. 


# Finding a string - using an if statement in a for loop
# For the first test
# shakes opens the desired file in the read format
shakes = open("test1improper.out.log", "r")

# The for loop will search line by line in shakes for the words "Total system energy", once found the line will be saved as a variable.
for line in shakes:
	if re.match("(.*)Improper angle energy(.*)",line):
		line1 = line
 
# Here extract number
num1 = []
# Use a for loop in order to go through each character in the line independently.
for t in line1.split():
# Save the number in the line to the variable num1 as a float
	try:
		num1.append(float(t))
# Otherwise, do nothing.
	except ValueError:
		pass
# Extract the only number in the list as a variable.
num1 = num1[0]



# This section is a list of all commented out things. This is because if you uncomment this section it will be easier to see why certain checks failed. Uncomment this to see the results from the log files, and the numbers which are being compared in the checks below. 
# The first four are the lines extracted from the log file. 
print line1
# The next four lines are the number extracted from the lines from the log file.
print num1

# Now check to see if it passes Test 5
# Check 1
#if num1 == -0.00 and num1 == num2:
#	c1 = 1
#	print "Check 1..."
#else: 
#	c1 = 0 
#	print "Check 1 fails"

# Now, we will see if it passes all the checks, which of course it will because there is only one. 
#if c1 == 1:
#	print "Pass Test 5: Improper Angle Energy Starting "
#else:
#	print "Test Fails - Check above."

