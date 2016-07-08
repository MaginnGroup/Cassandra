#!/usr/bin/env python

# This is Test 1 in a series of tests for a testSuite in order to check updates made to the Cassandra program. 
# Test 1 is an energy test which tests the starting energies for a molecule using the Lennard-Jones equation. This test does not perform any moves in the Cassandra program, it is simply an initial testing of energy. 
# Test 1 - Checks the starting energy for 4 different scenarios.
# Import Modules
import subprocess as sp #This module lets us run cassandra from python
import numpy as np #this module s the package for scientific computing in python
import random #This allows us to run random numbers
import os #idk what this one does yet
import re

#We will now create an input file (.inp file)
bold = '\033[1m' #Will make text bold
normal = '\033[0m' #Will make the next text normal(ie. unbold)

#This prints the starting line.
print "\n\n"+bold+"Test 1: LJ Starting Energy Commencing ... " + normal 
#Creates(opens if already existent) the file file.inp which will be our input file while running cassandra
input_inp = open("file.inp","w")
# Read in mcf and xyz files to use in inp generation
input_mcf = open("lj.mcf","w") #Creates mcf file and allows for edits
input_xyz = open("lj.xyz","w") #Creates an input xyz file (used for read_config)

#Write info for the MCF file. 
input_mcf.write("# Atom_Info\n1\n1   LJ   LJ   1.000   0.0   LJ   120.2722   1.000\n\n")
input_mcf.write("# Bond_Info\n0\n\n")
input_mcf.write("# Angle_Info\n0\n\n")
input_mcf.write("# Dihedral_Info\n0\n\n")
input_mcf.write("# Improper_Info\n0\n\n")
input_mcf.write("# Fragment_Info\n0\n\n")
input_mcf.write("# Fragment_Connectivity\n0\n\n")
input_mcf.write("# Intra_Scaling\n0. 0. 0. 0.\n0. 0. 0. 0.\n\n\n")
input_mcf.write("END")
input_mcf.close()

#Write info for the XYZ file. 
input_xyz.write("2\n\n") # This is the number of atoms in the simulation 
input_xyz.write("LJ  0.0  0.0   0.0\n") #Location of atom 1
input_xyz.write("LJ  1.0  0.0   0.0\n") #Location of atom 2
input_xyz.close()

#Write info into the file - this will create each section for the .inp file
# This input file is populated with numbers for an LJ simulation with argon
input_inp.write("# Run_Name\ntest1_check1.out\n!---------------\n\n")
input_inp.write("# Sim_Type\nNVT_MC\n!---------------\n\n")
input_inp.write("# Nbr_Species\n1\n!---------------\n\n")
input_inp.write("# VDW_Style\nLJ cut_tail 12.0\n!---------------\n\n")
input_inp.write("# Charge_Style\nNONE\n!---------------\n\n")
input_inp.write("# Intra_Scaling\n0.0 0.0 0.0 1.0\n!---------------vdw, coul lines for each species\n\n")
input_inp.write("# Mixing_Rule\nLB\n!---------------\n\n")
input_inp.write("# Seed_Info\n21498 489625\n!---------------\n\n")
input_inp.write("# Rcutoff_Low\n1.0\n!---------------\n\n")
input_inp.write("# Molecule_Files\nlj.mcf 2\n!----------------\n\n")
input_inp.write("# Box_Info\n1\nCUBIC\n100.0\n!---------------\n\n")
input_inp.write("# Temperature_Info\n300.0\n!---------------\n\n")
input_inp.write("# Move_Probability_Info\n\n")
input_inp.write("# Prob_Translation\n1.0\n1.00\n\n")
input_inp.write("# Done_Probability_Info\n!----------------\n\n")
input_inp.write("# Start_Type\nread_config 2 lj.xyz\n!---------------\n\n")
input_inp.write("# Run_Type\nEquilibration 100\n!---------------\n\n")
input_inp.write("# Average_Info\n1\n!---------------(0 = yes, 1 = no)\n\n")
input_inp.write("# Simulation_Length_Info\nUnits Steps\nProp_Freq 1\nCoord_Freq 1\nRun 0\n!---------------\n\n")
input_inp.write("# Property_Info 1\nEnergy_Total\nPressure\n!---------------\n\n")
input_inp.write("# Fragment_Files\n!---------------one line per fragment\n\n")
input_inp.write("# CBMC_Info\nkappa_ins 12\nkappa_rot 0\nkappa_dih 10\nrcut_cbmc 6.5\n!---------------\n\n")
input_inp.write("END")
input_inp.close()

# Read in input files
inp = open("file.inp").read() #This command reads in the input file
mcf = open("lj.mcf").read() #Reads in the mcf file
xyz = open("lj.xyz").read() #Reads in the orginal xyz file (used for read_config)

# The following line will print the input file when the hashtag is removed from in front of it.
#print str(inp)

# Run Cassandra Jobs
# Allows Cassandra run through Python
# This first subtest runs Cassandra with two molecules in a box of length 100, with the two molecules at a distance of sigma apart. Sigma in this case is 1.0 and epsilon is 120.0. The energy output here should match the Energy output of subtest 2, this is how we will check the accuracy of this test later on in this script.
proc = sp.Popen(["/afs/crc.nd.edu/x86_64_linux/c/cassandra/src/Cassandra_V1.2/Src/cassandra.exe " + "file.inp"], stdout=sp.PIPE, shell=True)
(out, err) = proc.communicate()

if err is not None:
	print("Error.Abort. ")

# Next, we will write over specific lines of the input file created above using a function, called replace_line that is created below. 
# This function takes three inputs: the name of the file where you would like to replace a line, the line number you would like to replace, and the text you want to replace the old text with. 
def replace_line(file_name, line_num, text):
	lines = open(file_name, 'r').readlines()
	lines[line_num] = text
	out = open(file_name, 'w')
	out.writelines(lines)
	out.close() # Closes the file so that the program doesn't explode. 

# Now that a function has been written, we can proceed by changing a line. 
# When indexing the line (ie the second input to the function, it is the line -1 that the vi open script says because python starts indexing at 0)
# This first change is for the second test that Brian, Ryan, and Eliseo want me to run (Doesn't sound as good as Ed, Ed, and Eddy) (sad.) 

#Changing the xyz input file 
# Change position of atom 1
replace_line('lj.xyz', 2, 'LJ -49.5  0.0  0.0 \n') 
# Change position of atom 2
replace_line('lj.xyz', 3, 'LJ  49.5  0.0  0.0 \n')

# Changing the input file:
# Changes the output name so we can run cassandra under a different name (and save those files too!)
replace_line('file.inp', 1, 'test1_check2.out\n') 

# Now, we will run Cassandra again, with the new numbers. 
# Run Cassandra Jobs - Again!
# Allows Cassandra run through Python
# The second subtest tests the same two molecules above, however they are now tested at a distance 0.5*sigma away from the wall of the box. This energy should match the above energy. 
proc = sp.Popen(["/afs/crc.nd.edu/x86_64_linux/c/cassandra/src/Cassandra_V1.2/Src/cassandra.exe " + "file.inp"], stdout=sp.PIPE, shell=True)
(out, err) = proc.communicate()

if err is not None:
	print("Error.Abort. ")

# Use the function from above (ie. replace_line)
# Now, we will perform the third subtest, which will test the starting energy at a minimum position, where the two molecules are a distance of (2^(1/6))*sigma apart. This should produce an energy of -1. 

# Changing the xyz input file 
# Change position of atom 1
replace_line('lj.xyz', 2, 'LJ 0.0       0.0  0.0 \n') 
# Change position of atom 2
replace_line('lj.xyz', 3, 'LJ 1.122462  0.0  0.0 \n')

# Changing the input file:
# Changes the output name so we can run cassandra under a different name (and save those files too!)
replace_line('file.inp', 1, 'test1_check3.out\n') 

# Cassandra will now be run again, for a third time
# Run Cassandra Job - Third subtest
# Allows Cassandra run through Python
proc = sp.Popen(["/afs/crc.nd.edu/x86_64_linux/c/cassandra/src/Cassandra_V1.2/Src/cassandra.exe " + "file.inp"], stdout=sp.PIPE, shell=True)
(out, err) = proc.communicate()

if err is not None:
	print("Error.Abort. ")

# The 4th and final subtest is the scenario where the two molecules are separated by a distance out of bounds of the Cassandra code and should return an energy of 0.0. In this case, we separate the two molecules by a distance of 5.0

# Changing the xyz input file 
# Change position of atom 1
replace_line('lj.xyz', 2, 'LJ 0.0  0.0  0.0 \n') 
# Change position of atom 2
replace_line('lj.xyz', 3, 'LJ 5.0  0.0  0.0 \n')

# Changing the input file:
# Changes the output name so we can run cassandra under a different name (and save those files too!)
replace_line('file.inp', 1, 'test1_check4.out\n') 
# For this edit we need to change the VDW_Style section as well. We will change it from a cut_tail 12.0 to a cut 2.5 since the position is 5.0 for the second atom. 
replace_line('file.inp', 13, 'LJ cut 2.5\n')

# Cassandra will now be run for a 4th time.
# Run Cassandra Jobs - Again!
# Allows Cassandra run through Python
proc = sp.Popen(["/afs/crc.nd.edu/x86_64_linux/c/cassandra/src/Cassandra_V1.2/Src/cassandra.exe " + "file.inp"], stdout=sp.PIPE, shell=True)
(out, err) = proc.communicate()

if err is not None:
	print("Error.Abort. ")

# For future reference :noh cancels all highlighting for when you accidently do that again

# Next, the four starting energies will be found in the output .log file for each test. This line will be extracted and saved in this python script so we can use it for comparison at the very end of the script. 
# Finding a string - using an if statement in a for loop
# For the first test
# shakes opens the desired file in the read format
shakes = open("test1_check1.out.log", "r")

# The for loop will search line by line in shakes for the words "Total system energy", once found the line will be saved as a variable.
for line in shakes:
	if re.match("(.*)Total system energy(.*)",line):
		line1 = line
 
# For the second test
# The same process is performed her as above, check commenting there if any questions.
shakes = open("test1_check2.out.log", "r")

for line in shakes:
	if re.match("(.*)Total system energy(.*)",line):
		line2 = line
  
# For the third test
# The same process is followed as in for test 1, check above for commentary. 
shakes = open("test1_check3.out.log", "r")

for line in shakes:
	if re.match("(.*)Total system energy(.*)",line):
		line3 = line
 
# For the fourth test
# The same process is followed as above, check test 1 for comments. 
shakes = open("test1_check4.out.log", "r")

for line in shakes:
	if re.match("(.*)Total system energy(.*)",line):
		line4 = line

# Next we will extract the numbers from each line of script extracted from the files above. 
# For number 1
# Set num1 equal to a blank string
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

# Now for number 2
# Set num2 equal to a blank list
num2 = []
# Use a for loop in order to go through each character in the line independently.
for t in line2.split():
# Save the number in the line to the list num2 
	try:
		num2.append(float(t))
	except ValueError:
		pass
# Extract the only number in the list as a variable
num2 = num2[0]

# Now for number 3
# Set num3 equal to a blank string 
num3 = []
# Use a for loop in order to go through each character in the line independently. 
for t in line3.split():
# Save the number in the line to the variable num3 as a float 
	try:
		num3.append(float(t))
	except ValueError:
		pass
# Extract the only number in the list as a variable
num3 = num3[0]

# Now for number 4
# Set num4 equal to a blank list
num4 = []
# Use a for loop in order to go through each character in the line independently. 
for t in line4.split():
# Save the number in the line to the variable num4 as a float
	try:
		num4.append(float(t))
	except ValueError:
		pass
# Extract the only number in the list as a variable
num4 = num4[0]

# This section is a list of all commented out things. This is because if you uncomment this section it will be easier to see why certain checks failed. Uncomment this to see the results from the log files, and the numbers which are being compared in the checks below. 
# The first four are the lines extracted from the log file. 
#print line1
#print line2
#print line3
#print line4
# The next four lines are the number extracted from the lines from the log file.
#print num1
#print num2
#print num3
#print num4

# Now we will compare and tell the user they passed and/or failed the test! 
# Check 1
if num1 == -0.00 and num1 == num2:
	c1 = 1
	print "Check 1..."
else: 
	c1 = 0 
	print "Check 1 fails"

# Now for check 2
if num2 == -0.00 and num2 == num2:
	c2 = 1
	print "Check 2..."
else: 
	c2 = 0
	print "Check 2 fails"

# Now for check three
if num3/100 == -1:
	c3 = 1
	print "Check 3..."
else:
	c3 =0
	print "Check 3 fails"

# Now for check four

if num4 == 0.0:
	c4 = 1
	print "Check 4..."
else:
	c4 = 0 
	print "Check 4 fails"  

# Now, we will see if Cassandra passes the entirety of test 1
if c1 == 1 and c2 == 1 and c3 == 1 and c4 == 1:
	print bold + "Pass Test 1: LJ Starting Energy" + normal
else:
	print bold + "Test Fails - Check above." + normal
