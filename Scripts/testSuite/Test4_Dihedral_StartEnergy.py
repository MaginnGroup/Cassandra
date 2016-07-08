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
import matplotlib.pyplot as pyplot

#We will now create an input file (.inp file)
bold = '\033[1m' #Will make text bold
normal = '\033[0m' #Will make the next text normal(ie. unbold)

#This prints the starting line.
print "\n\n"+bold+"Test 4: Dihedral Energy Commencing ... " + normal 
#Creates(opens if already existent) the file file.inp which will be our input file while running cassandra
# The input files are all angle. because it deals with the angle test.
input_inp = open("dihedral.inp","w")
# Read in mcf and xyz files to use in inp generation
input_mcf = open("dihedral.mcf","w") #Creates mcf file and allows for edits
input_xyz = open("dihedral.xyz","w") #Creates an input xyz file (used for read_config)

#Write info for the MCF file.
input_mcf.write("# Atom_Info\n4\n1   CH3_s1   C3   15.034   0.0   LJ   98.000   3.750\n2   CH2_s1   C2   14.027   0.0   LJ   46.000   3.950\n3   CH3_s1   C2   14.027   0.0   LJ   46.000   3.950\n4   CH3_s1   C3   15.034   0.0   LJ   98.000   3.750\n\n")
input_mcf.write("# Bond_Info\n3\n1   1   2   fixed   1.540\n2   2   3   fixed   1.540\n3   3   4   fixed   1.540\n\n")
input_mcf.write("# Angle_Info\n2\n1   1   2   3   harmonic  31250.0   113.5\n2   2   3   4   harmonic   31250.0   113.5\n\n")
input_mcf.write("# Dihedral_Info\n1\n1   1   2   3   4   OPLS    0.000   2.952   -0.567   6.579\n")
input_mcf.write("# Improper_Info\n0\n\n")
input_mcf.write("# Fragment_Info\n0\n\n")
input_mcf.write("# Fragment_Connectivity\n0\n\n")
input_mcf.write("# Intra_Scaling\n0. 0. 0. 0.\n0. 0. 0. 0.\n\n\n")
input_mcf.write("END")
input_mcf.close()

#Write info for the XYZ file. 
input_xyz.write("4\n\n") # This is the number of atoms in the simulation 
input_xyz.write("C1     -3.044  -1.167   -0.047 \n") #Location of atom 1
input_xyz.write("C2     -2.530  -0.441    1.211 \n") #Location of atom 2
input_xyz.write("C3     -0.994  -0.387    1.304 \n") #Location of atom 3
input_xyz.write("C4     -0.481   0.339    2.562 \n") #Location of atom 4
input_xyz.close()

#Write info into the file - this will create each section for the .inp file
# This input file is populated with numbers for an LJ simulation with argon
input_inp.write("# Run_Name\ntest4_OPLS_check1.out\n!---------------\n\n")
input_inp.write("# Sim_Type\nNVT_MC\n!---------------\n\n")
input_inp.write("# Nbr_Species\n1\n!---------------\n\n")
input_inp.write("# VDW_Style\nLJ cut_tail 14.0\n!---------------\n\n")
input_inp.write("# Charge_Style\nNONE\n!---------------\n\n")
input_inp.write("# Intra_Scaling\n0.0 0.0 0.0 1.0\n!---------------vdw, coul lines for each species\n\n")
input_inp.write("# Mixing_Rule\nLB\n!---------------\n\n")
input_inp.write("# Seed_Info\n21498 489625\n!---------------\n\n")
input_inp.write("# Rcutoff_Low\n1.0\n!---------------\n\n")
input_inp.write("# Molecule_Files\ndihedral.mcf 1\n!----------------\n\n")
input_inp.write("# Box_Info\n1\nCUBIC\n100.0\n!---------------\n\n")
input_inp.write("# Temperature_Info\n300.0\n!---------------\n\n")
input_inp.write("# Move_Probability_Info\n\n")
input_inp.write("# Prob_Translation\n1.0\n1.00\n\n")
input_inp.write("# Done_Probability_Info\n!----------------\n\n")
input_inp.write("# Start_Type\nread_config 1 dihedral.xyz\n!---------------\n\n")
input_inp.write("# Run_Type\nEquilibration 100\n!---------------\n\n")
input_inp.write("# Average_Info\n1\n!---------------(0 = yes, 1 = no)\n\n")
input_inp.write("# Simulation_Length_Info\nUnits Steps\nProp_Freq 1\nCoord_Freq 1\nRun 0\n!---------------\n\n")
input_inp.write("# Property_Info 1\nEnergy_Total\n!---------------\n\n")
input_inp.write("# Fragment_Files\n!---------------one line per fragment\n\n")
input_inp.write("# CBMC_Info\nkappa_ins 12\nkappa_rot 0\nkappa_dih 10\nrcut_cbmc 6.5\n!---------------\n\n")
input_inp.write("END")
input_inp.close()

# Read in input files
inp = open("dihedral.inp").read() #This command reads in the input file
mcf = open("dihedral.mcf").read() #Reads in the mcf file
xyz = open("dihedral.xyz").read() #Reads in the orginal xyz file (used for read_config)

# The following line will print the input file when the hashtag is removed from in front of it.
#print str(inp)

#OPLS

# Run Cassandra Jobs
# Allows Cassandra run through Python
# This first subtest runs Cassandra with two molecules in a box of length 100, with the two molecules at a distance of sigma apart. Sigma in this case is 1.0 and epsilon is 120.0. The energy output here should match the Energy output of subtest 2, this is how we will check the accuracy of this test later on in this script.
#print "Running"
proc = sp.Popen(["/afs/crc.nd.edu/x86_64_linux/c/cassandra/src/Cassandra_V1.2/Src/cassandra.exe " + "dihedral.inp"], stdout=sp.PIPE, shell=True)
(out, err) = proc.communicate()
#print "done"

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

# Now that a function has been written, we can proceed by changing a line. 
# When indexing the line (ie the second input to the function, it is the line -1 that the vi open script says because python starts indexing at 0)
# This first change is for the second test that Brian, Ryan, and Eliseo want me to run (Doesn't sound as good as Ed, Ed, and Eddy) (sad.) 

#Changing the xyz input file
# Change position of atom 1
replace_line('dihedral.xyz', 2, 'C1    -3.066  -0.379 -0.053 \n') 
# Change position of atom 2
replace_line('dihedral.xyz', 3, 'C2    -2.553   0.347  1.205 \n')
# Change position of atom 3
replace_line('dihedral.xyz', 4, 'C3    -1.018   0.408  1.310 \n')
# Change the position of atom 4 
replace_line('dihedral.xyz', 5, 'C4    -0.504   1.136  2.566 \n')

# Change the mcf file
replace_line('dihedral.mcf', 15, "1   1   2   3   harmonic 31250.0  114.0\n")
replace_line('dihedral.mcf', 16, "2   2   3   4   harmonic 31250.0  114.0\n")
   
# Changing the input file:
# Changes the output name so we can run cassandra under a different name (and save those files too!)
replace_line('dihedral.inp', 1, 'test4_OPLS_check2.out\n') 

# Now, we will run Cassandra again, with the new numbers. 
# Run Cassandra Jobs - Again!
# Allows Cassandra run through Python
# The second subtest tests the same two molecules above, however they are now tested at a distance 0.5*sigma away from the wall of the box. This energy should match the above energy. 
proc = sp.Popen(["/afs/crc.nd.edu/x86_64_linux/c/cassandra/src/Cassandra_V1.2/Src/cassandra.exe " + "dihedral.inp"], stdout=sp.PIPE, shell=True)
(out, err) = proc.communicate()

if err is not None:
	print("Error.Abort. ")

# Use the function from above (ie. replace_line)
# Now, we will perform the third subtest, which will test the starting energy at a minimum position, where the two molecules are a distance of (2^(1/6))*sigma apart. This should produce an energy of -1. 

#Changing the xyz input file
# Change position of atom 1
replace_line('dihedral.xyz', 2, 'C1     5.923    -0.508    0.000 \n') 
# Change position of atom 2
replace_line('dihedral.xyz', 3, 'C2     4.383    -0.508    0.000 \n')
# Change position of atom 3
replace_line('dihedral.xyz', 4, 'C3     3.744     0.893    0.003 \n')
# Change the position of atom 4 
replace_line('dihedral.xyz', 5, 'C4     2.204     0.893    0.002 \n')

# Change the mcf file
replace_line('dihedral.mcf', 15, "1   1   2   3   harmonic 31250.0  114.5\n")
replace_line('dihedral.mcf', 16, "2   2   3   4   harmonic 31250.0  114.5\n")
   
# Changing the input file:
# Changes the output name so we can run cassandra under a different name (and save those files too!)
replace_line('dihedral.inp', 1, 'test4_OPLS_check3.out\n') 

# Cassandra will now be run again, for a third time
# Run Cassandra Job - Third subtest
# Allows Cassandra run through Python
proc = sp.Popen(["/afs/crc.nd.edu/x86_64_linux/c/cassandra/src/Cassandra_V1.2/Src/cassandra.exe " + "dihedral.inp"], stdout=sp.PIPE, shell=True)
(out, err) = proc.communicate()

if err is not None:
	print("Error.Abort. ")

# For future reference :noh cancels all highlighting for when you accidently do that again

# Next, the three angle energies will be found in the output .log file for each test. This line will be extracted and saved in this python script so we can use it for comparison at the very end of the script. 
# Finding a string - using an if statement in a for loop
# For the first test
# shakes opens the desired file in the read format
shakes = open("test4_OPLS_check1.out.log", "r")

# The for loop will search line by line in shakes for the words "Total system energy", once found the line will be saved as a variable.
for line in shakes:
	if re.match("(.*)Dihedral angle energy(.*)",line):
		line1 = line
 
# For the second test
# The same process is performed her as above, check commenting there if any questions.
shakes = open("test4_OPLS_check2.out.log", "r")

for line in shakes:
	if re.match("(.*)Dihedral angle energy(.*)",line):
		line2 = line
  
# For the third test
# The same process is followed as in for test 1, check above for commentary. 
shakes = open("test4_OPLS_check3.out.log", "r")

for line in shakes:
	if re.match("(.*)Dihedral angle energy(.*)",line):
		line3 = line
 

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


# This section is a list of all commented out things. This is because if you uncomment this section it will be easier to see why certain checks failed. Uncomment this to see the results from the log files, and the numbers which are being compared in the checks below. 
# The first four are the lines extracted from the log file. 
#print line1
#print line2
#print line3
# The next four lines are the number extracted from the lines from the log file.
#print num1
#print num2
#print num3

print bold +  "Checking OPLS..." + normal
# Now we will compare and tell the user they passed and/or failed the test! 
# Check 1
if abs(0.00 - num1) <= 0.05:
	c1 = 1
	print "Check 1..."
else: 
	c1 = 0 
	print "Check 1 fails"

# Now for check 2
if abs(0.00 - num2) <= 0.05:
	c2 = 1
	print "Check 2..."
else: 
	c2 = 0
	print "Check 2 fails"

# Now for check three
if abs(0.00 - num3) <= 1.00:
	c3 = 1
	print "Check 3..."
else:
	c3 =0
	print "Check 3 fails"

# Now, we will see if Cassandra passes the entirety of test 1
if c1 == 1 and c2 == 1 and c3 == 1:
	t1 = 1
	print "Pass OPLS..."
else:
	t1 = 0
	print "OPLS Fails - Check above."


#HARMONIC
## Now, we will do the same thing but run it as a harmonic
replace_line('dihedral.mcf', 20 ,'1   1   2   3   4   harmonic  260.0   114.0 \n')
replace_line('dihedral.inp', 1, 'test4_harmonic_check1.out\n') 
replace_line('dihedral.mcf', 15, "1   1   2   3   harmonic 31250.0  113.5\n")
replace_line('dihedral.mcf', 16, "2   2   3   4   harmonic 31250.0  113.5\n")
   
#Changing the xyz input file
# Change position of atom 1
replace_line('dihedral.xyz', 2, 'C1    -3.044  -1.167 -0.047 \n') 
# Change position of atom 2
replace_line('dihedral.xyz', 3, 'C2    -2.530  -0.441  1.211 \n')
# Change position of atom 3
replace_line('dihedral.xyz', 4, 'C3    -0.994  -0.387  1.304 \n')
# Change the position of atom 4 
replace_line('dihedral.xyz', 5, 'C4    -0.481   0.339  2.562 \n')


# Run Cassandra Jobs
# Allows Cassandra run through Python
# This first subtest runs Cassandra with two molecules in a box of length 100, with the two molecules at a distance of sigma apart. Sigma in this case is 1.0 and epsilon is 120.0. The energy output here should match the Energy output of subtest 2, this is how we will check the accuracy of this test later on in this script.
#print "Running"
proc = sp.Popen(["/afs/crc.nd.edu/x86_64_linux/c/cassandra/src/Cassandra_V1.2/Src/cassandra.exe " + "dihedral.inp"], stdout=sp.PIPE, shell=True)
(out, err) = proc.communicate()
#print "done"

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

# Now that a function has been written, we can proceed by changing a line. 
# When indexing the line (ie the second input to the function, it is the line -1 that the vi open script says because python starts indexing at 0)
# This first change is for the second test that Brian, Ryan, and Eliseo want me to run (Doesn't sound as good as Ed, Ed, and Eddy) (sad.) 

#Changing the xyz input file
# Change position of atom 1
replace_line('dihedral.xyz', 2, 'C1    -3.066  -0.379 -0.053 \n') 
# Change position of atom 2
replace_line('dihedral.xyz', 3, 'C2    -2.553   0.347  1.205 \n')
# Change position of atom 3
replace_line('dihedral.xyz', 4, 'C3    -1.018   0.408  1.310 \n')
# Change the position of atom 4 
replace_line('dihedral.xyz', 5, 'C4    -0.504   1.136  2.566 \n')

# Change the mcf file
replace_line('dihedral.mcf', 15, "1   1   2   3   harmonic 31250.0  114.0\n")
replace_line('dihedral.mcf', 16, "2   2   3   4   harmonic 31250.0  114.0\n")
   
# Changing the input file:
# Changes the output name so we can run cassandra under a different name (and save those files too!)
replace_line('dihedral.inp', 1, 'test4_harmonic_check2.out\n') 

# Now, we will run Cassandra again, with the new numbers. 
# Run Cassandra Jobs - Again!
# Allows Cassandra run through Python
# The second subtest tests the same two molecules above, however they are now tested at a distance 0.5*sigma away from the wall of the box. This energy should match the above energy. 
proc = sp.Popen(["/afs/crc.nd.edu/x86_64_linux/c/cassandra/src/Cassandra_V1.2/Src/cassandra.exe " + "dihedral.inp"], stdout=sp.PIPE, shell=True)
(out, err) = proc.communicate()

if err is not None:
	print("Error.Abort. ")

# Use the function from above (ie. replace_line)
# Now, we will perform the third subtest, which will test the starting energy at a minimum position, where the two molecules are a distance of (2^(1/6))*sigma apart. This should produce an energy of -1. 

#Changing the xyz input file
# Change position of atom 1
replace_line('dihedral.xyz', 2, 'C1     5.923  -0.508   0.000 \n') 
# Change position of atom 2
replace_line('dihedral.xyz', 3, 'C2     4.383  -0.508   0.000 \n')
# Change position of atom 3
replace_line('dihedral.xyz', 4, 'C3     3.744   0.893   0.003 \n')
# Change the position of atom 4 
replace_line('dihedral.xyz', 5, 'C4     2.204   0.893   0.002 \n')


# Change the mcf file
replace_line('dihedral.mcf', 15, "1   1   2   3   harmonic 31250.0  114.5\n")
replace_line('dihedral.mcf', 16, "2   2   3   4   harmonic 31250.0  114.5\n")
   
# Changing the input file:
# Changes the output name so we can run cassandra under a different name (and save those files too!)
replace_line('dihedral.inp', 1, 'test4_harmonic_check3.out\n') 

# Cassandra will now be run again, for a third time
# Run Cassandra Job - Third subtest
# Allows Cassandra run through Python
proc = sp.Popen(["/afs/crc.nd.edu/x86_64_linux/c/cassandra/src/Cassandra_V1.2/Src/cassandra.exe " + "dihedral.inp"], stdout=sp.PIPE, shell=True)
(out, err) = proc.communicate()

if err is not None:
	print("Error.Abort. ")

# For future reference :noh cancels all highlighting for when you accidently do that again

# Next, the three angle energies will be found in the output .log file for each test. This line will be extracted and saved in this python script so we can use it for comparison at the very end of the script. 
# Finding a string - using an if statement in a for loop
# For the first test
# shakes opens the desired file in the read format
shakes = open("test4_harmonic_check1.out.log", "r")

# The for loop will search line by line in shakes for the words "Total system energy", once found the line will be saved as a variable.
for line in shakes:
	if re.match("(.*)Dihedral angle energy(.*)",line):
		line1 = line
 
# For the second test
# The same process is performed her as above, check commenting there if any questions.
shakes = open("test4_harmonic_check2.out.log", "r")

for line in shakes:
	if re.match("(.*)Dihedral angle energy(.*)",line):
		line2 = line
  
# For the third test
# The same process is followed as in for test 1, check above for commentary. 
shakes = open("test4_harmonic_check3.out.log", "r")

for line in shakes:
	if re.match("(.*)Dihedral angle energy(.*)",line):
		line3 = line


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


# This section is a list of all commented out things. This is because if you uncomment this section it will be easier to see why certain checks failed. Uncomment this to see the results from the log files, and the numbers which are being compared in the checks below. 
# The first four are the lines extracted from the log file. 
#print line1
#print line2
#print line3
# The next four lines are the number extracted from the lines from the log file.
#print num1
#print num2
#print num3

print bold +  "Checking harmonic..." + normal
# Now we will compare and tell the user they passed and/or failed the test! 
# Check 1
if abs(286.8 - num1) <= 1.000: 
	c1 = 1
	print "Check 1..."
else: 
	c1 = 0 
	print "Check 1 fails"

#Now for check 2
if abs(286 - num2) <= 1.000:
	c2 = 1
	print "Check 2..."
else: 
	c2 = 0
	print "Check 2 fails"

# Now for check three
if abs(286.5 - num3) <= 1.000:
	c3 = 1
	print "Check 3..."
else:
	c3 =0
	print "Check 3 fails"

# Now, we will see if Cassandra passes the entirety of test 1
if c1 == 1 and c2 == 1 and c3 == 1:
	t2 = 1
	print "Pass harmonic..."
else:
	t2 = 0
	print "harmonic fails - Check above."


## CHARMM
replace_line('dihedral.inp', 1, 'test4_CHARMM_check1.out\n') 
replace_line('dihedral.mcf', 20, '1   1   2   3   4   CHARMM   2.952   -0.567   6.579\n')
replace_line('dihedral.mcf', 15, "1   1   2   3   harmonic 31250.0  113.5\n")
replace_line('dihedral.mcf', 16, "2   2   3   4   harmonic 31250.0  113.5\n")


#Changing the xyz input file
# Change position of atom 1
replace_line('dihedral.xyz', 2, 'C1    -3.044  -1.167 -0.047 \n') 
# Change position of atom 2
replace_line('dihedral.xyz', 3, 'C2    -2.530  -0.441  1.211 \n')
# Change position of atom 3
replace_line('dihedral.xyz', 4, 'C3    -0.994  -0.387  1.304 \n')
# Change the position of atom 4 
replace_line('dihedral.xyz', 5, 'C4    -0.481   0.339  2.562 \n')

 
# Run Cassandra Jobs
# Allows Cassandra run through Python
# This first subtest runs Cassandra with two molecules in a box of length 100, with the two molecules at a distance of sigma apart. Sigma in this case is 1.0 and epsilon is 120.0. The energy output here should match the Energy output of subtest 2, this is how we will check the accuracy of this test later on in this script.
#print "Running"
proc = sp.Popen(["/afs/crc.nd.edu/x86_64_linux/c/cassandra/src/Cassandra_V1.2/Src/cassandra.exe " + "dihedral.inp"], stdout=sp.PIPE, shell=True)
(out, err) = proc.communicate()
#print "done"

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

# Now that a function has been written, we can proceed by changing a line. 
# When indexing the line (ie the second input to the function, it is the line -1 that the vi open script says because python starts indexing at 0)
# This first change is for the second test that Brian, Ryan, and Eliseo want me to run (Doesn't sound as good as Ed, Ed, and Eddy) (sad.) 

#Changing the xyz input file
# Change position of atom 1
replace_line('dihedral.xyz', 2, 'C1    -3.066  -0.379 -0.053 \n') 
# Change position of atom 2
replace_line('dihedral.xyz', 3, 'C2    -2.553   0.347  1.205 \n')
# Change position of atom 3
replace_line('dihedral.xyz', 4, 'C3    -1.018   0.408  1.310 \n')
# Change the position of atom 4 
replace_line('dihedral.xyz', 5, 'C4    -0.504   1.136  2.566 \n')

# Change the mcf file
replace_line('dihedral.mcf', 15, "1   1   2   3   harmonic 31250.0  114.0\n")
replace_line('dihedral.mcf', 16, "2   2   3   4   harmonic 31250.0  114.0\n")
   
# Changing the input file:
# Changes the output name so we can run cassandra under a different name (and save those files too!)
replace_line('dihedral.inp', 1, 'test4_CHARMM_check2.out\n') 

# Now, we will run Cassandra again, with the new numbers. 
# Run Cassandra Jobs - Again!
# Allows Cassandra run through Python
# The second subtest tests the same two molecules above, however they are now tested at a distance 0.5*sigma away from the wall of the box. This energy should match the above energy. 
proc = sp.Popen(["/afs/crc.nd.edu/x86_64_linux/c/cassandra/src/Cassandra_V1.2/Src/cassandra.exe " + "dihedral.inp"], stdout=sp.PIPE, shell=True)
(out, err) = proc.communicate()

if err is not None:
	print("Error.Abort. ")

# Use the function from above (ie. replace_line)
# Now, we will perform the third subtest, which will test the starting energy at a minimum position, where the two molecules are a distance of (2^(1/6))*sigma apart. This should produce an energy of -1. 
#Changing the xyz input file
# Change position of atom 1
replace_line('dihedral.xyz', 2, 'C1     5.923  -0.508  0.000 \n') 
# Change position of atom 2
replace_line('dihedral.xyz', 3, 'C2     4.383  -0.508  0.000 \n')
# Change position of atom 3
replace_line('dihedral.xyz', 4, 'C3     3.744   0.893  0.003 \n')
# Change the position of atom 4 
replace_line('dihedral.xyz', 5, 'C4     2.204   0.893  0.002 \n')

# Change the mcf file
replace_line('dihedral.mcf', 15, "1   1   2   3   harmonic 31250.0  114.5\n")
replace_line('dihedral.mcf', 16, "2   2   3   4   harmonic 31250.0  114.5\n")
   
# Changing the input file:
# Changes the output name so we can run cassandra under a different name (and save those files too!)
replace_line('dihedral.inp', 1, 'test4_CHARMM_check3.out\n') 

# Cassandra will now be run again, for a third time
# Run Cassandra Job - Third subtest
# Allows Cassandra run through Python
proc = sp.Popen(["/afs/crc.nd.edu/x86_64_linux/c/cassandra/src/Cassandra_V1.2/Src/cassandra.exe " + "dihedral.inp"], stdout=sp.PIPE, shell=True)
(out, err) = proc.communicate()

if err is not None:
	print("Error.Abort. ")

# For future reference :noh cancels all highlighting for when you accidently do that again

# Next, the three angle energies will be found in the output .log file for each test. This line will be extracted and saved in this python script so we can use it for comparison at the very end of the script. 
# Finding a string - using an if statement in a for loop
# For the first test
# shakes opens the desired file in the read format
shakes = open("test4_CHARMM_check1.out.log", "r")

# The for loop will search line by line in shakes for the words "Total system energy", once found the line will be saved as a variable.
for line in shakes:
	if re.match("(.*)Dihedral angle energy(.*)",line):
		line1 = line
 
# For the second test
# The same process is performed her as above, check commenting there if any questions.
shakes = open("test4_CHARMM_check2.out.log", "r")

for line in shakes:
	if re.match("(.*)Dihedral angle energy(.*)",line):
		line2 = line
  
# For the third test
# The same process is followed as in for test 1, check above for commentary. 
shakes = open("test4_CHARMM_check3.out.log", "r")

for line in shakes:
	if re.match("(.*)Dihedral angle energy(.*)",line):
		line3 = line
 

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


# This section is a list of all commented out things. This is because if you uncomment this section it will be easier to see why certain checks failed. Uncomment this to see the results from the log files, and the numbers which are being compared in the checks below. 
# The first four are the lines extracted from the log file. 
#print line1
#print line2
#print line3
# The next four lines are the number extracted from the lines from the log file.
#print num1
#print num2
#print num3


print bold +  "Checking CHARMM..." + normal
# Now we will compare and tell the user they passed and/or failed the test! 
# Check 1
if abs(200.8 - num1) <= 1.000:
	c1 = 1
	print "Check 1..."
else: 
	c1 = 0 
	print "Check 1 fails"

# Now for check 2
if abs(201 - num2) <= 1.00:
	c2 = 1
	print "Check 2..."
else: 
	c2 = 0
	print "Check 2 fails"

# Now for check three
if abs(200.5 - num3) <= 1.00:
	c3 = 1
	print "Check 3..."
else:
	c3 =0
	print "Check 3 fails"

# Now, we will see if Cassandra passes the entirety of test 1
if c1 == 1 and c2 == 1 and c3 == 1:
	t3 = 1
	print "Pass CHARMM..."
else:
	t3 = 0
	print "CHARMM Fails - Check above."

# Check the whole test 
if t1 == 1 and t2 == 1 and t3 == 1: 
	print bold + "Pass Test 4: Dihedral Starting Angle Energy" + normal
else: 
	print bold + "Test fails - check above"+ normal



