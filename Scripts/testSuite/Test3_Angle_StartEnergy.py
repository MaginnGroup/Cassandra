#!/usr/bin/env python

# This is Test 2 in a series of tests in order to check the functionality of Cassandra. 
# Test 2 - Checks the starting energy for 4 different scenarios for the Mie potential
# Import Modules
import subprocess as sp #This module lets us run cassandra from python
import numpy as np #this module s the package for scientific computing in python
import random #This allows us to run random numbers
import os #idk what this one does yet
import inspect
import re

#We will now create an input file (.inp file)
bold = '\033[1m' #Will make text bold
normal = '\033[0m' #Will make the next text normal(ie. unbold)

#This prints the starting line
print "\n\n"+bold+"Test 2: Mie Starting Energy Commencing... " + normal 
#Creates(opens if already existent) the file file.inp which will be our input file while running cassandra
input_inp = open("mie.inp","w")
# Read in mcf and xyz files to use in inp generation
input_mcf = open("mie.mcf","w") #Creates mcf file and allows for edits
input_xyz = open("mie.xyz","w") #Creates an input xyz file (used for read_config)

#Write info for the MCF file. 
input_mcf.write("# Atom_Info\n1\n1   Mie   Mie   1.000   0.0  Mie   120.2722   1.000   14.000   6.000\n\n")
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
input_inp.write("# Run_Name\ntest1mie.out\n!---------------\n\n")
input_inp.write("# Sim_Type\nNVT_MC\n!---------------\n\n")
input_inp.write("# Nbr_Species\n1\n!---------------\n\n")
input_inp.write("# VDW_Style\nMie cut_tail 7.0\n!---------------\n\n")
input_inp.write("# Charge_Style\nNONE\n!---------------\n\n")
input_inp.write("# Intra_Scaling\n0.0 0.0 0.0 1.0\n!---------------vdw, coul lines for each species\n\n")
input_inp.write("# Mixing_Rule\nLB\n!---------------\n\n")
input_inp.write("# Seed_Info\n55001389 1764321284\n!---------------\n\n")
input_inp.write("# Rcutoff_Low\n1.0\n!---------------\n\n")
input_inp.write("# Molecule_Files\nmie.mcf 2\n!----------------\n\n")
input_inp.write("# Box_Info\n1\nCUBIC\n100.0\n!---------------\n\n")
input_inp.write("# Temperature_Info\n300.0\n!---------------\n\n")
input_inp.write("# Move_Probability_Info\n\n")
input_inp.write("# Prob_Translation\n0.79\n0.5\n\n")
input_inp.write("# Done_Probability_Info\n!----------------\n\n")
input_inp.write("# Start_Type\nread_config 2 mie.xyz\n!---------------\n\n")
input_inp.write("# Run_Type\nEquilibration 100\n!---------------\n\n")
input_inp.write("# Average_Info\n1\n!---------------(0 = yes, 1 = no)\n\n")
input_inp.write("# Simulation_Length_Info\nUnits Steps\nProp_Freq 1\nCoord_Freq 1\nRun 0\n!---------------\n\n")
input_inp.write("# Property_Info 1\nEnergy_Total\nPressure\nVolume\nNmols\nDensity\n!---------------\n\n")
input_inp.write("# Fragment_Files\n!---------------one line per fragment\n\n")
input_inp.write("# CBMC_Info\nkappa_ins 12\nkappa_rot 0\nkappa_dih 10\nrcut_cbmc 6.5\n!---------------\n\n")
input_inp.write("END")
input_inp.close()

# Read in input files
inp = open("mie.inp").read() #This command reads in the input file
mcf = open("mie.mcf").read() #Reads in the mcf file
xyz = open("mie.xyz").read() #Reads in the orginal xyz file (used for read_config)

# The following line will print the input file when the hashtag is removed from in front of it.
#print str(inp)

# Run Cassandra Jobs
#print 'Running' 
# Allows Cassandra run through Python
proc = sp.Popen(["/afs/crc.nd.edu/x86_64_linux/c/cassandra/src/Cassandra_V1.2/Src/cassandra.exe " + "mie.inp"], stdout=sp.PIPE, shell=True)
(out, err) = proc.communicate()
#print 'done'

if err is not None:
	print("Error.Abort. ")

# Now, I will attempt to write over parts of the input file in order to run cassandra again and compare the two numbers (Warning: This attempt may fail.)

# Using a function that I create for ease...(supposedly...) (The riddled sarcasm in these comments are a result of me sitting inside all day, I'm much more of a sun person. Either they will have to go, or the next person can enjoy my amazing sense of humor. 

def replace_line(file_name, line_num, text):
	lines = open(file_name, 'r').readlines()
	lines[line_num] = text
	out = open(file_name, 'w')
	out.writelines(lines)
	out.close()

# Now that a function has been defined we can attempt to change a line
# When indexing the line (ie the second input to the function, it is the line -1 that the vi open script says because python starts indexing at 0)
# This first change is for the second test that Brian, Ryan, and Eliseo want me to run (Doesn't sound as good as Ed, Ed, and Eddy) (sad.) 
# Changing the xyz input file 
# Change position of atom 1
replace_line('mie.xyz', 2, 'Mie -49.5  0.0  0.0 \n') 
# Change position of atom 2
replace_line('mie.xyz', 3, 'Mie  49.5  0.0  0.0 \n')

# Changing the input file:
# Changes the output name so we can run cassandra under a different name (and save those files too!)
replace_line('mie.inp', 1, 'test2mie.out\n') 

# Now, we will run cassandra again, under the new name. So, we will get new output files and then hopefully we will be able to compare the two starting energies and they will match!! (If not, we fail and have to start again)
# Run Cassandra Jobs - Again!
#print 'Running x2' 
# Allows Cassandra run through Python
proc = sp.Popen(["/afs/crc.nd.edu/x86_64_linux/c/cassandra/src/Cassandra_V1.2/Src/cassandra.exe " + "mie.inp"], stdout=sp.PIPE, shell=True)
(out, err) = proc.communicate()
#print 'done 2'

if err is not None:
	print("Error.Abort. ")

# And that worked cuz I'm awesome. So, now we will work on the next two tests by doing the same thing but we will create more separate output files so we can still compare. I will most likely need to come back and comment this better. Oh well. so worth putting personality in my script.

# Use the function from above (ie. replace_line)
# Change #2 will test the third test --> at minimum position. (describe this better please Lizett. thanks. and go up and describe what the other positions are later. but right now I'm on a roll!)
# Changing the xyz input file 
# Change position of atom 1
replace_line('mie.xyz', 2, 'Mie  0.0       0.0  0.0 \n') 
# Change position of atom 2
replace_line('mie.xyz', 3, 'Mie  1.122462  0.0  0.0 \n')

# Changing the input file:
# Changes the output name so we can run cassandra under a different name (and save those files too!)
replace_line('mie.inp', 1, 'test3mie.out\n') 

# Now, we will run cassandra again, under the new name. So, we will get new output files and then hopefully we will be able to compare the two starting energies and they will match!! (If not, we fail and have to start again)
# Run Cassandra Jobs - Again!
#print 'Running x3' 
# Allows Cassandra run through Python
proc = sp.Popen(["/afs/crc.nd.edu/x86_64_linux/c/cassandra/src/Cassandra_V1.2/Src/cassandra.exe " + "mie.inp"], stdout=sp.PIPE, shell=True)
(out, err) = proc.communicate()
#print 'done 3'

if err is not None:
	print("Error.Abort. ")

# Sweet! That worked too. ok last one.
# The 4th test is the 3rd scenario given to me, out of bounds which should return an energy of affectively zero. (so we shall just copy and paste again and see what happens)
# Changing the xyz input file 
# Change position of atom 1
replace_line('mie.xyz', 2, 'Mie 0.0  0.0  0.0 \n') 
# Change position of atom 2
replace_line('mie.xyz', 3, 'Mie 5.0  0.0  0.0 \n')

# Changing the input file:
# Changes the output name so we can run cassandra under a different name (and save those files too!)
replace_line('mie.inp', 1, 'test4mie.out\n') 
# For this edit we need to change the VDW_Style section as well. We will change it from a cut_tail 12.0 to a cut 2.5 since the position is 5.0 for the second atom. 
replace_line('mie.inp', 13, 'Mie cut 2.5\n')

# Now, we will run cassandra again, under the new name. So, we will get new output files and then hopefully we will be able to compare the two starting energies and they will match!! (If not, we fail and have to start again)
# Run Cassandra Jobs - Again!
#print 'Running x4' 
# Allows Cassandra run through Python
proc = sp.Popen(["/afs/crc.nd.edu/x86_64_linux/c/cassandra/src/Cassandra_V1.2/Src/cassandra.exe " + "mie.inp"], stdout=sp.PIPE, shell=True)
(out, err) = proc.communicate()
#print 'done 4'

if err is not None:
	print("Error.Abort. ")

# That is the end of the Energy tests (there are 3 scenarios, cassandra was run a total of 4 times)
# AND THEY ALL WORK!! I'M A GENIUS! (okay I can probably cut down the line which allows cassandra to run in python to one variable but I'm so proud that it's fine for now. 

# For future reference :noh cancels all highlighting for when you accidently do that again

# And now we will find the 4 numbers we need and save them! YIPEE!

# Find a specific line of text and save it. This appears to be impossible. Just kidding. Evidently not impossible 
# We will do it for each of the four tests
shakes = open('test1mie.out.log', 'r')

for line in shakes:
	if re.match("(.*)Total system energy(.*)",line):
		line1 = line
#		print line1

# For the second test
shakes = open('test2mie.out.log', 'r')

for line in shakes:
	if re.match("(.*)Total system energy(.*)",line):
		line2 = line		
#		print line2

# For the third test
shakes = open('test3mie.out.log', 'r')

for line in shakes:
	if re.match("(.*)Total system energy(.*)",line):
		line3 = line
#		print line3

# For the fourth test
shakes = open('test4mie.out.log', 'r')

for line in shakes:
	if re.match("(.*)Total system energy(.*)",line):
		line4 = line
#		print line4

# Now we will extract and compare the number in order to Pass Test 2 (each of the 4 subtests will be referred to as a check)
num1 = []
for t in line1.split():
	try:
		num1.append(float(t))
	except ValueError:
		pass
num1 = num1[0]
#print num1

# Now we will extract number 2
num2 = []
for t in line2.split():
	try:
		num2.append(float(t))
	except ValueError:
		pass
num2 = num2[0]
#print num2

# Extract number 3
num3 = []
for t in line3.split():
	try:
		num3.append(float(t))
	except ValueError:
		pass
num3 = num3[0]
#print num3

# Extract number 4
num4 = []
for t in line4.split():
	try:
		num4.append(float(t))
	except ValueError:
		pass
num4 = num4[0]
#print num4

# Now we will compare and tell the user if they passed and/or failed Test 2
# Check 1
if num1 == -0.00 and num1 == num2:
	c1 = 1
	print "Check 1..."
else:
	c1 = 0
	print "Check 1 fails"

# Check 2
if num2 == -0.00 and num1 == num2:
	c2 = 1
	print "Check 2..."
else:
	c2 = 0
	print "Check 2 fails"

# Check 3
if num3 == -119.834:
	c3 = 1
	print "Check 3..."
else:
	c3 = 0
	print "Check 3 fails"

# Check 4
if num4 == 0.00:
	c4 = 1
	print "Check 4..."
else:
	c4 = 0
	print "Check 4 fails"

#Now, we will see if Test 2 passes 
if c1 == 1 and c2 == 1 and c3 == 1 and c4 == 1:
	print "Pass Test 2: Mie Starting Energy"
else: 
	print "Test Fails - Check above."
