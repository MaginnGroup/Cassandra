#!/usr/bin/env python

# This is Test 2 in a series of tests in order to check the functionality of Cassandra. 
#Test 2: This test tests 4 different configurations of two molecules in a box of size 100 Angstroms using the Mie potential equation. The first configuration is of two atoms a distance of sigma apart. The second is with each of the atoms a distance of 0.5*sigma away from the edge of opposite walls. The third test is a test where the two atoms are a distance of 2^(1/6)*sigma apart. The last test is of two atoms that are a distance of 5 angstroms apart which is out of bounds and thus the energy should be zero. 
#Note: The test uses a search to extract the energies from the log file created when running Cassandra based on the name of the energy. Furthermore, the extracted energies where than compared to the energies calculated by hand using the Mie potential equation. The exponents 14 and 6 were used.


# Test 2 - Checks the starting energy for 4 different scenarios for the Mie potential
# Import Modules
import subprocess as sp #This module lets us run cassandra from python
import numpy as np #this module s the package for scientific computing in python
import random #This allows us to run random numbers
import os #idk what this one does yet
import inspect
import re

# This command allows text to be bolded and unbolded in the command line. 
bold = '\033[1m' #Will make text bold
normal = '\033[0m' #Will make the next text normal(ie. unbold)

#This prints the starting line, telling the user which test is running
print "\n\n"+bold+"Test 2: Mie Starting Energy Commencing... " + normal 
#Creates(opens if already existent) the file mie.inp which will be our input file while running cassandra
input_inp = open("mie.inp","w")
# Read in mcf and xyz files to use in inp generation
input_mcf = open("mie.mcf","w") #Creates mcf file and allows for edits
input_xyz = open("mie.xyz","w") #Creates an input xyz file (used for read_config)

#Write info for the MCF file. 
# The atom info in this case coincides with the Mie potential. These numbers take the format of other mcf files, a sample mcf file was used for reference on structuring.
input_mcf.write("# Atom_Info\n1\n1   Mie   Mie   1.000   0.0  Mie   120.2722   1.000   14.000   6.000\n\n")
input_mcf.write("# Bond_Info\n0\n\n")
input_mcf.write("# Angle_Info\n0\n\n")
input_mcf.write("# Dihedral_Info\n0\n\n")
input_mcf.write("# Improper_Info\n0\n\n")
input_mcf.write("# Fragment_Info\n0\n\n")
input_mcf.write("# Fragment_Connectivity\n0\n\n")
input_mcf.write("# Intra_Scaling\n0. 0. 0. 0.\n0. 0. 0. 0.\n\n\n")
input_mcf.write("END")
# Closes the file to ensure that Cassandra does not explode.
input_mcf.close()

# Write info for the XYZ file. 
input_xyz.write("2\n\n") # This is the number of atoms in the simulation 
input_xyz.write("LJ  0.0  0.0   0.0\n") #Location of atom 1
input_xyz.write("LJ  1.0  0.0   0.0\n") #Location of atom 2
input_xyz.close()

# Write info into the file - this will create each section for the .inp file
# This input file is populated with numbers for an Mie simulation with an arbitrary atom
# This input decides what the output files will be named, it should be unique so as not to overwrite existing files. 
input_inp.write("# Run_Name\ntest2_check1.out\n!---------------\n\n")
# This line decides the simulation type. In this case, we are running an NVT Monte Carlo, that is constant mols, volume, and temperature.
input_inp.write("# Sim_Type\nNVT_MC\n!---------------\n\n")
# The number of species in the simulation
input_inp.write("# Nbr_Species\n1\n!---------------\n\n")
# This determines where to cut off the very inaccurate part of the spectrum
input_inp.write("# VDW_Style\nMie cut_tail 7.0\n!---------------\n\n")
# There is no charge, so NONE is imputted 
input_inp.write("# Charge_Style\nNONE\n!---------------\n\n")
input_inp.write("# Intra_Scaling\n0.0 0.0 0.0 1.0\n!---------------vdw, coul lines for each species\n\n")
input_inp.write("# Mixing_Rule\nLB\n!---------------\n\n")
# Seed info can be the generation of any random number 
input_inp.write("# Seed_Info\n55001389 1764321284\n!---------------\n\n")
input_inp.write("# Rcutoff_Low\n1.0\n!---------------\n\n")
# This line tells us from which mcf file we will be reading, and how many molecules should be read in.
input_inp.write("# Molecule_Files\nmie.mcf 2\n!----------------\n\n")
# Next, the size of the box is given in the following format: number of boxes, shape of the box, and the length of one side of the box. 
input_inp.write("# Box_Info\n1\nCUBIC\n100.0\n!---------------\n\n")
# The next line states the temperature at which the simulation will run. 
input_inp.write("# Temperature_Info\n300.0\n!---------------\n\n")
input_inp.write("# Move_Probability_Info\n\n")
input_inp.write("# Prob_Translation\n0.79\n0.5\n\n")
input_inp.write("# Done_Probability_Info\n!----------------\n\n")
# The start type is read_config so it will read in the xyz file in order to determine where the simulation should start
input_inp.write("# Start_Type\nread_config 2 mie.xyz\n!---------------\n\n")
input_inp.write("# Run_Type\nEquilibration 100\n!---------------\n\n")
input_inp.write("# Average_Info\n1\n!---------------(0 = yes, 1 = no)\n\n")
# The next is the line for the simulation length, how long the simulation will run. 
input_inp.write("# Simulation_Length_Info\nUnits Steps\nProp_Freq 1\nCoord_Freq 1\nRun 0\n!---------------\n\n")
# Property info, each line will be a column in the property (prp1) output file.
input_inp.write("# Property_Info 1\nEnergy_Total\nPressure\nVolume\nNmols\nDensity\n!---------------\n\n")
input_inp.write("# Fragment_Files\n!---------------one line per fragment\n\n")
input_inp.write("# CBMC_Info\nkappa_ins 12\nkappa_rot 0\nkappa_dih 10\nrcut_cbmc 6.5\n!---------------\n\n")
input_inp.write("END")
# Close the input file so that Cassandra doesn't explode
input_inp.close()

# Read in input files
inp = open("mie.inp").read() #This command reads in the input file
mcf = open("mie.mcf").read() #Reads in the mcf file
xyz = open("mie.xyz").read() #Reads in the orginal xyz file (used for read_config)

# The following line will print the input file when the hashtag is removed from in front of it.
#print str(inp)

# The following will run Cassandra from this python scripts. There will be four separate tests run (subtests) in order to check the Mie Starting Energies. 

# Run Cassandra Jobs 
# The first test will be two molecules at a distance sigma apart. 
# Allows Cassandra run through Python
proc = sp.Popen(["/afs/crc.nd.edu/x86_64_linux/c/cassandra/src/Cassandra_V1.2/Src/cassandra.exe " + "mie.inp"], stdout=sp.PIPE, shell=True)
(out, err) = proc.communicate()

if err is not None:
	print("Error.Abort. ")


# Function defined to replace a specific line of text with three inputs, 1) the name of the file where you would like to replace the line 2) the line number which you would like to replace (in python indexing starts at 0). 3) the text which you would like to have replace the old text.  
def replace_line(file_name, line_num, text):
	lines = open(file_name, 'r').readlines()
	lines[line_num] = text
	out = open(file_name, 'w')
	out.writelines(lines)
	out.close()

# This will set up the second subtest, where there are once again two molecules. Each molecule is positioned 1/2*sigma away from opposite walls of the box. 
# Changing the xyz input file 
# Change position of atom 1
replace_line('mie.xyz', 2, 'Mie -49.5  0.0  0.0 \n') 
# Change position of atom 2
replace_line('mie.xyz', 3, 'Mie  49.5  0.0  0.0 \n')

# Changing the input file:
# Changes the output name so we can run cassandra under a different name (and save those files too!)
replace_line('mie.inp', 1, 'test2_check2.out\n') 

# Now, we will run cassandra again, under the new name. So, we will get new output files and then hopefully we will be able to compare the two starting energies and they will match!! (If not, we fail and have to start again)
# Run Cassandra Jobs - Again!
# Allows Cassandra run through Python
proc = sp.Popen(["/afs/crc.nd.edu/x86_64_linux/c/cassandra/src/Cassandra_V1.2/Src/cassandra.exe " + "mie.inp"], stdout=sp.PIPE, shell=True)
(out, err) = proc.communicate()

if err is not None:
	print("Error.Abort. ")

# The function replace_line will be used to set up the third test. This test will find the energy of two molecules at a distance of 2^(1/6)*sigma. 
# Change the xyz file
# Change position of atom 1
replace_line('mie.xyz', 2, 'Mie  0.0       0.0  0.0 \n') 
# Change position of atom 2
replace_line('mie.xyz', 3, 'Mie  1.122462  0.0  0.0 \n')

# Changing the input file:
# Changes the output name so we can run cassandra under a different name (and save those files too!)
replace_line('mie.inp', 1, 'test2_check3.out\n') 

# Run Cassandra Jobs - Again!
# Allows Cassandra run through Python
proc = sp.Popen(["/afs/crc.nd.edu/x86_64_linux/c/cassandra/src/Cassandra_V1.2/Src/cassandra.exe " + "mie.inp"], stdout=sp.PIPE, shell=True)
(out, err) = proc.communicate()

if err is not None:
	print("Error.Abort. ")

# The 4th test is the 3rd scenario given to me, out of bounds which should return an energy of affectively zero. In this case, 5 away from the first molecule.
# Changing the xyz input file 
# Change position of atom 1
replace_line('mie.xyz', 2, 'Mie 0.0  0.0  0.0 \n') 
# Change position of atom 2
replace_line('mie.xyz', 3, 'Mie 5.0  0.0  0.0 \n')

# Changing the input file:
# Changes the output name so we can run cassandra under a different name (and save those files too!)
replace_line('mie.inp', 1, 'test2_check4.out\n') 
# For this edit we need to change the VDW_Style section as well. We will change it from a cut_tail 12.0 to a cut 2.5 since the position is 5.0 for the second atom. 
replace_line('mie.inp', 13, 'Mie cut 2.5\n')

# The last of the four subtests is testing a distance outside the bounds. In this case, at a distance 5 away.
# Run Cassandra Jobs - Again!
#print 'Running x4' 
# Allows Cassandra run through Python
proc = sp.Popen(["/afs/crc.nd.edu/x86_64_linux/c/cassandra/src/Cassandra_V1.2/Src/cassandra.exe " + "mie.inp"], stdout=sp.PIPE, shell=True)
(out, err) = proc.communicate()
#print 'done 4'

if err is not None:
	print("Error.Abort. ")

# And now we will find the 4 numbers we need and save them! YIPEE!

# Find a specific line of text and save it. This appears to be impossible. Just kidding. Evidently not impossible 
# We will do it for each of the four tests
# The following 4 lines of the code are designed to extract the line of code from the output file, which contain the starting energy for comparison.
shakes = open('test2_check1.out.log', 'r')

for line in shakes:
	if re.match("(.*)Total system energy(.*)",line):
		line1 = line

# For the second test
shakes = open('test2_check2.out.log', 'r')

for line in shakes:
	if re.match("(.*)Total system energy(.*)",line):
		line2 = line		

# For the third test
shakes = open('test2_check3.out.log', 'r')

for line in shakes:
	if re.match("(.*)Total system energy(.*)",line):
		line3 = line

# For the fourth test
shakes = open('test2_check4.out.log', 'r')

for line in shakes:
	if re.match("(.*)Total system energy(.*)",line):
		line4 = line

# Now we will extract and compare the number in order to Pass Test 2 (each of the 4 subtests will be referred to as a check)
# Now, we will etract the number from the line of string that we saved above. 
num1 = []
for t in line1.split():
	try:
		num1.append(float(t))
	except ValueError:
		pass
num1 = num1[0]

# Now we will extract number 2
num2 = []
for t in line2.split():
	try:
		num2.append(float(t))
	except ValueError:
		pass
num2 = num2[0]

# Extract number 3
num3 = []
for t in line3.split():
	try:
		num3.append(float(t))
	except ValueError:
		pass
num3 = num3[0]

# Extract number 4
num4 = []
for t in line4.split():
	try:
		num4.append(float(t))
	except ValueError:
		pass
num4 = num4[0]

# The following section of the code is the lines and extracetd numbers in a commented out form, this is for ease of review if an error occurs. Simply uncomment them to see the numbers generated by running the Cassandra script. 
#print line1
#print line2
#print line3
#print line4
# Now, the numbers will be printed as well 
#print num1
#print num2
#print num3
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
if num3 == -99.636:
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
	print bold + "Pass Test 2: Mie Starting Energy" + normal
else: 
	print bold + "Test Fails - Check above." + normal
