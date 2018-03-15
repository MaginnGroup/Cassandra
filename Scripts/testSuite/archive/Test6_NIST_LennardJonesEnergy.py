#Test 6: This test tests the LJ energy from the NIST website. There are two checks in the test. The 1st check tests with a cut_tail of 3*sigma. While the second check tests a cut of 4*sigma. Each of the two checks goes in more detail checking the total, vdw, and long range energies. It tests these 3 energies for the 4 different configurations from the NIST website which were downloaded and then edited for use in the Cassandra program.

#Note: The test uses a search to extract the energies from the log file created when running Cassandra based on the name of the energy. Furthermore, the extracted energies where than compared to the energies acquired from the NIST website. 


# Import modules:
import subprocess as sp 
import numpy as np
import random
import os
import re
import sys

# this test does something. 
#We will now create an input file (.inp file)
bold = '\033[1m' #Will make text bold
normal = '\033[0m' #Will make the next text normal(ie. unbold)

#This prints the starting line.
print "\n\n"+bold+"Test 6: NIST Lennard Jones Commencing ... " + normal 
#Creates(opens if already existent) the file file.inp which will be our input file while running cassandra
# The input files are all angle. because it deals with the angle test.
input_inp = open("nist.inp","w")
# Read in mcf and xyz files to use in inp generation
input_mcf = open("nist.mcf","w") #Creates mcf file and allows for edits
input_nist1 = open("nist1.xyz","w") #Creates an input xyz file (used for read_config)

# Define a variable for sigma
sigma = 3.54

#Write one line into xyz
input_nist1.write("800\n\n") # This is the number of atoms in the simulation
# Modify the numbers in lj_sample_config_periodic1.xyz (well an attempt to anyway)
with open("lj_sample_config_periodic1.xyz") as f: 
	for line_nbr, line in enumerate(f):
		if (line_nbr == 0):
			numbers_str = line.split()
			numbers_float = [float(x) for x in numbers_str]
		elif (line_nbr > 1):
			numbers_str = line.split()[1:len(line.split())]
			numbers_float = [float(x) for x in numbers_str]
			numbers_float = [num*sigma for num in numbers_float]
			columns = numbers_float
			print >> input_nist1, "LJ           " + (str(columns[0]) +  "          " + str(columns[1]) + "          " + str(columns[2])) 

input_nist1.close()

#Write info for the MCF file. 
input_mcf.write("# Atom_Info\n1\n1   LJ   LJ   1.000   0.0   LJ   120.2722   3.540\n\n")
input_mcf.write("# Bond_Info\n0\n\n")
input_mcf.write("# Angle_Info\n0\n\n")
input_mcf.write("# Dihedral_Info\n0\n\n")
input_mcf.write("# Improper_Info\n0\n\n")
input_mcf.write("# Fragment_Info\n0\n\n")
input_mcf.write("# Fragment_Connectivity\n0\n\n")
input_mcf.write("# Intra_Scaling\n0. 0. 0. 0.\n0. 0. 0. 0.\n\n\n")
input_mcf.write("END")
input_mcf.close()


# Write info into the file - this will create each section for the .inp file
# This input file is populated with numbers for an Mie simulation with an arbitrary atom
# This input decides what the output files will be named, it should be unique so as not to overwrite existing files. 
input_inp.write("# Run_Name\ntest6_check1_nist1.out\n!---------------\n\n")
# This line decides the simulation type. In this case, we are running an NVT Monte Carlo, that is constant mols, volume, and temperature.
input_inp.write("# Sim_Type\nNVT_MC\n!---------------\n\n")
# The number of species in the simulation
input_inp.write("# Nbr_Species\n1\n!---------------\n\n")
# This determines where to cut off the very inaccurate part of the spectrum
input_inp.write("# VDW_Style\nLJ cut_tail 10.620\n!---------------\n\n")
# There is no charge, so NONE is inputted 
input_inp.write("# Charge_Style\nNONE\n!---------------\n\n")
input_inp.write("# Intra_Scaling\n0.0 0.0 0.0 1.0\n!---------------vdw, coul lines for each species\n\n")
input_inp.write("# Mixing_Rule\nLB\n!---------------\n\n")
# Seed info can be the generation of any random number 
input_inp.write("# Seed_Info\n55001389 1764321284\n!---------------\n\n")
input_inp.write("# Rcutoff_Low\n0.0\n!---------------\n\n")
# This line tells us from which mcf file we will be reading, and how many molecules should be read in.
input_inp.write("# Molecule_Files\nnist.mcf 800\n!----------------\n\n")
# Next, the size of the box is given in the following format: number of boxes, shape of the box, and the length of one side of the box. 
input_inp.write("# Box_Info\n1\nCUBIC\n35.4\n!---------------\n\n")
# The next line states the temperature at which the simulation will run. 
input_inp.write("# Temperature_Info\n300.0\n!---------------\n\n")
input_inp.write("# Move_Probability_Info\n\n")
input_inp.write("# Prob_Translation\n0.79\n0.5\n\n")
input_inp.write("# Done_Probability_Info\n!----------------\n\n")
# The start type is read_config so it will read in the xyz file in order to determine where the simulation should start
input_inp.write("# Start_Type\nread_config 800 nist1.xyz\n!---------------\n\n")
input_inp.write("# Run_Type\nEquilibration 100\n!---------------\n\n")
input_inp.write("# Average_Info\n1\n!---------------(0 = yes, 1 = no)\n\n")
# The next is the line for the simulation length, how long the simulation will run. 
input_inp.write("# Simulation_Length_Info\nUnits Steps\nProp_Freq 1\nCoord_Freq 1\nRun 1\n!---------------\n\n")
# Property info, each line will be a column in the property (prp1) output file.
input_inp.write("# Property_Info 1\nEnergy_Total\nPressure\nVolume\nNmols\nDensity\n!---------------\n\n")
input_inp.write("# Fragment_Files\n!---------------one line per fragment\n\n")
input_inp.write("# CBMC_Info\nkappa_ins 12\nkappa_rot 0\nkappa_dih 10\nrcut_cbmc 6.5\n!---------------\n\n")
input_inp.write("END")
# Close the input file so that Cassandra doesn't explode
input_inp.close()

#open files for read
inp = open("nist.inp").read() #This command reads in the input file
mcf = open("nist.mcf").read() #Reads in the mcf file
xyz = open("nist1.xyz").read() #Reads in the orginal xyz file (used for read_config)


# Run Cassandra
#proc = sp.Popen(["/afs/crc.nd.edu/x86_64_linux/c/cassandra/src/Cassandra_V1.2/Src/cassandra.exe " + "nist.inp"], stdout=sp.PIPE, shell=True) 
proc = sp.Popen(["/afs/crc.nd.edu/user/l/lpink/Git/Cassandra/Src/cassandra.exe " + "nist.inp"], stdout=sp.PIPE, shell=True) 
(out, err) = proc.communicate()

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


# NIST 2
input_nist2 = open("nist2.xyz","w") #Creates an input xyz file (used for read_config)
#
# Define a variable for sigma
sigma = 3.54

#Write one line into xyz
input_nist2.write("200\n\n") # This is the number of atoms in the simulation
# Modify the numbers in lj_sample_config_periodic1.xyz (well an attempt to anyway)
with open("lj_sample_config_periodic2.xyz") as f: 
	for line_nbr, line in enumerate(f):
		if (line_nbr == 0):
			numbers_str = line.split()
			numbers_float = [float(x) for x in numbers_str]
		elif (line_nbr > 1):
			numbers_str = line.split()[1:len(line.split())]
			numbers_float = [float(x) for x in numbers_str]
			numbers_float = [num*sigma for num in numbers_float]
			columns = numbers_float
			print >> input_nist2, "LJ           " + (str(columns[0]) +  "          " + str(columns[1]) + "          " + str(columns[2])) 

input_nist2.close()

# We will now test the third set of data from the nist website
# Now, we will replace some lines for fun. 
replace_line("nist.inp", 1, "test6_check1_nist2.out\n")
replace_line("nist.inp", 60, "read_config 200 nist2.xyz\n")
replace_line("nist.inp", 37, "nist.mcf 200\n")
replace_line("nist.inp", 43, "28.32\n")

# open files for read
xyz = open("nist2.xyz").read() #Reads in the orginal xyz file (used for read_config)


# Run Cassandra
proc = sp.Popen(["/afs/crc.nd.edu/user/l/lpink/Git/Cassandra/Src/cassandra.exe " + "nist.inp"], stdout=sp.PIPE, shell=True) 
(out, err) = proc.communicate()

if err is not None: 
	print("Error.Abort. ")

# NIST 3
input_nist3 = open("nist3.xyz","w") #Creates an input xyz file (used for read_config)

# Define a variable for sigma
sigma = 3.54

#Write one line into xyz
input_nist3.write("400\n\n") # This is the number of atoms in the simulation
# Modify the numbers in lj_sample_config_periodic1.xyz (well an attempt to anyway)
with open("lj_sample_config_periodic3.xyz") as f: 
	for line_nbr, line in enumerate(f):
		if (line_nbr == 0):
			numbers_str = line.split()
			numbers_float = [float(x) for x in numbers_str]
		elif (line_nbr > 1):
			numbers_str = line.split()[1:len(line.split())]
			numbers_float = [float(x) for x in numbers_str]
			numbers_float = [num*sigma for num in numbers_float]
			columns = numbers_float
			print >> input_nist3, "LJ           " + (str(columns[0]) +  "          " + str(columns[1]) + "          " + str(columns[2])) 

input_nist3.close()

# We will now test the third set of data from the nist website
# Now, we will replace some lines for fun. 
replace_line("nist.inp", 1, "test6_check1_nist3.out\n")
replace_line("nist.inp", 60, "read_config 400 nist3.xyz\n")
replace_line("nist.inp", 37, "nist.mcf 400\n")
replace_line("nist.inp", 43, "35.4\n")



# open files for read
xyz = open("nist3.xyz").read() #Reads in the orginal xyz file (used for read_config)


# Run Cassandra
#proc = sp.Popen(["/afs/crc.nd.edu/x86_64_linux/c/cassandra/src/Cassandra_V1.2/Src/cassandra.exe " + "nist.inp"], stdout=sp.PIPE, shell=True) 
proc = sp.Popen(["/afs/crc.nd.edu/user/l/lpink/Git/Cassandra/Src/cassandra.exe " + "nist.inp"], stdout=sp.PIPE, shell=True) 
(out, err) = proc.communicate()

if err is not None: 
	print("Error.Abort. ")

# NOW FOR NIST 4
input_nist4 = open("nist4.xyz","w") #Creates an input xyz file (used for read_config)

# Define a variable for sigma
sigma = 3.54

#Write one line into xyz
input_nist4.write("30\n\n") # This is the number of atoms in the simulation
# Modify the numbers in lj_sample_config_periodic1.xyz (well an attempt to anyway)
with open("lj_sample_config_periodic4.xyz") as f: 
	for line_nbr, line in enumerate(f):
		if (line_nbr == 0):
			numbers_str = line.split()
			numbers_float = [float(x) for x in numbers_str]
		elif (line_nbr > 1):
			numbers_str = line.split()[1:len(line.split())]
			numbers_float = [float(x) for x in numbers_str]
			numbers_float = [num*sigma for num in numbers_float]
			columns = numbers_float
			print >> input_nist4, "LJ           " + (str(columns[0]) +  "          " + str(columns[1]) + "          " + str(columns[2])) 

input_nist4.close()

# now for nist 4
# Now, we will replace some lines for fun. 
replace_line("nist.inp", 1, "test6_check1_nist4.out\n")
replace_line("nist.inp", 60, "read_config 30 nist4.xyz\n")
replace_line("nist.inp", 37, "nist.mcf 30\n")
replace_line("nist.inp", 43, "28.32\n")


# open files for read
xyz = open("nist4.xyz").read() #Reads in the orginal xyz file (used for read_config)


# Run Cassandra
##proc = sp.Popen(["/afs/crc.nd.edu/x86_64_linux/c/cassandra/src/Cassandra_V1.2/Src/cassandra.exe " + "nist.inp"], stdout=sp.PIPE, shell=True) 
proc = sp.Popen(["/afs/crc.nd.edu/user/l/lpink/Git/Cassandra/Src/cassandra.exe " + "nist.inp"], stdout=sp.PIPE, shell=True) 
(out, err) = proc.communicate()

if err is not None: 
	print("Error.Abort. ")



# CHECK 2 - This is for a cut_tail with 4 sigma and again using all 4 nist information
# Now, we will replace some lines for fun. 
replace_line("nist.inp", 1, "test6_check2_nist1.out\n")
replace_line("nist.inp", 13, "LJ cut_tail 14.16\n")
replace_line("nist.inp", 60, "read_config 800 nist1.xyz\n")
replace_line("nist.inp", 37, "nist.mcf 800\n")
replace_line("nist.inp", 43, "35.4\n")




# Run Cassandra
#proc = sp.Popen(["/afs/crc.nd.edu/x86_64_linux/c/cassandra/src/Cassandra_V1.2/Src/cassandra.exe " + "nist.inp"], stdout=sp.PIPE, shell=True) 
proc = sp.Popen(["/afs/crc.nd.edu/user/l/lpink/Git/Cassandra/Src/cassandra.exe " + "nist.inp"], stdout=sp.PIPE, shell=True) 
(out, err) = proc.communicate()

if err is not None: 
	print("Error.Abort. ")

# Check 2 - Nist 2. 
replace_line("nist.inp", 1, "test6_check2_nist2.out\n")
replace_line("nist.inp", 60, "read_config 200 nist2.xyz\n")
replace_line("nist.inp", 37, "nist.mcf 200\n")
replace_line("nist.inp", 43, "28.32\n")


# Run Cassandra
##proc = sp.Popen(["/afs/crc.nd.edu/x86_64_linux/c/cassandra/src/Cassandra_V1.2/Src/cassandra.exe " + "nist.inp"], stdout=sp.PIPE, shell=True) 
proc = sp.Popen(["/afs/crc.nd.edu/user/l/lpink/Git/Cassandra/Src/cassandra.exe " + "nist.inp"], stdout=sp.PIPE, shell=True) 
(out, err) = proc.communicate()

if err is not None: 
	print("Error.Abort. ")

# Check 2 - Nist 3. 
replace_line("nist.inp", 1, "test6_check2_nist3.out\n")
replace_line("nist.inp", 60, "read_config 400 nist3.xyz\n")
replace_line("nist.inp", 37, "nist.mcf 400\n")
replace_line("nist.inp", 43, "35.4\n")


# Run Cassandra
#proc = sp.Popen(["/afs/crc.nd.edu/x86_64_linux/c/cassandra/src/Cassandra_V1.2/Src/cassandra.exe " + "nist.inp"], stdout=sp.PIPE, shell=True) 
proc = sp.Popen(["/afs/crc.nd.edu/user/l/lpink/Git/Cassandra/Src/cassandra.exe " + "nist.inp"], stdout=sp.PIPE, shell=True) 
(out, err) = proc.communicate()

if err is not None: 
	print("Error.Abort. ")

# Check 2 - Nist 4. 
replace_line("nist.inp", 1, "test6_check2_nist4.out\n")
replace_line("nist.inp", 60, "read_config 30 nist4.xyz\n")
replace_line("nist.inp", 37, "nist.mcf 30\n")
replace_line("nist.inp", 43, "28.32\n")


## Run Cassandra
##proc = sp.Popen(["/afs/crc.nd.edu/x86_64_linux/c/cassandra/src/Cassandra_V1.2/Src/cassandra.exe " + "nist.inp"], stdout=sp.PIPE, shell=True) 
proc = sp.Popen(["/afs/crc.nd.edu/user/l/lpink/Git/Cassandra/Src/cassandra.exe " + "nist.inp"], stdout=sp.PIPE, shell=True) 
(out, err) = proc.communicate()

if err is not None: 
	print("Error.Abort. ")

# Next, the four starting energies will be found in the output .log file for each test. This line will be extracted and saved in this python script so we can use it for comparison at the very end of the script. 
# Finding a string - using an if statement in a for loop
# For the first test
# shakes opens the desired file in the read format
shakes = open("test6_check1_nist1.out.log", "r")

# Now we will extract numbers --> CHECK 1 NIST 1
# The for loop will search line by line in shakes for the words "Total system energy", once found the line will be saved as a variable.
for line in shakes:
	if re.match("(.*)Total system energy(.*)",line):
		line_tot11 = line

# Extract vdw energy
shakes = open("test6_check1_nist1.out.log", "r")
for line in shakes:
	if re.match("(.*)Inter molecule vdw(.*)",line):
		line_vdw11 = line

#Extract long range correction energy
shakes = open("test6_check1_nist1.out.log", "r")
for line in shakes:
	if re.match("(.*)Long range correction(.*)",line):
		line_long11 = line


#print "1"
#print line_tot11
#print line_vdw11
#print line_long11

# CHECK 1 NIST 2
shakes = open("test6_check1_nist2.out.log", "r")
for line in shakes:
	if re.match("(.*)Total system energy(.*)",line):
		line_tot12 = line

# Extract vdw energy
shakes = open("test6_check1_nist2.out.log", "r")
for line in shakes:
	if re.match("(.*)Inter molecule vdw(.*)",line):
		line_vdw12 = line

#Extract long range correction energy
shakes = open("test6_check1_nist2.out.log", "r")
for line in shakes:
	if re.match("(.*)Long range correction(.*)",line):
		line_long12 = line


#print "2"
#print line_tot12
#print line_vdw12
#print line_long12

# CHECK 1 NIST 3
shakes = open("test6_check1_nist3.out.log", "r")
for line in shakes:
	if re.match("(.*)Total system energy(.*)",line):
		line_tot13 = line

# Extract vdw energy
shakes = open("test6_check1_nist3.out.log", "r")
for line in shakes:
	if re.match("(.*)Inter molecule vdw(.*)",line):
		line_vdw13 = line

#Extract long range correction energy
shakes = open("test6_check1_nist3.out.log", "r")
for line in shakes:
	if re.match("(.*)Long range correction(.*)",line):
		line_long13 = line


#print "3"
#print line_tot13
#print line_vdw13
#print line_long13

# CHECK 1 NIST 4
shakes = open("test6_check1_nist4.out.log", "r")
for line in shakes:
	if re.match("(.*)Total system energy(.*)",line):
		line_tot14 = line

# Extract vdw energy
shakes = open("test6_check1_nist4.out.log", "r")
for line in shakes:
	if re.match("(.*)Inter molecule vdw(.*)",line):
		line_vdw14 = line

#Extract long range correction energy
shakes = open("test6_check1_nist4.out.log", "r")
for line in shakes:
	if re.match("(.*)Long range correction(.*)",line):
		line_long14 = line


#print "4"
#print line_tot14
#print line_vdw14
#print line_long14

# CHECK 2 NIST 1
shakes = open("test6_check2_nist1.out.log", "r")
for line in shakes:
	if re.match("(.*)Total system energy(.*)",line):
		line_tot21 = line

# Extract vdw energy
shakes = open("test6_check2_nist1.out.log", "r")
for line in shakes:
	if re.match("(.*)Inter molecule vdw(.*)",line):
		line_vdw21 = line

#Extract long range correction energy
shakes = open("test6_check2_nist1.out.log", "r")
for line in shakes:
	if re.match("(.*)Long range correction(.*)",line):
		line_long21 = line


#print "5"
#print line_tot21
#print line_vdw21
#print line_long21

# CHECK 2 NIST 2
shakes = open("test6_check2_nist2.out.log", "r")
for line in shakes:
	if re.match("(.*)Total system energy(.*)",line):
		line_tot22 = line

# Extract vdw energy
shakes = open("test6_check2_nist2.out.log", "r")
for line in shakes:
	if re.match("(.*)Inter molecule vdw(.*)",line):
		line_vdw22 = line

#Extract long range correction energy
shakes = open("test6_check2_nist2.out.log", "r")
for line in shakes:
	if re.match("(.*)Long range correction(.*)",line):
		line_long22 = line


#print "6"
#print line_tot22
#print line_vdw22
#print line_long22

# CHECK 2 NIST 3
shakes = open("test6_check2_nist3.out.log", "r")
for line in shakes:
	if re.match("(.*)Total system energy(.*)",line):
		line_tot23 = line

# Extract vdw energy
shakes = open("test6_check2_nist3.out.log", "r")
for line in shakes:
	if re.match("(.*)Inter molecule vdw(.*)",line):
		line_vdw23 = line

#Extract long range correction energy
shakes = open("test6_check2_nist3.out.log", "r")
for line in shakes:
	if re.match("(.*)Long range correction(.*)",line):
		line_long23 = line


#print "7"
#print line_tot23
#print line_vdw23
#print line_long23

# CHECK 2 NIST 4
shakes = open("test6_check2_nist4.out.log", "r")
for line in shakes:
	if re.match("(.*)Total system energy(.*)",line):
		line_tot24 = line

# Extract vdw energy
shakes = open("test6_check2_nist4.out.log", "r")
for line in shakes:
	if re.match("(.*)Inter molecule vdw(.*)",line):
		line_vdw24 = line

#Extract long range correction energy
shakes = open("test6_check2_nist4.out.log", "r")
for line in shakes:
	if re.match("(.*)Long range correction(.*)",line):
		line_long24 = line


#print "8"
#print line_tot24
#print line_vdw24
#print line_long24

# NOW WE will extract the number!
# For check 1 nist 1 total
num111 = []
# Use a for loop in order to go through each character in the line independently.
for t in line_tot11.split():
	try:
		num111.append(float(t))
	except ValueError:
		pass
num111 = num111[0]

# for check 1 nist 1 vdw
num112 = []
for t in line_vdw11.split():
	try:
		num112.append(float(t))
	except ValueError:
		pass
num112 = num112[0]

# for check 1 nist 1 long
num113 = []
for t in line_long11.split():
	try:
		num113.append(float(t))
	except ValueError:
		pass
num113 = num113[0]

#print num111
#print num112
#print num113

# For check 1 nist 2 total
num121 = []
# Use a for loop in order to go through each character in the line independently.
for t in line_tot12.split():
# Save the number in the line to the variable num1 as a float
	try:
		num121.append(float(t))
# Otherwise, do nothing.
	except ValueError:
		pass
# Extract the only number in the list as a variable.
num121 = num121[0]

# for check 1 nist 2 vdw
num122 = []
# Use a for loop in order to go through each character in the line independently.
for t in line_vdw12.split():
# Save the number in the line to the variable num1 as a float
	try:
		num122.append(float(t))
# Otherwise, do nothing.
	except ValueError:
		pass
# Extract the only number in the list as a variable.
num122 = num122[0]

# for check 1 nist 2 long
num123 = []
# Use a for loop in order to go through each character in the line independently.
for t in line_long12.split():
# Save the number in the line to the variable num1 as a float
	try:
		num123.append(float(t))
# Otherwise, do nothing.
	except ValueError:
		pass
# Extract the only number in the list as a variable.
num123 = num123[0]

#print num121
#print num122
#print num123

# For check 1 nist 3 total
num131 = []
# Use a for loop in order to go through each character in the line independently.
for t in line_tot13.split():
	try:
		num131.append(float(t))
	except ValueError:
		pass
num131 = num131[0]

# for check 1 nist 3 vdw
num132 = []
for t in line_vdw13.split():
	try:
		num132.append(float(t))
	except ValueError:
		pass
num132 = num132[0]

# for check 1 nist 3 long
num133 = []
for t in line_long13.split():
	try:
		num133.append(float(t))
	except ValueError:
		pass
num133 = num133[0]

#print num131
#print num132
#print num133

# For check 1 nist 4 total
num141 = []
# Use a for loop in order to go through each character in the line independently.
for t in line_tot14.split():
	try:
		num141.append(float(t))
	except ValueError:
		pass
num141 = num141[0]

# for check 1 nist 4 vdw
num142 = []
for t in line_vdw14.split():
	try:
		num142.append(float(t))
	except ValueError:
		pass
num142 = num142[0]

# for check 1 nist 4 long
num143 = []
for t in line_long14.split():
	try:
		num143.append(float(t))
	except ValueError:
		pass
num143 = num143[0]

#print num141
#print num142
#print num143

# For check 2 nist 1 total
num211 = []
# Use a for loop in order to go through each character in the line independently.
for t in line_tot21.split():
	try:
		num211.append(float(t))
	except ValueError:
		pass
num211 = num211[0]

# for check 1 nist 1 vdw
num212 = []
for t in line_vdw21.split():
	try:
		num212.append(float(t))
	except ValueError:
		pass
num212 = num212[0]

# for check 1 nist 1 long
num213 = []
for t in line_long21.split():
	try:
		num213.append(float(t))
	except ValueError:
		pass
num213 = num213[0]

#print num211
#print num212
#print num213

# For check 1 nist 2 total
num221 = []
# Use a for loop in order to go through each character in the line independently.
for t in line_tot22.split():
# Save the number in the line to the variable num1 as a float
	try:
		num221.append(float(t))
# Otherwise, do nothing.
	except ValueError:
		pass
# Extract the only number in the list as a variable.
num221 = num221[0]

# for check 1 nist 2 vdw
num222 = []
# Use a for loop in order to go through each character in the line independently.
for t in line_vdw22.split():
# Save the number in the line to the variable num1 as a float
	try:
		num222.append(float(t))
# Otherwise, do nothing.
	except ValueError:
		pass
# Extract the only number in the list as a variable.
num222 = num222[0]

# for check 1 nist 2 long
num223 = []
# Use a for loop in order to go through each character in the line independently.
for t in line_long22.split():
# Save the number in the line to the variable num1 as a float
	try:
		num223.append(float(t))
# Otherwise, do nothing.
	except ValueError:
		pass
# Extract the only number in the list as a variable.
num223 = num223[0]

#print num221
#print num222
#print num223

# For check 1 nist 3 total
num231 = []
# Use a for loop in order to go through each character in the line independently.
for t in line_tot23.split():
	try:
		num231.append(float(t))
	except ValueError:
		pass
num231 = num231[0]

# for check 1 nist 3 vdw
num232 = []
for t in line_vdw23.split():
	try:
		num232.append(float(t))
	except ValueError:
		pass
num232 = num232[0]

# for check 1 nist 3 long
num233 = []
for t in line_long23.split():
	try:
		num233.append(float(t))
	except ValueError:
		pass
num233 = num233[0]

#print num231
#print num232
#print num233

# For check 1 nist 4 total
num241 = []
# Use a for loop in order to go through each character in the line independently.
for t in line_tot24.split():
	try:
		num241.append(float(t))
	except ValueError:
		pass
num241 = num241[0]

# for check 1 nist 4 vdw
num242 = []
for t in line_vdw24.split():
	try:
		num242.append(float(t))
	except ValueError:
		pass
num242 = num242[0]

# for check 1 nist 4 long
num243 = []
for t in line_long24.split():
	try:
		num243.append(float(t))
	except ValueError:
		pass
num243 = num243[0]

#print num241
#print num242
#print num243

# Lastly, we will compare all these numbers and test the Cassandra Code 
# Check 1
print "\n"+bold+"Check 1 (Cutoff Radius, 3 sigma): " + normal
# For NIST 1
print bold+"Configuration NIST 1: " + normal
# For Total
if (num111/100 + 4548.72) < 1:
	c1 = 1
	print "Total..."
else: 
	c1 = 0 
	print "Total Fails."

# For vdw
if (num112/100 + 4351.5) < 4:
	c2 = 1
	print "Vdw..."
else: 
	c2 = 0 
	print "Vdw Fails."

#for long
if (num113/100 + 198.49) < 1:
	c3 = 1
	print "Long range..."
else: 
	c3 = 0 
	print "Long Range Fails."

if c1 == 1 and c2 == 1 and c3 == 1:
	n1 = 1
	print bold + "NIST 1..." + normal
else:
	n1 = 0 
	print bold + "NIST 1 fails. See above." + normal

# For NIST 2
print bold+"Configuration NIST 2: " + normal
# For Total
if (num111/100 + 710.96) < 1:
	c1 = 1
	print "Total..."
else: 
	c1 = 0 
	print "Total Fails."

# For vdw
if (num112/100 + 690.0) < 5:
	c2 = 1
	print "Vdw..."
else: 
	c2 = 0 
	print "Vdw Fails."

#for long
if (num113/100 + 24.230) < 1:
	c3 = 1
	print "Long range..."
else: 
	c3 = 0 
	print "Long Range Fails."

if c1 == 1 and c2 == 1 and c3 == 1:
	n2 = 1
	print bold + "NIST 2..." + normal
else:
	n2 = 0 
	print bold + "NIST 2 fails. See above." + normal




# For NIST 3
print bold+"Configuration NIST 3: " + normal
# For Total
if (num131/100 + -1194.59) < 1:
	c1 = 1
	print "Total..."
else: 
	c1 = 0 
	print "Total Fails."

# For vdw
if (num132/100 + 1144.96 ) < 4:
	c2 = 1
	print "Vdw..."
else: 
	c2 = 0 
	print "Vdw Fails."

#for long
if (num133/100 + 49.622) < 1:
	c3 = 1
	print "Long range..."
else: 
	c3 = 0 
	print "Long Range Fails."

# pass nist 2 - check 1

if c1 == 1 and c2 == 1 and c3 == 1:
	n3 = 1
	print bold + "NIST 3..."+normal
else:
	n3 = 0 
	print bold + "NIST 3 fails. See above." + normal

#NIST 4
print bold+"Configuration NIST 4: " + normal
# For Total
if (num131/100 + -17.46) < 1:
	c1 = 1
	print "Total..."
else: 
	c1 = 0 
	print "Total Fails."

# For vdw
if (num132/100 + 16.92 ) < 1:
	c2 = 1
	print "Vdw..."
else: 
	c2 = 0 
	print "Vdw Fails."

#for long
if (num133/100 + 5.4517) < 1:
	c3 = 1
	print "Long range..."
else: 
	c3 = 0 
	print "Long Range Fails."

# pass nist 2 - check 1

if c1 == 1 and c2 == 1 and c3 == 1:
	n4 = 1
	print bold + "NIST 4..."+normal
else:
	n4 = 0 
	print bold + "NIST 4 fails. See above." + normal



#check all nist
if n1 == 1 and n2 == 1 and n3 == 1 and n4 == 1:
	C1 = 1
	print bold + "Check 1..." + normal
else:
	C1 = 0
	print bold + "check 1 fails. See above." + normal




# CHECK 2
# Lastly, we will compare all these numbers and test the Cassandra Code 
# Check 1
print "\n"+bold+"Check 2 (Cutoff Radius, 4 sigma): " + normal
# For NIST 1
print bold+"Configuration NIST 1: " + normal
# For Total
if (num211/100 + 4549.95) < 4:
	c1 = 1
	print "Total..."
else: 
	c1 = 0 
	print "Total Fails."

# For vdw
if (num212/100 + 4467.5) < 4:
	c2 = 1
	print "Vdw..."
else: 
	c2 = 0 
	print "Vdw Fails."

#for long
if (num213/100 + 83.769) < 1:
	c3 = 1
	print "Long range..."
else: 
	c3 = 0 
	print "Long Range Fails."

if c1 == 1 and c2 == 1 and c3 == 1:
	n1 = 1
	print bold + "NIST 1..." + normal
else:
	n1 = 0 
	print bold + "NIST 1 fails. See above." + normal

# For NIST 2
print bold+"Configuration NIST 2: " + normal
# For Total
if (num111/100 + 711.56) < 1:
	c1 = 1
	print "Total..."
else: 
	c1 = 0 
	print "Total Fails."

# For vdw
if (num112/100 + 704.60) < 5:
	c2 = 1
	print "Vdw..."
else: 
	c2 = 0 
	print "Vdw Fails."

#for long
if (num113/100 + 10.210) < 1:
	c3 = 1
	print "Long range..."
else: 
	c3 = 0 
	print "Long Range Fails."

if c1 == 1 and c2 == 1 and c3 == 1:
	n2 = 1
	print bold + "NIST 2..." + normal
else:
	n2 = 0 
	print bold + "NIST 2 fails. See above." + normal



# For NIST 3
print bold+"Configuration NIST 3: " + normal
# For Total
if (num231/100 + -1194.62) < 1:
	c1 = 1
	print "Total..."
else: 
	c1 = 0 
	print "Total Fails."

# For vdw
if (num232/100 + 1175.4) < 4:
	c2 = 1
	print "Vdw..."
else: 
	c2 = 0 
	print "Vdw Fails."

#for long
if (num233/100 + 20.942) < 1:
	c3 = 1
	print "Long range..."
else: 
	c3 = 0 
	print "Long Range Fails."

# pass nist 2 - check 1

if c1 == 1 and c2 == 1 and c3 == 1:
	n3 = 1
	print bold + "NIST 3..." + normal
else:
	n3 = 0 
	print bold + "NIST 3 fails. See above." + normal

#NIST 4
print bold+"Configuration NIST 4: " + normal
# For Total
if (num131/100 + -17.41) < 1:
	c1 = 1
	print "Total..."
else: 
	c1 = 0 
	print "Total Fails."

# For vdw
if (num132/100 + 17.060 ) < 1:
	c2 = 1
	print "Vdw..."
else: 
	c2 = 0 
	print "Vdw Fails."

#for long
if (num133/100 + 2.3008) < 1:
	c3 = 1
	print "Long range..."
else: 
	c3 = 0 
	print "Long Range Fails."

# pass nist 2 - check 1

if c1 == 1 and c2 == 1 and c3 == 1:
	n4 = 1
	print bold + "NIST 4..."+normal
else:
	n4 = 0 
	print bold + "NIST 4 fails. See above." + normal



#check all nist
if n1 == 1 and n2 == 1 and n3 == 1 and n4 == 1:
	C2 = 1
	print bold + "Check 2..." + normal
else:
	C2 = 0
	print bold + "check 2 fails. See above." + normal



# FINAL CHECK: 

if C1 == 1 and C2 == 1: 
	print bold + "Pass Test 6: NIST Lennard Jones Energy." + normal
else: 
	print bold + "Test 6 fails." + normal
