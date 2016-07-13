# This is Test 7 in a series of tests for a testSuite in order to check updates made to the Cassandra program. 
# Test 7 blah blah blah 
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
print "\n\n"+bold+"Test 7: Charge Energy (Water SPCE) Commencing ... " + normal 
#Creates(opens if already existent) the file file.inp which will be our input file while running cassandra
input_inp = open("water.inp","w")
# Read in mcf and xyz files to use in inp generation
input_mcf = open("water.mcf","w") #Creates mcf file and allows for edits
input_water1 = open("water1.xyz", "w")


# Write text file into an xyz file
#Write one line into xyz
input_water1.write("100\n\n") # This is the number of atoms in the simulation
# Modify the numbers in lj_sample_config_periodic1.xyz (well an attempt to anyway)
with open("spce_sample_config_periodic1.xyz") as f: 
	for line_nbr, line in enumerate(f):
		if (line_nbr == 0):
			numbers_str = line.split()
		elif (line_nbr > 1):
			numbers_str = line.split()[1:len(line.split())]
		#	numbers_float = [float(x) for x in numbers_str]
			numbers_float = [num for num in numbers_str]
			columns = numbers_float
			print >> input_water1, (str(columns[3]) +  "       " + str(columns[0]) +  "          " + str(columns[1]) + "          " + str(columns[2])) 


input_water1.close()


#Write info for the MCF file. 
input_mcf.write("# Atom_Info\n3\n1   O   O   15.999   -0.8476   LJ   78.1974311   3.16555789\n2   H   H   1.008   0.42380   LJ   0.000   0.000\n3   H   H   1.008   0.42380   LJ   0.000   0.000\n\n")
input_mcf.write("# Bond_Info\n2\n1   1   2   fixed   1.000\n2   1   3   fixed   1.000\n\n")
input_mcf.write("# Angle_Info\n1\n1   2   1   3   fixed   109.47\n\n")
input_mcf.write("# Dihedral_Info\n0\n\n")
input_mcf.write("# Improper_Info\n0\n\n")
input_mcf.write("# Intra_Scaling\n0.   0.   0.0000   1.\n0.   0.   0.0000   1.\n\n")
input_mcf.write("# Fragment_Info\n0\n\n")
input_mcf.write("# Fragment_Connectivity\n0\n\n")
input_mcf.write("END")
input_mcf.close()

# Write info into the file - this will create each section for the .inp file
# This input file is populated with numbers for an Mie simulation with an arbitrary atom
# This input decides what the output files will be named, it should be unique so as not to overwrite existing files. 
input_inp.write("# Run_Name\ntest7_water1.out\n!---------------\n\n")
# This line decides the simulation type. In this case, we are running an NVT Monte Carlo, that is constant mols, volume, and temperature.
input_inp.write("# Sim_Type\nNVT_MC\n!---------------\n\n")
# The number of species in the simulation
input_inp.write("# Nbr_Species\n1\n!---------------\n\n")
# This determines where to cut off the very inaccurate part of the spectrum
input_inp.write("# VDW_Style\nLJ cut_tail 10.000\n!---------------\n\n")
# There is no charge, so NONE is inputted 
input_inp.write("# Charge_Style\ncoul ewald 10.0 0.000001\n!---------------\n\n")
#input_inp.write("# Intra_Scaling\n0.0 0.0 0.0 1.0\n!---------------vdw, coul lines for each species\n\n")
input_inp.write("# Mixing_Rule\nLB\n!---------------\n\n")
# Seed info can be the generation of any random number 
input_inp.write("# Seed_Info\n21498  489625\n!---------------\n\n")
input_inp.write("# Rcutoff_Low\n0.85\n!---------------\n\n")
# This line tells us from which mcf file we will be reading, and how many molecules should be read in.
input_inp.write("# Molecule_Files\nwater.mcf 100\n!----------------\n\n")
# Next, the size of the box is given in the following format: number of boxes, shape of the box, and the length of one side of the box. 
input_inp.write("# Box_Info\n1\nCUBIC\n20.0\n!---------------\n\n")
# The next line states the temperature at which the simulation will run. 
input_inp.write("# Temperature_Info\n300.0\n!---------------\n\n")
input_inp.write("# Move_Probability_Info\n\n")
input_inp.write("# Prob_Translation\n90\n38.0\n\n")
input_inp.write("# Done_Probability_Info\n!----------------\n\n")
# The start type is read_config so it will read in the xyz file in order to determine where the simulation should start
input_inp.write("# Start_Type\nread_config 100  water1.xyz\n!---------------\n\n")
input_inp.write("# Run_Type\nEquilibration 100\n!---------------\n\n")
input_inp.write("# Average_Info\n1\n!---------------(0 = yes, 1 = no)\n\n")
# The next is the line for the simulation length, how long the simulation will run. 
input_inp.write("# Simulation_Length_Info\nUnits Steps\nProp_Freq 1\nCoord_Freq 1\nRun 1\n!---------------\n\n")
# Property info, each line will be a column in the property (prp1) output file.
input_inp.write("# Property_Info 1\nEnergy_Total\nPressure\nVolume\nNmols\nDensity\n!---------------\n\n")
input_inp.write("# Fragment_Files\n!---------------one line per fragment\n\n")
input_inp.write("# CBMC_Info\nkappa_ins 12\nrcut_cbmc 6.5\n!---------------\n\n")
input_inp.write("END")
# Close the input file so that Cassandra doesn't explode
input_inp.close()


#open files for read
inp = open("water.inp").read() #This command reads in the input file
mcf = open("water.mcf").read() #Reads in the mcf file
xyz = open("water1.xyz").read() #Reads in the orginal xyz file (used for read_config)


# Run Cassandra
#print "Runnning"
#proc = sp.Popen(["/afs/crc.nd.edu/x86_64_linux/c/cassandra/src/Cassandra_V1.2/Src/cassandra_intel_openMP.exe " + "nist.inp"], stdout=sp.PIPE, shell=True) 
proc = sp.Popen(["/afs/crc.nd.edu/user/l/lpink/Git/Cassandra/Src/cassandra_intel_openMP.exe " + "water.inp"], stdout=sp.PIPE, shell=True) 
(out, err) = proc.communicate()

if err is not None: 
	print("Error.Abort. ")

#print "done"

# Next, we will write over specific lines of the input file created above using a function, called replace_line that is created below. 
# This fucntion takes three inputs: the name of the file where you would like to replace a line, the line number you would like to replace, and the text you want to replace the old text with. 
def replace_line(file_name, line_num, text):
	lines = open(file_name, 'r').readlines()
	lines[line_num] = text
	out = open(file_name, 'w')
	out.writelines(lines)
	out.close() # Closes the file so that the program doesn't explode. 

# water 2
input_water2 = open("water2.xyz","w") #creates an input xyz file (used for read_config)
# write a new xyz file
input_water2.write("200\n\n") # this is the number of atoms in the simulation
# modify the numbers in lj_sample_config_periodic1.xyz (well an attempt to anyway)
with open("spce_sample_config_periodic2.xyz") as f: 
	for line_nbr, line in enumerate(f):
		if (line_nbr == 0):
			numbers_str = line.split()
		elif (line_nbr > 1):
			numbers_str = line.split()[1:len(line.split())]
			#numbers_float = [float(x) for x in numbers_str]
			numbers_float = [num for num in numbers_str]
			columns = numbers_float
			print >> input_water2, (str(columns[3]) + "           " + str(columns[0]) +  "          " + str(columns[1]) + "          " + str(columns[2])) 

input_water2.close()

# We will now test the second set of configurations from the water spce thing
# Now, we will replace some lines for fun. 
replace_line("water.inp", 1, "test7_water2.out\n")
replace_line("water.inp", 56, "read_config 200 water2.xyz\n")
replace_line("water.inp", 33, "water.mcf 200\n")
replace_line("water.inp", 39, "20\n")

# open files for read
xyz = open("water2.xyz").read() #Reads in the orginal xyz file (used for read_config)

# Run Cassandra
#print "Runnning 2"
#proc = sp.Popen(["/afs/crc.nd.edu/x86_64_linux/c/cassandra/src/Cassandra_V1.2/Src/cassandra_intel_openMP.exe " + "nist.inp"], stdout=sp.PIPE, shell=True) 
proc = sp.Popen(["/afs/crc.nd.edu/user/l/lpink/Git/Cassandra/Src/cassandra_intel_openMP.exe " + "water.inp"], stdout=sp.PIPE, shell=True) 
(out, err) = proc.communicate()

if err is not None: 
	print("Error.Abort. ")

#print "done 2"

# WATER 3
input_water3 = open("water3.xyz","w") #creates an input xyz file (used for read_config)
# write a new xyz file
input_water3.write("300\n\n") # this is the number of atoms in the simulation
# modify the numbers in lj_sample_config_periodic1.xyz (well an attempt to anyway)
with open("spce_sample_config_periodic3.xyz") as f: 
	for line_nbr, line in enumerate(f):
		if (line_nbr == 0):
			numbers_str = line.split()
		elif (line_nbr > 1):
			numbers_str = line.split()[1:len(line.split())]
			#numbers_float = [float(x) for x in numbers_str]
			numbers_float = [num for num in numbers_str]
			columns = numbers_float
			print >> input_water3, (str(columns[3]) + "           " + str(columns[0]) +  "          " + str(columns[1]) + "          " + str(columns[2])) 

input_water3.close()

# We will now test the second set of configurations from the water spce thing
# Now, we will replace some lines for fun. 
replace_line("water.inp", 1, "test7_water3.out\n")
replace_line("water.inp", 56, "read_config 300 water3.xyz\n")
replace_line("water.inp", 33, "water.mcf 300\n")
replace_line("water.inp", 39, "20\n")

# open files for read
xyz = open("water3.xyz").read() #Reads in the orginal xyz file (used for read_config)

# Run Cassandra
#print "Runnning 3"
#proc = sp.Popen(["/afs/crc.nd.edu/x86_64_linux/c/cassandra/src/Cassandra_V1.2/Src/cassandra_intel_openMP.exe " + "nist.inp"], stdout=sp.PIPE, shell=True) 
proc = sp.Popen(["/afs/crc.nd.edu/user/l/lpink/Git/Cassandra/Src/cassandra_intel_openMP.exe " + "water.inp"], stdout=sp.PIPE, shell=True) 
(out, err) = proc.communicate()

if err is not None: 
	print("Error.Abort. ")

#print "done 3"

# WATER 4
input_water4 = open("water4.xyz","w") #creates an input xyz file (used for read_config)
# write a new xyz file
input_water4.write("750\n\n") # this is the number of atoms in the simulation
# modify the numbers in lj_sample_config_periodic1.xyz (well an attempt to anyway)
with open("spce_sample_config_periodic4.xyz") as f: 
	for line_nbr, line in enumerate(f):
		if (line_nbr == 0):
			numbers_str = line.split()
		elif (line_nbr > 1):
			numbers_str = line.split()[1:len(line.split())]
			#numbers_float = [float(x) for x in numbers_str]
			numbers_float = [num for num in numbers_str]
			columns = numbers_float
			print >> input_water4, (str(columns[3]) + "           " + str(columns[0]) +  "          " + str(columns[1]) + "          " + str(columns[2])) 

input_water4.close()

# We will now test the second set of configurations from the water spce thing
# Now, we will replace some lines for fun. 
replace_line("water.inp", 1, "test7_water4.out\n")
replace_line("water.inp", 56, "read_config 750 water4.xyz\n")
replace_line("water.inp", 33, "water.mcf 750\n")
replace_line("water.inp", 39, "30\n")

# open files for read
xyz = open("water4.xyz").read() #Reads in the orginal xyz file (used for read_config)

# Run Cassandra
#print "Runnning 4"
#proc = sp.Popen(["/afs/crc.nd.edu/x86_64_linux/c/cassandra/src/Cassandra_V1.2/Src/cassandra_intel_openMP.exe " + "nist.inp"], stdout=sp.PIPE, shell=True) 
proc = sp.Popen(["/afs/crc.nd.edu/user/l/lpink/Git/Cassandra/Src/cassandra_intel_openMP.exe " + "water.inp"], stdout=sp.PIPE, shell=True) 
(out, err) = proc.communicate()

if err is not None: 
	print("Error.Abort. ")

#print "done 4"

# Next, the four starting energies will be found in the output .log file for each test. This line will be extracted and saved in this python script so we can use it for comparison at the very end of the script. 
# Finding a string - using an if statement in a for loop
# For the first test
# shakes opens the desired file in the read format

# WATER 1
shakes = open("test7_water1.out.log", "r")
for line in shakes:
	if re.match("(.*)Total system energy(.*)",line):
		line_tot1 = line
	elif re.match("(.*)Intra molecule q(.*)", line):
		line_intraQ1 = line
	elif re.match("(.*)Inter molecule vdw(.*)", line):
		line_vdw1 = line
	elif re.match("(.*)Long range correction(.*)", line):
		line_long1 = line
	elif re.match("(.*)Inter molecule q(.*)", line):
		line_interQ1 = line
	elif re.match("(.*)Reciprocal ewald(.*)", line):
		line_recip1 = line
	elif re.match("(.*)Self ewald(.*)", line):
		line_self1 = line

print "Water 1"
print line_tot1
print line_intraQ1
print line_vdw1
print line_long1
print line_interQ1
print line_recip1
print line_self1

# WATER 2
shakes = open("test7_water2.out.log", "r")
for line in shakes:
	if re.match("(.*)Total system energy(.*)",line):
		line_tot2 = line
	elif re.match("(.*)Intra molecule q(.*)", line):
		line_intraQ2 = line
	elif re.match("(.*)Inter molecule vdw(.*)", line):
		line_vdw2 = line
	elif re.match("(.*)Long range correction(.*)", line):
		line_long2 = line
	elif re.match("(.*)Inter molecule q(.*)", line):
		line_interQ2 = line
	elif re.match("(.*)Reciprocal ewald(.*)", line):
		line_recip2 = line
	elif re.match("(.*)Self ewald(.*)", line):
		line_self2 = line

print "Water 2"
print line_tot2
print line_intraQ2
print line_vdw2
print line_long2
print line_interQ2
print line_recip2
print line_self2

# WATER 3
shakes = open("test7_water3.out.log", "r")
for line in shakes:
	if re.match("(.*)Total system energy(.*)",line):
		line_tot3 = line
	elif re.match("(.*)Intra molecule q(.*)", line):
		line_intraQ3 = line
	elif re.match("(.*)Inter molecule vdw(.*)", line):
		line_vdw3 = line
	elif re.match("(.*)Long range correction(.*)", line):
		line_long3 = line
	elif re.match("(.*)Inter molecule q(.*)", line):
		line_interQ3 = line
	elif re.match("(.*)Reciprocal ewald(.*)", line):
		line_recip3 = line
	elif re.match("(.*)Self ewald(.*)", line):
		line_self3 = line

print "Water 3"
print line_tot3
print line_intraQ3
print line_vdw3
print line_long3
print line_interQ3
print line_recip3
print line_self3


# WATER 4
shakes = open("test7_water4.out.log", "r")
for line in shakes:
	if re.match("(.*)Total system energy(.*)",line):
		line_tot4 = line
	elif re.match("(.*)Intra molecule q(.*)", line):
		line_intraQ4 = line
	elif re.match("(.*)Inter molecule vdw(.*)", line):
		line_vdw4 = line
	elif re.match("(.*)Long range correction(.*)", line):
		line_long4 = line
	elif re.match("(.*)Inter molecule q(.*)", line):
		line_interQ4 = line
	elif re.match("(.*)Reciprocal ewald(.*)", line):
		line_recip4 = line
	elif re.match("(.*)Self ewald(.*)", line):
		line_self4 = line

print "Water 4"
print line_tot4
print line_intraQ4
print line_vdw4
print line_long4
print line_interQ4
print line_recip4
print line_self4

