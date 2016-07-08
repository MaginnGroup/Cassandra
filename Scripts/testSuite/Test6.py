
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
print "\n\n"+bold+"Test 6: Nist something Commencing ... " + normal 
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
print "running"
#proc = sp.Popen(["/afs/crc.nd.edu/x86_64_linux/c/cassandra/src/Cassandra_V1.2/Src/cassandra.exe " + "nist.inp"], stdout=sp.PIPE, shell=True) 
proc = sp.Popen(["/afs/crc.nd.edu/user/l/lpink/Git/Cassandra/Src/cassandra.exe " + "nist.inp"], stdout=sp.PIPE, shell=True) 
(out, err) = proc.communicate()

if err is not None: 
	print("Error.Abort. ")
print "done"

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

# open files for read
xyz = open("nist2.xyz").read() #Reads in the orginal xyz file (used for read_config)


# Run Cassandra
print "running 2"
proc = sp.Popen(["/afs/crc.nd.edu/user/l/lpink/Git/Cassandra/Src/cassandra.exe " + "nist.inp"], stdout=sp.PIPE, shell=True) 
(out, err) = proc.communicate()

if err is not None: 
	print("Error.Abort. ")
print "done 2"

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


# open files for read
xyz = open("nist3.xyz").read() #Reads in the orginal xyz file (used for read_config)


# Run Cassandra
print "running 3"
#proc = sp.Popen(["/afs/crc.nd.edu/x86_64_linux/c/cassandra/src/Cassandra_V1.2/Src/cassandra.exe " + "nist.inp"], stdout=sp.PIPE, shell=True) 
proc = sp.Popen(["/afs/crc.nd.edu/user/l/lpink/Git/Cassandra/Src/cassandra.exe " + "nist.inp"], stdout=sp.PIPE, shell=True) 
(out, err) = proc.communicate()

if err is not None: 
	print("Error.Abort. ")
print "done 3"

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

# open files for read
xyz = open("nist4.xyz").read() #Reads in the orginal xyz file (used for read_config)


# Run Cassandra
print "running"
##proc = sp.Popen(["/afs/crc.nd.edu/x86_64_linux/c/cassandra/src/Cassandra_V1.2/Src/cassandra.exe " + "nist.inp"], stdout=sp.PIPE, shell=True) 
proc = sp.Popen(["/afs/crc.nd.edu/user/l/lpink/Git/Cassandra/Src/cassandra.exe " + "nist.inp"], stdout=sp.PIPE, shell=True) 
(out, err) = proc.communicate()

if err is not None: 
	print("Error.Abort. ")
print "done"



# CHECK 2 - This is for a cut_tail with 4 sigma and again using all 4 nist information
# Now, we will replace some lines for fun. 
replace_line("nist.inp", 1, "test6_check2_nist1.out\n")
replace_line("nist.inp", 13, "LJ cut_tail 14.16\n")
replace_line("nist.inp", 60, "read_config 800 nist1.xyz\n")
replace_line("nist.inp", 37, "nist.mcf 800\n")



# Run Cassandra
print "running 5"
#proc = sp.Popen(["/afs/crc.nd.edu/x86_64_linux/c/cassandra/src/Cassandra_V1.2/Src/cassandra.exe " + "nist.inp"], stdout=sp.PIPE, shell=True) 
proc = sp.Popen(["/afs/crc.nd.edu/user/l/lpink/Git/Cassandra/Src/cassandra.exe " + "nist.inp"], stdout=sp.PIPE, shell=True) 
(out, err) = proc.communicate()

if err is not None: 
	print("Error.Abort. ")
print "done 5"

# Check 2 - Nist 2. 
replace_line("nist.inp", 1, "test6_check2_nist2.out\n")
replace_line("nist.inp", 60, "read_config 200 nist2.xyz\n")
replace_line("nist.inp", 37, "nist.mcf 200\n")

# Run Cassandra
print "running 6"
##proc = sp.Popen(["/afs/crc.nd.edu/x86_64_linux/c/cassandra/src/Cassandra_V1.2/Src/cassandra.exe " + "nist.inp"], stdout=sp.PIPE, shell=True) 
proc = sp.Popen(["/afs/crc.nd.edu/user/l/lpink/Git/Cassandra/Src/cassandra.exe " + "nist.inp"], stdout=sp.PIPE, shell=True) 
(out, err) = proc.communicate()

if err is not None: 
	print("Error.Abort. ")
print "done 6"

# Check 2 - Nist 3. 
replace_line("nist.inp", 1, "test6_check2_nist3.out\n")
replace_line("nist.inp", 60, "read_config 400 nist3.xyz\n")
replace_line("nist.inp", 37, "nist.mcf 400\n")

# Run Cassandra
print "running 7"
#proc = sp.Popen(["/afs/crc.nd.edu/x86_64_linux/c/cassandra/src/Cassandra_V1.2/Src/cassandra.exe " + "nist.inp"], stdout=sp.PIPE, shell=True) 
proc = sp.Popen(["/afs/crc.nd.edu/user/l/lpink/Git/Cassandra/Src/cassandra.exe " + "nist.inp"], stdout=sp.PIPE, shell=True) 
(out, err) = proc.communicate()

if err is not None: 
	print("Error.Abort. ")
print "done 7"

# Check 2 - Nist 4. 
replace_line("nist.inp", 1, "test6_check2_nist4.out\n")
replace_line("nist.inp", 60, "read_config 30 nist4.xyz\n")
replace_line("nist.inp", 37, "nist.mcf 30\n")

## Run Cassandra
print "running 8"
##proc = sp.Popen(["/afs/crc.nd.edu/x86_64_linux/c/cassandra/src/Cassandra_V1.2/Src/cassandra.exe " + "nist.inp"], stdout=sp.PIPE, shell=True) 
proc = sp.Popen(["/afs/crc.nd.edu/user/l/lpink/Git/Cassandra/Src/cassandra.exe " + "nist.inp"], stdout=sp.PIPE, shell=True) 
(out, err) = proc.communicate()

if err is not None: 
	print("Error.Abort. ")
print "done 8"


