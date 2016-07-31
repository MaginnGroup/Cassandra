# This is Test 7 in a series of tests for a testSuite in order to check updates made to the Cassandra program. 
#Test 7: This test tests the charge energy from the NIST website using water spce. This test computes the electrostatic energy. It tests this energy for the 4 different configurations from the NIST website which were downloaded and then edited for use in the Cassandra program.

#Note: The test uses a search to extract the energies from the log file created when running Cassandra based on the name of the energy. Furthermore, the extracted energies where than compared to the energies acquired from the NIST website. Also, the energies from the NIST website corresponds to the following Cassandra energies: Cself = Eself, Creciprical = Efourrier, Cintermoleculeq = Eintra + Ereal, Cinter_mol_vdw = Edisp, Clong_range = Elrc, Cintra = 0.0. Where the energies with C in the front are the Cassandra generated energies and the ones with E in the front are from the NIST website. 

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

#print "Water 1"
#print line_tot1
#print line_intraQ1
#print line_vdw1
#print line_long1
#print line_interQ1
#print line_recip1
#print line_self1

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

#print "Water 2"
#print line_tot2
#print line_intraQ2
#print line_vdw2
#print line_long2
#print line_interQ2
#print line_recip2
#print line_self2

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

#print "Water 3"
#print line_tot3
#print line_intraQ3
#print line_vdw3
#print line_long3
#print line_interQ3
#print line_recip3
#print line_self3


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

#print "Water 4"
#print line_tot4
#print line_intraQ4
#print line_vdw4
#print line_long4
#print line_interQ4
#print line_recip4
#print line_self4

# We will now extract the numbers 
# For water 1 - set empty variables
num_tot1 = []
num_intraQ1 = []
num_vdw1 = []
num_long1 = []
num_interQ1 = []
num_recip1 = []
num_self1 = []

# For water 2 - set empty variables
num_tot2 = []
num_intraQ2 = []
num_vdw2 = []
num_long2 = []
num_interQ2 = []
num_recip2 = []
num_self2 = []

# For water 3 - set empty variables
num_tot3 = []
num_intraQ3 = []
num_vdw3 = []
num_long3 = []
num_interQ3 = []
num_recip3 = []
num_self3 = []

# For water 4 - set empty variables
num_tot4 = []
num_intraQ4 = []
num_vdw4 = []
num_long4 = []
num_interQ4 = []
num_recip4 = []
num_self4 = []

# WATER 1
# Use a for loop in order to go through each character in the line independently.
# For the total (I have yet to figure out how to put this all on big loop like above
for t in line_tot1.split():
	try:
		num_tot1.append(float(t))
	except ValueError:
		pass
num_tot1 = num_tot1[0]

# for IntraQ
for t in line_intraQ1.split():

	try:
		num_intraQ1.append(float(t))
	except ValueError:
		pass
num_intraQ1 = num_intraQ1[0]

# for vdw
for t in line_vdw1.split():
	try:
		num_vdw1.append(float(t))
	except ValueError:
		pass
num_vdw1 = num_vdw1[0]

# for long
for t in line_long1.split():
	try:
		num_long1.append(float(t))
	except ValueError:
		pass
num_long1 = num_long1[0]

# for interQ
for t in line_interQ1.split():
	try:
		num_interQ1.append(float(t))
	except ValueError:
		pass
num_interQ1 = num_interQ1[0]

#for recip
for t in line_recip1.split():
	try:
		num_recip1.append(float(t))
	except ValueError:
		pass
num_recip1 = num_recip1[0]

#for self 
for t in line_self1.split():
	try:
		num_self1.append(float(t))
	except ValueError:
		pass
num_self1 = num_self1[0]

#print "Water 1"
#print num_tot1
#print num_intraQ1
#print num_vdw1
#print num_long1
#print num_interQ1
#print num_recip1
#print num_self1

# WATER 2
# Use a for loop in order to go through each character in the line independently.
# For the total (I have yet to figure out how to put this all on big loop like above
for t in line_tot2.split():
	try:
		num_tot2.append(float(t))
	except ValueError:
		pass
num_tot2 = num_tot2[0]

# for IntraQ
for t in line_intraQ2.split():

	try:
		num_intraQ2.append(float(t))
	except ValueError:
		pass
num_intraQ2 = num_intraQ2[0]

# for vdw
for t in line_vdw2.split():
	try:
		num_vdw2.append(float(t))
	except ValueError:
		pass
num_vdw2 = num_vdw2[0]

# for long
for t in line_long2.split():
	try:
		num_long2.append(float(t))
	except ValueError:
		pass
num_long2 = num_long2[0]

# for interQ
for t in line_interQ2.split():
	try:
		num_interQ2.append(float(t))
	except ValueError:
		pass
num_interQ2 = num_interQ2[0]

#for recip
for t in line_recip2.split():
	try:
		num_recip2.append(float(t))
	except ValueError:
		pass
num_recip2 = num_recip2[0]

#for self 
for t in line_self2.split():
	try:
		num_self2.append(float(t))
	except ValueError:
		pass
num_self2 = num_self2[0]

#print "Water 2"
print num_tot2
#print num_intraQ2
#print num_vdw2
#print num_long2
print num_interQ2
#print num_recip2
#print num_self2

# WATER 3
# Use a for loop in order to go through each character in the line independently.
# For the total (I have yet to figure out how to put this all on big loop like above
for t in line_tot3.split():
	try:
		num_tot3.append(float(t))
	except ValueError:
		pass
num_tot3 = num_tot3[0]

# for IntraQ
for t in line_intraQ3.split():

	try:
		num_intraQ3.append(float(t))
	except ValueError:
		pass
num_intraQ3 = num_intraQ3[0]

# for vdw
for t in line_vdw3.split():
	try:
		num_vdw3.append(float(t))
	except ValueError:
		pass
num_vdw3 = num_vdw3[0]

# for long
for t in line_long3.split():
	try:
		num_long3.append(float(t))
	except ValueError:
		pass
num_long3 = num_long3[0]

# for interQ
for t in line_interQ3.split():
	try:
		num_interQ3.append(float(t))
	except ValueError:
		pass
num_interQ3 = num_interQ3[0]

#for recip
for t in line_recip3.split():
	try:
		num_recip3.append(float(t))
	except ValueError:
		pass
num_recip3 = num_recip3[0]

#for self 
for t in line_self3.split():
	try:
		num_self3.append(float(t))
	except ValueError:
		pass
num_self3 = num_self3[0]

#print "Water 3"
print num_tot3
#print num_intraQ3
#print num_vdw3
#print num_long3
print num_interQ3
#print num_recip3
#print num_self3

# WATER 4
# Use a for loop in order to go through each character in the line independently.
# For the total (I have yet to figure out how to put this all on big loop like above
for t in line_tot4.split():
	try:
		num_tot4.append(float(t))
	except ValueError:
		pass
num_tot4 = num_tot4[0]

# for IntraQ
for t in line_intraQ4.split():
	try:
		num_intraQ4.append(float(t))
	except ValueError:
		pass
num_intraQ4 = num_intraQ4[0]

# for vdw
for t in line_vdw4.split():
	try:
		num_vdw4.append(float(t))
	except ValueError:
		pass
num_vdw4 = num_vdw4[0]

# for long
for t in line_long4.split():
	try:
		num_long4.append(float(t))
	except ValueError:
		pass
num_long4 = num_long4[0]

# for interQ
for t in line_interQ4.split():
	try:
		num_interQ4.append(float(t))
	except ValueError:
		pass
num_interQ4 = num_interQ4[0]

#for recip
for t in line_recip4.split():
	try:
		num_recip4.append(float(t))
	except ValueError:
		pass
num_recip4 = num_recip4[0]

#for self 
for t in line_self4.split():
	try:
		num_self4.append(float(t))
	except ValueError:
		pass
num_self4 = num_self4[0]

#print "Water 4"
print num_tot4
#print num_intraQ4
#print num_vdw4
#print num_long4
print num_interQ4
#print num_recip4
#print num_self4

# Lastly, we will compare all these numbers and test the Cassandra Code 
# WATER 1
print "\n"+bold+"Check 1 ( SPCE Water Configuration 1): " + normal
# For Total
if (num_tot1*.01/0.008314 + 4.88604) < 1:
	c1 = 1
	print "Total..."
else: 
	c1 = 0 
	print "Total Fails."

# For intraQ
if num_intraQ1 == 0:
	c2 = 1
	print "Intra Q..."
else: 
	c2 = 0 
	print "Intra Q Fails."


# For vdw
if ((num_vdw1*0.01/0.008314)/1000 - 99.5387) < 1:
	c3 = 1
	print "Vdw..."
else: 
	c3 = 0 
	print "Vdw Fails."

#for long
if ((num_long1*0.01/0.008314)/100 + 8.23715) < 1:
	c4 = 1
	print "Long range..."
else: 
	c4 = 0 
	print "Long Range Fails."

# For interQ
if (num_interQ1/100 + 4351.5) < 4:
	c2 = 1
	print "Inter Q..."
else: 
	c2 = 0
	print "Inter Q Fails."
print num_interQ1
print num_interQ1/((-558889+2809990)*.008314*100)

# For recip
if ((num_recip1*0.01/0.008314)/100 - 62.7000) < 1:
	c6 = 1
	print "Reciprocal..."
else: 
	c6 = 0 
	print "Reciprocal Fails."


# For self
if ((num_self1*0.01/0.008314)/100000 + 28.4469) < 1:
	c7 = 1
	print "Self Ewald..."
else: 
	c7 = 0 
	print "Self Ewald Fails."

#For water 1
if c2 == 1 and c3 == 1 and c4 == 1 and c6 == 1 and c7 == 1:
	w1 = 1
	print bold + "Check 1..." + normal
else:
	w1 = 0 
	print bold + "Check 1 fails. See above." + normal

# WATER 2
print "\n"+bold+"Check 2 ( SPCE Water Configuration 2): " + normal
# For Total
if (num_tot1*.01/0.008314 + 2 ) < 1:
	c1 = 1
	print "Total..."
else: 
	c1 = 0 
	print "Total Fails."

# For intraQ
if num_intraQ2 == 0:
	c2 = 1
	print "Intra Q..."
else: 
	c2 = 0 
	print "Intra Q Fails."


# For vdw
if ((num_vdw2*0.01/0.008314)/10000 - 19.3712) < 1:
	c3 = 1
	print "Vdw..."
else: 
	c3 = 0 
	print "Vdw Fails."

#for long
if ((num_long2*0.01/0.008314)/100 + 32.9486) < 1:
	c4 = 1
	print "Long range..."
else: 
	c4 = 0 
	print "Long Range Fails."

# For interQ
if (num_interQ2/100 + 4351.5) < 4:
	c2 = 1
	print "Inter Q..."
else: 
	c2 = 0 
	print "Inter Q Fails."

print num_interQ2
print num_interQ2/((-1192950+5619980)*.008314*100)

# For recip
if ((num_recip2*0.01/0.008314)/100 - 60.349) < 1:
	c6 = 1
	print "Reciprocal..."
else: 
	c6 = 0 
	print "Reciprocal Fails."


# For self
if ((num_self2*0.01/0.008314)/100000 + 56.8938) < 1:
	c7 = 1
	print "Self Ewald..."
else: 
	c7 = 0 
	print "Self Ewald Fails."


#For water 1
if c2 == 1 and c3 == 1 and c4 == 1 and c6 == 1 and c7 == 1:
	w2 = 1
	print bold + "Check 2..." + normal
else:
	w2 = 0 
	print bold + "Check 2 fails. See above." + normal

# WATER 3
print "\n"+bold+"Check 3 ( SPCE Water Configuration 3): " + normal
# For Total
if (num_tot1*.01/0.008314 + 2 ) < 1:
	c1 = 1
	print "Total..."
else: 
	c1 = 0 
	print "Total Fails."

# For intraQ
if num_intraQ3 == 0:
	c2 = 1
	print "Intra Q..."
else: 
	c2 = 0 
	print "Intra Q Fails."


# For vdw
if ((num_vdw3*0.01/0.008314)/10000 -35.4344) < 1:
	c3 = 1
	print "Vdw..."
else: 
	c3 = 0 
	print "Vdw Fails."

#for long
if ((num_long3*0.01/0.008314)/100 + 74.1343) < 1:
	c4 = 1
	print "Long range..."
else: 
	c4 = 0 
	print "Long Range Fails."

# For interQ
if (num_interQ3/100 + 4351.5) < 4:
	c2 = 1
	print "Inter Q..."
else: 
	c2 = 0 
	print "Inter Q Fails."

print num_interQ3
print num_interQ3/((-1962970+8429980)*.008314*100)

# For recip
if ((num_recip3*0.01/0.008314)/100 - 52.4461) < 1:
	c6 = 1
	print "Reciprocal..."
else: 
	c6 = 0 
	print "Reciprocal Fails."


# For self
if ((num_self3*0.01/0.008314)/100000 + 85.3407) < 1:
	c7 = 1
	print "Self Ewald..."
else: 
	c7 = 0 
	print "Self Ewald Fails."


#For water 1
if c2 == 1 and c3 == 1 and c4 == 1 and c6 == 1 and c7 == 1:
	w3 = 1
	print bold + "Check 3..." + normal
else:
	w3 = 0 
	print bold + "Check 3 fails. See above." + normal

# WATER 4
print "\n"+bold+"Check 4 ( SPCE Water Configuration 4): " + normal
# For Total
if (num_tot1*.01/0.008314 + 2 ) < 1:
	c1 = 1
	print "Total..."
else: 
	c1 = 0 
	print "Total Fails."

# For intraQ
if num_intraQ4 == 0:
	c2 = 1
	print "Intra Q..."
else: 
	c2 = 0 
	print "Intra Q Fails."


# For vdw
if ((num_vdw4*0.01/0.008314)/10000 - 44.8593) < 1:
	c3 = 1
	print "Vdw..."
else: 
	c3 = 0 
	print "Vdw Fails."

#for long
if ((num_long4*0.01/0.008314)/100 + 13.7286) < 1:
	c4 = 1
	print "Long range..."
else: 
	c4 = 0 
	print "Long Range Fails."

# For interQ
if (num_interQ4/100 + 4351.5) < 4:
	c2 = 1
	print "Inter Q..."
else: 
	c2 = 0 
	print "Inter Q Fails."


print num_interQ4
print num_interQ4/((-3572260+14148300)*.008314*100)

# For recip
if ((num_recip4*0.01/0.008314)/100 - 75.878) < 2:
	c6 = 1
	print "Reciprocal..."
else: 
	c6 = 0 
	print "Reciprocal Fails."


# For self
if ((num_self4*0.01/0.008314)/1000000 + 14.1483) < 1:
	c7 = 1
	print "Self Ewald..."
else: 
	c7 = 0 
	print "Self Ewald Fails."


#For water 4
if c2 == 1 and c3 == 1 and c4 == 1 and c6 == 1 and c7 == 1:
	w4 = 1
	print bold + "Check 4..." + normal
else:
	w4 = 0 
	print bold + "Check 4 fails. See above." + normal

#check all Water
if w1 == 1 and w2 == 1 and w3 == 1 and w4 == 1:
	print bold + "Pass Test 7: Charge (Water SPCE) Energy." + normal
else: 
	print bold + "Test 7 fails." + normal
