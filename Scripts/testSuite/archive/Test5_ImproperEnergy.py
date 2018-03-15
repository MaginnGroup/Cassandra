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
from itertools import *

#We will now create an input file (.inp file)
bold = '\033[1m' #Will make text bold
normal = '\033[0m' #Will make the next text normal(ie. unbold)

#This prints the starting line.
print "\n\n"+bold+"Test 5: Improper Angle Energy Commencing ... " + normal 


# Try to read in a certain number of lines of a file at a time.
#with open("trajectoryImproper.xyz", 'r') as points:
#f = open("trajectoryImproper.xyz", "r")

n = 1
k = 14

open("results_improper.txt", "w").close()


# We will attempt to do a loop here that does all of this like 20 times to start (actually needs to be done 1000 times)

n = 1

for n in range(1, 1002):
	#Creates(opens if already existent) the file file.inp which will be our input file while running cassandra
	# The input files are all angle. because it deals with the angle test.
	input_inp = open("improper.inp","w")
	# Read in mcf and xyz files to use in inp generation
	input_mcf = open("improper.mcf","w") #Creates mcf file and allows for edits
	#input_xyz = open("improper.xyz","w") #Creates an input xyz file (used for read_config)

	#Write info for the MCF file.
	input_mcf.write("# Atom_Info\n12\n1   C1_s1   C   12.0101   -0.115000   LJ   35.226   3.550   ring\n2   C2_s1   C   12.010   -0.115000   LJ   35.226   3.550   ring\n3   C3_s1   C   12.010   -0.115000   LJ   35.226   3.550   ring\n4   C4_s1   C   12.010   -0.115000   LJ   35.226   3.550   ring\n5   C5_s1   C   12.010   -0.115000   LJ   35.226    3.550   ring\n6   C6_s1   C   12.010   -0.115000   LJ   35.226   3.550   ring\n7   H1_s1   H   1.008   0.115000   LJ   15.097    2.420\n8   H2_s1   H   1.008   0.115000   LJ   15.097    2.420\n9   H3_s1   H   1.008    0.115000   LJ   15.097   2.420\n10   H4_s1   H   1.008   0.115000   LJ   15.097   2.420\n11   H5_s1   H   1.008   0.115000   LJ   15.097   2.420\n12   H6_s1   H   1.008   0.115000   LJ   15.097   2.420\n\n")
	input_mcf.write("# Bond_Info\n12\n1   1   2   fixed   1.394\n2   1   6   fixed   1.394\n3   1   7   fixed   1.080\n4   2   3   fixed   1.394\n5   2   8   fixed   1.080\n6   3   4   fixed   1.394\n7   3   9   fixed  1.080\n8   4   5  fixed   1.394\n9   4   10 fixed   1.080\n10   5   6   fixed   1.394\n11   5   11   fixed   1.080\n12   6   12   fixed   1.080\n\n")
	input_mcf.write("# Angle_Info\n18\n1 1 2 3 harmonic  40257.8  120.00\n2 1 2 8 harmonic  30193.4  120.00\n3 1 6 5 harmonic  40257.8  120.0\n4 1 6 12 harmonic  30193.4  120.00\n5 2 1 6 harmonic  40257.8  120.00\n6 2 1 7 harmonic 30193.4  120.00\n7 2 3 4 harmonic  40257.8  120.0\n8 2 3 9 harmonic  30193.4  120.0\n9 3 2 8 harmonic  30193.4  120.0\n10 3 4 5 harmonic  40257.8  120.00\n11 3 4 10 harmonic  30193.4  120.00\n12 4 3 9 harmonic  30193.4  120.0\n13 4 5 6 harmonic  40257.8  120.00\n14 4 5 11 harmonic  30193.4  120.00\n15 5 4 10 harmonic  30193.4  120.00\n16 5 6 12 harmonic  30193.4  120.00\n17 6 1 7 harmonic  30193.4  120.00\n18 6 5 11 harmonic  30193.4  120.00\n\n")
	input_mcf.write("# Dihedral_Info\n24\n1 6 1 2 3 CHARMM 12.9700774517  2  180.0\n2 7 1 2 3 CHARMM 17.5731271258  2  180.0\n3 6 1 2 8 CHARMM 17.5731271258  2  180.0\n4 7 1 2 8 CHARMM 10.0420779356  2  180.0\n5 2 1 6 5 CHARMM 12.9700774517  2  180.0\n6 7 1 6 5 CHARMM 17.5731271258  2  180.0\n7 2 1 6 12 CHARMM 17.5731271258  2  180.0\n8 7 1 6 12 CHARMM 10.0420779356  2  180.0\n9 1 2 3 4 CHARMM 12.9700774517  2  180.0\n10 8 2 3 4 CHARMM 17.5731271258  2.0  180.0\n11 1 2 3 9 CHARMM 17.5731271258  2  180.0\n12 8 2 3 9 CHARMM 10.0420779356  2  180.00\n13 2 3 4 5 CHARMM 12.9700774517  2  180.0\n14 9 3 4 5 CHARMM 17.5731271258  2  180.0\n15 2 3 4 10 CHARMM 17.5731271258  2  180.0\n16 9 3 4 10 CHARMM 10.0420779356  2  180.0\n17 3 4 5 6 CHARMM 12.9700774517  2  180.0\n18 10 4 5 6 CHARMM 17.5731271258  2  180.0\n19 3 4 5 11 CHARMM 17.5731271258  2.0  180.0\n20 10 4 5 11 CHARMM 17.5731271258  2.0  180.0\n21 4 5 6 1 CHARMM 12.9700774517  2  180.0\n22 11 5 6 1 CHARMM 17.5731271258  2 180.0\n23 4 5 6 12 CHARMM 17.5731271258  2  180.0\n24 11 5 6 12 CHARMM  10.0420779356  2  180.0\n\n")
	input_mcf.write("# Improper_Info\n6\n1 2 6 1 7 cvff 17.5731271258  -1 2\n2 1 3 2 8 cvff 17.5731271258   -1 2\n3 2 4 3 9 cvff 17.5731271258  -1 2\n4 3 5 4 10 cvff 17.5731271258  -1 2\n5 4 6 5 11 cvff 17.5731271258  -1 2\n6 1 5 6 12 cvff 17.5731271258  -1 2\n\n")
	input_mcf.write("# Fragment_Info\n0\n\n")
	input_mcf.write("# Fragment_Connectivity\n0\n\n")
	input_mcf.write("# Intra_Scaling\n0. 0. 0.0000 1.\n0. 0. 0.0000 1.\n\n\n")
	input_mcf.write("END")
	input_mcf.close()

	# Generate your xyz file again and again and again
	#for n in range(1,4):
	with open("trajectoryImproper.xyz", "r") as points:
			fail = open("fail.xyz", "w")
			data = iter(points)
			#print ''.join(list(islice(data,k-14,k)))
			#print 'new config'
			fail.write(''.join(list(islice(data, k-14, k))))
			k = k + 14
			fail.close()


	#Write info for the XYZ file. 
	#input_xyz.write("1\n\n") # This is the number of atoms in the simulation 
	#input_xyz.write("C      -0.697  -1.208  0.0005 \n") #Location of atom 1
	#input_xyz.write("C       0.698  -1.208  0.0015 \n") #Location of atom 2
	#input_xyz.write("C      1.395  0  0.0005 \n") #Location of atom 3
	#input_xyz.write("C      0.698  1.207  -0.0005 \n") #Location of atom 4
	#input_xyz.write("C      -0.697  1.207  -0.0005 \n") #Location of atom 5
	#input_xyz.write("C      -1.394  0  -0.0005 \n") #Location of atom 6
	#input_xyz.write("H      -1.239  -2.147  0.0005 \n") #Location of atom 7
	#input_xyz.write("H      1.24  -2.147  0.0015 \n") #Location of atom 8
	#input_xyz.write("H      2.479  -0.001  0.0005 \n") #Location of atom 9
	#input_xyz.write("H      1.24  2.146  -0.0005 \n") #Location of atom 10
	#input_xyz.write("H      -1.239  2.147  -0.0015 \n") #Location of atom 11
	#input_xyz.write("H      -2.479  0  -0.0005 \n") #Location of atom 12
	#input_xyz.close()

	#Write info into the file - this will create each section for the .inp file
	# This input file is populated with numbers for an LJ simulation with argon
	input_inp.write("# Run_Name\ntest5_improper.out\n!---------------\n\n")
	input_inp.write("# Sim_Type\nNVT_MC\n!---------------\n\n")
	input_inp.write("# Nbr_Species\n1\n!---------------\n\n")
	input_inp.write("# VDW_Style\nLJ cut_tail 14.0\n!---------------\n\n")
	input_inp.write("# Charge_Style\ncoul ewald 14.0 0.0000001\n!---------------\n\n")
	input_inp.write("# Intra_Scaling\n0.0 0.0 0.0 1.0\n!---------------vdw, coul lines for each species\n\n")
	input_inp.write("# Mixing_Rule\nLB\n!---------------\n\n")
	input_inp.write("# Seed_Info\n1211131639 1211131640\n!---------------\n\n")
	input_inp.write("# Rcutoff_Low\n1.0\n!---------------\n\n")
	input_inp.write("# Molecule_Files\nimproper.mcf 1\n!----------------\n\n")
	input_inp.write("# Box_Info\n1\nCUBIC\n50.0\n!---------------\n\n")
	input_inp.write("# Temperature_Info\n510.0\n!---------------\n\n")
	input_inp.write("# Move_Probability_Info\n\n")
	input_inp.write("# Prob_Translation\n0.40\n0.50\n14.0\n\n")
	input_inp.write("# Done_Probability_Info\n!----------------\n\n")
	input_inp.write("# Start_Type\nread_config 1 fail.xyz\n!---------------\n\n")
	input_inp.write("# Run_Type\nEquilibration 100\n!---------------\n\n")
	input_inp.write("# Average_Info\n1\n!---------------(0 = yes, 1 = no)\n\n")
	input_inp.write("# Simulation_Length_Info\nUnits Steps\nProp_Freq 1\nCoord_Freq 1\nRun 1\n!---------------\n\n")
	input_inp.write("# Property_Info 1\nEnergy_Total\nDensity\nNmols\nVolume\nPressure\n!---------------\n\n")
	input_inp.write("# Fragment_Files\n!---------------one line per fragment\n\n")
	input_inp.write("# CBMC_Info\nkappa_ins 12\nkappa_rot 0\nkappa_dih 12\nrcut_cbmc 6.5 6.5\n!---------------\n\n")
	input_inp.write("# Pair_Energy\nTRUE\n\n")
	input_inp.write("END")
	input_inp.close()

	# Read in input files
	inp = open("improper.inp").read() #This command reads in the input file
	mcf = open("improper.mcf").read() #Reads in the mcf file
	xyz = open("fail.xyz").read() #Reads in the orginal xyz file (used for read_config)
 
	# Run Cassandra
	#print "Runnning"
	#proc = sp.Popen(["/afs/crc.nd.edu/x86_64_linux/c/cassandra/src/Cassandra_V1.2/Src/cassandra_intel_openMP.exe " + "nist.inp"], stdout=sp.PIPE, shell=True) 
	proc = sp.Popen(["/afs/crc.nd.edu/user/l/lpink/Git/Cassandra/Src/cassandra_intel_openMP.exe " + "improper.inp"], stdout=sp.PIPE, shell=True) 
	(out, err) = proc.communicate()

	if err is not None: 
		print("Error.Abort. ")

	#print "done"

	# This fucntion takes three inputs: the name of the file where you would like to replace a line, the line number you would like to replace, and the text you want to replace the old text with. 
	def replace_line(file_name, line_num, text):
		lines = open(file_name, 'r').readlines()
		lines[line_num] = text
		out = open(file_name, 'w')
		out.writelines(lines)
		out.close() # Closes the file so that the program doesn't explode. 




	# Extract the line of improper number
	shakes = open("test5_improper.out.log", "r")
	for line in shakes:
		if re.match("(.*)Improper angle energy(.*)",line):
			line_im = line

	# Now, take just the number
	num = []
	for t in line_im.split():
		try:
			num.append(float(t))
		except ValueError:
			pass
	num = num[0]


	#print num
	#f = open("results_improper.txt", "w")
	num_str = str(num)

	with open("results_improper.txt", "a") as myfile:
		myfile.write(num_str + "\n")

	n += 1



results = open("givenResults.txt", "w")
with open("givenResults_improper.txt", "r") as inf:
	for line in inf:
		parts = line.split()
		if len(parts) > 1:
			print >> results, parts[1]

results.close()

# Now loop through for your checks, hopefully this works:
given_results = open("givenResults.txt", "r")
calc_results = open("results_improper.txt", "r")

given = list(given_results)
calc = list(calc_results)

#print given[0]
#print calc[0]

given_num = [float(i) for i in given]
calc_num = [float(i) for i in calc]

#print given_num[1]
#print calc_num[1]

# HERE IT GOES

check = open("check_test5.txt", "w")
n = 1
for n in range(0, 1001):
	if given_num[n] - calc_num[n] < 1: 
		diff = 1	
		n += 1
		print >> check, "pass"
	else:
		diff = 0
		n += 1
		print >> check, "fail"


shakes = open("check_test5.txt", "r")
for line in shakes: 
	if re.match("(.*)fail(.*)", line):
		print "Test 5 Fails. See check_test5.txt file."
	else: 
		c = 1

if c == 1:
	print "Pass..."	 
	print bold + "Pass Test 5: Improper Angle Energy" + normal


# Now, we will find a way to plot:

