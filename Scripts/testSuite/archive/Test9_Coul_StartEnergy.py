#!/usr/bin/env python

#*******************************************************************************
# SCRIPT:  Test9_Coul_StartEnergy.py
# VERSION: 2.0
# FEATURES: Compute the electrostatic energy of a dipole
#*******************************************************************************

#*******************************************************************************
# IMPORT MODULES
#*******************************************************************************
import subprocess as sp #This module lets us run cassandra from python
import numpy as np #this module s the package for scientific computing in python
import os, sys

#*******************************************************************************
# ARGUMENT PARSE
#*******************************************************************************

#*******************************************************************************
# VARIABLE DEFINITIONS
#*******************************************************************************
# Test description
test_no = 9
test_desc = "Coulombic energy of one dipole"

# Physical constants
# to calculate the analytic answer
elementary_charge = 1.602176565 # 10^-19, C
epsilon_0 = 8.854187817 # 10^-12, C^2 / N m^2
avogadro = 6.02214129 # 10^23, /mol
charge_factor = elementary_charge**2 * avogadro/(4*np.pi*epsilon_0) * 10000 # kJ A / mol

# Simulation parameters
nSpecies = 1
nMols = (1,)
nAtoms = (2,)
atomName = (('C','A'),)
atomCharge = ((1.0,-1.0),)
nBonds = (1,)
bondParms = (((1,2,"fixed",1.0),),)
box = 30

# Check parameters
# params for each simulation
numChecks = 2 # number of simulations to run
vdwStyle = ('lj minimum_image', 'none')
chargeStyle = ('coul minimum_image', 'coul ewald 10. 1e-5')
cassStr = (("Total system energy",),
           ("Total system energy", "Inter molecule q", "Reciprocal ewald", "Self ewald"))
analyticAnswer = [None] * numChecks #this list will hold the analytic answers
cassAnswer     = [None] * numChecks #this list will hold cassandra's answers
passCheck      = [None] * numChecks
errorTol = 1e-5

# Formatting variables
bold = '\033[1m' #Will make text bold
normal = '\033[0m' #Will make the next text normal(ie. unbold)

#*******************************************************************************
# FILE MANAGEMENT
#*******************************************************************************
cassDir = "/Users/rmullen2/dev/cassandra/Src/"
cassExe = "cassandra.dev"
cassRun = "test9.out"
inpName = "test9.inp"
xyzName = "test9.inp.xyz"
mcfName = ("test9.dipole.mcf",)


#*******************************************************************************
# FUNCTION DEFINITIONS
#*******************************************************************************
# This functon writes over specific lines of the input file created above using a function, called replace_line that is created below. 
# This function takes three inputs: the name of the file where you would like to replace a line, the line number you would like to replace, and the text you want to replace the old text with. 
def replace_line(file_name, line_num, text):
	lines = open(file_name, 'r').readlines()
	lines[line_num] = text
	out = open(file_name, 'w')
	out.writelines(lines)
	out.close() # Closes the file so that the program doesn't explode. 

# Error function
def erf(x):
	# constants
	a1 =  0.254829592
	a2 = -0.284496736
	a3 =  1.421413741
	a4 = -1.453152027
	a5 =  1.061405429
	p  =  0.3275911

	# Save the sign of x
	sign = 1
	if x < 0:
		sign = -1
	x = abs(x)

	# A & S 7.1.26
	t = 1.0/(1.0 + p*x)
	y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*np.exp(-x*x)

	return sign*y
#*******************************************************************************
# MAIN PROGRAM BEGINS HERE
#*******************************************************************************
#
# Step 1) Write input files
#
# Loop through the checks
# 
# Step 2) Calculate the correct answer
# Step 3) Run Cassandra to get its answer
# Step 4) Compare answers
#

#This prints the starting line.
print "\n\n"+bold+"Test " + str(test_no) +": " + test_desc + normal

# Step 1) Write input files
# 1.1) Write MCF for cation, anion
for s in range(nSpecies):
	mcf = open(mcfName[s],"w") #Creates mcf file and allows for edits
	mcf.write("# Atom_Info\n%d\n" % (nAtoms[s]))
	for i in range(nAtoms[s]):
		mcf.write("%d    %s    %s   1.0 %.1f   NONE\n" % (i+1,atomName[s][i],atomName[s][i],atomCharge[s][i]))
	mcf.write("\n# Bond_Info\n%d\n" % (nBonds[s]))
	for i in range(nBonds[s]):
		mcf.write("%d    %d    %d    %s    %9.3f\n" % (i+1,bondParms[s][i][0],bondParms[s][i][1],
                                                       bondParms[s][i][2],bondParms[s][i][3]))
	mcf.write("\n# Angle_Info\n0\n")
	mcf.write("\n# Dihedral_Info\n0\n")
	mcf.write("\n# Improper_Info\n0\n")
	mcf.write("\n# Fragment_Info\n0\n")
	mcf.write("\n# Fragment_Connectivity\n0\n")
	mcf.write("\n# Intra_Scaling\n0. 0. 0. 0.\n0. 0. 0. 0.\n")
	mcf.write("\nEND\n")
	mcf.close()

# Loop through checks
print "%-20s %-20s %18s %18s %18s %8s" % ("Property","ChargeStyle","Cassandra","Analytic","Relative_Err","Pass")
for i in range(numChecks):

	# Step 2) Calculate the correct answer
	if (chargeStyle[i] == 'coul minimum_image'):
		# Coulomb's law
		analyticAnswer[i] = [0.]
	elif ('ewald' in chargeStyle[i]):
		# Ewald summation
		ewald_tol = float(chargeStyle[i].split()[-1])
		rcut = float(chargeStyle[i].split()[2])
		alpha = np.sqrt(-np.log(ewald_tol)) / rcut
		k2_cut = (np.log(ewald_tol) / rcut / np.pi)**2
		nmax = int( -2.0 * np.log(ewald_tol) / rcut * box / (2 * np.pi)) + 1

		analyticAnswer[i] = [0.] * 4
		# Real
		analyticAnswer[i][1] = atomCharge[0][0] * atomCharge[0][1] * charge_factor * (- erf(alpha))
		# Reciprocal
		recip = 0.
		for nx in range(-nmax,nmax):
			kx = float(nx) / box
			real2 = (atomCharge[0][0] + atomCharge[0][1] * np.cos(2 * np.pi * kx))**2
			im2 = (atomCharge[0][1] * np.sin(2 * np.pi * kx))**2
			for ny in range(-nmax,nmax):
				ky = float(ny) / box
				for nz in range(-nmax,nmax):
					kz = float(nz) / box
					k2 = kx**2 + ky**2 + kz**2
					if (k2 == 0. or k2 > k2_cut):
						continue
					prefactor = 1.0/float(k2) * np.exp(-(np.pi/alpha)**2 * float(k2))
					recip = recip + prefactor * (real2 + im2)
		analyticAnswer[i][2] = recip * charge_factor / (2 * np.pi * box**3)
		# Self
		analyticAnswer[i][3] = - alpha / np.sqrt(np.pi) * (atomCharge[0][0]**2 + atomCharge[0][1]**2) * charge_factor
		# Total
		analyticAnswer[i][0] = sum(analyticAnswer[i][1:4])

	# Step 3) Run Cassandra to get its answer
	# 3.1) Write inp file

	# Combine file names and nmols
	mcfStr = ''
	for s in range(nSpecies):
		mcfStr = mcfStr + mcfName[s] + " " + str(nMols[s]) + '\n'

	# Write inp file
	inp = open(inpName,"w")
	inp.write("# Run_Name\n%s\n\n" % (cassRun))
	inp.write("# Sim_Type\nnvt\n\n")
	inp.write("# Nbr_Species\n%d\n\n" % (nSpecies))
	inp.write("# VDW_Style\n%s\n\n" % (vdwStyle[i]))
	inp.write("# Charge_Style\n%s\n\n" % (chargeStyle[i]))
	inp.write("# Seed_Info\n1 2\n\n")
	inp.write("# Rcutoff_Low\n0.1\n\n")
	inp.write("# Molecule_Files\n%s\n" % (mcfStr))
	inp.write("# Box_Info\n1\ncubic\n%.1f\n\n" % (box))
	inp.write("# Temperature_Info\n1.0\n\n")
	inp.write("# Move_Probability_Info\n\n")
	inp.write("# Prob_Translation\n1.0\n0.0 0.0\n\n")
	inp.write("# Done_Probability_Info\n\n")
	inp.write("# Start_Type\nread_config %s %s\n\n" % (' '.join(str(x) for x in nMols),xyzName))
	inp.write("# Run_Type\nEquilibration 100\n\n")
	inp.write("# Simulation_Length_Info\nunits steps\nprop_freq 1\ncoord_freq 1\nrun 0\n\n")
	inp.write("END\n")
	inp.close()

	# 3.1) Write xyz
	xyz = open(xyzName,"w")
	xyz.write("2\n\n") # This is the number of atoms in the simulation 
	xyz.write(" C  0.0  0.0  0.0\n") #Location of atom 1
	xyz.write(" A  1.0  0.0  0.0\n") #Location of atom 2
	xyz.close()

	# 3.2) Run Cassandra Jobs
	proc = sp.Popen([cassDir + cassExe + " " + inpName], stdout=sp.PIPE, shell=True)
	(out, err) = proc.communicate()

	if err is not None:
		print("Error.Abort.")

	# 3.3) Read logfile
	log = open(cassRun + ".log", "r")
	# search line by line in log for the words "Total system energy"
	nPrp = len(cassStr[i])
	cassAnswer[i] = [None] * nPrp
	for line in log:
		for j in range(nPrp):
			if (cassStr[i][j] in line):
				cassAnswer[i][j] = float(line.split()[-1])

	# Step 4) Compare answers
	for j in range(nPrp):
		if (analyticAnswer[i][j] == 0.):
			passCheck = cassAnswer[i][j] == 0.
			print "%-20s %-20s %18.6g %18.6g %18s %8s" % (cassStr[i][j],chargeStyle[i],
						cassAnswer[i][j],analyticAnswer[i][j],'',passCheck)
		else:
			errorRel = abs(cassAnswer[i][j] - analyticAnswer[i][j])/analyticAnswer[i][j]
			passCheck = abs(errorRel) <= errorTol
			print "%-20s %-20s %18.6g %18.6g %18.6g %8s" % (cassStr[i][j],chargeStyle[i],
						cassAnswer[i][j],analyticAnswer[i][j],errorRel,passCheck)
	print ""

# Clean up scratch files
os.system('rm ' + inpName)
os.system('rm ' + xyzName)
os.system('rm ' + ' '.join(mcfName))
os.system('rm ' + cassRun + '*')
