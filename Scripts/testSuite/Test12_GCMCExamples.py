#!/usr/bin/env python 
#*******************************************************************************
# SCRIPT:  Test9_NPTexamples.py
# VERSION: 1.0
# FEATURES: Test the results of a short simulation of the NPT examples
#*******************************************************************************

#*******************************************************************************
# IMPORT MODULES
#*******************************************************************************
import subprocess as sp #This module lets us run cassandra from python
import numpy as np #this module's the package for scientific computing in python
import random #This allows us to run random numbers
import os, sys
import re
from itertools import *
import argparse

#*******************************************************************************
# ARGUMENT PARSE
#*******************************************************************************

parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=
"""DESCRIPTION:
Runs the given test using the Cassandra executable specified.

EXAMPLES:
To run a test using a Cassandra executable inside of Cassandra/Src/:

	> python Test#_Description.py cassandra.exe
	> python Test#_Description.py cassandra_gfortran.exe

To run a test using a Cassandra executable elsewhere:

	> python Test#_Description.py /home/applications/cassandra.exe --absPath
	> python Test#_Description.py /home/applications/cassandra.exe -a


""")
parser.add_argument('cassandra_exe', 
                help="Cassandra executable [file name if inside Src/ directory or path with " +
                "indicated --absPath flag]")
parser.add_argument('--absPath','-a', action='store_true',
                help="Signals that the Cassandra executable is given as a path instead of " +
                "given as is in the Src driectory.")

args = parser.parse_args()


#*******************************************************************************
# VARIABLE DEFINITIONS
#*******************************************************************************
# Test description
test_no = 12
test_desc = "GCMC Examples"

# Check parameters
nChecks = 6 # number of simulations to run
title = ("Methane", "Methane_Butane", "Methane_Butane_Silicalite", "Methane_Silicalite", "Nitrogen", "Nitrogen_Silicalite")
analyticAnswer = [None] * nChecks # list to hold analytic answers
cassAnswer     = [None] * nChecks # list to hold cassandra's answers
passCheck      = [None] * nChecks # list to hold if cassandra passed each check
cassStr = [None]* nChecks
cassStr[0] = (("Energy_Total","Nmols","Volume"),)
cassStr[1] = (("Energy_Total","Nmols_1","Nmols_2","Pressure", "Mass_Density"),)
cassStr[2] = (("Energy_Total","Nmols_1","Nmols_2","Nmols_3","Volume", "Pressure"),)
cassStr[3] = (("Energy_Total","Nmols_1","Nmols_2"),)
cassStr[4] = (("Energy_Total","Nmols","Density", "Pressure", "Mass_Density"),)
cassStr[5] = (("Energy_Total","Nmols_1","Nmols_2", "Pressure"),)
endStep = (4000,50000,10000,4000,4000, 4000)
errorTol = 5e-2

 # Formatting variables
bold = '\033[1m' #Will make text bold
normal = '\033[0m' #Will make the next text normal(ie. unbold)

#*******************************************************************************
# FILE MANAGEMENT
#*******************************************************************************
testSuiteFolder = os.getcwd()
MainDir 	= testSuiteFolder[0:len(testSuiteFolder)-len('Scripts/testSuite')]
cassDir 	= MainDir + "Src/"
resourceDir = MainDir + "Scripts/testSuite/Resources/"
cassExe     = args.cassandra_exe
if args.absPath:
	cassDir = ""
resFolder	= (resourceDir+"exampleResources/GCMC/Methane/",
		resourceDir+"exampleResources/GCMC/Methane_Butane/",
		resourceDir+"exampleResources/GCMC/Methane_Butane_Silicalite/",
		resourceDir+"exampleResources/GCMC/Methane_Silicalite/",
		resourceDir+"exampleResources/GCMC/Nitrogen/",
		resourceDir+"exampleResources/GCMC/Nitrogen_Silicalite/")
cassRun 	= (resFolder[0]+"methane.out",resFolder[1]+"methane_butane.out",
		resFolder[2]+"methane_butane_Si.out",resFolder[3]+"methane_Si.out",
		resFolder[4]+"nitrogen.out",resFolder[5]+"nitrogen_Si.out")
inpName 	= ("methane.inp","methane_butane.inp","methane_butane_Si.inp",
		"methane_Si.out.inp","nitrogen.inp","nitrogem_Si.out")
resultFolder= (resourceDir+"exampleResults/GCMC/Methane/",
		resourceDir+"exampleResults/GCMC/Methane_Butane/",
		resourceDir+"exampleResults/GCMC/Methane_Butane_Silicalite/",
		resourceDir+"exampleResults/GCMC/Methane_Silicalite/",
		resourceDir+"exampleResults/GCMC/Nitrogen/",
		resourceDir+"exampleResults/GCMC/Nitrogen_Silicalite/")
prpResName  = (resultFolder[0]+"methane.out.prp",resultFolder[1]+"methane_butane.out.prp",
		resultFolder[2]+"methane_butane_Si.out.prp", resultFolder[3]+"methane_Si.out.prp", 
		resultFolder[4]+"nitrogen.out.prp", resultFolder[5]+"nitrogen_Si.out.prp")
#*******************************************************************************
# FUNCTION DEFINITIONS
#*******************************************************************************

#*******************************************************************************
# MAIN PROGRAM BEGINS HERE
#*******************************************************************************
# Step 1) Calculate the correct answer
# Step 2) Write input files
# Step 3) Run Cassandra to get its answer
# Step 4) Compare answers
#

#This prints the starting line.
print "\n\n"+bold+"Test " + str(test_no) +": " + test_desc + normal

FailCount = 0; # Counting number of failed tests


# Loop through checks
print "%-30s %-20s %18s %18s %18s %8s" % ("Title", "Property","Cassandra","Accepted","Relative_Err","Pass")

for i in range(nChecks):

	# Step 2) Write input files
	# 2.1) Write - only required to write on first check - mcf does not change
	# Step 3) Run Cassandra to get its answer

	# Change directory to the resource folder to run using given files
	os.chdir(resFolder[i])
	proc = sp.Popen([cassDir + cassExe + " " + inpName[i]], stdout=sp.PIPE, shell=True)
	(out, err) = proc.communicate()

	if err is not None:
		print("Error.Abort.")

	# Change directory back to test suite
	os.chdir(MainDir+"Scripts/testSuite/")

	nPrp = len(cassStr[i][0])

	# Convert prp files to a readable format with random access
	with open(cassRun[i] + ".prp") as Runproperties, open(prpResName[i]) as Resproperties:
		valuesRun = Runproperties.readlines()
		valuesAnswer = Resproperties.readlines()
		lines = len(valuesRun)
		cassrunResults=[None]*lines
		casspriorResults=[None]*lines
		for j,a in enumerate(valuesRun):
			cassrunResults[j] = a.split()
		for j,b in enumerate(valuesAnswer):
			casspriorResults[j] = b.split()

	# Determine number of rows until desired checking point
	nrows = len(cassrunResults)
	end = 110
	end1=end
	for rows in range(nrows):
		if ((rows > 2) and (end>rows) and (end==end1)) :
			if int(cassrunResults[rows][0]) == endStep[i]:
				end = rows

	# Exract answers for needed properties and desired check
	index=0;
	cassAnswer[i] = [None]*(nPrp)
	for c in range(len(cassrunResults[1])):
		if(cassrunResults[1][c] == cassStr[i][0][index]):
			cassAnswer[i][index] = float(cassrunResults[end][c-1])
			index = index+1
	index=0;
	analyticAnswer[i] = [None]*(nPrp)
	for c in range(len(casspriorResults[1])):
		if(casspriorResults[1][c] == cassStr[i][0][index]):
			analyticAnswer[i][index] = float(casspriorResults[end][c-1])
			index = index+1

	for j in range(nPrp):
		if (analyticAnswer[i][j] == 0.):
			passCheck = cassAnswer[i][j] == 0.
			if passCheck ==0:
				FailCount = FailCount+1;
			if (j == 0):
				print "%-30s %-20s %18.6g %18.6g %18s %8s" % (title[i],cassStr[i][0][j],
						cassAnswer[i][j],analyticAnswer[i][j],'',passCheck)
			else: 
				print "%-30s %-21s %17.6g %18.6g %18s %8s" % ('',cassStr[i][0][j],
						cassAnswer[i][j],analyticAnswer[i][j],'',passCheck)

		else:
			errorRel = abs(cassAnswer[i][j] - analyticAnswer[i][j])/analyticAnswer[i][j]
			passCheck = abs(errorRel) <= errorTol
			if passCheck ==0:
				FailCount = FailCount+1;
			if (j == 0):
				print "%-30s %-20s %18.6g %18.6g %18.6g %8s" % (title[i],cassStr[i][0][j],
						cassAnswer[i][j],analyticAnswer[i][j],errorRel,passCheck)
			else: 
				print "%-30s %-21s %17.6g %18.6g %18.6g %8s" % ('',cassStr[i][0][j],
						cassAnswer[i][j],analyticAnswer[i][j],errorRel,passCheck)

if (FailCount != 0):
	PassState = "False"
else:
	PassState = "True"

LastTest = open(MainDir+ 'Scripts/testSuite/testOutput/LastTest.txt',"w")
LastTest.write(PassState)
LastTest.close()

#*******************************************************************************
# DATA OUTPUT
#*******************************************************************************

#*******************************************************************************
# CLEAN UP SCRATCH FILES
#*******************************************************************************
#os.system('rm ' + inpName)
#os.system('rm ' + xyzName)
#os.system('rm ' + ' '.join(mcfName))
#os.system('rm ' + cassRun + '*')








