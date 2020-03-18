#!/usr/bin/env python 
#*******************************************************************************
# SCRIPT:  Test10_NVTexamples.py
# VERSION: 1.0
# FEATURES: Test the results of a short simulation of the NVT examples
#*******************************************************************************

#*******************************************************************************
# IMPORT MODULES
#*******************************************************************************
from __future__ import print_function
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
Runs the test suite using the Cassandra executable specified, including path.

EXAMPLES:
To run the test suite using a Cassandra executable inside of Cassandra/Src/:

	> python testSuite.py ../../Src/cassandra.exe
	> python testSuite.py ../../Src/cassandra_gfortran.exe

To run a test using a Cassandra executable elsewhere:

	> python Test#_Description.py /home/applications/cassandra.exe 
	> python Test#_Description.py /home/applications/cassandra_gfortran.exe

""")
parser.add_argument('cassandra_exe', 
                help="Cassandra executable, including path. To call an executable in the same"+
                "Cassandra package's Src folder, utilize ../../Src/cassandra.exe as the executable path.")

args = parser.parse_args()

#*******************************************************************************
# VARIABLE DEFINITIONS
#*******************************************************************************
# Test description
test_no = 11
test_desc = "GEMC Examples"

# Check parameters
nChecks = 5 # number of simulations to run
nBox = 2 # number of simulation boxes
title = ("2,2-dimethylhexane", "Cyclohexane", "Diethylether", "Isobutane", "Methane")
analyticAnswer = [None] * nChecks # list to hold analytic answers
cassAnswer     = [None] * nChecks # list to hold cassandra's answers
passCheck      = [None] * nChecks # list to hold if cassandra passed each check
cassStr = [None]* nChecks
cassPrint = [None]* nChecks
cassStr[0] = (("Energy_Total","Density","Nmols","Volume","Mass_Density"),)
cassStr[1] = (("Energy_Total","Density","Nmols","Volume","Pressure"),)
cassStr[2] = (("Energy_Total","Density","Nmols","Volume","Pressure"),)
cassStr[3] = (("Energy_Total","Density","Nmols","Volume"),)
cassStr[4] = (("Density","Nmols","Volume","Pressure"),)

cassPrint[0] = (("Energy_Total [kJ/mol-Ext]","Density [molec/A^3]","Nmols","Volume [A^3]","Mass_Density [kg/m^3]"),)
cassPrint[1] = (("Energy_Total [kJ/mol-Ext]","Density [molec/A^3]","Nmols","Volume [A^3]","Pressure [bar]"),)
cassPrint[2] = (("Energy_Total [kJ/mol-Ext]","Density [molec/A^3]","Nmols","Volume [A^3]","Pressure [bar]"),)
cassPrint[3] = (("Energy_Total [kJ/mol-Ext]","Density [molec/A^3]","Nmols","Volume [A^3]"),)
cassPrint[4] = (("Density [molec/A^3]","Nmols","Volume [A^3]","Pressure"),)

endStep = (2200,2200,2200,2200,2200)
errorTol = 5e-4

 # Formatting variables
bold = '\033[1m' #Will make text bold
normal = '\033[0m' #Will make the next text normal(ie. unbold)

#*******************************************************************************
# FILE MANAGEMENT
#*******************************************************************************
testSuiteFolder = os.getcwd()
MainDir 	= testSuiteFolder[0:len(testSuiteFolder)-len('Scripts/testSuite')]
resourceDir = MainDir + "Scripts/testSuite/Resources/"
cassExe     = args.cassandra_exe
resFolder	= (resourceDir+"exampleResources/GEMC/2,2-dimethylhexane/",
		resourceDir+"exampleResources/GEMC/Cyclohexane/",
		resourceDir+"exampleResources/GEMC/Diethylether/",
		resourceDir+"exampleResources/GEMC/Isobutane/",
		resourceDir+"exampleResources/GEMC/Methane/")
cassRun 	= ("gemc_dimethylhexane.out","gemc_cyclohexane.out","gemc_diethylether.out","gemc_isobutane.out","gemc_methane.out")
inpName 	= ("gemc_dimethylhexane.inp","gemc_cyclohexane.inp","gemc_diethylether.inp","gemc_isobutane.inp","gemc_methane.inp")
mcfRun 	= ("dimethylhexane.mcf","cyclohexane.mcf","diethylether.mcf","isobutane.mcf","methane.mcf")
species = ([1],[1],[1],[1],[1])
xyzFlag 	= (0,0,0,0,0)
xyzName 	=('','','','','')
chkFlag 	= (0,0,0,0,0)
chkName 	=('','','','','')
initFlag 	= (0,0,1,1,0)
initName 	=('','',"Init_Config","Init_Config",'')

resultFolder= (resourceDir+"exampleResults/GEMC/2,2-dimethylhexane/",
		resourceDir+"exampleResults/GEMC/Cyclohexane/",
		resourceDir+"exampleResults/GEMC/Diethylether/",
		resourceDir+"exampleResults/GEMC/Isobutane/",
		resourceDir+"exampleResults/GEMC/Methane/")
fileResultName  = (resultFolder[0]+"gemc_dimethylhexane.out",resultFolder[1]+"gemc_cyclohexane.out",resultFolder[2]+"gemc_diethylether.out",
		resultFolder[3]+"gemc_isobutane.out",resultFolder[4]+"gemc_methane.out")
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
print("\n\n"+bold+"Test " + str(test_no) +": " + test_desc + normal)

FailCount = 0; # Counting number of failed tests


# Loop through checks
print("%-30s %-35s %18s %18s %18s %8s" % ("Title", "Property","Cassandra","Accepted","Relative_Err","Pass"))

for i in range(nChecks):

	# Step 2) Write input files
	# 2.1) Write - only required to write on first check - mcf does not change
	# Step 3) Run Cassandra to get its answer

	os.system('cp '+resFolder[i]+inpName[i]+' '+inpName[i])
	os.system('cp '+resFolder[i]+mcfRun[i]+' '+mcfRun[i])
	for iSpec, jSpec in enumerate(species[i]):
		os.system('cp -r '+resFolder[i]+'species'+str(jSpec) + ' species'+str(jSpec))
	if xyzFlag[i]:
		os.system('cp '+resFolder[i]+xyzName[i]+' '+xyzName[i])
	if chkFlag[i]:
		os.system('cp '+resFolder[i]+chkName[i]+' '+chkName[i])
	if initFlag[i]:
		os.system('cp -r '+resFolder[i]+initName[i]+' ' + initName[i])
	
	proc = sp.Popen([cassExe + " " + inpName[i]], stdout=sp.PIPE, shell=True)
	(out, err) = proc.communicate()
	if err is not None:
		print("Error.Abort.")

	nPrp = len(cassStr[i][0])

	# Convert prp files to a readable format with random access
	cassrunResults=[None]*nBox
	casspriorResults=[None]*nBox
	with open(cassRun[i] + ".box1.prp") as Runproperties1, open(fileResultName[i] + ".box1.prp") as Resproperties1:
		valuesRun1 = Runproperties1.readlines()
		valuesAnswer1 = Resproperties1.readlines()
		for d in range(nBox):
			lines = len(valuesRun1)
			cassrunResults[d]=[None]*lines
			casspriorResults[d]=[None]*lines
		for j,a in enumerate(valuesRun1):
			cassrunResults[0][j] = a.split()
		for j,b in enumerate(valuesAnswer1):
			casspriorResults[0][j] = b.split()
	with open(cassRun[i] + ".box2.prp") as Runproperties2, open(fileResultName[i] + ".box2.prp") as Resproperties2:
		valuesRun2 = Runproperties2.readlines()
		valuesAnswer2 = Resproperties2.readlines()
		for j,a in enumerate(valuesRun2):
			cassrunResults[1][j] = a.split()
		for j,b in enumerate(valuesAnswer2):
			casspriorResults[1][j] = b.split()

	# Determine number of rows until desired checking point
	nrows = len(cassrunResults[0])
	end = 25
	for rows in range(nrows):
		if ((rows > 2) and (end>rows)) :
			if int(cassrunResults[0][rows][0]) == endStep[i]:
				end = rows

	# Exract answers for needed properties
	index=0;
	cassAnswer[i] = [None]*nBox
	for d in range(nBox):
		cassAnswer[i][d] = [None]*(nPrp)
	for c in range(len(cassrunResults[0][1])):
		if(cassrunResults[0][1][c] == cassStr[i][0][index]):
			cassAnswer[i][0][index] = float(cassrunResults[0][end][c-1])
			cassAnswer[i][1][index] = float(cassrunResults[1][end][c-1])
			index = index+1
	index=0;
	analyticAnswer[i] = [None]*nBox
	for d in range(nBox):
		analyticAnswer[i][d] = [None]*(nPrp)
	for c in range(len(casspriorResults[0][1])):
		if(casspriorResults[0][1][c] == cassStr[i][0][index]):
			analyticAnswer[i][0][index] = float(casspriorResults[0][end][c-1])
			analyticAnswer[i][1][index] = float(casspriorResults[1][end][c-1])

			index = index+1

	# Compare answers for each property
	
	for b in range(nBox):
		for j in range(nPrp):
			if (analyticAnswer[i][b][j] == 0.):
				passCheck = cassAnswer[i][b][j] == 0.
				if passCheck ==0:
					FailCount = FailCount+1;
				if (j == 0):
					print("%-30s %-35s %18.6g %18.6g %18s %8s" % (title[i]+" Box "+str(b+1),cassPrint[i][0][j],
							cassAnswer[i][b][j],analyticAnswer[i][b][j],'',passCheck))
				else: 
					print("%-30s %-36s %17.6g %18.6g %18s %8s" % ('',cassPrint[i][0][j],
							cassAnswer[i][b][j],analyticAnswer[i][b][j],'',passCheck))

			else:
				errorRel = abs(cassAnswer[i][b][j] - analyticAnswer[i][b][j])/analyticAnswer[i][b][j]
				passCheck = abs(errorRel) <= errorTol
				if passCheck ==0:
					FailCount = FailCount+1;
				if (j == 0):
					print("%-30s %-35s %18.6g %18.6g %18.6g %8s" % (title[i]+" Box "+str(b+1),cassPrint[i][0][j],
							cassAnswer[i][b][j],analyticAnswer[i][b][j],errorRel,passCheck))
				else: 
					print("%-30s %-36s %17.6g %18.6g %18.6g %8s" % ('',cassPrint[i][0][j],
							cassAnswer[i][b][j],analyticAnswer[i][b][j],errorRel,passCheck))

	os.system('rm '+inpName[i])
	os.system('rm '+mcfRun[i])
	for iSpec, jSpec in enumerate(species[i]):
		os.system('rm -r species'+str(jSpec))
	os.system('rm '+ cassRun[i]+ "*")
	if xyzFlag[i]:
		os.system('rm '+xyzName[i])
	if chkFlag[i]:
		os.system('rm '+chkName[i])
	if initFlag[i]:
		os.system('rm -r '+initName[i])

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
#for i in range(nChecks):
#	os.system('rm '+ "\""+ cassRun[i]+ ".prp\"")
#	os.system('rm '+ "\""+cassRun[i] + ".box1.prp\"")
#	os.system('rm '+ "\""+cassRun[i] + ".box2.prp\"")







