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
test_no = 10
test_desc = "NVT Examples"

# Check parameters
nChecks = 6 # number of simulations to run
title = ("2,2-dimethylhexane", "Diethylether - MIE", "Diethylether", "Water SPC", "Water TIP4P/2005", "Water TIP5P")
analyticAnswer = [None] * nChecks # list to hold analytic answers
cassAnswer     = [None] * nChecks # list to hold cassandra's answers
passCheck      = [None] * nChecks # list to hold if cassandra passed each check
cassStr = [None]* nChecks
cassPrint = [None]* nChecks
cassStr[0] = (("Energy_Total","Pressure"),)
cassStr[1] = (("Energy_Total","Pressure","Volume","Density"),)
cassStr[2] = (("Energy_Total","Pressure"),)
cassStr[3] = (("Energy_Total","Pressure"),)
cassStr[4] = (("Energy_Total","Pressure"),)
cassStr[5] = (("Energy_Total","Pressure"),)

cassPrint[0] = (("Energy_Total [kJ/mol-Ext]","Pressure [bar]"),)
cassPrint[1] = (("Energy_Total [kJ/mol-Ext]","Pressure [bar]","Volume [A^3]","Density [molec/A^3]"),)
cassPrint[2] = (("Energy_Total [kJ/mol-Ext]","Pressure [bar]"),)
cassPrint[3] = (("Energy_Total [kJ/mol-Ext]","Pressure [bar]"),)
cassPrint[4] = (("Energy_Total [kJ/mol-Ext]","Pressure [bar]"),)
cassPrint[5] = (("Energy_Total [kJ/mol-Ext]","Pressure [bar]"),)
endStep = (100,2200,100,100,100,100)
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
resFolder	= (resourceDir+"exampleResources/NVT/2,2-dimethylhexane/",
		resourceDir+"exampleResources/NVT/diethylether_mie/",
		resourceDir+"exampleResources/NVT/dimethylether/",
		resourceDir+"exampleResources/NVT/water_spc/",
		resourceDir+"exampleResources/NVT/water_tip4p2005/",
		resourceDir+"exampleResources/NVT/water_tip5p/")
cassRun 	= ("nvt.out","nvt_mie.out","nvt.out","nvt.out","nvt.out","nvt.out")
inpName 	= ("nvt.inp","nvt_mie.inp","nvt.inp","nvt.inp","nvt.inp","nvt.inp")
mcfRun 	= ("dmh.mcf","dee.mcf","dme.mcf","spc.mcf","tip4p.mcf","tip5p.mcf")
species = ([1],[1],[1],[1],[1],[1])
xyzFlag 	= (0,0,1,0,1,0)
xyzName 	=('','','nvt.inp.xyz','','nvt.inp.xyz','')
chkFlag 	= (1,0,0,1,0,0)
chkName 	=('nvt.inp.chk','','','nvt.inp.chk','','')


resultFolder= (resourceDir+"exampleResults/NVT/2,2-dimethylhexane/",
		resourceDir+"exampleResults/NVT/diethylether_mie/",
		resourceDir+"exampleResults/NVT/dimethylether/",
		resourceDir+"exampleResults/NVT/water_spc/",
		resourceDir+"exampleResults/NVT/water_tip4p2005/",
		resourceDir+"exampleResults/NVT/water_tip5p/")
prpResName  = (resultFolder[0]+"nvt.out.prp", resultFolder[1]+"nvt_mie.out.prp",
		resultFolder[2]+"nvt.out.prp",resultFolder[3]+"nvt.out.prp",
		resultFolder[4]+"nvt.out.prp",resultFolder[5]+"nvt.out.prp")
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

	proc = sp.Popen([cassExe + " " + inpName[i]], stdout=sp.PIPE, shell=True)
	(out, err) = proc.communicate()
	if err is not None:
		print("Error.Abort.")

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
	end = 15
	for rows in range(nrows):
		if ((rows > 2) and (end>rows)) :
			if int(cassrunResults[rows][0]) == endStep[i]:
				end = rows

	# Exract answers for needed properties
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
				print("%-30s %-35s %18.6g %18.6g %18s %8s" % (title[i],cassPrint[i][0][j],
						cassAnswer[i][j],analyticAnswer[i][j],'',passCheck))
			else: 
				print("%-30s %-36s %17.6g %18.6g %18s %8s" % ('',cassPrint[i][0][j],
						cassAnswer[i][j],analyticAnswer[i][j],'',passCheck))

		else:
			errorRel = abs(cassAnswer[i][j] - analyticAnswer[i][j])/analyticAnswer[i][j]
			passCheck = abs(errorRel) <= errorTol
			if passCheck ==0:
				FailCount = FailCount+1;
			if (j == 0):
				print("%-30s %-35s %18.6g %18.6g %18.6g %8s" % (title[i],cassPrint[i][0][j],
						cassAnswer[i][j],analyticAnswer[i][j],errorRel,passCheck))
			else: 
				print("%-30s %-36s %17.6g %18.6g %18.6g %8s" % ('',cassPrint[i][0][j],
						cassAnswer[i][j],analyticAnswer[i][j],errorRel,passCheck))

	os.system('rm '+inpName[i])
	os.system('rm '+mcfRun[i])
	for iSpec, jSpec in enumerate(species[i]):
		os.system('rm -r species'+str(jSpec))
	os.system('rm '+ cassRun[i]+ "*")
	if xyzFlag[i]:
		os.system('rm '+xyzName[i])
	if chkFlag[i]:
		os.system('rm '+chkName[i])

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








