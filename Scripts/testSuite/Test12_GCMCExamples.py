#!/usr/bin/env python 
#*******************************************************************************
# SCRIPT:  Test9_NPTexamples.py
# VERSION: 1.0
# FEATURES: Test the results of a short simulation of the NPT examples
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
test_no = 12
test_desc = "GCMC Examples"

# Check parameters
nChecks = 6 # number of simulations to run
title = ("Methane", "Methane_Butane", "Methane_Butane_Silicalite", "Methane_Silicalite", "Nitrogen", "Nitrogen_Silicalite")
analyticAnswer = [None] * nChecks # list to hold analytic answers
cassAnswer     = [None] * nChecks # list to hold cassandra's answers
passCheck      = [None] * nChecks # list to hold if cassandra passed each check
cassStr = [None]* nChecks
cassPrint = [None]* nChecks
cassStr[0] = (("Energy_Total","Nmols","Volume"),)
cassStr[1] = (("Energy_Total","Nmols_1","Nmols_2","Pressure", "Mass_Density"),)
cassStr[2] = (("Energy_Total","Nmols_1","Nmols_2","Nmols_3","Volume", "Pressure"),)
cassStr[3] = (("Energy_Total","Nmols_1","Nmols_2"),)
cassStr[4] = (("Energy_Total","Nmols","Density", "Pressure", "Mass_Density"),)
cassStr[5] = (("Energy_Total","Nmols_1","Nmols_2", "Pressure"),)

cassPrint[0] = (("Energy_Total [kJ/mol-Ext]","Nmols","Volume [A^3]"),)
cassPrint[1] = (("Energy_Total [kJ/mol-Ext]","Nmols_1","Nmols_2","Pressure [bar]", "Mass_Density [kg/m^3]"),)
cassPrint[2] = (("Energy_Total [kJ/mol-Ext]","Nmols_1","Nmols_2","Nmols_3","Volume [A^3]", "Pressure [bar]"),)
cassPrint[3] = (("Energy_Total [kJ/mol-Ext]","Nmols_1","Nmols_2"),)
cassPrint[4] = (("Energy_Total [kJ/mol-Ext]","Nmols","Density [molec/A^3]", "Pressure [bar]", "Mass_Density [kg/m^3]"),)
cassPrint[5] = (("Energy_Total [kJ/mol-Ext]","Nmols_1","Nmols_2", "Pressure [bar]"),)
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
resourceDir = MainDir + "Scripts/testSuite/Resources/"
cassExe     = args.cassandra_exe
resFolder	= (resourceDir+"exampleResources/GCMC/Methane/",
		resourceDir+"exampleResources/GCMC/Methane_Butane/",
		resourceDir+"exampleResources/GCMC/Methane_Butane_Silicalite/",
		resourceDir+"exampleResources/GCMC/Methane_Silicalite/",
		resourceDir+"exampleResources/GCMC/Nitrogen/",
		resourceDir+"exampleResources/GCMC/Nitrogen_Silicalite/")
cassRun 	= ("methane.out","methane_butane.out","methane_butane_Si.out",
		"methane_Si.out","nitrogen.out","nitrogen_Si.out")
inpName 	= ("methane.inp","methane_butane.inp","methane_butane_Si.inp",
		"methane_Si.inp","nitrogen.inp","nitrogen_Si.inp")
mcfRun 	= (["CH4.mcf"],["C4H10.mcf","CH4.mcf"],["C4H10.mcf","CH4.mcf","unitcell.mcf"],["CH4.mcf","SiO2.mcf"],
		["N23S.mcf"],["MFI.mcf","N23S.mcf"])
species = ([1],[1,2],[2,3],[2],[1],[2])
xyzFlag 	= (0,0,1,1,0,1)
xyzName 	=('','','unitcell.xyz','Si27ucEM.xyz','','MFI.xyz')
chkFlag 	= (0,0,0,0,0,0)
chkName 	=('','','','','','')
initFlag 	= (0,0,0,0,0,0)
initName 	=('','','','','','')


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
print("\n\n"+bold+"Test " + str(test_no) +": " + test_desc + normal)

FailCount = 0; # Counting number of failed tests


# Loop through checks
print("%-30s %-35s %18s %18s %18s %8s" % ("Title", "Property","Cassandra","Accepted","Relative_Err","Pass"))

for i in range(nChecks):

	# Step 2) Write input files
	# 2.1) Write - only required to write on first check - mcf does not change
	# Step 3) Run Cassandra to get its answer

	os.system('cp '+resFolder[i]+inpName[i]+' '+inpName[i])
	for iMCF, jMCF in enumerate(mcfRun[i]):
		os.system('cp '+resFolder[i]+jMCF+' '+jMCF)
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
	for iMCF, jMCF in enumerate(mcfRun[i]):
		os.system('rm '+jMCF)
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









