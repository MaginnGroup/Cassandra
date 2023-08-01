#!/usr/bin/env python 
#*******************************************************************************
# SCRIPT:  Test13_Widom.py
# VERSION: 1.0
# FEATURES: Test the results of a short Widom insertion simulation
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
test_no = 13
test_desc = "Widom Insertion & Trajectory Reader Examples"

num_threads = os.environ.get('OMP_NUM_THREADS')
if num_threads is not None:
    os.environ['OMP_NUM_THREADS'] = '1'

# Check parameters
nChecks = 3 # number of simulations to run
title = ("Diethylether", "Pentane", "Water SPCE")
analyticAnswer = [None] * nChecks # list to hold analytic answers
cassAnswer     = [None] * nChecks # list to hold cassandra's answers
passCheck      = [None] * nChecks # list to hold if cassandra passed each check
cassStr = [None]* nChecks
cassPrint = [None]* nChecks
cassStr[0] = (("Energy_Total","Pressure","Volume","Density"),)
cassStr[1] = (("Energy_Total","Pressure","Volume","Density", "Mass_Density"),)
cassStr[2] = (("Energy_Total","Pressure","Volume","Density", "Mass_Density"),)
cassPrint[0] = (("Energy_Total [kJ/mol-Ext]","Pressure [bar]","Volume [A^3]","Density [molec/A^3]", 
    "Shifted Chem. Potential [kJ/mol]", "Recommended Emax", "atompair rminsq index", "Avg widom_var", "Subgroup 1 <widom_var>", "Subgroup 2 <widom_var>", "Subgroup 3 <widom_var>"),)
cassPrint[1] = (("Energy_Total [kJ/mol-Ext]","Pressure [bar]","Volume [A^3]","Density [molec/A^3]", "Mass_Density [kg/m^3]", 
    "Shifted Chem. Potential [kJ/mol]", "Recommended Emax", "atompair rminsq index", "Avg widom_var", "Subgroup 1 <widom_var>", "Subgroup 2 <widom_var>", "Subgroup 3 <widom_var>"),)
cassPrint[2] = (("Energy_Total [kJ/mol-Ext]","Pressure [bar]","Volume [A^3]","Density [molec/A^3]", "Mass_Density [kg/m^3]", 
    "Shifted Chem. Potential [kJ/mol]", "Recommended Emax", "atompair rminsq index", "Avg widom_var", "Subgroup 1 <widom_var>", "Subgroup 2 <widom_var>", "Subgroup 3 <widom_var>"),)
endStep = (10,10,10)
errorTol = 5e-4

 # Formatting variables
bold = '\033[1m' #Will make text bold
normal = '\033[0m' #Will make the next text normal(ie. unbold)

#*******************************************************************************
# FILE MANAGEMENT
#*******************************************************************************
testSuiteFolder = os.getcwd()
MainDir         = testSuiteFolder[0:len(testSuiteFolder)-len('Scripts/testSuite')]
resourceDir = MainDir + "Scripts/testSuite/Resources/"
cassExe     = args.cassandra_exe
resFolder       = (resourceDir+"exampleResources/Widom/diethylether/",
                resourceDir+"exampleResources/Widom/pentane/",
                resourceDir+"exampleResources/Widom/water_spce/")
cassRun         = ("dee_widom.out","pentane_widom.out","spce_widom.out")
mcfRun  = ("dee.mcf","pentane.mcf","spce.mcf")
species = ([1],[1],[1])
inpName         = ("dee_widom.inp","pentane_widom.inp","spce_widom.inp")
xyzFlag         = (1,1,1)
hFileFlag       = (1,1,1)
xyzName         =('dee_traj.xyz','pentane_traj.xyz','spce_traj.xyz')
hFileName       =('dee_traj.H','pentane_traj.H','spce_traj.H')



resultFolder= (resourceDir+"exampleResults/Widom/diethylether/",
                resourceDir+"exampleResults/Widom/pentane/",
                resourceDir+"exampleResults/Widom/water_spce/")
prpResName  = (resultFolder[0]+"dee_widom.out.prp",resultFolder[1]+"pentane_widom.out.prp",resultFolder[2]+"spce_widom.out.prp")
logResName  = [resultFolder[i]+cassRun[i]+".log" for i in range(nChecks)]
wprpResName  = [resultFolder[i]+cassRun[i]+".spec1.wprp" for i in range(nChecks)]
wprp2ResName  = [resultFolder[i]+cassRun[i]+".spec1.wprp2" for i in range(nChecks)]
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
        if hFileFlag[i]:
            os.system('cp '+resFolder[i]+hFileName[i]+' '+hFileName[i])

        proc = sp.Popen([cassExe + " " + inpName[i]], stdout=sp.PIPE, shell=True)
        (out, err) = proc.communicate()
        if err is not None:
                print("Error.Abort.")

        nPrp = len(cassPrint[i][0])

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

        # Exract answers for needed properties and desired check
        index=0;
        cassAnswer[i] = [None]*(nPrp)
        for c in range(len(cassrunResults[1])):
                if(cassrunResults[1][c] == cassStr[i][0][index]):
                        cassAnswer[i][index] = float(cassrunResults[end][c-1])
                        index = index+1
        with open(cassRun[i] + ".log") as runLog:
            for line in runLog:
                if "chemical potential for species" in line:
                    cassAnswer[i][index] = float(line.split()[-2])
                elif "Recommended pair " in line:
                    cassAnswer[i][index] = float(line.split('=')[1].strip().split()[0])
                else:
                    continue
                index += 1

        with open(cassRun[i]+".rminsq2") as run_rminsq2:
            for line in run_rminsq2:
                lastline = line
            cassAnswer[i][index] = int(lastline.split()[-2].strip())
        index += 1
        #index += 1
        line_num = -1

        with open(cassRun[i] + ".spec1.wprp") as runwprp:
            for line in runwprp:
                lineSplit = line.split()
                if lineSplit[0] == str(endStep[i]):
                    cassAnswer[i][index] = float(lineSplit[1])
                    break
                line_num += 1

        index += 1
        with open(cassRun[i] + ".spec1.wprp2") as runwprp2:
            for iline in range(line_num):
                runwprp2.readline()
            lineSplit = runwprp2.readline().split()
            for subgroup in range(3):
                cassAnswer[i][index] = float(lineSplit[subgroup])
                index += 1

        index=0;
        analyticAnswer[i] = [None]*(nPrp)
        for c in range(len(casspriorResults[1])):
                if(casspriorResults[1][c] == cassStr[i][0][index]):
                        analyticAnswer[i][index] = float(casspriorResults[end][c-1])
                        index = index+1

        with open(logResName[i]) as resLog:
            for line in resLog:
                if "chemical potential for species" in line:
                    analyticAnswer[i][index] = float(line.split()[-2])
                elif "Recommended pair " in line:
                    analyticAnswer[i][index] = float(line.split('=')[1].strip().split()[0])
                else:
                    continue
                index += 1
        rminsqResName = logResName[i].replace("log","rminsq2")
        with open(rminsqResName) as res_rminsq2:
            for line in res_rminsq2:
                lastline = line
            analyticAnswer[i][index] = int(lastline.split()[-2].strip())
        index += 1
        line_num = -1

        with open(wprpResName[i]) as reswprp:
            for line in reswprp:
                lineSplit = line.split()
                if lineSplit[0] == str(endStep[i]):
                    analyticAnswer[i][index] = float(lineSplit[1])
                    break
                line_num += 1

        index += 1
        with open(wprp2ResName[i]) as reswprp2:
            for iline in range(line_num):
                reswprp2.readline()
            lineSplit = reswprp2.readline().split()
            for subgroup in range(3):
                analyticAnswer[i][index] = float(lineSplit[subgroup])
                index += 1

        for j in range(nPrp):
                if (analyticAnswer[i][j] == 0.):
                        passCheck = cassAnswer[i][j] == 0.
                        if passCheck == 0:
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
        if hFileFlag[i]:
            os.system('rm '+hFileName[i])


if (FailCount != 0):
        PassState = "False"
else:
        PassState = "True"
        
with open(MainDir+ 'Scripts/testSuite/testOutput/LastTest.txt',"w") as LastTest:
    LastTest.write(PassState)

if num_threads is not None:
    os.environ['OMP_NUM_THREADS'] = num_threads



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









