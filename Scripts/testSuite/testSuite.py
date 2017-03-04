#!/usr/bin/env python 
#*******************************************************************************
# SCRIPT:  testSuite.py
# VERSION: 2.0
# FEATURES: Test the results of all given tests scripts
#*******************************************************************************


#*******************************************************************************
# IMPORT MODULES
#*******************************************************************************
import subprocess
import os,sys
from testSuiteFunctions import checkLastTest
import argparse

#*******************************************************************************
# ARGUMENT PARSE
#*******************************************************************************

parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=
"""DESCRIPTION:
Runs the test suite using the Cassandra executable specified.

EXAMPLES:
To run the test suite using a Cassandra executable inside of Cassandra/Src/:

	> python testSuite.py cassandra.exe
	> python testSuite.py cassandra_gfortran.exe

To run the test suite using a Cassandra executable elsewhere:

	> python testSuite.py /home/applications/cassandra.exe --absPath
	> python testSuite.py /home/applications/cassandra.exe -a


""")
parser.add_argument('cassandra_exe', 
                help="Cassandra executable [file name if inside Src/ directory or path with " +
                "indicated --absPath flag]")
parser.add_argument('--absPath','-a', action='store_true',
                help="Signals that the Cassandra executable is given as a path instead of " +
                "given as is in the Src driectory.")

args = parser.parse_args()

#*******************************************************************************
# FILE MANAGEMENT
#*******************************************************************************

testSuiteFolder = os.getcwd()
MainDir = testSuiteFolder[0:len(testSuiteFolder)-len('Scripts/testSuite')]
TestOutput = MainDir + 'Scripts/testSuite/testOutput/LastTest.txt'
cassExe     = args.cassandra_exe

#*******************************************************************************
# MAIN PROGRAM BEGINS HERE
#*******************************************************************************

Passed = 0 # Total Passed

# Clear the current test output by running the check function with overwrite value as true:
checkLastTest(1,TestOutput)

# The following are the tests in the testSuite:
tests=("Test1_LJ_StartEnergy.py", "Test2_Mie_StartEnergy.py", "Test3_Angle_StartEnergy.py", 
	"Test4_Dihedral_StartEnergy.py", "Test5_ImproperEnergy.py","Test6_NIST_LennardJonesEnergy.py",
	"Test7_Water.py","Test8_Coul_StartEnergy.py", "Test9_NPTExamples.py",
	"Test10_NVTExamples.py", "Test11_GEMCExamples.py", "Test12_GCMCExamples.py")

for i,j in enumerate(tests):
	if i==len(tests)-1:
		overwrite=0
	else:
		overwrite=1
	if args.absPath:
		subprocess.call(["python", j, cassExe, "--absPath" ])
	else:
		subprocess.call(["python", j, cassExe])
	Passed += checkLastTest(overwrite,TestOutput)

print(" ")
print("---------------------------------------------------------------------------------------------------------------------")
print(" ")
print("Passed: " + str(Passed))
print("Failed: " + str(len(tests)-Passed))
