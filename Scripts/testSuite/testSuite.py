#!/usr/bin/env python 
#*******************************************************************************
# SCRIPT:  testSuite.py
# VERSION: 2.0
# FEATURES: Test the results of all given tests scripts
#*******************************************************************************


#*******************************************************************************
# IMPORT MODULES
#*******************************************************************************
from __future__ import print_function
import subprocess
import os,sys
from testSuiteFunctions import checkLastTest
import argparse
import warnings

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

To run the test suite using a Cassandra executable elsewhere:

	> python testSuite.py /home/applications/cassandra.exe 
	> python testSuite.py /home/applications/cassandra_gfortran.exe

""")
parser.add_argument('cassandra_exe', 
                help="Cassandra executable, including path. To call an executable in the same"+
                "Cassandra package's Src folder, utilize ../../Src/cassandra.exe as the executable path.")

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
	subprocess.call(["python", j, cassExe])
	Passed += checkLastTest(overwrite,TestOutput)

print(" ")
print("---------------------------------------------------------------------------------------------------------------------")
print(" ")
print(("Passed: " + str(Passed)))
print(("Failed: " + str(len(tests)-Passed)))

if Passed ==0:

	print("Significant error. Please check path to executable.")
	# Attempt to clean up scratch files
	os.system('rm -r species*')
	os.system('rm -r Init_Config')
	os.system('rm *.mcf')
	os.system('rm *.xyz')
	os.system('rm *.out*')
	os.system('rm *.inp')
	os.system('rm *.chk')
	os.system('rm *.mcf')


# Python2 deprecation warning
if (sys.version_info < (3,0)):
    warnings.showwarning("\nSupport for Python2 is deprecated in "
        "Cassandra and will be removed in a future release. Please "
        "consider switching to Python3.", DeprecationWarning,
        'testSuite.py', 97)

# Exit codes to signal pass/fail for AZP
if len(tests) - Passed > 0:
    exit(1)
else:
    exit(0)

