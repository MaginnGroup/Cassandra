#!/usr/bin/env python 
#*******************************************************************************
# SCRIPT:  gatherExamples.py
# VERSION: 1.0
# FEATURES: Gather all example .prps into the test folders 
#*******************************************************************************

#*******************************************************************************
# IMPORT MODULES
#*******************************************************************************
from __future__ import print_function
import subprocess
import os,sys
from testSuiteFunctions import checkLastTest
import argparse

#*******************************************************************************
# ARGUMENT PARSE
#*******************************************************************************

parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=
"""DESCRIPTION:
Path to the "Examples" directory that is desired to pull from. Typically located in the main Cassandra directory.

EXAMPLES:
To pull Examples from same Cassandra root as the testSuite

	> python testSuite.py ../../Examples

To pull Examples from another Cassandra root

	> python Test#_Description.py /home/applications/Cassandra/Examples

""")
parser.add_argument('example_folder', 
                help="Example folder to pull all .prps from.", type=bool)

args = parser.parse_args()

#*******************************************************************************
# PULL EXAMPLE RESULTS
#*******************************************************************************

print("\n- Copying .prp example results...\n")
for root, dirs, files in os.walk(args.example_folder, topdown=False):
	for name in files:
   		if name.split(".")[-1]=="prp":
			print((os.path.join(root, name)))
			os.system("cp "+ os.path.join(root, name)+ " " +"Resources/exampleResults/"+root.split("/")[-2]+"/"+root.split("/")[-1]+"/"+name)




