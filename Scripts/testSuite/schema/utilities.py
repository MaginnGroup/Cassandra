#!/usr/bin/env python
import sys
import json
import jsonschema
from jsonschema import validate
import numpy as np


def validateJsonFile(jsonInput, jsonSchema):
	jsonData = open(jsonInput,'r')
	schemaF = open(jsonSchema,'r')
	jsonInput = json.load(jsonData)
	jsonSchema = json.load(schemaF)
	jsonData.close()
	schemaF.close()
	
	try:
	        validate(jsonInput,jsonSchema)
	        sys.stdout.write("Json data structure successfully validated.\n")
	except jsonschema.exceptions.ValidationError as ve:
	        sys.stderr.write("Json data structure could not be validated.\n")
	        sys.stderr.write(str(ve) + "\n")

def json2input(jsonData, inputFile):
	inputF = open(inputFile,'w')
	
        inputF.write(''.join(["# Run_Name","\n",jsonData["runName"],"\n\n"]))
        inputF.write(''.join(["# Sim_Type","\n",jsonData["simType"],"\n\n"]))
        inputF.write(''.join(["# Nbr_Species","\n",str(jsonData["numSpecies"]),"\n\n"]))
        inputF.write(''.join(["# VDW_Style","\n",jsonData["dispersion"]["potential"]," ", \
                               jsonData["dispersion"]["cutType"]," ",
                               str(jsonData["dispersion"]["cutOff"]),"\n\n"]))

        if (jsonData["electrostatics"]["method"] == "none"):
        	inputF.write(''.join(["# Charge_Style","\n", \
                               jsonData["electrostatics"]["method"],"\n\n"]))
	else:
	        inputF.write(''.join(["# Charge_Style","\n","coul"," ", \
                               jsonData["electrostatics"]["method"]," ", \
                               str(jsonData["electrostatics"]["cutOff"])," ",
                               str(jsonData["electrostatics"]["parm1"]),"\n\n"]))

         
        inputF.write(''.join(["# Intra_Scaling","\n", \
                        ' '.join(map(str,jsonData["intraScaling"]["dispersion"])),"\n", \
                        ' '.join(map(str,jsonData["intraScaling"]["electrostatics"])),"\n\n",]))


        inputF.write(''.join(["# Mixing_Rule","\n",jsonData["mixingRule"],"\n\n"]))

        inputF.write(''.join(["# Seed_Info","\n", \
                        ' '.join(map(str,jsonData["seeds"])),"\n\n",]))

        inputF.write(''.join(["# Rcutoff_Low","\n", \
                        str(jsonData["cutOffLow"]),"\n\n",]))

	inputF.write(''.join(["# Molecule_Files","\n", \
                       ''.join(' '.join([str(a[0]),str(a[1]),'\n']) \
                        for a in zip(jsonData["mcfFiles"],map(str,jsonData["molNumber"]))),"\n\n"]))

        inputF.write("# Box_Info\n")
	nboxes = len(jsonData["box"])
        inputF.write(''.join([str(nboxes),"\n"]))
	for ibox in range(0, nboxes):
		orthogonalVectors=[]
		for index1,vector1 in enumerate(jsonData["box"][ibox]):
			for index2,vector2 in enumerate(jsonData["box"][ibox]):
				if (index1>=index2):
					continue
				if (np.dot(vector1, vector2)==0.0):
					orthogonalVectors.append(True)
				else:
					orthogonalVectors.append(False)
	
		if all(item == True for item in orthogonalVectors):
			vector1 = jsonData["box"][ibox][0]
			vector2 = jsonData["box"][ibox][1]
			vector3 = jsonData["box"][ibox][2]
			norm1 = np.linalg.norm(vector1)
			norm2 = np.linalg.norm(vector2)
 			norm3 = np.linalg.norm(vector3)
			if (norm1==norm2 and norm2==norm3):
				inputF.write("cubic\n")
				inputF.write(''.join([str(np.linalg.norm(vector1)),"\n\n"]))
		
			else: 
				inputF.write("orthorombic\n")
				inputF.write(' '.join([str(norm1),str(norm2),str(norm3),"\n\n"]))
		else:

			inputF.write("cell_matrix\n")
			vector1 = jsonData["box"][ibox][0]
			vector2 = jsonData["box"][ibox][1]
			vector3 = jsonData["box"][ibox][2]
			inputF.write(' '.join(str(a) for a in vector1))
			inputF.write('\n')
			inputF.write(' '.join(str(a) for a in vector2))
			inputF.write('\n')
			inputF.write(' '.join(str(a) for a in vector3))
			inputF.write('\n\n')

        inputF.write("# Temperature_Info\n300.0\n!---------------\n\n")
        inputF.write("# Move_Probability_Info\n\n")
        inputF.write("# Prob_Translation\n1.0\n1.00\n\n")
        inputF.write("# Done_Probability_Info\n!----------------\n\n")
        inputF.write("# Start_Type\nread_config 2 lj.xyz\n!---------------\n\n")
        inputF.write("# Run_Type\nEquilibration 100\n!---------------\n\n")
        inputF.write("# Average_Info\n1\n!---------------(0 = yes, 1 = no)\n\n")
        inputF.write("# Simulation_Length_Info\nUnits Steps\nProp_Freq 1\nCoord_Freq 1\nRun 0\n!---------------\n\n")
        inputF.write("# Property_Info 1\nEnergy_Total\nPressure\n!---------------\n\n")
        inputF.write("# Fragment_Files\n!---------------one line per fragment\n\n")
        inputF.write("# CBMC_Info\nkappa_ins 12\nkappa_rot 0\nkappa_dih 10\nrcut_cbmc 6.5\n!---------------\n\n")
        inputF.write("END")
        inputF.close()
        
        

def mie(n,m,epsilon, sigma, rij):

    factor1 = (n/(n-m))
    factor2 = (n/m)**(m/(n-m))
    factor3 = (sigma/rij)**n - (sigma/rij)**m
    return epsilon*factor1*factor2*factor3

def lj(epsilon,sigma,rij):

    factor1 = (sigma/rij)**12 - (sigma/rij)**6
    print factor1
    return epsilon * factor1

def harmonic(constant,nominal,rij):
    
    return constant*(rij-nominal)**2

