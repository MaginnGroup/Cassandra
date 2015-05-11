#!/usr/bin/env python
import sys, math, subprocess, os, linecache, argparse, re

#********************************************************************************
#   Cassandra - An open source atomistic Monte Carlo software package
#   developed at the University of Notre Dame.
#   http://cassandra.nd.edu
#   Prof. Edward Maginn <ed@nd.edu>
#   Copyright (2013) University of Notre Dame du Lac
#
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program.  If not, see <http://www.gnu.org/licenses/>.
#********************************************************************************

##################################################################################
#
# Cassandra script for generating fragment libraries automatically
#
# This script sets up a Cassandra simulation.
#
# Files required before running:
#	MCF files for each species
#	CAR files for each species
#	Master input file
#
# Input files produced:
# 
# /species?/fragments/species?_mcf_gen.inp
# /species?/frag?/frag?.inp
#
# The script overwrites the section of the input file where needed (i.e. # Frag_Info)
# with the aforementioned folder locations. 
#
# Usage: library_setup.py /path/cassandra.exe input_file.inp pdbfilespecies1.pdb pdfilespecies2.pdb ...
#
# Example: library_setup.py /home/Cassandra/cassandra.exe input_file.inp pdb1.pdb pdb2.pdb ...
##################################################################################
####FUNCTION DEFINITIONS####

def check_type_infilename(infilename):
        file = open(infilename,'r')
        for line in file:
                if '<molecule>' in line:
                        return 'cml'
        return 'pdb'

def is_number(a):
	try:
		float(a)
	except ValueError:
		return False
	return True

def cml_to_pdb(infilename):
	file = open(infilename,'r')
	cml_atom_info=[]
	cml_bond_info=[]
	for line_nbr, line in enumerate(file):
		if '<bondArray>' in line:
			cml_start_bonds = line_nbr + 1
		if '</bondArray>' in line:
			cml_end_bonds = line_nbr + 1
		if '<atomArray>' in line:
			cml_start_atom = line_nbr + 1
		if '</atomArray>' in line:
			cml_end_atom = line_nbr + 1

	for i in xrange(cml_start_atom+1, cml_end_atom):
		cml_atom_info.append(re.findall('"([^"]*)"',linecache.getline(infilename, i)))
			
	for i in xrange(cml_start_bonds+1, cml_end_bonds):
		cml_bond_info.append(re.findall('"([^"]*)"',linecache.getline(infilename, i))[0].split())

	temp=[]
	coordinates=[]
	for i,line in enumerate(cml_atom_info):
		for j,element in enumerate(line):
			if is_number(cml_atom_info[i][j]):
				temp.append(cml_atom_info[i][j])
		coordinates.append(temp)
		temp=[]

	filepdb = open(os.path.splitext(infilename)[0]+'.pdb','w')

	for line_nbr, line in enumerate(cml_atom_info):
		filepdb.write('HETATM       '+str(line_nbr+1)+'       '+line[1]+'       '+coordinates[line_nbr][0]+'     '+coordinates[line_nbr][1]+'     '+coordinates[line_nbr][2]+'     '+line[1]+'       '+line[5]+"\n")
		
	cml_aux_vector1 = []
	cml_aux_vector2 = []
	cml_aux_vector3 = []
	cml_aux_vector4 = []
	cml_aux_vector5 = []
	cml_aux_vector6 = []
	for line_nbr, line in enumerate(cml_atom_info):
		cml_current_atom = line[0]
		for line_bond_nbr, line_bond in enumerate(cml_bond_info):
			if cml_current_atom in line_bond:
				cml_aux_vector1.append(line_bond_nbr)
		cml_aux_vector2.append([cml_current_atom, cml_aux_vector1])
		cml_aux_vector1=[]

        #cml_aux_vector2 = [atom, [lines in cml_bond_info]]		
	
	for line in cml_aux_vector2:
		for x in line[1]:
			index=cml_bond_info[x].index(line[0])
			if index==0:
				otherindex=1
			if index==1:
				otherindex=0
			cml_aux_vector3.append(cml_bond_info[x][otherindex])
		cml_aux_vector4.append([line[0], cml_aux_vector3])
		cml_aux_vector3=[]

	#cml_aux_vector3 = [atom1 atom2]

				
	for line in cml_aux_vector4:
		for subline in line[1]:
			cml_match0=re.match(r"([a-z]+)([0-9]+)", subline, re.I)
			if cml_match0:
				cml_aux_vector5.append(cml_match0.groups()[-1])
		
		cml_match = re.match(r"([a-z]+)([0-9]+)", line[0], re.I)
		if cml_match:
			cml_split = cml_match.groups()
			atom_number = cml_split[-1]
		cml_aux_vector6.append([atom_number, cml_aux_vector5])
		cml_aux_vector5=[]

	for line in cml_aux_vector6:
		filepdb.write('CONECT     '+line[0]+'    ')
		for element in line[1]:
			filepdb.write(element+"   ")
		filepdb.write("\n")

	filepdb.close()
	return os.path.splitext(infilename)[0]+'.pdb'

####END FUNCTION DEFINITIONS####

parser = argparse.ArgumentParser("This script helps to set up a simulation by creating folders containing the fragment libraries and the MCF fragment files automatically")
parser.add_argument('cassandra_path', action="store")
parser.add_argument('inputfile', action="store")
parser.add_argument('pdbfiles', nargs=argparse.REMAINDER)
results = parser.parse_args()

pdbfile = results.pdbfiles


for index,eachfile in enumerate(pdbfile):
	if check_type_infilename(eachfile)=='cml':
		this_pdbfile = cml_to_pdb(eachfile)
		pdbfile[index]=this_pdbfile

cassandra_location = os.path.abspath(results.cassandra_path)
input_file = os.path.abspath(results.inputfile)

bold = '\033[1m'
normal = '\033[0m'

gcmc_flag = 0
chempot_flag = 0
fugacity_flag = 0

zbyomega = 1

print bold+"\n*********Cassandra setup*********\n"

print "Cassandra location: "+normal+cassandra_location
print bold+"Scanning input file..."

#Locate line number for the following keywords:
# Nbr_Species
# Molecule_Files
# Start_Type
# Fragment_Files
# Temperature

in_file = open(input_file,'r')
keyword_line = []

for line_number,line in enumerate(in_file):
	if not line.strip():
		continue
	line_parse = line.split()
	if line_parse[0]=='#':		
		if line_parse[1] == "Nbr_Species":
			nbr_species_line = line_number+1
		elif line_parse[1] == "Molecule_Files":
			molec_files_line = line_number+1
		elif line_parse[1] == "Start_Type":
			start_type_line = line_number+1
                elif line_parse[1] == "Fragment_Files":
			frag_files_line = line_number+1
		elif line_parse[1] == "Temperature_Info":
			temp_line = line_number+1
		elif line_parse[1] == "Sim_Type":
			sim_type_line = line_number + 1
		elif line_parse[1] == "Charge_Style":
			charge_style_line = line_number+1
		elif line_parse[1] == "Box_Info":
			box_info_line = line_number+1
		elif line_parse[1] == "Chemical_Potential_Info":
			chempot_flag = 1
		elif line_parse[1] == "Fugacity_Info":
			fugacity_flag = 1

in_file.close()

#Obtain Nbr_Boxes

nbr_boxes = int(linecache.getline(input_file, box_info_line+1))

#Obtain Charge_Style info

charge_style = []
for i in xrange(0,nbr_boxes):
	charge_style.append(linecache.getline(input_file, charge_style_line+1+i))


#Obtain simulation type
sim_type = linecache.getline(input_file,sim_type_line+1)
if "GCMC" in sim_type:
	gcmc_flag = 1
	if chempot_flag == 1 and fugacity_flag == 1:
		print bold + "Aborting. Both Chemical potential and Fugacity were specified\n"
		print bold + "within a GCMC simulation. Select only one of these."
		sys.exit()

	if fugacity_flag == 1:
		print bold+"Grand Canonical simulation found. Will look for Zig/Omega in log files."+normal



#Obtain number of species
nbr_species = int(linecache.getline(input_file,nbr_species_line+1))
print bold+"Number of species found: " + normal + str(nbr_species)

#Look for MCF files
mcf_files = []
for i in xrange(0,nbr_species):
	mcf_files.append(linecache.getline(input_file,molec_files_line+i+1).split()[0])
	print bold+"The MCF file number " + str(i+1) + " is: " + normal + mcf_files[i]

#Open the MCF files to figure out how many fragments each species has
line_where_fragment_info = []
for i in xrange(0,len(mcf_files)):
	current_mcf = open(mcf_files[i],'r')
	for line_number_mcf, line_mcf in enumerate(current_mcf):
		
		if not line_mcf.strip():
			continue
		line_mcf_parse = line_mcf.split()
		if len(line_mcf_parse)>1:
			if line_mcf_parse[1] == "Fragment_Info":
				line_where_fragment_info.append(line_number_mcf+1)
				current_mcf.close()
				break

#Open the MCF files to figure out where the atom info section is, and then see 
# if there is any atom labeled as "ring"

line_where_atom_info = []
for i in xrange(0,len(mcf_files)):
	current_mcf = open(mcf_files[i],'r')
	for line_number_mcf, line_mcf in enumerate(current_mcf):
	
		if not line_mcf.strip():
			continue
		line_mcf_parse = line_mcf.split()
		if len(line_mcf_parse)>1:
			if line_mcf_parse[1] == "Atom_Info":
				line_where_atom_info.append(line_number_mcf+1)
				current_mcf.close()
				break

atom_type_list = []
which_atoms_are_ring = []
temp_list = []
for i in xrange(0, nbr_species):
	number_atoms_this_is = int(linecache.getline(mcf_files[i],line_where_atom_info[i]+1).split()[0])
	for j in xrange(0,number_atoms_this_is):
		atom_type_list.append(linecache.getline(mcf_files[i],line_where_atom_info[i]+2+j).split()[1]+'_s'+str(i+1))
		length = len(linecache.getline(mcf_files[i],line_where_atom_info[i]+2+j).split())
		if length == 9:
			temp_list.append(str(j+1))
	which_atoms_are_ring.append(temp_list)
	temp_list=[]

#Find how many fragments there are and get a list of them.

number_fragments = []
fragment_list=[]
for i in xrange(0, nbr_species):
	number_fragments.append(int(linecache.getline(mcf_files[i],line_where_fragment_info[i]+1).split()[0]))
	print bold+"Species "+str(i+1)+" has "+str(number_fragments[i])+ " fragments"
	temp_list=[]
	for j in xrange(0,number_fragments[i]):
		temp_list.append(linecache.getline(mcf_files[i],line_where_fragment_info[i]+2+j).split()[2:])
	fragment_list.append(temp_list)

#Based on the list of atoms labeled as "ring", and the fragment list, find out which fragment
#is a ring. This will be done by stating that if three or more atoms labeled fragment are in
#a given fragment list, then this list is a ring.

temp_list=[]
what_fragments_are_ring=[]
for index,speciesi in enumerate(fragment_list):
	for index_fragment,fragment in enumerate(speciesi):
		intersection = set.intersection(set(fragment),set(which_atoms_are_ring[index]))
		if len(intersection)>2:
			temp_list.append(str(index_fragment))
	what_fragments_are_ring.append(temp_list)
	temp_list=[]

if len(what_fragments_are_ring)>0:
	print bold+"Molecules with rings found. These are:"
	for index,element in enumerate(what_fragments_are_ring):
		print "Species "+str(index+1)+" has "+str(len(element))+ " rings."

#We know what fragments are ring for each species. Copy only those pdbfiles for the species that have rings

files_ring_to_copy=[]
temp=[]
for ispecies,fragment_list_ring in enumerate(what_fragments_are_ring):
	if len(fragment_list_ring)>0:
		filename_mcf = os.path.splitext(os.path.basename(mcf_files[ispecies]))[0]
		for index,element in enumerate(pdbfile):
			if filename_mcf in element:
				temp.append(str(ispecies))
				temp.append(pdbfile[index])
				files_ring_to_copy.append(temp)
				temp=[]
#files_ring_to_copy = [[speciesring, locationpdb]]

#Look for temperature
temperature = float(linecache.getline(input_file,temp_line+1).split()[0])


#Rewrite the MCF files to append an '_s1' to the atom type
stride_2=0
for i in xrange(0,len(mcf_files)):
        current_mcf = open(mcf_files[i],'r')
	new_mcf_file = open(mcf_files[i]+"temp",'w')
	total_lines = len(current_mcf.readlines())
	current_mcf.close()
	number_of_atoms = int(linecache.getline(mcf_files[i],line_where_atom_info[i]+1))
	stride = 0
	for line_number in xrange(1,total_lines+1):
#		counter_atom_type_list = 0
	        if line_number + stride > line_where_atom_info[i]+1 and line_number + stride < line_where_atom_info[i]+1 + number_of_atoms:
			for j in xrange(0,number_of_atoms):
				this_line = linecache.getline(mcf_files[i],line_number+j+stride).split()
				this_line[1] = atom_type_list[j+stride_2]
				new_mcf_file.write('    '.join(this_line)+'\n')
			stride = stride+number_of_atoms	
			stride_2 = stride_2+number_of_atoms
			new_mcf_file.write('\n')
		else:
			new_mcf_file.write(linecache.getline(mcf_files[i],line_number+stride))
#                	new_mcf_file.write(linecache.getline(input_file,line_number))
	new_mcf_file.close()
	os.system("mv "+mcf_files[i]+'temp '+mcf_files[i])



################




#So far, we should have enough information to set up the whole simulation.
#First, let's create N number of folders, where N is the number of species
#Inside each folder, we'll include the following folders:
#/species?/fragments
#/species?/car
#/species?/frag?

for i in xrange(0, nbr_species):
	os.system("mkdir species"+str(i+1))
	os.system("mkdir species"+str(i+1)+"/fragments/")
	for j in xrange(0,number_fragments[i]):
		os.system("mkdir species"+str(i+1)+"/frag"+str(j+1))

#Now, create input files for fragment MCF generation
for i in xrange(0, nbr_species):
	for element in files_ring_to_copy:
		if str(i) in element[0]:
			os.system("cp "+element[1]+" species"+str(i+1)+"/fragments/molecule.pdb")

	print "\n\n"+bold+"Creating input MCF generation file for species "+str(i+1) +"..."+normal
	input_mcf_gen = open("species"+str(i+1)+"/fragments/species"+str(i+1)+"_mcf_gen.inp",'w')
	input_mcf_gen.write("# Run_Name\nspecies"+str(i+1)+"mcf_gen")
	input_mcf_gen.write("\n\n")
	input_mcf_gen.write("# Sim_Type\nMCF_Gen\n\n")
	input_mcf_gen.write("# Nbr_Species\n1\n\n")
	input_mcf_gen.write("# VDW_Style\nLJ cut_tail 12.0\n\n")
	input_mcf_gen.write("# Rcutoff_Low\n1.0\n\n")
	input_mcf_gen.write("# Mixing_Rule\nLB\n\n")
	input_mcf_gen.write("# Charge_Style\ncoul Ewald 12.0 0.000001\n\n")
	input_mcf_gen.write("# Intra_Scaling\n0.0 0.0 0.0 1.0\n0.0 0.0 0.0 1.0\n\n")
	input_mcf_gen.write("# Molecule_Files\n../../"+mcf_files[i]+" 50\n\n")
	input_mcf_gen.write("# Box_Info\n1\nCUBIC\n10.0 10.0 10.0\n\nEND")
	print "Done..."
	input_mcf_gen.close()
	print "Running Cassandra to generate MCF files..."+normal
	os.chdir('./species'+str(i+1)+'/fragments/')
	subprocess.call([cassandra_location,'species'+str(i+1)+'_mcf_gen.inp'])
	os.chdir('../../')

print bold+"Done..."

#Now, create fragment library generation files for each fragment of each species

for i in xrange(0, nbr_species):
	for j in xrange(0,number_fragments[i]):

		if str(j) in what_fragments_are_ring[i]:
			print bold+"Generating RING FRAGMENT library species "+str(i+1)+" fragment "+str(j+1)+normal
			input_frag = open("species"+str(i+1)+"/frag"+str(j+1)+"/frag"+str(j+1)+".inp",'w')
			input_frag.write("# Run_Name\nfrag"+str(j+1)+"\n\n")
			input_frag.write("# Sim_Type\nNVT_MC_Ring_Fragment\n\n")
			input_frag.write("# Nbr_Species\n1\n\n")
			input_frag.write("# VDW_Style\nLJ cut 14.0\n\n")
			input_frag.write("# Mixing_Rule\nLB\n\n")
			input_frag.write("# Rcutoff_Low\n1.0\n\n")


			for index,box_charge_style in enumerate(charge_style):
				if 'NONE' not in box_charge_style:
					input_frag.write("# Charge_Style\n"+charge_style[index]+"\n")
					break
				else:
					if index+1 == nbr_boxes:
						input_frag.write("# Charge_Style\nNONE\n\n")
						break


			input_frag.write("# Intra_Scaling\n0.0 0.0 0.0 1.0\n0.0 0.0 0.0 1.0\n\n")
			input_frag.write("# Box_Info\n1\nCUBIC\n50.0 50.0 50.0\n\n")
			input_frag.write("# Temperature_Info\n"+str(temperature)+"\n\n")
			input_frag.write("# Seed_Info\n706111630 70611631\n\n")
			input_frag.write("# Molecule_Files\n../fragments/frag_"+str(j+1)+"_1.mcf 1\n\n")
			input_frag.write("# Start_Type\nread_old\n../fragments/frag_"+str(j+1)+"_1.xyz 1\n\n")
			input_frag.write("# Move_Probability_Info\n# Prob_Ring\n1.0 35.0\n# Done_Probability_Info\n\n")
			input_frag.write("# Run_Type\nProduction 1000 10\n\n")
			input_frag.write("# Frequency_Info\nfreq_type    none\nNthermofreq  100\nNcoordfreq   5000\nMCsteps      1000000\n# Done_Frequency_Info\n\n")
			input_frag.write("# Property_Info 1\nEnergy_Total\n\n")
			input_frag.write("# File_Info\nfrag"+str(j+1)+".dat\n\n")
			input_frag.write("END")
			input_frag.close()
		        os.chdir('./species'+str(i+1)+'/frag'+str(j+1)+'/')
		        subprocess.call([cassandra_location,'frag'+str(j+1)+'.inp'])
			if gcmc_flag == 1:
				logfile = open('frag'+str(j+1)+'.log','r')
				for line in logfile:
					if 'Z/Omega' in line:
						zbyomega = zbyomega * float(line.split()[1])
				logfile.close()
			os.chdir('../../')
		else:
			print bold+"Generating fragment library species "+str(i+1)+" fragment "+str(j+1)+normal
			input_frag = open("species"+str(i+1)+"/frag"+str(j+1)+"/frag"+str(j+1)+".inp",'w')
			input_frag.write("# Run_Name\nfrag"+str(j+1)+"\n\n")
			input_frag.write("# Sim_Type\nNVT_MC_Fragment\n\n")
			input_frag.write("# Nbr_Species\n1\n\n")
			input_frag.write("# VDW_Style\nNONE\n\n")
			input_frag.write("# Rcutoff_Low\n0.0\n\n")
			input_frag.write("# Mixing_Rule\nLB\n\n")
			input_frag.write("# Charge_Style\nNONE\n\n")
			input_frag.write("# Intra_Scaling\n0.0 0.0 0.0 1.0\n\n")
			input_frag.write("# Molecule_Files\n../fragments/frag_"+str(j+1)+"_1.mcf 1\n\n")
			input_frag.write("# Box_Info\n1\nCUBIC\n20.0 20.0 20.0\n\n")
			input_frag.write("# Temperature_Info\n"+str(temperature)+"\n\n")
			input_frag.write("# Seed_Info\n706111630 70611631\n\n")
			input_frag.write("# Move_Probability_Info\n# Prob_Translation\n1.0\n0.2 10.0\n1.0\n# Done_Probability_Info\n\n")
			input_frag.write("# Start_Type\nread_old\n../fragments/frag_"+str(j+1)+"_1.xyz 1\n\n")
			input_frag.write("# Run_Type\nProduction 1000 10\n\n")
			input_frag.write("# Frequency_Info\nfreq_type    none\nNthermofreq  10\nNcoordfreq   90\nMCsteps      1100000\nNequilSteps  100000\n# Done_Frequency_Info\n\n")
			input_frag.write("# File_Info\nfrag"+str(j+1)+".dat\n\n")
			input_frag.write("END")
			input_frag.close()
#			raw_input("here"+str(i)+str(j))
		        os.chdir('./species'+str(i+1)+'/frag'+str(j+1)+'/')
		        subprocess.call([cassandra_location,'frag'+str(j+1)+'.inp'])
			if gcmc_flag == 1:
				logfile = open('frag'+str(j+1)+'.log','r')
				if fugacity_flag == 1:
					for line in logfile:
						if 'Z/Omega' in line:
							zbyomega = zbyomega * float(line.split()[1])
					logfile.close()
				os.chdir('../../')


#Go back to the master input file and rewrite the location of the fragment libraries.
in_file = open(input_file,'r')
new_file = open(input_file+"temp",'w')
total_lines = len(in_file.readlines())
in_file.close()
for line_number in xrange(1,total_lines):
	
	if line_number == frag_files_line:
		new_file.write("# Fragment_Files\n")
		total_frag_counter = 1
		for i in xrange(0, nbr_species):
		        for j in xrange(0,number_fragments[i]):
				new_file.write("species"+str(i+1)+"/frag"+str(j+1)+"/frag"+str(j+1)+".dat  " + str(total_frag_counter)+"\n")
				total_frag_counter = total_frag_counter + 1
			
		break			
	else:
		new_file.write(linecache.getline(input_file,line_number))

line_number = line_number + total_frag_counter

if gcmc_flag == 1 and fugacity_flag == 1:
	
	new_file.write("\n! DO NOT CHANGE THE SECTION ZIG BY OMEGA!")
	new_file.write("\n# Zig_By_Omega_Info\n")
	new_file.write(str(zbyomega)+"\n")

for line in xrange(line_number,total_lines+1):
	new_file.write(linecache.getline(input_file,line))

in_file.close()
new_file.close()

print bold+"Removing temporary input file"
os.system("rm "+ input_file+"; mv "+input_file+"temp "+input_file)
print "Finished"+normal
