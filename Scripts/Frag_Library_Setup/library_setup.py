#!/usr/bin/env python

#*******************************************************************************
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
#*******************************************************************************

#*******************************************************************************
#
# Script for generating fragment libraries automatically
#
# This script sets up a Cassandra simulation.
#
# Required files:
#	  Cassandra input file
#	  PDB files for each species
#	  MCF files for each species
#
# Fragment libraries are located in:
#   species?/frag?/frag?.dat
#
# where ? stands for the species id (i.e. species1, species2 ...)
#
# Note that the script overwrites the section of the input file where needed 
# (i.e. # Frag_Info) with the aforementioned folder locations. 
#
#*******************************************************************************

#*******************************************************************************
# IMPORT MODULES
#*******************************************************************************
import sys, math, subprocess, os, linecache, argparse, re

#*******************************************************************************
# ARGUMENT PARSE
#*******************************************************************************
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=
"""DESCRIPTION:
Creates a library of fragment conformations (e.g. frag1.dat) in the following directory tree:

  species?/fragments/species?_mcf_gen.inp
                     frag_?_?.mcf
                     frag_?_?.xyz
  species?/frag?/frag?.inp
                 frag?.dat

EXAMPLES:
To generate a fragment library of a water molecule:

    > library_setup.py /home/applications/cassandra.exe nvt.inp water.pdb

  or if the cassandra executable is in your PATH:

    > library_setup.py $(which cassandra.exe) nvt.inp water.pdb

To generate fragment libraries for a mixture of propane and butane:

    > library_setup.py cassandra.exe nvt.inp c3.pdb c4.pdb
""")
parser.add_argument('cassandra_path', 
                help="path to the cassandra executable")
parser.add_argument('input_file', 
                help="Parameters for the individual fragment MC runs are " +
                     "taken from INPUTFILE. The fragment .mcf files are " +
                     "generated from the molecular .mcf files given in " +
                     "INPUTFILE. The fragment library files (e.g. " +
                     "species1/frag1/frag1.dat) will be added to the " +
                     "# Frag_Info section of INPUTFILE.")
parser.add_argument('config_files', nargs='+',
                help="CONFIGFILES must be in either .pdb or .cml format. " +
                     "A .pdb file can be generated using Gaussview, " +
                     "while .cml files can be generated using Avogadro. " +
                     "CONFIGFILES provide the starting configuration for " +
                     "each fragment that has a ring. If the molecule is " +
                     "rigid, the coordinates are simply transferred from " +
                     "CONFIGFILE to the fragment library file.")
parser.add_argument('--nConfigs', '-n', type=int, default=100000,
                help="number of configurations to write to the fragment " +
                     "library")
parser.add_argument('--noFlags', action='store_true',
                help="Suppresses the _s? flags that are appended to each " +
                     "atom type by default in the .mcf")

args = parser.parse_args()

#*******************************************************************************
# FUNCTION DEFINITIONS
#*******************************************************************************

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
				temp.append(float(cml_atom_info[i][j]))
		coordinates.append(temp)
		temp=[]

	filepdb = open(os.path.splitext(infilename)[0]+'.pdb','w')

	for line_nbr, line in enumerate(cml_atom_info):
		filepdb.write(
			'%-6s%5d  %-3s              %8.3f%8.3f%8.3f%6.2f%6.2f          %2s   %2s\n' 
			% ('HETATM', line_nbr+1, line[1], 
				coordinates[line_nbr][0], 
				coordinates[line_nbr][1], 
				coordinates[line_nbr][2], 0., 0., line[1], line[5]) )
		
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

def read_pdb(pdb_file):
	"""arguments:
	pdb_file, string = filename of a pdb file
returns:
	atom_list, list of ints = atom numbers from pdb file
	atom_coords, dict = atom parameters: name, type, element, mass
"""
	# initialize variables
	atom_coords = {}
	pdb = open(pdb_file,'r')
	atom_nbr = 1
	for line in pdb:
		# read atom info
		this_line = line.split()
		if line[0:6]=='HETATM' or line[0:4]=='ATOM': 
			x = float(line[30:38].strip())
			y = float(line[38:46].strip())
			z = float(line[46:54].strip())
			atom_coords[atom_nbr] = (x,y,z)
			atom_nbr = atom_nbr + 1

	pdb.close()
	return atom_coords

#*******************************************************************************
# FILE MANAGEMENT
#*******************************************************************************
pdb_files = []
for index,each_file in enumerate(args.config_files):
	if check_type_infilename(each_file)=='cml':
		pdb_file = cml_to_pdb(each_file)
	else:
		pdb_file = each_file
	pdb_files.append(pdb_file)

cassandra_path = os.path.abspath(args.cassandra_path)
input_file = os.path.abspath(args.input_file)

#*******************************************************************************
# MAIN PROGRAM BEGINS HERE
#*******************************************************************************
#Initialize Variables
bold = '\033[1m'
normal = '\033[0m'
gcmc_flag = 0

print bold+"\n********************** Cassandra Setup *********************\n"+normal
print bold+"Cassandra location: "+normal+cassandra_path
print bold+"Scanning input file"+normal

#Locate line number for the following keywords:
# Nbr_Species
# Molecule_Files
# Start_Type
# Fragment_Files
# Temperature

in_file = open(input_file,'r')
keyword_line = []

# initialize optional keywords
intra_scaling_line = None
charge_style_line = None
mixing_rule_line = None

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
			sim_type_line = line_number+1
		elif line_parse[1] == "VDW_Style":
			vdw_style_line = line_number+1
		elif line_parse[1] == "Charge_Style":
			charge_style_line = line_number+1
		elif line_parse[1] == "Box_Info":
			box_info_line = line_number+1
		elif line_parse[1] == "Intra_Scaling":
			intra_scaling_line = line_number+1
		elif line_parse[1] == "Mixing_Rule":
			mixing_rule_line = line_number+1

in_file.close()

#Obtain Nbr_Boxes
nbr_boxes = int(linecache.getline(input_file, box_info_line+1))

#Obtain VDW_Style and Charge_Style info
if nbr_boxes == 1:
	vdw_style = linecache.getline(input_file, vdw_style_line+1).split()[0]
	if charge_style_line:
		charge_style = linecache.getline(input_file, charge_style_line+1).split()[0]
else:
	vdw_style = []
	charge_style = []
	for i in xrange(1,nbr_boxes+1):
		vdw_style.append(linecache.getline(input_file, vdw_style_line+i).split()[0])
		if charge_style_line:
			charge_style.append(linecache.getline(input_file, charge_style_line+i).split()[0])
	if vdw_style[1] != vdw_style[0]:
		vdw_style = raw_input("VDW_Style for boxes don't match. " +
		                      "Enter the VDW_Style to use (" + vdw_style[0] + "/" + 
		                      vdw_style[i-1] + "):")
	else:
		vdw_style = vdw_style[0]
	if vdw_style != 'LJ' and vdw_style != 'lj' and \
	   vdw_style != 'Mie' and vdw_style != 'mie':
		print "Only 'LJ' and 'Mie' VDW_Styles are supported"
		quit()
	if charge_style_line:
		if charge_style[1] != charge_style[0]:
			charge_style = raw_input("Charge_Style for boxes don't match. " +
										 "Enter the Charge_Style to use (" + charge_style[0] + "/" + 
										 charge_style[i-1] + "):")
		else:
			charge_style = charge_style[0]

#Obtain simulation type
sim_type = linecache.getline(input_file,sim_type_line+1)

#Obtain number of species
nbr_species = int(linecache.getline(input_file,nbr_species_line+1))
print bold+"  Number of species found: " + normal + str(nbr_species)

#Intra_Scaling
vdw_scaling = []
charge_scaling = []
for i in xrange(1,nbr_species+1):
	if intra_scaling_line:
		if charge_style_line and charge_style == 'coul':
			vdw_scaling.append(linecache.getline(input_file,intra_scaling_line+2*i-1))
			charge_scaling.append(linecache.getline(input_file,intra_scaling_line+2*i))
		else:
			vdw_scaling.append(linecache.getline(input_file,intra_scaling_line+i))

#Mixing_Rule
if mixing_rule_line:
	mixing_rule = linecache.getline(input_file,mixing_rule_line+1)
	if mixing_rule.strip() == "custom":
		i = 1
		while True:
			if linecache.getline(input_file,mixing_rule_line+1+i).strip() and \
			   linecache.getline(input_file,mixing_rule_line+1+i)[1] != '!':
				mixing_rule = mixing_rule + linecache.getline(input_file,mixing_rule_line+1+i)
				i += 1
			else:
				break

#Open PDB Files. Is there any files without the CONECT keyword (e.g. zeolite)? If so, 
#tag it. This species will be assumed to be comprised of rigid fragment.

pdb_without_conect = []
for this_species,each_pdb in enumerate(pdb_files):
	this_pdb = open(each_pdb,'r')
	conect_found = False	
	for line in this_pdb:
		if 'CONECT' in line:
			conect_found = True
			break
	if conect_found == False:
		pdb_without_conect.append(this_species)
	

#Look for MCF files
mcf_files = []
for i in xrange(0,nbr_species):
	mcf_files.append(linecache.getline(input_file,molec_files_line+i+1).split()[0])
	print bold+"  The MCF file number " + str(i+1) + " is: " + normal + mcf_files[i]

#Open the MCF files to record where atom, bond and fragment info is
line_where_atom_info = []
line_where_bond_info = []
line_where_fragment_info = []
for i in xrange(0,len(mcf_files)):
	current_mcf = open(mcf_files[i],'r')
	for line_number_mcf, line_mcf in enumerate(current_mcf):
		if not line_mcf.strip():
			continue
		line_mcf_parse = line_mcf.split()
		if len(line_mcf_parse)>1:
			if line_mcf_parse[1] == "Atom_Info":
				line_where_atom_info.append(line_number_mcf+1)
			elif line_mcf_parse[1] == "Bond_Info":
				line_where_bond_info.append(line_number_mcf+1)
			elif line_mcf_parse[1] == "Fragment_Info":
				line_where_fragment_info.append(line_number_mcf+1)

#is there is any atom labeled as "ring"?
atom_type_list = [0]*nbr_species
which_atoms_are_ring = []
temp_list = []
nbr_atoms = []
for i in xrange(0, nbr_species):
	nbr_atoms.append(int(linecache.getline(mcf_files[i],
	                     line_where_atom_info[i]+1).split()[0]))
	atom_type_list[i] = [0]*nbr_atoms[i]
	for j in xrange(0,nbr_atoms[i]):
		atom_type = linecache.getline(mcf_files[i],
		                              line_where_atom_info[i]+2+j).split()[1]
		#remove old '_s?' flags if present
		if atom_type[-3:-1] == '_s':
			atom_type = atom_type[:-3]
		#add '_s?" flags to atom_types, unless --noFlags given
		if not args.noFlags:
			atom_type = atom_type + '_s' + str(i+1)
		atom_type_list[i][j] = atom_type
		atom_in_ring = linecache.getline(mcf_files[i],
		               line_where_atom_info[i]+2+j).split()[-1] == 'ring'
		if atom_in_ring:
			temp_list.append(str(j+1))
	which_atoms_are_ring.append(temp_list)
	temp_list=[]

#Find how many fragments there are and get a list of them.
nbr_fragments = []
fragment_list = []
for i in xrange(0, nbr_species):
	nbr_fragments.append(int(linecache.getline(mcf_files[i],
	                            line_where_fragment_info[i]+1).split()[0]))
	print bold+"    Species "+str(i+1)+" has "+str(nbr_fragments[i])+ " fragments" + normal
	temp_list=[]
	for j in xrange(0,nbr_fragments[i]):
		temp_list.append(linecache.getline(mcf_files[i],
		                 line_where_fragment_info[i]+2+j).split()[2:])
	fragment_list.append(temp_list)

#Based on the list of atoms labeled as "ring", and the fragment list, find out 
#which fragment is a ring. This will be done by testing if three or more 
#atoms labeled ring are in a given fragment list, then this list is a ring.
temp_ring = []
temp_exoring = []
fragment_has_ring = [] # Boolean value for each i,j
ring_with_exoatoms = [] # List of ring fragments that have exoatoms for each i
molecules_have_ring = False
for i,species in enumerate(fragment_list):
	for j,fragment in enumerate(species):
		intersection = set.intersection(set(fragment),
		                                set(which_atoms_are_ring[i]))
		temp_ring.append(len(intersection) > 2)
		if temp_ring[j] and len(fragment) > len(intersection):
				temp_exoring.append(j)
	fragment_has_ring.append(temp_ring)
	ring_with_exoatoms.append(temp_exoring)
	temp_ring = []
	temp_exoring = []

if any([any(has_ring) for has_ring in fragment_has_ring]):
	print bold+"Molecules with rings found. These are:" + normal
	for i,species in enumerate(fragment_has_ring):
		print bold+"Species "+str(i+1)+" has "+str(sum(species))+ " rings."+normal

#We know what fragments are ring for each species. Copy only those config_files 
#for the species that have rings
files_ring_to_copy=[]
temp=[]
for ispecies,has_ring in enumerate(fragment_has_ring):
	if has_ring:
		filename_mcf = os.path.splitext(os.path.basename(mcf_files[ispecies]))[0]
		for index,element in enumerate(pdb_files):
			if filename_mcf in element:
				temp.append(str(ispecies))
				temp.append(pdb_files[index])
				files_ring_to_copy.append(temp)
				temp=[]
#files_ring_to_copy = [[speciesring, locationpdb]]

#Look for temperature
temperature = float(linecache.getline(input_file,temp_line+1).split()[0])

#Rewrite the MCF files to append an '_s1' to the atom type
for i in xrange(0,len(mcf_files)):
	current_mcf = open(mcf_files[i],'r')
	new_mcf_file = open(mcf_files[i]+"temp",'w')
	total_lines = len(current_mcf.readlines())
	current_mcf.close()
	stride = 0
	for line_number in xrange(1,total_lines+1):
		if line_number + stride > line_where_atom_info[i]+1 and \
		   line_number + stride <= line_where_atom_info[i]+1 + nbr_atoms[i]:
			for j in xrange(0,nbr_atoms[i]):
				this_line = linecache.getline(mcf_files[i],line_number+j+stride).split()
				this_line[1] = atom_type_list[i][j]
				new_mcf_file.write('    '.join(this_line)+'\n')
			stride = stride+nbr_atoms[i]
			new_mcf_file.write('\n')
		elif line_number + stride <= total_lines+1:
			new_mcf_file.write(linecache.getline(mcf_files[i],line_number+stride))
	new_mcf_file.close()
	os.system("mv "+mcf_files[i]+'temp '+mcf_files[i])



################




#So far, we should have enough information to set up the whole simulation.
#First, let's create N number of folders, where N is the number of species
#Inside each folder, we'll include the following folders:
#/species?/fragments
#/species?/frag?

for i in xrange(0, nbr_species):

	for j in xrange(0,nbr_fragments[i]):
		os.system("mkdir -p species"+str(i+1)+"/frag"+str(j+1))

	if i not in pdb_without_conect: 
		os.system("mkdir -p species"+str(i+1)+"/fragments/")

#Now, create input files for fragment MCF generation
for i in xrange(0, nbr_species):

	if i in pdb_without_conect:
		print "\n\n" + bold + "MCF generation file not created for species "+str(i+1)+normal
		if nbr_fragments[i] == 0:
			print bold+"No fragment configuration needed."+normal
		else:
			print bold+"Fragment configuration will be taken from PDB file."+normal
		continue

	if nbr_atoms[i] >= 3:
		for element in files_ring_to_copy:
			if str(i) in element[0]:
				os.system("cp "+element[1]+" species"+str(i+1)+"/fragments/molecule.pdb")

		print "\n\n"+bold+"Creating input MCF generation file for species "+str(i+1)+" "+normal
		input_mcf_gen = open("species"+str(i+1)+"/fragments/species"+str(i+1)+"_mcf_gen.inp",'w')
		input_mcf_gen.write("# Run_Name\nspecies"+str(i+1)+"_mcf_gen")
		input_mcf_gen.write("\n\n")
		input_mcf_gen.write("# Sim_Type\nMCF_Gen\n\n")
		input_mcf_gen.write("# Nbr_Species\n1\n\n")
		input_mcf_gen.write("# VDW_Style\n"+vdw_style+" minimum_image\n\n")
		input_mcf_gen.write("# Rcutoff_Low\n1.0\n\n")
		if mixing_rule_line:
			input_mcf_gen.write("# Mixing_Rule\n"+mixing_rule+"\n\n")
		if charge_style_line:
			if charge_style == 'coul':
				input_mcf_gen.write("# Charge_Style\n"+charge_style+" minimum_image\n\n")
				if intra_scaling_line:
					input_mcf_gen.write("# Intra_Scaling\n"+
													 vdw_scaling[i]+charge_scaling[i]+"\n")
			else:
				input_mcf_gen.write("# Charge_Style\n"+charge_style+"\n\n")
				if intra_scaling_line:
					input_mcf_gen.write("# Intra_Scaling\n"+vdw_scaling[i]+"\n")
		input_mcf_gen.write("# Molecule_Files\n../../"+mcf_files[i]+" 50\n\n")
		input_mcf_gen.write("# Box_Info\n1\nCUBIC\n30.0 30.0 30.0\n\nEND")
		input_mcf_gen.close()

		print bold+"Running Cassandra to generate MCF files"+normal
		sys.stdout.flush()
		os.chdir('./species'+str(i+1)+'/fragments/')
		subprocess.call([cassandra_path,'species'+str(i+1)+'_mcf_gen.inp'])
		os.chdir('../../')


#Test if fragment and/or ring is rigid
fragment_is_rigid = [] # Boolean entry for each i,j
ring_is_rigid = [] # frag_id's for each frag that has a rigid ring
for i in xrange(0, nbr_species):

	if nbr_atoms[i] < 3:

		temp_rigid = [True]
		temp_ring_rigid = []

	elif i in pdb_without_conect:
		
		temp_rigid = [True]
		temp_ring_rigid = []

	else:
		temp_rigid = []
		temp_ring_rigid = []
		for j in xrange(0,nbr_fragments[i]):
			#Read mcf to see how many angles are fixed
			frag_atoms_are_ring = []
			angle_parms = {}
			nbr_angles_fixed = 0
			frag_mcf= open("species" + str(i+1) + "/fragments/frag_" + str(j+1) + 
										 "_1.mcf",'r')
			line = frag_mcf.readline()
			while line:
				if "# Atom_Info" in line:
					nbr_frag_atoms = int(frag_mcf.readline().split()[0])
					for k in xrange(0,nbr_frag_atoms):
						line = frag_mcf.readline()
						atom_ID = int(line.split()[0])
						if line.split()[-1] == 'ring':
							frag_atoms_are_ring.append(atom_ID)
				elif "# Angle_Info" in line:
					nbr_angles = int(frag_mcf.readline().split()[0])
					for k in xrange(0,nbr_angles):
						line = frag_mcf.readline()
						angle_ID = tuple([int(x) for x in line.split()[1:4]])
						angle_type = line.split()[4]
						angle_parms[angle_ID] = angle_type
						if angle_type == 'fixed':
							nbr_angles_fixed += 1
				line = frag_mcf.readline()
			#If all angles are rigid, fragment is rigid
			temp_rigid.append(nbr_angles == nbr_angles_fixed)

			#Also check for a rigid ring (exoatoms may be flexible)
			if fragment_has_ring[i][j]:
				if temp_rigid[j]:
					#If the whole frag is rigid, the ring must be rigid
					temp_ring_rigid.append(j)
				elif any([angle_type == 'fixed' for angle_type in angle_parms.values()]):
					#Must have some 'fixed' angles to have a rigid ring
					nbr_ring_angles_fixed = 0
					for angle_ID in angle_parms:
						ring_angle = all([a in frag_atoms_are_ring for a in angle_ID])
						if ring_angle and angle_parms[angle_ID] == 'fixed':
							nbr_ring_angles_fixed += 1
					if nbr_ring_angles_fixed == len(frag_atoms_are_ring):
						temp_ring_rigid.append(j)

	fragment_is_rigid.append(temp_rigid)
	ring_is_rigid.append(temp_ring_rigid)

#Create input file for each fragment of each species

for i in xrange(0, nbr_species):
	for j in xrange(0,nbr_fragments[i]):
		if fragment_is_rigid[i][j]:
			#Read PDB
			atom_coords = read_pdb(pdb_files[i])
			#Rotate fragment
			#Write frag.dat
			output_frag = open("species"+str(i+1)+"/frag"+str(j+1)+"/frag"+str(j+1)+
												 ".dat",'w')
			output_frag.write('1\n') # Single conformation
			output_frag.write(str(temperature) + ' 0.0\n')
			for a in fragment_list[i][j]:
				output_frag.write('%s' % atom_type_list[i][int(a)-1])
				output_frag.write(' %8.3f %8.3f %8.3f\n' % atom_coords[int(a)])
			output_frag.close()
		else:
			if fragment_has_ring[i][j]:
				print bold+"Generating RING FRAGMENT library species "+str(i+1)+\
									 " fragment "+str(j+1)+normal
				input_frag = open("species"+str(i+1)+"/frag"+str(j+1)+"/frag"+str(j+1)+
													".inp",'w')
				input_frag.write("# Run_Name\nfrag"+str(j+1)+"\n\n")
				input_frag.write("# Sim_Type\nNVT_MC_Ring_Fragment\n\n")
				input_frag.write("# Nbr_Species\n1\n\n")
				input_frag.write("# VDW_Style\n"+vdw_style+" minimum_image\n\n")
				if mixing_rule_line:
					input_frag.write("# Mixing_Rule\n"+mixing_rule+"\n\n")
				input_frag.write("# Rcutoff_Low\n1.0\n\n")
				if charge_style_line:
					if charge_style == 'coul':
						input_frag.write("# Charge_Style\n"+charge_style+" minimum_image\n\n")
						if intra_scaling_line:
							input_frag.write("# Intra_Scaling\n"+
															 vdw_scaling[i]+charge_scaling[i]+"\n")
					else:
						input_frag.write("# Charge_Style\n"+charge_style+"\n\n")	
						if intra_scaling_line:
							input_frag.write("# Intra_Scaling\n"+vdw_scaling[i]+"\n")
				input_frag.write("# Box_Info\n1\nCUBIC\n50.0 50.0 50.0\n\n")
				input_frag.write("# Temperature_Info\n"+str(temperature)+"\n\n")
				input_frag.write("# Seed_Info\n706111630 70611631\n\n")
				input_frag.write("# Molecule_Files\n../fragments/frag_"+str(j+1)+
					               "_1.mcf 1\n\n")
				input_frag.write("# Start_Type\nread_config 1 ../fragments/frag_"+str(j+1)+
					               "_1.xyz 1\n\n")
				if j in ring_with_exoatoms[i]:
					if j in ring_is_rigid[i]:
						input_frag.write("# Move_Probability_Info\n" + 
														 "# Prob_Atom_Displacement\n1.0\n0.2 10.0\n\n"
														 "# Done_Probability_Info\n\n")
					else:
						input_frag.write("# Move_Probability_Info\n" + 
														 "# Prob_Ring\n0.5 35.0\n\n" + 
														 "# Prob_Atom_Displacement\n0.5\n0.2 10.0\n\n"
														 "# Done_Probability_Info\n\n")
				else:
					input_frag.write("# Move_Probability_Info\n# Prob_Ring\n1.0 35.0\n" + 
													 "# Done_Probability_Info\n\n")
				input_frag.write("# Run_Type\nProduction 1000 10\n\n")
				input_frag.write("# Simulation_Length_Info\nUnits    Steps\n"+
					               "prop_freq  100\ncoord_freq   5000\n"+
					               "run          "+str(100*args.nConfigs)+"\n\n")
				input_frag.write("# Property_Info 1\nEnergy_Total\n\n")
				input_frag.write("# File_Info\nfrag"+str(j+1)+".dat\n\n")
				input_frag.write("END")
				input_frag.close()
				os.chdir('./species'+str(i+1)+'/frag'+str(j+1)+'/')
				subprocess.call([cassandra_path,'frag'+str(j+1)+'.inp'])
				os.chdir('../../')
			else:
				print bold+"Generating fragment library species "+str(i+1)+\
					         " fragment " + str(j+1)+normal
				input_frag = open("species"+str(i+1)+"/frag"+str(j+1)+"/frag"+str(j+1)+
					                ".inp",'w')
				input_frag.write("# Run_Name\nfrag"+str(j+1)+"\n\n")
				input_frag.write("# Sim_Type\nNVT_MC_Fragment\n\n")
				input_frag.write("# Nbr_Species\n1\n\n")
				input_frag.write("# VDW_Style\n"+vdw_style+" minimum_image\n\n")
				input_frag.write("# Rcutoff_Low\n0.0\n\n")
				if mixing_rule_line:
					input_frag.write("# Mixing_Rule\n"+mixing_rule+"\n\n")
				if charge_style_line:
					if charge_style == 'coul':
						input_frag.write("# Charge_Style\n"+charge_style+" minimum_image\n\n")
						if intra_scaling_line:
							input_frag.write("# Intra_Scaling\n"+
															 vdw_scaling[i]+charge_scaling[i]+"\n")
					else:
						input_frag.write("# Charge_Style\n"+charge_style+"\n\n")
						if intra_scaling_line:
							input_frag.write("# Intra_Scaling\n"+vdw_scaling[i]+"\n")
				input_frag.write("# Molecule_Files\n../fragments/frag_"+str(j+1)+
												 "_1.mcf 1\n\n")
				input_frag.write("# Box_Info\n1\nCUBIC\n50.0 50.0 50.0\n\n")
				input_frag.write("# Temperature_Info\n"+str(temperature)+"\n\n")
				input_frag.write("# Seed_Info\n706111630 70611631\n\n")
				input_frag.write("# Move_Probability_Info\n# Prob_Translation\n1.0\n"+
												 "0.2 10.0\n1.0\n# Done_Probability_Info\n\n")
				input_frag.write("# Start_Type\nread_config 1 ../fragments/frag_"+str(j+1)+
												 "_1.xyz 1\n\n")
				input_frag.write("# Run_Type\nProduction 1000 10\n\n")
				input_frag.write("# Simulation_Length_Info\nUnits  Steps\n"+
				                 "prop_freq  10\ncoord_freq   90\n"+
				                 "run          "+str(11*args.nConfigs)+"\n"+
				                 "nequilsteps  "+str(args.nConfigs)+"\n\n")
				input_frag.write("# File_Info\nfrag"+str(j+1)+".dat\n\n")
				input_frag.write("END")
				input_frag.close()
#				raw_input("here"+str(i)+str(j))
				os.chdir('./species'+str(i+1)+'/frag'+str(j+1)+'/')
				subprocess.call([cassandra_path,'frag'+str(j+1)+'.inp'])
				os.chdir('../../')


#Go back to the master input file and rewrite the location of the fragment libraries.
in_file = open(input_file,'r')
new_file = open(input_file+"temp",'w')
total_lines = len(in_file.readlines())
in_file.close()
omit = False
for line_number in xrange(1,total_lines+1):
	
	if line_number == frag_files_line:
		new_file.write("# Fragment_Files\n")
		total_frag_counter = 0
		for i in xrange(0, nbr_species):
			for j in xrange(0,nbr_fragments[i]):
				total_frag_counter += 1
				new_file.write("species"+str(i+1)+"/frag"+str(j+1)+"/frag"+str(j+1)+
				               ".dat  " + str(total_frag_counter)+"\n")
		new_file.write('!------------------------------------------------------' + 
	 	               '---one line per fragment\n\n')
		omit = True
	elif not omit:
		new_file.write(linecache.getline(input_file,line_number))
	elif line_number > frag_files_line:
		if linecache.getline(input_file,line_number)[0] == '#' or \
		   'END' in linecache.getline(input_file,line_number):
			omit = False
			new_file.write(linecache.getline(input_file,line_number))

in_file.close()
new_file.close()

print bold+"Removing temporary input file"+normal
os.system("rm "+ input_file+"; mv "+input_file+"temp "+input_file)
print bold+"Finished"+normal
