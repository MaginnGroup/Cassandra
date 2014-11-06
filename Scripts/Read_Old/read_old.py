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

#********************************************************************************
#
# This script helps the user to restart a simulation using the read_old option.
# It creates N new xyz files, where N is the number of boxes. The location of
# these files is in the current directory. These xyz files are taken from the
# last configuration of the previous simulation.
#
# The maximum displacements are taken from the CHK file.
#
# A new inputfile will be created located in the current folder.
#
# Usage: read_old.py oldinputfile.inp newinputfile.inp
#
# A few notes:
#
#(1) The script will try to load the following files:
#        MCF file: it will use the path specified in oldinputfile.inp
#        CHK File: it will look for it in the current directory
#        XYZ File: it will look for it in the current directory
#        H File: it will look for it in the current directory
#
#(2) The script will create the following files whose name will be derived from old run name.
#        *XYZ_box?.xyz_new: Last configuration created from a previous simulation. 
#			    The headers of each configuration will be the number of molecules of each species.
#        newinputfile.inp: New input file containing new maximum displacements from the CHK file
#                          Start_Type set to read_old using the new *XYZ_box?.xyz_new files
#                          New box sizes taken from CHK File
#
#********************************************************************************




#!/usr/bin/env python
import sys, os, linecache, argparse
from numpy import dot

def truncate(number, digits):
	a=number.split('.')
	return ".".join([a[0],a[1][0:digits]])

####PARSING ARGUMENTS####

parser = argparse.ArgumentParser(description = "Helps to restart a simulation using the read_old keyword in Cassandra. The new xyz files will be created in the current directory. Maximum displacements taken from CHK file.")
parser.add_argument('inputfile', action="store", help="Original input file")
parser.add_argument('newfile', action="store", help="New input file created")
results = parser.parse_args()

input_file = results.inputfile
newinputfile = results.newfile

#########################

truncation_limit=3

in_file = open(input_file,'r')


#####LOCATE KEYWORDS IN INPUT FILE#####

keyword_line = []
for line_number,line in enumerate(in_file):
	if not line.strip():
		continue
	line_parse = line.split()
	if line_parse[0]=='#':		
		if line_parse[1] == "Run_Name":
			run_name_line = line_number+1
		elif line_parse[1] == "Sim_Type":
			sim_type_line = line_number+1
		elif line_parse[1] == "Nbr_Species":
			nbr_species_line = line_number+1
		elif line_parse[1] == "Molecule_Files":
			molec_files_line = line_number+1
		elif line_parse[1] == "Box_Info":
			box_info_line = line_number+1


in_file.close()

run_name = linecache.getline(input_file,run_name_line+1)[0:-1]
sim_type = linecache.getline(input_file,sim_type_line+1)[0:-1]

####Obtain Box Info####

nboxes = int(linecache.getline(input_file, box_info_line+1))

box_type=[]
line_box_info_begin_end=[box_info_line]
i=2
box_count=0
while 1: 
	if 'CUBIC' in linecache.getline(input_file,box_info_line+i) or 'ORTHOGONAL' in linecache.getline(input_file,box_info_line+i) or 'CELL_MATRIX' in linecache.getline(input_file,box_info_line+i):
		box_type.append(linecache.getline(input_file,box_info_line+i).split()[0])
		box_count=box_count+1
	if box_count == nboxes:
		line_box_info_begin_end.append(box_info_line+i+1)
		break
	i = i + 1


nbr_species = int(linecache.getline(input_file,nbr_species_line+1))

xyz_file=[]
h_file=[]
mcf_files = []

for i in xrange(nboxes):
        xyz_file.append(run_name + "_box"+str(i+1)+".xyz")
        h_file.append(run_name + "_H.box"+str(i+1))

for i in xrange(0,nbr_species):
	mcf_files.append(linecache.getline(input_file,molec_files_line+i+1).split()[0])





####Obtain number of atoms per species from MCF Files####

atom_per_species = []
for i in xrange(0,len(mcf_files)):
	current_mcf = open(mcf_files[i],'r')
	for line_number_mcf, line_mcf in enumerate(current_mcf):
		
		if not line_mcf.strip():
			continue
		line_mcf_parse = line_mcf.split()
		if len(line_mcf_parse)>1:
			if line_mcf_parse[1] == "Atom_Info":
				atom_per_species.append(int(linecache.getline(mcf_files[i], line_number_mcf+2)))
				current_mcf.close()
				break

####Extract XYZ files####

#line_number will contain the line where the last config starts

line_number=[]
number_molecules=[]
total_lines=[]
out_xyz_file=[]

for j in xrange(nboxes):

	os.system('tail -'+str(1+nbr_species)+' '+h_file[j] +' > temp.dat'+str(j))

	file = open('temp.dat'+str(j),'r')
	total_lines_temp = len(file.readlines())
	file.close

	for i in xrange(total_lines_temp):
		if i==0: continue
		number_molecules.append(int(linecache.getline('temp.dat'+str(j),i+1).split()[1]))

	number_of_lines_last_snapshot = dot(number_molecules, atom_per_species)

	os.system('tail -'+str(number_of_lines_last_snapshot+2)+' '+xyz_file[j] +' > tempxyz_new'+str(j))

	file_old = open('tempxyz_new'+str(j),'r')

	new_xyz = open(xyz_file[j] + "_new",'w')
	for line_number, line in enumerate(file_old):
		if line_number == 0:
			for x in number_molecules:
				new_xyz.write(str(x)+" ")
			new_xyz.write('\n')	
			continue

		if not line.strip():
			continue
		new_xyz.write(line)

	file_old.close()
	new_xyz.close()
	number_molecules=[]
	os.system('rm temp.dat'+str(j)+'; rm tempxyz_new'+str(j))
	out_xyz_file.append(xyz_file[j] + '_new')





####READ CHECKPOINT FILE#####

checkfile = run_name + ".chk"

file = open(checkfile,'r')
for number, line in enumerate(file):
	if line ==' ******** Box info ***********\n':
		break
		
#number is one line before Box info tag

boxlength=[]
maxvol=[]
number=number+5
for i in xrange(nboxes):
	if "CELL_MATRIX" in box_type[i]:
		boxlength.append([linecache.getline(checkfile,number+0), linecache.getline(checkfile,number+1),linecache.getline(checkfile,number+2)])
		number = number + 6
		maxvol.append(linecache.getline(checkfile,number).split()[0])
		number = number + 4

	elif "ORTHOGONAL" in box_type[i]:
		boxlength.append([linecache.getline(checkfile,number+0).split()[0],linecache.getline(checkfile,number+1).split()[1],linecache.getline(checkfile,number+2).split()[2]])
		number = number + 6
		maxvol.append(linecache.getline(checkfile,number).split()[0])
		number = number + 4

	elif "CUBIC" in box_type[i]:
		boxlength.append(linecache.getline(checkfile,number).split()[0])
		number = number + 6
		maxvol.append(linecache.getline(checkfile,number).split()[0])
        	number = number + 4

maxtrans_rot_dih = []
line_number = 4
for i in xrange(nboxes):
	for j in xrange(nbr_species):
		maxtrans_rot_dih.append(linecache.getline(checkfile,line_number).split())
		line_number = line_number + 3
	line_number = line_number+1

file.close()


####CREATE NEW INPUT FILE####

file_new = open(newinputfile,'w')
file_old = open(input_file,'r')
total_lines_old = len(file_old.readlines())

stride = 0
for i in xrange(total_lines_old):
	thisline = linecache.getline(input_file,i+stride)


	if "Box_Info" in thisline:

		file_new.write('# Box_Info\n')
		file_new.write(str(nboxes) + '\n')
		for j in xrange(nboxes):
			
			if "CELL_MATRIX" in box_type[j]:
				file_new.write(box_type[j]+"\n")
				for k in boxlength[j]:
					file_new.write(str(k))
				file_new.write("\n")
				#stride = 2*nboxes+3+stride

			elif "ORTHOGONAL" in box_type[j]:
				file_new.write(box_type[j]+"\n")
				for k in boxlength[j]:
					file_new.write(str(k)+" ")
				file_new.write("\n\n")
				#stride = 2*nboxes+3+stride	


			elif "CUBIC" in box_type[j]:
				file_new.write(box_type[j]+"\n")
				file_new.write(boxlength[j]+'  '+boxlength[j]+'  '+boxlength[j]+'\n\n')
		stride = line_box_info_begin_end[1]-line_box_info_begin_end[0]+1 +stride


	if "Prob_Translation" in thisline:

                file_new.write('# Prob_Translation\n')
                file_new.write(linecache.getline(input_file,i+1+stride))
		for k in xrange(nboxes*nbr_species):
        	        file_new.write(truncate(maxtrans_rot_dih[k][0],truncation_limit)+' ')
			if (k+1)%nbr_species == 0:
				file_new.write('\n')
		stride = 2*nboxes+ stride+1
		file_new.write('\n')


	if "Prob_Rotation" in thisline:

                file_new.write('# Prob_Rotation\n')

                file_new.write(linecache.getline(input_file,i+1+stride))
                for k in xrange(nboxes*nbr_species):
                        file_new.write(truncate(str(float(maxtrans_rot_dih[k][1])*180/3.1416),truncation_limit)+' ')
                        if (k+1)%nbr_species == 0:
                                file_new.write('\n')

		stride = 2*nboxes+ stride+1
                file_new.write('\n')


	if "Prob_Volume" in thisline:
                file_new.write('# Prob_Volume\n')
                file_new.write(linecache.getline(input_file,i+1+stride))
                for k in xrange(nboxes):
                        file_new.write(truncate(maxvol[k],truncation_limit)+' \n')

		stride = nboxes+2 + stride
                file_new.write('\n')
                

	if "Start_Type" in thisline:
		file_new.write('# Start_Type\n')
                file_new.write('read_old\n')
                for k in xrange(nboxes):
                        file_new.write(out_xyz_file[k]+' \n')

		stride = 2+nboxes+stride
                file_new.write('\n')

	else:
		insertedline = linecache.getline(input_file,i+stride)

		file_new.write(insertedline)

file_new.close()

