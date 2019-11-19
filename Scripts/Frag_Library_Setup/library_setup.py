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
#     Cassandra input file
#     PDB files for each species
#     MCF files for each species
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

import sys
import os
import subprocess
import argparse
from linecache import getline
import re
import math

def main():

    args = get_args()

    cassandra_path = os.path.abspath(args.cassandra_path)
    input_file = os.path.abspath(args.input_file)

    pdb_files = get_pdb_files(args.config_files):

    output = ( color.BOLD +
               "\n******************* Cassandra Setup *******************\n"
               "Cassandra location: " + color.END + cassandra_path +
               "\n" + color.BOLD +
               "Scanning input file\n" + color.END
             )
    print(output)

    keyword_line = parse_input_keywords(input_file)

    # Extract some info from input (.inp) file
    nbr_boxes = int(getline(input_file,keyword_line['Box_Info']+1)
    sim_type = getline(input_file,keyword_line['Sim_Type']+1)
    nbr_species = int(getline(input_file,keyword_line['Nbr_Species']+1))
    output = ( color.BOLD + "  Number of species found: " +
               color.END + str(nbr_species)
             )
    print(output)

    # Extract more info from input file (>1 line so move to fxns)
    mcf_files = parse_molecule_filenames(input_file,keyword_line,nbr_species)
    temperature = parse_temperature(input_file,keyword_line,nbr_boxes)
    vdw_style = parse_vdw_style(input_file,keyword_line,nbr_boxes)
    charge_style = parse_charge_style(input_file,keyword_line,nbr_boxes)
    mixing_rule = parse_mixing_rule(input_file,keyword_line)
    vdw_scaling, charge_scaling = parse_intra_scaling(input_file,
                                                      keyword_line,
                                                      nbr_boxes,
                                                      charge_style)

    # Create a dictionary with many of the inp data make some later
    # function calls a tad cleaner. None of this data will be changed
    # following the initial read above
    inp_data = { 'nbr_boxes' : nbr_boxes,
                 'sim_type' : sim_type,
                 'nbr_species' : nbr_species,
                 'temperature' : temperature,
                 'vdw_style' : vdw_style,
                 'charge_style' : charge_style,
                 'mixing_rule' : mixing_rule,
                 'vdw_scaling' : vdw_scaling,
                 'charge_scaling' : charge_scaling
               }
    
    # Find species that are fixed...assuming fixed if pdb does not
    # have CONECT records --> need to find a better way of doing this
    rigid_species = check_pdb_conect(pdb_files)

    # Parse MCF files for keywords
    mcf_keyword_line = parse_mcf_keywords(mcf_files)

    # Parse atom info section of each mcf file
    atom_types, ring_atoms = parse_mcf_atom_info(mcf_files,mcf_keyword_line,
                                                 args.noFlags)
    nbr_atoms = [len(sp_atom_types) for sp_atom_types in atom_types]

    # Parse fragment info section of each mcf file
    fragment_list = parse_mcf_frag_info(mcf_files,mcf_keyword_line)
    nbr_fragments = [len(sp_frag_list) for sp_frag_list in fragment_list]
    # Append 's_idx' to the atom types of the mcf files
    append_speciesidx(mcf_files)

    frag_ring_status = id_ring_fragments(fragment_list,ring_atoms)
    # TODO: Write this.
    # I don't think I need ring_is_rigid --> should be able to combine
    # frag_rigid_status and frag_ring_status together
    frag_rigid_status = id_rigid_fragments(fragment_list,ring_atoms)

    #####################################################################
    #####################################################################

    # Let the fun begin. Create folder structure
    create_directory_structure(nbr_species,nbr_fragments)
    # Calls to Cassandra to generate fragment MCF files
    create_fragment_mcf_files(inp_data,nbr_atoms,nbr_fragments,
                              rigid_species,frag_ring_status,
                              mcf_files,pdb_files,cassandra_path)


def frag_rigid_status():
    """Identify the rigid fragments"""
    frag_rigid_status = []
    for ispecies in range(nbr_species):
        if nbr_atoms[ispecies] < 3:
            frag_rigid_status.append([True])
        elif ispecies in rigid_species:
            frag_rigid_status.append([True])
        else:
            for jfrag in range(nbr_fragments[ispecies]):
                nbr_angles_fixed = 0
                filen = ("species" + str(ispecies+1) + "/fragments/frag_" +
                         str(jspecies+1) + "_1.mcf")
                with open(filen) as frag_mcf:
                # RSD TODO:: WORKING HERE NOW

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
        new_file.write(getline(input_file,line_number))
    elif line_number > frag_files_line:
        if getline(input_file,line_number)[0] == '#' or \
           'END' in getline(input_file,line_number):
            omit = False
            new_file.write(getline(input_file,line_number))

in_file.close()
new_file.close()

print bold+"Removing temporary input file"+normal
os.system("rm "+ input_file+"; mv "+input_file+"temp "+input_file)
print bold+"Finished"+normal

def parse_input_keywords(input_file):
    """Return a dictionary of line numbers for each required keyword"""

    required_keywords = ['Nbr_Species',
                         'Molecule_Files',
                         'Start_Type',
                         'Fragment_Files',
                         'Temperature_Info',
                         'Sim_Type',
                         'VDW_Style',
                         'Box_Info' ]

    optional_keywords = [ 'Intra_Scaling',
                          'Charge_Style',
                          'Mixing_Rule' ]

    # Extract keyword line numbers
    keyword_line = {}
    with open(input_file) as f:
        for line_number, line in enumerate(f):
            line_parse = line.strip().split()
            try:
                if line_parse[0] == '#':
                    keyword = line_parse[1]
                    if keyword in required_keywords:
                        keyword_line[keyword] = line_number+1
                    elif keyword in optional_keywords:
                        keyword_line[keyword] = line_number+1

    # Check that all required keywords are present
    if not all(keyword in keyword_line for keyword in required_keywords):
        raise ValueError("Input file does not have all required "
                "keywords. Required keywords are: {}".format(
                    required_keywords))

    return keyword_line

def parse_molecule_filenames(input_file,keyword_line,nbr_species):
    """Parse the input file and extract the MCF file names"""

    mcf_files = []
    for i in range(0,nbr_species):
        mcf_file = getline(input_file,
                         keyword_line['Molecule_Files'+i+1).strip()
        if mcf_file == '':
            raise ValueError("Must specify one MCF file per species.")
        try:
            if mcf_file[0] == "!"
        except IndexError:
            raise ValueError("Must specify one MCF file per species.")
        mcf_files.append(mcf_file)
        output = ( color.BOLD +
                   "   The MCF file number " + str(i+1) + "is: " +
                   color.END + mcf_files[i]
                 )
        print(output)

    return mcf_files

def parse_temperature(input_file,keyword_line,nbr_boxes):
    """Parse the input file and extract the temperature"""
    temperature_list = []
    for i in range(1,nbr_boxes+1):
        try:
            temperature_list.append(float(getline(input_file,
                                keyword_line['Temperature_Info']+i).strip()))
        except ValueError:
            raise ValueError("One temperature must be specified for "
                    "each simulation box. See manual for details.")
        if temperature_list[i-1] != temperature_list[0]:
            raise ValueError("Temperature of box {} does not match "
                    "the temperature of box {}.".format(1,i))

    return temp_list[0]

def parse_vdw_style(input_file,keyword_line,nbr_boxes):
    """Parse the input file and return the VDW style"""

    vdw_style = []
    for i in range(1,nbr_boxes+1):
        vdw_style.append(getline(input_file,
                         keyword_line['VDW_Style']+i).split()[0])
        if vdw_style[i-1] != vdw_style[0]:
            raise ValueError("VDW_Style for boxes {} and {} does not "
                            "match. Please fix your input file.".format(
                                1,i))

    vdw_style = vdw_style[0]
    if vdw_style.lower() != 'lj' and vdw_style.lower() != 'mie':
        raise ValueError("Only 'LJ' and 'Mie' VDW_Styles are supported")

    return vdw_style

def parse_charge_style(input_file,keyword_line,nbr_boxes):
    """Parse the input file and return the charge style"""

    if keyword_line['Charge_Style'] is None:
        return None

    charge_style = []
    for i in range(1,nbr_boxes+1):
        charge_style.append(getline(intput_file,
                            keyword_line['Charge_Style']+i).split()[0])
        if charge_style[i-1] != charge_style[0]:
            raise ValueError("Charge_Style for boxes {} and {} does not "
                        "match. Please fix your input file.".format(1,i))
    charge_style = charge_style[0]
    if charge_style.lower() != 'coul' and charge_style.lower() != 'none':
        raise ValueError("Only 'Coul' and 'None' Charge_Styles are supported")

    return charge_style

def parse_mixing_rule(input_file,keyword_line):
    """Parse the input file and return the mixing rules section"""

    mixing_line = keyword_line['Mixing_Rule']
    if mixing_line is None:
        return None
    else:
        mixing_rule = getline(input_file,mixing_line+1)
        if mixing_rule.strip() == "custom":
            i = 1
            while True:
                if getline(input_file,mixing_line+1+i).strip() and \
                   getline(input_file,mixing_line+1+i)[1] != '!':
                    mixing_rule += getline(input_file,mixing_line+1+i)
                    i += 1
                else:
                    break
    return mixing_rule

def parse_intra_scaling(input_file,keyword_line,nbr_boxes,charge_style):
    """Parse the input file and return the VDW and charge intra scaling"""

    vdw_scaling = []
    charge_scaling = []
    for i in range(1,nbr_species+1):
        if keyword_line['Intra_Scaling']:
            if charge_style.lower() == 'coul':
                vdw_scaling.append(getline(input_file,
                                    keyword_line['Intra_Scaling']+2*i-1))
                charge_scaling.append(getline(input_file,
                                    keyword_line['Intra_Scaling']+2*i))
            else:
                vdw_scaling.append(getline(input_file,
                                   keyword_line['Intra_Scaling']+i))

    return vdw_scaling, charge_scaling

def parse_mcf_keywords(mcf_files):

    keywords = ['Atom_Info','Fragment_Info']
    mcf_keyword_line = { keyword : [] for keyword in keywords }

    for i, mcf_file in enumerate(mcf_files):
        with open(mcf_file) as mcf:
            for line_number, line in enumerate(mcf):
                line_parse = line.strip().split()
                try:
                    if line_parse[0] == '#':
                        keyword = line_parse[1]
                        if keyword in keywords:
                            mcf_keyword_line[keyword].append(line_number+1)

    for key,val in mcf_keyword_line.items():
        if len(val) != len(mcf_files):
            raise ValueError("Not all MCF files have the required "
                "keywords: {}".format(keywords))

    return mcf_keyword_line

def parse_mcf_atom_info(mcf_files,mcf_keyword_line,noflags):
    """Parse all MCF files and return list of atom types and ring atoms."""

    atom_line = mcf_keyword_line['Atom_Info']
    atom_types = []
    ring_atoms = []

    for ispecies,mcf_file in enumerate(mcf_files):
        nbr_atoms = int(getline(mcf_file,atom_line[ispecies]+1).split()[0])
        atom_types[ispecies].append([])
        ring_atoms = [ispecies].append([])
        for atidx in range(0,nbr_atoms):
            line = getline(mcf_file,atom_line[ispecies]+2+atidx).split()
            atom_type = line[1]
            #remove old '_s?' flags if present
            if atom_type[-3:-1] == '_s':
                atom_type = atom_type[:-3]
            #add '_s?" flags to atom_types, unless --noFlags given
            if not args.noFlags:
                atom_type = atom_type + '_s' + str(i+1)
            atom_types[ispecies].append(atom_type)
            # Check for ring flag
            if line[-1] == 'ring':
                species_ring_atoms.append(str(atidx+1))

        ring_atoms.append(species_ring_atoms)

    return atom_types, ring_atoms

def parse_mcf_fragment_info(mcf_files,mcf_keyword_line):
    """Parse all MCF files and return list of fragments for each species"""

    frag_line = mcf_keyword_line['Fragment_Info']
    for ispecies,mcf_file in enumerate(mcf_files):
        nbr_fragments = int(getline(mcf_file,frag_line[ispecies]+1).split()[0])
        output = ( color.BOLD +
                   "    Species "+ str(i+1) + " has " + str(nbr_fragments) +
                   " fragments" + color.END
                 )
        print(output)
        sp_frag_list = []
        for idx in range(nbr_fragments):
            sp_frag_list.append(getline(mcf_file,
                                    frag_line[ispecies]+2+idx).split()[2:])

        fragment_list.append(sp_frag_list)

    return fragment_list

def id_ring_fragments(fragment_list,ring_atoms):
    """Identify the ring status of each fragment as 'exo' 'ring' or False.

    'exo' indicates a fragment with ring and exoring atoms
    'ring' indicates a fragment with only ring atoms
    False indicates a fragment with no ring atoms

    The intersection of the atoms labeled 'ring' and the fragment lists
    for each species is used to identify ring fragments. The intersection
    must have at least three overlapping atoms for the fragment to be
    a ring fragment. A comparison between the size of the fragment and
    the number of intersecting ring atoms is used to identify the
    presence or absence of exoring atoms.
    """
    frag_ring_status = []
    for ispecies,sp_frag_list in enumerate(fragment_list):
        frag_ring_status.append([])
        for j,fragment in enumerate(sp_frag_list):
            intersection = set.intersection(set(fragment),
                                            set(ring_atoms[ispecies]))
            if len(intersection) > 2 and len(fragment) > len(intersection):
                frag_ring_status[ispecies].append('exo')
            elif len(intersection) > 2:
                frag_ring_status[ispecies].append('ring')
            else:
                frag_ring_status[ispecies].append(False)

    if any([any(has_ring) for has_ring in frag_ring_status]):
        output = ( color.BOLD +
                  "Molecules with rings found. These are: " +
                   color.END
                 )
        print(output)
        for ispecies,sp_frag_ring_status in enumerate(frag_ring_status):
            output = ( color.BOLD +
                       "Species " + str(ispecies+1) "has " +
                       sum([status != False for status in sp_frag_ring_status]) +
                       " rings." + color.END
                     )
            print(output)

    return frag_ring_status

def append_speciesidx(mcf_files,atom_types):

    for ispecies,mcf_file in enumerate(mcf_files):
        nbr_atoms = len(atom_types[ispecies])
        with open(mcf_file) as mcf:
            total_lines = len(mcf.readlines())
        with open('temp.mcf') as new_mcf:
            atidx = 0
            for line_number in range(1,total_lines+1):
                if line_number > atom_line[ispecies]+1 and \
                   line_number <= atom_line[ispecies]+1 + nbr_atoms:
                    line = getline(mcf_file,line_number).split()
                    line[1] = atom_types[ispecies][atidx]
                    atidx+=1
                    new_mcf.write('    '.join(line)+'\n')
                else:
                    new_mcf.write(getline(mcf_file,line_number))
        os.system("mv "+'temp.mcf '+mcf_file)

def create_directory_structure(nbr_species,nbr_fragments):
    """Setup the directory structure for the fragment library

    Create N number of folders, where N is the number of species
    Inside each folder, we'll include the following folders:
    /species?/fragments
    /species?/frag?
    """

    for ispecies in range(nbr_species):
        for frag_idx in range(nbr_fragments[ispecies]):
            command = ("mkdir -p species" + str(ispecies+1) +
                       "/frag"+ str(frag_idx+1))
            os.system(command)
        if ispecies not in rigid_species:
            command = ("mkdir -p species" + str(i+1) + "/fragments/")
            os.system(command)

def check_type_infilename(infilename):
    """Return the filetype (CML or PDB). Assumes PDB if not CML."""

    with open(infilename) as f:
        for line in f:
            if '<molecule>' in line:
                return 'cml'
    return 'pdb'

def get_pdb_files(config_files):
    """Return a list of pdb files from the config files.

    Converts any cml files to pdb files as required.
    """
    pdb_files = []
    for index,each_file in enumerate(config_files):
        if check_type_infilename(each_file)=='cml':
            pdb_file = cml_to_pdb(each_file)
        else:
            pdb_file = each_file
        pdb_files.append(pdb_file)

    return pdb_files

def check_pdb_conect(pdb_files):
    """Return species number that do not have conect records.

    These files are assumed to be rigid fragments (e.g., zeolite)
    """

    rigid_species = []
    for this_species,each_pdb in enumerate(pdb_files):
        with open(each_pdb) as pdb:
            conect_found = False
            for line in this_pdb:
                if 'CONECT' in line:
                    conect_found = True
                    break
            if conect_found == False:
                rigid_species.append(this_species)

    return rigid_species


def is_number(a):
    try:
        float(a)
    except ValueError:
        return False
    return True

def cml_to_pdb(infilename):

    # ID starting and ending lines of atom/bond sections
    with open(infilename) as f:
        for line_nbr, line in enumerate(f):
            if '<bondArray>' in line:
                cml_start_bonds = line_nbr + 1
            if '</bondArray>' in line:
                cml_end_bonds = line_nbr + 1
            if '<atomArray>' in line:
                cml_start_atom = line_nbr + 1
            if '</atomArray>' in line:
                cml_end_atom = line_nbr + 1

    cml_atom_info=[]
    for i in xrange(cml_start_atom+1, cml_end_atom):
        cml_atom_info.append(re.findall('"([^"]*)"',getline(infilename, i)))

    cml_bond_info=[]
    for i in xrange(cml_start_bonds+1, cml_end_bonds):
        cml_bond_info.append(re.findall('"([^"]*)"',getline(infilename, i))[0].split())

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

def create_fragment_mcf_files(inp_data,nbr_atoms,nbr_fragments,
                              rigid_species,frag_ring_status,
                              mcf_files,pdb_files,cassandra_path):
    """Create MCF files for each fragment
    """
    nbr_species = inp_data['nbr_species']

    for ispecies in range(nbr_species):
        # Extract the mcf file for later
        mcf_file = mcf_file[ispecies]

        # First check if rigid species. If it is we do nothing here.
        if ispecies in rigid_species:
            output = ("\n\n" + color.BOLD +
                      "MCF generation file not created for species " +
                      str(ispecies+1) + color.END + "\n")
            if nbr_fragments[ispecies] == 0:
                output += (color.BOLD + "No fragment configuration needed." +
                          color.END )
            else:
                output += (color.BOLD +
                          "Fragment configuration will be taken "
                          "from the PDB file." + color.END)
            print(output)
            continue

        # All species that are not rigid and include more than two atoms
        if nbr_atoms[ispecies] >= 3:
            output = "\n\n" + color.BOLD +
                     "Creating input file for MCF generation file for " +
                     "species " + str(ispecies+1) + " " + color.END)
            print(output)

            # Write a cassandra input file
            write_input_mcfgen(ispecies,inp_data,mcf_file)

            # If ring fragment we also need to copy the pdb file here
            if sum([ status != False
                     for status in frag_ring_status[ispecies]]) > 0:
                pdb_search_name = os.path.splitext(os.path.basename(
                                                   mcf_file)) + '.pdb'
                found_pdb = False
                for pdb_path in pdb_files:
                    pdb_name = os.path.basename(pdb_path)
                    if pdb_name == pdb_search_name:
                        command = ("cp " + pdb_path + " species" + 
                                    str(i+1) + "/fragments/molecule.pdb")
                        os.system(command)
                        found_pdb = True
                        break
                if found_pdb == False:
                    raise ValueError("Unable to find PDB file "
                            "corresponding to {}".format(mcf_file))

            output = (color.BOLD +
                      "Running Cassandra to generate MCF files" +
                      color.END)
            print(output)
            sys.stdout.flush()


            # Now we can actually call Cassandra (only to generate an MCF)
            os.chdir('./species'+str(i+1)+'/fragments/')
            subprocess.call([cassandra_path,'species'+str(i+1)+'_mcf_gen.inp'])
            os.chdir('../../')



def write_input_mcfgen(ispecies,inp_data,mcf_file):

    filename = ('./species' + str(ispecies+1) + '/fragments/' +
                'species' + str(ispecies+1) + '_mcf_gen.inp' )

    inp_content = file_templates.inp_base_content.format(
                              run_name = "species" + str(i+1) + "_mcfgen",
                              sim_type = "MCF_Gen",
                              mcf_file = "../../" + mcf_file,
                              vdw_style = vdw_style,
                              rcutoff_low = rcutoff_low)

    if keyword_lines['Mixing_Rule'] is not None:
        inp_content += "# Mixing_Rule\n"
        inp_content += inp_data['mixing_rule']
        inp_content += "\n\n"
    if keyword_lines['Charge_Style'] is not None:
        inp_content += "# Charge_Style\n"
        inp_content += inp_data['charge_style']
        if charge_style.lower() == 'coul':
            inp_content += " minimum image"
        inp_content += "\n\n"
    if keyword_lines['Intra_Scaling'] is not None:
        inp_content += "# Intra_Scaling\n"
        inp_content += inp_data['vdw_scaling'][ispecies]
        inp_content += "\n"
        if charge_style.lower() == 'coul':
            inp_content += inp_data['charge_scaling'][ispecies]
            inp_content += "\n"
        inp_content += "\n"

    inp_content += "END"

    with open(filename,'w') as finp:
        finp.write(inp_content)



def read_pdb(pdb_file):
    """ Read a PDB file and return the atomic coordinates.

    Assumes all lines with coordinate information start with
    'HETATM' or 'ATOM'. Assumes PDB file has coordinates written
    in correct columns x = [30:38], y = [48:46], z = [46:54].

    Parameters:
    -----------
    pdb_file : string
        path to pdb file

    Returns:
    --------
    atom_coords : dict
        Atom coordinates. Key is the atom number (starting from 1).
        Value is a tuple, (x,y,z) coordinates.
    """
    atom_coords = {}
    atom_nbr = 1
    with open(pdb_file) as pdb:
        for line in pdb:
            # read atom info
            this_line = line.split()
            if line[0:6]=='HETATM' or line[0:4]=='ATOM':
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                atom_coords[atom_nbr] = (x,y,z)
                atom_nbr += 1

    return atom_coords

def parse_args():
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

    return args

class color:
    PURPLE = '\033[95m'
    CYAN = '\033[96m'
    DARKCYAN = '\033[36m'
    BLUE = '\033[94m'
    GREEN = '\033[92m'
    YELLOW = '\033[93m'
    RED = '\033[91m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'
    END = '\033[0m'

class file_templates:
    inp_base_content = """
# Run_Name
{run_name}

# Sim_Type
{sim_type}

# Nbr_Species
1

# Molecule_Files
{mcf_file} 1

# VDW_Style
{vdw_style} minimum_image

# Rcutoff_Low
{rcutoff_low}

# Box_Info
1
CUBIC
50.0 50.0 50.0

"""

    inp_lib_content = """
"""


