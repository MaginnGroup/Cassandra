#****************************************************************************
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
#****************************************************************************

#****************************************************************************
# SCRIPT: mcfgen.py
# VERSION: 1.1
# NEW FEATURES: read GROMACS forcefield files (.itp or .top)
#
# VERSION: 1.0
# ORIGINAL FEATURES: generate a molecular connectivity file (.mcf) from a 
#   configuration file (.pdb or .cml) and a custom forcefield file (.ff)
#****************************************************************************

import sys
import os
import argparse
import linecache
import networkx as nx
import re

def main():

    args = parse_args()

    # Check the input files
    configFile, configFileType, basename = check_config_file(args.configFile)
    mcfFile = check_mcf_file(args.mcfFile, basename, args.ffTemplate)
    ffFile, ffFileType = check_ff_file(args.ffFile, basename, args.ffTemplate)

    # Always start by reading PDB file
    print("Reading Modified PDB File...", end='')
    atom_list, atom_parms, atomtype_parms, molecule_graph = read_pdb(configFile)
    print("success\n")

    # RSD TODO: Plug and play parts of MCF writer from mBuild.
    # Next analyze bonds, angles, dihedrals, fragments, and rings
    bond_list = id_bonds(molecule_graph)
    print("Created bond_list...",end='')
    angle_list = id_angles(molecule_graph)
    print("angle_list...",end='')
    dihedral_list = id_dihedrals(molecule_graph)
    print("dihedral_list...")
    in_ring, frag_list, frag_conn = id_rings_fragments(molecule_graph)
    print("Identified fragments and rings...\n\nSummary:")
    print("----------------------------------------")

    print("{:4d} {:40s}".format(len(bond_list),'bonds'))
    print("{:4d} {:40s}".format(len(angle_list),'angles'))
    print("{:4d} {:40s}".format(len(dihedral_list),'dihedrals'))
    print("{:4d} {:40s}".format(len(frag_list),'fragments'))
    print("----------------------------------------")

    # Two uses of mcfgen: generate blank FF template or generate MCF file
    if args.ffTemplate:
        # Generate FF template
        write_ff_template(ffFile, atom_list, bond_list, angle_list,
                dihedral_list, atom_parms, atomtype_parms)
        print('Finished')
    else:
        # Generate MCF file
        # Read parms
        if ffFileType == 'native':
            atom_parms, atomtype_parms, bond_parms, angle_parms, dihedral_parms, improper_parms, scaling_1_4 = \
                read_native_ff(ffFile, atom_parms, atomtype_parms)
            if dihedral_style == 'none':
                dihedral_list = []
        elif ffFileType == 'gromacs':
            atom_parms, atomtype_parms, bond_parms, angle_parms, dihedral_parms, scaling_1_4 = \
                read_gromacs_ff_wrapper(ffFile, atom_parms, atomtype_parms)
            improper_parms = {}

        # Check parms
        atom_parms = check_atom_parms(atom_list, atom_parms, atomtype_parms)
        bond_parms = check_bond_parms(bond_list, bond_parms, atom_parms)
        angle_parms = check_angle_parms(angle_list, angle_parms, atom_parms)
        dihedral_parms = check_dihedral_parms(dihedral_list, dihedral_parms,
                atom_parms)

        # Got all the parms we need? write MCF.
        with open(mcfFile,'w') as mcf:
            header = ( '!****************************************'
                       '*****************************************\n'
                       '!Molecular connectivity file for {} \n'
                       '!****************************************'
                       '*****************************************'
                       '\n\n\n'.format(configFile)
                     )
            mcf.write(header)

            write_mcf_atoms(mcf, atom_list, atom_parms, in_ring)
            write_mcf_bonds(mcf, atom_list, bond_list, bond_parms)
            write_mcf_angles(mcf, atom_list, angle_list, angle_parms)
            write_mcf_dihedrals(mcf, atom_list, dihedral_list, dihedral_parms)
            write_mcf_impropers(mcf, atom_list, improper_parms)
            write_mcf_fragments(mcf, atom_list, frag_list, frag_conn)
            write_mcf_intra_scaling(mcf, scaling_1_4)

    # Clean up PDB file if input was CML
    if configFileType == 'cml':
        os.system("rm " + configFile)

#****************************************************************************
#                          END OF MAIN FUNCTION
#****************************************************************************

def id_rings_fragments(molecule_graph):
    """Identifies the rings and fragments of the molecule

    See Shah and Maginn, JCP, 135, 134121, 2011, doi:10.1063/1.3644939
    for a description of fragments

    Parameters
    ----------
    molecule_graph : nx.Graph
        Networkx graph with connectivity of atoms

    Returns
    ---------
    in_ring : list
        True for each atom in a ring
    frag_list : list
        Atom ids belonging to each fragment
    frag_conn : list
        Fragment ids of connected fragments

    """

    # Create a few structures
    frag_list = []
    frag_conn = []
    in_ring = { pdb_idx : False for pdb_idx in molecule_graph.nodes() }

    # First create a neighbor list for each atom
    neigh_dict = { node : list(neigh for neigh in
                               molecule_graph.neighbors(node))
                          for node in molecule_graph.nodes() }

    # Check not null graph (single particle PDB)
    #if molecule_graph.number_of_nodes() == 0:
    #    return in_ring, [], []

    # Check if entire molecule is connected. Warn if not.
    if nx.is_connected(molecule_graph) == False:
        print("WARNING: Not all components of the molecule are connected. "
              "MCF files are generally for a single molecule and thus "
              "everything should be connected through bonds.")

    all_rings = nx.cycle_basis(molecule_graph)
    for ring in all_rings:
        for idx in ring:
            in_ring[idx] = True

    # ID fused rings
    fused_rings = []
    rings_to_remove = []
    for i in range(len(all_rings)):
        ring1 = all_rings[i]
        for j in range(i+1, len(all_rings)):
            ring2 = all_rings[j]
            shared_atoms = list(set(ring1) & set(ring2))
            if len(shared_atoms) == 2:
                fused_rings.append(list(set(ring1+ring2)))
                rings_to_remove.append(ring1)
                rings_to_remove.append(ring2)
    for ring in rings_to_remove:
        all_rings.remove(ring)
    all_rings = all_rings + fused_rings
    # ID fragments which contain a ring
    for ring in all_rings:
        adjacentatoms = []
        for idx in ring:
            if len(neigh_dict[idx]) > 2:
                adjacentatoms.append(list(set(neigh_dict[idx])-set(ring)))
        tmp=filter(None, adjacentatoms)
        adjacentatoms = [x for sublist in tmp for x in sublist]
        frag_list.append(ring+adjacentatoms)
    # Now ID the other fragments
    for idx in neigh_dict:
        if len(neigh_dict[idx]) > 1:
            if in_ring[idx] == True:
                continue
            else:
                frag_list.append([idx]+neigh_dict[idx])
    # Now find connectivity (shared bonds)
    for i in range(len(frag_list)):
        frag1 = frag_list[i]
        for j in range(i+1, len(frag_list)):
            frag2 = frag_list[j]
            shared_atoms = list(set(frag1) & set(frag2))
            if len(shared_atoms) == 2:
                frag_conn.append([i, j])
            elif len(shared_atoms) > 2:
                warnings.warn('Fragments share more than two atoms...'
                      'something may be going awry unless there are'
                      'fused rings in your system. See below for details.')
                print('Fragment 1 atoms:')
                print(frag1)
                print('Fragment 2 atoms:')
                print(frag2)

    return in_ring, frag_list, frag_conn

def id_bonds(molecule_graph):
    """Identify the bonds in a molecule

    Parameters
    ----------
    molecule_graph : nx.Graph()
        Networkx graph describing molecule connectivity

    Returns
    -------
    bond_list : list
        List of tuples, one for each unique bond

    """
    bond_list = []
    for atom1, atom2 in molecule_graph.edges():
        bond_list.append((atom1,atom2))

    return bond_list

def id_angles(molecule_graph):
    """Identify the angles in a molecule

    Parameters
    ----------
    molecule_graph : nx.Graph()
        Networkx graph describing molecule connectivity

    Returns
    -------
    angle_list : list
        List of tuples, one for each unique angle

    """
    angle_list = []
    for atom2 in molecule_graph.nodes():
        neighs = [atom for atom in molecule_graph.neighbors(atom2)]
        if len(neighs) > 1:
            for i in range(len(neighs)):
                for j in range(i+1,len(neighs)):
                    angle_list.append((neighs[i],atom2,neighs[j]))

    return angle_list

def id_dihedrals(molecule_graph):
    """Identify the dihedrals in a molecule

    Parameters
    ----------
    molecule_graph : nx.Graph()
        Networkx graph describing molecule connectivity

    Returns
    -------
    dihedral_list : list
        List of tuples, one for each unique dihedral

    """
    dihedral_list = []

    for atom1 in molecule_graph.nodes():
        for atom2 in molecule_graph.neighbors(atom1):
            for atom3 in molecule_graph.neighbors(atom2):
                if atom1 == atom3:
                    continue
                for atom4 in molecule_graph.neighbors(atom3):
                    if atom1 == atom4 or atom2 == atom4:
                        continue
                    ijkl = (atom1,atom2,atom3,atom4)
                    lkji = ijkl[::-1]
                    if ijkl not in dihedral_list and lkji not in dihedral_list:
                        dihedral_list.append(ijkl)

    return dihedral_list

def write_ff_template(ffFile, atom_list, bond_list, angle_list, dihedral_list,
        atom_parms, atomtype_parms):

    print("\n\n*********Force field template file generation*********\n")

    global vdw_style
    global dihedral_style

    # Extract the number of atomtypes
    n_atomtypes = len(atomtype_parms)

    # Get vdw_style and dihedral_style from user
    vdw_style = input("Enter the VDW type (LJ/Mie):")
    if len(dihedral_list) > 0:
        dihedral_style = input("Enter the dihedral functional form "
		               "(CHARMM/OPLS/harmonic/none): ")
    else:
        dihedral_style = "NONE"
    if not dihedral_style or dihedral_style == "NONE" or dihedral_style == "none":
        dihedral_style = "none"

    # Identify different types of angles and dihedrals based on atom types
    # Each row in dihedarlList_atomType has each dihedral expressed as atomType
    # Each row in angleList_atomType has each angle expressed as atomType

    with open(ffFile,'w') as ff:
        # Write list of atomtypes
        ff.write("atomtypes\n"+str(n_atomtypes)+ "\n\n")
        ff.write("begin atom-atomtype\n")
        for ff_idx, pdb_idx in enumerate(atom_list):
            atype = atom_parms[pdb_idx]['type']
            ff.write(str(ff_idx+1) + " " + atype + "\n")
        ff.write("end atom-atomtype\n\n")

        ff.write("vdwtype "+vdw_style+"\n")
        ff.write("dihedraltype "+dihedral_style+"\n")
        if dihedral_style != 'none':
            ff.write("scaling_1_4_vdw \n")
            ff.write("scaling_1_4_charge \n")
        ff.write("\n")

        # Write nonbonded entries
        for atype in atomtype_parms.keys():
            ff.write("nonbonded\n")
            ff.write(atype + "\n")
            ff.write("Sigma \n")
            ff.write("Epsilon \n")
            if vdw_style == 'Mie':
                ff.write("Repulsive_Exponent \n")
                ff.write("Dispersive_Exponent \n")
            ff.write("atom_type_charge \n\n")

        # Id unique bond types
        bond_types = set()
        for bond in bond_list:
            bond_type = (atom_parms[bond[0]]['type'],
                         atom_parms[bond[1]]['type'])
            if bond_type[::-1] not in bond_types:
                bond_types.add(bond_type)

        # Write bond entries
        for bond_type in bond_types:
            ff.write("bonds\n")
            ff.write(bond_type[0] + " " + bond_type[1] + "\n")
            ff.write("Length \n")
            ff.write("Constant \n\n")

        # Id unique angle types
        angle_types = set()
        for angle in angle_list:
            angle_type = (atom_parms[angle[0]]['type'],
                          atom_parms[angle[1]]['type'],
                          atom_parms[angle[2]]['type'])
            if angle_type[::-1] not in angle_types:
                angle_types.add(angle_type)

        # Write angle entries
        for angle_type in angle_types:
            ff.write("angles\n")
            ff.write(angle_type[0] + " " + angle_type[1] + " " +
                     angle_type[2] + "\n")
            ff.write("Angle \n")
            ff.write("Constant \n\n")

        # Id unique dihedral types
        dihedral_types = set()
        for dihedral in dihedral_list:
            dihedral_type = (atom_parms[dihedral[0]]['type'],
                             atom_parms[dihedral[1]]['type'],
                             atom_parms[dihedral[2]]['type'],
                             atom_parms[dihedral[3]]['type'])
            if dihedral_type[::-1] not in dihedral_types:
                dihedral_types.add(dihedral_type)

        # Write dihedral entries
        if dihedral_style != "none":
            for dihedral_type in dihedral_types:
                ff.write("dihedrals\n")
                ff.write(dihedral_type[0] + " " + dihedral_type[1] + " " +
                         dihedral_type[2] + " " + dihedral_type[3] + "\n")
                if dihedral_style == "CHARMM":
                    ff.write("K \n")
                    ff.write("n \n")
                    ff.write("Gamma \n\n")
                elif dihedral_style == "OPLS":
                    ff.write("a0 \n")
                    ff.write("a1 \n")
                    ff.write("a2 \n")
                    ff.write("a3 \n\n")
                elif dihedral_style == "harmonic":
                    ff.write("Angle \n")
                    ff.write("Constant \n\n")

        # Write charge entries
        for ff_idx,pdb_idx in enumerate(atom_list):
            ff.write("charge\n")
            ff.write(str(ff_idx+1) + " \n\n")

        ff.write("end\n")

def read_pdb(pdb_file):

    atom_list = []
    atom_parms = {}
    atomtype_parms = {}
    molecule_graph = nx.Graph()
    repeatedIndex = False

    with open(pdb_file) as pdb:

        for line in pdb:
            this_line = line.strip().split()
            # Read atom info
            if line[0:6]=='HETATM' or line[0:4]=='ATOM':
                # Extract pdb index and atom type (custom)
                idx = int(line[6:11].strip())
                atype = this_line[-1]
                atom_list.append(idx)
                molecule_graph.add_node(idx)
                # Check if this index already exists. Add if not
                if idx not in atom_parms:
                    atom_parms[idx] = {}
                    atom_parms[idx]['type'] = atype
                else:
                    repeatedIndex = True
                    if atom_parms[idx]['type'] != atype:
                        raise ValueError("PDB contains a repeated index "
                                "with different atom types: atom {} "
                                "cannot have types {} and {}".format(
                                idx, atom_parms[idx]['type'], atype))

                # Extract the element
                element = line[76:78].strip().title()
                if not element: # element is blank
                    element = 'X'
                # Check if this atomtype already exists. Add if not
                if atype not in atomtype_parms.keys():
                    atomtype_parms[atype] = {}
                    atomtype_parms[atype]['element'] = element
                    atomtype_parms[atype]['mass'] = get_mass(element)
                elif atomtype_parms[atype]['element'] != element:
                    raise ValueError("PDB contains a repeated type with "
                          "different elements: atom type {} cannot be "
                          "elements {} and {}".format(atype,
                              atomtype_parms[atype], element))
            # Read bond info
            if line[0:6] == 'CONECT':
                if repeatedIndex:
                    raise ValueError("PDB contains a repeated index. "
                        "Cannot determine bond connectivity.")
                idx = this_line[1]
                for jdx in this_line[2:]:
                    molecule_graph.add_edge(int(idx),int(jdx))

    return atom_list, atom_parms, atomtype_parms, molecule_graph

def read_native_ff(ffFile, atom_parms, atomtype_parms):
    """arguments:
    ffFile, string = filename of the forcefield file
    atomParms, dict = i: {name:, type:, element:, mass:}
returns:
    atomParms, dict = i: {name:, type:, element:, mass:, vdw:}
    bond_parms, dict = (i,j): (type, length)
    angle_parms, dict = (i,j,k): (type, [ktheta,] angle)
    dihedral_parms, dict = (i,j,k,l): (type, parms)
    improperParms, dict = (i,j,k,l): (type, kpsi, angle)
"""

    bond_parms = {}
    angle_parms = {}
    dihedral_parms = {}
    improper_parms = {}
    scaling_1_4 = {}
    global dihedral_style

    with open(ffFile) as ff:
        line = ff.readline()
        while line:
            if 'dihedraltype' in line:
                dihedral_style = line.split()[1]
                if dihedral_style.lower() == 'none':
                    dihedral_style = 'none'
                    scaling_1_4['vdw'] = 0.
                    scaling_1_4['charge'] = 0.
            elif 'vdwtype' in line:
                vdw_style = line.split()[1]
            elif 'scaling_1_4_vdw' in line:
                scaling_1_4['vdw'] = float(line.split()[1])
            elif 'scaling_1_4_charge' in line:
                scaling_1_4['charge'] = float(line.split()[1])
            elif 'nonbonded' in line:
                atype = ff.readline().strip()
                sigma = float(ff.readline().split()[1])
                epsilon = float(ff.readline().split()[1])
                if vdw_style == 'Mie':
                    repulsive_exponent = float(ff.readline().split()[1])
                    dispersive_exponent = float(ff.readline().split()[1])
                try:
                    atom_type_charge = ff.readline().split()[1] # store as string
                except:
                    atom_type_charge = None

                if atype not in atomtype_parms:
                    atomtype_parms[atype] = {}
                atomtype_parms[atype]['bondtype'] = atype # for compatibility with gromacs
                if vdw_style == 'LJ':
                    atomtype_parms[atype]['vdw'] = (vdw_style, epsilon, sigma)
                elif vdw_style == 'Mie':
                    atomtype_parms[atype]['vdw'] = (vdw_style, epsilon, sigma, repulsive_exponent, dispersive_exponent)
                atomtype_parms[atype]['charge'] = (atom_type_charge)
            elif 'bonds' in line:
                bond_type = tuple(ff.readline().split())
                distance = float(ff.readline().split()[1])
                bond_style = ff.readline().split()[1]
                if bond_style != 'fixed':
                    raise ValueError('Cassandra does not support ' + bond_style + 'bonds.')
                bond_parms[bond_type] = (bond_style, distance)
            elif 'angles' in line:
                angle_type = tuple(ff.readline().split())
                theta = float(ff.readline().split()[1])
                ktheta = ff.readline().split()[1]
                if ktheta == 'fixed':
                    angle_parms[angle_type] = ('fixed', theta)
                else:
                    ktheta = float(ktheta)
                    angle_parms[angle_type] = ('harmonic', ktheta, theta)
            elif 'dihedrals' in line:
                dihedral_type = tuple(ff.readline().split()) #atomType
                if dihedral_style == 'CHARMM':
                    a0 = float(ff.readline().split()[1])
                    a1 = float(ff.readline().split()[1])
                    delta = float(ff.readline().split()[1])
                    dihedral_parms[dihedral_type] = (dihedral_style, a0, a1, delta)
                elif dihedral_style == 'OPLS':
                    c0 = float(ff.readline().split()[1])
                    c1 = float(ff.readline().split()[1])
                    c2 = float(ff.readline().split()[1])
                    c3 = float(ff.readline().split()[1])
                    dihedral_parms[dihedral_type] = (dihedral_style, c0, c1, c2, c3)
                elif dihedral_style == 'harmonic':
                    phi = float(ff.readline().split()[1])
                    kphi = float(ff.readline().split()[1])
                    dihedral_parms[dihedral_type] = (dihedral_style, kphi, phi)
                elif dihedral_style == 'none':
                    dihedral_parms[dihedral_type] = (dihedral_style,)
                    scaling_1_4['vdw'] = 0.
                    scaling_1_4['charge'] = 0.
            elif 'impropers' in line:
                index = tuple([int(i) for i in ff.readline().split()]) #atomNumber
                psi = float(ff.readline().split()[1])
                kpsi = float(ff.readline().split()[1])
                improper_parms[index] = ('harmonic',kpsi,psi)
            elif 'charge' in line:
                data = ff.readline().split()
                index = int(data[0]) #atomNumber
                if len(data)>1:
                    atom_parms[index]['charge'] = data[1] #store as string

            # else, if the information will be provided by atom type and corrected by checkParms(), do nothing
            line = ff.readline()

    return atom_parms, atomtype_parms, bond_parms, angle_parms, dihedral_parms, improper_parms, scaling_1_4

def read_gromacs_ff_wrapper(ffFile, atom_parms, atomtype_parms):

    bond_parms = {}
    angle_parms = {}
    dihedral_parms = {}
    scaling_1_4 = {}

    vdw_style = None
    comb_rule = None

    atom_parms, atomtype_parms, bond_parms, angle_parms, dihedral_parms, scaling_1_4 = \
            _gmx_read_ff(ffFile, atom_parms, atomtype_parms, bond_parms,
                            angle_parms, dihedral_parms, vdw_style,
                            comb_rule)

    return atom_parms, atomtype_parms, bond_parms, angle_parms, dihedral_parms

def _gmx_read_ff(ffFile, atom_parms, atomtype_parms, bond_parms,
                    angle_parms, dihedral_parms, vdw_style, comb_rule):
    """arguments:
    ffFile, string = filename of the forcefield file
    atomParms, dict = i: {name:, type:, element:, mass:}
    bond_parms, dict = (i,j): (type, length)
    angle_parms, dict = (i,j,k): (type, [ktheta,] angle)
    dihedral_parms, dict = (i,j,k,l): (type, parms)
optional arguments:
    vdw_style, string = 'LJ'
    comboRule, string = identifies the LJ parameters provided
returns:
    atomParms, dict = i: {name:, type:, element:, mass:, vdw:}
    bond_parms, dict = (i,j): (type, length)
    angle_parms, dict = (i,j,k): (type, [ktheta,] angle)
    dihedral_parms, dict = (i,j,k,l): (type, parms)
"""

    IG_CONSTANT = 0.008314

    with open(ffFile) as ff:
        line = ff.readline()
        while line:
            # Each line can start with (1) # include statement,
            # (2) section header, (3) section data, (4) comment,
            # or (5) nothing. Each if is a different possibility.
            if line.startswith('#include'):
                includeFile = os.path.expanduser(line.split()[1].strip('"'))
                if not os.path.isfile(includeFile):
                    parentDir = os.path.dirname(ffFile)
                    includeFile = os.path.join(parentDir, includeFile)
                if os.path.isfile(includeFile):
                    atom_parms, atomtype_parms, bond_parms, angle_parms, \
                            dihedral_parms, scaling_1_4 = _gmx_read_ff(
                            includeFile, atom_parms, atomtype_parms,
                            bond_parms, angle_parms, dihedral_parms,
                            scaling_1_4, vdw_style, comb_rule)
                else:
                    print( 'WARNING: Topology file {} not found. '
                           'Continuing without reading file'.format(
                            includeFile))
            # Extract the section header
            elif line.startswith('['):
                section = line.strip()
                line = ff.readline()
                # Section ends with a blank line
                while line and not line.isspace():
                    if not line.startswith(';') and not line.startswith('#'):
                        data = line.split()
                        if 'default' in section:
                            # Look for function type 1, to make
                            # sure this is Lennard-Jones
                            if data[0] == '1':
                                vdw_style = 'LJ'
                            else:
                                raise ValueError('Nonbonded forcefield is '
                                    'not LJ. Cassandra only supports LJ '
                                    'interactions at this time.')
                            comb_rule = data[1]
                            scaling_1_4['vdw'] = float(data[3])
                            scaling_1_4['charge'] = float(data[4])
                        # Look for atom_parms
                        elif 'atoms' in section:
                            index = int(data[0]) #atomNumber
                            if index not in atom_parms:
                                atom_parms[index] = {}
                            atom_parms[index]['charge'] = data[6] #store as string
                        elif 'atomtypes' in section:
                            base = None
                            for i in range(len(data)):
                                if data[i] == 'A' or data[i] == 'D':
                                    base=i
                                    break
                            if base is None:
                                raise Error('ptype is not A or D. Cassandra only supports point particles.\n')
                            index = data[0] #atomType
                            if index not in atomtype_parms:
                                atomtype_parms[index] = {}
                            atomtype_parms[index]['bondtype'] = data[0] if data[1].isdigit else data[1]
                            atomtype_parms[index]['mass'] = float(data[base-2])
                            atomtype_parms[index]['charge'] = data[base-1] #store as string
                            if comb_rule == '1':
                                C6 = float(data[base+1])
                                C12 = float(data[base+2])
                                sigma = ((C12 / C6)**(1./6.)) * 10
                                epsilon = C6**2 / 4. / C12 / IG_CONSTANT
                            elif comb_rule == '2' or comb_rule == '3':
                                sigma = float(data[base+1]) * 10
                                epsilon = float(data[base+2]) / IG_CONSTANT
                            atomtype_parms[index]['vdw'] = (vdw_style, epsilon, sigma)
                        # Look for bond_parms
                        elif 'bond' in section:
                            if 'type' in section:
                                index = (data[0], data[1]) #atomTypes
                            else:
                                index = (int(data[0]), int(data[1])) #atomNumbers
                            if not any([data[2]=='1', data[2]=='2', data[2]=='3',
                                          data[2]=='4']):
                                raise Error('Cassandra only supports fixed bonds at this time.')
                            if len(data) > 3:
                                b0 = float(data[3]) * 10
                                bond_parms[index] = ('fixed', b0)
                        # Look for angle_parms
                        elif 'angle' in section:
                            if 'type' in section:
                                index = (data[0], data[1], data[2])
                            else:
                                index = (int(data[0]), int(data[1]), int(data[2]))
                            if data[3]!='1':
                                raise Error('Cassandra supports fixed or harmonic angles.')
                            if len(data) > 4:
                                theta = float(data[4])
                                ktheta = float(data[5]) / 2. / IG_CONSTANT
                                angle_parms[index] = ('harmonic', ktheta, theta)
                        # Look for dihedral_parms
                        elif 'dihedral' in section:
                            if 'type' in section:
                                index = (data[0], data[1], data[2], data[3])
                            else:
                                index = (int(data[0]), int(data[1]), int(data[2]), int(data[3]))
                            if len(data) > 5:
                                if data[4]=='1' or data[4]=='4' or data[4]=='9': #CHARMM
                                    delta = float(data[5])
                                    a0 = float(data[6])
                                    a1 = float(data[7])
                                    dihedral_parms[index] = ('CHARMM', a0, a1, delta)
                                elif data[4] == '2': #harmonic
                                    phi = float(data[5])
                                    kphi = float(data[6]) / 2. / IG_CONSTANT
                                    dihedral_parms[index] = ('harmonic', kphi, phi)
                                elif data[4] == '3': #Ryckaert-Bellemans
                                    c0 = float(data[5])
                                    c1 = float(data[6])
                                    c2 = float(data[7])
                                    c3 = float(data[8])
                                    c4 = float(data[9])
                                    c5 = float(data[10])
                                    a0 = c0 + c1 + c2 + c3
                                    a1 = - c1 - 3. * c3 / 4.
                                    a2 = - c2 / 2.
                                    a3 = - c3 / 4.
                                    if not c4 == 0. and not c5 == 0.:
                                        raise ValueError('Can only convert '
                                             'Ryckaert-Bellemans dihedrals '
                                             'to OPLS if c4==0 and c5==0.\n')
                                    dihedral_parms[index] = ('OPLS', a0, a1, a2, a3)
                                elif data[4] == '5': #OPLS
                                    a0 = 0.
                                    a1 = float(data[5]) / 2.
                                    a2 = float(data[6]) / 2.
                                    a3 = float(data[7]) / 2.
                                    dihedral_parms[index] = ('OPLS', a0, a1, a2, a3)
                    line = ff.readline()
            line = ff.readline()

    ff.close()

    return atom_parms, atomtype_parms, bond_parms, angle_parms, dihedral_parms, scaling_1_4

def check_atom_parms(atom_list, atom_parms, atomtype_parms):

    for pdb_idx in atom_list:
        for parm in ['vdw', 'charge', 'element', 'mass', 'bondtype']:
            if parm not in atom_parms[pdb_idx]:
                atype = atom_parms[pdb_idx]['type']
                if parm in atomtype_parms[atype]:
                    atom_parms[pdb_idx][parm] = atomtype_parms[atype][parm]
                else:
                    raise ValueError('{} parms for atom {} not '
                        'found.'.format(parm,pdb_idx))

    return atom_parms

def check_bond_parms(bond_list, bond_parms, atom_parms):

    for bond in bond_list:
        ij = bond
        ji = ij[::-1]
        if ij not in bond_parms and ji not in bond_parms:
            ijType = tuple([atom_parms[i]['bondtype'] for i in ij])
            jiType = ijType[::-1]
            if ijType in bond_parms:
                bond_parms[ij] = bond_parms[ijType]
            elif jiType in bond_parms:
                bond_parms[ij] = bond_parms[jiType]
            else:
                raise ValueError('Bond parms for atoms {} not '
                    'found.'.format(bond))

    return bond_parms

def check_angle_parms(angle_list, angle_parms, atom_parms):

    for angle in angle_list:
        ijk = angle
        kji = ijk[::-1]
        if ijk not in angle_parms and kji not in angle_parms:
            ijkType = tuple([atom_parms[i]['bondtype'] for i in ijk])
            kjiType = ijkType[::-1]
            if ijkType in angle_parms:
                angle_parms[ijk] = angle_parms[ijkType]
            elif kjiType in angle_parms:
                angle_parms[ijk] = angle_parms[kjiType]
            else:
                raise ValueError('Angle parms for atoms {} not '
                    'found.'.format(angle))

    return angle_parms

def check_dihedral_parms(dihedral_list, dihedral_parms, atom_parms):

    for dihedral in dihedral_list:
        ijkl = dihedral
        lkji = ijkl[::-1]
        if ijkl not in dihedral_parms and lkji not in dihedral_parms:
            ijklType = tuple([atom_parms[i]['bondtype']for i in ijkl])
            lkjiType = ijklType[::-1]
            jkDefault = tuple(['X', ijklType[1], ijklType[2], 'X'])
            kjDefault = tuple(['X', lkjiType[1], lkjiType[2], 'X'])
            if ijklType in dihedral_parms:
                dihedral_parms[ijkl] = dihedral_parms[ijklType]
            elif lkjiType in dihedral_parms:
                dihedral_parms[ijkl] = dihedral_parms[lkjiType]
            elif jkDefault in dihedral_parms:
                dihedral_parms[ijkl] = dihedral_parms[jkDefault]
            elif kjDefault in dihedral_parms:
                dihedral_parms[ijkl] = dihedral_parms[kjDefault]
            else:
                raise ValueError('Dihedral parms for atoms {} not '
                    'found.'.format(dihedral))

    return dihedral_parms

def write_mcf_atoms(mcf, atom_list, atom_parms, in_ring):
    """Write the atoms in the system.
    """

    header = ('!Atom Format\n'
              '!index type element mass charge vdw_style parameters\n'
              '!vdw_style="LJ", parms=epsilon sigma\n'
              '!vdw_style="Mie", parms=epsilon sigma '
              'repulsion_exponent dispersion_exponent\n'
              '\n# Atom_Info\n'
              )
    mcf.write(header)
    mcf.write(str(len(atom_list))+'\n')
    for mcf_idx, pdb_idx in enumerate(atom_list):
        mcf.write('{:<4d}'.format(mcf_idx+1))
        mcf.write('  {:<6s}'.format(atom_parms[pdb_idx]['type']))
        mcf.write('  {:<2s}'.format(atom_parms[pdb_idx]['element']))
        mcf.write('  {:7.3f}'.format(atom_parms[pdb_idx]['mass']))
        mcf.write('  {:s}'.format(atom_parms[pdb_idx]['charge']))
        mcf.write('  {:3s}'.format(atom_parms[pdb_idx]['vdw'][0]))
        for i in range(1,len(atom_parms[pdb_idx]['vdw'])):
            mcf.write('  {:8.3f}'.format(atom_parms[pdb_idx]['vdw'][i]))
        if in_ring[pdb_idx] == True:
            mcf.write('  ring')
        mcf.write('\n')

def write_mcf_bonds(mcf, atom_list, bond_list, bond_parms):
    """Write the bonds in the system.
    """

    header = ( '\n!Bond Format\n'
               '!index i j type parameters\n' 
               '!type="fixed", parms=bondLength\n'
               '\n# Bond_Info\n'
             )
    mcf.write(header)
    mcf.write('{:d}\n'.format(len(bond_list)))

    for ctr, bond in enumerate(bond_list):
        if bond not in bond_parms:
            bond = bond[::-1]
        ij = tuple([atom_list.index(idx)+1 for idx in bond]) #mcf indices
        mcf.write('{:<4d}'.format(ctr+1))
        mcf.write('  {:<4d}  {:<4d}'.format(ij[0],ij[1]))
        mcf.write('  {:5s}  {:8.3f}'.format(bond_parms[bond][0],
                                            bond_parms[bond][1]))
        mcf.write('\n')

def write_mcf_angles(mcf, atom_list, angle_list, angle_parms):
    """Write the angles in the system.
    """

    header = ( '\n!Angle Format\n'
               '!index i j k type parameters\n'
               '!type="fixed", parms=equilibrium_angle\n'
               '!type="harmonic", parms=force_constant equilibrium_angle\n'
               '\n# Angle_Info\n'
             )
    mcf.write(header)
    mcf.write('{:d}\n'.format(len(angle_list)))

    for ctr,angle in enumerate(angle_list):
        if angle not in angle_parms:
            angle = angle[::-1]
        ijk = tuple([atom_list.index(idx)+1 for idx in angle])
        mcf.write('{:<4d}'.format(ctr+1))
        mcf.write('  {:<4d}  {:<4d}  {:<4d}'.format(ijk[0],ijk[1],ijk[2]))
        if angle_parms[angle][0] == 'harmonic':
            mcf.write('  {:<8s}  {:8.1f}  {:8.2f}'.format(angle_parms[angle][0],
                                                          angle_parms[angle][1],
                                                          angle_parms[angle][2]))
        elif angle_parms[angle][0] == 'fixed':
            mcf.write('  {:<8s}  {:8.2f}'.format(angle_parms[angle][0],
                                                 angle_parms[angle][1]))
        mcf.write('\n')

def write_mcf_dihedrals(mcf, atom_list, dihedral_list, dihedral_parms):
    """Write the dihedrals in the system.
    """

    header = ( '\n!Dihedral Format\n'
               '!index i j k l type parameters\n'
               '!type="none"\n'
               '!type="CHARMM", parms=a0 a1 delta\n'
               '!type="OPLS", parms=c0 c1 c2 c3\n'
               '!type="harmonic", parms=force_constant equilibrium_dihedral\n'
               '\n# Dihedral_Info\n'
             )
    mcf.write(header)
    mcf.write('{:d}\n'.format(len(dihedral_list)))

    for ctr, dihedral in enumerate(dihedral_list):
        if dihedral not in dihedral_parms:
            dihedral = dihedral[::-1]
        ijkl = tuple([atom_list.index(idx)+1 for idx in dihedral])
        mcf.write('{:<4d}'.format(ctr+1))
        mcf.write('  {:<4d}  {:<4d}  {:<4d}  {:<4d}'.format(ijkl[0],ijkl[1],
                                                            ijkl[1],ijkl[2]))
        if dihedral_parms[dihedral][0] == 'CHARMM':
            mcf.write('  {:<8s}  {:8.3f}  {:8.1f}  {:8.1f}'.format(
                                                dihedral_parms[dihedral][0],
                                                dihedral_parms[dihedral][1],
                                                dihedral_parms[dihedral][2],
                                                dihedral_parms[dihedral][3]))
        elif dihedral_parms[dihedral][0] == 'OPLS':
            mcf.write('  {:<8s}  {:8.3f}  {:8.3f}  {:8.3f}  {:8.3f}'.format(
                                                dihedral_parms[dihedral][0],
                                                dihedral_parms[dihedral][1],
                                                dihedral_parms[dihedral][2],
                                                dihedral_parms[dihedral][3],
                                                dihedral_parms[dihedral][4]))
        elif dihedral_parms[dihedral][0] == 'harmonic':
            mcf.write('  {:<8s}  {:8.1f}  {:8.2f}'.format(
                                                dihedral_parms[dihedral][0],
                                                dihedral_parms[dihedral][1],
                                                dihedral_parms[dihedral][2]))
        elif dihedral_parms[dihedral][0] == 'none':
            mcf.write('  {:<8s}'.format(dihedral_parms[dihedral]))
        mcf.write('\n')

def write_mcf_impropers(mcf, atom_list, improper_parms):
    """Write the impropers in the system.
    """

    header = ( '\n!Improper Format\n'
               '!index i j k l type parameters\n'
               '!type="harmonic", parms=force_constant equilibrium_improper\n'
               '\n# Improper_Info\n'
             )
    mcf.write(header)
    mcf.write('{:d}\n'.format(len(improper_parms)))

    for improper_idx, improper in enumerate(improper_parms.items()):
        ijkl = tuple([atom_list.index(idx)+1 for idx in improper])
        mcf.write('{:<4d}'.format(improper_idx+1))
        mcf.write('  {:<4d}  {:<4d}  {:<4d}  {:<4d}'.format(ijkl))
        mcf.write('  {:<8s}  {:8.1f}  {:8.2f}'.format(improper_parms[improper]))
        mcf.write('\n')

def write_mcf_fragments(mcf, atom_list, frag_list, frag_conn):
    """Write the fragments in the molecule.
    """

    header = ('\n!Fragment Format\n'
              '!index number_of_atoms_in_fragment branch_point other_atoms\n'
              '\n# Fragment_Info\n'
             )
    mcf.write(header)

    # Special cases first
    if len(frag_list) == 0:
        if len(atom_list) == 1:
            mcf.write('1\n')
            mcf.write('1 1 1\n')
        elif len(atom_list) == 2:
            mcf.write('1\n')
            mcf.write('1 2 1 2\n')
        else:
            warnings.warn('More than two atoms present but '
                          'no fragments identified.')
            mcf.write('0\n')
    else:
        mcf.write('{:d}\n'.format(len(frag_list)))
        for i, frag in enumerate(frag_list):
            mcf.write('{:d}    {:d}'.format(i+1, len(frag)))
            for idx in frag:
                mcf.write('    {:d}'.format(idx))
            mcf.write('\n')

    mcf.write('\n\n# Fragment_Connectivity\n')
    mcf.write('{:d}\n'.format(len(frag_conn)))
    for i, conn in enumerate(frag_conn):
        mcf.write('{:d}    {:d}    {:d}\n'.format(i+1, conn[0]+1, conn[1]+1))

def write_mcf_intra_scaling(mcf,scaling_1_4):
    """Write the intramolecular scaling in the molecule.
    """
    header = ( '\n!Intra Scaling\n'
               '!vdw_scaling    1-2 1-3 1-4 1-N\n'
               '!charge_scaling 1-2 1-3 1-4 1-N\n'
               '\n# Intra_Scaling\n'
             )

    mcf.write(header)
    mcf.write('0. 0. {:.4f} 1.\n'.format(scaling_1_4['vdw']))
    mcf.write('0. 0. {:.4f} 1.\n'.format(scaling_1_4['charge']))

def check_config_file(configFile):

    basename = os.path.splitext(os.path.basename(configFile))[0]
    configFileType = check_type_configfile(configFile)
    if configFileType == 'cml':
        cml_to_pdb(configFile)
        configFile = basename + '.pdb'

    print("\nThe configFile format was determined "
          "to be {}\n".format(configFileType))

    if not os.path.isfile(configFile):
        raise FileNotFoundError('Could not find file '
                '{}'.format(configFile))

    return configFile, configFileType, basename

def check_mcf_file(mcfFile,basename,ffTemplate):

    if mcfFile is None:
        mcfFile = basename + '.mcf'

    if not ffTemplate:
        if os.path.isfile(mcfFile):
            overwrite = input('MCF file {} already exists. Would you like '
                              'to overwrite the existing file?'
                              '(y/n):'.format(mcfFile))
            if overwrite.lower() == 'n' or overwrite.lower() == 'no':
                exit()
            elif overwrite.lower() == 'y' or overwrite.lower() == 'yes':
                bak_ctr = 1
                while os.path.isfile(mcfFile + '.bak' + str(bak_ctr)):
                    bak_ctr += 1
                    if bak_ctr > 100:
                        raise ValueError('There seems to be a problem.')
                print('Just in case we are creating a backup of '
                      '{} at {}\n'.format(mcfFile, mcfFile +
                      '.bak' + str(bak_ctr)))
                os.system('mv ' + mcfFile + ' ' + mcfFile + '.bak' + str(bak_ctr))
            else:
                raise ValueError('Invalid selection.')

    return mcfFile

def check_ff_file(ffFile,basename,ffTemplate):

    if ffFile is None:
        ffFile = basename + '.ff'
        ffFileType = 'native'
    else:
        ffFileExt = os.path.splitext(ffFile)[1]
        if ffFileExt == '.ff':
            ffFileType = 'native'
        elif ffFileExt == '.itp' or ffFileExt == '.top':
            ffFileType = 'gromacs'
        else:
            raise ValueError("Unsupported forcefield (ff) file "
                    "type: {}".format(ffFileType))

    if ffTemplate:
        if os.path.isfile(ffFile):
            overwrite = input('FF file {} already exists. Would you like '
                              'to overwrite the existing file?'
                              '(y/n):'.format(ffFile))
            if overwrite.lower() == 'n' or overwrite.lower() == 'no':
                exit()
            elif overwrite.lower() == 'y' or overwrite.lower() == 'yes':
                bak_ctr = 1
                while os.path.isfile(ffFile + '.bak' + str(bak_ctr)):
                    bak_ctr += 1
                    if bak_ctr > 100:
                        raise ValueError('There seems to be a problem.')
                print('Just in case we are creating a backup of '
                      '{} at {}\n'.format(ffFile, ffFile +
                      '.bak' + str(bak_ctr)))
                os.system('mv ' + ffFile + ' ' + ffFile + '.bak' + str(bak_ctr))
            else:
                raise ValueError('Invalid selection.')
    else:
        if not os.path.isfile(ffFile):
            raise FileNotFoundError('Could not find file '
                    '{}'.format(ffFile))

    return ffFile, ffFileType

def check_type_configfile(infilename):
    file = open(infilename,'r')
    for line in file:
        if '<molecule>' in line:
            return 'cml'
    return 'pdb'

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
        cml_atom_info.append(re.findall('"([^"]*)"',
                             linecache.getline(infilename, i)))

    for i in xrange(cml_start_bonds+1, cml_end_bonds):
        cml_bond_info.append(re.findall('"([^"]*)"',
                             linecache.getline(infilename, i))[0].split())


    filePdb = open(basename+'.pdb','w')

    for line_nbr, line in enumerate(cml_atom_info):
        pdbFormat='%-6s%5d %-4s%60s%2s  %s\n'
        recordName = 'HETATM'
        atomNumber = line_nbr + 1
        atomElement = line[1]
        atomName = atomElement + str(atomNumber)
        atomType = line[5]
        filePdb.write(pdbFormat % (recordName, atomNumber, atomName, '',
                                    atomElement, atomType))

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
        filePdb.write('CONECT     '+line[0]+'    ')
        for element in line[1]:
            filePdb.write(element+"   ")
        filePdb.write("\n")

    filePdb.close()

def get_mass(element):

    mass_dict = {
        'H': 1.0079,
        'He': 4.0026,
        'Li': 6.941,
        'Be': 9.0122,
        'B': 10.811,
        'C': 12.0107,
        'N': 14.0067,
        'O': 15.9994,
        'F': 18.9984,
        'Ne': 20.1797,
        'Na': 22.9897,
        'Mg': 24.305,
        'Al': 26.9815,
        'Si': 28.0855,
        'P': 30.9738,
        'S': 32.065,
        'Cl': 35.453,
        'Ar': 39.948,
        'K': 39.0983,
        'Ca': 40.078,
        'Sc': 44.9559,
        'Ti': 47.867,
        'V': 50.9415,
        'Cr': 51.9961,
        'Mn': 54.938,
        'Fe': 55.845,
        'Co': 58.9332,
        'Ni': 58.6934,
        'Cu': 63.546,
        'Zn': 65.39,
        'Ga': 69.723,
        'Ge': 72.64,
        'As': 74.9216,
        'Se': 78.96,
        'Br': 79.904,
        'Kr': 83.8,
        'Rb': 85.4678,
        'Sr': 87.62,
        'Y': 88.9059,
        'Zr': 91.224,
        'Nb': 92.9064,
        'Mo': 95.94,
        'Tc': 98,
        'Ru': 101.07,
        'Rh': 102.9055,
        'Pd': 106.42,
        'Ag': 107.8682,
        'Cd': 112.411,
        'In': 114.818,
        'Sn': 118.71,
        'Sb': 121.76,
        'Te': 127.6,
        'I': 126.9045,
        'Xe': 131.293,
        'Cs': 132.9055,
        'Ba': 137.327,
        'La': 138.9055,
        'Ce': 140.116,
        'Pr': 140.9077,
        'Nd': 144.24,
        'Pm': 145,
        'Sm': 150.36,
        'Eu': 151.964,
        'Gd': 157.25,
        'Tb': 158.9253,
        'Dy': 162.5,
        'Ho': 164.9303,
        'Er': 167.259,
        'Tm': 168.9342,
        'Yb': 173.04,
        'Le': 174.967,
        'Hf': 178.49,
        'Ta': 180.9479,
        'W': 183.84,
        'Re': 186.207,
        'Os': 190.23,
        'Ir': 192.217,
        'Pt': 195.078,
        'Au': 196.9665,
        'Hg': 200.59,
        'Tl': 204.3833,
        'Pb': 207.2,
        'Bi': 208.9804,
        'Po': 209,
        'At': 210,
        'Rn': 222,
        'Fr': 223,
        'Ra': 226,
        'Ac': 227,
        'Th': 232.0381,
        'Pa': 231.0359,
        'U': 238.0289,
        'Np': 237,
        'Pu': 244,
        'Am': 243,
        'Cm': 247,
        'Bk': 247,
        'Cf': 251,
        'Es': 252,
        'Fm': 257,
        'Md': 258,
        'No': 259,
        'Lr': 262,
        'Rf': 261,
        'Db': 262,
        'Sg': 266,
        'Bh': 264,
        'Hs': 277,
        'Mt': 268,
        'C1': 13.0186,
        'C2': 14.0265,
        'C3': 15.0344,
        'C4': 16.0423,
        'CH': 13.0186,
        'CH2': 14.0265,
        'CH3': 15.0344,
        'CH4': 16.0423
    }

    return mass_dict[element]

def parse_args():
    parser = argparse.ArgumentParser(
    formatter_class=argparse.RawDescriptionHelpFormatter, description=
    """DESCRIPTION:
    To generate a Molecular Connectivity File (.mcf), you will need both a
    configuration file and a file with the force field parameters.

    EXAMPLES
    To generate a .mcf from forcefield parameters in the literature:

      1) Run
           > python mcfgen.py molecule.pdb --ffTemplate
         to generate an intermediate file molecule.ff.
      2) Fill out molecule.ff with the parameters found in the literature.
      3) Run
           > python mcfgen.py molecule.pdb

    To generate a .mcf from a GROMACS forcefield file

      1) Run
           > python mcfgen.py molecule.pdb -f molecule.itp

    NOTES:
      1) The Protein Data Bank format (.pdb) defines a fixed 80-character format for
         HETATM and ATOM records. Entries are not necessarily separated by white
         space. The following information is read from each HETATM/ATOM record in
         the supplied .pdb file:

         COLUMNS        DATA TYPE       CONTENTS
         ---------------------------------------------------------------------------
          7 - 11        Integer(5)      Atom serial number.

         31 - 38        Real(8.3)       Orthogonal coordinates for X in Angstroms.

         39 - 46        Real(8.3)       Orthogonal coordinates for Y in Angstroms.

         47 - 54        Real(8.3)       Orthogonal coordinates for Z in Angstroms.

         77 - 78        RString(2)      Element symbol, right-justified.
         ---------------------------------------------------------------------------

         The element name is used to look up the atomic mass. If the element field
         is blank or contains an unknown element, the user will be prompted for the
         atomic mass, unless the --massDefault option is used, in which case the
         unknown element will given a mass of MASSDEFAULT.

         In addition to the standard element symbols, the following pseudo-atoms
         are recognized:

         PSEUDO-ATOM         SYMBOL      MASS
         ---------------------------------------------------------------------------
         Methane, CH4        C4          16.0423
         Methyl, -CH3        C3          15.0344
         Methanediyl, -CH2-  C2          14.0265
         Methanetriyl, -CH<  C1          13.0186
         ---------------------------------------------------------------------------

         The atom type for each HETATM/ATOM record is taken from the last white-
         space separated column. This may be the element name, or the user may need
         to append the atom type to the end of the record.

         After the HETATM/ATOM records, molecular connectivity is determined from a
         series of CONECT records. If no CONECT records are given for a poly-atomic
         molecule, the molecule will be listed with zero fragments. A molecule with
         zero fragments cannot be inserted, deleted or regrown in Cassandra. CONECT
         records have the following format:

         COLUMNS         DATA TYPE        CONTENTS
         ---------------------------------------------------------------------------
          7 - 11         Integer(5)       Atom serial number

         12 - 16         Integer(5)       Serial number of bonded atom

         17 - 21         Integer(5)       Serial number of bonded atom

         22 - 26         Integer(5)       Serial number of bonded atom

         27 - 31         Integer(5)       Serial number of bonded atom
         ---------------------------------------------------------------------------

         If the atom is bonded to more than 4 other atoms, the pattern is continued.
         The CONECT record is parsed as white-space separated data.

      2) This script does not currently support multiple dihedrals for the same 4
         atoms. If mulitple parameters are given for the same dihedral sequence,
         only the last parameters will be used.

      3) Improper definitions are not currently read from Gromacs forcefield files.
    """)
    parser.add_argument('configFile',
                        help="""CONFIGFILE must be in either .pdb or .cml format. """ +
                             """A .pdb file can be generated using Gaussview, """ +
                             """while .cml files can be generated using Avogadro.""")
    parser.add_argument('--ffTemplate', '-t', action='store_true',
                        help="""Generate a blank force field file template.""")
    parser.add_argument('--ffFile', '-f',
                        help="""The default FFFILE is molecule.ff, a custom format """ +
                             """for this script. Alternatively, the forcefile """ +
                             """parms can be supplied in GROMACS format (.itp).""")
    parser.add_argument('--mcfFile', '-m',
                        help="""The default MCFFILE is molecule.mcf.""")
    parser.add_argument('--massDefault', '-d', type=float,
                        help="""Provide a default mass for unknown elements.""")
    parser.add_argument('--verbose', '-v', action='store_true',
                        help="""Print an atom-by-atom summary of the topology scan""")

    args = parser.parse_args()

    return args

if __name__ == "__main__":
    main()
