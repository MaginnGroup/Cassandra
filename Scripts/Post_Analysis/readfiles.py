# Cassandra Classes
# This script contains all the functions required to read in the input/mcf files.
# It gets called from the highest level scripts.

from molecule_info import *

def read_atom_info(mcf,natoms):
	i = 0
	molecule = molecule_info()
	for line in mcf:
	
		atom_name = line.split()[1] #atom name
		atom_type = line.split()[2] #atom type
		mass = line.split()[3] #charge
		charge = line.split()[4]
		eps = line.split()[6] #epsilon
		sigma = line.split()[7] #sigma
		
		molecule.add_atom_info(natoms,atom_name, atom_type, sigma, eps, mass, charge)
		molecule.calc_MW(mass)	
		i+=1
		if (i == natoms):
			return molecule

def read_angle_info(mcf,molecule,nangles):
	i = 0

	if (nangles == 0):
		return molecule	
	else:
		for line in mcf:
			atom1_ndx = line.split()[1]
			atom2_ndx = line.split()[2]
			atom3_ndx = line.split()[3]	
			molecule.add_angle_info(atom1_ndx, atom2_ndx, atom3_ndx)
			i+=1
			if (i == nangles):
				return molecule

def read_dihedral_info(mcf,molecule,ndihedrals):
	i = 0
	if (ndihedrals == 0):
		return molecule

	else:
		for line in mcf:
			atom1_ndx = line.split()[1]
			atom2_ndx = line.split()[2]
			atom3_ndx = line.split()[3]
			atom4_ndx = line.split()[4]
			molecule.add_dihedral_info(atom1_ndx, atom2_ndx, atom3_ndx, atom4_ndx)
			i+=1
			if (i == ndihedrals):
				return molecule

def read_mcf_file(mcffile):
	with open(mcffile,'r') as mcf:
		for line in mcf:
			if "# Atom_Info" in line:
				line = mcf.next()
				natoms = int(line.split()[0])
				molecule = read_atom_info(mcf,natoms)
				mcffile = mcffile.strip('./')
			elif "# Angle_Info" in line:
				line = mcf.next()
				nangles = int(line.split()[0])
				molecule.nangles = nangles
				molecule = read_angle_info(mcf,molecule,nangles)
			elif "# Dihedral_Info" in line:
				line = mcf.next()
				ndihedrals = int(line.split()[0])
				molecule.ndihedrals = ndihedrals
				molecule = read_dihedral_info(mcf,molecule,ndihedrals)

	molecule.add_molecule_name(mcffile.replace('.mcf',''))
	return molecule

def read_H_file(Hfile,nspecies):
	Lx_array = []
	Ly_array = []
	Lz_array = []
	nmolecules_array =[]
	with open(Hfile,'r') as H:
		for line in H:
			line = H.next()
			Lx_array.append(float(line.split()[0]))
			line = H.next()
			Ly_array.append(float(line.split()[1]))
			line = H.next()
			Lz_array.append(float(line.split()[2]))
			line = H.next()
			line = H.next()
			for i in range(nspecies):
				line = H.next()
				nmolecules_array.append(line.split()[1])
		return (Lx_array,Ly_array,Lz_array,nmolecules_array)


def read_inp_file(inpfile):
	nmolecules_array = []
	with open(inpfile,'r') as inp:
		for line in inp:
			if "# Sim_Type" in line:
				line = inp.next()
				sim_type = line.split()[0]
	with open(inpfile,'r') as inp:			
		for line in inp:
			if "# Nbr_Species" in line:
				line = inp.next()
				nspecies = int(line.split()[0])
	with open(inpfile,'r') as inp:
		for line in inp:
			if "# Molecule_Files" in line:
				print nspecies
				for i in range(nspecies):
					line = inp.next()
					nmolecules_array.append(int(line.split()[1]))
	
	return (sim_type,nspecies,nmolecules_array)
				
