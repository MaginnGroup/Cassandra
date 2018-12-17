#molecule_info.py


# This library assumes an MCF file has been read and parsed.
# We will store the molecule info into this class.


import numpy as np

class molecule_info:
	def __init__(self):
		self.natoms = 0
		self.atom_name = []
		self.atom_type = []
		self.sigma = []
		self.eps = []
		self.mass = []
		self.charge = []
		self.name = '' 
		self.MW = 0
		self.com = 0
		self.nangles = 0
		self.angles = []
		self.ndihedrals = 0
		self.dihedrals = []

	def add_atom_info(self,number_atoms, atom_name,atom_type,sigma,eps,mass,charge):
		self.natoms =number_atoms
		self.atom_name.append(atom_name)
		self.atom_type.append(atom_type)
		self.sigma.append(sigma)
		self.eps.append(eps)
		self.mass.append(mass)
		self.charge.append(charge)

	def calc_MW(self,mass):
		self.MW = self.MW +float(mass)

	def add_molecule_name(self,name):
		self.name = name


	def add_angle_info(self,atom1_ndx,atom2_ndx,atom3_ndx):
		self.angles.append([atom1_ndx,atom2_ndx,atom3_ndx])

	def add_dihedral_info(self,atom1_ndx, atom2_ndx, atom3_ndx, atom4_ndx):
		self.dihedrals.append([atom1_ndx,atom2_ndx,atom3_ndx,atom4_ndx])

