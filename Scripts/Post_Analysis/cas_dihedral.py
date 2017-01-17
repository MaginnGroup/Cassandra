#!add path to python here

import numpy as np
import argparse
from argparse import RawTextHelpFormatter
import matplotlib.pyplot as plt

# F2PY Shared Objects
import cas_palib
# Cassandra Classes
from molecule_info import *
import readfiles 

#DEFINE FUNCTIONS

def obtain_dihedral_indices(output_ndx):
#!Function will determine species and atom index given a cumulative atom type index
#!Used for atom-atom RDF
	dihedral_ndx = 0
	for i in range(nspecies):
		for j in range(molecule_info[i].ndihedrals):
			if output_ndx <= dihedral_ndx:
				species_ndx = i
				atom1_ndx = j
				return(i,j)
			dihedral_ndx+=1


#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
parser = argparse.ArgumentParser('''*** CASSANDRA Dihedral Distribution Function Analysis Tool ***

Example Usage:

	python cas_dihedral.py -f test.xyz -m s1.mcf s2.mcf (-deg 1) (-o dihedral.xvg)

''',formatter_class=RawTextHelpFormatter)
parser.add_argument('-f',action='store',help ='Trajectory File (xyz format)')
parser.add_argument('-m',nargs = '+',action = 'store', help='MCF file(s) for each species. Must be in correct order.')
parser.add_argument('-H',action='store',help='Optional: H Matrix File (specify only if file name without extension is different from xyz file name)',default = 'default')
parser.add_argument('-b',action='store',help='Optional: Beginning frame (default = 1)', default = 1, type=int)
parser.add_argument('-e',action='store',help='Optional: Ending frame (default=last frame)', default=0,type=int)
parser.add_argument('-deg',action='store',help='Optional: degree per bin (default=1 degree/bin)',default=1,type=float)
parser.add_argument('-o',action='store',help='Optional: Output file name',default='dihedral.xvg')

#DEFINE PARSER ARGS
args = parser.parse_args()
nspecies = len(args.m)
binwidth = args.deg
nslices = int(360.0/(args.deg))
xyzfile = args.f

if args.H == 'default':
	Hfile = xyzfile[:-4]+'.H'
else:
	Hfile = args.H
begin_frame = args.b-1
end_frame = args.e-1
outfile = args.o
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#INITIALIZE
mcffiles = []
molecule_info = []
natoms=[]
max_natoms = 0
max_nmols = 0
species_ndx = 0
atom1_ndx = 0
atom2_ndx = 0
atom3_ndx = 0
atom4_ndx = 0
atom1_type = "temp"
atom2_type = "temp"
atom3_type = "temp"
atom4_type = "temp"
dihedral = np.zeros(nslices+1,np.float)

#READ MCF FILES
#create an object of molecule_info class for each species containing:
#        self.natoms = int
#        self.atom_name = array
#        self.atom_type = array
#        self.sigma = array
#        self.eps = array
#        self.mass = array
#        self.charge = array
#        self.name = char
#        self.MW = float

for i in range(nspecies):
	mcffile = args.m[i]
	mcffiles.append(mcffile)
	molecule_info.append(readfiles.read_mcf_file(mcffile))


# READ H-MATRIX FILE
#Read in the Hmatrix file to obtain the box lengths and the number of molecules
#for each snapshot. 
#Note: Matrices are flattened
H = readfiles.read_H_file(Hfile,nspecies)
Lx = H[0]
Ly = H[1]
Lz = H[2]

#Store number of frames
nframes = int(len(H[3])/nspecies)
if end_frame == -1:
	end_frame = nframes

if (begin_frame > nframes) or (end_frame > nframes):
	print "Error: Beginning or ending frame is larger than the total number of frames ("+str(nframes)+")"
	quit()
elif (begin_frame < 0) or (end_frame < 0):
	print "Error: Specified frames are outside of bounds"
	quit()

#Take the binwidth from the the smallest Length
Lmin = np.min([np.min(Lx[begin_frame:end_frame]),np.min(Ly[begin_frame:end_frame]),
		np.min(Lz[begin_frame:end_frame])])


#Construct matrix that contains number of each species for each frame
nmolecules = np.reshape(H[3],(nframes,nspecies))


# (OPTIONAL) READ INPUT FILE
#out = readfiles.read_inp_file(inpfile)
#molecule = readmcf.read_mcf_file(mcffile)

#Check MCF Files are consistent with first frame


#Store number of atoms of each species
for i in range(nspecies):
	natoms.append(molecule_info[i].natoms)

#Construct matrix (# species by max # atomtype) of masses of each atom type

#Find max number of atoms in a species
for i in range(nspecies):
	if max_natoms < molecule_info[i].natoms:
		max_natoms = molecule_info[i].natoms

#Find max number of molecules in a species
for i in range(nframes):
	for j in range(nspecies):
		if max_nmols < nmolecules[i][j]:
			max_nmols = nmolecules[i][j]


#Initialize mass for each atom in the system
atype_mass = np.zeros((max_natoms,nspecies))

#Loop through species and fill masses of each atom type
for i in range(nspecies):
	for j in range(max_natoms):
		try:
			atype_mass[j][i] = float(molecule_info[i].mass[j]) #convert to kg
		except:
			pass


ia = 0
print "Provide the index of the dihedral you would like to analyze:\n"
for i in range(nspecies):
	print "Species",str(i+1)+":\n"
	if molecule_info[i].ndihedrals == 0:
		print "Species",str(i+1),"has no dihedrals."

	for j in range(molecule_info[i].ndihedrals):
		#temporarily declare atom indices and names 
		#(these will be overwritten once the user provides an index for the output)
		atom1_ndx = int(molecule_info[i].dihedrals[j][0])-1
		atom1_type = molecule_info[i].atom_type[atom1_ndx]
		atom2_ndx = int(molecule_info[i].dihedrals[j][1])-1
		atom2_type = molecule_info[i].atom_type[atom2_ndx]
		atom3_ndx = int(molecule_info[i].dihedrals[j][2])-1
		atom3_type = molecule_info[i].atom_type[atom3_ndx]
		atom4_ndx = int(molecule_info[i].dihedrals[j][3])-1
		atom4_type = molecule_info[i].atom_type[atom4_ndx]
		print "("+str(ia+1)+"):",atom1_type,atom2_type, atom3_type, atom4_type,\
			" (atom #: "+str(atom1_ndx+1)+"-"+str(atom2_ndx+1)+"-"+\
			str(atom3_ndx+1)+"-"+str(atom4_ndx+1)+")"
		ia += 1
	print "\n"

#obtain index for output from user
output_ndx = int(raw_input())-1

#need to determine species and dihedral index
species_ndx, dihedral_ndx = obtain_dihedral_indices(output_ndx)

#store atom indexes and names based on user input
atom1_ndx = int(molecule_info[species_ndx].dihedrals[dihedral_ndx][0])-1
atom1_type = molecule_info[species_ndx].atom_type[atom1_ndx]
atom2_ndx = int(molecule_info[species_ndx].dihedrals[dihedral_ndx][1])-1
atom2_type = molecule_info[species_ndx].atom_type[atom2_ndx]
atom3_ndx = int(molecule_info[species_ndx].dihedrals[dihedral_ndx][2])-1
atom3_type = molecule_info[species_ndx].atom_type[atom3_ndx]
atom4_ndx = int(molecule_info[species_ndx].dihedrals[dihedral_ndx][3])-1
atom4_type = molecule_info[species_ndx].atom_type[atom4_ndx]

print "\n"
print "Dihedral",str(output_ndx+1),"(species",str(species_ndx+1)+") has been selected.\n"
print "Atoms",atom1_type+",",atom2_type+",",atom3_type+",","and",atom4_type,"will be analyzed.\n\n"

#OBTAIN THE INDICES OF SPECIES OR ATOM PAIR OUTPUT


##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
dihedral = cas_palib.cas_dihedral(xyzfile,nspecies, nslices,
						nframes, begin_frame, end_frame,
						Lx,Ly, Lz, 
						species_ndx,  
						atom1_ndx, atom2_ndx, atom3_ndx, atom4_ndx,
						atom1_type, atom2_type, atom3_type, atom4_type,
						atype_mass, max_natoms,max_nmols,
						natoms, nmolecules,dihedral)

bins = np.arange(0,binwidth*(nslices+1),binwidth)-180


outfmt = '%12.4f%12.4f\n'
soutfmt = '%12.4f'
if (outfile == 'dihedral.xvg'):
	outfile = outfile[:-4]+'_'+atom1_type+'_'+atom2_type+'_'+atom3_type+'_'+atom4_type+'.xvg'

outfilew = open(outfile,'w')

for i in range(len(bins)):
	outfilew.write(outfmt%(bins[i],dihedral[i]))
	print soutfmt%bins[i],soutfmt%dihedral[i]
outfilew.close()



