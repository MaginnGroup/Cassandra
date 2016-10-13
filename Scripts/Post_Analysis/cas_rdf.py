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

def obtain_atom_indices(output_ndx):
#!Function will determine species and atom index given a cumulative atom type index
#!Used for atom-atom RDF
	atom_ndx = 0
	for i in range(nspecies):
		for j in range(molecule_info[i].natoms):
			if output_ndx <= atom_ndx:
				return(i,j)
			atom_ndx+=1


#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
parser = argparse.ArgumentParser('''*** CASSANDRA Radial Distribution Function Analysis Tool ***

Example Usage:

	python cas_rdf.py -f test.xyz -m s1.mcf s2.mcf (-bins 100) (-o rdf.xvg)

''',formatter_class=RawTextHelpFormatter)
parser.add_argument('-f',action='store',help ='Trajectory File (xyz format)')
parser.add_argument('-m',nargs = '+',action = 'store', help='MCF file(s) for each species. Must be in correct order.')
parser.add_argument('-H',action='store',help='Optional: H Matrix File (specify only if file name without extension is different from xyz file name)',default = 'default')
parser.add_argument('-b',action='store',help='Optional: Beginning frame (default = 1)', default = 1, type=int)
parser.add_argument('-e',action='store',help='Optional: Ending frame (default=last frame)', default=0,type=int)
parser.add_argument('-o',action='store',help='Optional: Output file name',default='rdf.xvg')
parser.add_argument('-bins',action='store',help='Optional: Number of bins',default=100,type=int)
parser.add_argument('-com',action='store_true',help='Optional (Boolean): Output COM-COM rdf profile (default is atom-atom)', default=False)

#DEFINE PARSER ARGS
args = parser.parse_args()
nspecies = len(args.m)
xyzfile = args.f
com_flag = args.com
nbins = args.bins

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
species1_ndx = 0
species2_ndx = 0
atype1_ndx = 0
atype2_ndx = 0
atype1_name = "temp"
atype2_name = "temp"
rdf = np.zeros(nbins,np.float)

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
Lx_array = H[0]
Ly_array = H[1]
Lz_array = H[2]

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
Lmin = np.min([np.min(Lx_array[begin_frame:end_frame]),np.min(Ly_array[begin_frame:end_frame]),
		np.min(Lz_array[begin_frame:end_frame])])

binwidth = Lmin/2.0/nbins

#Construct matrix that contains number of each species for each frame
nmolecules = np.reshape(H[3],(nframes,nspecies))


# (OPTIONAL) READ INPUT FILE
#out = readfiles.read_inp_file(inpfile)
#molecule = readmcf.read_mcf_file(mcffile)


#Store number of atoms of each species
for i in range(nspecies):
	natoms.append(molecule_info[i].natoms)


#Check MCF Files are consistent with first frame



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

#Now we will fill in the indices with the mass. If no molecule is present,
#the mass should equal 0.

#Initialize mass for each atom in the system
atype_mass = np.zeros((max_natoms,nspecies))

#Loop through species and fill masses of each atom type
for i in range(nspecies):
	for j in range(max_natoms):
		try:
			atype_mass[j][i] = float(molecule_info[i].mass[j]) #convert to kg
		except:
			pass




#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#OBTAIN THE INDICES OF SPECIES OR ATOM PAIR OUTPUT

#For COM RDF
if com_flag == True:
	print """
Center-of-mass RDF specified.
Provide index of species 1:
	 """
	for i in range(nspecies):
		print "("+str(i+1)+"): ", molecule_info[i].name 
	print "\n"
	species1_ndx = int(raw_input())-1
	
	print """
Provide index of species 2:
 """
	for i in range(nspecies):
		print "("+str(i+1)+"): ", molecule_info[i].name 
	print "\n"
	species2_ndx = int(raw_input())-1
	
	outname1 = molecule_info[species1_ndx].name	
	outname2 = molecule_info[species2_ndx].name
	print '\n'+"Obtaining center-of-mass RDF for",outname1,"and", outname2
	print "\n"

#FOR ATOM RDF
else:
	print """
Atom-Atom RDF specified.
Provide index of atomtype 1:
"""
	atom_ndx = 0
	for i in range(nspecies):
		for j in range(molecule_info[i].natoms):
			print "("+str(atom_ndx+1)+"):", molecule_info[i].atom_type[j]
			atom_ndx +=1
	print "\n"
	output1_ndx = int(raw_input())-1
	
	#Determine species1 index and atom1 index from user input index
	species1_ndx, atype1_ndx = obtain_atom_indices(output1_ndx)

	print """
Provide index of atomtype 2:
 """
	atom_ndx = 0	
	for i in range(nspecies):
		for j in range(molecule_info[i].natoms):
			print "("+str(atom_ndx+1)+"):", molecule_info[i].atom_type[j]
			atom_ndx +=1
	print "\n"
	output2_ndx = int(raw_input())-1
	
	#Determine species2 index and atom2 index from user input index
	species2_ndx, atype2_ndx = obtain_atom_indices(output2_ndx)
	

	atype1_name = molecule_info[species1_ndx].atom_type[atype1_ndx]
	atype2_name = molecule_info[species2_ndx].atom_type[atype2_ndx]
	outname1 = atype1_name
	outname2 = atype2_name
	print '\n'+"Obtaining atom-atom RDF for",outname1,"and", outname2
	print "\n"



#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
rdf = cas_palib.cas_rdf(xyzfile,
						nbins,
						binwidth,
						nspecies,
						nframes,
						begin_frame,
						end_frame,
						Lx_array,Ly_array, Lz_array, Lmin,
						species1_ndx, species2_ndx,
						atype1_ndx, atype2_ndx, 
						atype1_name, atype2_name,
						atype_mass, max_natoms,max_nmols,
						natoms, nmolecules,
						com_flag,
						rdf)



distance = np.arange(0,binwidth*nbins,binwidth)

#print len(distance), len(rdf)

outfmt = '%12.4f%12.4f\n'
soutfmt = '%12.4f'
if (outfile == 'rdf.xvg'):
	outfile = outfile[:-4]+'_'+outname1+'_'+outname2+'.xvg'

outfilew = open(outfile,'w')

for i in range(len(distance)):
	outfilew.write(outfmt%(distance[i],rdf[i]))
	print soutfmt%distance[i],soutfmt%rdf[i]
outfilew.close()
#PLOT RESULTS INTO FIGURE
#plt.figure(1)
#plt.plot(distance,rdf)
#plt.savefig(outfile[:-4]+'.png')




