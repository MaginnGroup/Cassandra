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
# SCRIPT: mcfgen.py
# VERSION: 1.1
# NEW FEATURES: read GROMACS forcefield files (.itp or .top)
#
# VERSION: 1.0
# ORIGINAL FEATURES: generate a molecular connectivity file (.mcf) from a 
#   configuration file (.pdb or .cml) and a custom forcefield file (.ff)
#*******************************************************************************

#*******************************************************************************
# IMPORT MODULES
#*******************************************************************************
import sys, os, argparse, linecache, re

#*******************************************************************************
#ARGUMENT PARSE
#*******************************************************************************
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=
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
			 
NOTE: The mass of pseudoatoms must be manually entered into the .mcf file.
""")
parser.add_argument('configFile', 
                help="""CONFIGFILE must be in either .pdb or .cml 
format. A .pdb file can be generated using Gaussview, 
while .cml files can be generated using Avogadro. 
""")
parser.add_argument('--ffTemplate', action='store_true',
                help="""Generate a blank force field file template.
""")
parser.add_argument('--ffFile', '-f', nargs=1, 
                help="""The default FFFILE is molecule.ff, a custom format for 
this script. Alternatively, the forcefile parms can 
be supplied in GROMACS format (.itp).
""")
parser.add_argument('--mcfFile', '-m', nargs=1, 
                help="""The default MCFFILE is molecule.mcf.
""")


args = parser.parse_args()


configFile=args.configFile
basename = os.path.splitext(os.path.basename(configFile))[0]
ffTemplate = args.ffTemplate

if args.ffFile:
	ffFile = args.ffFile[0]
	ffFileExt = os.path.splitext(ffFile)[1]
	if ffFileExt == '.ff':
		ffFileType = 'native'
	elif ffFileExt == '.itp' or ffFileExt == '.top':
		ffFileType = 'gromacs'
else:
	ffFile = basename + '.ff'
	ffFileType = 'native'

if args.mcfFile:
	mcfFile = args.mcfFile[0]
else:
	mcfFile = basename + '.mcf'

#*******************************************************************************
# VARIABLE DEFINITIONS
#*******************************************************************************
Rg = 0.008314 #universal gas constant, kJ/mol
periodicTable={'H': 1.0079,
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
'C4': 16.0423}

#*******************************************************************************
# FUNCTION DEFINITIONS
#*******************************************************************************

def check_type_infilename(infilename):
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


	filePdb = open(infilename+'.pdb','w')

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
		filePdb.write('CONECT     '+line[0]+'    ')
		for element in line[1]:
			filePdb.write(element+"   ")
		filePdb.write("\n")

	filePdb.close()

def lookup(atomlook):
	tempfile = open ('temporary.temp', 'r')
	for line in tempfile:
		linestring = line.split()
		if linestring[0]==atomlook:
			line=line.split()
			return line
	tempfile.close()

def initialize(infilename):

	#This function will set the starting atom for the subsequent ring scan.
	#If an a
	#It can handle the following special cases:
	#1) Argon: If no CONECT keyword is found, it will automatically assume it's argon.
        #2) Cyclic UA model: In this case, a 'ghost' molecule will be used as a starting point.


	global cyclic_ua_atom
	conect_found = False
	atom =''
	ifile = open(infilename, 'r')
	tempfile = open('temporary.temp', 'w')
	for line in ifile:
		linestring = line.split()
		if not line.strip():
			continue
		if linestring[0]=='CONECT':
			conect_found = True
			tempfile.write(line[6:len(line)])
			numconnect=len(linestring)-2
			if numconnect==1:
				atom=linestring[1]
				cyclic_ua_atom = False
	ifile.close()

	if conect_found == False: #This means that this is Argon
		atom = ''
		cyclic_ua_atom = False
		tempfile.close()
		return atom

	if cyclic_ua_atom == True and atom == '' : 
	#This means molecule has at least one cyclic united atom model
		tempfile.close()
		atom = '0'

		#Rewrite tempfile to include a dummy atom 0

		tempfile2=open('temporary.temp2','w')
		tempfile1=open('temporary.temp','r')
		for line in tempfile1:
			this_line_new = line.split()
			if this_line_new[0] == '1':
				this_line_new.append('0')
			for each_atom in this_line_new:
				tempfile2.write(each_atom+'  ')
			tempfile2.write('\n')
		tempfile2.write('0 1')

		tempfile2.close()
		tempfile1.close()
		os.system('mv temporary.temp2 temporary.temp')

	return atom


def isRing(scannedList, atom):
	for i in xrange(0,len(scannedList)):
		if scannedList[i] == atom:
			return True


def cleanGenList(branchList):
	for i in xrange(0,len(genScannedList)):
		if genScannedList[i] == branchList[0]:
			for j in xrange(0,len(branchList)):
				del(genScannedList[-1])
			return


def isolateRing(theList):
	location=[]
	seen=set()
	seen_add=seen.add
	seenTwice=set(x for x in theList if x in seen or seen_add(x))
	repeatedAtom=list(seenTwice)

	for index, item in enumerate(theList):
		if item == repeatedAtom[0]:
			location.append(index)
	ringList.append(theList[location[0]:location[1]])
	#logfile.write("The ring list is: " + str(ringList) + "\n")
	cleanGenList(ringList[-1])
	
		
def isAlreadyInRing(atom):
	for row in ringList:
		for i in range(0,len(row)):
			if row[i]==atom:
				return True
	
def bondID():

	#This function will populate the variable "bondList"

	tempfile = open('temporary.temp','r')
	oldatom=[]
	for line in tempfile:
		atomlist1=line.split()
		newAtom=atomlist1[0]
		for i in range(1, len(atomlist1)):
			if atomlist1[i] in oldatom:
				continue
			bondList.append(atomlist1[0] + " " + atomlist1[i])
			oldatom.append(newAtom)
	tempfile.close()


def scan(newAtom,oldLine,scanning):

        #This function does a recursive scan on the molecule. The main product of this function
	#is the list ringList, which contains atoms that form a ring.

	scannedList=[]
	#branchalreadyscanned=[0]
	while 1==1:


		#logfile.write("Atom " + newAtom + "\n")

		if isAlreadyInRing(newAtom)==True:
			#logfile.write("Atom " + newAtom + " was already included in a ring\n\n")
			break

		if isRing(genScannedList, newAtom)==True:
			#logfile.write("A ring was found because atom " + newAtom + 
			#              " was already scanned\n\n")
			genScannedList.append(newAtom)
			return "EndRing", scannedList

		scannedList.append(newAtom)
		genScannedList.append(newAtom)
		newLine=lookup(newAtom)
		totalBondedAtoms=len(newLine)-1

		if totalBondedAtoms>2:  #branch point?

			#logfile.write("A branch point was found in atom " + newAtom + "\n")
			for j in xrange(1,totalBondedAtoms+1):
				if newLine[j]==oldLine[0]: 
				#or newLine[j]==branchalreadyscanned[0]: 
				#If old atom in list is scanned
					continue
				newAtom=newLine[j]
				endType, branchList = scan(newAtom, newLine, True)
				#branchalreadyscanned=lookup(newAtom)
				if endType=="EndPoint":
					cleanGenList(branchList)
				if endType=="EndRing":
					isolateRing(genScannedList)
				#if endType=="EndScan":
			break
	
		for i in xrange(1,totalBondedAtoms+1):

			if scanning==True and totalBondedAtoms==1:
				#logfile.write("End point found at " + newAtom + "\n")
				return "EndPoint", scannedList

			if scanning==False and totalBondedAtoms==1: #Initial step
				scanning = True

			if len(oldLine)==0:
				newAtom=newLine[i]
				oldLine=newLine
				continue

			if newLine[i]==oldLine[0]: #If old atom in list is scanned
				continue
	
			newAtom=newLine[i]
			oldLine=newLine
			break

	return "EndScan", scannedList


def fragID():
	adjacentatoms=[]
	tempfile=open('temporary.temp','r')
	for eachring in ringList:
		for eachatom in eachring:
			vector=lookup(eachatom)
			if len(vector)>3:
				adjacentatoms.append(list(set(vector)-set(eachring)))

		adjacentatoms1=filter(None,adjacentatoms)
		adjacentatoms = [x for sublist in adjacentatoms1 for x in sublist]
		fragList.append(eachring+adjacentatoms)
		adjacentatoms=[]

	for line in tempfile:
		inaring=False
		linestring=line.split()
		if len(linestring)>2:
			anchoratom=linestring[0]
		else:
			continue

		for row in ringList:
			if anchoratom in set(row):
				inaring=True

		if inaring==True:
			continue
		else:
			fragList.append(linestring)

	tempfile.close()

def angleID():

	#This function will populate the list angleList

	tempfile=open('temporary.temp','r')
	for line in tempfile:
		atomlist=line.split()
		if len(atomlist)<3:
			continue
		apex=atomlist[0]
		del(atomlist[0])
		for i in range(0,len(atomlist)-1):
			for j in range(i+1,len(atomlist)):
				angleList.append(atomlist[i] + " " + apex + " " + atomlist[j])

	removedoublecounting("angles")



	
def dihedralID():

	#This function populates the variable dihedralList

	tempfile=open('temporary.temp','r')

	for line in tempfile:
		atomlist1=line.split()
		atom1=atomlist1[0]
		for i in range(1,len(atomlist1)):
			atomlist2=lookup(atomlist1[i])
			atom2=atomlist2[0]
			if atom2==atom1:
				continue
			for j in range(1,len(atomlist2)):
				atomlist3=lookup(atomlist2[j])
				atom3=atomlist3[0]
				if atom3==atom2 or atom3==atom1:
					continue
				for k in range(1,len(atomlist3)):
					atomlist4=lookup(atomlist3[k])
					atom4=atomlist4[0]
					if atom4==atom3 or atom4==atom2 or atom4==atom1:
						continue
					dihedralList.append(atom1 + " " +atom2 + " " +atom3 + " " +atom4)

	removedoublecounting("dihedrals")
	
				
def atomInfo(infilename, outfilename): 
	
	ifile = open(infilename, 'r')
	ofile = open(outfilename,'w')
	ofile.write("!**********************************************************************************\n")
	ofile.write("!Molecular connectivity file for " + infilename + "\n")
	ofile.write("!**********************************************************************************\n")

	ofile.write("\n# Atom_Info\n")

	matrix = forcefield("nonbonded","read")
	matrixcharge = forcefield("charges","read")

	ofile.write(str(len(listofnames))+"\n")
	for line in ifile:
		if not line.strip():
			continue
		if line[0:6]=='HETATM' or line[0:4]=='ATOM':
			atomnumber = line[6:11].strip()
			atomname = line[12:16].strip()
			element = line[76:78].strip()
			if atomname == '':
				atomname = element+ atomnumber
			atomtype = line[80:].strip()
			atommass = str(periodicTable[element])
			vdwtype = "LJ"

			for rowmatrix in matrix:
				for item in rowmatrix[0]:
					if int(item)==int(atomnumber)-1:
						sigma = rowmatrix[2]
						epsilon = rowmatrix[3]
						atom_type_charge = rowmatrix[4]
						for index,line_charge in enumerate(matrixcharge):
							if atomnumber == line_charge[0]:
								charge = line_charge[1] 

			in_a_ring = False
			for eachring in ringList:
				if atomNumber in eachring:
					in_a_ring = True
			if cyclic_ua_atom == True and len(ringList)==0: #This means if there is only one cyclic united atom molecule
				in_a_ring=True

			ofile.write(atomnumber + spacing + atomname + spacing + atomelement + spacing + atommass + spacing + charge + spacing + vdwtype + spacing + epsilon + spacing + sigma)
			if in_a_ring:
				ofile.write(spacing + "ring")
			ofile.write("\n")
	
	ifile.close()
	ofile.close()
	#carfile.close()

def bondInfo(outfilename):
	ofile = open(outfilename,'a')
	ofile.write("\n# Bond_Info\n")
	ofile.write(str(len(bondList))+"\n")
	
	matrix = forcefield("bonded","read")
	

	for index, row in enumerate(bondList):
		bond=row.split()
		for rowmatrix in matrix:
			for item in rowmatrix[0]:
				if index == item:
					length=rowmatrix[2]
					constant = rowmatrix[3]
		ofile.write(str(index+1) + spacing + bond[0] + spacing + bond[1] + spacing + "fixed" + spacing + length+"\n")

	ofile.close()


def angleInfo(outfilename):

	ofile = open(outfilename,'a')

	ofile.write("\n# Angle_Info\n")
	ofile.write(str(len(angleList))+"\n")

	matrix = forcefield("angles","read") #matrix=[repeatedlines, angle type, force constant, angle)

	for index, row in enumerate(angleList):
		angle=row.split()
		for rowmatrix in matrix:
			for item in rowmatrix[0]:
				if index == item:
					angleconstant=rowmatrix[3]
					if "fixed" in  angleconstant:
						angletype = "fixed"
						angleconstant = ""
					else:
						angletype = "harmonic"
					eqmangle=rowmatrix[2]

		ofile.write(str(index+1) + spacing + angle[0] + spacing + angle[1] + spacing + angle[2] + spacing + angletype + spacing + angleconstant + spacing + eqmangle + "\n")

	ofile.close()



def dihedralInfo(outfilename):
	ofile = open(outfilename,'a')

	ofile.write("\n# Dihedral_Info\n")
	ofile.write(str(len(dihedralList))+"\n")

	matrix=forcefield("dihedrals","read")   #matrix=[repeatedlines, dihedral type, K, n, gamma)

	for index, row in enumerate(dihedralList):
		dihedral=row.split()
		for rowmatrix in matrix:
			for item in rowmatrix[0]:
				if index == item:
					if dihedralType == "CHARMM": 
						K = rowmatrix[2]
						n = rowmatrix[3]
						gamma = rowmatrix[4]
						ofile.write(str(index+1) + spacing + dihedral[0] + spacing + dihedral[1] + spacing + dihedral[2] + spacing + dihedral[3] + spacing + dihedralType + spacing + K + spacing + n + spacing + gamma + "\n")
					elif dihedralType == "OPLS":
						a0 = rowmatrix[2]
						a1 = rowmatrix[3]
						a2 = rowmatrix[4]
						a3 = rowmatrix[5]
						ofile.write(str(index+1) + spacing + dihedral[0] + spacing + dihedral[1] + spacing + dihedral[2] + spacing + dihedral[3] + spacing + dihedralType + spacing + a0 + spacing + a1 + spacing + a2 + spacing + a3 + "\n")
					elif dihedralType == "harmonic":
						k = rowmatrix[2]
						phi = rowmatrix[3]
						ofile.write(str(index+1) + spacing + dihedral[0] + spacing + dihedral[1] + spacing + dihedral[2] + spacing + dihedral[3] + spacing + dihedralType + spacing + k + spacing + phi+ "\n")



				

	ofile.close()
		

def fragInfo(mcfFile):
	"""arguments:
	mcfFile, string, the .mcf file to write to
returns:
	none
"""
	mcf = open(mcfFile,'a')
	mcf.write("\n!Fragment Format\n")
	mcf.write('!index number_of_atoms_in_fragment branch_point other_atoms\n')
	mcf.write("\n# Fragment_Info\n")

	global fragList

	if cyclic_ua_atom == True and len(ringList)==0: 
	#This means if there is only one cyclic united atom molecule
		fragList=[]
		fragList.append([str(i) for i in atomList])

	if len(fragList)==0:
		if len(atomList)==1:
			mcf.write('1\n')
			mcf.write('1 1 1\n')
		elif len(atomList)==2:
			mcf.write('1\n')
			mcf.write('1 2 1 2\n')

		# the below elif statement allows for zeolites
		elif len(atomList)>2:
			mcf.write('1\n')
			atomList_indices = []
			counter = 0
			# count the atoms
			for item in atomList:
				counter = counter + 1
				atomList_indices.append(counter)
			str_to_add = " "
			str_type_list = [str(item) for item in atomList_indices]
			new_str_list = str_to_add.join(str_type_list)
			mcf.write("1 %s %s" %(str(atomList_indices[-1]),new_str_list))
	else:
		mcf.write(str(len(fragList))+"\n")
		for index, row in enumerate(fragList):
			mcf.write(str(index+1) + spacing + str(len(row)))
			for item in row:
				mcf.write(spacing+item)
			mcf.write("\n")
	mcf.close()


def tablelookup(element, thingtolook):

	for i in xrange(0,108):
		line = pt[i]
		if line[3]==element:
			if thingtolook=="mass":
				return line[1]
			elif thingtolook=="name":
				return line[2]
			elif thingtolook=="atomic_number":
				return line[0]


def removedoublecounting(thingtoclean):

	if thingtoclean=="dihedrals":
		linesrepeated=[]
		cleanlist=[]
		for index,row in enumerate(dihedralList):
			set1=set(row.split())
			for i in range(index+1, len(dihedralList)):
				set2=set(dihedralList[i].split())
				if set1==set2:
					linesrepeated.append(i)
	
		for index,row in enumerate(dihedralList):
			if index in linesrepeated:
				continue
			cleanlist.append(row)

		del(dihedralList[0:len(dihedralList)])	
		for i in range(0,len(cleanlist)):
			dihedralList.append(cleanlist[i])

	if thingtoclean=="angles":
		linesrepeated=[]
		cleanlist=[]
		for index,row in enumerate(angleList):
			set1=set(row.split())
			for i in range(index+1, len(angleList)):
				set2=set(angleList[i].split())
				if set1==set2:
					linesrepeated.append(i)
	
		for index,row in enumerate(angleList):
			if index in linesrepeated:
				continue
			cleanlist.append(row)

		del(angleList[0:len(angleList)])	
		for i in range(0,len(cleanlist)):
			angleList.append(cleanlist[i])

	if thingtoclean=="fragConn":
		linesrepeated=[]
		cleanlist=[]
		for index,row in enumerate(fragConn):
			set1=set(row.split())
			for i in range(index+1, len(fragConn)):
				set2=set(fragConn[i].split())
				if set1==set2:
					linesrepeated.append(i)
	
		for index,row in enumerate(fragConn):
			if index in linesrepeated:
				continue
			cleanlist.append(row)
		
		del(fragConn[0:len(fragConn)])	
		for i in range(0,len(cleanlist)):
			fragConn.append(cleanlist[i])

def fragConnectivity():

	#This function will populate the list fragConn
	group=[]
	j=0

	if cyclic_ua_atom == True and len(ringList)==0:
		fragConn.append("0")
		return

	for index1, row in enumerate(fragList):
		for j in range(0,len(row)-1):
			for i in range(j+1,len(row)):
				element1=row[j]
				element2=row[i]
				group.append(element1)
				group.append(element2)
				for index2, otherrow in enumerate(fragList):
					#if otherrow in oldrows:
					#	continue
					a=set(group)
					b=set(otherrow)
					if a<=b and index1!=index2:
						fragConn.append(str(index1+1) + spacing + str(index2+1))
				group=[]

	removedoublecounting("fragConn")
			
			
def fragConnectivityInfo(mcfFile):
	mcf=open(mcfFile,'a')
	mcf.write('\n\n# Fragment_Connectivity\n')
	mcf.write(str(len(fragConn))+'\n')
	for index, row in enumerate(fragConn):
		mcf.write(str(index+1) + spacing + row)
		mcf.write('\n')
	mcf.write('\n\nEND\n')
	mcf.close()

def improperInfo(mcfFile):
	mcf=open(mcfFile,'a')
	mcf.write("\n# Improper_Info\n")
	mcf.write("0\n")
	mcf.close()

def ffFileGeneration(infilename,outfilename):

	print "\n\n*********Force field template file generation*********\n"

	print "Reading Modified PDB File..."
	ifile = open(infilename,'r')
	atomtypenames=[]
	for line in ifile:
		if not line.strip():
			continue
		linestring = line.split()
		if linestring[0]=='HETATM' or linestring[0] == 'ATOM':
			atomNumber = linestring[1]
			name = linestring[-1]
			listofnames.append(atomNumber + " " + name)
			atomtypenames.append(name)
	ifile.close()
	atomtypenames = list(set(atomtypenames))
	numAtomTypes = str(len(atomtypenames))

	global dihedralType
	if len(dihedralList) > 0:
		dihedralType = raw_input("Enter the dihedral functional form " + 
		                         "(CHARMM/OPLS/harmonic): ")
	else:
		dihedralType = "NONE" 
		#This is just a default. It will not affect molecules with no dihedrals.

	#Identify different types of angles and dihedrals based on atom types
	#Each row in dihedarlList_atomType has each dihedral expressed as atomType
	#Each row in angleList_atomType has each angle expressed as atomType


	for i, rowdihedral in enumerate(dihedralList):
		rowdihedralstr=rowdihedral.split()
		dihedralList_atomType.append(lookupinlist(listofnames, rowdihedralstr))

	for i, rowangle in enumerate(angleList):
		rowanglestr=rowangle.split()
		angleList_atomType.append(lookupinlist(listofnames, rowanglestr))

	for i, rowbond in enumerate(bondList):
		rowbondstr = rowbond.split()
		bondList_atomType.append(lookupinlist(listofnames, rowbondstr))

	ff=open(ffFile,'a')
	ff.write("atomtypes\n"+numAtomTypes+ "\n\n")
	ff.write("begin atom-atomtype\n")
	for line in listofnames:
		ff.write(line + "\n")
	ff.write("end atom-atomtype\n\n")

	ff.write("dihedraltype "+dihedralType+"\n\n")
	ff.close()

	forcefield("nonbonded","write")
	forcefield("bonded","write")
	forcefield("angles","write")
	forcefield("dihedrals","write")
	forcefield("charges","write")

def forcefield(degreeoffreedom,action):

	#Answer will be a matrix whose rows, given by linematrix, will be:
	#[list of atoms repeated, type of atom, sigma, epsilon]

	global dihedralType

	answer=[]
	linematrix=[]

	# If the script will read from the ff file, the list that relates atom name
	# and atom type is still unknown. 
	# The following lines read the atom-atomtype section in the file, and save it 
	# in listofnames. 
	# listofnames=[atomNumber atomType]

	if action == "read":
		
		listofnames[:] = []
		lineindex = 0
		flag = 0
		ff = open(ffFile,'r')
		for linea in ff:
			lineindex=lineindex+1
			if "begin atom-atomtype" in linea:
				flag = 1
				continue

			if "end atom-atomtype" in linea:
				flag = 0
				break	

			if flag ==1 :
				listofnames.append(linea[0:-1])


	if degreeoffreedom == "nonbonded":
	# The following lines will read listofnames and identify which atoms numbers 
	# share the same atom type. One variable is produced: Linesrepeated, which 
	# contains those atoms that share the same atom type.
		if action == "read":
			linesrepeated=[]
			rowalreadyscanned=[]

			for index,row in enumerate(listofnames):
				rowvector1=row.split()
				flag=0
				for oldsets in rowalreadyscanned:
					oldset=oldsets.split()
					if oldset[1]==rowvector1[1]:
						flag=1
					
				if flag==1: continue
	
				for i in range(index+1, len(listofnames)):
					rowvector2=listofnames[i].split()
					if rowvector1[1]==rowvector2[1]:
						linesrepeated.append(i)
						
				linesrepeated.append(index)
				rowalreadyscanned.append(row)

				linematrix.append(linesrepeated)

				# Using the forcefield file, we need to build a matrix with the 
				# following form:
				#   MATRIX == #[list of atoms repeated, type of atom, sigma, epsilon, charge]
				# To do that, we need to find where the "nonbonded" labels are are in 
				# the ff file so we can get sigma, epsilon, and charge

				linewherenonbonded=[]
				ff=open(ffFile,'r')
				for index, line in enumerate(ff):
					if "nonbonded" in line:
						linewherenonbonded.append(index+1)
				ff.close()
			
				# Now, the following four lines will contain the elements of the matrix

				for i in linewherenonbonded:
					vartemporal = linecache.getline(ffFile, i+1).split()[0]
					if vartemporal==rowvector1[1]:
						linematrix.append(linecache.getline(ffFile, i+1)[0:-1])
						linematrix.append(linecache.getline(ffFile, i+2)[6:-1])
						linematrix.append(linecache.getline(ffFile, i+3)[8:-1])
						linematrix.append(linecache.getline(ffFile, i+4)[17:-1])
				answer.append(linematrix)
				linematrix=[]		
				linesrepeated=[]

		elif action == "write":

			linesrepeated=[]
			rowalreadyscanned=[]

			for index,row in enumerate(listofnames):
				rowvector1=row.split()
				flag=0
				for oldsets in rowalreadyscanned:
					oldset=oldsets.split()
					if oldset[1]==rowvector1[1]:
						flag=1
					
				if flag==1: continue
	
				for i in range(index+1, len(listofnames)):
					rowvector2=listofnames[i].split()
					if rowvector1[1]==rowvector2[1]:
						linesrepeated.append(i)
						
				linesrepeated.append(index)
				rowalreadyscanned.append(row)


	
				ff=open(ffFile,'a')
				ff.write("nonbonded\n")
				
				ff.write(row.split()[1] + "\n")
				ff.write("Sigma \n")
				ff.write("Epsilon \n")
				ff.write("atom_type_charge \n\n")
				linesrepeated=[]
				ff.close()	


	if degreeoffreedom=="bonded":

		if action == "read":

			for i, rowbond in enumerate(bondList):
				rowbondstr = rowbond.split()
				bondList_atomType.append(lookupinlist(listofnames, rowbondstr))

			linesrepeated=[]
			rowalreadyscanned=[]

			for index,row in enumerate(bondList_atomType):
				flag=0
				for oldsets in rowalreadyscanned:
					if oldsets==row:
						flag=1
					else:
						reverse = [j for j in reversed(oldsets)]
						if row==reverse:
							flag=1
				if flag==1: continue

				for i in range(index+1, len(bondList_atomType)):
					if row==bondList_atomType[i]:
						linesrepeated.append(i)
					else:
						reverse = [j for j in reversed(bondList_atomType[i])]
						if row==reverse:
							linesrepeated.append(i)
	
				linesrepeated.append(index)
				rowalreadyscanned.append(row)

				linematrix.append(linesrepeated)

				#We need to find where the "bonds" labels are are in the ff
				#file so we can get K and angle

				linewherebonded=[]
				ff=open(ffFile,'r')
				for index, line in enumerate(ff):
					if "bonds" in line:
						linewherebonded.append(index+1)
				ff.close()
			
				#Now, the following two lines will contain the elements of the matrix

				for i in linewherebonded:

					#Read from ffFile which bondtype we'll work with (e.g. A-A, A-B)
					vartemporal = linecache.getline(ffFile, i+1)[0:-1].split()

					if set(vartemporal)==set(row):
						linematrix.append(linecache.getline(ffFile, i+1)[0:-1])
						linematrix.append(linecache.getline(ffFile, i+2)[6:-1])
						linematrix.append(linecache.getline(ffFile, i+3)[8:-1])

				answer.append(linematrix)
				linematrix=[]		
				linesrepeated=[]


		if action=="write":

			linesrepeated=[]
			rowalreadyscanned=[]

			for index,row in enumerate(bondList_atomType):
				flag=0
				for oldsets in rowalreadyscanned:
					if oldsets==row:
						flag=1
					else:
						reverse = [j for j in reversed(oldsets)]
						if row==reverse:
							flag=1
				if flag==1: continue

				for i in range(index+1, len(bondList_atomType)):
					if row==bondList_atomType[i]:
						linesrepeated.append(i)
					else:
						reverse = [j for j in reversed(bondList_atomType[i])]
						if row==reverse:
							linesrepeated.append(i)
	
				linesrepeated.append(index)
				rowalreadyscanned.append(row)

				ff=open(ffFile,'a')
				ff.write("bonds\n")
				ff.write(str(row[0]) + " " + str(row[1]) + "\n")
				ff.write("Length \n")
				ff.write("Constant \n\n")
				linesrepeated=[]
				ff.close()


	if degreeoffreedom=="angles":

		if action == "read":
			for i, rowangle in enumerate(angleList):
				rowanglestr=rowangle.split()
				angleList_atomType.append(lookupinlist(listofnames, rowanglestr))

			linesrepeated=[]
			rowalreadyscanned=[]

			for index,row in enumerate(angleList_atomType):
				flag=0
				for oldsets in rowalreadyscanned:
					if oldsets==row:
						flag=1
					else:
						reverse = [j for j in reversed(oldsets)]
						if row==reverse:
							flag=1
				if flag==1: continue

				for i in range(index+1, len(angleList_atomType)):
					if row==angleList_atomType[i]:
						linesrepeated.append(i)
					else:
						reverse = [j for j in reversed(angleList_atomType[i])]
						if row==reverse:
							linesrepeated.append(i)
	
				linesrepeated.append(index)
				rowalreadyscanned.append(row)


				linematrix.append(linesrepeated)

				#We need to find where the "angles" labels are are in the ff
				#file so we can get K and angle

				linewhereangle=[]
				ff=open(ffFile,'r')
				for index, line in enumerate(ff):
					if "angles" in line:
						linewhereangle.append(index+1)
				ff.close()
			
				#Now, the following two lines will contain the elements of the matrix

				for i in linewhereangle:

					#Read from ffFile which angle type we'll work with (e.g. A-A-B, A-B-A)
					vartemporal = linecache.getline(ffFile, i+1)[0:-1].split()

					if vartemporal==row:
						linematrix.append(linecache.getline(ffFile, i+1)[0:-1])
						linematrix.append(linecache.getline(ffFile, i+2)[6:-1])
						linematrix.append(linecache.getline(ffFile, i+3)[8:-1])

				answer.append(linematrix)
				linematrix=[]		
				linesrepeated=[]


		if action=="write":

			linesrepeated=[]
			rowalreadyscanned=[]

			for index,row in enumerate(angleList_atomType):
				flag=0
				for oldsets in rowalreadyscanned:
					if oldsets==row:
						flag=1
					else:
						reverse = [j for j in reversed(oldsets)]
						if row==reverse:
							flag=1
				if flag==1: continue

				for i in range(index+1, len(angleList_atomType)):
					if row==angleList_atomType[i]:
						linesrepeated.append(i)
					else:
						reverse = [j for j in reversed(angleList_atomType[i])]
						if row==reverse:
							linesrepeated.append(i)
	
				linesrepeated.append(index)
				rowalreadyscanned.append(row)


				ff=open(ffFile,'a')
				ff.write("angles\n")
				ff.write(row[0] + " " + row[1] + " " +row[2] + "\n")
				ff.write("Angle \n")
				ff.write("Constant \n\n")
				linesrepeated=[]
				ff.close()


	if degreeoffreedom=="dihedrals":

		if action=="read":

			for i, rowdihedral in enumerate(dihedralList):
				rowdihedralstr=rowdihedral.split()
				dihedralList_atomType.append(lookupinlist(listofnames, rowdihedralstr))
		
			linesrepeated=[]
			rowalreadyscanned=[]

			for index,row in enumerate(dihedralList_atomType):
				flag=0
				for oldsets in rowalreadyscanned:
					if oldsets==row:
						flag=1
					else:
						reverse = [j for j in reversed(oldsets)]
						if row==reverse:
							flag=1
				if flag==1: continue

				for i in range(index+1, len(dihedralList_atomType)):
					if row==dihedralList_atomType[i]:
						linesrepeated.append(i)
					else:
						reverse = [j for j in reversed(dihedralList_atomType[i])]
						if row==reverse:
							linesrepeated.append(i)
	
				linesrepeated.append(index)
				rowalreadyscanned.append(row)


				linematrix.append(linesrepeated)

				#We need to find where the "dihedrals" labels are are in the ff
				#file so we can get K and angle

				linewheredihedral=[]
				ff=open(ffFile,'r')
				for index, line in enumerate(ff):
					if "dihedrals" in line:
						linewheredihedral.append(index+1)
					if line[0:12] == "dihedraltype":
						dihedralType = line.split()[1]
				ff.close()
			

				for i in linewheredihedral:

					vartemporal = linecache.getline(ffFile, i+1)[0:-1].split()

					if vartemporal==row:
						if dihedralType == "CHARMM":
							linematrix.append(linecache.getline(ffFile, i+1)[0:-1])
							linematrix.append(linecache.getline(ffFile, i+2)[1:-1])
							linematrix.append(linecache.getline(ffFile, i+3)[1:-1])
							linematrix.append(linecache.getline(ffFile, i+4)[5:-1])
						elif dihedralType == "OPLS":
							linematrix.append(linecache.getline(ffFile, i+1)[0:-1])
							linematrix.append(linecache.getline(ffFile, i+2)[2:-1])
							linematrix.append(linecache.getline(ffFile, i+3)[2:-1])
							linematrix.append(linecache.getline(ffFile, i+4)[2:-1])
							linematrix.append(linecache.getline(ffFile, i+5)[2:-1])
						elif dihedralType == "harmonic":
							linematrix.append(linecache.getline(ffFile, i+1)[0:-1])
							linematrix.append(linecache.getline(ffFile, i+2)[2:-1])
							linematrix.append(linecache.getline(ffFile, i+3)[3:-1])


				answer.append(linematrix)
				linematrix=[]		
				linesrepeated=[]


		if action=="write":


			linesrepeated=[]
			rowalreadyscanned=[]

			for index,row in enumerate(dihedralList_atomType):
				flag=0
				for oldsets in rowalreadyscanned:
					if oldsets==row:
						flag=1
					else:
						reverse = [j for j in reversed(oldsets)]
						if row==reverse:
							flag=1
				if flag==1: continue

				for i in range(index+1, len(dihedralList_atomType)):
					if row==dihedralList_atomType[i]:
						linesrepeated.append(i)
					else:
						reverse = [j for j in reversed(dihedralList_atomType[i])]
						if row==reverse:
							linesrepeated.append(i)
	
				linesrepeated.append(index)
				rowalreadyscanned.append(row)



				ff=open(ffFile,'a')
				ff.write("dihedrals\n")
				ff.write(row[0] + " " + row[1] + " " +row[2] + " " + row[3] + "\n")

#				global dihedralType

				if dihedralType == "CHARMM":
					ff.write("K \n")
					ff.write("n \n")
					ff.write("Gamma \n\n")
					linesrepeated=[]
					ff.close()
				elif dihedralType == "OPLS":
					ff.write("a0 \n")
					ff.write("a1 \n")
					ff.write("a2 \n")
					ff.write("a3 \n\n")
					linesrepeated=[]
					ff.close()
				elif dihedralType == "harmonic":
					ff.write("K \n")
					ff.write("phi \n\n")
					linesrepeated=[]
					ff.close()


	if degreeoffreedom=="charges":

		if action == "write":
		
			for index,row in enumerate(listofnames):
				ff=open(ffFile,'a')
				ff.write("charge\n")
				ff.write(row.split()[0] + " \n\n")
	

			ff.write("end")
			ff.close()

		if action == "read":

			ff=open(ffFile,'r')
			for index, line in enumerate(ff):
				if "charge" in line:
					answer.append(linecache.getline(ffFile, index+2).split())
			ff.close()
		
			


	if action=="read": return answer


def lookupinlist(namelist, degreeoffreedomlist):
	row = []
	for element in degreeoffreedomlist:
		for item in namelist:
			itemstr=item.split()
			if itemstr[0] == element:
				row.append(itemstr[1])
				break
	return row

def readPdb(pdbFile):
	"""arguments:
	pdbFile, string = filename of a pdb file
returns:
	atomList, list of ints = atom numbers from pdb file
	atomParms, dict = atom parameters: name, type, element, mass
"""
	# initialize variables
	atomList = []
	atomParms = {}
	bondList = []
	pdb = open(pdbFile,'r')

	for line in pdb:
		# read atom info
		this_line = line.split()
		if line[0:6]=='HETATM' or line[0:4]=='ATOM': 
			i = int(line[6:11].strip())
			atomList.append(i)
			atomParms[i] = {}
			atomParms[i]['name'] = line[12:16].strip()
			atomParms[i]['element'] = line[76:78].strip().title()
			if atomParms[i]['name'] == '':
				atomParms[i]['name'] = atomParms[i]['element'] + i
			atomParms[i]['type'] = line.split()[-1].strip()
			atomParms[i]['mass'] = periodicTable[atomParms[i]['element']]
		# read bond info
		if "CONECT" in this_line:
			lineList = line.split()
			i = lineList[1]
			for j in lineList[2:]:
				if not any([True for bond in bondList if bond[0]==int(j)]):
					bondList.append((int(i), int(j)))

	pdb.close()
	return atomList, atomParms

class Error(Exception):

	"""Base class for exceptions in this module."""
	pass

def readNative(ffFile, atomParms):
	"""arguments:
	ffFile, string = filename of the forcefield file
	atomParms, dict = i: {name:, type:, element:, mass:}
returns:
	atomParms, dict = i: {name:, type:, element:, mass:, vdw:}
	bondParms, dict = (i,j): (type, length)	
	angleParms, dict = (i,j,k): (type, [ktheta,] angle)
	dihedralParms, dict = (i,j,k,l): (type, parms)
"""
	
	bondParms = {}
	angleParms = {}
	dihedralParms = {}

	ff = open(ffFile,'r')

	line = ff.readline()
	while line:
		if 'dihedraltype' in line:
			dihedralType = line.split()[1]
		elif 'nonbonded' in line:
			index = ff.readline().strip() #atomType
			sigma = float(ff.readline().split()[1])
			epsilon = float(ff.readline().split()[1])
			try:
				atom_type_charge = float(ff.readline().split()[1])
			except:
				atom_type_charge = 'None'

			if index not in atomParms: 
				atomParms[index] = {}
			atomParms[index]['vdw'] = ('LJ', epsilon, sigma)
			atomParms[index]['charge'] = (atom_type_charge)
		elif 'bonds' in line:
			index = tuple(ff.readline().split()) #atomType
			distance = float(ff.readline().split()[1])
			bondType = ff.readline().split()[1]
			if bondType != 'fixed':
				raise Error('Cassandra does not support ' + bondType + 'bonds.')
			bondParms[index] = (bondType, distance)
		elif 'angles' in line:
			index = tuple(ff.readline().split()) #atomType
			theta = float(ff.readline().split()[1])
			ktheta = ff.readline().split()[1]
			if ktheta == 'fixed':
				angleParms[index] = ('fixed', theta)
			else:
				ktheta = float(ktheta)
				angleParms[index] = ('harmonic', ktheta, theta)
		elif 'dihedrals' in line:
			index = tuple(ff.readline().split()) #atomType
			if dihedralType == 'CHARMM':
				a0 = float(ff.readline().split()[1])
				a1 = float(ff.readline().split()[1])
				delta = float(ff.readline().split()[1])
				dihedralParms[index] = (dihedralType, a0, a1, delta)
			elif dihedralType == 'OPLS':
				c0 = float(ff.readline().split()[1])
				c1 = float(ff.readline().split()[1])
				c2 = float(ff.readline().split()[1])
				c3 = float(ff.readline().split()[1])
				dihedralParms[index] = (dihedralType, c0, c1, c2, c3)
			elif dihedralType == 'harmonic':
				kphi = float(ff.readline().split()[1])
				phi = float(ff.readline().split()[1])
				dihedralParms[index] = (dihedralType, kphi, phi)
			elif dihedralType == 'none':
				dihedralParms[index] = (dihedralType,)
		elif 'charge' in line:
			data = ff.readline().split()
			index = int(data[0]) #atomNumber
			if len(data)>1:
				atomParms[index]['charge'] = float(data[1])
				
			# else, if the information will be provided by atom type and corrected by checkParms(), do nothing
				
		line = ff.readline()

	return atomParms, bondParms, angleParms, dihedralParms

def readGromacs(ffFile, atomParms, bondParms, angleParms, dihedralParms, 
                vdwType = None, comboRule = None):
	"""arguments:
	ffFile, string = filename of the forcefield file
	atomParms, dict = i: {name:, type:, element:, mass:}
	bondParms, dict = (i,j): (type, length)	
	angleParms, dict = (i,j,k): (type, [ktheta,] angle)
	dihedralParms, dict = (i,j,k,l): (type, parms)
optional arguments:
	vdwType, string = 'LJ'
	comboRule, string = identifies the LJ parameters provided
returns:
	atomParms, dict = i: {name:, type:, element:, mass:, vdw:}
	bondParms, dict = (i,j): (type, length)	
	angleParms, dict = (i,j,k): (type, [ktheta,] angle)
	dihedralParms, dict = (i,j,k,l): (type, parms)
"""

	ff = open(ffFile,'r')

	line = ff.readline()
	while line:
		if line.startswith('#include'):
			includeFile = os.path.expanduser(line.split()[1].strip('"'))
			if not os.path.isfile(includeFile):
				parentDir = os.path.dirname(ffFile)
				includeFile = os.path.join(parentDir, includeFile)
			if os.path.isfile(includeFile):
				atomParms, bondParms, angleParms, dihedralParms = readGromacs(
				includeFile, atomParms, bondParms, angleParms, dihedralParms, 
				vdwType, comboRule)
			else:
				print 'WARNING: Topology file ' + includeFile + ' not found. ' + \
				      'Continuing without reading file.'
		elif line.startswith('['):
			section = line.strip() #store the section header
			line = ff.readline()
			while line and not line.isspace(): #section ends with a blank line
				if not line.startswith(';'):
					data = line.split()
					if 'default' in section:
						# Look for function type 1, to make sure this is Lennard-Jones
						if data[0] == '1':
							vdwType = 'LJ'
						else:
							raise Error('Nonbonded forcefield is not LJ. Cassandra only ' + 
													'supports LJ interactions at this time.')
						comboRule = data[1]
					# Look for atomParms
					elif 'atoms' in section:
						index = int(data[0]) #atomNumber
						if index not in atomParms:
							atomParms[index] = {}
						atomParms[index]['charge'] = float(data[6])
					elif 'atomtypes' in section:
						if data[4] != 'A':
							raise Error('ptype is not A. Cassandra only supports point ' + 
													'particles.')
						index = data[0].upper() #atomType
						if index not in atomParms:
							atomParms[index] = {}
						atomParms[index]['mass'] = float(data[2])
						atomParms[index]['charge'] = float(data[3])
						if comboRule == '1':
							C6 = float(data[5])
							C12 = float(data[6])
							sigma = ((C12 / C6)**(1/6.)) * 10
							epsilon = C6**2 / 4 / C12 / Rg
						elif comboRule == '2' or comboRule == '3':
							sigma = float(data[5]) * 10
							epsilon = float(data[6]) / Rg
						atomParms[index]['vdw'] = (vdwType, epsilon, sigma)
					# Look for bondParms
					elif 'bond' in section:
						if 'type' in section:
							index = (data[0].upper(), data[1].upper()) #atomTypes
						else:
							index = (int(data[0]), int(data[1])) #atomNumbers
						if not any([data[2]=='1', data[2]=='2', data[2]=='3', 
							          data[2]=='4']):
							raise Error('Cassandra only supports fixed bonds at this time.')
						b0 = float(data[3]) * 10
						bondParms[index] = ('fixed', b0)
					# Look for angleParms
					elif 'angle' in section:
						if 'type' in section:
							index = (data[0].upper(), data[1].upper(), data[2].upper())
						else:
							index = (int(data[0]), int(data[1]), int(data[2]))
						if data[3]!='1':
							raise Error('Cassandra supports fixed or harmonic angles.')
						theta = float(data[4])
						ktheta = float(data[5]) / 2. / Rg
						angleParms[index] = ('harmonic', ktheta, theta)
					# Look for dihedralParms
					elif 'dihedral' in section:
						if 'type' in section:
							index = (data[0].upper(), data[1].upper(), data[2].upper(), 
							         data[3].upper())
						else:
							index = (int(data[0]), int(data[1]), int(data[2]), int(data[3]))
						if data[4]=='1' or data[4]=='4' or data[4]=='9': #CHARMM
							delta = float(data[5])
							a0 = float(data[6])
							a1 = float(data[7])
							dihedralParms[index] = ('CHARMM', a0, a1, delta)
						elif data[4] == '2': #harmonic
							phi = float(data[5])
							kphi = float(data[6]) / 2. / Rg
							dihedralParms[index] = ('harmonic', kphi, phi)
						elif data[4] == '3': #Ryckaert-Bellemans
							c0 = float(data[5])
							c1 = float(data[6])
							c2 = float(data[7])
							c3 = float(data[8])
							c4 = float(data[9])
							a0 = c0 + c2 / 2. + c4/ 2.
							a1 = c1 + 3. * c3 / 4.
							a2 = c2 / 2. + c4 / 2.
							a3 = c3 / 4.
							if not c4 == 0. and c5 == 0.:
								raise Error('Can only convert Ryckaert-Bellemans dihedrals ' + 
								            'to OPLS if c4==0 and c5==0.\n')
							dihedralParms[index] = ('OPLS', a0, a1, a2, a3)
						elif data[4] == '5': #OPLS
							a0 = 0.
							a1 = float(data[5]) / 2.
							a2 = float(data[6]) / 2.
							a3 = float(data[7]) / 2.
							dihedralParms[index] = ('OPLS', a0, a1, a2, a3)
						else:
							raise Error('Cassandra supports OPLS, CHARMM and harmonic ' + 
													'dihedrals.')
				line = ff.readline()
		line = ff.readline()

	ff.close()

	return atomParms, bondParms, angleParms, dihedralParms

def checkParms(atomList, bondList, angleList, dihedralList, 
               atomParms, bondParms, angleParms, dihedralParms):
	"""arguments:
	atomList, list =
	bondList, list =
	angleList, list =
	dihedralList, list =
	atomParms, dict =
	bondParms, dict =
	angleParms, dict =
	dihedralParms, dict =
returns:
	atomParms, dict =
	bondParms, dict =
	angleParms, dict =
	dihedralParms, dict =
"""

	# Check that we have all the parms we need
	for i in atomList:
		for parm in ['vdw', 'charge']:
			if parm not in atomParms[i]:
				iType = atomParms[i]['type']
				if parm in atomParms[iType]:
					atomParms[i][parm] = atomParms[iType][parm]
				else:
					raise Error(parm + ' parms for atom ' + str(i) + ' not found.')
	for bond in bondList:
		ij = tuple([int(i) for i in bond.split()])
		ji = ij[::-1]
		if ij not in bondParms and ji not in bondParms:
			ijType = tuple([atomParms[i]['type'] for i in ij])
			jiType = ijType[::-1]
			if ijType in bondParms:
				bondParms[ij] = bondParms[ijType]
			elif jiType in bondParms:
				bondParms[ij] = bondParms[jiType]
			else:
				raise Error('Bond parms for atoms ' + bond + ' not found.')
	for angle in angleList:
		ijk = tuple([int(i) for i in angle.split()])
		kji = ijk[::-1]
		if ijk not in angleParms and kji not in angleParms:
			ijkType = tuple([atomParms[i]['type'] for i in ijk])
			kjiType = ijkType[::-1]
			if ijkType in angleParms:
				angleParms[ijk] = angleParms[ijkType]
			elif kjiType in angleParms:
				angleParms[ijk] = angleParms[kjiType]
			else:
				raise Error('Angle parms for atoms ' + angle + ' not found.')
	for dihedral in dihedralList:
		ijkl = tuple([int(i) for i in dihedral.split()])
		lkji = ijkl[::-1]
		if ijkl not in dihedralParms and lkji not in dihedralParms:
			ijklType = tuple([atomParms[i]['type']for i in ijkl])
			lkjiType = ijklType[::-1]
			if ijklType in dihedralParms:
				dihedralParms[ijkl] = dihedralParms[ijklType]
			elif lkjiType in dihedralParms:
				dihedralParms[ijkl] = dihedralParms[lkjiType]
			else:
				raise Error('Dihedral parms for atoms ' + dihedral + ' not found.')	
			
	return atomParms, bondParms, angleParms, dihedralParms
	print "Revised atomParms."
	print atomParms

def writeMcf(configFile, mcfFile, 
             atomList, bondList, angleList, dihedralList, ringList,
             atomParms, bondParms, angleParms, dihedralParms):
	"""arguments:
	configFile, string = configuration file
	mcfFile, string = molecular connectivity file
	atomList, list = 
	bondList, list = 
	angleList, list = 
	dihedralList, list = 
	ringList, list = 
	atomParms, dict = 
	bondParms, dict = 
	angleParms, dict = 
	dihedralParms, dict =
returns:
	none
"""

	mcf = open(mcfFile, 'w')	

	mcf.write('!***************************************' + 
	          '****************************************\n')
	mcf.write('!Molecular connectivity file for ' + configFile + '\n')
	mcf.write('!***************************************' + 
	          '****************************************\n')

	mcf.write('!Atom Format\n')
	mcf.write('!index type element mass charge type parameters\n' + 
	          '!type="LJ", parms=epsilon sigma\n')
	mcf.write('\n# Atom_Info\n')
	mcf.write(str(len(atomList))+'\n')
	for i in atomList:
		mcf.write('%-4d'     % (i))
		mcf.write('  %-6s'   % (atomParms[i]['name']))
		mcf.write('  %-2s'   % (atomParms[i]['element']))
		mcf.write('  %7.3f'  % (atomParms[i]['mass']))
		mcf.write('  %6.3f' % (atomParms[i]['charge']))
		mcf.write('  %2s  %8.3f  %8.3f' % atomParms[i]['vdw'])
		for ring in ringList:
			if str(i) in ring: mcf.write('  ring')
		mcf.write('\n')
	
	if bondList:
		mcf.write('\n!Bond Format\n')
		mcf.write('!index i j type parameters\n' + 
							'!type="fixed", parms=bondLength\n')
	mcf.write('\n# Bond_Info\n')
	mcf.write(str(len(bondList)) + '\n')
	nBond = 0
	for bond in bondList:
		ij = tuple([int(i) for i in bond.split()])
		if ij not in bondParms: ij = ij[::-1]
		nBond = nBond + 1
		mcf.write('%-4d' % (nBond))
		mcf.write('  %-4d  %-4d' % ij)
		mcf.write('  %5s  %8.3f' % bondParms[ij])
		mcf.write('\n')

	if angleList:
		mcf.write('\n!Angle Format\n')
		mcf.write('!index i j k type parameters\n' +
							'!type="fixed", parms=equilibrium_angle\n' + 
							'!type="harmonic", parms=force_constant equilibrium_angle\n')
	mcf.write('\n# Angle_Info\n')
	mcf.write(str(len(angleList)) + '\n')
	nAngle = 0
	for angle in angleList:
		ijk = tuple([int(i) for i in angle.split()])
		if ijk not in angleParms: ijk = ijk[::-1]
		nAngle = nAngle + 1
		mcf.write('%-4d' % (nAngle))
		mcf.write('  %-4d  %-4d  %-4d' % ijk)
		if angleParms[ijk][0] == 'harmonic':
			mcf.write('  %-8s  %8.1f  %8.2f' % angleParms[ijk])
		elif angleParms[ijk][0] == 'fixed':
			mcf.write('  %-8s  %8.2f' % angleParms[ijk])
		mcf.write('\n')

	if dihedralList:
		mcf.write('\n!Dihedral Format\n')
		mcf.write('!index i j k l type parameters\n' + 
							'!type="none"\n' + 
							'!type="CHARMM", parms=a0 a1 delta\n' + 
							'!type="OPLS", parms=c0 c1 c2 c3\n' + 
							'!type="harmonic", parms=force_constant equilibrium_dihedral\n')
	mcf.write('\n# Dihedral_Info\n')
	mcf.write(str(len(dihedralList)) + '\n')
	nDihedral = 0
	for dihedral in dihedralList:
		ijkl = tuple([int(i) for i in dihedral.split()])
		if ijkl not in dihedralParms: ijkl = ijkl[::-1]
		nDihedral = nDihedral + 1
		mcf.write('%-4d' % (nDihedral))
		mcf.write('  %-4d  %-4d  %-4d  %-4d' % ijkl)
		if dihedralParms[ijkl][0] == 'CHARMM':
			mcf.write('  %-8s  %8.3f  %8.1f  %8.1f' % dihedralParms[ijkl])
		elif dihedralParms[ijkl][0] == 'OPLS':
			mcf.write('  %-8s  %8.3f  %8.3f  %8.3f  %8.3f' % dihedralParms[ijkl])
		elif dihedralParms[ijkl][0] == 'harmonic':
			mcf.write('  %-8s  %8.1f  %8.2f' % dihedralParms[ijkl])
		elif dihedralParms[ijkl][0] == 'none':
			mcf.write('  %-8s' % dihedralParms[ijkl])
		mcf.write('\n')

	mcf.close()
	improperInfo(mcfFile)
	fragInfo(mcfFile)
	fragConnectivityInfo(mcfFile)


#*******************************************************************************
# FILE MANAGEMENT
#*******************************************************************************
infilename_type = check_type_infilename(configFile)
if infilename_type == 'cml':
	cml_to_pdb(configFile)
        configFile = configFile + '.pdb'

#*******************************************************************************
# MAIN PROGRAM BEGINS HERE
#*******************************************************************************
#INITIALIZE LISTS AND VARIABLES
genScannedList=[]
ringList=[]
fragList=[]
angleList=[]
dihedralList=[]
bondList=[]
bondList_atomType = []
angleList_atomType = []
dihedralList_atomType = []
spacing="    "
fragConn=[]
listofnames=[]
cyclic_ua_atom = True

#RING SCAN
initialatom=initialize(configFile)
if initialatom == '':
	tempfile = open('temporary.temp','w')
	tempfile.write('1')
	tempfile.close()
	initialatom = '1'
	
scan_result = scan(initialatom,[],False)
if cyclic_ua_atom==True:
	#Clean "ghost atom" from temporary_file
	tempfile2=open('temporary.temp2','w')
	tempfile1=open('temporary.temp','r')
	for line in tempfile1:
		this_line_new = line.split()
		if this_line_new[-1] == '0':		
			for each_atom in this_line_new[0:-1]:
				tempfile2.write(each_atom+"  ")
			tempfile2.write("\n")
		elif '0 1' in line:
			continue
		else:
			tempfile2.write(line)
	tempfile2.close()
	tempfile1.close()
	os.system('mv temporary.temp2 temporary.temp')

#IDENTIFY BONDS, FRAGMENTS, ANGLES, DIHEDRALS
bondID()
fragID()
angleID()
dihedralID()
fragConnectivity()

#PRESENT SCAN SUMMARY TO USER
print "\n\n*********Generation of Topology File*********\n"
print "Summary: "
print "There are " + str(len(bondList)) + " bonds identified." 

if cyclic_ua_atom == True and len(scan_result[1]) > 2:
	print "Cyclic united atom molecule with no branches"
else:
	print "There are " + str(len(ringList)) + " rings identified. These rings are: "

	for row in ringList:
		print row

	print "There are " + str(len(fragList)) + " fragments identified"
print "There are " + str(len(angleList)) + " angles identified"
print "There are " + str(len(dihedralList)) + " dihedrals identified"
if ffTemplate:
	#GENERATE BLANK FORCEFIELD TEMPLATE
	if os.path.isfile(ffFile) and not os.path.isfile(ffFile + '.BAK'):
		os.system('mv ' + ffFile + ' ' + ffFile + '.BAK')
	ffFileGeneration(configFile,ffFile)
	print 'Finished'
else:
	# GENERATE MCF FILE
	atomList, atomParms, = readPdb(configFile)
	# Read parms
	if ffFileType == 'native':
		atomParms, bondParms, angleParms, dihedralParms = \
			readNative(ffFile, atomParms)
	elif ffFileType == 'gromacs':
		atomParms, bondParms, angleParms, dihedralParms = \
			readGromacs(ffFile, atomParms, {}, {}, {})
	
	# Check parms
	# The *Parms dictionaries may contain entries with atomType indices.
	# We need to make sure all the atomNumber indices are populated.
	atomParms, bondParms, angleParms, dihedralParms = \
		checkParms(atomList,  bondList,  angleList,  dihedralList,
		           atomParms, bondParms, angleParms, dihedralParms)

	# Got all the parms we need? write Mcf.
	writeMcf(configFile, mcfFile, 
	         atomList, bondList, angleList, dihedralList, ringList,
	         atomParms, bondParms, angleParms, dihedralParms)

os.system("rm temporary.temp")
if infilename_type == 'cml':
	os.system("rm " + configFile)
