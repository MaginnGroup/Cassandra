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
# This script helps the user to set up a Molecular Connectivity File (MCF). 
# The steps for using are as follows:
#
# 1) Generate a PDB of CML file
# 2) Run
#      > mcfgen.py -fffile pdborcml.xxx molecule.ff
# 3) Fill out molecule.ff with literature force field values. Look after units
#    and functional force field forms.
# 4) Run
#      > mcfgen.py -cassandra pdborcml.xxx molecule.mcf molecule.ff
# 5) Check molecule.mcf
#
#********************************************************************************

#!/usr/bin/env python
import sys, os, argparse, linecache, re

#VARIABLES
pt=[['1','1.0079','Hydrogen','H'], 
['2','4.0026','Helium','He'], 
['3','6.941','Lithium','Li'], 
['4','9.0122','Beryllium','Be'], 
['5','10.811','Boron','B'], 
['6','12.0107','Carbon','C'], 
['7','14.0067','Nitrogen','N'], 
['8','15.9994','Oxygen','O'], 
['9','18.9984','Fluorine','F'], 
['10','20.1797','Neon','Ne'], 
['11','22.9897','Sodium','Na'], 
['12','24.305','Magnesium','Mg'], 
['13','26.9815','Aluminum','Al'], 
['14','28.0855','Silicon','Si'], 
['15','30.9738','Phosphorus','P'], 
['16','32.065','Sulfur','S'], 
['17','35.453','Chlorine','Cl'], 
['18','39.948','Argon','Ar'], 
['19','39.0983','Potassium','K'], 
['20','40.078','Calcium','Ca'], 
['21','44.9559','Scandium','Sc'], 
['22','47.867','Titanium','Ti'], 
['23','50.9415','Vanadium','V'], 
['24','51.9961','Chromium','Cr'], 
['25','54.938','Manganese','Mn'], 
['26','55.845','Iron','Fe'], 
['27','58.9332','Cobalt','Co'], 
['28','58.6934','Nickel','Ni'], 
['29','63.546','Copper','Cu'], 
['30','65.39','Zinc','Zn'], 
['31','69.723','Gallium','Ga'], 
['32','72.64','Germanium','Ge'], 
['33','74.9216','Arsenic','As'], 
['34','78.96','Selenium','Se'], 
['35','79.904','Bromine','Br'], 
['36','83.8','Krypton','Kr'], 
['37','85.4678','Rubidium','Rb'], 
['38','87.62','Strontium','Sr'], 
['39','88.9059','Yttrium','Y'], 
['40','91.224','Zirconium','Zr'], 
['41','92.9064','Niobium','Nb'], 
['42','95.94','Molybdenum','Mo'], 
['43','98','Technetium','Tc'], 
['44','101.07','Ruthenium','Ru'], 
['45','102.9055','Rhodium','Rh'], 
['46','106.42','Palladium','Pd'], 
['47','107.8682','Silver','Ag'], 
['48','112.411','Cadmium','Cd'], 
['49','114.818','Indium','In'], 
['50','118.71','Tin','Sn'], 
['51','121.76','Antimony','Sb'], 
['52','127.6','Tellurium','Te'], 
['53','126.9045','Iodine','I'], 
['54','131.293','Xenon','Xe'], 
['55','132.9055','Cesium','Cs'], 
['56','137.327','Barium','Ba'], 
['57','138.9055','Lanthanum','La'], 
['58','140.116','Cerium','Ce'], 
['59','140.9077','Praseodymium','Pr'], 
['60','144.24','Neodymium','Nd'], 
['61','145','Promethium','Pm'], 
['62','150.36','Samarium','Sm'], 
['63','151.964','Europium','Eu'], 
['64','157.25','Gadolinium','Gd'], 
['65','158.9253','Terbium','Tb'], 
['66','162.5','Dysprosium','Dy'], 
['67','164.9303','Holmium','Ho'], 
['68','167.259','Erbium','Er'], 
['69','168.9342','Thulium','Tm'], 
['70','173.04','Ytterbium','Yb'], 
['71','174.967','Lutetium','Lu'], 
['72','178.49','Hafnium','Hf'], 
['73','180.9479','Tantalum','Ta'], 
['74','183.84','Tungsten','W'], 
['75','186.207','Rhenium','Re'], 
['76','190.23','Osmium','Os'], 
['77','192.217','Iridium','Ir'], 
['78','195.078','Platinum','Pt'], 
['79','196.9665','Gold','Au'], 
['80','200.59','Mercury','Hg'], 
['81','204.3833','Thallium','Tl'], 
['82','207.2','Lead','Pb'], 
['83','208.9804','Bismuth','Bi'], 
['84','209','Polonium','Po'], 
['85','210','Astatine','At'], 
['86','222','Radon','Rn'], 
['87','223','Francium','Fr'], 
['88','226','Radium','Ra'], 
['89','227','Actinium','Ac'], 
['90','232.0381','Thorium','Th'], 
['91','231.0359','Protactinium','Pa'], 
['92','238.0289','Uranium','U'], 
['93','237','Neptunium','Np'], 
['94','244','Plutonium','Pu'], 
['95','243','Americium','Am'], 
['96','247','Curium','Cm'], 
['97','247','Berkelium','Bk'], 
['98','251','Californium','Cf'], 
['99','252','Einsteinium','Es'], 
['100','257','Fermium','Fm'], 
['101','258','Mendelevium','Md'], 
['102','259','Nobelium','No'], 
['103','262','Lawrencium','Lr'], 
['104','261','Rutherfordium','Rf'], 
['105','262','Dubnium','Db'], 
['106','266','Seaborgium','Sg'], 
['107','264','Bohrium','Bh'], 
['108','277','Hassium','Hs'], 
['109','268','Meitnerium','Mt']]


   
#logfile=open('logfile.log','w')
genscannedlist=[]
ringlist=[]
fraglist=[]
anglelist=[]
dihedrallist=[]
bondlist=[]
spacing="    "
frag_conn=[]
dihedrallist_atomtypes=[]
anglelist_atomtypes=[]
bondlist_atomtypes=[]
listofnames=[]
maffa_list_atom=[]
lammps = False
cassandra = False
ff_flag = False

#FUNCTION DEFINITIONS

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
		cml_atom_info.append(re.findall('"([^"]*)"',linecache.getline(infilename, i)))
			
	for i in xrange(cml_start_bonds+1, cml_end_bonds):
		cml_bond_info.append(re.findall('"([^"]*)"',linecache.getline(infilename, i))[0].split())


	filepdb = open(infilename+'.pdb','w')

	for line_nbr, line in enumerate(cml_atom_info):
		filepdb.write('HETATM       '+str(line_nbr+1)+'       '+line[1]+'       x       y       z       '+line[1]+'       '+line[5]+"\n")
		
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

def lookup(atomlook):
    tempfile = open ('temporary.temp', 'r')
    for line in tempfile:
        linestring = line.split()
        if linestring[0]==atomlook:
           line=line.split()
           return line
    tempfile.close()

def initialize(infilename,outfilename):
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

	if cyclic_ua_atom == True and atom=='' : #This means molecule has at least one cyclic united atom model
		tempfile.close()
		atom = '0'

		#Rewrite tempfile to include a dummy atom 0

		tempfile2=open('temporary.temp2','w')
		tempfile1=open('temporary.temp','r')
		for line in tempfile1:
			this_line_new = line.split()
			if this_line_new[0]=='1':
				this_line_new.append('0')
			for each_atom in this_line_new:
				tempfile2.write(each_atom+'  ')
			tempfile2.write('\n')
		tempfile2.write('0 1')

		tempfile2.close()
		tempfile1.close()
		os.system('mv temporary.temp2 temporary.temp')


	return atom



def isring(scannedlist, atom):
	for i in xrange(0,len(scannedlist)):
		if scannedlist[i]==atom:
			return True

def cleangenlist(branchlist):
	for i in xrange(0,len(genscannedlist)):
		if genscannedlist[i]==branchlist[0]:
			for j in xrange(0,len(branchlist)):
				del(genscannedlist[-1])
			return

def isolatering(thelist):
	location=[]
	seen=set()
	seen_add=seen.add
	seen_twice=set(x for x in thelist if x in seen or seen_add(x))
	repeatedatom=list(seen_twice)

	for index, item in enumerate(thelist):
		if item==repeatedatom[0]:
			location.append(index)
	ringlist.append(thelist[location[0]:location[1]])
	#logfile.write("The ring list is: " + str(ringlist) + "\n")
	cleangenlist(ringlist[-1])
	
		
def isalreadyinring(atom):
	for row in ringlist:
		for i in range(0,len(row)):
			if row[i]==atom:
				return True
	
def bond_identification():

	tempfile = open('temporary.temp','r')
	oldatom=[]
	for line in tempfile:
		atomlist1=line.split()
		newatom=atomlist1[0]
		for i in range(1, len(atomlist1)):
			if atomlist1[i] in oldatom:
				continue
			bondlist.append(atomlist1[0] + " " + atomlist1[i]+"\n")
			oldatom.append(newatom)
	tempfile.close()


def scan(newatom,oldline,scanning):
        scannedlist=[]
	#branchalreadyscanned=[0]
	while 1==1:


		#logfile.write("Atom " + newatom + "\n")

		if isalreadyinring(newatom)==True:

			#logfile.write("Atom " + newatom + " was already included in a ring\n\n")
			break

		if isring(genscannedlist, newatom)==True:

			#logfile.write("A ring was found because atom " + newatom + " was already scanned\n\n")
			genscannedlist.append(newatom)
			return "EndRing", scannedlist

		scannedlist.append(newatom)
		genscannedlist.append(newatom)
		newline=lookup(newatom)
		totalbondedatoms=len(newline)-1

		if totalbondedatoms>2:  #branch point?

			#logfile.write("A branch point was found in atom " + newatom + "\n")
			for j in xrange(1,totalbondedatoms+1):
			        if newline[j]==oldline[0]: #or newline[j]==branchalreadyscanned[0]: #If old atom in list is scanned
					continue
				newatom=newline[j]
				endtype, branchlist = scan(newatom, newline, True)
				#branchalreadyscanned=lookup(newatom)
				if endtype=="EndPoint":
					cleangenlist(branchlist)
				if endtype=="EndRing":
					isolatering(genscannedlist)
				#if endtype=="EndScan":
					
				
			break
	
   	        for i in xrange(1,totalbondedatoms+1):

			if scanning==True and totalbondedatoms==1:
				#logfile.write("End point found at " + newatom + "\n")
				return "EndPoint", scannedlist

			if scanning==False and totalbondedatoms==1: #Initial step
				scanning = True

			if len(oldline)==0:
				newatom=newline[i]
				oldline=newline
				continue

			if newline[i]==oldline[0]: #If old atom in list is scanned
				continue
	
			newatom=newline[i]
			oldline=newline
			break

	return "EndScan", scannedlist


def fragidentification():

	adjacentatoms=[]
	tempfile=open('temporary.temp','r')
	for eachring in ringlist:
		for eachatom in eachring:
			vector=lookup(eachatom)
			if len(vector)>3:
				adjacentatoms.append(list(set(vector)-set(eachring)))

		adjacentatoms1=filter(None,adjacentatoms)
		adjacentatoms = [x for sublist in adjacentatoms1 for x in sublist]
		fraglist.append(eachring+adjacentatoms)
		adjacentatoms=[]

	for line in tempfile:
		inaring=False
		linestring=line.split()
		if len(linestring)>2:
			anchoratom=linestring[0]
		else:
			continue

		for row in ringlist:
			if anchoratom in set(row):
				inaring=True
                       
		if inaring==True:                                                             
			continue
		else:
			fraglist.append(linestring)

	tempfile.close()

def angleidentification():
	tempfile=open('temporary.temp','r')
	for line in tempfile:
		atomlist=line.split()
		if len(atomlist)<3:
			continue
		apex=atomlist[0]
		del(atomlist[0])
		for i in range(0,len(atomlist)-1):
			for j in range(i+1,len(atomlist)):
				anglelist.append(atomlist[i] + " " + apex + " " + atomlist[j])

	removedoublecounting("angles")



	
def dihedralidentification():
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
					dihedrallist.append(atom1 + " " +atom2 + " " +atom3 + " " +atom4)

	removedoublecounting("dihedrals")
	
				
def atom_info(infilename, outfilename): 
#THE CAR FILE IS GENERATED HERE!
	
	ifile = open(infilename, 'r')
	ofile = open(outfilename,'w')
#	carfile = open(outfilename + ".car",'w')
	ofile.write("!**********************************************************************************\n")
	ofile.write("!Molecular connectivity file for " + infilename + "\n")
	ofile.write("!**********************************************************************************\n")

#	ofile.write("\n# Species_Type\nSORBATE\n\n# 1st_Fragment_Ins_Style\nCOM\n")
#	ofile.write("\n# 1st_Fragment_Ins_Style\nCOM\n")
	ofile.write("\n# Atom_Info\n")

	matrix = forcefield("nonbonded","read")
	matrixcharge = forcefield("charges","read")

	ofile.write(str(len(listofnames))+"\n")
	
	for line in ifile:
		linestring = line.split()
	        if not line.strip():
        	        continue
		if linestring[0]=='HETATM' or linestring[0]=='ATOM':
			atomnumber=linestring[1]
        		atomname=linestring[2]+linestring[1]
			element = linestring[2]
			weight = tablelookup(element,"weight")
			vdwtype = "LJ"

			#carfile.write(atomname + spacing + "0.0" + spacing + "0.0" + spacing + "0.0\n")
			
			for rowmatrix in matrix:
				for item in rowmatrix[0]:
					if int(item)==int(atomnumber)-1:
						sigma = rowmatrix[2]
						epsilon = rowmatrix[3]
						for index,line_charge in enumerate(matrixcharge):
							if atomnumber == line_charge[0]:
								charge = line_charge[1] 

			in_a_ring = False
			for eachring in ringlist:
				if atomnumber in eachring:
					in_a_ring = True
#tlacaelel
        		if cyclic_ua_atom == True and len(ringlist)==0: #This means if there is only one cyclic united atom molecule
				in_a_ring=True

			if in_a_ring:
				ofile.write(atomnumber + spacing + atomname + spacing + element + spacing + weight + spacing + charge +spacing +  vdwtype + spacing + epsilon + spacing + sigma + spacing + "ring\n")
			else:
				ofile.write(atomnumber + spacing + atomname + spacing + element + spacing + weight + spacing + charge +spacing +  vdwtype + spacing + epsilon + spacing + sigma + "\n")
	
	ifile.close()
	ofile.close()
	#carfile.close()


def bond_info(outfilename):



	ofile = open(outfilename,'a')
	ofile.write("\n# Bond_Info\n")
	ofile.write(str(len(bondlist))+"\n")
	
	matrix = forcefield("bonded","read")
	

	for index, row in enumerate(bondlist):
		bond=row.split()
		for rowmatrix in matrix:
			for item in rowmatrix[0]:
				if index == item:
					length=rowmatrix[2]
					constant = rowmatrix[3]
		ofile.write(str(index+1) + spacing + bond[0] + spacing + bond[1] + spacing + "fixed" + spacing + length+"\n")

	ofile.close()


def angle_info(outfilename):

	ofile = open(outfilename,'a')

	ofile.write("\n# Angle_Info\n")
	ofile.write(str(len(anglelist))+"\n")

	matrix = forcefield("angles","read") #matrix=[repeatedlines, angle type, force constant, angle)

	for index, row in enumerate(anglelist):
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



def dihedral_info(outfilename):

	ofile = open(outfilename,'a')

	ofile.write("\n# Dihedral_Info\n")
	ofile.write(str(len(dihedrallist))+"\n")

	matrix=forcefield("dihedrals","read")   #matrix=[repeatedlines, dihedral type, K, n, gamma)

	for index, row in enumerate(dihedrallist):
		dihedral=row.split()
		for rowmatrix in matrix:
			for item in rowmatrix[0]:
				if index == item:
					if dihedraltype == "CHARMM": 
						K = rowmatrix[2]
						n = rowmatrix[3]
						gamma = rowmatrix[4]
						ofile.write(str(index+1) + spacing + dihedral[0] + spacing + dihedral[1] + spacing + dihedral[2] + spacing + dihedral[3] + spacing + dihedraltype + spacing + K + spacing + n + spacing + gamma + "\n")
					elif dihedraltype == "OPLS":
						a0 = rowmatrix[2]
						a1 = rowmatrix[3]
						a2 = rowmatrix[4]
						a3 = rowmatrix[5]
						ofile.write(str(index+1) + spacing + dihedral[0] + spacing + dihedral[1] + spacing + dihedral[2] + spacing + dihedral[3] + spacing + dihedraltype + spacing + a0 + spacing + a1 + spacing + a2 + spacing + a3 + "\n")
					elif dihedraltype == "harmonic":
						k = rowmatrix[2]
						phi = rowmatrix[3]
						ofile.write(str(index+1) + spacing + dihedral[0] + spacing + dihedral[1] + spacing + dihedral[2] + spacing + dihedral[3] + spacing + dihedraltype + spacing + k + spacing + phi+ "\n")



				

	ofile.close()
		

def fragment_info(outfilename):
	ofile=open(outfilename,'a')
	ofile.write("\n# Fragment_Info\n")

	global fraglist

#tlacaelel
        if cyclic_ua_atom == True and len(ringlist)==0: #This means if there is only one cyclic united atom molecule
		fraglist=[]
		vector = []
		for i in xrange(len(listofnames)):
			vector.append(str(i+1))
		fraglist.append(vector)

	if len(fraglist)==0 and len(listofnames)==1:
		ofile.write("1\n1 1 1")
		return
	elif len(fraglist)==0 and len(listofnames)==2:
		ofile.write("1\n1 2 1 2")
		return


	ofile.write(str(len(fraglist))+"\n")
	for index, row in enumerate(fraglist):
		ofile.write(str(index+1) + spacing + str(len(row)))
		for item in row:
			ofile.write(spacing+item)
		ofile.write("\n")
	ofile.close()





def tablelookup(element, thingtolook):

	for i in xrange(0,108):
		line = pt[i]
		if line[3]==element:
			if thingtolook=="weight":
				return line[1]
			elif thingtolook=="name":
				return line[2]
			elif thingtolook=="atomic_number":
				return line[0]
		

def removedoublecounting(thingtoclean):

	if thingtoclean=="dihedrals":
		linesrepeated=[]
		cleanlist=[]
		for index,row in enumerate(dihedrallist):
			set1=set(row.split())
			for i in range(index+1, len(dihedrallist)):
				set2=set(dihedrallist[i].split())
				if set1==set2:
					linesrepeated.append(i)
	
		for index,row in enumerate(dihedrallist):
			if index in linesrepeated:
				continue
			cleanlist.append(row)

		del(dihedrallist[0:len(dihedrallist)])	
		for i in range(0,len(cleanlist)):
			dihedrallist.append(cleanlist[i])

	if thingtoclean=="angles":
		linesrepeated=[]
		cleanlist=[]
		for index,row in enumerate(anglelist):
			set1=set(row.split())
			for i in range(index+1, len(anglelist)):
				set2=set(anglelist[i].split())
				if set1==set2:
					linesrepeated.append(i)
	
		for index,row in enumerate(anglelist):
			if index in linesrepeated:
				continue
			cleanlist.append(row)

		del(anglelist[0:len(anglelist)])	
		for i in range(0,len(cleanlist)):
			anglelist.append(cleanlist[i])

	if thingtoclean=="frag_conn":
		linesrepeated=[]
		cleanlist=[]
		for index,row in enumerate(frag_conn):
			set1=set(row.split())
			for i in range(index+1, len(frag_conn)):
				set2=set(frag_conn[i].split())
				if set1==set2:
					linesrepeated.append(i)
	
		for index,row in enumerate(frag_conn):
			if index in linesrepeated:
				continue
			cleanlist.append(row)
		
		del(frag_conn[0:len(frag_conn)])	
		for i in range(0,len(cleanlist)):
			frag_conn.append(cleanlist[i])

def fragment_connectivity():
	group=[]
	j=0

#tlacaelel
        if cyclic_ua_atom == True and len(ringlist)==0:
		frag_conn.append("0")
		return

	for index1, row in enumerate(fraglist):
		for j in range(0,len(row)-1):
			for i in range(j+1,len(row)):
				element1=row[j]
				element2=row[i]
				group.append(element1)
				group.append(element2)
				for index2, otherrow in enumerate(fraglist):
					#if otherrow in oldrows:
					#	continue
					a=set(group)
					b=set(otherrow)
					if a<=b and index1!=index2:
						frag_conn.append(str(index1+1) + spacing + str(index2+1))
				group=[]

	removedoublecounting("frag_conn")

			
			
def fragment_connectivity_info(outfilename):
	ofile=open(outfilename,'a')
	ofile.write("\n# Fragment_Connectivity\n")
	ofile.write(str(len(frag_conn))+"\n")
	for index, row in enumerate(frag_conn):
		ofile.write(str(index+1) + spacing + row)
		ofile.write("\n")
	ofile.write("\n\nEND")
	ofile.close()

def improper_info(outfilename):
	ofile=open(outfilename,'a')
	ofile.write("\n# Improper_Info\n")
	ofile.write("0\n")
	ofile.close()

def fffile_generation(infilename,outfilename):

	print "\n\n*********Force field template file generation*********\n"

	print "Reading Modified PDB File..."
	ifile = open(infilename,'r')
	atomtypenames=[]
	for line in ifile:
	        if not line.strip():
	                continue
		linestring = line.split()
		if linestring[0]=='HETATM' or linestring[0] == 'ATOM':
			atomnumber = linestring[1]
			name = linestring[-1]
			listofnames.append(atomnumber + " " + name)
			atomtypenames.append(name)
	ifile.close()
	atomtypenames = list(set(atomtypenames))
	numatomtypes = str(len(atomtypenames))

	global dihedraltype
	if len(dihedrallist) > 0:
	        dihedraltype = raw_input("Enter the dihedral functional form (CHARMM/OPLS/harmonic): ")
	else:
		dihedraltype = "NONE" #This is just a default. It will not affect molecules with no dihedrals.

	#IDENTIFY DIFFERENT TYPES OF ANGLES AND DIHEDRALS BASED ON ATOM TYPES
	#EACH ROW IN DIHEDRALLIST_ATOMTYPES HAS EACH DIHEDRAL EXPRESSED IN ATOM TYPES TERMS
	#EACH ROW IN ANGLELIST_ATOMTYPES HAS EACH ANGLE EXPRESSED IN ATOM TYPES TERMS


	for i, rowdihedral in enumerate(dihedrallist):
		rowdihedralstr=rowdihedral.split()
		dihedrallist_atomtypes.append(lookupinlist(listofnames, rowdihedralstr))

	for i, rowangle in enumerate(anglelist):
		rowanglestr=rowangle.split()
		anglelist_atomtypes.append(lookupinlist(listofnames, rowanglestr))

	for i, rowbond in enumerate(bondlist):
		rowbondstr = rowbond.split()
		bondlist_atomtypes.append(lookupinlist(listofnames, rowbondstr))

	ff_file=open(outfilename,'a')
	ff_file.write("atomtypes\n"+numatomtypes+ "\n\n")
	ff_file.write("begin atom-atomtype\n")
	for line in listofnames:
		ff_file.write(line + "\n")
	ff_file.write("end atom-atomtype\n\n")

	ff_file.write("dihedraltype "+dihedraltype+"\n\n")

	ff_file.close()

	forcefield("nonbonded","write")
	forcefield("bonded","write")
	forcefield("angles","write")
	forcefield("dihedrals","write")
	forcefield("charges","write")


def forcefield(degreeoffreedom,action):

	#Answer will be a matrix whose rows, given by linematrix, will be:
	#[list of atoms repeated, type of atom, sigma, epsilon]

	global dihedraltype

	answer=[]
	linematrix=[]

        #If the script will read from the ff file, the list that relates atom name and atom type is still unknown. 
	#The following lines read the atom-atomtype section in the file, and save it in listofnames. 
	#listofnames=[atomnumber atomtype]

	if action == "read":
		
		listofnames[:]=[]
		lineindex=0
		flag=0
		ff_file=open(fffile,'r')
		for linea in ff_file:
			lineindex=lineindex+1
			if "begin atom-atomtype" in linea:
				flag = 1
				continue

			if "end atom-atomtype" in linea:
				flag=0
				break	

			if flag==1:
				listofnames.append(linea[0:-1])


	if degreeoffreedom=="nonbonded":
	#The following lines will read listofnames and identify which atoms numbers share the same atom type. 
	#One variable is produced: Linesrepeated, which contains those atoms that share the same atom type.
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

				#Using the forcefield file, we need to build a matrix with the following form
				#MATRIX == #[list of atoms repeated, type of atom, sigma, epsilon]
				#To do that, wee need to find where the "nonbonded" labels are are in the ff
				#file so we can get sigma and epsilon

				linewherenonbonded=[]
				ff_file=open(fffile,'r')
				for index, line in enumerate(ff_file):
					if "nonbonded" in line:
						linewherenonbonded.append(index+1)
				ff_file.close()

			
				#Now, the following three lines will contain the elements of the matrix

				for i in linewherenonbonded:
					vartemporal = linecache.getline(fffile, i+1).split()[0]
					if vartemporal==rowvector1[1]:
						linematrix.append(linecache.getline(fffile, i+1)[0:-1])
						linematrix.append(linecache.getline(fffile, i+2)[6:-1])
						linematrix.append(linecache.getline(fffile, i+3)[8:-1])

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


	
				ff_file=open(outfilename,'a')
				ff_file.write("nonbonded\n")
				
				ff_file.write(row.split()[1] + "\n")
				ff_file.write("Sigma \n")
				ff_file.write("Epsilon \n\n")
				linesrepeated=[]
				ff_file.close()	


		


	if degreeoffreedom=="bonded":

		if action == "read":

			for i, rowbond in enumerate(bondlist):
				rowbondstr = rowbond.split()
				bondlist_atomtypes.append(lookupinlist(listofnames, rowbondstr))

			linesrepeated=[]
			rowalreadyscanned=[]

			for index,row in enumerate(bondlist_atomtypes):
				flag=0
				for oldsets in rowalreadyscanned:
					if oldsets==row:
						flag=1
					else:
						reverse = [j for j in reversed(oldsets)]
						if row==reverse:
							flag=1
				if flag==1: continue

				for i in range(index+1, len(bondlist_atomtypes)):
					if row==bondlist_atomtypes[i]:
						linesrepeated.append(i)
					else:
						reverse = [j for j in reversed(bondlist_atomtypes[i])]
						if row==reverse:
							linesrepeated.append(i)
	
				linesrepeated.append(index)
				rowalreadyscanned.append(row)

				linematrix.append(linesrepeated)

				#We need to find where the "bonds" labels are are in the ff
				#file so we can get K and angle

				linewherebonded=[]
				ff_file=open(fffile,'r')
				for index, line in enumerate(ff_file):
					if "bonds" in line:
						linewherebonded.append(index+1)
				ff_file.close()
			
				#Now, the following two lines will contain the elements of the matrix

				for i in linewherebonded:

					#Read from fffile which bondtype we'll work with (e.g. A-A, A-B)
					vartemporal = linecache.getline(fffile, i+1)[0:-1].split()

					if set(vartemporal)==set(row):
						linematrix.append(linecache.getline(fffile, i+1)[0:-1])
						linematrix.append(linecache.getline(fffile, i+2)[6:-1])
						linematrix.append(linecache.getline(fffile, i+3)[8:-1])

				answer.append(linematrix)
				linematrix=[]		
				linesrepeated=[]


		if action=="write":

			linesrepeated=[]
			rowalreadyscanned=[]

			for index,row in enumerate(bondlist_atomtypes):
				flag=0
				for oldsets in rowalreadyscanned:
					if oldsets==row:
						flag=1
					else:
						reverse = [j for j in reversed(oldsets)]
						if row==reverse:
							flag=1
				if flag==1: continue

				for i in range(index+1, len(bondlist_atomtypes)):
					if row==bondlist_atomtypes[i]:
						linesrepeated.append(i)
					else:
						reverse = [j for j in reversed(bondlist_atomtypes[i])]
						if row==reverse:
							linesrepeated.append(i)
	
				linesrepeated.append(index)
				rowalreadyscanned.append(row)

				ff_file=open(outfilename,'a')
				ff_file.write("bonds\n")
				ff_file.write(str(row[0]) + " " + str(row[1]) + "\n")
				ff_file.write("Length \n")
				ff_file.write("Constant \n\n")
				linesrepeated=[]
				ff_file.close()


	if degreeoffreedom=="angles":

		if action == "read":
			for i, rowangle in enumerate(anglelist):
				rowanglestr=rowangle.split()
				anglelist_atomtypes.append(lookupinlist(listofnames, rowanglestr))

			linesrepeated=[]
			rowalreadyscanned=[]

			for index,row in enumerate(anglelist_atomtypes):
				flag=0
				for oldsets in rowalreadyscanned:
					if oldsets==row:
						flag=1
					else:
						reverse = [j for j in reversed(oldsets)]
						if row==reverse:
							flag=1
				if flag==1: continue

				for i in range(index+1, len(anglelist_atomtypes)):
					if row==anglelist_atomtypes[i]:
						linesrepeated.append(i)
					else:
						reverse = [j for j in reversed(anglelist_atomtypes[i])]
						if row==reverse:
							linesrepeated.append(i)
	
				linesrepeated.append(index)
				rowalreadyscanned.append(row)


				linematrix.append(linesrepeated)

				#We need to find where the "angles" labels are are in the ff
				#file so we can get K and angle

				linewhereangle=[]
				ff_file=open(fffile,'r')
				for index, line in enumerate(ff_file):
					if "angles" in line:
						linewhereangle.append(index+1)
				ff_file.close()
			
				#Now, the following two lines will contain the elements of the matrix

				for i in linewhereangle:

					#Read from fffile which angle type we'll work with (e.g. A-A-B, A-B-A)
					vartemporal = linecache.getline(fffile, i+1)[0:-1].split()

					if vartemporal==row:
						linematrix.append(linecache.getline(fffile, i+1)[0:-1])
						linematrix.append(linecache.getline(fffile, i+2)[6:-1])
						linematrix.append(linecache.getline(fffile, i+3)[8:-1])

				answer.append(linematrix)
				linematrix=[]		
				linesrepeated=[]


		if action=="write":

			linesrepeated=[]
			rowalreadyscanned=[]

			for index,row in enumerate(anglelist_atomtypes):
				flag=0
				for oldsets in rowalreadyscanned:
					if oldsets==row:
						flag=1
					else:
						reverse = [j for j in reversed(oldsets)]
						if row==reverse:
							flag=1
				if flag==1: continue

				for i in range(index+1, len(anglelist_atomtypes)):
					if row==anglelist_atomtypes[i]:
						linesrepeated.append(i)
					else:
						reverse = [j for j in reversed(anglelist_atomtypes[i])]
						if row==reverse:
							linesrepeated.append(i)
	
				linesrepeated.append(index)
				rowalreadyscanned.append(row)


				ff_file=open(outfilename,'a')
				ff_file.write("angles\n")
				ff_file.write(row[0] + " " + row[1] + " " +row[2] + "\n")
				ff_file.write("Angle \n")
				ff_file.write("Constant \n\n")
				linesrepeated=[]
				ff_file.close()


	if degreeoffreedom=="dihedrals":

		if action=="read":

			for i, rowdihedral in enumerate(dihedrallist):
				rowdihedralstr=rowdihedral.split()
				dihedrallist_atomtypes.append(lookupinlist(listofnames, rowdihedralstr))
		
			linesrepeated=[]
			rowalreadyscanned=[]

			for index,row in enumerate(dihedrallist_atomtypes):
				flag=0
				for oldsets in rowalreadyscanned:
					if oldsets==row:
						flag=1
					else:
						reverse = [j for j in reversed(oldsets)]
						if row==reverse:
							flag=1
				if flag==1: continue

				for i in range(index+1, len(dihedrallist_atomtypes)):
					if row==dihedrallist_atomtypes[i]:
						linesrepeated.append(i)
					else:
						reverse = [j for j in reversed(dihedrallist_atomtypes[i])]
						if row==reverse:
							linesrepeated.append(i)
	
				linesrepeated.append(index)
				rowalreadyscanned.append(row)


				linematrix.append(linesrepeated)

				#We need to find where the "dihedrals" labels are are in the ff
				#file so we can get K and angle

				linewheredihedral=[]
				ff_file=open(fffile,'r')
				for index, line in enumerate(ff_file):
					if "dihedrals" in line:
						linewheredihedral.append(index+1)
					if line[0:12] == "dihedraltype":
						dihedraltype = line.split()[1]
				ff_file.close()
			

				for i in linewheredihedral:

					vartemporal = linecache.getline(fffile, i+1)[0:-1].split()

					if vartemporal==row:
						if dihedraltype == "CHARMM":
							linematrix.append(linecache.getline(fffile, i+1)[0:-1])
							linematrix.append(linecache.getline(fffile, i+2)[1:-1])
							linematrix.append(linecache.getline(fffile, i+3)[1:-1])
							linematrix.append(linecache.getline(fffile, i+4)[5:-1])
						elif dihedraltype == "OPLS":
                                                        linematrix.append(linecache.getline(fffile, i+1)[0:-1])
                                                        linematrix.append(linecache.getline(fffile, i+2)[2:-1])
                                                        linematrix.append(linecache.getline(fffile, i+3)[2:-1])
                                                        linematrix.append(linecache.getline(fffile, i+4)[2:-1])
                                                        linematrix.append(linecache.getline(fffile, i+5)[2:-1])
						elif dihedraltype == "harmonic":
                                                        linematrix.append(linecache.getline(fffile, i+1)[0:-1])
                                                        linematrix.append(linecache.getline(fffile, i+2)[2:-1])
                                                        linematrix.append(linecache.getline(fffile, i+3)[3:-1])
 

				answer.append(linematrix)
				linematrix=[]		
				linesrepeated=[]


		if action=="write":


			linesrepeated=[]
			rowalreadyscanned=[]

			for index,row in enumerate(dihedrallist_atomtypes):
				flag=0
				for oldsets in rowalreadyscanned:
					if oldsets==row:
						flag=1
					else:
						reverse = [j for j in reversed(oldsets)]
						if row==reverse:
							flag=1
				if flag==1: continue

				for i in range(index+1, len(dihedrallist_atomtypes)):
					if row==dihedrallist_atomtypes[i]:
						linesrepeated.append(i)
					else:
						reverse = [j for j in reversed(dihedrallist_atomtypes[i])]
						if row==reverse:
							linesrepeated.append(i)
	
				linesrepeated.append(index)
				rowalreadyscanned.append(row)



				ff_file=open(outfilename,'a')
				ff_file.write("dihedrals\n")
				ff_file.write(row[0] + " " + row[1] + " " +row[2] + " " + row[3] + "\n")

#				global dihedraltype

				if dihedraltype == "CHARMM":
					ff_file.write("K \n")
					ff_file.write("n \n")
					ff_file.write("Gamma \n\n")
					linesrepeated=[]
					ff_file.close()
				elif dihedraltype == "OPLS":
					ff_file.write("a0 \n")
					ff_file.write("a1 \n")
					ff_file.write("a2 \n")
					ff_file.write("a3 \n\n")
					linesrepeated=[]
					ff_file.close()
				elif dihedraltype == "harmonic":
                                        ff_file.write("K \n")
                                        ff_file.write("phi \n\n")
					linesrepeated=[]
					ff_file.close()


	if degreeoffreedom=="charges":

		if action == "write":
			for index,row in enumerate(listofnames):
				ff_file=open(outfilename,'a')
				ff_file.write("charge\n")
				ff_file.write(row.split()[0] + " \n\n")
	

			ff_file.write("end")
			ff_file.close()

		if action == "read":

			ff_file=open(fffile,'r')
			for index, line in enumerate(ff_file):
				if "charge" in line:
					answer.append(linecache.getline(fffile, index+2).split())
			ff_file.close()
		
			


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


#ARGUMENT PARSE

parser = argparse.ArgumentParser("This script helps to generate a Molecular Connectivity File (MCF) using a PDB or CML file. An intermediate .ff file must be filled out with the parameters found in the literature")
parser.add_argument('pdbfile', action="store")
#parser.add_argument('outfile', action="store")
#parser.add_argument('forcefile', nargs='?', action="store")

parser.add_argument('-cassandra', help="Generate MCF file for Cassandra", action="store_true")
parser.add_argument('-fffile', help="Generate a standard force field file", action="store_true")

results = parser.parse_args()


if results.pdbfile==None:
	print "Error: No PDB specified. Order: mcfgen.py [flag] pdbfile.pdb"
	exit()

basename = os.path.splitext(os.path.basename(results.pdbfile))[0]

infilename=results.pdbfile

if results.cassandra:
	cassandra=True
	outfilename=basename+".mcf"
	fffile = basename+".ff"

#	if results.forcefile==None:
#		print "Error: No Force File specified. Order: pdb2mcf.py pdbfile mcffile forcefile -cassandra"
#		exit()
#	else:
#		fffile = results.forcefile

elif results.fffile:
	ff_flag = True
	outfilename=basename+".ff"
#	fffile = results.forcefile



infilename_type = check_type_infilename(infilename)
if infilename_type == 'cml':
	cml_to_pdb(infilename)
        infilename = infilename+'.pdb'




cyclic_ua_atom = True
	
#MAIN PROGRAM BEGINS HERE

#RING SCAN
#logfile.write("******Ring scan begins******\n")
initialatom=initialize(infilename, outfilename)
if initialatom == '':
	tempfile = open('temporary.temp','w')
	tempfile.write('1')
	tempfile.close()
	initialatom = '1'
	

#logfile.write("Scan will start in atom " + str(initialatom) + "\n")
scan_result = scan(initialatom,[],False)


#logfile.write("******Ring scan ends******\n")
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


#BOND_IDENTIFICATION
bond_identification()

#FRAGMENT IDENTIFICATION
#logfile.write("******Identification of fragments*******\n")
fragidentification()
#logfile.write("There are "+str(len(fraglist))+" fragments in this molecule\n")
#for row in fraglist:
	#logfile.write(str(row) + "\n")
#logfile.write("******End of identification of fragments*******\n")

#ANGLE IDENTIFICATION
#logfile.write("******Identification of angles*******\n")
angleidentification()
#for row in anglelist:
	#logfile.write(str(row) + "\n")
#logfile.write("******End of angle identification*******\n")


#DIHEDRAL IDENTIFICATION
#logfile.write("******Dihedral identification******\n")
dihedralidentification()
#logfile.write("There are " + str(len(dihedrallist)) + " dihedrals\n")
#for row in dihedrallist:
	#logfile.write(row + "\n")
#logfile.write("******End of Dihedral identification******\n")

#FRAGMENT_CONNECTIVITY
fragment_connectivity()

#PRESENT SCAN SUMMARY TO USER


print "\n\n*********Generation of Topology File*********\n"
print "Summary: "
print "There are " + str(len(bondlist)) + " bonds identified." 

if cyclic_ua_atom == True and len(scan_result[1]) > 2:
 	print "Cyclic united atom molecule with no branches"
else:
	print "There are " + str(len(ringlist)) + " rings identified. These rings are: "

	for row in ringlist:
		print row

	print "There are " + str(len(fraglist)) + " fragments identified"
print "There are " + str(len(anglelist)) + " angles identified"
print "There are " + str(len(dihedrallist)) + " dihedrals identified"

#FORCE FIELD FILE/MCF/LAMMPS GENERATION

if ff_flag: 
	os.system("rm " + outfilename)
	fffile_generation(infilename,outfilename)
	print "Finished"

if cassandra:

	#ATOM_INFO SECTION
	#logfile.write("******Writing Atom_Info section******\n")
	atom_info(infilename, outfilename)
	#logfile.write("******End of Atom_Info section******\n")

	#BOND_INFO SECTION
	#logfile.write("******Writing Bond_Info section******\n")
	bond_info(outfilename)
	#for row in bondlist:
	#	logfile.write(row + "\n")
	#logfile.write("******End of Writing Bond_Info section******\n")

	#ANGLE_INFO SECTION
	angle_info(outfilename)

	#DIHEDRAL_INFO SECTION
	dihedral_info(outfilename)

	#FRAGMENT_INFO SECTION
	fragment_info(outfilename)

	#IMPROPER_INFO SECTION
	improper_info(outfilename)

	#FRAGMENT_CONNECTIVITY SECTION
	fragment_connectivity_info(outfilename)
	#logfile.close()
os.system("rm temporary.temp")
if infilename_type == 'cml':
        os.system("rm " + infilename)

