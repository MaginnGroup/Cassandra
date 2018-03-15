#!/usr/bin/env python

#*******************************************************************************
# PROPERTIES
#*******************************************************************************
import math
import numpy as np
import random

def LJ(r,epsilon,sigma):
	"""r, float = distance between two atoms
	epsilon, float = interaction energy
	sigma, float = vdw diameter"""
	sigByR6 = (sigma/r)**6
	sigByR12 = sigByR6**2
	u = 4 * epsilon * (sigByR12 - sigByR6)
	return u

def Mie(r,epsilon,sigma, n, m):
	"""r, float = distance between two atoms
	epsilon, float = interaction energy
	sigma, float = vdw diameter
	n, float = larger exponent
	m, float = smaller exponent"""
	sigByrn = (sigma/r)**n
	sigByrm = (sigma/r)**m
	u = ((n/(n-m))*((n/m)**(m/(n-m)))) * epsilon * ( (sigByrn - sigByrm) )
	return u

def MieCutAndShift(r,epsilon,sigma, n, m, cut):
	"""r, float = distance between two atoms
	epsilon, float = interaction energy
	sigma, float = vdw diameter
	n, float = larger exponent
	m, float = smaller exponent
	cut, float = r to cut potential at"""
	sigByrn = (sigma/r)**n
	sigByrm = (sigma/r)**m
	sigByrnCut = (sigma/cut)**n
	sigByrmCut = (sigma/cut)**m
	u = ((n/(n-m))*((n/m)**(m/(n-m)))) * epsilon * ( (sigByrn - sigByrm) - (sigByrnCut - sigByrmCut) )
	return u

def AngleEnergyHarmonic(force,nomangle, theta):
	"""force, float = bond angle force constant - K0 - units of K/rad2
	nomangle, float = nominal bond angle - Theta0 - units of degrees
	theta, float = actual bond angle - units of degrees
	Outputs energy in kJ/mol, rounded to the 3rd digit"""
	nomangle = nomangle * math.pi / 180
	theta = theta * math.pi / 180
	E = round(force*((theta-nomangle)**2)*0.008314,3)
	return E

def DihedEnergyOLPS(a0, a1, a2, a3, theta):
	"""
	theta, float = actual bond angle - units of degrees"""
	theta = theta * math.pi / 180
	E = round(a0 + a1*(1+math.cos(theta)) + a2*(1-math.cos(2*theta)) + a3*(1+math.cos(3*theta)),3)
	return E

def DihedEnergyCHARMM(a0, a1, delta, theta):
	"""
	theta, float = actual bond angle - units of degrees"""
	theta = theta * math.pi / 180
	delta = delta * math.pi / 180
	E =  round(a0*(1+math.cos(a1*theta-delta)),3)
	return E

def DihedEnergyHarmonic(force, nomangle, theta):
	"""force, float = bond angle force constant - K0 - units of K/rad2
	nomangle, float = nominal bond angle - Theta0 - units of degrees
	theta, float = actual bond angle - units of degrees"""
	nomangle = nomangle * math.pi / 180
	theta = theta * math.pi / 180
	E = round(force*((theta-nomangle)**2)*0.008314,3)
	return E

def Distance(x1,y1,z1, x2, y2, z2):
	"""xyz coordinates of two points"""
	d = math.sqrt((x2-x1)**2+(y2-y1)**2+(z2-z1)**2)
	return d

def xyzFromDistanceRandom(P1, d):
	"""P1 = vector of floats (x,y,z) for point 1
	d, float = distance desired"""
	P1 = list(P1)
	P1 = np.array([P1[0],P1[1],P1[2]])
	theta = random.random()*(2*math.pi)
	phi = random.random()*(math.pi)
	dvec = np.array([d*math.cos(theta)*math.sin(phi),d*math.sin(theta)*math.sin(phi),d*math.cos(phi)])
	P2 = P1 + dvec
	P2 = (P2[0],P2[1],P2[2])
	return P2

def xyzFromAngleRandom(P1, P2, d, angle):
	"""P1 = vector of floats (x,y,z) for point 1
	P2 = vector of floats (x,y,z) for point 2
	angle, float = angle desired
	d, float = distance between point 2 and point 3
	c, vector of floats for center of possible circle
	h, distance from c to P2
	n, unit vector from c to point 2 (parrallel to h)
	u, unit vector from c to a point on circle (Perpindicular to h)
	P3 = random point with desired angle
	returns xyz coordinates of third point with the inputted angle"""
	P1 = list(P1)
	P2 = list(P2)
	P1 = np.array([P1[0],P1[1],P1[2]])
	P2 = np.array([P2[0],P2[1],P2[2]])
	dist12 = Distance(P1[0],P1[1],P1[2],P2[0],P2[1],P2[2])
	theta = math.pi-(math.pi*angle/180)
	h = d*math.cos(theta)
	r = d*math.sin(theta)
	c = (P2-P1)*h/dist12+ P2
	n = ((c-P2)/(h))

	if (n[0] == 0):
		u = np.array([1,0,0])
	elif (n[1] == 0):
		u = np.array([0,1,0])
	elif (n[2] == 0):
		u = np.array([0,0,1])
	else:
		u = np.array([-2/(n[0]),1/(n[1]),1/(n[2])])
		u = u / (Distance(u[0],u[1],u[2],0,0,0))

	phi = random.random()*(2*math.pi)
	P3 = c + r*math.cos(phi)*u + r*math.sin(phi)*(np.cross(n,u)) 
	P3 = (P3[0],P3[1],P3[2])
	return P3

def xyzFromDihedral(P1, P2, P3, d, angle34, dihedangle):
	"""P1 = vector of floats (x,y,z) for point 1
	P2 = vector of floats (x,y,z) for point 2
	P3 = vector of floats (x,y,z) for point 3
	angle34, float = angle desired between P2, P3, P4
	dihedangle, float = angle desired for dihedral
	d, float = distance between point 3 and point 4
	c1, vector of floats for center of circle that created P1
	c4, vector of floats for center of possible circle for P4
	h, distance from c4 to point 3
	n, unit vector from c4 to point 3 (parrallel to h)
	u, unit vector from c4 to a point on circle c4 that has dihedral angle of 0
	P4 = random point with desired angle
	returns xyz coordinates of fourth point with the inputted dihedral angle"""
	P1 = list(P1)
	P2 = list(P2)
	P3 = list(P3)
	P1 = np.array([P1[0],P1[1],P1[2]])
	P2 = np.array([P2[0],P2[1],P2[2]])
	P3 = np.array([P3[0],P3[1],P3[2]])
	dist12 = Distance(P1[0],P1[1],P1[2],P2[0],P2[1],P2[2])
	dist23 = Distance(P2[0],P2[1],P2[2],P3[0],P3[1],P3[2])
	theta1 = math.acos(np.dot(P2-P1,P3-P2)/(dist12*dist23))
	h1 = dist12*math.cos(theta1)
	c1 = (P2-P3)*h1/dist23 + P2
	u = (P1-c1)/(dist12*math.sin(theta1))

	theta4 = math.pi-(math.pi*angle34/180)
	h4 = d*math.cos(theta4)
	r4 = d*math.sin(theta4)
	c4 = (P3-P2)*h4/dist23+ P3
	n = ((c4-P3)/(h4))

	phi = (math.pi*dihedangle/180)
	P3 = c4 + r4*math.cos(phi)*u + r4*math.sin(phi)*(np.cross(n,u)) 
	P3 = (P3[0],P3[1],P3[2])
	return P3

def replace_line(file_name, line_num, text):
	# This functon writes over specific lines of the input file created above using a function, called replace_line that is created below. 
	# This function takes three inputs: the name of the file where you would like to replace a line, the line number you would like to replace, and the text you want to replace the old text with. 
	lines = open(file_name, 'r').readlines()
	lines[line_num] = text
	out = open(file_name, 'w')
	out.writelines(lines)
	out.close()

def checkLastTest(overwrite,file_name):
	# This function checks the results of the previous test and then overwrites for the next test to populate the results if overwrite = true
	with open(file_name) as results:
		BoolResult = results.readlines()
		if (BoolResult[0] == 'True'):
			Passed=1
		else:
			Passed=0
	if overwrite:
		LastTest = open(file_name,"w")
		LastTest.write(' ')
		LastTest.close()
	return Passed

def erf(x):
	# Error function
	# constants
	a1 =  0.254829592
	a2 = -0.284496736
	a3 =  1.421413741
	a4 = -1.453152027
	a5 =  1.061405429
	p  =  0.3275911

	# Save the sign of x
	sign = 1
	if x < 0:
		sign = -1
	x = abs(x)

	# A & S 7.1.26
	t = 1.0/(1.0 + p*x)
	y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*np.exp(-x*x)

	return sign*y
