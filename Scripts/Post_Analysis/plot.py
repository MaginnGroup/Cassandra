import numpy as np
import matplotlib.pyplot as plt
import argparse
import numpy as np
from argparse import RawTextHelpFormatter


"""
***************************************************************************************************
This script can be used to visualize your output property files. It  will print out the block
average properties and errors. For multiple property files, the output columns must be the
same property or this script will not work properly.

Requirements:
python 2.7.x
matplotlib 1.x

Example Usage:

python plot.py RUN_NAME.box1.prp1 RUN_NAME.box2.prp1 -skip NUM_LINES -setx -block NUM_BLOCKS  -lg 0

Help:
python plot.py -h

****************************************************************************************************
"""


parser = argparse.ArgumentParser("Plots data files for identical output columns", formatter_class=RawTextHelpFormatter)
parser.add_argument('file',nargs = '+', action='store',help='enter filename')
parser.add_argument('-skip', action ='store', help='Number lines to skip', default=1)
parser.add_argument('-block', action='store', help='Number of blocks for block averaging', default = 5)
parser.add_argument('-setx', action ='store_true', help='Change x-axis column from MC Steps')
parser.add_argument('-legend', action ='store',help='Legend positions:\n0. upper left (default)\n'
	+ '1. upper right\n2. lower left\n3. lower right\n', default = 0)
parser.add_argument('-marker', action='store',help='Example Options: \n'
	+'\'o\' - circles\n\'^\' - up triangles\n\'<\' - left triangle\n\'>\' - right triangle\n' 
	+ '\'x\' - X\'s\n\'.\' - dotted\n\'s\' - squares\n' 
	+'(see Matplotlib Documentation for more options)',default = '-')
 
args = parser.parse_args()
block = int(args.block)
file1 = open(args.file[0],'r')
file1.readline()
line1 = file1.readline()
data_label= line1.split()
line2 = file1.readline()
x_axis_flag = args.setx
marker = args.marker
x_axis = 1
stride = -1
conv_x = 1
conv_y = 1

#determine output property
print '\n'
for i in range(len(data_label)):
	if (data_label[i] == '#' or data_label[i] == 'MC_STEP' 
			or data_label[i] == 'MC_SWEEP' 
			or data_label[i] == '# MC_STEP'):
		stride = stride+1
	else:
		print str(i-stride)+'. ' + str(data_label[i])


print '\n Enter index of property  you would like to output:'+'\n'
y_axis = int(raw_input())+stride


#determine output units
output_options = ['Energy_Total','Enthalpy','Energy_Intra','Energy_Elec','Energy_LJ','Chemical_Potential',
				'Nmols','Subensemble','Volume','Pressure','Temperature','Density']
output_unit = ['(kJ/mol)-Ext','(kJ/mol)-Ext','(kJ/mol)-Ext','(kJ/mol)-Ext','(kJ/mol)-Ext','(kJ/mol)',
				 ' ',' ','(A^3)','(bar)','(K)','(kg/m^3)']
for i, option in enumerate(output_options):
	if data_label[y_axis] == option:
		output_ndx = i
		break
	else:
		output_ndx = 7 #empty

print '\n'+'Output:', data_label[y_axis], output_unit[output_ndx],'\n'
file1.close()

#if output is Density, convert value to kg/m^3
if data_label[y_axis]=="Density":
    print 'Input Molecular Weight (kg/mol) (i.e. 0.01802 for water):'
    MW = raw_input()
    conv_y = float(MW)/6.0221413E23*1E30

#set x-axis to different property
if x_axis_flag:
	print '\n'+'Enter index of x-axis:' +'\n'
	x_axis = int(raw_input())+stride
	outputx_ndx = np.where(np.in1d(output_options,data_label[x_axis]))
	if data_label[x_axis]=="Density":
    		print 'Input Molecular Weight (kg/mol) (i.e. 0.01802 for water):'
    		MW = raw_input()
    		conv_x = float(MW)/6.0221413E23*1E30

#print and plot
data = np.zeros(len(args.file))
legend = ['' for x in range(len(args.file))]
plt.figure(1)

for i in range(len(args.file)):
    legend[i] = args.file[i]

fmtfile = '%'+str(len(legend[0]))+'s'
fmtout= '%15s'
fmt1 = '%15.7f'
fmt2 = '%15.7f'
print '*'*(len(legend[0])+35)+'\n'
print fmtfile%('File Name'),fmtout%('Average'), fmtout%('Std Dev.')
print '*'*(len(legend[0])+35)+'\n'

for i in range(len(data)):
    data =  np.loadtxt(str(args.file[i]),skiprows=int(args.skip))
    
    x_data=data[:,x_axis-1]*conv_x
    y_data=data[:,y_axis-1]*conv_y
    y_data_block = np.split(y_data[len(y_data)%block:],block)
    avg_y_data_block =  np.average(y_data_block,1)
    print legend[i], fmt1%np.average(avg_y_data_block), fmt2%np.std(avg_y_data_block,ddof=1)
    plt.plot(x_data,y_data,marker)
    plt.hold(True)

print '*'*(len(legend[0])+35)+'\n'

if x_axis_flag:
	plt.xlabel(data_label[x_axis]+' '+output_unit[outputx_ndx])
else:
	plt.xlabel(data_label[1])

plt.ylabel(data_label[y_axis]+' '+output_unit[output_ndx])


#plot legend options
legend_pos = int(args.legend)

if legend_pos == 0:
	pos = 'upper left'
elif legend_pos == 1:
	pos = 'upper right'
elif legend_pos == 2:
	pos = 'lower left'
elif legend_pos == 3:
	pos = 'lower right'

plt.legend(legend,loc=pos,fontsize='small',numpoints=1)

plt.show()
