# This is the whole testSuite. So there are many different modules to this testSuite. Each Test will test something different described above the code which allows that test to run. Each individual test is saves as a separate python script. The individual python scripts are not added to the testSuite until they are fully functional. 
# First, the subprocess module needs to be imported into this python script. 
import subprocess

# The following are the tests in the testSuite:

# Test 1: LJ Starting Energy:
#	This test is designed to... 
subprocess.call(["python", "Test1_LJ_StartEnergy.py"])


# Test 2: Mie_Starting Energy:
#	This this is designed to... 
subprocess.call(["python", "Test2_Mie_StartEnergy.py"])

# And that worked, so pretty much, I'm a genius. No biggie.

# Test 3: Angle Starting Energy: 
# 	This test is designed to... 
subprocess.call(["python", "Test3_Angle_StartEnergy.py"])

# Test 4: Dihedral Starting Energy:
	# This test is designed to... 
subprocess.call(["python", "Test4_Dihedral_StartEnergy.py"])

# Test 5: Improper Energy:
subprocess.call(["python", "Test5_ImproperEnergy.py"])

# Test 6: NIST Lennard Jones Energy:
	# This test is designed to... 
subprocess.call(["python", "Test6_NIST_LennardJonesEnergy.py"])

# Test 7: Charge Energy (Water SPCE): 
	# This test is designed to...
subprocess.call(["python", "Test7_ChargeEnergy.py"])
