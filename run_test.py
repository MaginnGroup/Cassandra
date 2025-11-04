# In run_test.py
import sys
import os

# --- Setup Python Path ---
# This tells Python "look in the 'python' folder for modules to import"
# We do this so we can find our new 'api' module.
current_dir = os.path.dirname(os.path.abspath(__file__))
python_dir = os.path.join(current_dir, "python")
sys.path.append(python_dir)

# --- Import Your Module ---
# Now this import will work
try:
    import api
except ImportError:
    print(f"Error: Could not find the 'api.py' module in {python_dir}")
    print("Please make sure it exists.")
    sys.exit(1)

print("Successfully imported 'api' module.")

# --- Define Paths ---
# We'll test with the SPC water example
# We build the paths relative to this script's location.
example_inp = os.path.join(current_dir, "Examples/NVT/water_spc/nvt.inp")

# IMPORTANT: You must have a compiled 'cassandra.exe' for this to work.
# We'll assume it's in the 'Src' directory.
cassandra_exe = os.path.join(current_dir, "Src/cassandra_gfortran.exe")

# --- Check if Files Exist ---
if not os.path.exists(example_inp):
    print(f"Error: Example input file not found at {example_inp}")
    sys.exit(1)

if not os.path.exists(cassandra_exe):
    print(f"Error: Cassandra executable not found at {cassandra_exe}")
    print("Please compile Cassandra first by running 'make' in the 'Src' directory.")
    sys.exit(1)

print(f"Found Cassandra executable: {cassandra_exe}")
print(f"Found example input: {example_inp}")

# --- Run the Test ---
print("\n=== CALLING api.run_cassandra() ===")
try:
    api.run_cassandra(
        input_file_path=example_inp,
        cassandra_exe_path=cassandra_exe
    )
    print("=== TEST SCRIPT FINISHED ===")

except Exception as e:
    print(f"An error occurred while running api.run_cassandra: {e}")