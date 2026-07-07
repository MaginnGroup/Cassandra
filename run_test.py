# Smoke test for the Cassandra Python wrapper.
# See python/README.md for architecture notes and usage examples.
import sys
import os

current_dir = os.path.dirname(os.path.abspath(__file__))
python_dir = os.path.join(current_dir, "python")
sys.path.append(python_dir)

try:
    import api
except ImportError:
    print(f"Error: Could not find the 'api.py' module in {python_dir}")
    print("Please make sure it exists.")
    sys.exit(1)

example_inp = os.path.join(current_dir, "Examples/NVT/water_spc/nvt.inp")
cassandra_exe = os.path.join(current_dir, "Src/cassandra_gfortran.exe")

print(f"Found example input: {example_inp}")
print(f"Using Cassandra executable: {cassandra_exe}")

print("\n=== CALLING api.run_cassandra() ===")
try:
    api.run_cassandra(
        input_file_path=example_inp,
        cassandra_exe_path=cassandra_exe,
    )
    print("=== TEST SCRIPT FINISHED ===")
except FileNotFoundError as exc:
    print(f"Error: {exc}")
    sys.exit(1)
except api.CassandraRunError as exc:
    print(f"Cassandra run failed: {exc}")
    sys.exit(1)
