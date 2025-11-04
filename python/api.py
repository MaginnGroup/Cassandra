# In python/api.py
import subprocess
import os
import shlex

def run_cassandra(input_file_path: str, cassandra_exe_path: str = "cassandra.exe"):
    """
    Runs the Cassandra executable as a subprocess.

    Args:
        input_file_path: The path to the .inp file (e.g., "Examples/NVT/water_spc/nvt.inp").
        cassandra_exe_path: The path to the compiled cassandra executable.
    """
    print(f"--- Running Cassandra for {input_file_path} ---")

    # Get the FULL, ABSOLUTE path to the input file and executable
    inp_file_abs = os.path.abspath(input_file_path)
    exe_path_abs = os.path.abspath(cassandra_exe_path)

    # Get the directory where the input file lives.
    # We will run Cassandra from this directory so all output
    # files (log, prp, etc.) end up in the right place.
    run_directory = os.path.dirname(inp_file_abs)
    
    # This is the command Cassandra expects:
    # /path/to/executable /path/to/input.inp
    command = shlex.split(f"{exe_path_abs} {inp_file_abs}")

    print(f"Running in directory: {run_directory}")
    print(f"Running command: {command}")

    # Use subprocess.run to call the executable
    result = subprocess.run(
        command,
        cwd=run_directory,         # 'cwd' = run from this directory
        capture_output=True,     # Capture stdout/stderr
        text=True                  # As text, not bytes
    )

    # Print the output from Cassandra
    print("--- CASSANDRA STDOUT ---")
    print(result.stdout)
    
    if result.stderr:
        print("--- CASSANDRA STDERR ---")
        print(result.stderr)

    print("--- Cassandra run complete ---")
    return result