<!-- end list -->

```python
# In python/api.py
import subprocess
import os
import shlex

def run_cassandra(input_file_path: str, cassandra_exe_path: str = "cassandra.exe"):
    """
    Runs the Cassandra executable as a subprocess.

    Args:
        input_file_path: The path to the .inp file (e.g., "Examples/run.inp").
        cassandra_exe_path: The path to the compiled cassandra executable.
    """
    print(f"--- Running Cassandra for {input_file_path} ---")

    # Get the directory where the input file lives
    run_directory = os.path.dirname(input_file_path)
    
    # Get just the name of the input file
    inp_file = os.path.basename(input_file_path)

    # We must run Cassandra in the same directory as the input file
    # so it can find all the other files (like .mol, .dat)
    
    # Note: shlex.split is safer than a simple list for commands
    command = shlex.split(f"{cassandra_exe_path} -i {inp_file}")

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
```