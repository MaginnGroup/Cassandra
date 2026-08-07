# Molecule to PDB

This command-line tool creates an explicit-hydrogen, three-dimensional PDB file
from either a SMILES string or a chemical name. RDKit generates and optimizes
the conformers. Chemical-name lookup uses the PubChem PUG REST service and thus
requires an internet connection. The output contains a complete, symmetric
`CONECT` section suitable for Cassandra's `mcfgen.py`: every atom has its own
connectivity record, and each bond appears in both directions. 

## Where this tool fits in the setup workflow

You need one pdb file for each species you wish to simulate. You can obtain these in different ways. This tool will generate one for you using a molecule name or SMILES string. It relies on RDKit to do this. 

SMILES/name
  → molecule_to_pdb.py → PDB (+ CONECT)
  → hand-edit atom types onto each HETATM/ATOM line
  → mcfgen.py --ffTemplate → .ff (fill by hand; Cassandra units)
  → mcfgen.py → .mcf
  → library_setup.py → fragment libraries
  → Cassandra .inp

Run these tools from the directory containing the molecule files 
(same convention as mcfgen.py).

Full returning-user checklist:
[docs/mcf-setup-workflow.md](../../docs/mcf-setup-workflow.md).

## Short hand-edit required after running this script (and after collecting any generic pdb)

*NOTE:* you must manually edit the pdb after it is created by adding force field atom types to the end of each atom line. PDB files only know about elements not atom types. See Examples/NVT/diethylether_mie/dee.pdb to see what this looks like. 

# After molecule_to_pdb (element only):
HETATM    1  C  ...                      C

# After hand-edit for mcfgen (atom type "CT" at end of line):
HETATM    1  C  ...                      C  CT

## Setup

This project uses the existing `cassandra-dev` Conda environment. Activate the Conda environment and install RDKit from `conda-forge`:

```bash
conda activate cassandra-dev
conda install -c conda-forge rdkit
```

RDKit is then part of `cassandra-dev` and is available whenever that environment is active. 
Prefer conda-forge into cassandra-dev. requirements.txt is only a pip fallback if you are not using Conda.

Confirm that the shell is using the intended Python interpreter and that RDKit imports successfully:

```bash
echo "$CONDA_DEFAULT_ENV"
which python
python --version
python -c "import rdkit; print('RDKit:', rdkit.__version__)"
```

The environment name should be `cassandra-dev`, and `which python` should
normally report:

```text
/opt/anaconda3/envs/cassandra-dev/bin/python
```

If both `(cassandra-dev)` and `(.venv)` appear in the shell prompt, deactivate
the `.venv` environment before continuing:

```bash
deactivate
conda activate cassandra-dev
```

## Examples

From SMILES:

```bash
python molecule_to_pdb.py --smiles 'CCO' --output ethanol.pdb
```

From a chemical name:

```bash
python molecule_to_pdb.py --name '1,1,1,2,2-pentafluoroethane' --output r125.pdb
```

Sample additional conformers for a flexible molecule:

```bash
python molecule_to_pdb.py --name ibuprofen --conformers 50 --output ibuprofen.pdb
```

The name and SMILES options are deliberately separate. This avoids silently
interpreting a short chemical name as a valid SMILES string.

## Molecular-simulation warning

The PDB contains atom identities, connectivity records, explicit hydrogens, and
one optimized geometry. It does **not** contain force-field atom types, partial
charges, bonded parameters, or validated protonation/stereochemistry choices.
Check those separately before using the structure in a molecular simulation. 
