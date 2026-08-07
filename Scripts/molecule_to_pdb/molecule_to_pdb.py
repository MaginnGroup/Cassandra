#!/usr/bin/env python3
"""Create a 3D PDB file from a SMILES string or chemical name.

The program follows this sequence:

1. Obtain a SMILES representation directly from the user or from PubChem.
2. Ask RDKit to convert the SMILES text into a molecular graph.
3. Add explicit hydrogen atoms.
4. Generate several possible three-dimensional conformations.
5. Optimize each conformation and select the lowest-energy one.
6. Write the selected coordinates and connectivity to a PDB file.

This is a coordinate-generation tool, not a complete force-field assignment
program. A simulation package will still need atom types, partial charges, and
the appropriate bonded and nonbonded parameters.

Written by: Edward Maginn
Date: 2026-08-07
History:
2026-08-07: Initial version
"""

from __future__ import annotations

import argparse
import json
import re
import sys
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

# RDKit is the chemistry toolkit that performs structure parsing, 3D embedding,
# geometry optimization, and PDB output. Keeping the import in a try/except lets
# us give students a useful installation message instead of a traceback.
try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
except ImportError as exc:  # pragma: no cover - dependency guard
    raise SystemExit(
        "RDKit is required. Install dependencies with:\n"
        "  python -m pip install -r requirements.txt"
    ) from exc


# PubChem provides a web API called PUG REST. We use it only for chemical-name
# lookup; providing --smiles works without contacting PubChem.
PUBCHEM_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"


def smiles_from_name(name: str, timeout: float = 20.0) -> str:
    """Resolve a chemical name to a SMILES string using PubChem PUG REST.

    A URL cannot contain arbitrary spaces and punctuation, so ``quote`` converts
    the supplied name into a URL-safe representation. PubChem returns JSON, a
    structured text format that Python converts into dictionaries and lists.
    """
    encoded_name = urllib.parse.quote(name, safe="")
    url = (
        f"{PUBCHEM_BASE}/compound/name/{encoded_name}/property/"
        "ConnectivitySMILES,CanonicalSMILES/JSON"
    )
    # Supplying a User-Agent is good practice when accessing a public web API.
    request = urllib.request.Request(
        url,
        headers={"User-Agent": "molecule-to-pdb/1.0 (Python urllib)"},
    )
    # Network operations can fail for several reasons. A 404 specifically means
    # that PubChem did not recognize the requested name; other failures may be
    # server errors, network problems, or timeouts.
    try:
        with urllib.request.urlopen(request, timeout=timeout) as response:
            payload = json.load(response)
    except urllib.error.HTTPError as exc:
        if exc.code == 404:
            raise ValueError(f"PubChem found no compound named {name!r}.") from exc
        raise RuntimeError(f"PubChem lookup failed with HTTP {exc.code}.") from exc
    except urllib.error.URLError as exc:
        raise RuntimeError(f"Could not contact PubChem: {exc.reason}") from exc

    # The requested properties are stored in the first compound record in this
    # nested JSON structure. A defensive check gives a clearer error if PubChem
    # ever returns an incomplete or differently shaped response.
    try:
        properties = payload["PropertyTable"]["Properties"][0]
    except (KeyError, IndexError, TypeError) as exc:
        raise RuntimeError("PubChem returned an unexpected response.") from exc

    # Current PubChem responses normally use ConnectivitySMILES. Accepting the
    # older CanonicalSMILES key makes the program more tolerant of API changes.
    for key in ("ConnectivitySMILES", "CanonicalSMILES"):
        if properties.get(key):
            return str(properties[key])
    raise RuntimeError("PubChem response did not contain a usable SMILES string.")


def molecule_from_smiles(smiles: str) -> "Chem.Mol":
    """Convert SMILES text into an RDKit molecular graph with hydrogens.

    ``MolFromSmiles`` creates atoms and bonds but normally leaves most hydrogen
    atoms implicit. Explicit hydrogens are needed when we want them to appear as
    atoms with their own Cartesian coordinates in the PDB file.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES string: {smiles!r}")
    return Chem.AddHs(mol)


def generate_3d(
    mol: "Chem.Mol", *, seed: int = 20260806, conformers: int = 10
) -> tuple["Chem.Mol", str, float | None]:
    """Embed conformers, optimize them, and retain the lowest-energy result.

    A conformer is one possible three-dimensional arrangement of the same
    molecular graph. Flexible molecules can have many conformers, so sampling
    several candidates is generally better than generating only one.

    Returns the selected molecule, the force-field name, and its energy. The
    energy is useful for ranking conformers of this molecule; it should not be
    interpreted as an absolute thermodynamic energy.
    """

    # ETKDG combines experimental torsion preferences with distance geometry.
    # A fixed random seed makes classroom examples reproducible: running the
    # same command should produce the same initial conformers.
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    params.useRandomCoords = False
    conformer_ids = list(AllChem.EmbedMultipleConfs(mol, conformers, params))
    if not conformer_ids:
        raise RuntimeError("RDKit could not generate a 3D conformer.")

    # Each conformer is optimized to move it toward a nearby local minimum.
    # MMFF94 is preferred when parameters exist for every atom. UFF covers a
    # wider range of elements and serves as a fallback for many molecules.
    force_field = "none"
    energies: dict[int, float] = {}
    if AllChem.MMFFHasAllMoleculeParams(mol):
        force_field = "MMFF94"
        results = AllChem.MMFFOptimizeMoleculeConfs(mol, numThreads=0)
    else:
        try:
            force_field = "UFF"
            results = AllChem.UFFOptimizeMoleculeConfs(mol, numThreads=0)
        except (RuntimeError, ValueError):
            # The embedded coordinates are still usable if UFF is unavailable,
            # but there will be no optimized energy with which to rank them.
            force_field = "none"
            results = []

    # RDKit returns (convergence_status, energy) for each optimized conformer.
    # We retain the energy here; advanced versions of the exercise could also
    # warn when convergence_status is nonzero.
    for conf_id, result in zip(conformer_ids, results):
        _, energy = result
        energies[conf_id] = float(energy)

    # Missing energies are treated as infinity. If every energy is missing,
    # Python's min() keeps the first embedded conformer.
    best_id = min(conformer_ids, key=lambda cid: energies.get(cid, float("inf")))
    best_energy = energies.get(best_id)

    # The molecule currently holds all sampled conformers. Make a copy and keep
    # only the selected conformer so the PDB writer outputs one coordinate set.
    best = Chem.Mol(mol)
    selected = Chem.Conformer(mol.GetConformer(best_id))
    best.RemoveAllConformers()
    best.AddConformer(selected, assignId=True)
    return best, force_field, best_energy


def safe_stem(text: str) -> str:
    """Turn user-supplied text into a conservative default filename."""

    # Replace spaces and shell-sensitive punctuation, limit the length, and use
    # a generic name if nothing remains after cleaning.
    stem = re.sub(r"[^A-Za-z0-9_.-]+", "_", text.strip()).strip("._")
    return stem[:80] or "molecule"


def write_pdb(mol: "Chem.Mol", output: Path, title: str, smiles: str) -> None:
    """Write one RDKit molecule to a Cassandra-compatible PDB file.

    RDKit normally writes a compact CONECT section in which a bond may be listed
    in only one direction. Cassandra's ``mcfgen.py`` instead expects every atom
    to have a connectivity record. We therefore replace RDKit's CONECT records
    with a complete, symmetric adjacency list.
    """

    # parents=True allows an output such as structures/r125.pdb even when the
    # structures directory does not exist yet.
    output.parent.mkdir(parents=True, exist_ok=True)

    # RDKit places the _Name property in the PDB COMPND record. We also preserve
    # the exact SMILES used to create the structure in a human-readable remark.
    mol.SetProp("_Name", title)
    pdb_block = Chem.MolToPDBBlock(mol)

    # Retain RDKit's coordinate and metadata records, but remove its CONECT and
    # END records. A new connectivity section and one final END are added below.
    coordinate_lines = [
        line
        for line in pdb_block.splitlines()
        if not line.startswith("CONECT") and line.strip() != "END"
    ]

    # RDKit atom indices start at zero, whereas PDB atom serial numbers start at
    # one. GetNeighbors() reads the molecular graph created from the SMILES. A
    # line is written even for an atom with no neighbors because mcfgen.py uses
    # these records as a complete atom-by-atom adjacency table.
    connectivity_lines = []
    for atom in mol.GetAtoms():
        serial = atom.GetIdx() + 1
        neighbor_serials = sorted(
            neighbor.GetIdx() + 1 for neighbor in atom.GetNeighbors()
        )
        conect = f"CONECT{serial:5d}"
        conect += "".join(f"{neighbor:5d}" for neighbor in neighbor_serials)
        connectivity_lines.append(conect)

    remarks = f"REMARK  1 GENERATED FROM SMILES: {smiles}\n"
    complete_pdb = "\n".join(coordinate_lines + connectivity_lines + ["END", ""])
    output.write_text(remarks + complete_pdb, encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    """Define the command-line interface and return its argument parser."""

    parser = argparse.ArgumentParser(
        description="Generate an explicit-hydrogen 3D PDB from SMILES or a chemical name."
    )
    # A mutually exclusive group means the user must supply exactly one source:
    # either --smiles or --name, but never both in the same command.
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--smiles", help="SMILES string (quote it in the shell)")
    source.add_argument("--name", help="chemical name to resolve through PubChem")
    parser.add_argument("-o", "--output", type=Path, help="output PDB path")
    parser.add_argument(
        "--conformers", type=int, default=10, help="conformers to sample (default: 10)"
    )
    parser.add_argument("--seed", type=int, default=20260806, help="random seed")
    parser.add_argument("--timeout", type=float, default=20.0, help="name lookup timeout")
    return parser


def main() -> int:
    """Run the complete conversion workflow and return a shell exit status."""

    args = build_parser().parse_args()
    if args.conformers < 1:
        raise SystemExit("--conformers must be at least 1")

    try:
        # Name input must first be translated to SMILES. Direct SMILES input can
        # proceed immediately and does not require an internet connection.
        if args.name:
            smiles = smiles_from_name(args.name, timeout=args.timeout)
            title = args.name
            print(f"Resolved {args.name!r} to SMILES: {smiles}")
        else:
            smiles = args.smiles
            title = "SMILES molecule"

        # At this point both input routes have produced the same representation,
        # so the remainder of the workflow is identical for names and SMILES.
        mol = molecule_from_smiles(smiles)
        mol, force_field, energy = generate_3d(
            mol, seed=args.seed, conformers=args.conformers
        )
        # Use the requested path when present. Otherwise, a chemical name gives
        # a descriptive filename, while direct SMILES input uses molecule.pdb.
        default_name = safe_stem(args.name or "molecule") + ".pdb"
        output = args.output or Path(default_name)
        write_pdb(mol, output, title, smiles)
    # Expected user, chemistry, and network errors are printed concisely. A zero
    # return value means success to the shell; one indicates failure.
    except (ValueError, RuntimeError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1

    energy_text = "unavailable" if energy is None else f"{energy:.6f}"
    print(f"Wrote {output.resolve()}")
    print(f"Atoms: {mol.GetNumAtoms()} (hydrogens included)")
    print(f"Geometry optimization: {force_field}; energy: {energy_text}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
