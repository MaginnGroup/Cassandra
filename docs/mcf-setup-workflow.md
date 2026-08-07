# MCF Setup Workflow (V2 checklist)

Returning-user checklist for building a new molecule from SMILES/name through
a runnable Cassandra input. This is **not** a full manual — it points to the
detailed Script READMEs.

**Related:** [docs index](README.md) ·
[molecule_to_pdb](../Scripts/molecule_to_pdb/README.md) ·
[mcfgen](../Scripts/MCF_Generation/README) ·
[library_setup](../Scripts/Frag_Library_Setup/README) ·
[smoke test](smoke-test.md)

---

## Pipeline overview

```text
SMILES/name
  → molecule_to_pdb.py → PDB (+ CONECT)
  → hand-edit atom types onto each HETATM/ATOM line
  → mcfgen.py --ffTemplate → .ff (fill by hand; Cassandra units)
  → mcfgen.py → .mcf
  → library_setup.py → fragment libraries
  → Cassandra .inp
```

Work in the directory that will hold the molecule files (PDB, FF, MCF, inp).
Activate `cassandra-dev` first.

```bash
conda activate cassandra-dev
cd /path/to/your/molecule/project   # e.g. ~/CassandraV2/projects/hfc125
```

---

## Step 1 — Generate a PDB

```bash
python ~/CassandraV2/Cassandra/Scripts/molecule_to_pdb/molecule_to_pdb.py \
  --name '1,1,1,2,2-pentafluoroethane' --output r125.pdb
```

Or from SMILES:

```bash
python ~/CassandraV2/Cassandra/Scripts/molecule_to_pdb/molecule_to_pdb.py \
  --smiles 'CCO' --output ethanol.pdb
```

`--name` needs internet (PubChem). Details:
[Scripts/molecule_to_pdb/README.md](../Scripts/molecule_to_pdb/README.md).

**Check:** PDB has explicit hydrogens and a full `CONECT` section.

---

## Step 2 — Add atom types and inspect CONECT

`mcfgen.py` needs a force-field **atom type** at the end of each atom line.
`molecule_to_pdb` only writes the element.

```text
# After molecule_to_pdb:
HETATM    1  C  ...                      C

# After hand-edit (example type CT):
HETATM    1  C  ...                      C  CT
```

Reference style: `Examples/NVT/diethylether_mie/dee.pdb`.

Also inspect `CONECT`: every bond should appear in both directions; highly
branched atoms may need merged CONECT lines (see MCF_Generation README note on
PF6-style connectivity).

**Check:** every atom line has a type; CONECT looks complete and symmetric.

---

## Step 3 — Build FF template, fill parameters, write MCF

From the molecule directory (PDB present; same basename for FF):

```bash
python ~/CassandraV2/Cassandra/Scripts/MCF_Generation/mcfgen.py r125.pdb --ffTemplate
```

Edit `r125.ff` by hand with literature parameters. **Cassandra units:**

| Parameter | Units |
|-----------|--------|
| LJ ε | K |
| LJ σ | Å |
| Bond length (fixed) | Å |
| Angle θ₀ | degrees |
| Angle \(K_\theta\) | **K/rad²** (energy/\(k_B\); no ½ in \(E = K(\theta-\theta_0)^2\)) |
| OPLS / CHARMM dihedral coeffs | kJ/mol |
| Harmonic dihedral / improper \(K\) | K/rad² |

If a published angle \(K\) is in kJ/mol/rad²:

\[
K_{\text{Cassandra}} = K_{\text{published}} / R
\quad (R \approx 8.314 \times 10^{-3}\ \mathrm{kJ\,mol^{-1}\,K^{-1}})
\]

Confirm whether the published form includes \(1/2\); Cassandra does **not**.

Then generate the MCF:

```bash
python ~/CassandraV2/Cassandra/Scripts/MCF_Generation/mcfgen.py r125.pdb
```

Details and caveats: [Scripts/MCF_Generation/README](../Scripts/MCF_Generation/README).

**Check:** `r125.mcf` exists; fragment section looks reasonable; charges sum ~0.

---

## Step 4 — Fragment library setup

Needs: MCF, a draft `.inp`, and the PDB(s) in the same directory.

```bash
python ~/CassandraV2/Cassandra/Scripts/Frag_Library_Setup/library_setup.py \
  ~/CassandraV2/Cassandra/Src/cassandra_gfortran.exe \
  r125-equil.inp \
  r125.pdb
```

This creates `species1/frag*/frag*.dat` (and related inputs) and updates the
`# Fragment_Files` section of the inp.

Details: [Scripts/Frag_Library_Setup/README](../Scripts/Frag_Library_Setup/README).

**Check:** fragment `.dat` files exist; `# Fragment_Files` in the inp points to them.

---

## Step 5 — Finish the `.inp` and run

Use an existing example inp as a template (NPT/NVT/GCMC as needed). Key items:

- `# Molecule_Files` — MCF and max molecules
- `# Start_Type` — e.g. `make_config N`
- `# Fragment_Files` — from library_setup
- `# Rcutoff_Low` — hardcore overlap cutoff (Å)

**Fluorinated molecules (e.g. R125):** geminal F···H 1–3 distances are often
~1.9 Å. `# Rcutoff_Low 2.0` can abort after `make_config` with
`Atomic overlap in the configuration` even when the box is fine. Use **1.0**
or **0.85** unless you know all nonbonded 1–3 distances are longer.

Then run:

```bash
~/CassandraV2/Cassandra/Src/cassandra_gfortran.exe r125-equil.inp
```

Or via the Python wrapper once paths are set (see [python/README.md](../python/README.md)).

**Check:** log reaches `Cassandra simulation complete`; `.prp` / `.xyz` written.

---

## Quick “did I forget?” list

- [ ] `cassandra-dev` active; RDKit installed if using molecule_to_pdb
- [ ] Atom types on every PDB atom line
- [ ] FF filled in Cassandra units (angle \(K\) in K/rad²)
- [ ] Fragment libraries generated
- [ ] `# Rcutoff_Low` safe for this chemistry
- [ ] Do not `git add` simulation outputs blindly

---

## Related docs

- [docs/README.md](README.md) — V2 docs index
- [smoke-test.md](smoke-test.md) — running the Python smoke test
- [V2_DEVELOPMENT_ENVIRONMENT.md](V2_DEVELOPMENT_ENVIRONMENT.md) — conda / compile
