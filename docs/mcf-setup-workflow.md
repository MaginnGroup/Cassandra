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

OPLS dihedrals are converted to Ryckaert–Bellemans form when the MCF is read —
see [dihedral-opls-rb.md](dihedral-opls-rb.md).

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

## Important: PDB geometry ≠ MCF geometry

This is easy to miss and worth understanding before you trust a run.

### What the PDB is for

`molecule_to_pdb` (PubChem / RDKit from SMILES) produces a **reasonable 3-D
sketch**: connectivity (`CONECT`) plus Cartesian coordinates that look like
the molecule. Those coordinates are **not** required to match your force-field
bond lengths or equilibrium angles.

What Cassandra needs from the PDB for setup is mainly:

1. **Atom order** (must match the MCF `# Atom_Info` order after `mcfgen`)
2. **Connectivity** (`CONECT` → bonds in the MCF)
3. A starting guess for fragment generation

The **MCF** is the geometry authority: fixed bond lengths \(r_0\), harmonic
angle \(\theta_0\) / \(K_\theta\), dihedrals, charges, LJ parameters.

### Worked example — R125 (`r125.pdb` vs `r125.mcf`)

For 1,1,1,2,2-pentafluoroethane (project `hfc125`), the SMILES-generated PDB
and the literature-parameter MCF disagree. Atom order matches (C1, C2, F1–F3,
F4, F5, H1 → MCF atoms 1–8).

**Bonds (Å)** — PDB measured vs MCF fixed \(r_0\):

| # | Atoms | \(r\) (PDB) | \(r_0\) (MCF) | \(\Delta\) |
|---|-------|-------------|---------------|------------|
| 1 | C–C (1–2) | 1.521 | 1.538 | −0.017 |
| 2 | C–F (1–6) | 1.361 | 1.350 | +0.011 |
| 3 | C–F (1–7) | 1.360 | 1.350 | +0.010 |
| 4 | C–H (1–8) | 1.093 | 1.096 | −0.003 |
| 5 | C–F (2–3) | 1.357 | 1.350 | +0.007 |
| 6 | C–F (2–4) | 1.356 | 1.350 | +0.006 |
| 7 | C–F (2–5) | 1.352 | 1.350 | +0.002 |

**Angles (°)** — PDB measured vs MCF \(\theta_0\) (all harmonic in this MCF):

| # | Atoms | \(\theta\) (PDB) | \(\theta_0\) (MCF) | \(\Delta\) |
|---|-------|------------------|--------------------|------------|
| 1 | C–C–F (2–1–6) | 111.40 | 109.24 | +2.16 |
| 2 | C–C–F (2–1–7) | 111.38 | 109.24 | +2.14 |
| 3 | C–C–H (2–1–8) | 111.15 | 110.32 | +0.83 |
| 4 | F–C–F (6–1–7) | 106.30 | 107.36 | −1.06 |
| 5 | F–C–H (6–1–8) | 108.22 | 108.79 | −0.57 |
| 6 | F–C–H (7–1–8) | 108.21 | 108.79 | −0.58 |
| 7 | C–C–F (1–2–3) | 111.18 | 109.24 | +1.94 |
| 8 | C–C–F (1–2–4) | 111.17 | 109.24 | +1.93 |
| 9 | C–C–F (1–2–5) | 112.65 | 109.24 | +3.41 |
| 10 | F–C–F (3–2–4) | 106.63 | 107.36 | −0.73 |
| 11 | F–C–F (3–2–5) | 107.46 | 107.36 | +0.10 |
| 12 | F–C–F (4–2–5) | 107.46 | 107.36 | +0.10 |

So the starter PDB can be off by ~0.02 Å in bonds and a few degrees in angles.
That is **normal** for a SMILES/PubChem geometry.

### How Cassandra “corrects” this

Production configurations are **not** copies of the PDB Cartesians.

1. **`library_setup.py`** builds fragment libraries using the MCF (and
   Cassandra) so fragment geometries respect **fixed bond lengths** in the
   MCF.
2. **`# Start_Type make_config`** assembles the initial box from those
   fragments / MCF topology. Molecules in the run therefore have bonds at
   the MCF \(r_0\) values (constraints), while harmonic angles fluctuate
   around \(\theta_0\) under Monte Carlo.

You can verify after a short run with the diagnostics:

```bash
# Bonds should sit on r0 (near machine precision aside from XYZ output rounding)
python ~/CassandraV2/Cassandra/diagnostics/bond_distribution.py \
  run.out.xyz species.mcf --bond 1 --no-show

# Angles: sampled PDF vs Boltzmann at T (mean near theta0)
python ~/CassandraV2/Cassandra/diagnostics/angle_distribution.py \
  run.out.xyz species.mcf --inp run.inp --angle 1 --no-show
```

For the R125 equilibration movie, all seven fixed bonds matched MCF \(r_0\)
to ~\(10^{-8}\)–\(10^{-13}\) Å even though the original `r125.pdb` did not —
exactly the behavior above.

**Teaching takeaway:** do not expect the PDB to equal the MCF geometry; expect
the **simulation** (after fragment setup + `make_config`) to equal the MCF
constraints. If `bond_distribution` shows large \(|r - r_0|\`, something is
wrong (wrong MCF, atom-order mismatch, or species layout) — not “the PDB was
approximate.”

### Packmol / `read_config` — do wrong bonds “heal”?

Some workflows build the initial box with **Packmol** (or another packer) and
start Cassandra with `# Start_Type read_config`. That can work, but
intramolecular geometry is **not** automatically repaired the way
`make_config` is.

Abbreviated rules (fuller discussion planned for the user-guide / Read the
Docs):

- Cassandra treats MCF bonds as **fixed**. On energy evaluation it **checks**
  \(|r - r_0|\) against a tolerance (default **0.01 Å**) and **aborts** if a
  bond is broken — it does **not** rescale coordinates to \(r_0\).
- Translate, rotate, angle, dihedral, and volume moves **preserve** existing
  bond lengths (rigid-body / COM shifts). Harmonic **angles** and **dihedrals**
  can relax under their MC moves; fixed **bonds** do not.
- Fragment **regrowth** / insert / `make_config` rebuild molecules from the
  fragment library (MCF \(r_0\)). That fixes bonds only for molecules that
  actually get rebuilt — not every molecule in the box by default.
- So a Packmol box with SMILES-like monomers (bonds off by ~0.01–0.02 Å, as in
  the R125 PDB above) may **fail at startup** or, if under tolerance, **keep**
  slightly wrong bonds for the whole run unless those molecules are regrown.

**Practical advice:** pack **MCF-correct** monomers (or prefer
`make_config` + fragment libraries). Do not assume “regrowth will fix
everything over time.”

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
- [ ] Understand PDB ≠ MCF geometry; MCF + fragments/`make_config` set bonds
- [ ] Fragment libraries generated
- [ ] `# Rcutoff_Low` safe for this chemistry
- [ ] Optional: `bond_distribution` / `angle_distribution` on a short run
- [ ] Do not `git add` simulation outputs blindly

---

## Related docs

- [docs/README.md](README.md) — V2 docs index
- [dihedral-opls-rb.md](dihedral-opls-rb.md) — OPLS → RB conversion at MCF read
- [smoke-test.md](smoke-test.md) — running the Python smoke test
- [diagnostics/README.md](../diagnostics/README.md) — bond / angle / dihedral checks
- [V2_DEVELOPMENT_ENVIRONMENT.md](V2_DEVELOPMENT_ENVIRONMENT.md) — conda / compile
