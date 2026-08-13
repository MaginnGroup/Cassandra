# OPLS dihedrals and Ryckaert–Bellemans (RB) conversion

How Cassandra treats **OPLS** proper dihedrals when reading an MCF, and why
energy evaluation usually sees an **RB torsion** instead.

**Related:** [MCF setup workflow](mcf-setup-workflow.md) ·
`Get_Dihedral_Info` in `Src/input_routines.f90` ·
`Compute_Molecule_Dihedral_Energy` in `Src/energy_routines.f90`

---

## What you write in the MCF

Under `# Dihedral_Info`, an OPLS line looks like:

```text
# Dihedral_Info
N
i  a1 a2 a3 a4  OPLS  a0  a1  a2  a3
```

with `a0`…`a3` in **kJ/mol**. The analytic form is

\[
E_{\mathrm{OPLS}}(\phi)
  = a_0
  + a_1\bigl(1+\cos\phi\bigr)
  + a_2\bigl(1-\cos 2\phi\bigr)
  + a_3\bigl(1+\cos 3\phi\bigr).
\]

You can also enter type `RB` / `Ryckaert-Bellemans` with up to six \(C_k\)
coefficients (kJ/mol) directly.

---

## What Cassandra does at read time

In `Get_Dihedral_Info`:

1. Convert \(a_0\ldots a_3\) to internal (atomic) energy units.
2. If all coefficients are ~0, mark the dihedral `none`.
3. Otherwise mark it **RB-formatted** and map to

\[
E_{\mathrm{RB}}(\phi)=\sum_{k=0}^{5} C_k\,\cos^k\phi
\]

with \(C_4=C_5=0\) for a pure OPLS term and

| Coefficient | From OPLS |
|-------------|-----------|
| \(C_0\) | \(a_2 + (a_0+a_1+a_2+a_3)\) |
| \(C_1\) | \(a_1 - 3 a_3\) |
| \(C_2\) | \(-2 a_2\) |
| \(C_3\) | \(4 a_3\) |

4. Store the torsion as type `"RB torsion"` (`int_rb_torsion`).
5. **Combine** with any existing RB entry on the same four atoms (see below).

Some CHARMM dihedrals with integer \(n\) and \(\delta \in \{0^\circ,180^\circ\}\)
are also rewritten as RB polynomials and can combine the same way.

---

## Combining multiple terms on the same atoms

Force fields (and MCFs) often list **several** contributions that belong to
the **same** torsion — for example two OPLS lines, an OPLS line plus an RB
line, or several CHARMM multiplicities that Cassandra could map to RB.

Because

\[
E_{\mathrm{RB}}=\sum_k C_k\cos^k\phi
\]

is linear in the \(C_k\), adding two potentials is the same as adding their
coefficient vectors.

When a new RB-formatted term is read, `Get_Dihedral_Info`:

1. Builds its \(C_k\) vector.
2. Looks through the RB list already built for this species.
3. Declares a **match** if the four atom indices are identical **or**
   reversed \((a_1,a_2,a_3,a_4)\) vs \((a_4,a_3,a_2,a_1)\).
   Reversed order is allowed because \(\cos(-\phi)=\cos\phi\), so the RB
   energy does not depend on the sign of \(\phi\).
4. On a match: `rb_c = rb_c + C_new` (element-wise) and no new list slot.
5. If none match: append a new `"RB torsion"` entry.

**Not combined this way:** harmonic dihedrals, and CHARMM terms that could
not be rewritten as RB (non-integer \(n\), or \(\delta\) not \(0^\circ\)/\(180^\circ\)).
Those remain separate entries evaluated in the second pass of
`Compute_Molecule_Dihedral_Energy`.

**Example:** an MCF with two `OPLS` lines on atoms 2–3–4–5 (or 5–4–3–2)
becomes **one** RB torsion whose \(C_k\) are the sum of the two mapped
coefficient sets. The energy loop therefore evaluates one polynomial for
that unique atom quartet.

---

## What the energy routine sees

`Compute_Molecule_Dihedral_Energy` first loops `1 … ndihedrals_rb` and
evaluates the RB polynomial via `DOT_PRODUCT(cos^k(phi), rb_c)`.

The `CASE(int_opls)` branch in that routine is a **legacy fallback**. After a
normal MCF read, OPLS terms have already been converted, so that case is
almost never hit. Keep writing `OPLS` in the MCF when that is your force-field
source — the conversion is automatic and intentional.

---

## Why convert?

- One polynomial evaluation per unique torsion (after combining).
- Same code path for native RB and OPLS-derived terms.
- \(\cos(-\phi)=\cos\phi\), so reversed atom order does not change the RB
  energy (handy when combining duplicates).
