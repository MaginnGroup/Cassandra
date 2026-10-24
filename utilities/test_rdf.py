#!/usr/bin/env python3
#*******************************************************************************
#  test_rdf.py
#
#  Synthetic boxes for the RDF: a known pair distance, a molecule that
#  straddles the periodic boundary, and the ideal-gas coordination integral.
#  The 243 K movie is only a smoke test that the reader finishes.
#
#  Author:  Edward J. Maginn (EJM), University of Notre Dame
#  Written: 10/24/26
#
#*******************************************************************************
"""Check RDF peaks, unwrapped centers of mass, and the coordination integral."""

from __future__ import annotations

import sys
import tempfile
from pathlib import Path

import numpy as np

_REPO_ROOT = Path(__file__).resolve().parent.parent
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from utilities.rdf import (
    _frame_g,
    _pair_histogram,
    bin_edges,
    compute_rdf,
    first_minimum_after_peak,
    running_coordination,
    cutoff_label,
    distance_axis_label,
    distance_for_axis,
    plot_title,
    resolve_legend,
    selection_tag,
    write_plot,
)


NVT243 = Path("/Users/edwardmaginn/CassandraV2/projects/hfc125/NVT-243")
NVT243_MCF = Path(
    "/Users/edwardmaginn/CassandraV2/projects/hfc125/NPT-Runs/NPT-243/r125.mcf"
)

ONE_ATOM_MCF = """# Atom_Info
1
1 C C 12.0 0.0

# Bond_Info
0

# Angle_Info
0
"""

DIATOMIC_MCF = """# Atom_Info
2
1 C C 12.0 0.0
2 C C 12.0 0.0

# Bond_Info
1
1 1 2 fixed 0.4

# Angle_Info
0
"""


def _write_h(path: Path, length: float, nmol: int) -> None:
    volume = length ** 3
    path.write_text(
        "\n".join(
            [
                f"   {volume:.8f}",
                f"   {length:.8f}  0.0  0.0",
                f"   0.0  {length:.8f}  0.0",
                f"   0.0  0.0  {length:.8f}",
                "",
                "   1",
                f"   1   {nmol}",
                "",
            ]
        ),
        encoding="utf-8",
    )


def _write_xyz(path: Path, rows: list[tuple[str, float, float, float]]) -> None:
    lines = [f"{len(rows):8d}", " MC_STEP: 1"]
    for element, x, y, z in rows:
        lines.append(f"{element:<2s} {x:18.8f} {y:18.8f} {z:18.8f}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _peak_bin(result) -> int:
    return int(np.argmax(result.g))


def main() -> int:
    edges = bin_edges(5.0, 0.05)
    g_ideal = np.ones(len(edges) - 1)
    rho = 0.02
    integrated = running_coordination(g_ideal, edges, rho)[-1]
    analytic = (4.0 / 3.0) * np.pi * (5.0 ** 3) * rho
    if abs(integrated - analytic) > 1.0e-9:
        raise SystemExit(
            f"ideal-gas coordination {integrated} is not {analytic}"
        )

    centers = np.array([0.25, 0.75, 1.25, 1.75, 2.25, 2.75])
    curve = np.array([0.2, 2.0, 3.0, 1.4, 0.7, 1.1])
    minimum = first_minimum_after_peak(centers, curve)
    if minimum is None or abs(minimum - 2.25) > 1.0e-12:
        raise SystemExit(f"first minimum was {minimum}, expected 2.25")
    shoulder_centers = np.arange(0.05, 1.15, 0.1)
    shoulder_curve = np.array([0.2, 1.2, 2.5, 1.7, 1.4, 1.42, 1.3, 0.9, 0.7, 0.75, 1.1])
    shoulder_minimum = first_minimum_after_peak(shoulder_centers, shoulder_curve)
    if shoulder_minimum is None or abs(shoulder_minimum - 0.85) > 1.0e-12:
        raise SystemExit(
            "first-shell minimum should be the lowest bin between the "
            f"second and third g(r)=1 crossings, got {shoulder_minimum}"
        )

    cell = np.diag([10.0, 10.0, 10.0])
    short_edges = np.array([0.0, 1.0])
    a_points = np.array(
        [[0.0, 0.0, 0.0], [0.2, 0.0, 0.0], [9.8, 0.0, 0.0], [9.6, 0.0, 0.0]]
    )
    a_molecules = np.array([0, 0, 1, 1])
    b_points = np.array([[0.1, 0.0, 0.0], [9.7, 0.0, 0.0]])
    b_molecules = np.array([0, 1])
    aa_inter = _pair_histogram(
        a_points, a_molecules, a_points, a_molecules, cell, short_edges,
        like=True, exclude_intramolecular=True,
    )
    aa_all = _pair_histogram(
        a_points, a_molecules, a_points, a_molecules, cell, short_edges,
        like=True, exclude_intramolecular=False,
    )
    ab_inter = _pair_histogram(
        a_points, a_molecules, b_points, b_molecules, cell, short_edges,
        like=False, exclude_intramolecular=True,
    )
    ab_all = _pair_histogram(
        a_points, a_molecules, b_points, b_molecules, cell, short_edges,
        like=False, exclude_intramolecular=False,
    )
    if aa_inter[0] != 4 or aa_all[0] != 6 or ab_inter[0] != 4 or ab_all[0] != 8:
        raise SystemExit("site-pair histogram did not count PBC or intramolecular pairs correctly")
    shells = np.array([1.0])
    g_aa, rho_aa = _frame_g(
        aa_inter, shells, 10.0, nmol=2, sites_i=2, sites_j=2,
        like=True, exclude_intramolecular=True,
    )
    g_ab, rho_ab = _frame_g(
        ab_inter, shells, 10.0, nmol=2, sites_i=2, sites_j=1,
        like=False, exclude_intramolecular=True,
    )
    if g_aa[0] != 10.0 or rho_aa != 0.2 or g_ab[0] != 10.0 or rho_ab != 0.1:
        raise SystemExit("intermolecular site-pair normalization is not exact")
    g_ab_all, rho_ab_all = _frame_g(
        ab_all, shells, 10.0, nmol=2, sites_i=2, sites_j=1,
        like=False, exclude_intramolecular=False,
    )
    if g_ab_all[0] != 10.0 or rho_ab_all != 0.2:
        raise SystemExit("all-site normalization is not exact")

    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        (root / "one.mcf").write_text(ONE_ATOM_MCF, encoding="utf-8")
        _write_xyz(
            root / "pair.xyz",
            [("C", 0.0, 0.0, 0.0), ("C", 3.2, 0.0, 0.0)],
        )
        _write_h(root / "pair.H", 10.0, 2)
        pair = compute_rdf(
            root / "pair.xyz",
            root / "one.mcf",
            root / "pair.H",
            mode="com",
            bin_width=0.5,
            r_max=5.0,
            r_cut=4.0,
        )
        peak = _peak_bin(pair)
        if not (3.0 <= pair.r_center[peak] < 3.5):
            raise SystemExit(f"pair peak is at {pair.r_center[peak]} A, expected 3.0–3.5")
        if np.count_nonzero(pair.g > 0.0) != 1:
            raise SystemExit("the two-particle RDF has more than one occupied bin")
        if selection_tag("com", element_i=None, element_j=None, atom_type_i=None, atom_type_j=None, atom_index_i=None, atom_index_j=None) != "com":
            raise SystemExit("COM file tag should be com")
        if selection_tag("atom", element_i="F", element_j=None, atom_type_i=None, atom_type_j=None, atom_index_i=None, atom_index_j=None) != "F-F":
            raise SystemExit("like-atom file tag should repeat the selection")
        if plot_title("COM-COM RDF, species 1", "cn") != "COM-COM CN, species 1":
            raise SystemExit("CN title should name the coordination number")
        if plot_title("COM-COM RDF, species 1", "rdf", "HFC-125 COM") != "HFC-125 COM":
            raise SystemExit("a requested title should replace the default")
        atom_note = "Atom-atom RDF, species 1: atom type 'c3_1_s1' - atom type 'c3_1_s1' (same selection)"
        if plot_title(atom_note, "rdf") != atom_note[: -len(" (same selection)")]:
            raise SystemExit("the default atom title should omit the selection note")
        if resolve_legend("rdf", None) != "none" or resolve_legend("cn", None) != "none":
            raise SystemExit("RDF and CN figures should have no legend by default")
        if resolve_legend("both", None) != "with-cutoff":
            raise SystemExit("the combined figure should name the curves and the cutoff")
        if float(distance_for_axis(10.0, "nm")) != 1.0 or distance_axis_label("A") != "r (Å)":
            raise SystemExit("the plot axis should default to nm and draw angstroms as Å")
        if cutoff_label(6.0, "nm") != "r_cut 0.6 nm" or cutoff_label(6.0, "A") != "r_cut 6 Å":
            raise SystemExit("the cutoff legend should use the plot unit")
        plot_path = root / "pair_rdf_cn.png"
        write_plot(pair, plot_path, "both")
        write_plot(pair, root / "pair_cn.png", "cn")
        if not plot_path.is_file() or plot_path.stat().st_size < 1000:
            raise SystemExit("combined RDF plot was not written")

        (root / "dimer.mcf").write_text(DIATOMIC_MCF, encoding="utf-8")
        _write_xyz(
            root / "split.xyz",
            [
                ("C", 0.2, 0.0, 0.0),
                ("C", 9.8, 0.0, 0.0),
                ("C", 3.0, 0.0, 0.0),
                ("C", 3.4, 0.0, 0.0),
            ],
        )
        _write_h(root / "split.H", 10.0, 2)
        split = compute_rdf(
            root / "split.xyz",
            root / "dimer.mcf",
            root / "split.H",
            mode="com",
            bin_width=0.2,
            r_max=6.0,
            r_cut=6.0,
        )
        peak = _peak_bin(split)
        if abs(split.r_center[peak] - 3.2) > 0.15:
            raise SystemExit(
                f"PBC-split COM peak is at {split.r_center[peak]} A, expected 3.2"
            )
        false_com = np.flatnonzero(np.abs(split.r_center - 1.8) < 0.15)
        if false_com.size and np.any(split.g[false_com] > 0.0):
            raise SystemExit("wrapped coordinates produced the middle-of-box COM distance")

        atoms = compute_rdf(
            root / "split.xyz",
            root / "dimer.mcf",
            root / "split.H",
            mode="atom",
            element_i="C",
            bin_width=0.2,
            r_max=6.0,
            r_cut=6.0,
        )
        if "excluded" not in atoms.exclusion_note:
            raise SystemExit("atom–atom RDF did not report the intramolecular exclusion")

    if NVT243.is_dir() and NVT243_MCF.is_file():
        smoke = compute_rdf(
            NVT243 / "R125-T243-NVT.xyz",
            NVT243_MCF,
            NVT243 / "R125-T243-NVT.H",
            mode="com",
            bin_width=0.2,
            r_cut=6.0,
        )
        if smoke.n_frames != 10 or not np.all(np.isfinite(smoke.g)):
            raise SystemExit("243 K COM smoke run did not produce a finite 10-frame table")
        if not np.isfinite(smoke.coordination):
            raise SystemExit("243 K coordination number is not finite")
        fluorine = compute_rdf(
            NVT243 / "R125-T243-NVT.xyz",
            NVT243_MCF,
            NVT243 / "R125-T243-NVT.H",
            mode="atom",
            element_i="F",
            bin_width=0.2,
            r_cut=6.0,
        )
        if fluorine.n_frames != 10 or not np.all(np.isfinite(fluorine.g)):
            raise SystemExit("243 K fluorine RDF did not produce a finite table")
        text = fluorine.table()
        if "mean_shell_count" not in text or "coordination number" not in text:
            raise SystemExit("RDF table is missing the coordination columns")
    else:
        print("NVT-243 smoke test skipped; the project trajectory is not on this machine")

    print("rdf tests passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
