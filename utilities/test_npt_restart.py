#!/usr/bin/env python3
#*******************************************************************************
#  test_npt_restart.py
#
#  Molecular centers of mass are mapped from L = 10 A to L = 8 A without
#  scaling bonds, angles, or dihedrals.  The generated NVT input must drop
#  the volume move. The GCMC input must keep a placeholder chemical potential.
#
#  Author:  Edward J. Maginn (EJM), University of Notre Dame
#  Written: 10/16/26
#*******************************************************************************
"""Check NPT-to-NVT scaling and the GCMC template."""

from __future__ import annotations

import json
import sys
import tempfile
from pathlib import Path

import numpy as np

_REPO_ROOT = Path(__file__).resolve().parent.parent
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from diagnostics.io.hfile import read_h
from diagnostics.io.mcf import read_mcf
from diagnostics.io.xyz import read_xyz
from utilities.npt_restart import _mic_delta, _unwrapped_molecule, build_restart


INP = """! leftover comment that Cassandra ignores
# Run_Name
demo-npt

# Sim_Type
npt_mc

# Nbr_Species
1

# Molecule_Files
demo.mcf 1

# Box_Info
1
cubic
10.0

# Temperature_Info
300.0

# Pressure_Info
2.0

# Move_Probability_Info

# Prob_Translation
0.40
0.30

# Prob_Rotation
0.40
30.0

# Prob_Regrowth
0.195
1.0

# Prob_Volume
0.005
800

# Done_Probability_Info

# Start_Type
read_config 1 old.xyz

# Fragment_Files
species1/frag1/frag1.dat  1

END
"""


def _write_inputs(root: Path) -> tuple[Path, Path, Path, Path]:
    inp = root / "demo.inp"
    inp.write_text(INP, encoding="utf-8")
    (root / "demo.mcf").write_text(
        """# Atom_Info
1
1 C C 12.011 0.0

# Bond_Info
0

# Angle_Info
0
""",
        encoding="utf-8",
    )
    xyz = root / "demo.restart.xyz"
    xyz.write_text(
        "       1\n RESTART MC_STEP: 10\nC             1.00000000        2.00000000        3.00000000\n",
        encoding="utf-8",
    )
    hfile = root / "demo.restart.H"
    hfile.write_text(
        "\n".join(
            [
                "   1000.00000000000000",
                "   10.00000000000000  0.00000000000000  0.00000000000000",
                "   0.00000000000000  10.00000000000000  0.00000000000000",
                "   0.00000000000000  0.00000000000000  10.00000000000000",
                "",
                "   1",
                "   1   1",
                "",
            ]
        ),
        encoding="utf-8",
    )
    averages = root / "averages.json"
    averages.write_text(
        json.dumps({"volume_average": 512.0, "cubic_length_A": 8.0}) + "\n",
        encoding="utf-8",
    )
    return inp, xyz, hfile, averages


BOUNDARY_MCF = """# Atom_Info
4
1 C C 12.0 0.0
2 C C 12.0 0.0
3 O O 16.0 0.0
4 H H 1.0 0.0

# Bond_Info
3
1 1 2 fixed 0.4
2 2 3 fixed 0.4123105626
3 3 4 fixed 0.3741657387

# Angle_Info
2
1 1 2 3 fixed 165.9637565
2 2 3 4 fixed 157.2149227

# Dihedral_Info
1
1 1 2 3 4 none
"""


def _write_boundary_inputs(root: Path) -> tuple[Path, Path, Path, Path]:
    inp = root / "boundary.inp"
    inp.write_text(INP.replace("demo.mcf 1", "boundary.mcf 1"), encoding="utf-8")
    (root / "boundary.mcf").write_text(BOUNDARY_MCF, encoding="utf-8")
    xyz = root / "boundary.restart.xyz"
    xyz.write_text(
        """       4
 RESTART MC_STEP: 10
C             9.70000000        0.00000000        0.00000000
C             0.10000000        0.00000000        0.00000000
O             0.50000000        0.10000000        0.00000000
H             0.80000000        0.30000000        0.20000000
""",
        encoding="utf-8",
    )
    hfile = root / "boundary.restart.H"
    hfile.write_text(
        "\n".join(
            [
                "   1000.00000000000000",
                "   10.00000000000000  0.00000000000000  0.00000000000000",
                "   0.00000000000000  10.00000000000000  0.00000000000000",
                "   0.00000000000000  0.00000000000000  10.00000000000000",
                "",
                "   1",
                "   1   1",
                "",
            ]
        ),
        encoding="utf-8",
    )
    averages = root / "boundary.averages.json"
    averages.write_text(
        json.dumps({"volume_average": 512.0, "cubic_length_A": 8.0}) + "\n",
        encoding="utf-8",
    )
    return inp, xyz, hfile, averages


def _internal_vectors(coords: object, cell: object) -> tuple[object, object, object]:
    """Bond vectors, an angle cosine, and a dihedral cosine for the test chain."""
    import numpy as np

    xyz = np.asarray(coords, dtype=float)
    h = np.asarray(cell, dtype=float)
    b01 = _mic_delta(xyz[1] - xyz[0], h)
    b12 = _mic_delta(xyz[2] - xyz[1], h)
    b23 = _mic_delta(xyz[3] - xyz[2], h)
    angle_cos = float(np.dot(-b01, b12) / np.linalg.norm(b01) / np.linalg.norm(b12))
    n1 = np.cross(b01, b12)
    n2 = np.cross(b12, b23)
    dihedral_cos = float(np.dot(n1, n2) / np.linalg.norm(n1) / np.linalg.norm(n2))
    return np.vstack((b01, b12, b23)), angle_cos, dihedral_cos


def _section_values(text: str, header: str) -> list[str]:
    """Data lines following ``header`` until the next section."""
    lines = text.splitlines()
    collecting = False
    values: list[str] = []
    for line in lines:
        stripped = line.strip()
        if stripped.startswith("#"):
            if collecting:
                break
            collecting = stripped.split()[0:2] == header.split()
            continue
        if collecting and stripped and not stripped.startswith("!"):
            values.append(stripped)
    return values


def main() -> int:
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        inp, xyz, hfile, averages = _write_inputs(root)
        out = root / "nvt-out"
        build_restart(
            inp=inp,
            out_dir=out,
            ensemble="nvt",
            prp=None,
            averages=averages,
            restart_xyz=xyz,
            restart_h=hfile,
            discard_sweeps=None,
            run_name="demo-nvt",
        )
        written = (out / "demo-nvt.inp").read_text(encoding="utf-8")
        if "npt_mc" in written or "# Pressure_Info" in written or "# Prob_Volume" in written:
            raise SystemExit("NVT input still contains NPT-only sections")
        if "nvt_mc" not in written:
            raise SystemExit("NVT sim type missing")
        if "read_config 1 demo-nvt.xyz" not in written:
            raise SystemExit(f"unexpected start line:\n{written}")
        if abs(float(_section_values(written, "# Box_Info")[2]) - 8.0) > 1e-8:
            raise SystemExit("box length was not set to 8 A")
        probs = []
        for header in ("# Prob_Translation", "# Prob_Rotation", "# Prob_Regrowth"):
            probs.append(float(_section_values(written, header)[0]))
        if abs(sum(probs) - 1.0) > 1e-6:
            raise SystemExit(f"NVT probabilities sum to {sum(probs)}, not 1")

        scaled = read_xyz(out / "demo-nvt.xyz")
        expected = [0.8, 1.6, 2.4]
        if scaled.n_frames != 1 or max(abs(a - b) for a, b in zip(scaled.coords[0, 0], expected)) > 1e-8:
            raise SystemExit(f"coordinates were not scaled: {scaled.coords}")
        cell = read_h(out / "demo-nvt.H").frames[0]
        if abs(cell.volume - 512.0) > 1e-6 or abs(cell.cell[0, 0] - 8.0) > 1e-8:
            raise SystemExit(f"scaled cell is wrong: volume {cell.volume}, H {cell.cell}")

        boundary_inp, boundary_xyz, boundary_h, boundary_averages = _write_boundary_inputs(root)
        boundary_out = root / "boundary-out"
        build_restart(
            inp=boundary_inp,
            out_dir=boundary_out,
            ensemble="nvt",
            prp=None,
            averages=boundary_averages,
            restart_xyz=boundary_xyz,
            restart_h=boundary_h,
            discard_sweeps=None,
            run_name="boundary-nvt",
        )
        old_cell = read_h(boundary_h).frames[0].cell
        new_cell = read_h(boundary_out / "boundary-nvt.H").frames[0].cell
        old_coords = read_xyz(boundary_xyz).coords[0]
        new_coords = read_xyz(boundary_out / "boundary-nvt.xyz").coords[0]
        old_vectors, old_angle, old_dihedral = _internal_vectors(old_coords, old_cell)
        new_vectors, new_angle, new_dihedral = _internal_vectors(new_coords, new_cell)
        if not np.allclose(old_vectors, new_vectors, rtol=0.0, atol=2.0e-8):
            raise SystemExit("PBC-spanning molecule had its bond vectors changed")
        if abs(old_angle - new_angle) > 2.0e-8 or abs(old_dihedral - new_dihedral) > 2.0e-8:
            raise SystemExit("PBC-spanning molecule had its angle or dihedral changed")
        mcf = read_mcf(root / "boundary.mcf")
        old_contiguous = _unwrapped_molecule(old_coords, mcf, old_cell)
        new_contiguous = _unwrapped_molecule(new_coords, mcf, new_cell)
        masses = np.array([atom.mass for atom in mcf.atoms])
        old_fractional_com = np.linalg.solve(old_cell, np.average(old_contiguous, axis=0, weights=masses))
        new_fractional_com = np.linalg.solve(new_cell, np.average(new_contiguous, axis=0, weights=masses))
        if not np.allclose(old_fractional_com % 1.0, new_fractional_com % 1.0, rtol=0.0, atol=2.0e-8):
            raise SystemExit("molecular COM fractional coordinate changed")

        disconnected = root / "disconnected.mcf"
        disconnected.write_text(
            """# Atom_Info
2
1 C C 12.0 0.0
2 H H 1.0 0.0

# Bond_Info
0

# Angle_Info
0
""",
            encoding="utf-8",
        )
        try:
            _unwrapped_molecule(
                np.array([[1.0, 1.0, 1.0], [9.0, 1.0, 1.0]]),
                read_mcf(disconnected),
                old_cell,
            )
        except ValueError as exc:
            if "needs bonds" not in str(exc):
                raise
        else:
            raise SystemExit("disconnected multi-atom topology was accepted")

        bad_count_h = root / "bad-count.restart.H"
        bad_count_h.write_text(
            hfile.read_text(encoding="utf-8").replace("   1   1", "   1   2"),
            encoding="utf-8",
        )
        try:
            build_restart(
                inp=inp,
                out_dir=root / "bad-count",
                ensemble="nvt",
                prp=None,
                averages=averages,
                restart_xyz=xyz,
                restart_h=bad_count_h,
                discard_sweeps=None,
                run_name="bad-count",
            )
        except ValueError as exc:
            if "atom count disagrees" not in str(exc):
                raise
        else:
            raise SystemExit("restart .H / XYZ count mismatch was accepted")

        movie = root / "movie.xyz"
        movie.write_text(xyz.read_text(encoding="utf-8") * 2, encoding="utf-8")
        try:
            build_restart(
                inp=inp,
                out_dir=root / "bad",
                ensemble="nvt",
                prp=None,
                averages=averages,
                restart_xyz=movie,
                restart_h=hfile,
                discard_sweeps=None,
                run_name="nope",
            )
        except ValueError as exc:
            if "frames" not in str(exc):
                raise
        else:
            raise SystemExit("multi-frame movie was accepted")

        gcmc_out = root / "gcmc-out"
        build_restart(
            inp=inp,
            out_dir=gcmc_out,
            ensemble="gcmc",
            prp=None,
            averages=averages,
            restart_xyz=xyz,
            restart_h=hfile,
            discard_sweeps=None,
            run_name="demo-gcmc",
        )
        gcmc = (gcmc_out / "demo-gcmc.inp").read_text(encoding="utf-8")
        if "# Chemical_Potential_Info" not in gcmc or "# Pressure_Info" in gcmc:
            raise SystemExit("GCMC template did not replace the pressure section")
        if _section_values(gcmc, "# Chemical_Potential_Info") != ["0.0"]:
            raise SystemExit("chemical potential placeholder missing")
        if _section_values(gcmc, "# Prob_Insertion")[0] != "0.150000":
            raise SystemExit("insertion probability placeholder missing")
        if _section_values(gcmc, "# Prob_Deletion") != ["0.150000"]:
            raise SystemExit("deletion probability placeholder missing")

    print("npt_restart synthetic tests passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
