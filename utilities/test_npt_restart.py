#!/usr/bin/env python3
#*******************************************************************************
#  test_npt_restart.py
#
#  A one-atom cubic box is scaled from L = 10 A to the average L = 8 A.
#  The generated NVT input must drop the volume move.  The GCMC input must
#  keep a placeholder chemical potential instead of inventing one.
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

_REPO_ROOT = Path(__file__).resolve().parent.parent
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from diagnostics.io.hfile import read_h
from diagnostics.io.xyz import read_xyz
from utilities.npt_restart import build_restart


INP = """! leftover comment that Cassandra ignores
# Run_Name
demo-npt

# Sim_Type
npt_mc

# Nbr_Species
1

# Molecule_Files
demo.mcf 5

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
read_config 5 old.xyz

# Fragment_Files
species1/frag1/frag1.dat  1

END
"""


def _write_inputs(root: Path) -> tuple[Path, Path, Path, Path]:
    inp = root / "demo.inp"
    inp.write_text(INP, encoding="utf-8")
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
                "   1   5",
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
        if "read_config 5 demo-nvt.xyz" not in written:
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
