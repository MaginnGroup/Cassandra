#!/usr/bin/env python
""""""
import argparse
import io
import numpy as np
from pathlib import Path


def lammpstrjconvert(
    lammpstrjpath, n_list, fstr="%f", Hpath=None, xyzpath=None, getframes=None
):
    """Convert LAMMMPS custom dump file to .xyz and .H files to be read by Cassandra.

    Inputs:
        lammpstrjpath: path-like object, such as pathlib.Path or string, containing
            the path to the LAMMPS dump file to read and convert.
        n_list: list of integers containing the number of molecules of each species
            in the order in which the species are listed in the Cassandra input file.
        fstr: format string designating the format with which to write
            the coordinate floats in the .xyz file.  The default is "%f".
        Hpath: path-like object, such as pathlib.Path or string, containing
            the path to which to write the .H file.  If None (the default), this is
            lammpstrjpath with the parent directories and ".lammpstrj" suffix (if present)
            stripped and ".H" appended.
        xyzpath: path-like object, such as pathlib.Path or string, containing
            the path to which to write the .H file.  If None (the default), this is
            lammpstrjpath with the parent directories and ".lammpstrj" suffix (if present)
            stripped and ".xyz" appended.
        getframes: sequence or 1-D array of integers designating the zero-based indices of
            specific frames to write.  If None (the default), all frames are written.
            If empty, the function effectively does nothing.  The frames are written in
            the same order with which they are listed in getframes, which is not
            necessarily in ascending order.
    """
    ltpath = Path(lammpstrjpath)
    full_fstr = "X " + fstr + " " + fstr + " " + fstr
    colnames = ("id", "xu", "yu", "zu")

    if Hpath is None:
        Hpath = Path(ltpath.stem + ".H")
    else:
        Hpath = Path(Hpath)
    if xyzpath is None:
        xyzpath = Path(ltpath.stem + ".xyz")
    else:
        xyzpath = Path(xyzpath)

    nonsorted = False

    if getframes is None:
        frame_array = []
    elif len(getframes):
        frame_array = np.array(getframes)
        nonsorted = any(
            [frame_array[i] >= frame_array[i + 1] for i in range(len(frame_array) - 1)]
        )
    else:
        return

    with ltpath.open() as ltfile, Hpath.open("w") as Hfile, xyzpath.open(
        "w"
    ) as xyzfile:
        eofreached = False
        # Define nested functions findheading and convert_frame
        def findheading(tgt, lineadvance=True):
            eof_flag = False
            tgt_hit = False
            while not (tgt_hit):
                if lineadvance:
                    thislinestr = ltfile.readline()
                else:
                    thislinestr = ltfile.readline(len(tgt))
                if len(tgt) <= len(thislinestr):
                    tgt_hit = thislinestr[0 : len(tgt)] == tgt
                elif not thislinestr:
                    eof_flag = True
                    tgt_hit = True
            return eof_flag

        def convert_frame():
            xy = 0.0
            xz = 0.0
            yz = 0.0
            timestep = int(ltfile.readline().strip())
            eofreached = findheading("ITEM: NUMBER OF ATOMS")
            n_atoms = int(ltfile.readline().strip())
            eofreached = findheading("ITEM: BOX BOUNDS", False)
            boxbounds_list = ltfile.readline().strip().split()
            extent_x_str = ltfile.readline().strip().split()
            extent_y_str = ltfile.readline().strip().split()
            extent_z_str = ltfile.readline().strip().split()
            xlo = float(extent_x_str[0])
            xhi = float(extent_x_str[1])
            ylo = float(extent_y_str[0])
            yhi = float(extent_y_str[1])
            zlo = float(extent_z_str[0])
            zhi = float(extent_z_str[1])
            if boxbounds_list[0] == "xy":
                triclinic_flag = True
                xy = float(extent_x_str[2])
                xz = float(extent_y_str[2])
                yz = float(extent_z_str[2])
                xlo -= min([0.0, xy, xz, xy + xz])
                xhi -= max([0.0, xy, xz, xy + xz])
                ylo -= min([0.0, yz])
                yhi -= max([0.0, yz])
            xx = xhi - xlo
            yy = yhi - ylo
            zz = zhi - zlo
            a = [xx, 0, 0]
            b = [xy, yy, 0]
            c = [xz, yz, zz]
            lmat = np.array([a, b, c])
            volume = np.inner(a, np.cross(b, c))
            # boxes always have origin at (0,0,0) in Cassandra, but not always in lammps
            box_center = np.sum(lmat, axis=1) * 0.5 + np.array([xlo, ylo, zlo])
            eofreached = findheading("ITEM: ATOMS ", False)
            coldict = {
                colname: i
                for i, colname in enumerate(ltfile.readline().strip().split())
            }  # gives indices of columns
            xyz_buffer = io.StringIO()
            for i in range(n_atoms):
                xyz_buffer.write(ltfile.readline())
            xyz_buffer.seek(0)
            xyz = np.loadtxt(
                xyz_buffer, usecols=[coldict[colname] for colname in colnames]
            )
            xyz_buffer.close()
            xyz = xyz[np.argsort(xyz[:, 0])]
            xyz[:, 1:] -= box_center  # center the box at (0, 0, 0)
            nspecies = len(n_list)
            Hfile.write("{:^26.17g}\n".format(volume))
            Hfile.write(
                "{:^26.17g}{:^26.17g}{:^26.17g}\n".format(
                    lmat[0, 0], lmat[0, 1], lmat[0, 2]
                )
            )
            Hfile.write(
                "{:^26.17g}{:^26.17g}{:^26.17g}\n".format(
                    lmat[1, 0], lmat[1, 1], lmat[1, 2]
                )
            )
            Hfile.write(
                "{:^26.17g}{:^26.17g}{:^26.17g}\n\n".format(
                    lmat[2, 0], lmat[2, 1], lmat[2, 2]
                )
            )
            Hfile.write("{:>12d}\n".format(nspecies))
            for i in range(nspecies):
                Hfile.write("{:>12d}{:>12d}\n".format(i + 1, n_list[i]))
            # write xyz file
            xyzfile.write("{:>12d}\n".format(n_atoms))
            xyzfile.write(" TIMESTEP: {:>11d}\n".format(timestep))
            np.savetxt(xyzfile, xyz[:, 1:], fmt=full_fstr)

        # End of nested function definitions
        if nonsorted:
            streampos = {}
            framerange = np.arange(max(frame_array + 1))
            for iframe in framerange:
                eofreached = findheading("ITEM: TIMESTEP")
                if eofreached:
                    raise ValueError(
                        "Frame indices specified in getframes exceed the maximum frame index "
                        + str(iframe - 1)
                    )
                if iframe in frame_array:
                    streampos[iframe] = ltfile.tell()
            for iframe in frame_array:
                ltfile.seek(streampos[iframe])
                convert_frame()
        else:
            eofreached = findheading("ITEM: TIMESTEP")
            iframe = 0
            while not eofreached:
                if getframes is None or iframe in frame_array:
                    convert_frame()
                eofreached = findheading("ITEM: TIMESTEP")
                if getframes is not None:
                    if iframe == frame_array[-1]:
                        eofreached = True
                    elif eofreached:
                        raise ValueError(
                            "Frame indices specified in getframes exceed the maximum frame index "
                            + str(iframe)
                        )
                iframe += 1


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--format", default="%f")
    parser.add_argument("fname")
    parser.add_argument("nmols", nargs="+", type=int)
    args = parser.parse_args()
    lammpstrjconvert(lammpstrjpath=args.fname, n_list=args.nmols, fstr=args.format)
