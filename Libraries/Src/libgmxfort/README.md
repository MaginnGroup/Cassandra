# libgmxfort

**Repository update 10/03/2020: This project has been archived and is read-only. I can no
longer support users with the library. I hope it was useful to you in your analysis.
Please fork and keep developing it if you have found it useful.**

[![DOI](https://zenodo.org/badge/71169051.svg)](https://zenodo.org/badge/latestdoi/71169051)

This is a modern Fortran library for reading in an analyzing
[GROMACS](http://www.gromacs.org/) compressed
trajectory files (.xtc) and index files (.ndx). Some helpful utilities functions
are also included like periodic boundary condition and distance functions. It
uses an object-oriented philosophy for reading in, storing, and retrieving
simulation data for analysis.

Also check out [dcdfort](https://github.com/wesbarnett/dcdfort).

## Requirements

[xdrfile](https://github.com/wesbarnett/libxdrfile)>=2.1.2 is required.

## Compilation

After cloning the repository, or extracting the release tarball, cd into the
repository. Then:

```bash
mkdir build
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=/usr/local
make
```

## Testing

To test your build, do:

```bash
make test
```

If the tests do not pass, please file an issue.

## Installation

The following will install the library to the location specified by
`DCMAKE_INSTALL_PREFIX`.

```bash
make install
```

## Usage

Compile your program with `-lgmxfort`. You may also need to use `-I` to point to
where the modules files are even with all of the right environment variables set
(by default at `/usr/local/include`). 

### Linking other cmake projects

A file is included to easily link other cmake projects to the gmxfort
installation. Use `find_package ( gmxfort )` and the variables
`gmxfort_INCLUDE_DIRS` and `gmxfort_LIBRARIES`.

### pkg-config

A pkg-config file is included, so that it can
be used in your program compilations. You may need to set `PKG_CONFIG_PATH` to
find the file (by default in the directory `/usr/local/lib/pkgconfig`)

### API

To use the library always put `use gmxfort_trajectory` first in order to use the
`Trajectory` class and `use gmxfort_utils` in order to use any of the other
utilities.  There is an example in the `example` folder on how to do this.
Additionally [sasa](https://github.com/wesbarnett/sasa) uses libgmxfort.

Typically you will open a trajectory file (and optionally a corresponding index
file). Then you will read in the entire trajectory file at once, or you can read
it in in chunks. Then you should close the trajectory file when done.

The simplest way to use this library is to construct a `Trajectory` object and
then use the `read()` method:

```fortran
use gmxfort_trajectory
implicit none
type(Trajectory) :: trj
call trj%read("traj.xtc")
```

The `read()` method opens the xtc file, reads in all information, and then
closes it. The `trj` object in this example now stores all of the coordinates and
information from the .xtc file.

If you have a corresponding index file, you can add a second argument to
`open`:

```fortran
call trj%read("traj.xtc", "index.ndx")
```

Now information regarding the index groups is stored in memory.

If you want to read in the trajectory file in frame-by-frame use `read_next()`
instead of `read()`. To use this, you must additionally open and close the xtc
file on your own. By default it reads in one frame:

```fortran
integer :: n
call trj%open("traj.xtc", "index.ndx")
n = trj%read_next()
call trj%close()
```

To read in more than one, specify an argument. The following reads in 10 frames:

```fortran
n = trj%read_next(10)
```

`read_next()` returns the number of frames actually read in. It is a function,
and not a subroutine. This is useful for using it with a `do while` loop. For
example:

```fortran
use gmxfort_trajectory

implicit none

type(Trajectory) :: trj
integer :: i, n

call trj%open("traj.xtc", "index.ndx")

n = trj%read_next(10)
do while (n > 0)
    do i = 1, n
        ! do some things with the frames read in
    end do
    n = trj%read_next(10)
end do

call trj%close()
```

After calling `read()` or `read_next()` every atom's coordinates are accessible
via the `x()` method. For example, to get the coordinates of the first atom in
the first frame you would do the following. The frame is the first argument and
the atom number is the second argument. 

```fortran
real :: myatom(3)
! ...
myatom = trj%x(1, 1)
```

**Note**: Fortran uses one-based indexing, and that convention is retained here.

If you read in an index file, you can get atom coordinates in relationship to
that. The following gets the fifth atom in index group `C` in the `10`th frame:

```fortran
myatom = trj%x(10, 5, "C")
```

**Note**: If you have more than one group in your index file with the same name,
this will simply use the first group with that name. It's best not to repeat
group names in your index file. The library will give you a warning if it finds
that an index name is duplicated:

    LIBGMXFORT WARNING: Index group OW specified more than once in index file.

Note that when you use `x()` you will still have to give it the frame number as
the first argument even if you only read in one frame with `read_next()`.  You
can always get the number of frames in a trajectory file object with the
`nframes` member:

```fortran
integer :: n
! ...
n = trj%nframes
```

You can also get the number of atoms with the `natoms()` method:

```fortran
integer :: n
! ...
n = trj%natoms()
```

If you want to know how many atoms are in an index group include the group name
as an argument. In this example the group name is "C":

```fortran
n = trj%natoms("C")
```

To get the box coordinates, use `box`. The following gets the box of the `2`nd
frame:

```fortran
real :: mybox(3,3)
! ...
mybox = trj%box(2)
```

You can also get the simulation time and step corresponding with a frame you
read in, using `time` and `step`, respectively. The following get the time associated
with the first frame read in:

```fortran
real :: mytime
! ...
mytime = trj%time(1)
```

And now the step for the same:

```fortran
integer :: mystep
! ...
mystep = trj%step(1)
```

As shown above, the most common use of this library is to use `read()` or
`read_next()` to save all atom locations and then use getters like `x()` and
`natoms()` to get information about them by specifying an index group as an
argument. To save memory, you can save just a specific index group with
`read()`:


```fortran
trj%read(xtcfile, ndxfile, "C")
```

If you do this, you only have access to the group above, and you should never 
pass an index group name to getters like `x()`, since only one group is
available. If you do specify a group in a getter after already specifying it in
`read()` or `read_next()`, you will get this error:

    LIBGMXFORT ERROR: Do not specify an index group in x() when already specifying an index group with read() or read_next().

There are several functions and subroutines in the `gmxfort_utils` module,
including periodic boundary and distance calculations. Check out the source file
for what is available.

## License

libgmxfort

Copyright (C) 2016-2019 James W. Barnett

This program is free software; you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation; either version 2 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
Street, Fifth Floor, Boston, MA 02110-1301 USA.
