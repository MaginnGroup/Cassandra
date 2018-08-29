<a href="https://cassandra.nd.edu">
<p align="center">
  <img src="https://cassandra.nd.edu/images/visual-identity/cassandra_logo_square-full.png" width="500" title="Cassandra Logo">
</p>
</a>

![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)
[![PRs Welcome](https://img.shields.io/badge/PRs-welcome-brightgreen.svg?style=flat-square)](http://makeapullrequest.com)

## Contents
1. [Overview](#overview)
2. [Features](#features)
3. [How to Get and Compile Cassandra](#obtaining)
4. [Getting Started](#starting)
5. [Documentation](https://cassandra.nd.edu/index.php/documentation)
6. [Forum](https://cassandra.nd.edu/index.php/community)
7. [Contribute]
8. [Citation](#citation)


## <a name="overview"></a>Overview

[Cassandra](https://cassandra.nd.edu/) is an open source Monte Carlo package developed in the [Maginn group](http://sites.nd.edu/maginn-group/) at the [University of Notre Dame](http://www.nd.edu) to perform atomistic simulations of molecules composed of rings, chains, or both. Cassandra is suited to compute the thermodynamic properties of fluids and phase equilibria. It handles a standard "Class I"-type force field having fixed bond lengths. Cassandra uses OpenMP parallelization and comes with a number of scripts, utilities and examples to help with simulation setup. It is released under the GNU General Public License.

## <a name="features"></a> Features

The following features are supported in version 1.2:

1. Ensembles
    - Canonical
    - Isothermal-Isobaric
    - Grand-canonical
    - Constant volume Gibbs
    - Constant pressure Gibbs

2. Moves
    - Translation
    - Rotation
    - Box volume
    - Molecule regrowth
    - Molecule insertion
    - Molecule deletion
    - Molecule swap
    - Dihedral change
    - Angle change
    - Atom displacement
    - Ring flip

3. Potentials
    - LJ and Mie
    - Fixed bond lengths
    - Harmonic and fixed angles
    - OPLS, CHARMM and harmonic proper dihedrals
    - Harmonic improper dihedrals
    - Point partial charges

4. Other
    - Cut-off schemes: cut, cut and shift, cut and switch, standard long-tail correction
    - Electorstatics: Ewald or damped shifted force methods

## <a name="obtaining"></a> How to Get and Compile Cassandra

You can get a free copy of Cassandra by [downloading the tarball](https://cassandra.nd.edu/index.php/download). Alternatively, you can use GitHub to get a copy of the bleeding edge version:

```
> git clone
https://github.com/MaginnGroup/Cassandra.git
```

The ```/Src/``` directory contains the Makefiles that you can use to compile the code. Makefiles contain the compilation options and set of directives used to automate the build. At present, Makefiles for the
Intel Fortran Compiler (12.1), gfortran i4.4.7 20120313 (Red Hat 4.4.7-4) and
Portland group compiler (PGI) (14.6) are included in the distribution.

To compile Cassandra, first remove any object files using the 'clean' command

```> make clean```

The following table provides examples of the commands you need to use to compile the source code as provided:

|   Build type    | Intel Compiler      |    GNU Fortran   |
|---    |---    |---    |
|  Debug     |   ```make -f Makefile```    | ```make -f Makefile.gfortran```   |
|  OpenMP    |   ```make -f Makefile.intel.openMP```     |   ```make -f Makefile.gfortran.openMP```    |

Each of these commands will produce an executable with the name specified in the field EXEC NAME in the relevant Makefile.

You can take these Makefiles as templates
creating your own customized Makefile. Depending on the architecture
of the machine you are using, you will need to change compilation
options and flags.  Note that to change compiler options appropriate
for your environment, the F90FLAGS line in the Makefiles can be
edited. You can also modify the optimization options to
improve the speed of the code or operate in debug mode.

Before running the openMP enabled executable, the environment variable
OMP_NUM_THREADS will have to be set to the number of cores you want to
run the simulation on. For example, for 12-core simulation, the
following command is used in tcsh.

```> setenv OMP_NUM_THREADS 12```

## <a name="started"></a> Getting Started

**Examples**. This directory contains examples of short Cassandra simulations in
NVT, NPT, grand canonical and Gibbs ensembles for a number of systems molecules with varying degrees
of conformational complexity (LJ particles, branch points, ring moieties) and
those requiring computation of electrostatic interactions are included.

**Documentation**. Chapter 3 of the [documentation](https://cassandra.nd.edu/index.php/documentation/) provides an overview of the files and workflow to setup a Cassandra simulation. Additionally, the [workshop materials](https://cassandra.nd.edu/index.php/documentation) contain more examples and slides that have been used for teaching.

## <a name="citation"></a> Citation

 J. K. Shah, E. Marin‐Rimoldi, R. G. Mullen, B. P. Keene, S. Khan, A. S. Paluch, N. Rai, L. L. Romanielo, T. W. Rosch, B. Yoo, E. J. Maginn. J. Comput. Chem. 2017, 38, 1727–1739. [DOI: 10.1002/jcc.24807](https://onlinelibrary.wiley.com/doi/abs/10.1002/jcc.24807)
