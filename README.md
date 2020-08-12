<a href="https://cassandra.nd.edu">
<p align="center">
  <img src="https://cassandra.nd.edu/images/visual-identity/cassandra_logo_square-full.png" width="500" title="Cassandra Logo">
</p>
</a>

![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)
[![Build Status](https://dev.azure.com/MaginnGroup/Cassandra/_apis/build/status/MaginnGroup.Cassandra?branchName=master)](https://dev.azure.com/MaginnGroup/Cassandra/_build/latest?definitionId=2&branchName=master)
[![PRs Welcome](https://img.shields.io/badge/PRs-welcome-brightgreen.svg?style=flat-square)](http://makeapullrequest.com)
[![Chat on Slack](https://img.shields.io/badge/chat-on_slack-green.svg?longCache=true&style=flat&logo=slack)](https://join.slack.com/t/cassandra-nd/shared_invite/enQtNDU0NDIyMDc5ODEyLWE5Yzc1YmQ3MzdhNWVmYTkyNjY0YWQ3ZDJmMjQ5ZTliYTc3MTNlZTg2MzgxOTdjY2Y1ZTNiMjVhMmE4NGEzYWM)

## Contents
1. [Overview](#overview)
2. [Features](#features)
3. [Installation](#installation)
4. [Getting Started](#starting)
5. [Documentation](https://github.com/MaginnGroup/Cassandra/releases/latest/download/user_guide.pdf)
6. [Forum](https://cassandra.nd.edu/index.php/community)
7. [Contribute](#contributing)
8. [Citation](#citation)


## <a name="overview"></a>Overview

[Cassandra](https://cassandra.nd.edu/) is an open source Monte Carlo package developed in the [Maginn group](http://sites.nd.edu/maginn-group/) at the [University of Notre Dame](http://www.nd.edu) to perform atomistic simulations of molecules composed of rings, chains, or both. Cassandra is suited to compute the thermodynamic properties of fluids and phase equilibria. It handles a standard "Class I"-type force field having fixed bond lengths. Cassandra uses OpenMP parallelization and comes with a number of scripts, utilities and examples to help with simulation setup. It is released under the GNU General Public License.

A Python interface to Cassandra is under development as well. Further details are
available [here](https://mosdef-cassandra.readthedocs.io/en/latest/).

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

## <a name="installation"></a> Installation

We offer the option to install from source or through `conda-forge`.
The simplest way to get Cassandra is through `conda-forge`.

### Installation with conda

If you already have [conda](https://docs.conda.io/en/latest/miniconda.html)
installed, you can install Cassandra into a new `conda` environment
with the following command:

```
> conda create --name cassandra -c conda-forge cassandra
```

Once you activate the environment, Cassandra should be accesible on your `PATH`.
You can test with the following:

```
> conda activate cassandra
> which cassandra.exe
> which library_setup.py
> which mcfgen.py
```

If the installation was successful, you should see the location of each of
those files. The `conda-forge` installation uses `gfortran` with
OPENMP parallelization. The number of parallel threads can be controlled
by setting the `OMP_NUM_THREADS` environment variable. For example, in `bash`:

```
export OMP_NUM_THREADS=4
```

Cassandra parallelizes reasonably up to 8-12 threads.

### Installation from source

The source code tarball for the latest release is available through our [releases
page](https://github.com/maginngroup/cassandra/releases/latest). Alternatively,
you can clone the GitHub repository to get the bleeding edge version:

```
> git clone
https://github.com/MaginnGroup/Cassandra.git
```

In either case, the ```/Src/``` directory contains the Makefiles that you can use to compile the code.
Makefiles contain the compilation options and set of directives used to automate the build.
At present, Makefiles for the Intel Fortran Compiler, gfortran and
Portland group compiler (PGI) are included in the distribution.

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
OMP_NUM_THREADS will have to be set to the number of threads you want to
run the simulation on. For example, for 12-thread simulation, the
following command is used in tcsh:

```> setenv OMP_NUM_THREADS 12```

or alternatively in bash:

```> export OMP_NUM_THREADS=12```

## <a name="started"></a> Getting Started

**/Examples**. This directory contains examples of short Cassandra simulations in
NVT, NPT, grand canonical and Gibbs ensembles for a number of systems molecules with varying degrees
of conformational complexity (LJ particles, branch points, ring moieties) and
those requiring computation of electrostatic interactions are included.

**/Documentation**. Chapter 3 of the [documentation](https://github.com/MaginnGroup/Cassandra/releases/latest/download/user_guide.pdf) provides an overview of the files and workflow to setup a Cassandra simulation. Additionally, the [workshop materials](https://cassandra.nd.edu/index.php/documentation) contain more examples and slides that have been used for teaching.

## <a name="contributing"></a> Contributing

Depending on whether you are a user or a developer, you can contribute to Cassandra in the following ways:

1. **User**
    - Ask usage questions and participating in discussions using the [forum](https://cassandra.nd.edu/index.php/community).
    - Request new features using the [forum](https://cassandra.nd.edu/index.php/community).
    - Reporting bugs using GitHub issues (preferred) or [forum](https://cassandra.nd.edu/index.php/community).

2. **Developer**
    - Contribute code or documentation using GitHub pull request workflow. See a primer on Cassandra Software Development on GitHub. See [CONTRIBUTING](https://github.com/MaginnGroup/Cassandra/blob/master/CONTRIBUTING.md) for more details on the contributing process.
    - Ask development questions (i.e. you are trying to fix a bug or contribute with a new feature) using [Slack](https://cassandra-nd.slack.com/messages/general/). 

## <a name="citation"></a> Citation

 J. K. Shah, E. Marin‐Rimoldi, R. G. Mullen, B. P. Keene, S. Khan, A. S. Paluch, N. Rai, L. L. Romanielo, T. W. Rosch, B. Yoo, E. J. Maginn. J. Comput. Chem. 2017, 38, 1727–1739. [DOI: 10.1002/jcc.24807](https://onlinelibrary.wiley.com/doi/abs/10.1002/jcc.24807)
