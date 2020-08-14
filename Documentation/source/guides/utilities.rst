Utilities
=========

.. _sec:mcfgen:

Generate a Molecular Connectivity File
--------------------------------------

The script ``mcfgen.py`` is a tool that aims to ease the setup of molecular
connectivity files from scratch (see the :ref:`sec:mcf_file` section to learn
more about MCFs), as the generation of these files by hand can be error prone.
In this section, a pentane MCF will be generated to demonstrate the use of this
tool. The Transferable Potentials for Phase Equilibria (TraPPE) force field will
be used to represent the pentane molecular interactions. This force field
involves a pairwise-additive 12-6 Lennard-Jones potential to represent the
dispersion-repulsion interactions. Additionally, bond angles and dihedral angles
are represented through harmonic and OPLS functional forms, respectively. Bond
lengths are kept constant. The force field mathematical expression becomes

.. math::

   \begin{aligned}
   U = \sum_{angles} & K_\theta(\theta-\theta_0)^2 + \\
   \sum_{dihedrals} & \frac{1}{2}K_1[1+cos(\phi)]+\frac{1}{2}K_2[1-cos(2\phi)] + \frac{1}{2}K_3[1+cos(3\phi)]+\frac{1}{2}K_4[1-cos(4\phi)] + \\
   \sum_{i} \sum_{i>j} & 4 \epsilon_{ij} \left [  \left ( \frac {\sigma_{ij}} { r_{ij} }\right )^{12} - \left ( \frac {\sigma_{ij}} { r_{ij} }\right )^{6}\ \right ]\end{aligned}

First, generate (or obtain) a PDB file or a CML file. To generate a PDB or CML
file, software such as Gaussview or Avogadro can be used.  Alternatively, PDB
files can be downloaded from the internet (e.g., `<https://www.rcsb.org>`_).
In this example, a pentane PDB file using the program Gaussview v5.08 was
generated as shown below.

.. code-block::

    REMARK   1 File created by GaussView 5.0.8
    ATOM      1  C1  PENL    1       2.142   1.395  -8.932  1.00  0.00          C
    ATOM      2  C2  PENL    1       3.631   1.416  -8.537  1.00  0.00          C
    ATOM      3  C3  PENL    1       4.203  -0.012  -8.612  1.00  0.00          C
    ATOM      4  C4  PENL    1       5.691   0.009  -8.218  1.00  0.00          C
    ATOM      5  C5  PENL    1       5.691   0.009  -8.218  1.00  0.00          C
    TER
    CONECT    1    2
    CONECT    2    1    3
    CONECT    3    2    4
    CONECT    4    3    5
    CONECT    5    4

Append a column containing the atom types:

.. code-block::

    REMARK   1 File created by GaussView 5.0.8
    ATOM      1  C1  PENL    1       2.142   1.395  -8.932  1.00  0.00          C3  CH3
    ATOM      2  C2  PENL    1       3.631   1.416  -8.537  1.00  0.00          C2  CH2
    ATOM      3  C3  PENL    1       4.203  -0.012  -8.612  1.00  0.00          C2  CH2
    ATOM      4  C4  PENL    1       5.691   0.009  -8.218  1.00  0.00          C2  CH2
    ATOM      5  C5  PENL    1       5.691   0.009  -8.218  1.00  0.00          C3  CH3
    TER
    CONECT    1    2
    CONECT    2    1    3
    CONECT    3    2    4
    CONECT    4    3    5
    CONECT    5    4

Avogadro v1.1.1 can also be used to generate CML files. Below is an
example of a CML file generated using Avogadro.

.. code-block::

    <molecule>
     <atomArray>
      <atom id="a1" elementType="C" x3="-0.367789" y3="-0.161907" z3="0.206019"/>
      <atom id="a2" elementType="C" x3="-1.354811" y3="-1.178938" z3="-0.372241"/>
      <atom id="a3" elementType="C" x3="-2.735586" y3="-0.597632" z3="-0.678858"/>
      <atom id="a4" elementType="C" x3="-3.435276" y3="0.007943" z3="0.526735"/>
      <atom id="a5" elementType="C" x3="1.027694" y3="-0.340782" z3="-0.372648"/>
     </atomArray>
     <bondArray>
      <bond atomRefs2="a1 a2" order="1"/>
      <bond atomRefs2="a3 a4" order="1"/>
      <bond atomRefs2="a3 a2" order="1"/>
      <bond atomRefs2="a5 a1" order="1"/>
     </bondArray>
    </molecule>

Modify the pentane united atom CML file. Note that the atom type is
appended as a last column between quotation marks.

.. code-block::

    <molecule>
     <atomArray>
      <atom id="a1" elementType="C" x3="-0.367789" y3="-0.161907" z3="0.206019"/> "CH2"
      <atom id="a2" elementType="C" x3="-1.354811" y3="-1.178938" z3="-0.372241"/> "CH2"
      <atom id="a3" elementType="C" x3="-2.735586" y3="-0.597632" z3="-0.678858"/> "CH2"
      <atom id="a4" elementType="C" x3="-3.435276" y3="0.007943" z3="0.526735"/> "CH3"
      <atom id="a5" elementType="C" x3="1.027694" y3="-0.340782" z3="-0.372648"/> "CH3"
     </atomArray>
     <bondArray>
      <bond atomRefs2="a1 a2" order="1"/>
      <bond atomRefs2="a3 a4" order="1"/>
      <bond atomRefs2="a3 a2" order="1"/>
      <bond atomRefs2="a5 a1" order="1"/>
     </bondArray>
    </molecule>

In the terminal, run the following command:

.. code-block:: bash

    python mcfgen.py pentane.pdb –ffTemplate

This command will create an .ff file. The first three sections of the FF file
are displayed next. Do not modify these.

.. code-block::

    atomtypes
    2

    begin atom-atomtype
    1 CH3
    2 CH2
    3 CH2
    4 CH2
    5 CH3
    end atom-atomtype

    dihedraltype OPLS

The force field parameters for non-bonded (not shown), bonds, angle, dihedral
(not shown) and coulombic interactions (not shown) must be entered next to the
corresponding keyword. For example, the angle type CH3 CH2 CH2 has an angle of
114.0. This value must be placed next to the “Angle” keyword.

.. code-block::

    bonds
    CH2 CH2
    Length 1.54
    Constant fixed

    angles
    CH3 CH2 CH2
    Angle 114.0
    Constant 31250.0


For more examples of filled ff files, please refer to the examples
contained in the ``/Scripts/MCF_Generation/`` directory. Using the filled
.ff file, run:

.. code-block:: bash

    python mcfgen.py pentane.pdb

Check the file newly created pentane.mcf for any possible errors. This example
can be found in the directory ``/Scripts/MCF_Generation/PDB/``

Note that if an MCF for a rigid solid is being created, this last step
must include the ``--solid`` flag, as

.. code-block:: bash

    python mcfgen.py zeolite.pdb --solid


.. _sec:libgen:

Generate Library of Fragment Configurations
-------------------------------------------

The goal of the script ``library_setup.py`` is to automate the generation of
fragment libraries. As a starting point, the script requires the simulation
input file, and the MCF and PDB files for each of the species. To run this
script, type

.. code-block:: bash

    python library_setup.py $PATH$/cassandra.exe input_file.inp pdbfilespecies1.pdb pdfilespecies2.pdb ...

This script will create the necessary files to create the fragment libraries. It
will also run Cassandra to generate these libraries, whose location will be at
``/species?/frag?/frag?.inp``, where ’?’ refers to the species number, for
example, species 1, species 2 etc. Note that the script overwrites the section
of the input file where needed (i.e. ``# Fragment_Files``) with the aforementioned
directory locations.
