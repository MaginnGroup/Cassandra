Simulating Rigid Solids and Surfaces
====================================

Simulations involving a rigid solid or surface can be performed in constant
volume ensembles (i.e., NVT and GCMC). For example, an adsorption isotherm can
be computed with a GCMC simulation in which fluid molecules are inserted into a
crystalline solid. In addition to the files described in :ref:`ch:input_files`,
the following files are required to run a simulation with a rigid solid or
surface:

-  a molecular connectivity file with force field parameters for each atom in
   the solid (``*.mcf``)

-  a fragment library file listing the coordinates of each atom in the solid
   (``*.dat``)

-  a configuration file with the initial coordinates of the all atoms in the
   system (``*.xyz``)

The MCF and fragment library file can be created using the scripts discussed in
:ref:`sec:mcfgen` and :ref:`sec:libgen`.  Each of these scripts requires a
starting PDB configuration file. The input file also requires specific
parameters as given in the following section. Sample input scripts for GCMC
simulations of various fluids in silicalite are included in the Examples/GCMC
directory of the Cassandra distribution.

Input file
----------

The input file should be completed as described in :ref:`sec:input_file` but
with the following parameters:

-  Under the keyword ``# Prob_Translation``, the translation width for the solid
   is 0.

-  Under the keyword ``# Prob_Rotation``, the rotation width for the solid is 0.

-  Under the keyword ``# Prob_Regrowth``, the regrowth probability for the solid
   is 0.

-  The volume dimensions under the keyword ``# Box_Info`` must match the crystal
   dimensions.

-  Under the keyword ``# Start_Type``, the ``read_config`` or ``add_to_config``
   options must be used.

GCMC simulations require the following additional parameters:

-  Under the keyword ``# Prob_Insertion`` insertion method for the solid
   is ``none``.

-  Under the keyword ``# Chemical_Potential_Info``, no chemical
   potential is required for the solid.

.. _sec:solid_pdb:

PDB file
--------

A PDB configuration file for a zeolite can be created in the following
ways, among others:

-  Manually, with atomic coordinates from the literature. For example,
   the atomic coordinates of silicalite are published
   `here <https://doi.org/10.1021/j150615a020>`_.

-  From a Crystallographic Information File (CIF), which can be
   downloaded from the
   `Database of Zeolite Structures <http://www.iza-structure.org/databases>`_.
   From the home page, click
   on the menu "All codes" and select your zeolite. The website
   will display structural information about the zeolite and will have a
   link to download a CIF. The CIF contains information about the
   zeolite structure such as cell parameters, space groups, T and O atom
   coordinates. A CIF can be converted into a PDB file using either
   Mercury or VESTA, both of which are available to download for free.
   For example, using VESTA:

   #. From the File menu, click Open. Then download the CIF (e.g.
      MFI.cif)

   #. From the Objects menu, click Boundary. Enter the desired number of
      replicas along each axis. To output a single unit cell, enter -0
      to 1 in the x, y and z ranges. To output a 2x2x2 crystal, enter -1
      to 1 in the x, y and z ranges.

   #. From the File menu, click Export Data. Enter a filename ending
      with .pdb (e.g. MFI.pdb)

.. _sec:solid_mcf:

Molecular connectivity file
---------------------------

Since the solid structure will be rigid, no bond distances, angle
parameters or dihedral parameters are needed in the MCF. The PDB file
for the rigid solid does not list CONECT information, so the
``mcfgen.py`` script will not include bond, angle, or dihedral sections
in the force field template (``*.ff``) or MCF. The number of fragments will
be zero. Only nonbonded parameters are needed.

.. _sec:fragment file:

Fragment library files
----------------------

The ``library_setup.py`` script will not create a fragment library since
the MCF lists zero fragments.

.. _sec:solid_xyz:

Configuration xyz file
----------------------

A simulation is initiated from an xyz file using the ``read_config`` or
``add_to_config`` options. :ref:`sec:start_type` details the ``read_config`` and
``add_to_config`` options.
