.. highlight:: none

.. _ch:input_files:

Files Required to Run Cassandra
===============================

.. _sec:input_file:

Simulation Input File
---------------------

This is a required file that is given as an argument to the Cassandra
executable. Example input files for each ensemble are provided in the
Examples directory that can be modified for new simulations. The input
file is divided into sections. Each section begins with a section header
that starts with a ``#``, e.g. ``# Run_Name``, and ends with a blank line.
Section ``# Move_Probability_Info`` is an exception and terminates with
``# Done_Move_Probability_Info``, because subsections e.g. ``#
Prob_Translation`` are separated by blank lines. Comment lines begin with
``!`` and are ignored. Sections in the input file can be listed in any
order, but the order and format of keywords and parameters given in each
section are important unless otherwise noted below. Previously, some
keywords were capitalized, e.g. ``CUBIC``, some contained an initial
capital, e.g. ``Units``, and some were all lowercase, e.g. ``kappa_ins``. New
in version 1.2, all keywords are supported in lowercase text; each word
in a section header must still begin with an initial capital.

.. _sec:run_name:

Run Name
~~~~~~~~

| ``# Run_Name``
| *Character*

The run name is specified on the next line following the keyword. This
name is used as a prefix for all the files produced by the simulation.
For example::

    # Run_Name
    dee.out

| Cassandra will then use ``dee.out`` as prefix for all output files
  created.

Simulation Type
~~~~~~~~~~~~~~~

| ``# Sim_Type``
| *Character*

Sets the ensemble (and thus the suite of moves) of a Cassandra
simulation. The following ensembles are supported:

-  ``nvt`` or ``nvt_mc`` (canonical ensemble)
-  ``npt`` or ``npt_mc`` (isothermal-isobaric ensemble)
-  ``gcmc`` (grand canonical ensemble)
-  ``gemc`` (Gibbs ensemble)
-  ``gemc_npt`` (Multi-species Gibbs ensemble)
-  ``nvt_min`` (canonical ensemble, only moves which lower the energy are
   accepted)
-  ``fragment`` or ``nvt_mc_fragment`` (canonical ensemble simulation of a
   fragment)
-  ``ring_fragment`` or ``nvt_mc_ring_fragment`` (canonical ensemble
   simulation of a ring fragment)

For example,::

    # Sim_Type
    npt

will run a Monte Carlo simulation in the isothermal-isobaric ensemble in which
the number of molecules of each species :math:`N`, the pressure :math:`P` and
temperature :math:`T` are held constant.

.. note::
    Simulation types ``fragment`` and ``ring_fragment`` are used only for
    generating a fragment library. 

Number of species
~~~~~~~~~~~~~~~~~

| ``# Nbr_Species``
| *Integer*

Total number of species in the simulation. For ionic systems, each ion
is counted as a separate species. For example, for a mixture of two
species, use the following::

    # Nbr_Species
    2

.. _sec:vdw_style:

VDW Style
~~~~~~~~~

| ``# VDW_Style``
| *Character(i,1)* [*Character(i,2) Real(i,3) Real(i,4)/Logical(i,4)*]

This keyword specifies the functional form of repulsion dispersion
interactions to be used and if tail corrections are added for box
:math:`i`. One line is required for each box. *Character(i,1)*
specifies the van der Waals model and can be ``lj`` for a
Lennard-Jones 12-6 potential, ``mie`` for a Mie potential, or ``none``
to turn off all repulsion-dispersion interactions. *Character(i,2)*
and *Real(i,3)* are required for ``lj`` or ``mie``. *Character(i,2)*
specifies how the Lennard-Jones potential is truncated. Options are
``cut``, ``cut_tail``, ``cut_switch``, or ``cut_shift``. Refer to
Chapter [Chapter:Force Field] for the functional forms. The other
parameters *Real(i,3)* and *Real(i,4)/Logical(i,4)* depend on the
selection of *Character(i,2)* as described below:

-  | ``cut``: This option cuts the potential at the distance specified by
     *Real(i,3)*. The fourth parameter is omitted. For example, to simulate one
     box with a 14 Å cutoff specify the following:

    .. code-block:: none
        
        # VDW_Style
        lj cut 14.0

   | Similarly, for a two box simulations such as used in the Gibbs
     ensemble where both boxes have a 14 Å cutoff, use the following:
    
    .. code-block:: none

        # VDW_Style
        lj cut 14.0
        lj cut 14.0

-  | ``cut_tail``: This options cuts the potential off at a distance
     corresponding to *Real(i,3)* and applies analytic tail corrections
     to the energy and pressure. An optional fourth argument
     *Logical(i,4)* can be set to ``true``, in which case *Real(i,3)* is
     ignored and the cutoff distance is always set to half of the
     simulation box length. The cutoff will change during the course of
     the simulation when attempting volume moves. This option is
     provided to enable reproduction of literature simulations that use
     a cut off distance of half the simulation box length, but its use
     is discouraged.

   | For example, to simulate one box with a 14 Å cutoff using tail
     corrections, specify the following:

    .. code-block:: none

        # VDW_Style
        lj cut_tail 14.0

   | For a two box simulation where the first box has a 14 Å cutoff and
     the second one has a 20 Å cutoff, use the following:

    .. code-block:: none

        # VDW_Style
        lj cut_tail 14.0
        lj cut_tail 20.0

-  | ``cut_switch``: This option cuts the potential off and smoothly
     brings the potential to zero using a spline. The potential is
     cutoff and the spline turned on at a distance specified by
     *Real(i,3)* (:math:`r_{on}` in Eq [Eq:cut\_switch]) and the
     potential goes to zero at a distance specified by *Real(i,4)*
     (:math:`r_{off}` in Eq [Eq:cut\_switch]).

   | A one box simulation using the ``cut_switch`` option would
     be specified as follows:

    .. code-block:: none

        # VDW_Style
        lj cut_switch 12.0 14.0

   | In this case, the Lennard-Jones potential would end at 12.0 Å and
     be smoothly taken to zero at 14.0 Å. :math:`r_{on} < r_{off}` or
     *Real(i,3)* :math:`<` *Real(i,4)*.

-  | ``cut_shift``: This option cuts the potential off at a distance
     specified by *Real(i,3)* and shifts the entire potential so that at
     this distance the potential is zero. The fourth parameter
     *Real(i,4)/Logical(i,4)* is omitted. The functional form of this
     potential is given in eq [Eq:cut\_shift].

   | To perform a two box simulation with a ``cut_shift`` option in which
     both boxes have a 10.5 Å cutoff, use the following:

    .. code-block:: none

        # VDW_Style
        lj cut_shift 10.5
        lj cut_shift 10.5

.. note:: 
    For all options, cutoff distances must be less than or equal to
    the shortest edge length of a simulation box.

Charge Style
~~~~~~~~~~~~

| ``# Charge_Style``
| *Character(i,1)* [*Character(i,2) Real(i,3) Real(i,4)*]

Cassandra allows the use of fixed partial charges on atomic centers
using a Coulomb potential of the form given in Eq [Eq:Coulomb]. If
this section is missing from the input file, the electrostatic energy
of the simulation will not be computed. If you do not wish to use a
Coulomb potential for box *i*, set *Character(i,1)* to ``none``. If
``none`` is selected for *Character(i,1)* then *Character(i,2)*,
*Real(i,3)* and *Real(i,4)* are omitted.

For example,

.. code-block:: none

    # Charge_Style
    none

should be used if you have no partial charges and are simulating a
single box (or the section can just be omitted).

To compute the electrostatic energy for box *i*, this section must be
included and *Character(i,1)* set to ``coul``. For this option,
*Character(i,2)* can be set to ``ewald`` if you want to use an Ewald
sum to compute Coulombic interactions, ``dsf`` if you want to use the
`Damped Shifted Force method <https://doi.org/10.1063/1.2206581>`_
by Fennell *et al.*, or it can be set to ``cut``,
in which case the Coulombic interactions will be cut off and the long
range interactions ignored. For the Ewald option, *Real(i,3)* is the
real space cutoff distance and *Real(i,4)* specifies the accuracy of
the Ewald summation. A reasonable value for the accuracy is
:math:`10^{-5}`. Note that the number of reciprocal vectors for the
Ewald summation is determined in the code based on the accuracy
parameter. For more details, see the
`paper by Fincham <https://doi.org/10.1080/08927029408022180>`_.

For example,

.. code-block:: none

    # Charge_Style
    coul ewald 12.0 1E-5

will use the Ewald sum for a single box. The real space cutoff will be
12 Å and the accuracy will be :math:`10^{-5}`. If you have two boxes,
like in a Gibbs ensemble calculation, then you could use the
following:

.. code-block:: none

    # Charge_Style
    coul ewald 12.0 1E-5
    coul ewald 30.0 1E-5

This will use an Ewald sum for both boxes. In the first box, the real
space cutoff will be 12 Å while in the second box a larger cutoff of
30 Å will be used.

.. note::
    When performing Gibbs ensemble simulations of vapor-liquid equilibria, the
    vapor box is often much larger than the liquid box. In this case, you will
    want to use a longer real space cutoff for the larger vapor box to avoid
    using too many reciprocal space vectors.

.. note::
    Also note that the real space cutoffs must always be less than or equal to
    half of the shortest edge length of a simulation box.

If you wish to use the Damped Shifted Force method, the entry
*Real(i,3)* is the electrostatic energy cutoff distance and
*Real(i,4)* is an optional entry to specify the damping parameter. If
not specified, Cassandra will set this value algorithmically from the
cutoff radius. For example,

.. code-block:: none

    # Charge_Style
    coul dsf 12.0 0.20

will use the Damped Shifted Force method for a single box. The
electrostatic energy cutoff will be set to 12 Å and the damping
parameter will be set to 0.20, which is a reasonable value for typical
liquid phase simulations.

.. note::

    If the cutoff in ``VDW_Style`` is set to half of the simulation box length,
    any cutoff distance specified in the ``Charge_Style`` section will default to
    the half of the simulation box length. In the case of Ewald summation,
    however, the accuracy will be the same as *Real(i,4)*.

Mixing Rule
~~~~~~~~~~~

| ``# Mixing_Rule``
| *Character*

Sets the method by which van der Waals interactions between unlike atoms are
calculated. Acceptable options are ``lb`` for Lorentz-Berthelot, ``geometric``
for geometric mixing rule and ``custom`` for allowing the user to provide
specific values. To use either ``lb`` or ``geometric`` keywords with the Mie
potential, all atomtypes must have the same repulsive and dispersive exponents.
If this section is missing, ``lb`` is used as default.

To illustrate the use of the ``custom`` option, consider a mixture of methane
(species 1) and butane (species 2) united atom models using a Lennard-Jones
potential. Methane has a single atomtype, CH4. Butane has two atomtypes:
pseudoatoms 1 and 4 are type CH3, pseudoatoms 2 and 3 are type CH2. The cross
interaction table is as follows:

.. code-block:: none

    # Mixing_Rule
    custom
    CH4 CH3 120.49 3.75
    CH4 CH2 82.51 3.83
    CH3 CH2 67.14 3.85

The order in which atom types are listed is unimportant, but the atom
types must match exactly the types given in each MCF. The
Lennard-Jones potential requires two parameters: an energy parameter
with units K, and a collision diameter with units Å. The Mie potential
requires four parameters: an energy parameter with units K, a
collision diameter with units Å, a repulsive exponent, and a
dispersive exponent.

.. _sec:seeds:

Starting Seed
~~~~~~~~~~~~~

| ``# Seed_Info``
| *Integer(1) Integer(2)*

Inputs for the starting random number seeds for the simulation.  Cassandra uses
a random number generator
`proposed by L’Ecuyer <https://doi.org/10.1090/S0025-5718-99-01039-X>`_,
which takes five seeds to calculate a random number, out of which
three are defined internally while two *Integer(1)*
and *Integer(2)* are supplied by the user.

As an example,

.. code-block:: none

    # Seed_Info
    1244432 8263662

.. note::
    Note that two independent simulations can be run using the same input
    information if different seeds are used. If two simulations having exactly
    the same input information and the same seeds are run, the results will be
    identical.

.. note::
    When a ‘checkpoint’ file is used to restart a simulation (see ``# Start_Type``
    below), the user supplied seeds will be overwritten by those present in the
    checkpoint file. If ``# Start_Type`` is ``read_config``, then the seeds
    specified in the input file are used.


Minimum Cutoff
~~~~~~~~~~~~~~

| ``# Rcutoff_Low``
| *Real*

Sets the minimum allowable distance in Å between two atoms. Any MC move
bringing two sites closer than this distance will be immediately rejected. It
avoids numerical problems associated with random moves that happen to place
atoms very close to one another such that they will have unphysically strong
repulsion or attraction. This distance must be less than the intramolecular
distance of all atoms in a species which are not bonded to one another. For
models that use dummy sites without explicitly defining bonds between dummy and
atomic sites of the molecules (for example, the TIP4P water model), it is
important that the minimum distance is set to be less than the shortest
distance between any two sites on the molecule. For most systems, 1 Å seems to
work OK, but for models with dummy sites, a shorter value may be required.

Pair Energy Storage
~~~~~~~~~~~~~~~~~~~

| ``# Pair_Energy``
| *Logical*

Cassandra can use a time saving feature in which the energies between
molecules are stored and used during energy evaluations after a move,
thereby saving a loop over all molecules. This requires more memory,
but it can be faster. The default is to not use this feature. If you
wish to use this, set *Logical* to ``true``.

Molecule Files
~~~~~~~~~~~~~~

| ``# Molecule_Files``
| *Character(i,1) Integer(i,2)*

This specifies the name of the molecular connectivity file (MCF) and
the maximum total number of molecules of a given species specified by
this MCF. A separate line is required for each species present in the
simulation. *Character(i,1)* is the name of the MCF for species *i*.
*Integer(i,2)* is the maximum number of molecules expected for the
species.

For example,

.. code-block:: none

    # Molecule_Files 
    butane.mcf 100
    hexane.mcf 20
    octane.mcf 5

specifies that there are three different species, and the MCFs state
the names of the files where information on the three species can be
found. Species 1 is butane, species 2 is hexane and species 3 is
octane. There can be a maximum of 100 butane molecules, 20 hexane
molecules and 5 octane molecules in the total system. The maximum
number of molecules specified here will be used to allocate memory for
each species, so do not use larger numbers than are needed.

Simulation Box
~~~~~~~~~~~~~~

| ``# Box\_Info``
| *Integer(1)*
| *Character(i)*
| *Real(i,1)* [*Real(i,2) Real(i,3)*]

This section sets parameters for the simulation boxes. *Integer(1)*
specifies the total number of boxes in the simulation. Gibbs ensemble
simulations must have two boxes. *Character(i)* is the shape of the
:math:`i`\ th simulation box. The supported keywords are ``cubic``,
``orthogonal``, and ``cell_matrix``.

If *Character(i)* is ``cubic``, *Real(i,1)* is the length of the box
edges in Å. Information for additional boxes is provided in an
analogous fashion and is separated from the previous box by a blank
line. For a two box simulation, box information is given as:

.. code-block:: none

    # Box_Info 2
    cubic
    30.0

    cubic
    60.0

This will construct a 30 x 30 x 30 Å cube and the second a 60 x 60 x
60 Å cube.

The options ``orthogonal`` and ``cell_matrix`` are only supported for
constant volume simulations (i.e. NVT or GCMC) which only have 1 box.
If *Character(1)* is ``orthogonal``, *Real(1,1) Real(1,2) Real(1,3)*
are the length, width and height that define the simulation box. For
example,

.. code-block:: none

    # Box_Info 1
    orthogonal
    30.0 35.0 40.0

This will create a simulation box with dimensions 30.0 x 35.0 x 40.0
Å.

A non-orthogonal box is created by setting *Character(1)* to
``cell_matrix``. In this case, three basis vectors are needed to
define the simulation box. Each vector is entered as a column of a 3x3
matrix. For example,

.. code-block:: none

    # Box_Info 1
    cell_matrix
    30  0  0
    0  35  0
    0   2 40

defines a simulation box with basis vectors (30, 0, 0), (0, 35, 2) and
(0, 0, 40).

Temperature
~~~~~~~~~~~

| ``# Temperature_Info``
| *Real(i)*

*Real(i)* is the temperature in Kelvin for box :math:`i`. For GEMC,
the temperature of box 2 will be read from a second line:

.. code-block:: none

    # Temperature_Info
    300.0
    300.0

Pressure
~~~~~~~~

| ``# Pressure_Info``
| *Real(i)*

*Real(i)* is the pressure setpoint in bar for box :math:`i`. For GEMC,
the pressure of box 2 will be read from a second line:

.. code-block:: none

    # Pressure_Info
    1.0
    1.0

If the simulation type does not require an input pressure (e.g., NVT),
this section will be ignored.

Chemical Potential
~~~~~~~~~~~~~~~~~~

| ``# Chemical\_Potential\_Info``
| *Real(1) ... Real(n)*

where *n* is the number of insertable species and *Real(i)* is the
chemical potential setpoint (shifted by a species-specific constant)
of insertable species *i* in kJ/mol. Each chemical potential will be
assigned in the order species appear in the ``Molecule_Files``
section. For species with insertion method none, the chemical
potential can be listed as none or omitted. This section is only read
for grand canonical simulations. See Eq. ([eq:muShift]) for more
information. For example, the adsorption of methane (species 2) in a
zeolite (species 1) can be computed by inserting methane molecules
into a box with a zeolite crystal. In this example, only one chemical
potential (for methane) is required and the following are equivalent:

.. code-block:: none

    # Chemical_Potential_Info
    -35.0

.. code-block:: none

    # Chemical_Potential_Info
    none -35.0

.. warning::

    Specifying the chemical potential as ``0.0`` is **not** the same as
    ``none``. 

Move Probabilities
~~~~~~~~~~~~~~~~~~

| ``# Move_Probability_Info``
| ``[subsections]``
| ``# Done_Probability_Info``

This section specifies the probabilities associated with different
types of MC moves to be performed during the simulation. The section
begins with the header ``# Move_Probability_Info`` and is terminated by
the footer ``# Done_Probability_Info``. All the move probability
subsections must be between the section header and footer.

.. note::

    If the move probabilities do not sum to 1.0, then the probability of
    each move will be divided by the total.

Translation
^^^^^^^^^^^

| ``# Prob_Translation``
| *Real(1)*
| *Real(i,1) ... Real(i,n)* \*One line required for each box :math:`i`

where :math:`n` is the number of species. *Real(1)* is the probability
of performing a center of mass translation move. *Real(i,j)* is the
maximum displacement in Å of species :math:`j` in box :math:`i`. This
subsection is optional in all ensembles.

For example, if you have three species and two boxes, you could
specify the translation probability as:

.. code-block:: none

    # Prob_Translation
    0.25
    2.0 2.5 1.0
    12.0 12.0 12.0

This will tell Cassandra to attempt center of mass translations 25% of
the total moves. For box 1, the maximum displacement will be 2.0 Å for
species 1, 2.5 Å for species 2, and 1.0 Å for species 3. For box 2,
the maximum displacement for all species is 12.0 Å.
For a simulation that involves solid frameworks, set the maximum
displacement of the solid species to zero. Every molecule in the
simulation with a maximum displacement greater than zero has an equal
chance of being moved.

Rotation
^^^^^^^^

| ``# Prob_Rotation``
| *Real(1)*
| *Real(i,1) ... Real(i,n)* \*One line required for each box :math:`i`

where :math:`n` is the number of species. The probability of performing a
rotation move is specified by *Real(1)* while *Real(i,j)* denotes the maximum
rotation for species :math:`j` in box :math:`i` in degrees about the x, y or
z-axis. The axis will be chosen with uniform probability. This subsection is
optional for all ensembles.

For example, if you are simulating a single species in two boxes, you could
specify the rotational probability as:

.. code-block:: none

    # Prob_Rotation
    0.25
    30.0 180.0

Twenty-five percent of the attempted moves will be rotations.  Molecules in box
1 will be rotated a maximum of 30 around the x, y, or z-axis. Molecules in box
2 will be rotated a maximum of 180 around the x, y, or z-axis.

If all species are point particles (such as single-site Lennard-Jones
particles), this section should be omitted. For a multi-species system, set
*Real(i,j)* to zero for point particles and solid frameworks.

Linear molecules are a special case. A molecule is identified as
linear if all angles in the MCF are fixed at 180. If a linear molecule
were aligned with the axis of rotation, then the molecular orientation
would not be changed. Therefore, linear molecules are rotated by
choosing a random unit vector with uniform probability without regard
to the molecule’s current orientation or the maximum rotation. As with
non-linear molecules, if *Real(i,j)* is zero, no molecules of species
:math:`j` will be rotated.

For a single box simulation of a non-linear molecule (species 1), a
linear molecule (species 2), and a point particle (species 3), you
could specify:

.. code-block:: none

    # Prob_Rotation
    0.25
    30.0 10.0 0.0

Molecules of species 1 will be rotated a maximum of 30 around the x, y
or z-axis, molecules of species 2 will be rotated by choosing a random
unit vector, and the point particles will not be rotated.

Angle
^^^^^

| ``# Prob\_Angle``
| *Real(1)*

A molecule will be selected at random and its angle will be perturbed based on
its Boltzmann weighted distribution. The probability of attempting this move is
the only required input. It is specified by *Real(1)*. 

For example,

.. code-block:: none

    # Prob_Angle
    0.3 

tells Cassandra to attempt angle moves 30% of the total moves for all molecules
containing angles within a given box.

.. note:: 

    Note that this move is rarely needed since the fragment
    libraries should already provide efficient sampling of angles. This
    move, however, may improve sampling of angles for large molecules in
    the case where parts of its fragments are rarely regrown by a regrowth
    move.


Dihedral
^^^^^^^^

| ``# Prob_Dihedral``
| *Real(1)*
| *Real(1) ... Real(n)*

The probability of performing a dihedral move is specified by
*Real(1)* while *Real(n)* denotes the maximum width of a dihedral
angle displacement for each species. The maximum width is given in
degrees. 

For example,

.. code-block:: none

    # Prob_Dihedral
    0.3
    20 0.0

tells Cassandra to attempt dihedral moves 30% of the total moves for all
molecules containing dihedrals within a given box. The maximum dihedral width
will be 20 for species 1 and 0.0 for species 2.  Since the maximum dihedral
width of species 2 is set to 0.0 in both boxes, no dihedral moves will be
attempted on species 2. Note that a single max dihedral width is provided, even
if species 1 may contain many dihedrals. This is also true for simulations with
more than one box. Also note that the same max dihedral width is used for
systems containing more than one box.

.. note::
    Note that this move is rarely needed since the regrowth moves
    should already provide efficient sampling of dihedrals. This move,
    however, may improve sampling of dihedrals for large molecules in the
    case where the parts of its fragments are rarely regrown (albeit a
    small maximum width is provided).

Regrowth
^^^^^^^^

| ``# Prob_Regrowth``
| *Real(1)*
| *Real(2,1) ... Real(2,n)*

where :math:`n` is the number of species. A regrowth move consists of deleting
part of the molecule randomly and then regrowing the deleted part via
configurational bias algorithm. This can result in relatively substantial
conformational changes for the molecule, but the cost of this move is higher
than that of a simple translation or rotation. The probability of attempting a
regrowth move is specified by *Real(1)* while *Real(2,i)* specifies the
relative probability of performing this move on species :math:`i`. The relative
probabilities must sum to 1 otherwise Cassandra will quit with an error. This
subsection is optional for all ensembles.

For example, if simulating 70 molecules of species 1 and 30 molecules of
species 2, you could specify the following:

.. code-block:: none

    # Prob_Regrowth
    0.3
    0.7 0.3

Thirty percent of the attempted moves will be regrowth moves. Seventy percent
of the regrowth moves will be attempted on a molecule of species 1 and the
balance of regrowth moves on a molecule of species 2.

.. note::
 
    *Real(2,i)* should be set to zero for monatomic, linear, or rigid
    species, including solid frameworks.

Volume
^^^^^^

| ``# Prob_Volume``
| *Real(1)*
| *Real(2)*
| [\ *Real(3)*]

*Real(1)* is the relative probability of attempting a box volume
change. Since volume changes are computationally expensive, this
probability should normally not exceed 0.05 and values from 0.01-0.03
are typical. *Real(2)* is the maximum volume displacement in
Å\ :sup:`3` for box 1. *Real(3)* is the maximum volume displacement
in Å\ :sup:`3` for box 2, and is only required for GEMC-NPT
simulations. The attempted change in box volume is selected from a
uniform distribution. This subsection is required for NPT, GEMC-NPT
and GEMC-NVT simulations.

For example, if you are simulating a liquid with a single box in the NPT
ensemble, the following:

.. code-block:: none

    # Prob_Volume
    0.02
    300

tells Cassandra to attempt volume moves 2% of the total moves. The box volume
would be changed by random amounts ranging from -300 Å\ :sup:`3` to +300 Å\
:sup:`3`. For a liquid box 20 Å per side, this would result in a maximum box
edge length change of about 0.25 Å, which is a reasonable value. Larger volume
changes should be used for vapor boxes. If you wish to perform a GEMC-NPT
simulation, you might specify the following:

.. code-block:: none

    # Prob_Volume
    0.02
    300
    5000

This tells Cassandra to attempt volume moves 2% of the total moves. The first
box volume (assumed here to be smaller and of higher density, such as would
occur if it were the liquid box) would be changed by random amounts ranging
from -300 Å\ :math:`^3` to +300 Å\ :math:`^3`. The second box volume would be
changed by random amounts ranging from -5000 Å\ :math:`^3` to +5000 Å\
:math:`^3`. As with all move probabilities, you can experiment with making
larger or smaller moves. Note that if the ``# Run_Type`` is ``equilibration``,
Cassandra will attempt to optimize the magnitude of the volume change to
achieve about 50% acceptance rates.

.. note::

    The volume perturbation move is only supported for cubic boxes.

Insertion and Deletion Moves
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

| ``# Prob_Insertion``
| *Real(1)*
| *Character(2,1) ... Character(2,n)*

where :math:`n` is the number of species. *Real(1)* sets the probability of
attempting insetion moves. *Character(2,i)* is the insertion method and can be
either ``cbmc`` or ``none``. If ``cbmc``, species :math:`i` will be inserted by
assembling its fragments using configurational bias Monte Carlo. If ``none``,
species :math:`i` will not be inserted or deleted. This subsection is required
for GCMC simulations.

If there is more than one insertable species, each is chosen for an insertion
attempt with equal probability. For example, if you are performing a GCMC
simulation with two species that can be inserted, you might specify the
following:

.. code-block:: none

    # Prob_Insertion
    0.1
    cbmc cbmc

This tells Cassandra to attempt insertions 10% of the total moves
and both species will be inserted using CBMC. If only species 1 is to
be inserted or deleted, use:

.. code-block:: none

    # Prob_Insertion
    0.1
    cbmc none


| ``# Prob_Deletion``
| *Real(1)*

*Real(1)* is the probability of attempting to delete a molecule during a
simulation, and must match the insertion probability to satisfy microscopic
reversibility. The molecule to delete is selected by first choosing a species
with uniform probability, and then choosing a molecule of that species with
uniform probability. If a species has the insertion method ``none``, no attempt
is made to delete it. This subsection is required for GCMC simulations.

| ``# Prob_Swap``
| *Real(1)*
| *Character(2,1) ... Character(2,n)*
| [\ *prob\_swap\_species Real(3,1) ... Real(3,n)*]
| [\ *prob\_swap\_from\_box Real(4,1) ... Real(4,i)*]

where :math:`n` is the number of species and :math:`i` is the number of boxes.
*Real(1)* is the probability of attempting to transfer a molecule from one box
to another. Similar to the ``# Prob_Insertion`` subsection, *Character(2,i)* is
the insertion method and can be ``cbmc`` or ``none``. If ``cbmc``, species
:math:`i` will be inserted by assembling its fragments using configurational
bias Monte Carlo. If ``none``, species :math:`i` will not be transferred
between boxes.  This subsection is required for a GEMC simulation.


For example, while performing a GEMC simulation for three species the first two
of which are exchanged while the third is not, specify the following:

.. code-block:: none

    # Prob_Swap
    0.1
    cbmc cbmc none

This tells Cassandra to attempt swap moves 10% of the total moves. Attempts
will be made to transfer species 1 and 2 between available boxes while
molecules of species 3 will remain in the boxes they are present in at the
start of the simulation.

By default, a molecule is chosen for the attempted swap with uniform
probability (amongst swappable molecules). As a result, if one species has a
much higher mole fraction in the system (e.g. if calculating methane solubility
in liquid water), then most attempted swaps will be of the more abundant
species. This behavior can be changed by using the optional keywords
``prob_swap_species`` and ``prob_swap_from_box``.

The keyword ``prob_swap_species`` must be given with :math:`n` options:
*Real(3,j)* is the probability of selecting species :math:`j`. The keyword
prob\_swap\_from\_box must be given with :math:`i` options: *Real(4,j)* is the
probability of selecting a molecule from box :math:`j`. For example, to select
a molecule of species 1 for 90% of attempted swaps and to select box 2 as the
donor box for 75% of attempted swaps, use:

.. code-block:: none

    # Prob_Swap
    0.1
    cbmc cbmc none
    prob_swap_species 0.9 0.1 0.0
    prob_swap_from_box 0.25 0.75

The probability of selecting a species with insertion method ``none`` must be 0.

Ring Flip Move
^^^^^^^^^^^^^^
| ``# Prob\_Ring``
| *Real(1) Real(2)*

This subsection is used when flip moves are to be attempted to sample bond
angles and dihedral angles in a ring fragment. For more details on this move
see `Shah and Maginn <https://doi.org/10.1063/1.3644939>`_.
The relative probability of attempting
a flip move is specified by *Real(1)* while the maximum angular displacement in
degrees for the move is given by *Real(2)*. For example, if the flip is to be
attempted 30% of the time and the maximum angular displacement for the move is
20 specify the following:

.. code-block:: none

    # Prob_Ring
    0.30 20.0

.. note::

    Note that this subsection is used only in input files that generate
    configuration libraries of ring moieties. The input file of the actual
    simulation would involve the ``# Prob_Regrowth`` keyword.

.. _sec:start_type:

Start Type
~~~~~~~~~~

| ``# Start_Type``
| *Character(1) [options]*

This section specifies whether Cassandra generates an initial
configuration or uses a previously generated configuration to start a
simulation. *Character(1)* can be one of four keywords: ``make_config``
``read_config``, ``add_to_config``, or ``checkpoint``. The specifications
for *[options]* depends on *Character(1)* as described below:


-  | ``make_config`` will generate an initial configuration using a
     configurational biased scheme. The number of molecules of each
     species is specified as follows:

   | ``make_config`` *Integer(1) ... Integer(n)*
   | where *n* is the number of species and *Integer(i)* is the number
     of molecules of species :math:`i` to insert into the box. This
     keyword can be repeated for each box. For example, to generate an
     initial configuration with 100 molecules of species 1 and 75
     molecules of species 2:

     .. code-block:: none

        # Start_Type
        make_config 100 75

   | If the simulation also has a second box with 25 molecules of
     species 2 only:

     .. code-block:: none

        # Start_Type
        make_config 100 75
        make_config   0 25

-  | ``read_config`` will use the coordinates from a ``.xyz`` file. For
     example, a configuration generated at one temperature may be used
     to initiate a simulation at another temperature. After ``read_config``,
     the number of molecules of each species must be given, followed by
     the ``.xyz`` filename:

   | ``read_config`` *Integer(1) ... Integer(n) Character(1)*
   | where *n* is the number of species, *Integer(i)* is the number of
     molecules of species :math:`i` to read from file *Character(1)*.
     This keyword can be repeated for each box. For example, to start a
     simulation using a configuration of 50 molecules each of species 1
     and 2:

     .. code-block:: none

        # Start_Type
        read_config 50 50 liquid.xyz

   | If the simulation also has a second box with 10 molecules of
     species 1 and 90 molecules of species 2:

     .. code-block:: none

        # Start_Type
        read_config 50 50 liquid.xyz
        read_config 10 90 vapor.xyz

   | The ``.xyz`` files must have the following format:

     .. code-block:: none

        <number of atoms>
        comment line
        <element> <x> <y> <z>
        ...

-  | ``add_to_config`` will read the coordinates from an .xyz file,
     but then insert additional molecules. After ``add_to_config`` specify
     the number of molecules of each species to be read, followed by the
     .xyz filename, followed by the number of molecules of each species
     to be added:

   | ``add_to_config`` *Integer(1) ... Integer(n) Character(1)
     Integer(n+1) ... Integer(2n)*
   | where *n* is the number of species, *Integer(1)* through
     *Integer(n)* are the number of molecules of each species to read
     from file *Character(1)*, and *Integer(n+1)* through *Integer(2n)*
     are the number of molecules of each species to add to the
     configuration. This keyword can be repeated for each box. For
     example, to start a simulation by reading in a zeolite (speces 1)
     configuration and adding 30 molecules of methane (species 2):

     .. code-block:: none

        # Start_Type
        add_to_config 1 0 MFI.xyz 0 30

   | where the file ``MFI.xyz`` contains the coordinates of a unit cell
     of MFI silicalite.

-  | ``checkpoint`` this keyword is used to restart a simulation from
     a checkpoint file. During the course of a simulation, Cassandra
     periodically generates a checkpoint file (``*.chk``) containing
     information about the total number of translation, rotation and
     volume moves along with the random number seeds and the coordinates
     of each molecule and its box number at the time the file is
     written. Cassandra provides the capability of restarting from this
     state point in the event that a simulation crashes or running a
     production simulation from an equilibrated configuration. For this
     purpose, in addition to the checkpoint keyword, additional
     information in the form of the name of the checkpoint file
     *Character(1)* is required in the following format:

   | ``checkpoint`` *Character(1)*
   | For example, to continue simulations from a checkpoint file
     ``methane_vle_T148.chk``, you might specify:

    .. code-block:: none

        # Start_Type
        checkpoint methane_vle_T148.chk

    .. note::

        Note that when a checkpoint file is used to restart a simulation,
        the seeds for random number generation supplied by the user will be
        overwritten by those present in the checkpoint file. By contrast,
        if ``# Start_Type`` is ``read_config``, then the seeds specified
        in the input file are used.

.. note::
    Unless starting from a checkpoint file, input files for a multi-box
    simulation must have one line for each box in the ``Start_Type``
    section. Each line can start with a different keyword. For example, a
    GEMC simulation of a water(1)-methane(2) mixture can begin from an
    equilibrated water box and a new vapor box:
    .. code-block::

        # Start_Type
        read_config 100 0 water.xyz
        make_config  50  50


Run Type
~~~~~~~~

| ``# Run_Type``
| *Character(1)* *Integer(1)* [*Integer(2)*]

This section is used to specify whether a given simulation is an equilibration
or a production run. For an equilibration run, the maximum translational,
rotational, torsional and volume widths (for an NPT or a GEMC simulation) are
adjusted to achieve 50% acceptance rates. During a production run, the maximum
displacement width for different moves are held constant.

Depending on the type of the simulation, *Character(1)* can be set to either
``equilibration`` or ``production``. For an ``equilibration`` run, *Integer(1)*
denotes the number of MC steps performed for a given thermal move before the
corresponding maximum displacement width is updated. *Integer(2)* is the number
of MC volume moves after which the volume displacement width is updated. This
number is optional if no volume moves are performed during a simulation (for
example in an NVT or a GCMC simulation). When the run type is set to
``production``, the MC moves refer to the frequency at which the acceptance
ratios for various moves will be computed and output to the log file. These
acceptance rates should be checked to make sure proper sampling is achieved.

For an NPT equilibration run in which the widths of the thermal move are to be
updated after 1000 MC moves and maximum volume displacements after 100 volume
moves, specify the following:

.. code-block:: none

    # Run_Type
    equilibration 100 10

For an NVT production run in which the acceptance ratios of various thermal
moves are printed to the log file after every 250 MC steps of a given thermal
move, use the following:

.. code-block:: none

    # Run_Type
    production 250

Simulation Length
~~~~~~~~~~~~~~~~~

| ``# Simulation_Length_Info``
| *units Character(1)*
| *prop\_freq Integer(2)*
| *coord\_freq Integer(3)*
| *run Integer(4)*
| [\ *steps\_per\_sweep Integer(5)*]
| [\ *block\_averages Integer(6)*]

This section specifies the frequency at which thermodynamic properties and
coordinates are output to a file. The ``units`` keyword determines the method
by which the simulation is terminated and data is output.  *Character(1)* can
be minutes, steps, or sweeps. Thermodynamic quantities are output every
*Integer(2)* units, coordinates are written to the disk every *Integer(3)*
units and the simulation will stop after *Integer(4)* units.

If *Character(1)* is minutes, then the simulation runs for a specified time.
For example, to run a simulation for 60 minutes with thermodynamic properties
written every minute and coordinates output every 10 minutes, use:

.. code-block:: none

    # Simulation_Length_Info
    units minutes
    prop_freq 1
    coord_freq 10
    run 60

If *Character(1)* is steps, the simulation runs for a specified number of MC
steps. An MC step is defined as a single MC move, regardless of type and
independent of system size. To run a simulation of 50,000 steps such that
thermodynamic quantities are printed every 100 MC steps and coordinates are
output every 10,000 steps, use:

.. code-block:: none

    # Simulation_Length_Info
    units steps
    prop_freq 100
    coord_freq 10000
    run 50000

If *Character(1)* is sweeps, the simulation runs for a specified number of MC
sweeps. The number of MC steps per sweep can be defined using the optional
keyword ``steps_per_sweep``. The default ``steps_per_sweep`` value is the sum
of the weights of each move type. A sweep is typically defined as the number of
MC moves needed for every move to be attempted for every molecule.

For example, in a water box of 100 molecules in the NPT ensemble, a sweep would
be 201 moves-100 translations, 100 rotations and 1 volume change. To run a
simulation of 1,000 sweeps with thermodynamic quantities are printed every 100
sweeps and coordinates are output every 100 sweeps, use the following:

.. code-block:: none

    # Simulation_Length_Info
    units sweeps
    prop_freq 100
    coord_freq 100
    run 1000
    steps_per_sweep 201

The optional keyword ``block_avg_freq`` switches the thermodynamic output from
instantaneous values to block average values, where *Integer(6)* is the number
of units per block. The number of blocks is given by *Integer(4)*/*Integer(6)*
and the number of data points per block is *Integer(6)*/*Integer(2)*. For
example, during a run of 1,000,000 steps, with properties computed every 100
steps and averaged every 100,000 steps, specify:

.. code-block:: none

    # Simulation_Length_Info
    units steps
    run 1000000
    block_avg_freq 100000
    prop_freq 100
    coord_freq 100

This simulation will output 10 averages, and each average will be computed from
1000 data points.

Property Output
~~~~~~~~~~~~~~~

| ``# Property_Info`` *Integer(i)*
| *Character(j) \*One line for each property :math:`j`*

This section provides information on the properties that are output.
More than one section is allowed for multiple boxes. In this case,
each section is separated by a blank line. *Integer(i)* is the
identity of the box for which the properties are desired.
*Character(i,j)* is the property that is to be output. Each property
is specified on a separate line. The supported keywords are:

* ``energy_total``: Total energy of the system (extensive) in kJ/mol
* ``energy_lj``: Lennard-Jones energy of the sytem in kJ/mol
* ``energy_elec``: Electrostatic energy of the sytem in kJ/mol
* ``energy_intra``: Total intramolecular energy of the system including
  bonded and non-bonded interactions in kJ/mol
* ``enthalpy``: Enthalpy of the system (extensive) kJ/mol. The enthalpy
  is computed using the pressure setpoint for isobaric simulations and
  the computed pressure for all other ensembles.
* ``pressure``: Pressure of the system in bar
* ``volume``: Volume of the system in Å\ :sup:`3`
* ``nmols``: Number of molecules of each species
* ``density``: Density of each species in #/Å\ :sup:`3`
* ``mass_density``: Density of the system in kg/m\ :sup:`3`

For example, if you would like total energy, volume and pressure of a one box
system to be written, you may specify the following:

.. code-block:: none

    # Property_Info 1
    energy_total
    volume
    pressure

For a GEMC-NVT simulation, total energy and density of all the species in box 1
and total energy, density of all the species in box 2 along with the pressure
may be output using the following format:

.. code-block:: none

    # Property_Info 1
    energy_total
    density

    # Property_Info 2
    energy_total
    density
    pressure

Fragment Files
~~~~~~~~~~~~~~

| ``# Fragment_Files``
| *Character(i)* *Integer(i)* \*One line for each fragment :math:`i`

In this section, information about the fragment library is specified.
*Character(i)* gives the location of the fragment library of fragment
:math:`i`; *Integer(i)* is the corresponding integer id specifying the
type of the fragment.

.. note::

    This section is automatically generated by ``library_setup.py``. However,
    if there is a need to change this section, follow the example given below.

For a simulation involving two species of which the first one contains three
distinct fragments and species 2 has two identical fragments, this section
might look like:

.. code-block:: none

    # Fragment_Files
    frag_1_1.dat 1
    frag_2_1.dat 2
    frag_3_1.dat 3
    frag_1_2.dat 4
    frag_1_2.dat 4

This tells Cassandra to use the files ``frag_1_1.dat``, ``frag_2_1.dat`` and
``frag_3_1.dat`` for the three fragments of species 1. Since species 2 has two
identical fragment, Cassandra will use the same fragment library ``frag_1_2.dat``
for these fragments.

Verbosity in log file
~~~~~~~~~~~~~~~~~~~~~

| ``# Verbose_Logfile``
| *Logical*

This optional section is used to control the level of detail about the
simulation setup that is output to the log file. Controlling this can be useful
for development purposes. If this section is missing, *Logical* is set to
``false`` by default. Supported options for *Logical* are ``true`` or
``false``.

File Info
~~~~~~~~~

| ``# File_Info``
| *Character*

This section is used only while generating a fragment library.  Cassandra will
use the filename specified in *Character* to store different conformations of
the fragment being simulated.

.. note::

    This section is automatically handled
    by ``library_setup.py``. However, if the user wishes to modify this part,
    use the following template:

    .. code-block:: none

        # File_Info
        frag.dat

    This tells Cassandra to store the fragment library in the file named ``frag.dat``.

CBMC parameters
~~~~~~~~~~~~~~~

| ``# CBMC_Info``
| ``kappa_ins`` *Integer(1)*
| ``kappa_dih`` *Integer(2)*
| ``rcut_cbmc`` *Real(3,1)* [*Real(3,2)*]

Cassandra utilizes a configurational bias methodology based on
`sampling a library of fragment conformations <https://doi.org/10.1063/1.3644939>`_.
This section sets a number of parameters required for biased insertion/deletion (refer to
the sections ``# Prob_Insertion``, ``# Prob_Deletion`` and ``# Prob_Swap`` and
configurational regrowth (``# Prob_Regrowth``).

This section is only required if molecules are regrown, inserted and/or
deleted.  Keyword ``kappa_ins`` is required if the section ``# Start_Type`` is
given with keyword ``make_config`` or ``add_to_config``, or if the section ``#
Sim_Type`` is ``gcmc``, ``gemc`` or ``gemc_npt``.

Keyword ``kappa_ins`` is required if section ``# CBMC_Info`` is required.
For a biased insertion, a fragment is chosen to insert first in proportion to
the number of atoms in fragment. For example, to insert a united-atom molecule
of ethylbenzene, the ring fragment has 7 pseudoatoms while the other has 3. The
ring fragment will be inserted first with a probability of 0.7.  By contrast,
to insert a united-atom molecule of dodecane, all ten fragments have 3
pseudoatoms and so one is chosen with uniform probability. After choosing a
Boltzmann-distributed conformation and an orientation with uniform probability,
*Integer(1)* trial positions are generated for the center-of-mass of the
fragment. One of the trial positions is then selected randomly based on the
Boltzmann weight of the energy of the trial position.

Keyword ``kappa_dih`` is required if any species composed of multiple
fragments is inserted/deleted or regrown. Additional fragments are
added to the growing molecule using *Integer(2)* trial dihedral angles
that connect the new fragment to the existing part of molecule.

Keyword ``rcut_cbmc`` is required if section ``# CBMC_Info`` is required.
For all the trials, energy of the partially grown molecule with itself
and surrounding molecules is to be calculated. For this purpose, a
short cutoff is used. *Real(4,i)* specifies the cutoff distance in
Å for box :math:`i`. A short cutoff is fast, but might miss some
overlaps. You can experiment with this value to optimize it for your
system.

For a GEMC simulation in which 12 candidate positions are generated
for biased insertion/deletion, 10 trials for biased dihedral angle
selection and the cutoff for biasing energy calculation is set to 5.0
Å in box 1 and 6.5 Å in box 2, this section would look like:

.. code-block:: none

    # CBMC_Info
    kappa_ins 12
    kappa_dih 10
    rcut_cbmc 5.0 6.5


.. _sec:mcf_file:

Molecular Connectivity File
---------------------------

A Molecular Connectivity File (MCF) defines the information related to
bonds, angles, dihedrals, impropers fragments and non bonded
interactions for a given species. One MCF is required for each species
present in the system. The information contained in this file involves
the force field parameters, atoms participating in each of the
interactions and the functional form used in each potential
contribution. The keywords are preceeded by a ``#`` and comments follow
a ``!``. Similarly to the input file, the order of the keywords is not
important. A complete list of the keywords is provided below.

.. note::

    MCFs are generated by the script ``mcfgen.py`` automatically.  The
    following description is provided for the users who wish to modify the MCF
    or manually write the MCF.


.. warning::

    Parameters for all of the following keywords must be separated by spaces
    only. Do not use the tab character.


Atom Info
~~~~~~~~~

| ``# Atom_Info``
| *Integer(1)*
| *Integer(2) Character(3)\*6 Character(4)\*2 Real(5) Real(6) Character(7)\*20 Optional\_Parms Character(fin)*

This keyword specifies the information for non-bonded interactions. It
is a required keyword in the MCF. If not specified, the code will
abort. The inputs are specified below:

-  *Integer(1)*: Total number of atoms in the species.

-  *Integer(2)*: Atom index.

-  *Character(3)\*6*: Atom type up to 6 characters. This string of
   characters should be unique for each interaction site in the system,
   i.e. do not use the same atom type for two atoms in the same (or
   different) species unless the (pseudo)atoms have the same atom types.

-  *Character(4)\*2*: Atom element name up to 2 characters.

-  *Real(5)*: Mass of the atom in amu. Note that for united atom models,
   this would be the mass of the entire pseudoatom.

-  *Real(6)*: Charge on the atom.

-  *Character(7)*: The functional form for van der Waals (vdW)
   interactions. Options are ``LJ`` for Lennard-Jones, ``Mie`` for the Mie
   potential, or ``NONE`` if the atom type does not have vdW interactions.
   ``LJ`` and ``Mie`` cannot be used in the same simulation. This must match
   what is given for ``# VDW_Style`` (:ref:`sec:vdw_style`) in
   the input file.

-  *Character(fin)*: The final entry on the line is ``ring`` only if the
   atom is part of a ring fragment. Note that a ring fragment is defined
   as those atoms that belong to the ring (e.g. in cyclohexane, all the
   six carbons) and any atom directly bonded to these ring atoms (e.g.
   in cyclohexane, all the hydrogens). In other words, all of the ring
   and exoring atoms are given the ring flag. For atoms that are not
   part of rings, leave this field blank.

Additional parameters are required for LJ and Mie potentials. For LJ,

-  *Real(8)*: The energy parameter in K.

-  *Real(9)*: Collision diameter (:math:`\sigma`) in Å.

For Mie,

-  *Real(8)*: The energy parameter in K.

-  *Real(9)*: Collision diameter (:math:`\sigma`) in Å.

-  *Real(10)*: The repulsive exponent.

-  *Real(11)*: The dispersive exponent.

.. note::

    For single-fragment species, the branch point atom
    is listed as the first atom.

For example, for a united atom pentane model:

.. code-block:: none

    # Atom_Info
    5
    1 CH3_s1 C 15.0107 0.0 LJ 98.0 3.75
    2 CH2_s1 C 14.0107 0.0 LJ 46.0 3.95
    3 CH2_s1 C 14.0107 0.0 LJ 46.0 3.95
    4 CH2_s1 C 14.0107 0.0 LJ 46.0 3.95
    5 CH3_s1 C 15.0107 0.0 LJ 98.0 3.75

The number below the keyword ``# Atom_Info`` specifies a species with
5 interaction sites, consistent with a united atom pentane model. The
first column specifies the atom ID of each of the pseudo atoms. The
second and third columns provide the atom type and atom name,
respectively. The fourth column represents the atomic mass of each
pseudoatom. Note that the mass of ``CH3_s1`` is 15.0107 for this united
atom model, as it involves a carbon and three hydrogen atoms. The same
applies for all other interaction sites. The fifth column contains the
partial charges placed on each of these pseudoatoms. The sixth,
seventh and eighth columns contain the repulsion-dispersion functional
form, the energy parameter and the collision diameter respectively. In
this case, the usual Lennard-Jones functional form is used. Note that
none of these atoms used the flag ``ring``, as no rings are present in
this molecule.

For a molecule containing rings, for example cyclohexane:

.. code-block:: none

    # Atom_Info
    6
    1 CH_s1 C 13.0107 0.0 LJ 52.5 3.91 ring
    2 CH_s1 C 13.0107 0.0 LJ 52.5 3.91 ring
    3 CH_s1 C 13.0107 0.0 LJ 52.5 3.91 ring
    4 CH_s1 C 13.0107 0.0 LJ 52.5 3.91 ring
    5 CH_s1 C 13.0107 0.0 LJ 52.5 3.91 ring
    6 CH_s1 C 13.0107 0.0 LJ 52.5 3.91 ring

.. note::

    The flag ‘ring’ was appended as the last column for each site in this
    cyclic molecule.

For the SPC/E water model:

.. code-block:: none

    # Atom_Info
    3
    1 O1_s1 O 16.00 -0.8476 LJ 78.20 3.1656
    2 H2_s1 H 1.000 0.4238 NONE
    3 H3_s1 H 1.000 0.4238 NONE

.. note::

    This is a molecule with a single fragment, so the branch point atom is
    the first atom in the list.

For a single-site model of CO2 using the Mie potential:

.. code-block:: none

    # Atom_Info
    1
    1 CO2 C 44.00 0.0 Mie 361.69 3.741 23.0 6.66

where the last two parameters are the repulsive and dispersive
exponents, respectively.

Bond Info
~~~~~~~~~

| ``# Bond_Info``
| *Integer(1)*
| *Integer(i,2) Integer(i,3) Integer(i,4) Character(i,5) Real(i,6) Real(i,7)*

This section provides information on the number of bonds in a molecule
and atoms involved in each bond along with its type. It is a required
keyword in the MCF. If not specified, the code will abort. The inputs
are specified below:

-  *Integer(1)*: Total number of bonds in the species. From the next
   line onwards, the bonds are listed sequentially and information for
   each bond is included on a separate line.

-  *Integer(i,2)*: Index of the :math:`i^{th}` bond.

-  *Integer(i,3) Integer(i,4)*: IDs of the atoms participating in the
   bond.

-  *Character(i,5)*: Type of the bond. At present only ‘fixed’ is
   permitted.

-  *Real(i,6)*: Specifies the bond length for a particular bond in Å.

.. note::

    At present, Cassandra simulations can be carried out only
    for fixed bond length systems.


For example, for the water model SPC/E, the ``# Bond_Info`` section is
the following:

.. code-block:: none

    # Bond_Info
    2
    1 1 2 fixed 1.0
    2 1 3 fixed 1.0

In the above example, two bonds are specified whose fixed length is
set to 1.0 Å.

Angle Info
~~~~~~~~~~

| ``# Angle_Info``
| *Integer(1)*
| *Integer(i,2) Integer(i,3) Integer(i,4) Integer(i,5) Character(i,6) Real(i,7) Real(i,8)*

The section lists the information on the angles in the species. It is
a required keyword in the MCF. If not specified, the code will abort.

-  *Integer(1)*: Number of angles in the species.

-  *Integer(i,2)*: Index of the :math:`i^{th}` angle.

-  *Integer(i,3) Integer(i,4) Integer(i,5)*: IDs of the atoms
   participating in the :math:`i^{th}` angle. Note that *Integer(i,4)*
   is the ID of the central atom.

-  *Character(i,6)*: Type of the angle. Currently, Cassandra supports
   ‘fixed’ and ‘harmonic’ (see :ref:`sec:ff_angles`) angles.
   For the ‘fixed’ option, *Real(i,7)* is the value of the angle and
   *Real(i,8)* is ignored by the code if specified. In the case of
   ’harmonic’ potential type, *Real(i,7)* specifies the harmonic force
   constant (:math:`K/rad^2`) while *Real(i,8)* is the nominal bond
   angle (in degrees).

For example, for a united atom pentane molecule with flexible angles,
this section is the following:

.. code-block:: none

    # Angle_Info
    3
    1 1 2 3 harmonic 31250.0 114.0
    2 2 3 4 harmonic 31250.0 114.0
    3 3 4 5 harmonic 31250.0 114.0

In the above example, the three angles between the pseudoatoms found
in the pentane model are specified. The three angles have an harmonic
potential, whose force constant is equal and is set to 31250.0
K/rad\ :sup:`2`. Finally, the equilibrium angle for these angles is
114.0°.

An example for SPC/E water model with fixed angles is:

.. code-block:: none

    # Angle_Info
    1
    1 2 1 3 fixed 109.47

This model has only one angle that is set to 109.47°.
No force constant is provided as the angle is fixed.

Dihedral Info
~~~~~~~~~~~~~

| ``# Dihderal_Info``
| *Integer(1)*
| *Integer(i,2) Integer(i,3) Integer(i,4) Integer(i,5) Integer(i,6) Character(i,7) Real(i,8) Real(i,9) Real(i,10) Real(i,11)*

This section of the MCF lists the number of dihedral angles and
associated information for a given species. It is a required keyword
in the MCF. If not specified, the code will abort.

-  *Integer(1)*: Lists the number of dihedral angles.

-  *Integer(i,2)*: Index of the :math:`i^{th}` dihedral angle.

-  *Integer(i,3): Integer(i,6)* - IDs of the atoms in the :math:`i^{th}`
   dihedral angle.

-  *Character(i,7)* : Dihedral potential type. Acceptable options are ``OPLS``,
   ``CHARMM``, ``harmonic`` and ``none``. If ``OPLS`` dihedral potential type is
   selected, then the real numbers *Real(i,8) - Real(i,11)* are the coefficients
   in the Fourier series (see :ref:`sec:ff_dihedrals`). The units are in kJ/mol. For
   the ``CHARMM`` dihedral potential type, three additional parameters are
   specified: :math:`a_0, a_1` and :math:`\delta` (see :ref:`sec:ff_dihedrals`). If
   ``harmonic`` dihedral potential type is used, then two additional parameters,
   :math:`K_{phi}` and :math:`\phi_0` (see :ref:`sec:ff_dihedrals`), are
   specified. For the ``none`` dihedral potential type, no additional parameters
   are necessary.

For example, for a united atom pentane molecule using an OPLS dihedral
potential type, the dihedrals are specified as follows:

.. code-block:: none

    # Dihedral_Info
    2
    1 1 2 3 4 OPLS 0.0 2.95188 -0.5670 6.5794
    2 2 3 4 5 OPLS 0.0 2.95188 -0.5670 6.5794

In this model two dihedral angles are specified by atoms 1,2,3,4 and
2,3,4,5. This model uses an OPLS functional form and thus four
parameters are provided after the OPLS flag.

Intramolecular Scaling
~~~~~~~~~~~~~~~~~~~~~~

| ``# Intra_Scaling``
| *Real(i,1) Real(i,2) Real(i,3) Real(i,4)*
| *Real(i,5) Real(i,6) Real(i,7) Real(i,8)*

This section sets the intramolecular scaling for 1-2, 1-3, 1-4 and 1-N
interactions within a given species. 1-2 means interactions between
atom 1 and another atom 2 directly bonded to it, 1-3 means
interactions between atom 1 and other atoms 3 separated from atom 1 by
exactly two bonds, etc. The first line corresponds to the VDW scaling:
*Real(i,1) Real(i,2) Real(i,3) Real(i,4)* apply to 1-2, 1-3, 1-4 and
1-N interactions, where N corresponds to all atoms separated from atom
1 by more than three bonds. The second line corresponds to the Coulomb
scaling: *Real(i,5) Real(i,6) Real(i,7) Real(i,8)* apply to 1-2, 1-3,
1-4 and 1-N interactions. Note that intramolecular scaling applies to
all the boxes in the simulation.

For example,

.. code-block:: none

    # Intra_Scaling
    0.0 0.0 0.5 1.0
    0.0 0.0 0.5 1.0

turns off 1-2 and 1-3 interactions, scales the VDW and
Coulombic interactions for 1-4 atoms by 50%, and uses full
interactions for all other atom pairs in the species.

.. note::

    If the ``# Intra_Scaling`` section is missing from the MCF, it will be
    looked for in the input file. If provided, the values in the MCF file
    will always override any values provided in the input file.

Fragment Info
~~~~~~~~~~~~~

| ``# Fragment_Info``
| *Integer(1)*
| *Integer(i,2) Integer(i,3) Integer(i,4) Integer(i,5) ...
  Integer(i,2+Integer(i,3))*

This section defines the total number of fragments in a given species.
It is an optional keyword. However, if the species is composed of
fragments, then this section must be specified. The inputs are
specified below:

-  *Integer(1)*: Total number of fragments.

-  *Integer(i,2)*: Index of the :math:`i^{th}` fragment.

-  *Integer(i,3)*: Number of atoms in the :math:`i^{th}` fragment.

-  *Integer(i,4) ... Integer(i,2+integer(i,3))*: List of the atom IDs in
   the fragment. The first atom ID is that for the branch point atom.
   .. warning::

        Atom ordering for the remaining atoms must match the order of atoms
        in the fragment library files.

For example, for a pentane united atom model:

.. code-block:: none

    # Fragment_Info
    3
    1 3 2 1 3
    2 3 3 2 4
    3 3 4 3 5

This specifies three fragments. Each of these fragments has three atoms. The
first atom specified for each of the fragments is the branch point atom.

Fragment Connectivity
~~~~~~~~~~~~~~~~~~~~~

| ``# Fragment_Connectivity``
| *Integer(1)*
| *Integer(i,2) Integer(i,3) Integer(i,4)*

The section lists the fragment connectivity - which fragment is bonded
to which other fragment. It is a required keyword if
``Fragment_Info`` is specified.

-  *Integer(1)*: total number of fragment connections.

-  *Integer(i,2)*: index of the :math:`i^{th}` fragment connectivity.

-  *Integer(i,3) Integer(i,4)*: fragment IDs participating in the
   connectivty.

For example, for a pentane united atom model:

.. code-block:: none

    # Fragment_Connectivity
    2
    1 1 2
    2 2 3

In this example, there are three fragments, therefore, two fragment
connectivities must be specified. Note that fragment 1 is connected to fragment
2 and fragment 2 is connected to fragment 3.
