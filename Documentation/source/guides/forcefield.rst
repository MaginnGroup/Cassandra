Force Field
===========

Bonds
-----

Cassandra is designed assuming all bond lengths are fixed. If you wish
to utilize a force field developed with flexible bond lengths, we
recommend that you either use the nominal or “equilibrium” bond lengths
of the force field as the fixed bond lengths specified for a Cassandra
simulation or carry out an energy minimization of the molecule with a
package that treats flexible bond lengths and utilize the bond lengths
obtained from the minimization.

.. table:: Cassandra units for bonds
    :widths: auto
    :align: center

    =========== ========= ========
    Parameter   Symbol    Units
    =========== ========= ========
    Bond length :math:`l` Å
    =========== ========= ========

.. _sec:ff_angles:

Angles
------

Cassandra supports two types of bond angles:

-  ``fixed`` : The angle declared as fixed is not perturbed during the
   course of the simulation.

-  ``harmonic`` : The bond angle energy is calculated as

   .. math::

      E_\theta = K_\theta (\theta - \theta_0)^2
      \label{Eq:angle_potential}

where the user must specify :math:`K_\theta` and :math:`\theta_0`.
Note that a factor of :math:`1/2` is **not used** in the energy
calculation of a bond angle. Make sure you know how the force
constant is defined in any force field you use.

.. table:: Cassandra units for bond angles
    :widths: auto
    :align: center

    ========================== =================  ==================
    Parameter                       Symbol            Units
    ========================== =================  ==================
    Nominal bond angle          :math:`\theta_0`   degrees
    Bond angle force constant   :math:`K_\theta`   K/rad\ :sup:`2`
    ========================== =================  ==================


.. _sec:ff_dihedrals:

Dihedrals
---------

Cassandra can handle four different types of dihedral angles:

-  ``OPLS``: The functional form of the dihedral potential is

   .. math::

      E_\phi = a_0 + a_1\, \left ( 1 + \cos(\phi)  \right ) + a_2 \, \left ( 1 -
        \cos(2\phi)\right ) + a_3 \, \left ( 1 + \cos (3\phi)\right )

where :math:`a_0`, :math:`a_1`, :math:`a_2` and :math:`a_3` are
specified by the user.

-  ``CHARMM``: The functional form of the potential is

   .. math::

      E_\phi = a_0  (1 + \cos (a_1\phi - \delta))

where :math:`a_0`, :math:`a_1` and :math:`\delta` are specified by
the user.

-  ``harmonic``: The dihedral potential is of the form:

   .. math::

      E_\phi = K_\phi (\phi - \phi_0)^2

where :math:`K_\phi` and :math:`\phi_0` are specified by the user.

-  ``none`` : There is no dihedral potential between the given atoms.

.. table:: Cassandra units for proper dihedrals
    :widths: auto
    :align: center

    =========================== ==================================================== ================ 
     Functional Form              Parameter                                               Units
    =========================== ==================================================== ================ 
     OPLS                        :math:`a_0`, :math:`a_1`, :math:`a_2`, :math:`a_3`   kJ/mol
     CHARMM                      :math:`a_0`                                          kJ/mol
     CHARMM                      :math:`a_1`                                          dimensionless
     CHARMM                      :math:`\delta`                                       degrees
     harmonic                    :math:`K_\phi`                                       K/rad\ :sup:`2`
     harmonic                    :math:`\phi_0`                                       degrees
    =========================== ==================================================== ================ 

Impropers
---------

Improper energy calculations can be carried out with the following two
options:

-  ``none``: The improper energy is set to zero for the improper angle.

-  ``harmonic``: The following functional form is used to calculate the
   energy due to an improper angle

   .. math:: E_\psi = K_\psi \left ( \psi - \psi_0 \right )^2

where :math:`K_\psi` and :math:`\psi_0` are specified by the user.

.. table:: Cassandra units for impropers
    :widths: auto
    :align: center

    ================= ================ ================ 
     Parameter         Symbol            Units
    ================= ================ ================ 
     Force constant    :math:`K_\psi`   K/rad\ :sup:`2`
     Improper          :math:`\psi_0`   degrees
    ================= ================ ================ 

Nonbonded
---------

The nonbonded interactions between two atoms :math:`i` and :math:`j` are
due to repulsion-dispersion interactions and electrostatic interactions
(if any).

Repulsion-Dispersion Interactions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The repulsion-dispersion interactions can take one of the following
forms:

-  Lennard-Jones 12-6 potential (LJ):

   .. math:: {\cal V}(r_{ij})= 4 \epsilon_{ij} \left [  \left ( \frac {\sigma_{ij}} { r_{ij} }\right )^{12} - \left ( \frac {\sigma_{ij}} { r_{ij} }\right )^{6}\ \right ]

where :math:`\epsilon_{ij}` and :math:`\sigma_{ij}` are the energy
and size parameters set by the user. For unlike interactions,
different combining rules can be used, as described elsewhere. Note
that this option only evaluates the energy up to a specified cutoff
distance. As described below, analytic tail corrections to the
pressure and energy can be specified to account for the finite cutoff
distance.

-  Cut and shift potential:

   .. math::

      {\cal V}(r_{ij})= 4 \epsilon_{ij} \left [  \left ( \frac {\sigma_{ij}} { r_{ij} }\right )^{12} - \left ( \frac {\sigma_{ij}} { r_{ij} }\right )^{6}\ \right ] -  4 \epsilon_{ij} \left [  \left ( \frac {\sigma_{ij}} { r_{cut}}\right )^{12} - \left ( \frac {\sigma_{ij}} { r_{cut} }\right )^{6}\ \right ]
      \label{Eq:cut_shift}

where :math:`\epsilon_{ij}` and :math:`\sigma_{ij}` are the energy
and size parameters set by the user and :math:`r_{cut}` is the cutoff
distance. This option forces the potential energy to be zero at the
cutoff distance. For unlike interactions, different combining rules
can be used, as described elsewhere.

-  Cut and switch potential:

   .. math::

      {\cal V}(r_{ij})= 4 \epsilon_{ij} \left [  \left ( \frac {\sigma_{ij}} { r_{ij} }\right )^{12} - \left ( \frac {\sigma_{ij}} { r_{ij} }\right )^{6}\ \right ] f
       \label{Eq:cut_switch}

The factor :math:`f` takes the following values:

   .. math::

      \begin{aligned}
          f = 
          \begin{cases}
          
              1.0 \, \, \, &  r_ {ij}  \le r_{on} \\
              \frac { (r_{off}^2 - r_{ij}^2)^2 (r_{off}^2 - 3r_{on}^2 + 2r_{ij}^2)} {\left ( r_{off}^2 - r_{on}^2 \right )^3}  \, \, \,  & r_{on} < r_{ij} < r_{off}\\
              0.0 \, \, \, & r_{ij} \ge r_{off} 
              
          \end{cases}\end{aligned}

where :math:`\epsilon_{ij}` and :math:`\sigma_{ij}` are the energy
and size parameters set by the user. This option smoothly forces the
potential to go to zero at a distance :math:`r_{off}`, and begins
altering the potential at a distance of :math:`r_{on}`. Both of these
parameters must be specified by the user. For unlike interactions,
different combining rules can be used, as described elsewhere.

-  Mie potential (generalized form of LJ):

   .. math::

      {\cal V}(r_{ij})=  \left ( \frac{n}{n-m} \right ) \left ( \frac {n}{m} \right )^{\frac{m}{n-m}}\epsilon_{ij} \left [  \left ( \frac {\sigma_{ij}} { r_{ij} }\right )^{n} - \left ( \frac {\sigma_{ij}} { r_{ij} }\right )^{m}\ \right  ] 
       \label{Eq:mie}

where :math:`\epsilon_{ij}` and :math:`\sigma_{ij}` are the energy
and size parameters and :math:`n` and :math:`m` are the repulsive and
attractive exponents set by the user. This option allows for the use
of a generalized LJ potential (for LJ, :math:`n` = 12 and :math:`m` =
6). Note that this option only evaluates the energy up to a specified
cutoff distance. Both n and m can take on separate integer or float
values set by the user. For unlike interactions, different combining
rules can be used, as described elsewhere.

-  Mie cut and shift potential:

   .. math::

      {\cal V}(r_{ij})=  \left ( \frac{n}{n-m} \right ) \left ( \frac {n}{m} \right )^{\frac{m}{n-m}}\epsilon_{ij} \left [  \left ( \frac {\sigma_{ij}} { r_{ij} }\right )^{n} - \left ( \frac {\sigma_{ij}} { r_{ij} }\right )^{m}\ \right  ] -  \left ( \frac{n}{n-m} \right ) \left ( \frac {n}{m} \right )^{\frac{m}{n-m}}\epsilon_{ij} \left [  \left ( \frac {\sigma_{ij}} { r_{cut}}\right )^{n} - \left ( \frac {\sigma_{ij}} { r_{cut} }\right )^{m}\ \right]  
       \label{Eq:mie_cut_shift}

where :math:`\epsilon_{ij}` and :math:`\sigma_{ij}` are the energy
and size parameters and :math:`n` and :math:`m` are the repulsive and
attractive exponents set by the user. This option forces the
potential energy to be zero at the cutoff distance (i.e. setting
:math:`n` = 12 and :math:`m` = 6 provides the same potential as the
LJ cut and shift option). For unlike interactions, different
combining rules can be used, as described elsewhere.

-  Tail corrections: If the Lennard-Jones potential is used, standard
   Lennard-Jones tail corrections are used to approximate the long range
   dispersion interactions

.. table:: Cassandra units for repulsion-dispersion interactions
    :widths: auto
    :align: center

    =================== ====================== ================ 
     Parameter            Symbol                 Units
    =================== ====================== ================ 
     Energy parameter    :math:`\epsilon/k_B`     K
     Collision diameter  :math:`\sigma`           Å
    =================== ====================== ================ 
    
Electrostatics
~~~~~~~~~~~~~~

Electrostatic interactions are given by Coulomb’s law

.. math::

   {\cal V}_{elec} (r_{ij}) = \frac{1}{4\pi\epsilon_0} \frac {q_i q_j} {r_{ij}}.
   \label{Eq:Coulomb}

where :math:`q_i` and :math:`q_j` are the partial charges set by the
user, which are placed on atomic positions given by :math:`r_i` and
:math:`r_j`. In a simulation, the electrostatic interactions are
calculated using either an Ewald summation, the
`Damped Shifted Force <https://doi.org/10.1063/1.2206581>`_,
or a direct summation using the minimum image convention.
Note that the total energy that is printed out
in the property file is extensive. Consequently, to obtain intensive
energies, the printed energies must divided by the total number of
molecules in the system.


.. table:: Cassandra units for coulombic interactions
    :widths: auto
    :align: center

    ========= ========= =======
    Parameter  Symbol    Units
    ========= ========= =======
    Charge    :math:`q`   e
    ========= ========= =======

Summary of Cassandra units
--------------------------

.. table:: Summary of Cassandra units for input parameters
    :widths: auto
    :align: center

    =========================== ==================================================== ===================
     Item                         Parameter                                               Units
    =========================== ==================================================== ===================
    Bonds                       :math:`l`                                             Å
    Bond angles                 :math:`\theta_0`                                      degrees
    Bond angles                 :math:`K_\theta`                                      K/rad\ :sup:`2`
    OPLS dihedrals              :math:`a_0`, :math:`a_1`, :math:`a_2`, :math:`a_3`    kJ/mol
    CHARMM dihedrals            :math:`a_0`                                           kJ/mol
    CHARMM dihedrals            :math:`a_1`                                           dimensionless
    CHARMM dihedrals            :math:`\delta`                                        degrees
    Harmonic dihedrals          :math:`K_\phi`                                        K/rad\ :sup:`2`
    Harmonic dihedrals          :math:`\phi_0`                                        degrees
    Impropers                   :math:`K_\psi`                                        K/rad\ :sup:`2`
    Impropers                   :math:`\psi_0`                                        degrees

    Simulation box length                                                             Å
    Distances                                                                         Å
    Volume                                                                            Å\ :sup:`3`
    Rotational width                                                                  degrees
    Temperature                                                                       K
    Pressure                                                                          bar
    Chemical potential                                                                kJ/mol
    Energy                                                                            kJ/mol 
    =========================== ==================================================== ===================


