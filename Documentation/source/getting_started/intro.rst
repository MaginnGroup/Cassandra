
Introduction
~~~~~~~~~~~~

Cassandra is an open source Monte Carlo package capable of simulating any number
of molecules composed of rings, chains, or both. It can be used to simulate
compounds such as small organic molecules, oligomers, aqueous solutions and
ionic liquids. It handles standard ''Class I''-type force fields having fixed
bond lengths, harmonic bond angles and improper angles, a CHARMM or OPLS-style
dihedral potential, a Lennard-Jones 12-6 potential and fixed partial charges. It
does *not* treat flexible bond lengths. Cassandra uses OpenMP
parallelization and comes with a number of scripts, utilities and examples to
help with simulation setup.

Cassandra is capable of simulating systems in the following ensembles: 

* Canonical (NVT) 
* Isothermal-isobaric (NPT) 
* Grand canonical (\ :math:`\mu`\ VT) 
* Constant volume Gibbs (NVT-Gibbs) 
* Constant pressure Gibbs (NPT- Gibbs)

..  Distribution
    ============

    Cassandra is distributed as a gzipped tar file ``Cassandra_V*.tar.gz``. You can
    unpack the distribution by running the command::

        tar -xzf Cassandra_V1.2.tar.gz

    Upon successful unpacking of the archive file, the ``Cassandra_V1.2`` directory
    will have a number of subdirectories. Please refer to the ``README`` file in the
    main Cassandra directory for a detailed information on each of the
    subdirectories.

    * ``Documentation`` - contains the user guide
    * ``Examples`` - contains example input files and short simulations of various
      systems in the above ensembles.
    * ``MCF`` - molecular connectivity files for a number of molecules. These can
      be used as the basis for generating your own MCF files for molecules of 
      interest.
    * ``Scripts`` - useful scripts to set up simulation input files.
    * ``Src`` - Cassandra source code.
