.. _sec:trajectory_reader:

Trajectory Reader
=================

Pregenerated trajectories can be read by Cassandra for thermodynamic property calculations and/or Widom test particle insertions.
This can be especially useful for computing chemical potentials for molecular dynamics trajectories, as the trajectory can be generated and written out by a molecular dynamics engine, potentially with more efficient sampling than can be achieved with Cassandra.  The resulting trajectory, once converted to ``.xyz`` and ``.H`` files like those output by Cassandra, can be read frame-by-frame into Cassandra, which can perform Widom insertions on the trajectory frames as it would during a NVT or NPT Cassandra simulation.

Input file
----------

In order to use the trajectory reader, the following features must be included in the input file:

- Under the keyword ``# Sim_Type``, the simulation type must be ``pregen``, rather than a thermodynamic ensemble.

- The section ``# Pregen_Info`` must be included as described in :ref:`sec:input_file`.

- While the ``# Box_Info`` section must be included as for other types of simulations, the box size included there
  will be internally overwritten for each frame by the box size in the pregenerated trajectory ``.H`` file.
  A box size must still be included and the box shape must accurately describe the actual box shape for the pregenerated trajectory.
  For example, if the box size in the input file is ``cubic`` but the box of the pregenerated trajectory is triclinic, 
  the results will not be accurate.

- Under the keyword ``Simulation_Length_Info``, the units must be steps, which are equivalent to simulation frames when the trajectory reader
  is used.  Writing a simulation length that exceeds the number of frames in the pregenerated trajectory does not cause an error, 
  but Cassandra will write a warning in the log file and consider the simulation completed upon trying to read a frame that doesn't exist 
  in the pregenerated trajectory files.  The warning message will include the number of frames Cassandra was able to read.

Unlike for other types of simulations, the sections ``# Move_Probability_Info``, ``# Start_Type``, and ``# Run_Type`` are not required in
the input file and will be ignored if they are present, since they do not apply.

Pregenerated Trajectory Files
-----------------------------

A pregenerated trajectory must be supplied to Cassandra through ``.H`` and ``.xyz`` files, which are described in :ref:`sec:output_files`.
Cassandra does not use the contents of the elements column of the ``.xyz`` file, so it does not need to be accurate, 
though it does need to be present.  The atoms of a molecule must be listed together in the same order as in the corresponding MCF file, 
and all of the molecules of species 1 must be listed first, then species 2, and so on.

A utility script for converting LAMMPS dump files to ``.H`` and ``.xyz`` files is included with Cassandra and documented
in :ref:`sec:lammpstrjconvert`.
