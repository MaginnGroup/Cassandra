.. highlight:: none

.. _ch:hma:

Harmonically Mapped Averaging (HMA)
===================================

.. _sec:hma:

Background
----------

Harmonically mapped averaging (HMA) is a technique that improves the precision
of measured themodynamic properites of crystalline systems.  Beyond precision,
other benefits include smaller potential-truncation effects, finite-size
effects, faster equilibration and shorter decorrelation time.

Limitations
-----------

HMA can only be used in the NVT ensemble and when the atoms only vibrate around
their lattice sites (without diffusing).  HMA can be used in Cassandra with any
pair potential, but not with Ewald summation.

Truncation
----------

In order to supress property fluctuations due to atoms moving into and out of
the truncation radius, the potential in Cassandra is truncated based on the
lattice site distance instead of the distance between the atoms.  This will lead
to a small difference in the measured properties, but the difference will
disappear in the limit of long truncation.

Example
-------

An simuation that uses HMA is included in `Examples/NVT/HMA/ <https://github.com/MaginnGroup/Cassandra/tree/master/Examples/NVT/HMA>`_ in the source tree.  The example simulation runs a Lennard-Jones FCC crystal near its melting point.  Energy and pressure are computed using both conventional and HMA and reported in the property output.

References
----------

Sabry G. Moustafa, Andrew J. Schultz, and David A. Kofke, Very fast averaging of thermal properties of crystals by molecular simulation, `Phys. Rev. E [92], 043303 (2015) <https://dx.doi.org/10.1103/PhysRevE.92.043303>`_
