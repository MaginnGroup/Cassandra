|License|
|Citing|
|Version|
|Azure|
|GitHub|

.. |License| image:: https://img.shields.io/badge/license-GPL--3.0-green
   :target: reference/license.html
.. |Citing| image:: https://img.shields.io/badge/cite-Cassandra-blue
   :target: reference/citing.html
.. |Version| image:: https://img.shields.io/conda/vn/conda-forge/cassandra
   :target: https://anaconda.org/conda-forge/cassandra
.. |Azure| image:: https://dev.azure.com/MaginnGroup/Cassandra/_apis/build/status/MaginnGroup.Cassandra?branchName=master
.. |GitHub| image:: https://img.shields.io/badge/contribute_on-GitHub-lightgrey
   :target: https://github.com/MaginnGroup/Cassandra


Overview
~~~~~~~~

Cassandra is an open source Monte Carlo software package developed in the
`Maginn group <http://sites.nd.edu/maginn-group/>`_ at the
University of Notre Dame. It is designed to perform atomistic simulations
of molecules composed of rings, chains, or both.

.. warning::

    ReadTheDocs documentation for Cassandra is currently under
    construction. The `PDF reference manual
    <https://github.com/MaginnGroup/Cassandra/releases/latest/download/user_guide.pdf>`_
    is still considered the authoritative source during our transition to
    ReadTheDocs. This message will be removed once the transition to
    ReadTheDocs is complete.

Resources
~~~~~~~~~

* :doc:`Installation guide <getting_started/install>`: Instructions for installing Cassandra
* :ref:`basics`: The basic workflow to use Cassandra
* `MoSDeF Cassandra <https://mosdef-cassandra.readthedocs.io>`_ : A full-fledged
  Python wrapper for Cassandra
* `GitHub repository <https://github.com/MaginnGroup/Cassandra>`_: View the source code, contribute, and raise issues
* `Workshop Materials
  <https://cassandra.nd.edu/images/code/cassandra_workshop_materials_June2016.tar.gz>`_: Notes on statistical
  mechanics and additional tutorials from a June 2016 Cassandra Workshop


.. Cassandra is suited to compute the thermodynamic properties of fluids
   and phase equilibria. It handles a standard "Class I"-type force field
   having fixed bond lengths. Cassandra uses OpenMP parallelization and
   comes with a number of scripts, utilities and examples to help with
   simulation setup. It is released under the GNU General Public License.


Citation
~~~~~~~~

Please :doc:`cite our publication <reference/citing>` if you use this software as part of your research.

Installation
~~~~~~~~~~~~

Complete installation instructions are :doc:`here <getting_started/install>`.
A conda installation is available:

.. code-block:: bash

    conda create --name mc -c conda-forge cassandra

Credits
~~~~~~~
Development of Cassandra was supported by the National Science Foundation
under grant number ACI-1339785. Any opinions, findings, and conclusions or
recommendations expressed in this material are those of the author(s) and do
not necessarily reflect the views of the National Science Foundation.

Complete acknowledgements can be found :doc:`here <reference/acknowledgements>`.
