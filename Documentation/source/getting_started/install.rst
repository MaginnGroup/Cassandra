Installation
============

We recommend the conda installation for beginning users. The conda installation
will install the Cassandra executable and also provide the ``library_setup.py``
and ``mcfgen.py`` auxillary scripts. If you wish to contribute to Cassandra
or have access to an Intel compiler, you may wish to install from source.

.. note::

    The Intel compiler offers substantial performance improvements
    compared to ``gfortran``. If you are running a large number of
    production calculations and have access to the Intel compiler
    you may want to take the time to install from source.


Installing with conda
~~~~~~~~~~~~~~~~~~~~~

If you already have
`conda <https://docs.conda.io/en/latest/miniconda.html>`_ installed,
you can create a new conda environment and install Cassandra with
a single command:

.. code-block:: bash

    conda create --name mc -c conda-forge cassandra

The command creates a new conda environment (``mc``) and installs
``cassandra``. The ``-c`` flag specifies the conda channels that
is used to install cassandra. To use the environment,
run ``conda activate mc``.

After activating the environment, you can test your installation
by checking that the required executables are on your ``PATH``
with the following commands:

.. code-block:: bash

    which cassandra.exe
    which mcfgen.py
    which library_setup.py

The version of ``cassandra`` installed through ``conda`` uses
OpenMP parallelization. The number of parallel threads is controlled
through the ``OMP_NUM_THREADS`` environment variable. For example, to
use eight threads with a bash terminal you would run:

.. code-block:: bash

    export OMP_NUM_THREADS=8

Installing from source
~~~~~~~~~~~~~~~~~~~~~~

Cassandra may alternatively be installed from source. There are two
methods for obtaining the source code: (1) downloading the
tarball of the latest release from our `GitHub releases page
<https://github.com/MaginnGroup/Cassandra/releases/latest/>`_, or (2)
cloning the GitHub repository. The command to clone the repository
is:

.. code-block:: bash

    git clone https://github.com/maginngroup/cassandra.git


If you download the tarball from the GitHub releases page,
you will need to unpack it:

.. code-block:: bash

    tar -xzvf Cassandra-1.2.4.tar.gz

In either case, after obtaining the source code, go into the
``Src`` directory and run the following:

.. code-block:: bash

    make -f Makefile.gfortran.openMP
    cd ../
    mkdir bin/
    mv Src/cassandra_gfortran_openMP.exe ./bin/cassandra.exe
    cp Scripts/Frag_Library_Setup/library_setup.py ./bin/.
    cp Scripts/MCF_Generation/mcfgen.py ./bin/.

.. note::
    There are several different Makefiles in the Src directory.
    The Makefiles with the ``.openMP`` extension use
    OpenMP parallelization.

Finally, if you wish, you can add ``Cassandra-1.2.4/bin``
to your ``PATH``:

.. code-block:: bash

    export PATH=path_to_install/Cassandra-1.2.4/bin:${PATH}

Unless you add the preceding line to your ``.bashrc`` you will need to
run it every time you open a new terminal window.
