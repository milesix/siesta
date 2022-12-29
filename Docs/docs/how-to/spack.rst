:sequential_nav: next

..  _building_with_spack:

Building Siesta with spack
==========================

`Spack <https://spack.io>`_ is a package manager targeting research
software and supercomputers, although it proves useful even for
software management on a personal machine. Please consult Spack
documentation for installation instructions.

It is possible to install SIESTA with an **experimental spack fork**
created by Vladimir Dikan (not yet with the core spack repository),
which can be obtained through::

    git clone -b siesta-develop https://github.com/vdikan/spack.git

Its installation and configuration is as described in the standard
spack docs (v.15.4-16.0).

After configuration of compilers for Spack, the installation of SIESTA
is done in principle with a single spec command. For example, the command::

    spack install siesta@master +utils ^openmpi +cxx +cxx_exceptions

installs a parallel version of siesta from GitLab's @master branch with
mpi provided by openmpi, with C++ support, as well as key siesta
utilities like ``denchar``, ``tbtrans`` and ``stm``.

In order to make the executables available either run::

    spack load -r siesta@master    # "@" precedes the version installed by spack

or use Spack's built-in `modulefiles generation mechanism <https://spack.readthedocs.io/en/latest/module_file_support.html>`_.

Available Spack versions of SIESTA
-----------------------------------------

At the moment there are a few SIESTA versions visible for the
experimental Spack package, namely:

- **@master** - for the master branch of the project on GitLab

- **@psml** - for the branch with PSML pseudopotentials support (downloaded from Git)

- **@elsi** - for the branch with ELSI+PEXSI support (downloaded from Git). Requires MPI built with cxx, as in the example above.

- **@4.1-b4** -for the stable version hosted on Launchpad



  








