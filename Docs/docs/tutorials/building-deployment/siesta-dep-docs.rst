:sequential_nav: next

..  _tutorial-deployment:

SIESTA Deployment Options
=========================

    :Author: Vladimir Dikan <vdikan@icmab.es> (ICMAB-CSIC)

.. contents:: Table of Contents
   :depth: 3

As you might know (or will realize during the school), SIESTA has
quite a number of capabilities and operation modes. And and as many
other HPC codes, it relies on quite a number of dependencies and build
options, that sometimes give their users hard times with configuration
of their research environments.

Below is a review of some aspects of compilation and deployment of
SIESTA. Some pre-configured options are discussed, followed by an
overview of SIESTA's general Makefile template and dependencies for a
manual compilation. Finally, a couple of scripted provisioning
systems are mentioned to give users an idea of how building if
SIESTA can be streamlined for different environments.

Before continuing please, remember to first and foremost address the
`SIESTA manual <https://siesta-project.org/SIESTA_MATERIAL/Docs/Manuals/siesta-MaX-1.0-3.pdf>`_ and in particular the "Compilation" section for
instructions on installation of SIESTA.

Also please keep in mind that this is not an exhaustive overview of all
possible installation options for SIESTA. After the recent migration
of its codebase to the `GitLab repository <https://gitlab.com/siesta-project/siesta>`_, other installation and
packaging systems are developed further to provide deployment options
for SIESTA. Stay updated with the SIESTA project announcements!

1 Ready-to-Use Options
----------------------

For personal use, either just to try siesta or to run standard and not
very complex jobs, there are virtualized "plug-and-play" variants
available:

- ``QuantumMobile`` virtual machine with a number of computational
  codes pre-installed, along with pre-configured for them ``AiiDA``
  instance.

- Conda installation: See :ref:`this how-to<building_with_conda>`.

- Containers for different platforms such as ``Docker`` and
  ``Singularity``

1.1 QuantumMobile VM
~~~~~~~~~~~~~~~~~~~~

The `QuantumMobile <https://quantum-mobile.readthedocs.io>`_ virtual machine is set up with pre-installed SIESTA,
as well as with a number of other computational codes curated by
MARVEL and MaX European initiatives. It may be the easiest
installation option for beginners. In fact, QuantumMobile will be used
during the SIESTA school for its students' accounts.

.. note::
   It's recommended to use `VirtualBox 6.1 <https://www.virtualbox.org/>`_ with QuantumMobile.

[Screenshots not included in this version]

Download the image from `QuantumMobile Releases <https://quantum-mobile.readthedocs.io/en/latest/releases/index.html>`_ page, launch VirtualBox
and select the ``.ova`` file in the Import Appliance menu.

After launching the VM one can see that ``siesta`` is available in the
terminal emulator, along with some of its useful utilities like
``tbtrans``, ``vibra``, ``gnubands``, ``denchar`` etc.

Moreover, after launching AiiDA activation command (as the splash
screen suggests) we obtain AiiDA's ``verdi`` shell with siesta
code already pre-configured to be used with AiiDA.

The header of SIESTA output shows the code configuration options and
dependencies with which the code was built inside the VM.

1.2 Containers
~~~~~~~~~~~~~~

When a more lightweight variant of virtual environment is preferred,
the containers such as `Docker <https://www.docker.com/>`_
and `Singularity <https://sylabs.io/singularity/>`_
can be used. We provide
experimental containers for SIESTA, and plan to auto-generate them for
future releases.

1.2.1 Docker
^^^^^^^^^^^^

Docker images for SIESTA can be found on this `Docker Hub page <https://hub.docker.com/r/vdikan/siesta-dist/tags?page=1&ordering=last_updated>`_.
Sample usage:

.. code:: sh

    docker pull vdikan/siesta-dist:master
    docker run --interactive --tty -w /app -v "$(pwd):/app" vdikan/siesta-dist:master

Assuming one has all required input in the current local directory,
the command above will launch a container for **siesta@master**
development branch, mount said directory as a volume under ``/app`` path
and pass control of the terminal inside the container to the user.
Executables from the Siesta package will be visible and should work
while inside the container, and after shutdown the local directory
mounted as ``/app`` should persist together with all the generated output.

1.2.2 Singularity
^^^^^^^^^^^^^^^^^

The Singularity project is an alternative to Docker aimed primarily
for scientists and their supercomputers. At a glance, using
Singularity is more convenient because it helps running containerized
software in a manner more familiar to researchers (that is, running
binaries within a scope of some directory path).

It is arguably most convenient to convert existing Docker images to
Singularity:

.. code:: sh

    singularity build siesta-dist.sif docker://vdikan/siesta-dist:master

Will produce a Singularity container image from the Docker image from
the previous section, and

.. code:: sh

    singularity run siesta-dist.sif

will launch it in a way similar to the Docker container.

.. note::
   There is no need to specifically designate and mount a
   volume with the current folder so that the results of calculations are
   saved: Singularity `does so automatically <https://sylabs.io/guides/3.0/user-guide/bind_paths_and_mounts.html>`_.

2 Source Code Compilation
-------------------------

The detailed instructions how to build SIESTA from source code are
written in the `SIESTA manual <https://siesta-project.org/SIESTA_MATERIAL/Docs/Manuals/siesta-MaX-1.0-3.pdf>`_.
Here some details are clarified.

2.1 Staging the build directory
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

It's recommended to obtain stable versions of SIESTA from the Siesta
`Project releases <https://gitlab.com/siesta-project/siesta/-/releases>`_ page, e.g.:

.. code:: sh

      wget https://gitlab.com/siesta-project/siesta/-/archive/v4.1.5/siesta-v4.1.5.tar.gz
      tar -xzf siesta-v4.1.5.tar.gz
      cd siesta-v4.1.5/
      ls

    AUTHORS COPYING Docs Examples NOTICE.txt Obj Pseudo README.md ReleaseNotes.md Src Tests Tutorials Util version.info

The source code is placed in ``Src/`` directory, but as the manual
states, compilation inside ``Src/`` is not allowed and will cause an
error:

.. code:: sh

      cd Src/
      make

    ** You can no longer build SIESTA in Src.
    ** Go to the Obj directory and see the README file.
    Makefile:62: recipe for target 'what' failed
    make: *** [what] Error 1

The correct way to build SIESTA is to use a separate build directory.
The default such directory is present in the archive and is called
``Obj/``, but users may have any number of build directories on the same
level, with different configurations of SIESTA.

First, the build directory should be staged by calling a setup script:

.. code:: sh

      cd Obj/
      sh ../Src/obj_setup.sh

    *** Compilation setup done.
    *** Remember to copy an arch.make file into the directory.
    *** These files are template arch.make files:
    ***    gfortran.make (for gfortran compiler)
    ***    intel.make (for intel compiler)
    ***    DOCUMENTED-TEMPLATE.make (requires customization)

Then an ``arch.make`` file should be placed with options that customize
the build. The output of previous command lists a few templates
located in ``Obj/``. New releases also come with a number of well-structured
and self-documented ``arch.make`` templates located in
``Obj/ARCH_EXPERIMENTAL/``

.. code:: sh

       tree Obj/

    Obj/
    ├── ARCH-EXPERIMENTAL
    │   ├── cte-gcc.mk
    │   ├── cte-ibm.mk
    │   ├── gcc-modules.mk
    │   ├── master-raw.make
    │   ├── mn-intel.mk
    │   └── README
    ├── DOCUMENTED-TEMPLATE.make
    ├── gfortran.make
    ├── intel.make
    └── README

It is recommended to use these new templates to compile SIESTA.

.. code:: sh

    cp ARCH-EXPERIMENTAL/master-raw.make ./arch.make
    $EDITOR arch.make

2.2 ``arch.make`` configuration
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The structure of the new ``arch.make`` templates is very
self-explanatory. Users should edit only the header part before this
line:

.. code:: sh

    #--------------------------------------------------------
    # Nothing should need to be changed below
    #--------------------------------------------------------

First section specifies internal flags for SIESTA dependencies. Set
them to ``1`` for the required options, e.g. ``WITH_MPI=1`` The rest leave blank.

The next section defines linking symbols for mandatory dependencies as well
for some external libraries marked as required. When requesting any of
them, uncomment the corresponding line and put the appropriate library
symbols. Please, contact your system administrator to find out the
correct library locations on a shared system.

Example: ``netcdf`` dependency switched  on with separate
``netcdf-fortran`` interface:

.. code:: sh

    WITH_NETCDF=1
    WITH_SEPARATE_NETCDF_FORTRAN=1
    WITH_NCDF=1
    ...
    NETCDF_ROOT=$(NETCDF_HOME)  # /path/to/netcdf-c
    NETCDF_FORTRAN_ROOT=$(NETCDF_HOME)  # /path/to/netcdf-fortran

Finally, specify compiler options that aren't visible or do differ
from the shell environment (alternative ways possible, see the inline
documentation commentaries).

.. code:: sh

    FC_PARALLEL=mpif90
    FC_SERIAL=gfortran
    FPP = $(FC_SERIAL) -E -P -x c
    FFLAGS = -O2
    FFLAGS_DEBUG= -g -O0
    RANLIB=echo
    ...

2.3 Siesta Dependencies
~~~~~~~~~~~~~~~~~~~~~~~

The dependencies for SIESTA can be divided in two groups: *common*
ones that are widely used for many projects in HPC and computational
research, and *specific* for SIESTA project: libraries used in siesta
program (almost) exclusively, written by authors of SIESTA. The
majority of SIESTA-specific dependencies are grouped on the
corresponding section on the project page: `https://gitlab.com/siesta-project/libraries <https://gitlab.com/siesta-project/libraries>`_

.. note::
   SIESTA code already ships with minimum necessary
   libraries. In principle, even use of ``MPI`` is optional: SIESTA can be
   built and run in serial mpde. The only requirement is installed
   ``SCALAPACK`` library in case when ``MPI`` is requested.

The following table lists some external libraries that can be
linked with SIESTA and extend its capabilities:

.. table::

    +------------------------------------------------------------------+----------------------------------------------------------------------+
    | Common                                                           | SIESTA-specific                                                      |
    +==================================================================+======================================================================+
    | MPI, BLAS, LAPACK                                                | `libfdf <https://gitlab.com/siesta-project/libraries/libfdf>`_       |
    +------------------------------------------------------------------+----------------------------------------------------------------------+
    | `ELSI <https://wordpress.elsi-interchange.org/>`_ + ext. solvers | `xmlf90 <https://gitlab.com/siesta-project/libraries/xmlf90>`_       |
    +------------------------------------------------------------------+----------------------------------------------------------------------+
    | `NetCDF <https://www.unidata.ucar.edu/software/netcdf/>`_        | `libPSML <https://gitlab.com/siesta-project/libraries/libpsml>`_     |
    +------------------------------------------------------------------+----------------------------------------------------------------------+
    | `libXC <https://www.tddft.org/programs/libxc/>`_                 | `libGridXC <https://gitlab.com/siesta-project/libraries/libgridxc>`_ |
    +------------------------------------------------------------------+----------------------------------------------------------------------+
    | `flook <https://github.com/ElectronicStructureLibrary/flook>`_   | `flos <https://github.com/siesta-project/flos>`_                     |
    +------------------------------------------------------------------+----------------------------------------------------------------------+

2.4 Building Utilities
~~~~~~~~~~~~~~~~~~~~~~

Many useful features of SIESTA are extracted as separate utilities in
the ``Util/`` library. One can build them separately or perform a batch
compilation with ``Util/build_all.sh``.

.. note::
   Utilities rely on the ``arch.make`` configuration for SIESTA.

Build them after successful compilation of the ``siesta`` executable in
``Obj/``, and/or edit the corresponding ``Makefile``-s accordingly
(paying attention to environment variables pointing to ``siesta``, e.g.
``OBJDIR``). Then proceed with:

.. code:: sh

    cd Util && ./build_all.sh

3 Scripted Installations
------------------------

Accurate management of dependencies and build environment for SIESTA
(or any comparably complex scientific code) can be cumbersome. That is
why several scripted systems exist that can aid setting up of said
environment.

Those systems and packages differ in complexity and operation quality,
being either a set of ``Make``-scripts, ``Git`` actions, or a set of
packages for software managers such as ``Spack``. Users and
system administrators have a choice of options in cases where a
`1 Ready-to-Use Options`_ solution does not suite.

A couple of variants are described in this section.

3.1 Siesta-Install-Scripts
~~~~~~~~~~~~~~~~~~~~~~~~~~

A lightweight set of scripts for manual installation of siesta and its
dependencies. Authored by **Alberto Garcia** and **Xe Hu**.

The `project repository <https://gitlab.com/mailhexu/siesta-install-scripts>`_ contains branches adapted for different operating
systems and HPC clusters; e.g. users of ``Debian/Ubuntu`` can use
`this branch <https://gitlab.com/mailhexu/siesta-install-scripts/-/tree/ubuntu>`_. For documentation address the ``README`` files and inline
documentation comments.

.. code:: sh

    git clone https://gitlab.com/mailhexu/siesta-install-scripts.git
    cd siesta-install-scripts/
    # 1. address README-s in subdirectories
    # 2. inspect ./Tarballs/download.sh
    ./Tarballs/download.sh
    cp cp Config/gnu/* Tarballs/  # copy Makefile for GNU target platform
    cp Config/siesta.common.arch.make Tarballs/
    # 3. inspect and adapt Makefiles
    # 4. inspect and edit ./do_all.sh
    ./do_all.sh

3.2 Spack Siesta package
~~~~~~~~~~~~~~~~~~~~~~~~

It is possible to install SIESTA with experimental
`Spack <https://spack.io/>`_ package.
Spack is a package manager targeting research software and
supercomputers, although it proves useful even for software management
on a personal machine. Please consult Spack documentation for
installation instructions.

.. warning::
   In order to get Spack siesta package for now one needs to clone
   this `dedicated branch of spack <https://github.com/vdikan/spack/tree/siesta-develop>`_
   (called ``siesta-develop``), not the core spack repository! The rest of
   installation described in the docs is valid for said branch as for
   core spack v.15.4-16.0

We are working towards providing a single one-liner command to install
SIESTA releases with a single one-line specification for Spack.
At the moment you can build experimental versions of SIESTA with Spack
obtained through:

.. code:: sh

    git clone -b siesta-develop https://github.com/vdikan/spack.git

After configuration of compilers for Spack, the installation of SIESTA
is done in principle with a single spec command, e.g.:

.. code:: sh

    spack install siesta@master +utils ^openmpi +cxx +cxx_exceptions

installs parallel version of siesta from GitLab's @master branch with
mpi provided by openmpi, with C++ support, as well as key siesta
utilities like ``denchar``, ``tbtrans`` and ``stm``.

In order to make the executables available either run:

.. code:: sh

    spack load -r siesta@master    # "@" precedes the version installed by spack

or address Spack's built-in `modulefiles generation mechanism <https://spack.readthedocs.io/en/latest/module_file_support.html>`_.

3.2.1 Available Spack versions of SIESTA:
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

At the moment there are few SIESTA versions visible for the
experimental Spack package, namely:

- **@master** - for the master branch of the project on GitLab

- **@psml** - for the branch with PSML pseudopotentials support (downloaded from Git)

- **@elsi** - for the branch with ELSI+PEXSI support (downloaded from Git). Requires MPI built with cxx, as in the example above.

- **@4.1-b4** -for the stable version hosted on Launchpad
