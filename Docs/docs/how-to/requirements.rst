.. sectionauthor:: Yann Pouillon <y.pouillon@simuneatomistics.com>

Minimum requirements
====================

In this section, we describe what the minimum requirements to build SIESTA
are, so that a computer expert could jump to the build itself with this
information at hand. In :doc:`build-prep-env`, we detail how to install these
requirements.

Hardware requirements
---------------------

SIESTA can work on a broad variety of computer architectures, ranging from
Raspberry PI to massively parallel supercomputers. To make it work, you will
need at least:

- 1 Gb of RAM, even if small atomic structures can even be accessed with less.
- Around 2 Gb of disk storage, in particular if you want to run the SIESTA
  test suite.

SIESTA will likely work on any CPU from various vendors. Please note, however,
that most tests and benchmarks performed during development are using Intel
and ARM architectures.

For all versions of SIESTA
--------------------------

To build any version of SIESTA, you need at least:

- A C compiler.
- A C++ compiler (might not be strictly necessary, still highly recommended).
- A Fortran compiler.
- The GNU Make utility.

In addition, if you want to run a parallel version of the SIESTA executable,
you will have to install a MPI distribution.

To build the SIESTA user manual, you will also need a LaTeX distribution.

Differences between SIESTA 4.x and other versions
-------------------------------------------------

Up to SIESTA 4.x, the source code came with all the necessary dependencies to
build a serial executable without having to install anything more than the
minimum requirements. However, from versions 5.x and above, as well as for
MaX-related versions, the source code is distributed without external
libraries. This is an application of the `separation of concerns`_ software
design principle.

As a consequence, a set of libraries has to be installed before the latter
versions of SIESTA can be compiled. They include:

- LibFDF
- LibGridXC
- LibPSML
- XMLF90

In turn, it is highly recommended to install LibXC and use it within
LibGridXC, in order to be able to access the 600+ exchange-correlation
functionals provided by LibXC from SIESTA.

.. _`separation of concerns`: https://en.wikipedia.org/wiki/Separation_of_concerns
