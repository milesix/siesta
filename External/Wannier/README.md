Preparing Wannier90 for SIESTA calculations
===========================================


Motivations
-----------

Wannier90 is mostly meant to be used as a standalone program, reading
its parameters from input files and writing results to output
files. In these conditions, performing a wannierisation with SIESTA
would mean prepare the input data, interrupt the SIESTA calculation,
run Wannier90 for each manifold of interest, and then continue with
the SIESTA calculation. Such a procedure would invlove several manual
steps, which would make it particularly error-prone.

In order to streamline the whole process and have Wannier90 run
embedded in SIESTA, it is necessary to build a modified version of its
source code, wrapping the main wannier90 program as a subroutine (see
below for more details). These changes allow Wannier90 to receive its
input directly from SIESTA and to be run several times with different
input parameters without having to stop the program. They have been
kept minimalistic, in order not to interfere with the internals of
Wannier90 and make it easier to keep SIESTA synchronised with the
evolution of the source code.


Using the modified Wannier90 in Siesta (CMake only)
---------------------------------------------------

Just define the WANNIER90_PACKAGE environment variable to point to
a pristine wannier90-3.1.0.tar.gz package (possibly remote):

    export WANNIER90_PACKAGE=/path/to/wannier90-3.1.0.tar.gz
    export WANNIER90_PACKAGE=https://github.com/wannier-developers/wannier90/archive/v3.1.0.tar.gz

and set -DWITH_WANNIER90=ON in the CMake invocation line for
Siesta. The build system will unpack and patch the wannier90
distribution automatically, compile the appropriate version
(serial or mpi) of the library, and link it to Siesta.

Simple wrapper for wannier90
----------------------------

Within the context of the Siesta project, a small patch for
wannier90-3.1.0 has been created to enable a simpler interface between
the codes. The main change is the wrapping of the ``wannier_prog``
program to turn it into a subroutine that can be called directly from
another program. That is the only entry point of the API, as
implemented in the ``wannier90_m`` module.

In addition, a simple CMake building system has been implemented to facilitate
the automatic compilation of the wrapper.

The ``wannier90_wrapper`` subroutine accepts as arguments:

* The ``seedname`` or file prefix.
* The MPI communicator to be used (only in the MPI version)
* (Optional) nnkp_mode: a flag to request only the generation of a .nnkp file.
* (Optional) dryrun_mode: a flag to perform a simple dry-run to check the input file.
* (Optional) Arrays to hold the information about k-point neighbors
* (Optional) Arrays to hold the information about the unitary matrices

The last two are convenience arguments to enable client codes to extract the needed information
directly, without having to read the .nnkp and .chk files, respectively.

All other interaction with wannier90 is through files:

* ``seedname.win`` (provided by client) is the input file
* ``seedname.amn`` (provided by client) 
* ``seedname.mmn`` (provided by client)
* ``seedname.{eig_ext}`` (provided by client)
* ``seedname.wout`` (generated by wannier90)
*  ... plus other plotting and auxiliary files generated by wannier90.

The standard command-line interaction with wannier90 has been bypassed.

In MPI operation, no initialization is done by wannier90, and the communicator to be used is
passed explicitly. Internally, the communicator is named ``mpi_comm_w90``. The ``mpif.h`` file
is read once in the auxiliary module w90_mpi, which is used by relevant pieces of the code.

Some minor changes have been needed to make sure that all variables are deallocated at the end
of the wannier90 run, so that the wrapper can be called repeatedly without errors.

### Compilation of the wrapper with CMake

(Note that the standard makefile-based building system with **NOT** work for the wrapper.)

Basic incantation (items in brackets are optional, needed for MPI operation and to help finding
an appropriate Lapack library):

```
  cmake -S. -B _build [ -DWITH_MPI=ON ] [-DLAPACK_LIBRARY=   ] [ -DCMAKE_INSTALL_PREFIX=/path/to/inst ]
  cmake --build _build

  # (optional)  cd _build ; ctest ; cd ..

  # (optional) cmake --install _build
```

  
