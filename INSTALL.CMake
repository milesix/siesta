Installation of Siesta with CMake

This is a work-in-progress implementation. Please report problems and suggestions.

It needs cmake >= 3.17 (3.14 if libgridxc tests are disabled), and the
ninja (>=1.10) backend. Both cmake and ninja can be installed easily
in most systems, in particular with conda.

The CMake approach facilitates the handling of the external libraries
that must be installed before Siesta can be compiled. Depending on the
needs and experience of the user, several modes of operation are
available.

A) Installation on-the-fly, using git submodules for the external
libraries.

- The xmlf90, libpsml, and libgridxc (required) libraries are
  available as submodules in the ExtLibs subdirectory, set to the
  appropriate branches and/or tags for use in the current version of
  Siesta. To use them, if you got your sources from the git repo, you
  just need to do (just once, after you have fetched/cloned this
  Siesta branch):

     git submodule update --init --recursive <path>

  where <path> is the path to the submodule. Leave out '<path>' to
  enable all submodules.

  (If you plan to compile Siesta in a system without internet access
   (e.g. MareNostrum), you can use the 'stage_submodules.sh' script
   in the top level of the distribution to gather all the source files
   for later uploading)

  (If you got a full-source package with the ExtLibs subdirectories already
   populated, you do not need to do anything at this point).

  The basic invocation will then be:
  
  cmake -S. -GNinja -B_build -DCMAKE_INSTALL_PREFIX=/path/to/installation
  cmake --build _build
    
  (Optionally followed by  '-- -j 1 -v'  for more info)
  (Optionally followed by  '-- target'  to build a specific target)

  cmake --install _build

If MPI is found on the system, it will be used by default (WITH_MPI=ON).

If NetCDF is found on the system, it will be used by default (WITH_NETCDF=ON), and
so will the NCDF library (See below for a fuller option list).

The system will try to find a LAPACK library automatically, but it will honor
the setting of the LAPACK_LIBRARY variable, and/or a directory path LAPACK_LIBRARY_DIR
in which to search.

  cmake [ other options ]  -DLAPACK_LIBRARY="openblas"
  cmake [ other options ]  -DLAPACK_LIBRARY_DIR="/opt/lapack"

A search for a BLAS library will be carried out before searching for LAPACK (honoring
if present the analogous variables BLAS_LIBRARY and BLAS_LIBRARY_DIR).
Targets BLAS::BLAS and LAPACK::LAPACK will be defined, and the latter made to depend on
the first. See the appropriate Config/cmake/Find* files for more details, including the
possible need to set BLAS_LIBRARY or LAPACK_LIBRARY to "NONE" if either of them is
redundant.
  
If WITH_MPI=ON, the system will search also for a Scalapack library (honoring
if present the variables SCALAPACK_LIBRARY and SCALAPACK_LIBRARY_DIR). Please watch out for
MPI compatibility issues, as they are NOT currently checked.

(The CMake support for the linear-algebra library search has been
borrowed from the DFTB+ project).

To have libgridxc compiled with libxc support, change the invocation to:

   cmake -S. -GNinja -B_build -DWITH_LIBXC=ON \
        -DCMAKE_PREFIX_PATH="$LIBXC_ROOT" \
        -DCMAKE_INSTALL_PREFIX=/path/to/installation

   (Libxc must be pre-installed on the system, in the path pointed to by $LIBXC_ROOT)

B) The same mechanism would work without submodules, as long as there
   is network access through the https protocol. In this case (the
   "fetch" scenario), cmake will download the source for the external
   libraries and compile them. The same flags as above apply.

C) Pre-installed external libraries can be discovered using cmake packages or pkg-config.
   
   The basic incantation is:

   cmake -S. -GNinja -B_build \
        -DCMAKE_PREFIX_PATH="$XMLF90_ROOT;$PSML_ROOT;$GRIDXC_ROOT;$LIBXC_ROOT" \
        -DCMAKE_INSTALL_PREFIX=/path/to/installation

   where the *_ROOT variables point to the appropriate installation roots of the libraries.

   Note that when a given library's installation cannot be found, the
   system will fall back to the "submodule" and "fetch" scenarios, in
   that order.

Developers and experienced users can experiment with different
versions of the libraries by setting different values for the *_ROOT
variables, or by selectively leaving out *_ROOT variables from
CMAKE_PREFIX_PATH and populating the submodule directories appropriately.

* Accepted options:

      - WITH_UTILS: Compiles programs in Util besides TBtrans (currently
                    just a few representative examples) (ON by default)
		    NOTE: You need to set and pass $FFTW_ROOT for Util/STM/ol-stm.

      - WITH_MPI:   It is ON by default if MPI is found on the system.
                    To deactivate, use -DWITH_MPI=OFF

      - WITH_NETCDF: It is ON by default if NetCDF is found on the system.
                     At the same time, it will set WITH_NCDF to ON.
		     WITH_NCDF_PARALLEL will remain OFF by default.

      - WITH_NCDF:   It is ON by default if NetCDF is found on the system.
                     To deactivate, use -DWITH_NCDF=OFF

      - WITH_NCDF_PARALLEL: It is OFF by default. Only meaningful if NetCDF is ON.
                            To enable, use -DWITH_NCDF_PARALLEL=ON


      - WITH_LIBXC: When compiling GRIDXC, use the optional libxc library. The latter
                    must be pre-compiled, and incorporated into the build by adding
		    the appropriate $LIBXC_ROOT to the CMAKE_PREFIX_PATH.

      - WITH_ELPA:  See config/cmake/search_for_elpa.cmake. Due to the
                    idiosynchracies of the ELPA build system, some preparation
		    is needed.
		    
      - WITH_LUA:   The environmental variable FLOOK_ROOT must be set and
                    point to an installation of the flook library.


 Quirks, workarounds, and missing bits:
 
    * Version, compiler, and library information is now designed to be
      incorporated into the Siesta executable through a 'configure'd"
      version-info template, but the variables needed are not fully
      implemented. The version is taken from the project declaration.

    * The handling of external libraries should be improved:

      - Lapack and Scalapack need the linker flags in environment
        variables.

      - ELPA is handled through its pkgconfig files,
        which might not be well constructed (see the modules
        in config/cmake)

      - The flook library needs its path to be set in FLOOK_ROOT, and
        the linking string is hard-coded, as there is no proper
        pkgconfig file yet.

      - GRIDXC cannot yet be incorporated from a pre-compiled package. It
        must be compiled from the submodule source in ExtLibs/libgridxc or
        by fetching the source on-the-fly.
        Some more work is needed to interface correctly to versions of
        libgridxc installed with the autotools build system.
        
      - When using MPI, a single MPI version of GRIDXC is compiled, and it
        is used also for those targets that would only need a serial version.
        This will work in most systems, but in some others (e.g. IBM's mpi)
        one would have to use 'mpirun ...' for serial utilities
        linking to GRIDXC (gen-basis, for example).

      - Some utility programs are still not enabled for building with
        CMake. (One can use the WITH_COMPACT_LOG=1 option in the hybrid
	approach (see INSTALL) to ease the burden of preparing lists
        of source files.)


 