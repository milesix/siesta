Installation of Siesta with CMake

This is a work-in-progress implementation. Please report problems and suggestions.

It needs cmake >= 3.14, and (if used) the ninja (>=1.10) backend.
Both cmake and ninja can be installed easily in most systems, in
particular with conda.

The CMake approach facilitates the handling of the external libraries
that must be installed before Siesta can be compiled. Depending on the
needs and experience of the user, several modes of operation are
available.

_Experimental_) SPACK packages are available in Config/spack_package_defs

   After setting up spack and defining compilers, etc, a user can simply install
   a new repo with the 'siesta-project' namespace:

     spack repo add /path/to/spack_package_defs

   and issue commands such as:

     spack install xmlf90
     spack info siesta
     spack spec siesta +mpi +netcdf +libxc +elpa
     spack install siesta -mpi build_type=Debug

   Note that the spack builtin repo *might have* other Siesta-related
   recipes prepared in the past by other members of the community.
   By making sure that the 'siesta-project' repo is listed first in the
   spack repository chain, those can be avoided. Check:

     spack repo list
     
   At this point, the spack recipes pull sources from development
   branches of Siesta, without a well-defined source id. This is
   a temporary situation during the final stages of development of
   the CMake framework and the update of the dependency libraries.
   
A) Installation on-the-fly, using git submodules for the external
libraries.

- The xmlf90, libpsml, and libgridxc (required) libraries are
  available as submodules in the External subdirectory, set to the
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

  (If you got a full-source package with the External subdirectories already
   populated, you do not need to do anything at this point).

  The basic invocation will then be: (-GNinja is optional)
  
  cmake -S. -GNinja -B_build -DCMAKE_INSTALL_PREFIX=/path/to/installation
  cmake --build _build
    
  (Optionally followed by  '-- -j 1 -v'  for more info)
  (Optionally followed by  '-- target'  to build a specific target)

  cmake --install _build

If MPI is found on the system, it will be used by default (WITH_MPI=ON).

If NetCDF is found on the system, it will be used by default (WITH_NETCDF=ON), and
so will the NCDF library (See below for a fuller option list).

The system will try to find a LAPACK library automatically, using a
custom module. The search can be helped by the setting of the
LAPACK_LIBRARY variable, and/or a directory path LAPACK_LIBRARY_DIR in
which to search.

  cmake [ other options ]  -DLAPACK_LIBRARY="openblas"
  cmake [ other options ]  -DLAPACK_LIBRARY_DIR="/opt/lapack"

The search will also take into account the setting of the variables
BLAS_LIBRARIES and LAPACK_LIBRARIES, which are handled by the standard
CMake 'findLAPACK' module.

NOTE: In systems in which Lapack is implicitly linked (e.g., Cray systems),
it might be necessary to set LAPACK_LIBRARY to "NONE". Newer versions of
CMake detect this automatically.
  
If WITH_MPI=ON, the system will search also for a Scalapack library
(honoring if present the variables SCALAPACK_LIBRARY and
SCALAPACK_LIBRARY_DIR). Please watch out for MPI compatibility issues,
as they are NOT currently checked. The same "implicit linking" note
above applies.

(The CMake support for the linear-algebra library search has been
borrowed from the DFTB+ and SIRIUS projects.)

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

      - WITH_MPI:   It is ON by default. To deactivate, use -DWITH_MPI=OFF

      - WITH_OPENMP: It is OFF by default. To activate, use -DWITH_OPENMP=ON.
                     This option affects only Siesta and TBTrans for now.
		     Check the setting of OMP_NUM_THREADS for execution.

      - WITH_NETCDF: It is ON by default.
                     At the same time, it will set WITH_NCDF to ON.
		     WITH_NCDF_PARALLEL will remain OFF by default.

      - WITH_NCDF:   Its default value is determined by that of WITH_NETCDF.
                     To activate/deactivate independently, use -DWITH_NCDF=ON/OFF

      - WITH_NCDF_PARALLEL: It is OFF by default. Only meaningful if NetCDF is ON.
                            To enable, use -DWITH_NCDF_PARALLEL=ON

      - WITH_LIBXC: When compiling GRIDXC, use the optional libxc library. The latter
                    must be pre-compiled, and incorporated into the build by adding
		    the appropriate $LIBXC_ROOT to the CMAKE_PREFIX_PATH.

      - WITH_GRID_SP: Use single-precision (SP) for grid arrays.
                      It needs to find a SP libgridxc.

      - WITH_ELPA:  See config/cmake/search_for_elpa.cmake. Due to the
                    idiosyncrasies of the ELPA build system, some preparation
		    might be needed.
		    
      - WITH_FFTW:  (ON by default) If not set, Util/STM/ol-stm will
                    not be compiled.

      - WITH_LUA: The flook library can be compiled on-the-fly, using
                  the code in a git submodule in
                  External/Lua-Engine/flook, or in a tarball pointed
                  to by FLOOK_PACKAGE. Pre-compiled installations can
                  be used by setting the environment variable
                  FLOOK_ROOT or a CMake variable (-DFLOOK_ROOT=/path/to/flook),
		  or by adding $FLOOK_ROOT to CMAKE_PREFIX_PATH.

      - WITH_PSOLVER: The psolver library can be compiled on-the-fly,
                  using the code in git submodules in External/BigDFT,
                  Pre-compiled installations can be used by setting
                  the environment variable PSOLVER_ROOT or a CMake
                  variable (-DPSOLVER_ROOT=/path/to/psolver), or by adding
                  $PSOLVER_ROOT to CMAKE_PREFIX_PATH.

* Handling of compiler flags

  For Intel and GNU compilers, a set of default flags for the
  different build_types (Debug, Check, Release, RelWithDebInfo, None) are
  hard-coded in the file Config/cmake/flags.cmake.

  Those settings can be changed in the command line by setting either
  of

       Fortran_FLAGS            (will affect all build_types)
       Fortran_FLAGS_DEBUG      (this in the sense of "not optimized", with debug symbols)
       Fortran_FLAGS_CHECK      (with extra checks on array bounds, pointers, etc)
       Fortran_FLAGS_RELEASE
       Fortran_FLAGS_RELWITHDEBINFO   (release with debug symbols)
       Fortran_FLAGS_NONE      

* Toolchain files

  As an experimental feature, some pre-packaged "toolchain" files can
  be found in the Config/cmake/toolchains directory. They contain
  compiler flags and library search hints.  They can be used as:

     cmake -C /path/to/toolchain/file ......

  or as

     cmake -DCMAKE_TOOLCHAIN_FILE=/path/to/toolchain/file

  These files can also be used to set options. Some variables can be
  in the toolchain file, and some others given in the command line:

     cmake -C /path/to/toolchain/file -DWITH_MPI=OFF 

  If a variable is set *both* in the toolchain file and in the command
  line, the command-line value takes precedence.
  
 * Quirks, workarounds, and missing bits:

    * The handling of external libraries should be improved:

      - ELPA is handled through its pkgconfig files,
        which might not be well constructed (see the modules
        in config/cmake)

      - When using MPI, a single MPI version of GRIDXC is compiled (or
        found) and it is used also for the utility programs that would
        only need a serial version (since they call only the 'atomxc'
        routine).  This will work in most systems, but in some others
        (e.g. IBM's mpi) one would have to use 'mpirun ...' for serial
        utilities linking to GRIDXC (gen-basis, for example).

    * Tests of the installation are very basic for now

       Go into the _build directory and type:

       ctest  [options]

      If the required external libraries have been compiled as part of
      the current CMake invocation, installation tests for them will
      also be executed.
      
(Thanks to DFTB+ and other projects for ideas and support files)

 