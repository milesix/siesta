Installation of Siesta with CMake

This is a work-in-progress implementation. Most features work, but fine
details might not. Please report problems and suggestions.

It needs cmake >= 3.14, and the ninja (>=1.10) backend. Both can be installed
easily in most systems, in particular with conda.

External libraries are discovered using pkg-config. The basic incantation is:

  cmake -S. -GNinja -B_build \
        -DCMAKE_PREFIX_PATH="$XMLF90_ROOT;$PSML_ROOT;$GRIDXC_ROOT;$LIBXC_ROOT" \
        -DCMAKE_INSTALL_PREFIX=/path/to/installation
    
  cmake --build _build
    
  (Optionally followed by  '-- -j 1 -v'  for more info)
  (Optionally followed by  '-- target'  to build a specific target)

  cmake --install _build

Note the explicit mention of the paths to the installation roots of
the libraries. To install those libraries you can use the scripts
provided in the Siesta distribution (in Docs/install_*.sh)

If MPI is found on the system, it will be used by default (WITH_MPI=ON).

If NetCDF is found on the system, it will be used by default (WITH_NETCDF=ON), and
so will the NCDF library (See below for a fuller option list).

In addition to the above, the environment variable LAPACK_LIBS must be
set to the right linker string. If WITH_MPI=ON, the SCALAPACK_LIBS
environment variable must be set appropriately (please watch out for
MPI compatibility issues).


    * Accepted options:

      - WITH_UTILS: Compiles programs in Util besides TBtrans (currently
                    just a few representative examples)

      - WITH_MPI:   It is ON by default if MPI is found on the system.
                    To deactivate, use -DWITH_MPI=OFF

      - WITH_NETCDF: It is ON by default if NetCDF is found on the system.
                     At the same time, it will set WITH_NCDF to ON.
		     WITH_NCDF_PARALLEL will remain OFF by default.

      - WITH_NCDF:   It is ON by default if NetCDF is found on the system.
                     To deactivate, use -DWITH_NCDF=OFF

      - WITH_NCDF_PARALLEL: It is OFF by default. Only meaningful if NetCDF is ON.
                            To enable, use -DWITH_NCDF_PARALLEL=ON


      - WITH_ELPA:  See config/cmake/search_for_elpa.cmake. Due to the
                    idiosynchracies of the ELPA build system, some preparation
		    is needed.
		    
      - WITH_LUA:   The environmental variable FLOOK_ROOT must be set and
                    point to an installation of the flook library.


 Quirks, workarounds, and missing bits:
 
    * Version, compiler, and library information is now designed to
      be incorporated through a 'configure'd" version-info template,
      but the variables needed are not fully implemented. The version
      is taken from the project declaration, not from version.info.

    * The handling of external libraries should be improved:

      - Lapack and Scalapack need the linker flags in environment
        variables.

      - ELPA is handled through its pkgconfig files,
        which might not be well constructed (see the modules
        in config/cmake)

      - The flook library needs its path to be set in FLOOK_ROOT.

      - For GRIDXC, the search will try 'libgridxc_dp_mpi.pc' if
        WITH_MPI=ON, or 'libgridxc_dp.pc' for the serial case. This
        will work in most cases, but in some systems (e.g. IBM's mpi)
        one would have to use 'mpirun ...' for serial utilities
        linking to GRIDXC (gen-basis, for example).

      - Very soon it should be possible to use CMake to build external
        libraries as sub-projects. (See INSTALL for a hybrid approach
        using CMake and makefiles).

      - Some utility programs are still not enabled for building with
        CMake. (One can use the WITH_COMPACT_LOG=1 option in the hybrid
	approach (see INSTALL) to ease the burden of preparing lists
        of source files.)


 