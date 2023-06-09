#
# Vega supercomputer at Maribor, Slovenia, with ELSI + external PEXSI 
#
# Using the gompic easyBuild toolchain (GCC, OpenBLAS, OpenMPI, CUDA)
#
# Set up the modules as:
#
# ml gompic/2020b
# ml OpenBLAS/0.3.12-GCC-10.2.0
# ml ScaLAPACK/2.1.0-gompic-2020b
# ml netCDF-Fortran/4.5.3-gompic-2020b
# ml CMake/3.18.4-GCCcore-10.2.0
# ml FFTW/3.3.8-gompic-2020b
#
#---- If OpenBLAS is threaded, use OpenMP, but set OMP_NUM_THREADS to 1 in most cases.
#
set(WITH_OPENMP "ON" CACHE BOOL "with OpenMP")
#
set(WITH_ELSI "ON" CACHE BOOL "with ELSI")
#
# This is temporarily needed because the discovery of an external PEXSI in the ELSI config file is
# not complete.
#
set(ELSI_HAS_EXTERNAL_PEXSI "ON" CACHE BOOL "the linked ELSI uses an external PEXSI library")

set(LAPACK_LIBRARY
      "-L /cvmfs/sling.si/modules/el7/software/OpenBLAS/0.3.12-GCC-10.2.0/lib -lopenblas -lpthread -lm -ldl"
        CACHE STRING "lapack library chosen")

set(SCALAPACK_LIBRARY
     "-L/cvmfs/sling.si/modules/el7/software/ScaLAPACK/2.1.0-gompic-2020b/lib -lscalapack"
     CACHE STRING "scalapack library chosen")
     
set(CMAKE_BUILD_TYPE "RelWithDebInfo" CACHE STRING "build_type")

#
# You will need to specify the paths to the installation of ELSI and PEXSI
# in the CMAKE_PREFIX_PATH variable, either here on in the command line.
# Assuming those are in the environment variables $ELSI_ROOT and $PEXSI_ROOT:
#
set(CMAKE_PREFIX_PATH "$ENV{ELSI_ROOT};$ENV{PEXSI_ROOT}" CACHE STRING "paths to ELSI and PEXSI")
