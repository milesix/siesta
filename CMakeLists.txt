
cmake_minimum_required(VERSION 3.14)

project(
  "siesta"
  LANGUAGES Fortran C
  VERSION "4.1.6"
  DESCRIPTION "A great DFT code"
)

option(WITH_MPI "Whether Siesta should support MPI-parallelism" FALSE)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Debug")
endif()

# Follow GNU conventions for installing directories
include(GNUInstallDirs)


# General configuration information
add_subdirectory("config")

# -- Internal library dependencies
if (WITH_MPI)
 add_subdirectory("Src/MPI")
endif()

add_subdirectory("Src/fdf")
add_subdirectory("Src/wxml")
add_subdirectory("Src/SiestaXC")
add_subdirectory("Src/xmlparser")

# -- Siesta executable
add_subdirectory("Src")

# -- Utilities

add_subdirectory("Util")

# add the testsuite
#enable_testing()
#add_subdirectory("Testers")
