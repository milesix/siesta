
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

# Collect sources

set(srcs)
add_subdirectory("Src/SiestaXC")
add_subdirectory("Src/fdf")
add_subdirectory("Src/xmlparser")
add_subdirectory("Src/wxml")
add_subdirectory("Src")

add_executable(${PROJECT_NAME}-siesta "${srcs}")

set_target_properties(
  ${PROJECT_NAME}-siesta
  PROPERTIES
  VERSION "${PROJECT_VERSION}"
  OUTPUT_NAME "siesta"
)

# find_package("lapack" REQUIRED)
# find_package("blas" REQUIRED)
# set(
#   lib-deps
#   "blas::blas"
#   "lapack::lapack"
# )

# target_link_libraries(
#   ${PROJECT_NAME}-siesta
#   PUBLIC
#   "${lib-deps}"
# )

set(LIBS "-lveclibfort")
TARGET_LINK_LIBRARIES(
  ${PROJECT_NAME}-siesta
  PRIVATE
  ${LIBS})
  
install(
  TARGETS ${PROJECT_NAME}-siesta
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  )


# Package license files
install(
  FILES
  "COPYING"
  DESTINATION "${CMAKE_INSTALL_DATADIR}/licenses/${PROJECT_NAME}"
)

# add the testsuite
#enable_testing()
#add_subdirectory("Testers")
