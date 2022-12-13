#
# -- Search for Psolver using pkg-config data
# 
find_package(PkgConfig REQUIRED)
pkg_check_modules(PSOLVER psolver)

if (PSOLVER_FOUND)

   message(DEBUG "PSOLVER Libdir: ${PSOLVER_LIBDIR}")
   message(DEBUG "Psolver Libraries: ${PSOLVER_LIBRARIES}")
   message(STATUS "Psolver Link libraries: ${PSOLVER_LINK_LIBRARIES}")

   message(STATUS " ----- ")
   message(STATUS "   ** Make sure that there are no stray linear algebra libraries **")
   message(STATUS "   ** in the above list of libraries                             **")
   message(STATUS " ----- ")
 
   message(DEBUG "PSOLVER_INCLUDEDIR: ${PSOLVER_INCLUDEDIR}")
   message(STATUS "PSOLVER_INCLUDE_DIRS: ${PSOLVER_INCLUDE_DIRS}")

   add_library(psolver::psolver INTERFACE IMPORTED GLOBAL)
   target_link_libraries(psolver::psolver INTERFACE  ${PSOLVER_LINK_LIBRARIES})
   target_include_directories(psolver::psolver INTERFACE ${PSOLVER_INCLUDE_DIRS})

else()

#
# Compile sources on-the-fly
#

include(ExternalProject)

set(config_flags)  # parameters desired for ./configure of Autotools

# Where things will be installed for internal use
set(install_prefix "${CMAKE_CURRENT_BINARY_DIR}/install")

message(STATUS " --- Entering External/BigDFT")

message(STATUS "     -- LAPACK_LIBRARY: ${LAPACK_LIBRARY}")
if(WITH_MPI)
  message(STATUS "     -- MPI_Fortran_COMPILER: ${MPI_Fortran_COMPILER}")
  set(Fortran_COMPILER "${MPI_Fortran_COMPILER}")
else()
  set(Fortran_COMPILER "${CMAKE_Fortran_COMPILER}")
endif()

if(WITH_OPENMP)
  set(config_flags "--with-openmp") 
else()
  set(config_flags "--with-openmp=no") 
endif()

#
# Complain early if basic tools are not found
#
find_program(MAKE_EXECUTABLE NAMES gmake make mingw32-make REQUIRED)
find_program(AUTORECONF_EXECUTABLE NAMES autoreconf REQUIRED)

#
# --- yaml:  basic C compilation
#            Avoid using shared libraries
#
ExternalProject_Add(yaml
 URL ${CMAKE_CURRENT_LIST_DIR}/yaml
 UPDATE_DISCONNECTED true
 CONFIGURE_HANDLED_BY_BUILD true
 CONFIGURE_COMMAND <SOURCE_DIR>/configure
                   CC=${CMAKE_C_COMPILER}
		   --enable-shared=no --prefix=${install_prefix}
 BUILD_COMMAND ${MAKE_EXECUTABLE} -j
 INSTALL_COMMAND ${MAKE_EXECUTABLE} -j install
 TEST_COMMAND ""
)

ExternalProject_Add_Step(yaml
  autoreconf
  COMMAND ${AUTORECONF_EXECUTABLE} -ifv
  DEPENDEES download
  DEPENDERS configure
  WORKING_DIRECTORY <SOURCE_DIR>
)

#
# --- futile:  Use MPI if enabled, and re-use LAPACK variable
#              Avoid using shared libraries
#              OpenMP used if enabled (see setting of config_flags above)
#              Fortran flags as set at the end of the flags.cmake file

ExternalProject_Add(futile
 URL ${CMAKE_CURRENT_LIST_DIR}/futile/
 PATCH_COMMAND patch -p1 -i "${CMAKE_SOURCE_DIR}/External/BigDFT/Patches/futile.patch"
 UPDATE_DISCONNECTED true
 CONFIGURE_HANDLED_BY_BUILD true
 CONFIGURE_COMMAND <SOURCE_DIR>/configure
                   FC=${Fortran_COMPILER}
                   FCFLAGS=${Fortran_FLAGS_CURRENT}
		   --with-yaml-path=${install_prefix}
		   --with-ext-linalg=${LAPACK_LIBRARY}
		   --enable-shared=no
		   ${config_flags} --prefix=${install_prefix}
 BUILD_COMMAND ${MAKE_EXECUTABLE} -j
 INSTALL_COMMAND ${MAKE_EXECUTABLE} -j install
 TEST_COMMAND ""
)

ExternalProject_Add_Step(futile
  autoreconf
  COMMAND ${AUTORECONF_EXECUTABLE} -ifv
  DEPENDEES download
  DEPENDERS configure
  WORKING_DIRECTORY <SOURCE_DIR>
)

add_dependencies(futile yaml)

#
# ---   atlab: Use previous .pc files to provide information

ExternalProject_Add(atlab
 URL ${CMAKE_CURRENT_LIST_DIR}/atlab
 UPDATE_DISCONNECTED true
 CONFIGURE_HANDLED_BY_BUILD true
 CONFIGURE_COMMAND FC=${Fortran_COMPILER}
                   PKG_CONFIG_PATH=${install_prefix}/lib/pkgconfig
                   <SOURCE_DIR>/configure  FCFLAGS=${Fortran_FLAGS_CURRENT}
		   ${config_flags} --prefix=${install_prefix}
 BUILD_COMMAND ${MAKE_EXECUTABLE} -j
 INSTALL_COMMAND ${MAKE_EXECUTABLE} -j install
 TEST_COMMAND ""
)

ExternalProject_Add_Step(atlab
  autoreconf
  COMMAND ${AUTORECONF_EXECUTABLE} -ifv
  DEPENDEES download
  DEPENDERS configure
  WORKING_DIRECTORY <SOURCE_DIR>
)

add_dependencies(atlab futile)

#
# ---   Psolver: Use previous .pc files to provide information
#

ExternalProject_Add(psolver
 URL ${CMAKE_CURRENT_LIST_DIR}/psolver
 PATCH_COMMAND patch -p1 -i "${CMAKE_SOURCE_DIR}/External/BigDFT/Patches/psolver.patch"
 UPDATE_DISCONNECTED true
 CONFIGURE_HANDLED_BY_BUILD true
 CONFIGURE_COMMAND <SOURCE_DIR>/configure FC=${Fortran_COMPILER} FCFLAGS=${Fortran_FLAGS_CURRENT}
                   PKG_CONFIG_PATH=${install_prefix}/lib/pkgconfig
		   ${config_flags} --prefix=${install_prefix}
 BUILD_COMMAND ${MAKE_EXECUTABLE} -j
 INSTALL_COMMAND ${MAKE_EXECUTABLE} -j install
 TEST_COMMAND ""
)

ExternalProject_Add_Step(psolver
  autoreconf
  COMMAND ${AUTORECONF_EXECUTABLE} -ifv
  DEPENDEES download
  DEPENDERS configure
  WORKING_DIRECTORY <SOURCE_DIR>
)

add_dependencies(psolver atlab)

#
# Note that the library names are hardwired here
# The installation in ${install_prefix} provides a .pc
# file, but we cannot use it... (or I do not know how)
     # (Maybe use a custom command to run pkg-config and pick up the right libraries)
 
 set(PSOLVER_LIBRARIES "-L${install_prefix}/lib -lPSolver-1 -latlab-1 -lfutile-1 ${LAPACK_LIBRARY} -L${install_prefix}/lib -lyaml -ldl")

 add_library(psolver::psolver INTERFACE IMPORTED GLOBAL)

 file(MAKE_DIRECTORY "${install_prefix}/include")

 SET_TARGET_PROPERTIES(psolver::psolver PROPERTIES
      INTERFACE_LINK_LIBRARIES       ${PSOLVER_LIBRARIES}
      INTERFACE_INCLUDE_DIRECTORIES "${install_prefix}/include"
 )
 
# This is needed to trigger the compilation of the external projects

add_dependencies(psolver::psolver psolver)

endif()

