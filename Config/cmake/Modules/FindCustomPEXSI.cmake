# find PEXSI via brute force for now
#
# This is a kludge to work around some deficiencies
# in target installation in the pexsi-2.1 Cmake build system

# The environment variable PEXSI_ROOT must point to the
# root of the pexsi-2.1 installation

include(FindPackageHandleStandardArgs)

# In pre-release 2.1 the (static) libraries are all
# in ${PEXSI_ROOT}/lib
# ... but we provide extra hooks for other places
# for the subordinate libraries
#

find_library(PEXSI_LIBRARIES
  NAMES pexsi
  PATH_SUFFIXES lib
  HINTS
  ENV PEXSI_ROOT
  DOC "PEXSI libraries list")

find_library(PARMETIS_LIBRARIES
  NAMES parmetis
  PATH_SUFFIXES lib
  HINTS
  ENV PEXSI_ROOT
  ENV PARMETIS_ROOT
  DOC "Parmetis libraries list")

find_library(METIS_LIBRARIES
  NAMES metis
  PATH_SUFFIXES lib
  HINTS
  ENV PEXSI_ROOT
  ENV METIS_ROOT
  DOC "Metis libraries list")

find_library(SUPERLU_DIST_LIBRARIES
  NAMES superlu_dist
  PATH_SUFFIXES lib
  HINTS
  ENV PEXSI_ROOT
  ENV SUPERLU_DIST_ROOT
  DOC "SuperLU_dist libraries list")

find_path(PEXSI_INCLUDE_DIR
  NAMES f_ppexsi_interface.mod
  PATH_SUFFIXES include
  HINTS
  ENV PEXSI_ROOT)

find_package_handle_standard_args(CustomPEXSI "DEFAULT_MSG"
                                  PEXSI_LIBRARIES PARMETIS_LIBRARIES METIS_LIBRARIES
				  SUPERLU_DIST_LIBRARIES PEXSI_INCLUDE_DIR)


message("PEXSI_INCLUDE_DIR: ${PEXSI_INCLUDE_DIR}")
message("PEXSI_LIBRARIES: ${PEXSI_LIBRARIES}")
message("PARMETIS_LIBRARIES: ${PARMETIS_LIBRARIES}")
message("METIS_LIBRARIES: ${PARMETIS_LIBRARIES}")
message("SUPERLU_LIBRARIES: ${SUPERLU_DIST_LIBRARIES}")


if(CustomPEXSI_FOUND AND NOT TARGET PEXSI::PEXSI)


# Instead of 'STATIC' (the current setting for pre-release 2.1),
# we might want to use 'SHARED' if appropriate (to be implemented)

  add_library(pexsi STATIC IMPORTED)
  set_target_properties(pexsi PROPERTIES
  			      IMPORTED_LOCATION "${PEXSI_LIBRARIES}"
                              IMPORTED_LINK_INTERFACE_LANGUAGES "CXX")
  add_library(parmetis STATIC IMPORTED)
  set_target_properties(parmetis PROPERTIES
  			      IMPORTED_LOCATION "${PARMETIS_LIBRARIES}"
                              IMPORTED_LINK_INTERFACE_LANGUAGES "CXX")
  add_library(metis STATIC IMPORTED)
  set_target_properties(metis PROPERTIES
  			      IMPORTED_LOCATION "${METIS_LIBRARIES}"
                              IMPORTED_LINK_INTERFACE_LANGUAGES "CXX")
  add_library(superlu_dist STATIC IMPORTED)
  set_target_properties(superlu_dist PROPERTIES
  			      IMPORTED_LOCATION "${SUPERLU_DIST_LIBRARIES}"
                              IMPORTED_LINK_INTERFACE_LANGUAGES "CXX")

  
  add_library(PEXSI::PEXSI INTERFACE IMPORTED)

  # I need to add by hand "-lstdc++"
  # This is non-optimal, since other platforms might have a different
  # name for the C++ library

  target_link_libraries(PEXSI::PEXSI
                          INTERFACE
                          pexsi
			  parmetis
			  metis
			  superlu_dist)
			  
  set_target_properties(PEXSI::PEXSI
                        PROPERTIES
                        INTERFACE_INCLUDE_DIRECTORIES "${PEXSI_INCLUDE_DIR}")

  # This does not seem to be needed, but YMMV
  #  target_link_libraries(PEXSI::PEXSI INTERFACE MPI::MPI_CXX)
  
endif()

