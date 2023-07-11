# find PEXSI via brute force for now

include(FindPackageHandleStandardArgs)

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
  DOC "Parmetis libraries list")

find_library(METIS_LIBRARIES
  NAMES metis
  PATH_SUFFIXES lib
  HINTS
  ENV PEXSI_ROOT
  DOC "Metis libraries list")

find_library(SUPERLU_LIBRARIES
  NAMES superlu_dist
  PATH_SUFFIXES lib
  HINTS
  ENV PEXSI_ROOT
  DOC "SuperLU_dist libraries list")

find_path(PEXSI_INCLUDE_DIR
  NAMES f_ppexsi_interface.mod
  PATH_SUFFIXES include
  HINTS
  ENV PEXSI_ROOT)

find_package_handle_standard_args(CustomPEXSI "DEFAULT_MSG"
                                  PEXSI_LIBRARIES PARMETIS_LIBRARIES METIS_LIBRARIES
				  SUPERLU_LIBRARIES PEXSI_INCLUDE_DIR)


message("PEXSI_INCLUDE_DIR: ${PEXSI_INCLUDE_DIR}")

message("PEXSI_LIBRARIES: ${PEXSI_LIBRARIES}")
message("PARMETIS_LIBRARIES: ${PARMETIS_LIBRARIES}")

if(CustomPEXSI_FOUND AND NOT TARGET PEXSI::PEXSI)
  add_library(PEXSI::PEXSI INTERFACE IMPORTED)
  set_target_properties(PEXSI::PEXSI PROPERTIES
                                     INTERFACE_INCLUDE_DIRECTORIES "${PEXSI_INCLUDE_DIR}"
				     IMPORTED_LINK_INTERFACE_LANGUAGES "CXX"
                                     INTERFACE_LINK_LIBRARIES
		 "${PEXSI_LIBRARIES};${PARMETIS_LIBRARIES};${METIS_LIBRARIES};${SUPERLU_LIBRARIES};-lstdc++")
endif()

