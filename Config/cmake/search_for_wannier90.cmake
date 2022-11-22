#
# Search for old-style, manually compiled, libwannier90 installations.
# Note that the preferred way to install is through automatic patching
# of a pristine wannier90 package or a pristine wannier90 source directory.
#
# We rely on the environmental variable WANNIER90_ROOT
#
set(WANNIER90_ROOT "$ENV{WANNIER90_ROOT}" CACHE FILEPATH "wannier90 installation path")

 if("${WANNIER90_ROOT}" STREQUAL "" )

    message(SEND_ERROR "Variable WANNIER90_ROOT not set")
     
 else()
 
  message("Using wannier90 installation in ${WANNIER90_ROOT}")
  set(WANNIER90_LIBS "-L${WANNIER90_ROOT}/lib -lwannier")
  set(WANNIER90_INCLUDE_DIRS "${WANNIER90_ROOT}/include")
  message("WANNIER90_LIBS: ${WANNIER90_LIBS}")
  message("WANNIER90_INCLUDE_DIRS: ${WANNIER90_INCLUDE_DIRS}")

  add_library(libwannier90::libwannier90 INTERFACE IMPORTED GLOBAL)
  target_link_libraries(libwannier90::libwannier90
               INTERFACE "${WANNIER90_LIBS}")
  target_include_directories(libwannier90::libwannier90
               INTERFACE
              "${WANNIER90_INCLUDE_DIRS}")

 endif()
