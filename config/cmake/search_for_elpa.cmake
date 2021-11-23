#
#
# -- Search for ELPA using pkg-config data
# (This, with further refinemens,
#  should go in a module in the cmake subdir)
# 
# There are several non-standard issues:
#
# 1. The pkgconfig filename is typicall of the form
#
#     elpa-2020.05.001.rc1.pc
#
# instead of simply 'elpa.pc'.
#
# 2. The include directory entry in the pkgconfig file
#    is missing the 'modules' leaf component.
#
# The most expedient solution for (1) is to put a symbolic
# link to the desired version somewhere in the PKG_CONFIG_PATH. For example:
#
#   ln -sf /path/to/original/.pc/file $HOME/lib/pkgconfig/elpa.pc
#
#   export PKG_CONFIG_PATH=$HOME/lib/pkgconfig:${PKG_CONFIG_PATH}
#
# (Note: CMAKE_PREFIX_PATH can also be used, but without the
#  '/lib/pkgconfig' part.
#
# A solution for (2) is encoded below.
#
# 
 find_package(PkgConfig REQUIRED)
 pkg_check_modules(ELPA REQUIRED elpa)

 message(STATUS "ELPA Libdir: ${ELPA_LIBDIR}")
 message(STATUS "Elpa Libraries: ${ELPA_LIBRARIES}")
 message(STATUS "Elpa Link libraries: ${ELPA_LINK_LIBRARIES}")
 message(STATUS "ELPA_INCLUDEDIR: ${ELPA_INCLUDEDIR}")
 message(STATUS "ELPA_INCLUDE_DIRS: ${ELPA_INCLUDE_DIRS}")
 #
 # Fix non-standard setting
 #
 set(ELPA_INCFLAGS "${ELPA_INCLUDE_DIRS}/modules")
 message(STATUS "ELPA_INCFLAGS to be used: ${ELPA_INCFLAGS")
