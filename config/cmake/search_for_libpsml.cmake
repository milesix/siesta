#
# Note that the project name is 'psml', not 'libpsml',
# in versions of the library with cmake support.
#
# The library target name coming out of this module
# is nevertheless libpsml
#
find_package(psml QUIET)
if (psml_FOUND)
  message(STATUS "Found libpsml cmake package")
  set_target_properties(psml::psml-lib PROPERTIES IMPORTED_GLOBAL TRUE)
  add_library(libpsml ALIAS psml::psml-lib)

else()

  find_package(PkgConfig QUIET)
  pkg_check_modules(libpsml IMPORTED_TARGET GLOBAL libpsml>=1.1.10)

  if(libpsml_FOUND)
    message(STATUS "Found libpsml through pkgconfig")
    add_library(libpsml ALIAS PkgConfig::libpsml)
  elseif(DOWNLOAD_FALLBACK)
    include(FetchContent)
    FetchContent_Declare(libpsml
        GIT_REPOSITORY https://gitlab.com/siesta-project/libraries/libpsml
        GIT_TAG cmake
    )
    message(STATUS "... Downloading libpsml with git into build hierarchy")
    FetchContent_MakeAvailable(libpsml)
    message(STATUS "... Adding libpsml target as a dependency")
  else()
    message(FATAL_ERROR "Cannot find libpsml... enable download or use a submodule")
  endif()
endif()
