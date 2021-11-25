find_package(psml QUIET)
if (libpsml_FOUND)
  message(STATUS "Found libpsml cmake package")
  add_library(libpsml ALIAS psml::libpsml-lib)

else()

  find_package(PkgConfig QUIET)
  pkg_check_modules(libpsml IMPORTED_TARGET GLOBAL libpsml>=1.1.9)

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
