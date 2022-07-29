find_package(gridxc QUIET)
if (gridxc_FOUND)
  message(STATUS "Found gridxc cmake package")
  # This promotes the target to GLOBAL, so that an alias can be
  # defined (see https://stackoverflow.com/questions/45401212/how-to-make-imported-target-global-afterwards)
  set_target_properties(gridxc::gridxc-lib PROPERTIES IMPORTED_GLOBAL TRUE)
  add_library(gridxc ALIAS gridxc::gridxc-lib)

  # If the package contains libxc support, linking is automatic

else()

  find_package(PkgConfig QUIET)
  if (WITH_MPI)
    pkg_check_modules(gridxc IMPORTED_TARGET GLOBAL libgridxc_dp_mpi>=0.10.0)
  else()
    pkg_check_modules(gridxc IMPORTED_TARGET GLOBAL libgridxc_dp>=0.10.0)
  endif()  

  if(gridxc_FOUND)
    message(STATUS "Found gridxc through pkgconfig")
    if (WITH_LIBXC)
      target_link_libraries(PkgConfig::gridxc INTERFACE ${LIBXC_LINK_LIBRARIES})
    endif()
    add_library(gridxc ALIAS PkgConfig::gridxc)

  elseif(DOWNLOAD_FALLBACK)
    include(FetchContent)
    #
    # Since we are compiling gridxc, we request LIBXC support
    # This assumes that we can find libxc (through pkgconfig)
    #
    if(WITH_LIBXC)
      set(WITH_LIBXC CACHE INTERNAL "ON")
    endif()
    #
    FetchContent_Declare(gridxc
        GIT_REPOSITORY https://gitlab.com/siesta-project/libraries/libgridxc
        GIT_TAG cmake
    )
    message(STATUS "... Downloading gridxc with git into build hierarchy")
    FetchContent_MakeAvailable(gridxc)
    message(STATUS "... Adding gridxc target as a dependency")

  else()
    message(FATAL_ERROR "Cannot find gridxc... enable download or use a submodule")
  endif()
endif()
