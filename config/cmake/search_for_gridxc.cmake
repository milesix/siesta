#find_package(gridxc QUIET)
if (gridxc_FOUND)
  message(STATUS "Found gridxc cmake package")
  add_library(gridxc ALIAS gridxc::gridxc-lib)
  #  target_link_libraries(${atom_target} PRIVATE
   #                     gridxc::gridxc-lib)
  # If the package contains libxc support, linking is automatic
  # 
else()

  find_package(PkgConfig QUIET)
  pkg_check_modules(gridxc IMPORTED_TARGET GLOBAL gridxc>=0.11.0)

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
    # I cannot see why this does not work...
    # target_link_libraries(${atom_target} PRIVATE gridxc::gridxc-lib)
    # ... whereas this does
    #target_link_libraries(${atom_target} PRIVATE gridxc)
    ##set(GRIDXC_COMPILED "True")
  else()
    message(FATAL_ERROR "Cannot find gridxc... enable download or use a submodule")
  endif()
endif()
