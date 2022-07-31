find_package(xmlf90 QUIET)
if (xmlf90_FOUND)
  message(STATUS "Found xmlf90 cmake package")
  add_library(xmlf90 ALIAS xmlf90::xmlf90-lib)
  #  target_link_libraries(${atom_target} PRIVATE
  #                     xmlf90::xmlf90-lib)

else()

  find_package(PkgConfig QUIET)
  pkg_check_modules(xmlf90 IMPORTED_TARGET GLOBAL xmlf90>=1.5.4)

  if(xmlf90_FOUND)
    message(STATUS "Found xmlf90 through pkgconfig")
    add_library(xmlf90 ALIAS PkgConfig::xmlf90)
    #target_link_libraries(${atom_target} PRIVATE xmlf90)
  elseif(DOWNLOAD_FALLBACK)
    include(FetchContent)
    FetchContent_Declare(xmlf90
        GIT_REPOSITORY https://gitlab.com/siesta-project/libraries/xmlf90
        GIT_TAG cmake
    )
    message(STATUS "... Downloading xmlf90 with git into build hierarchy")
    FetchContent_MakeAvailable(xmlf90)
    message(STATUS "... Adding xmlf90 target as a dependency")
    # I cannot see why this does not work...
    # target_link_libraries(${atom_target} PRIVATE xmlf90::xmlf90-lib)
    # ... whereas this does
   # target_link_libraries(${atom_target} PRIVATE xmlf90)
  else()
    message(FATAL_ERROR "Cannot find xmlf90... enable download or use a submodule")
  endif()
endif()
