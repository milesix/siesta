#
# Inspired by the DFTB+ project
# Stops the code if the source and the build folders are identical.
#
function(siesta_util_ensure_out_of_source_build)

  get_filename_component(srcdir "${CMAKE_CURRENT_SOURCE_DIR}" REALPATH)
  get_filename_component(bindir "${CMAKE_CURRENT_BINARY_DIR}" REALPATH)

  if("${srcdir}" STREQUAL "${bindir}")
    message(FATAL_ERROR
      "It is not allowed to configure and build this project from its source folder. Please, create a \
separate build directory and invoke CMake from that directory.")
  endif()

endfunction()



#[==========================[
Function to dynamically add a directory as a build-dependency
if the directory exists and has some files in it
The number of files should probably be tweaked, well.

This routine accepts a set of arguments to handle how sub-directories
are added.

  DIRECTORY : required
    which directory to dynamically add
  OPTION : required
    name of the option that is exposed to the user
  NITEMS : required
    number of items that should be in DIRECTORY before it is eligeble for
    adding as a sub-project.
  HELP : optional
    help text exposed for the OPTION, defaults to a description of the OPTION and DIRECTORY
  DEFAULT : optional
    a default value of OPTION in case it is eligeble for addition.
    Defaults to TRUE.

#]==========================]
function(add_subdirectory_option)
  set(options "")
  set(oneValueArgs DIRECTORY OPTION HELP NITEMS DEFAULT)
  set(multiValueArgs "")
  cmake_parse_arguments(_asop "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

  if(NOT DEFINED _asop_DIRECTORY)
    message(FATAL_ERROR "add_subdirectory_option requires DIRECTORY argument")
  endif()
  if(NOT DEFINED _asop_OPTION)
    message(FATAL_ERROR "add_subdirectory_option requires OPTION argument")
  endif()
  if(NOT DEFINED _asop_NITEMS)
    set(_asop_NFILES 1)
    message(VERBOSE "add_subdirectory_option set NITEMS to 1 (more than 1 file+directory in DIRECTORY)")
  endif()
  if(NOT DEFINED _asop_DEFAULT)
    set(_asop_DEFAULT TRUE)
    message(VERBOSE "add_subdirectory_option set DEFAULT to TRUE if the directory exists")
  endif()
  if(NOT DEFINED _asop_HELP)
    set(_asop_HELP
      "Include support for option ${_asop_OPTION} with a default value of ${_asop_DEFAULT} "
      "if the folder ${_asop_DIRECTORY} exists with >=${_asop_NITEMS} files+directories present.")
  endif()

  message(VERBOSE "Checking directory ${_asop_DIRECTORY} for content")

  file(GLOB _result
    LIST_DIRECTORIES TRUE
    "${_asop_DIRECTORY}/*"
    )
  list(LENGTH _result _result_len)
  message(VERBOSE "Directory ${_asop_DIRECTORY} contains ${_result_len} files and directories")
  if(_result_len GREATER _asop_NITEMS)
    option(${_asop_OPTION} "${_asop_HELP}" ${_asop_DEFAULT})
  else()
    set(${_asop_OPTION} FALSE CACHE BOOL "${_asop_HELP}" FORCE)
  endif()
  if ( ${${_asop_OPTION}} )
    message(STATUS "Adding support with ${_asop_OPTION}=TRUE")
    list(APPEND CMAKE_MESSAGE_INDENT "  ")
    add_subdirectory("${_asop_DIRECTORY}")
    list(POP_BACK CMAKE_MESSAGE_INDENT)
  endif()
endfunction()
