#
# Inspired by the DFTB+ project
# Stops the code if the source and the build folders are identical.
#
function(ensure_out_of_source_build)

  get_filename_component(srcdir "${CMAKE_CURRENT_SOURCE_DIR}" REALPATH)
  get_filename_component(bindir "${CMAKE_CURRENT_BINARY_DIR}" REALPATH)

  if("${srcdir}" STREQUAL "${bindir}")
    message(FATAL_ERROR
      "It is not allowed to configure and build this project from its source folder. Please, create a \
separate build directory and invoke CMake from that directory.")
  endif()

endfunction()

