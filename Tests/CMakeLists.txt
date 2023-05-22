file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/install_tests")

#
#  Function to set up tests
#
# It accepts a few arguments
# NAME : directory of test (optional, an unnamed argument will be accepted as NAME)
# PSEUDO_DIR : files to be copied in as pseudos
# MPI_NPROC : number of MPI processors
# OMP_NPROC : number of OpenMP processors
# EXTRA_FILES : extra files to be copied in to the run directory
# SERIAL : to skip MPI (regardless of other arguments)
# This function uses these external variables:
#  WITH_MPI
#  WITH_OPENMP
#  num_ranks
function(siesta_test)
  set(options SERIAL)
  set(oneValueArgs NAME MPI_NPROC PSEUDO_DIR)
  set(multiValueArgs EXTRA_FILES)
  cmake_parse_arguments(_stest "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

  # Check for bad arguments

  # A single unparsed-argument will be the equivalent of NAME
  # in case it has not been passed.
  if(DEFINED _stest_UNPARSED_ARGUMENTS)
    if(NOT DEFINED _stest_NAME)
      # get name from the first argument
      list(POP_BACK _stest_UNPARSED_ARGUMENTS _stest_NAME)
    endif()
    list(LENGTH _stest_UNPARSED_ARGUMENTS _stest_len)
    if(_stest_len GREATER 0)
      message(FATAL_ERROR "Unparsed arguments in siesta_test test=${_stest_NAME}"
        "Arguments are:\n"
        "  ${_stest_UNPARSED_ARGUMENTS}")
    endif()
  endif()


  # Keyword only options
  # In this case serial may be used to overwrite MPI and OpenMP tests
  if(NOT DEFINED _stest_SERIAL)
    set(_stest_SERIAL FALSE)
  endif()

  # Name of test *must* be defined
  if(NOT DEFINED _stest_NAME)
    message(FATAL_ERROR "siesta_test missing test directory argument: NAME")
  endif()

  message(VERBOSE "Adding new test ${_stest_NAME}")
  list(APPEND CMAKE_MESSAGE_INDENT "  ")

  # Check for number of ranks. It will be reduced to {num_ranks} if too high.
  if(DEFINED _stest_MPI_NPROC)
    if(WITH_MPI)
      if(_stest_MPI_NPROC GREATER num_ranks)
        message(WARNING "Reducing number of MPI processors in test ${_stest_NAME} to ${num_ranks}<${_stest_MPI_NPROC}")
        set(_stest_MPI_NPROC ${num_ranks})
      endif()
    endif()
    message(VERBOSE "using MPI with ${_stest_MPI_NPROC} processors")
  else()
    # set default value
    set(_stest_MPI_NPROC ${num_ranks})
  endif()

  # Check for OpenMP
  # It will understand this variable as the maximum value
  if(DEFINED ENV{OMP_NUM_THREADS})
    # set the OMP_NPROC to this value if undefined
    if(DEFINED _stest_OMP_NPROC)
      if(_stest_OMP_NPROC GREATER ENV{OMP_NUM_THREADS})
        message(WARNING "Reducing number of OpenMP ranks in test ${_stest_NAME} to ${ENV{OMP_NUM_THREADS}}<${_stest_OMP_NPROC}")
        set(_stest_OMP_NPROC ${ENV{OMP_NUM_THREADS}})
      endif()
    else()
      set(_stest_OMP_NPROC 1)
    endif()
    message(VERBOSE "using OpenMP with ${_stest_OMP_NPROC} processors")
  else()
    # set default value of 1
    set(_stest_OMP_NPROC 1)
  endif()

  if(NOT DEFINED _stest_PSEUDO_DIR)
    set(_stest_PSEUDO_DIR "${CMAKE_CURRENT_SOURCE_DIR}/Pseudos")
  endif()
  message(VERBOSE "test(${_stest_NAME}) will be using pseudos from: ${_stest_PSEUDO_DIR}")

  # Build up the variables
  set(_stest_test ${_stest_NAME}) # name of test
  set(_stest_env "SIESTA_PS_PATH=${_stest_PSEUDO_DIR}") # env-vars to add
  set(_stest_pre "")
  set(_stest_post "")
  if( WITH_MPI AND (NOT _stest_SERIAL) )
    set(_stest_pre
      ${MPIEXEC_NUMPROC_FLAG} ${_stest_MPI_NPROC} ${MPIEXEC_PREFLAGS}
      )
    set(_stest_post
      ${MPIEXEC_POSTFLAGS}
      )
    set(_stest_test "${_stest_test}_mpi_np${_stest_MPI_NPROC}")
  endif()

  if( WITH_OPENMP AND (NOT _stest_SERIAL) )
    list(APPEND _stest_env "OMP_NUM_THREADS=${_stest_OMP_NPROC}")
    set(_stest_test "${_stest_test}_omp_np${_stest_OMP_NPROC}")
  endif()

  # Create the working directory
  # Since we might have multiple tests for the same directory,
  # we will need to parse the options before doing anything
  set(_stest_wdir "${CMAKE_CURRENT_BINARY_DIR}/${_stest_test}")

  # Now we can prepare the test directory:
  file(MAKE_DIRECTORY "${_stest_wdir}")
  file(COPY "${_stest_NAME}/${_stest_NAME}.fdf"
    DESTINATION "${_stest_wdir}")

  # Copy in extra files
  if(DEFINED _stest_EXTRA_FILES)
    message(VERBOSE "test(${_stest_NAME}) copying extra files: ${_stest_EXTRA_FILES}")
    foreach(_extra IN LISTS _stest_EXTRA_FILES)
      file(COPY "${_stest_NAME}/${_extra}" DESTINATION "${_stest_wdir}")
    endforeach()
  endif()

  # Now add the test

if( WITH_MPI AND (NOT _stest_SERIAL) )
  add_test(
    NAME "siesta-${_stest_test}"
    COMMAND ${MPIEXEC_EXECUTABLE} ${_stest_pre} $<TARGET_FILE:siesta-siesta> -out OUT ${_stest_NAME}.fdf ${_stest_post}
    WORKING_DIRECTORY "${_stest_wdir}"  )
else()
  add_test(
    NAME "siesta-${_stest_test}"
    COMMAND $<TARGET_FILE:siesta-siesta> -out OUT ${_stest_NAME}.fdf ${_stest_post}
    WORKING_DIRECTORY "${_stest_wdir}"  )
endif()

  # Add environment variables to the test
  foreach(_env IN LISTS _stest_env)
    set_tests_properties(siesta-${_stest_test}
      PROPERTIES
      ENVIRONMENT "${_env}")
  endforeach()

  list(POP_BACK CMAKE_MESSAGE_INDENT)

endfunction()


siesta_test(h2o)
#siesta_test(ch4)
siesta_test(wannier
  EXTRA_FILES wannier.nnkp wannier.win)

siesta_test(si_pdos_kgrid)

if(WITH_LUA)
  siesta_test(lua_h2o
    EXTRA_FILES siesta.lua)
endif()

if(WITH_DFTD3)
  siesta_test(NAME dftd3)
endif()
