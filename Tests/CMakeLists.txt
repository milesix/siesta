file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/install_tests")

#
#  Function to set up tests
#
# It accepts a few arguments
# NAME : directory of test
# PSEUDO_DIR : files to be copied in as pseudos
# MPI_NPROC : number of MPI processors
# EXTRA_FILES : extra files to be copied in to the run directory
# SERIAL : to skip MPI (regardless of other arguments)
function(siesta_test)
  set(options SERIAL)
  set(oneValueArgs NAME MPI_NPROC PSEUDO_DIR)
  set(multiValueArgs EXTRA_FILES)
  cmake_parse_arguments(_stest "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

  # Check for bad arguments
  if(DEFINED _stest_UNPARSED_ARGUMENTS)
    message(FATAL_ERROR "Unparsed arguments in siesta_test test=${_stest_NAME}")
  endif()

  # Options
  if(NOT DEFINED _stest_SERIAL)
    set(_stest_SERIAL FALSE)
  endif()

  if(NOT DEFINED _stest_NAME)
    message(FATAL_ERROR "siesta_test missing test directory argument: NAME")
  endif()
  if(DEFINED _stest_MPI_NPROC)
    if(WITH_MPI)
      if(_stest_MPI_NPROC GREATER num_ranks)
        message(WARNING "Reducing number of ranks in test ${_stest_NAME} to ${num_ranks}<${_stest_MPI_NPROC}")
        set(_stest_MPI_NPROC ${num_ranks})
      endif()
    endif()
  else()
    # set default value
    set(_stest_MPI_NPROC ${num_ranks})
  endif()

  # Create the working directory
  set(_stest_wdir "${CMAKE_BINARY_DIR}/install_tests/${_stest_NAME}")
  file(MAKE_DIRECTORY "${_stest_wdir}")

  file(COPY "${_stest_NAME}/${_stest_NAME}.fdf"
    DESTINATION "${_stest_wdir}")

  if(NOT DEFINED _stest_PSEUDO_DIR)
    set(_stest_PSEUDO_DIR "${CMAKE_SOURCE_DIR}/Tests/Pseudos")
    message(VERBOSE "test(${_stest_NAME}) will be using pseudos from: ${_stest_PSEUDO_DIR}")
  endif()

  if(DEFINED _stest_EXTRA_FILES)
    message(STATUS "test(${_stest_NAME}) copying extra files: ${_stest_EXTRA_FILES}")
    foreach(_extra IN LISTS _stest_EXTRA_FILES)
      file(COPY "${_stest_NAME}/${_extra}" DESTINATION "${_stest_wdir}")
    endforeach()
  endif()

  if(WITH_MPI AND (NOT _stest_SERIAL))
    set(mpiexec_args
      ${MPIEXEC_NUMPROC_FLAG} ${_stest_MPI_NPROC}
      ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:siesta-siesta> -out OUT ${_stest_NAME}.fdf ${MPIEXEC_POSTFLAGS})
    add_test(
      NAME siesta-${_stest_NAME}_mpi_np_${_stest_MPI_NPROC}
      COMMAND ${MPIEXEC_EXECUTABLE} ${mpiexec_args}
      WORKING_DIRECTORY "${_stest_wdir}")
	     
    set_tests_properties(siesta-${_stest_NAME}_mpi_np_${_stest_MPI_NPROC}
      PROPERTIES
      ENVIRONMENT "SIESTA_PS_PATH=${_stest_PSEUDO_DIR}")

  else()
 
    add_test(
      NAME siesta-${_stest_NAME}
      COMMAND $<TARGET_FILE:siesta-siesta> -out OUT ${_stest_NAME}.fdf
      WORKING_DIRECTORY "${_stest_wdir}" )

    set_tests_properties(siesta-${_stest_NAME}
      PROPERTIES
      ENVIRONMENT "SIESTA_PS_PATH=${_stest_PSEUDO_DIR}")

  endif()

endfunction()



siesta_test(NAME h2o)
#siesta_test(NAME ch4)
siesta_test(NAME wannier
  EXTRA_FILES wannier.nnkp wannier.win)

setup_test(NAME si_pdos_kgrid)

if(WITH_LUA)
  siesta_test(NAME lua_h2o
    EXTRA_FILES siesta.lua)
endif()

if(WITH_DFTD3)
  siesta_test(NAME dftd3)
endif()
