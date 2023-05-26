file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/install_tests")

#
#  Function to set up tests
#
function(setup_test name extrafiles)
 set(wdir "${CMAKE_BINARY_DIR}/install_tests/${name}")
 file(MAKE_DIRECTORY "${wdir}")
 
 file(COPY "${name}/${name}.fdf"  DESTINATION "${wdir}")

 if(extrafiles)
   separate_arguments(extrafiles)
   foreach(_extra IN LISTS extrafiles)
     file(COPY ${name}/${_extra} DESTINATION "${wdir}")
   endforeach()
 endif()

 if(WITH_MPI)
   set(mpiexec_args
          ${MPIEXEC_NUMPROC_FLAG} ${num_ranks}
          ${MPIEXEC_PREFLAGS}
	  $<TARGET_FILE:siesta-siesta> -out OUT ${name}.fdf
          ${MPIEXEC_POSTFLAGS})
	  
   add_test( NAME    siesta-${name}_mpi_np_${num_ranks}
             COMMAND  ${MPIEXEC_EXECUTABLE} ${mpiexec_args}
             WORKING_DIRECTORY "${wdir}")
	     
   set_tests_properties(siesta-${name}_mpi_np_${num_ranks}
                       PROPERTIES
                       ENVIRONMENT "SIESTA_PS_PATH=${CMAKE_SOURCE_DIR}/Tests/Pseudos")
 		       
    
 else()
 
   add_test( NAME    siesta-${name}
             COMMAND $<TARGET_FILE:siesta-siesta> -out OUT ${name}.fdf
             WORKING_DIRECTORY "${wdir}" )

   set_tests_properties(siesta-${name}
                       PROPERTIES
                       ENVIRONMENT "SIESTA_PS_PATH=${CMAKE_SOURCE_DIR}/Tests/Pseudos")
		       
 endif(WITH_MPI)

endfunction()




setup_test("h2o" "")
#setup_test("ch4" "")
#setup_test("si_pdos_kgrid" "")
#setup_test("wannier" "wannier.nnkp wannier.win")
setup_test("gaas_fatbands" "")
setup_test("si_fatbands" "")


if(WITH_LUA)
  setup_test("lua_h2o" "siesta.lua")
endif()

if(WITH_DFTD3)
  setup_test("dftd3" "")
endif()
