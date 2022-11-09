file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/install_tests")

function(setup_test name pseudos extrafiles)
 set(wdir "${CMAKE_BINARY_DIR}/install_tests/${name}")
 file(MAKE_DIRECTORY "${wdir}")
 
 file(COPY "${name}/${name}.fdf"  DESTINATION "${wdir}")

 separate_arguments(pseudos)
 foreach(_ps IN LISTS pseudos)
   file(COPY Pseudos/${_ps} DESTINATION "${wdir}")
 endforeach()

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
	  
   add_test( NAME    siesta-${name}_mpi
             COMMAND ${MPIEXEC_EXECUTABLE} ${mpiexec_args}
             WORKING_DIRECTORY "${wdir}" )

 else()
 
   add_test( NAME    siesta-${name}
             COMMAND siesta-siesta -out OUT ${name}.fdf
             WORKING_DIRECTORY "${wdir}" )

 endif(WITH_MPI)

endfunction()

#
#
#
setup_test("h2o" "H.psf O.psf" "")
#setup_test("ch4" "H.pbe.psf C.pbe.psf" "")
setup_test("wannier" "Si.psf" "wannier.nnkp wannier.win")


if(WITH_LUA)
  setup_test("lua_h2o" "O.psf H.psf" "siesta.lua")
endif()
