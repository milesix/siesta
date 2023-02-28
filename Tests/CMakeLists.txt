file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/install_tests")

#
#  Function to set up tests
#
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
	  
   add_test( NAME    siesta-${name}_mpi_np_${num_ranks}
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
setup_test("si_pdos_kgrid" "Si.psf" "")
setup_test("pdos_gamma_nc_fe" "Fe_soc.psml" "")
setup_test("pdos_gamma_soc_fe" "Fe_soc.psml" "")
setup_test("wannier" "Si.psf" "wannier.nnkp wannier.win")


if(WITH_LUA)
  setup_test("lua_h2o" "O.psf H.psf" "siesta.lua")
endif()

if(WITH_DFTD3)
  setup_test("dftd3" "C.psf N.psf O.psf H.psf" "")
endif()
