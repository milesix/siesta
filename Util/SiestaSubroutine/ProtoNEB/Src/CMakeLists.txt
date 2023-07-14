#
set(top_srcdir "${PROJECT_SOURCE_DIR}/Src")

if (WITH_MPI)
#
# ProtoNEB: Uses the MPI dispatch
#
add_executable(protoNEB
    protoNEB.F90
 
   ${top_srcdir}/reinit_m.F90
   ${top_srcdir}/siesta_cmlsubs.F90
   ${top_srcdir}/siesta_init.F
   ${top_srcdir}/m_io_yaml.F90
   ${top_srcdir}/siesta_end.F
   ${top_srcdir}/fsiesta_mpi.F90
   ${top_srcdir}/init_output.f90
)

target_link_libraries(protoNEB ${PROJECT_NAME}-libsiesta)

install(
  TARGETS protoNEB
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  )

endif(WITH_MPI)




  
