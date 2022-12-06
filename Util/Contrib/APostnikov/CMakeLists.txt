add_library(aux_contrib_ap OBJECT
       makebox.f  fillbox.f inver3.f  hit.f 
       opnout.f  test_ani.f test_md.f test_mdc.f wraxsf1.f wraxsf2.f
       read_xv.f  intpl04.f displa.f read_ev.f itochar.f w_arrow.f
       w_movie.f
)

#
add_executable(eig2bxsf
   eig2bxsf.f
   $<TARGET_OBJECTS:aux_contrib_ap>
)

add_executable(xv2xsf
   xv2xsf.f
   $<TARGET_OBJECTS:aux_contrib_ap>
)

add_executable(md2axsf
   md2axsf.f
   $<TARGET_OBJECTS:aux_contrib_ap>
)

add_executable(rho2xsf
   rho2xsf.f
   $<TARGET_OBJECTS:aux_contrib_ap>
)

add_executable(vib2xsf
   vib2xsf.f
   $<TARGET_OBJECTS:aux_contrib_ap>
)

add_executable(fmpdos fmpdos.f)


install(TARGETS
        eig2bxsf xv2xsf md2axsf rho2xsf vib2xsf fmpdos
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} )
