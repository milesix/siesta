set(top_srcdir "${CMAKE_SOURCE_DIR}/Src")

set(sources
    ts2ts.f90 )

list(APPEND sources

   ${top_srcdir}/precision.F
   ${top_srcdir}/parallel.F
   ${top_srcdir}/alloc.F90
   ${top_srcdir}/units.f90
   ${top_srcdir}/m_io.f
)

add_executable(ts2ts ${sources})
target_link_libraries(ts2ts
          PRIVATE
	  ${PROJECT_NAME}-libsys
	  libfdf::libfdf)

install(TARGETS ts2ts
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} )
