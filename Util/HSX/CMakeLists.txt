#
# Go the extra mile and treat this as a library
# that users can link to in their own programs.
#
add_library(hsx hsx_m.f90)
set_target_properties(hsx
  PROPERTIES
  OUTPUT_NAME "hsx"
  VERSION "${PROJECT_VERSION}"
  Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/include"
)

target_include_directories(
  hsx
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

install(
  TARGETS hsx
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  )
install(
  DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/include/"
  DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
)

#
# Now the utilities
# 
add_executable( hsx2hs hsx2hs.f90)
add_executable( hs2hsx hs2hsx.f90)

target_link_libraries(hsx2hs hsx)
target_link_libraries(hs2hsx hsx)

install(
  TARGETS hsx2hs hs2hsx
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  )

