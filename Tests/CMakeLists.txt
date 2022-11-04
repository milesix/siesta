file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/install_tests")

if(WITH_LUA)
 file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/install_tests/lua-h2o")

 file(COPY lua_h2o/lua_h2o.fdf lua_h2o/siesta.lua DESTINATION "${CMAKE_BINARY_DIR}/install_tests/lua-h2o")
 file(COPY Pseudos/O.psf Pseudos/H.psf DESTINATION "${CMAKE_BINARY_DIR}/install_tests/lua-h2o")

 add_test( NAME    siesta-lua-h2o
           COMMAND siesta-siesta -out OUT lua_h2o.fdf 
           WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/install_tests/lua-h2o" )
endif()
