#
# Two simple tests to check the installation
# They re-use files from the Examples directory, copying them to
# the build directory to avoid polluting the source tree.
#
# The same basic mechanism could be used with other material (e.g.
# that in Tests/)
#
file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/install_tests")
file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/install_tests/MgO")
file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/install_tests/H2O")

file(COPY MgO/MgO.fdf DESTINATION "${CMAKE_BINARY_DIR}/install_tests/MgO")
file(COPY Vps/O.psf Vps/Mg.psf DESTINATION "${CMAKE_BINARY_DIR}/install_tests/MgO")

file(COPY H2O/h2o.fdf DESTINATION "${CMAKE_BINARY_DIR}/install_tests/H2O")
file(COPY Vps/O.psf Vps/H.psf DESTINATION "${CMAKE_BINARY_DIR}/install_tests/H2O")

add_test(  NAME    siesta-mgo_test
           COMMAND siesta-siesta -out OUT MgO.fdf
           WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/install_tests/MgO" )

add_test(  NAME    siesta-h2o_test
           COMMAND siesta-siesta -out OUT h2o.fdf
           WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/install_tests/H2O" )
