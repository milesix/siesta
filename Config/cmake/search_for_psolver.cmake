#
#
# -- Search for Psolver using pkg-config data
# 
 find_package(PkgConfig REQUIRED)
 pkg_check_modules(PSOLVER REQUIRED psolver)

 message(DEBUG "PSOLVER Libdir: ${PSOLVER_LIBDIR}")
 message(DEBUG "Psolver Libraries: ${PSOLVER_LIBRARIES}")
 message(STATUS "Psolver Link libraries: ${PSOLVER_LINK_LIBRARIES}")

 message(STATUS " ----- ")
 message(STATUS "   ** Make sure that there are no stray linear algebra libraries **")
 message(STATUS "   ** in the above list of libraries                             **")
 message(STATUS " ----- ")
 
 message(DEBUG "PSOLVER_INCLUDEDIR: ${PSOLVER_INCLUDEDIR}")
 message(STATUS "PSOLVER_INCLUDE_DIRS: ${PSOLVER_INCLUDE_DIRS}")

add_library(psolver::psolver INTERFACE IMPORTED GLOBAL)
 target_link_libraries(psolver::psolver INTERFACE  ${PSOLVER_LINK_LIBRARIES})
 target_include_directories(psolver::psolver INTERFACE ${PSOLVER_INCLUDE_DIRS})
