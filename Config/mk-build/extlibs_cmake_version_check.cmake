if(${CMAKE_VERSION} VERSION_LESS "3.17.0")
    message( FATAL_ERROR "On-the-fly dependency compilation needs CMake >= 3.17")
endif()
