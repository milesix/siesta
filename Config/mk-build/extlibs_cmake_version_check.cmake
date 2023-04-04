if(${CMAKE_VERSION} VERSION_LESS "3.14.0")
    message( FATAL_ERROR "On-the-fly dependency compilation needs CMake >= 3.14")
endif()
