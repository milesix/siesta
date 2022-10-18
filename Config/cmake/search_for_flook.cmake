#
# We rely on the environmental variable FLOOK_ROOT
#
set(FLOOK_ROOT "$ENV{FLOOK_ROOT}" CACHE FILEPATH "flook installation path")

 if("${FLOOK_ROOT}" STREQUAL "" )

    message(SEND_ERROR "Variable FLOOK_ROOT not set")
     
 else()
 
  message("Using flook installation in ${FLOOK_ROOT}")
  set(FLOOK_LIBS "-L${FLOOK_ROOT}/lib -lflookall -ldl")
  set(FLOOK_INCLUDE_DIRS "${FLOOK_ROOT}/include")
  message("FLOOK_LIBS: ${FLOOK_LIBS}")
  message("FLOOK_INCLUDE_DIRS: ${FLOOK_INCLUDE_DIRS}")

  add_library(flook::flook INTERFACE IMPORTED GLOBAL)
  target_link_libraries(flook::flook
               INTERFACE "${FLOOK_LIBS}")
  target_include_directories(flook::flook
               INTERFACE
              "${FLOOK_INCLUDE_DIRS}")

 endif()
