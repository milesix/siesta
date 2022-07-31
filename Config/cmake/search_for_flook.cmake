#
# We rely on the environmental variable FLOOK_ROOT
#
 if("$ENV{FLOOK_ROOT}" STREQUAL "")
   message(SEND_ERROR "Cannot find FLOOK_ROOT")
 else()
  set(FLOOK_ROOT $ENV{FLOOK_ROOT})
  message("Found FLOOK_ROOT: ${FLOOK_ROOT}")
  set(FLOOK_LIBS "-L${FLOOK_ROOT}/lib -lflookall -ldl")
  set(FLOOK_INCLUDE_DIRS "${FLOOK_ROOT}/include")
  message("FLOOK_LIBS: ${FLOOK_LIBS}")
  message("FLOOK_INCLUDE_DIRS: ${FLOOK_INCLUDE_DIRS}")
 endif()
