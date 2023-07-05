add_library(
  ${PROJECT_NAME}-libunits
  STATIC
  units.F90
  units_legacy.f90
  units_codata2018.f90
  units_common.f90
  )

# This is duplicated in ${PROJECT_SOURCE_DIR}/CMakeLists.txt
# The reason is that we might remove the legacy code in all codes
# and fully resort to only handling this in units* files in this library.
# That would be the longterm goal.

# Lower case the units
string(TOLOWER "${WITH_UNITS_CONVENTION}" siesta_units_convention)
if("${siesta_units_convention}" STREQUAL "original" OR
   "${siesta_units_convention}" STREQUAL "legacy")
  target_compile_definitions(${PROJECT_NAME}-units
    PRIVATE
    SIESTA__UNITS_ORIGINAL)
elseif("${siesta_units_convention}" STREQUAL "codata2018")
  # pass - do nothing - the default
else()
  message(FATAL_ERROR "Unknown argument for WITH_UNITS_CONVENTION; not one of [original|legacy,codata2018]")
endif()

# Ensure module files are added to the inclusion
target_include_directories(
  ${PROJECT_NAME}-libunits
  INTERFACE "${CMAKE_CURRENT_BINARY_DIR}")
