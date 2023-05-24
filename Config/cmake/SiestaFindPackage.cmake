# This file is a modification of the one with the
# name "sdftd3_find_package" in s-dft3. The original terms follow:
# -----------------------------------------
# This file is part of s-dftd3.
# SPDX-Identifier: LGPL-3.0-or-later
#
# s-dftd3 is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# s-dftd3 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with s-dftd3.  If not, see <https://www.gnu.org/licenses/>.

# This changed version streamlines the effects and uses a function rather
# than a macro.
# It enables a bit more things and ensures expressivity in the options.

# Handling of subproject dependencies
function(Siesta_find_package)
  # Define options
  set(options QUIET)
  set(oneValueArgs
    # the package name:
    # Used for variable lookups {NAME}_FIND_METHODS etc.
    NAME
    URL REVISION # for git repo's
    SUBMODULE_DIR)
  set(multiValueArgs)
  cmake_parse_arguments(_f "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

  # check arguments
  if(NOT DEFINED _f_NAME)
    # Check that the package as an unnamed argument
    list(LENGTH _f_UNPARSED_ARGUMENTS nargs)
    if(nargs GREATER 0)
      list(POP_FRONT _f_UNPARSED_ARGUMENTS _f_NAME)
    else()
      message(FATAL_ERROR "Siesta_find_package missing argument: NAME")
    endif()
  endif()

  set(pkg "${_f_NAME}")
  # convert package name to lower and upper case
  string(TOLOWER "${pkg}" pkg_lc)
  string(TOUPPER "${pkg}" pkg_uc)

  macro(mymsg check msg)
    if(NOT _f_QUIET)
      message(${check} "${msg}")
    endif()
  endmacro()

  mymsg(CHECK_START "Searching for ${pkg}")
  # always add space (for sub-stuff)
  list(APPEND CMAKE_MESSAGE_INDENT "| ")

  # Allow users to pass options via flags
  #  *_URL
  #  *_REVISION
  #  *_SOURCE_DIR
  foreach(n IN ITEMS ${pkg_lc} ${pkg_uc} ${pkg})

    # Also allow externally setting the URL and REVISION
    if(DEFINED "${n}_URL")
      set(_f_URL "${${n}_URL}")
    endif()

    if(DEFINED "${n}_REVISION")
      set(_f_REVISION "${${n}_REVISION}")
    endif()

    # The submodule directory is a SOURCE_DIR from the outside
    # world. Hence it should look up the SOURCE_DIR
    if(DEFINED "${n}_SOURCE_DIR")
      set(_f_SUBMODULE_DIR "${${n}_SOURCE_DIR}")
    endif()
  endforeach()

  # Determine the allowed find-methods
  set(all_f_methods "cmake" "pkgconf" "subproject" "submodule" "fetch")
  set(allowed_f_methods "cmake" "pkgconf")
  if(DEFINED _f_SUBMODULE_DIR)
    list(APPEND allowed_f_methods "subproject" "submodule")
  endif()
  if(DEFINED _f_URL AND DEFINED _f_REVISION)
    list(APPEND allowed_f_methods "fetch")
  elseif(DEFINED _f_URL)
    message(WARNING "Siesta_find_package: ${pkg} will not allow fetch as the URL is supplied, but no revision, please use -D${pkg}_REVISION= to specify revision")
  elseif(DEFINED _f_REVISION)
    message(WARNING "Siesta_find_package: ${pkg} will not allow fetch as the REVISION is supplied, but no URL, please use -D${pkg}_URL= to specify URL")
  endif()

  # default the searched methods
  set(f_methods "${allowed_f_methods}")
  if(DEFINED "${PROJECT_NAME}_FIND_METHOD")
    set(f_methods "${${PROJECT_NAME}_FIND_METHOD}")
  endif()

  foreach(n IN ITEMS ${pkg_lc} ${pkg_uc} ${pkg})
    # Parse the option of the *find-methods*
    if(DEFINED "${n}_FIND_METHOD")
      set(f_methods "${${n}_FIND_METHOD}")
    endif()
  endforeach()

  # Strip away elements that are not in the allowed-f-methods
  # If there is no submodule, it shouldn't be allowed.
  list(REMOVE_DUPLICATES f_methods)
  # lower
  string(TOLOWER "${f_methods}" f_methods)
  foreach(method IN LISTS f_methods)
    if(NOT method IN_LIST allowed_f_methods)
      list(REMOVE_ITEM f_methods ${method})
    endif()
  endforeach()
  if("subproject" IN_LIST f_methods AND "submodule" IN_LIST f_methods)
    list(REMOVE_ITEM f_methods "subproject")
  endif()

  if(TARGET "${pkg}::${pkg}")
    # skip all methods
    set(f_methods)
  endif()

  #debug
  #message("METHODS = ${f_methods} | ${allowed_f_methods}")

  foreach(method IN ITEMS ${f_methods})
    # assert that method is in all_f_methods
    if(NOT method IN_LIST all_f_methods)
      message(FATAL_ERROR "Siesta_find_package: find-method ${method} is not in allowed methods ${all_f_methods}")
    endif()

    if("${method}" STREQUAL "cmake")
      mymsg(CHECK_START "CMake package lookup")

      find_package("${pkg}" CONFIG)
      if("${pkg}_FOUND")
        mymsg(CHECK_PASS "found")
        break()
      elseif(NOT _f_QUIET)
        mymsg(CHECK_FAIL "not found")
      endif()
    endif()

    if("${method}" STREQUAL "pkgconf")
      find_package(PkgConfig QUIET REQUIRED)

      mymsg(CHECK_START "pkg-config package lookup")

      pkg_check_modules("${pkg_uc}" QUIET "${pkg}")
      if("${${pkg_uc}_FOUND}")
        mymsg(CHECK_PASS "found")

        add_library("${pkg}::${pkg}" INTERFACE IMPORTED GLOBAL)
        target_link_libraries(
          "${pkg}::${pkg}"
          INTERFACE
          "${${pkg_uc}_LINK_LIBRARIES}"
        )
        target_include_directories(
          "${pkg}::${pkg}"
          INTERFACE
          "${${pkg_uc}_INCLUDE_DIRS}"
        )
        break()
      else()
        mymsg(CHECK_FAIL "not found")
      endif()
    endif()

    if("${method}" STREQUAL "subproject" OR "${method}" STREQUAL "submodule")
      mymsg(CHECK_START "sub-project in folder: ${_f_SUBMODULE_DIR}")

      set("${pkg_uc}_SOURCE_DIR" "${_f_SUBMODULE_DIR}")
      if(EXISTS "${${pkg_uc}_SOURCE_DIR}/CMakeLists.txt")
        mymsg(CHECK_PASS "found")

        set("${pkg_uc}_BINARY_DIR" "${_f_SUBMODULE_DIR}")
        add_subdirectory(
          "${${pkg_uc}_SOURCE_DIR}"
          "${${pkg_uc}_BINARY_DIR}"
        )

        add_library("${pkg}::${pkg}" INTERFACE IMPORTED GLOBAL)
        target_link_libraries("${pkg}::${pkg}" INTERFACE "${pkg}")

        # We need the module directory in the subproject before we finish the configure stage
        if(NOT EXISTS "${${pkg_uc}_BINARY_DIR}/include")
          make_directory("${${pkg_uc}_BINARY_DIR}/include")
        endif()

        break()
      else()
        mymsg(CHECK_FAIL "not found")
      endif()
    endif()

    if("${method}" STREQUAL "fetch")
      mymsg(CHECK_START "fetching from ${_f_URL} (${pkg})")

      include(FetchContent)
      FetchContent_Declare(
        "${pkg_lc}"
        GIT_REPOSITORY "${_f_URL}"
        GIT_TAG "${_f_REVISION}"
      )
      FetchContent_MakeAvailable("${pkg_lc}")

      add_library("${pkg}::${pkg}" INTERFACE IMPORTED GLOBAL)
      target_link_libraries("${pkg}::${pkg}" INTERFACE "${pkg}")

      # We need the module directory in the subproject before we finish the configure stage
      FetchContent_GetProperties("${pkg_lc}" BINARY_DIR "${pkg_uc}_BINARY_DIR")

      mymsg(CHECK_PASS "fetched")

      if(NOT EXISTS "${${pkg_uc}_BINARY_DIR}/include")
        make_directory("${${pkg_uc}_BINARY_DIR}/include")
      endif()

      break()
    endif()

  endforeach()

  list(POP_BACK CMAKE_MESSAGE_INDENT)
  if(NOT TARGET "${pkg}::${pkg}")
    mymsg(CHECK_FAIL "not found")
    message(FATAL_ERROR "Could not find dependency ${pkg}")
  endif()

  mymsg(CHECK_PASS "found")

  # notify about the found-variables
  set("${pkg}_FOUND" TRUE PARENT_SCOPE)
  set("${pkg_uc}_FOUND" TRUE PARENT_SCOPE)

endfunction()

