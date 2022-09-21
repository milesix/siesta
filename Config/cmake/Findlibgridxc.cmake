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

set(_lib "libgridxc")
set(_pkg "LIBGRIDXC")
set(_url "https://gitlab.com/siesta-project/libraries/libgridxc")
set(_tag "cmake-master")

# For now, allow only compilation as a subproject or by fetching,
# or as a pre-compiled cmake-ready package
set(LIBGRIDXC_FIND_METHOD "cmake" "subproject" "fetch")

if(NOT DEFINED "${_pkg}_FIND_METHOD")
  if(DEFINED "${PROJECT_NAME}-dependency-method")
    set("${_pkg}_FIND_METHOD" "${${PROJECT_NAME}-dependency-method}")
  else()
    set("${_pkg}_FIND_METHOD" "cmake" "pkgconf" "subproject" "fetch")
  endif()
  set("_${_pkg}_FIND_METHOD")
endif()

include("${CMAKE_CURRENT_LIST_DIR}/s-dftd3-utils.cmake")

sdftd3_find_package("${_lib}" "${${_pkg}_FIND_METHOD}" "${_url}" "${_tag}")

if (WITH_LIBXC)
  target_link_libraries(libgridxc::libgridxc INTERFACE libxc::XC_Fortran)
endif()

if(DEFINED "_${_pkg}_FIND_METHOD")
  unset("${_pkg}_FIND_METHOD")
  unset("_${_pkg}_FIND_METHOD")
endif()
unset(_lib)
unset(_pkg)
unset(_url)
unset(_tag)
