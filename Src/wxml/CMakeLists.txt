set(dir "${CMAKE_CURRENT_SOURCE_DIR}")

list(
  APPEND srcs

  "${dir}/flib_wcml.f90"
  "${dir}/flib_wstml.f90"
  "${dir}/flib_wxml.f90"
  "${dir}/m_wcml_coma.f90"
  "${dir}/m_wxml_array_str.f90"
  "${dir}/m_wxml_buffer.f90"
  "${dir}/m_wxml_core.f90"
  "${dir}/m_wxml_dictionary.f90"
  "${dir}/m_wxml_elstack.f90"
  "${dir}/m_wxml_error.f90"
  "${dir}/m_wxml_escape.f90"
  "${dir}/m_wxml_overloads.f90"
  "${dir}/m_wxml_text.F90"

)

set(srcs "${srcs}" PARENT_SCOPE)
