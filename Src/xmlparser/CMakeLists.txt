set(dir "${CMAKE_CURRENT_SOURCE_DIR}")

list(
  APPEND srcs

  "${dir}/flib_sax.f90"
  "${dir}/m_buffer.f90"
  "${dir}/m_charset.f90"
  "${dir}/m_converters.f90"
  "${dir}/m_debug.f90"
  "${dir}/m_dictionary.f90"
  "${dir}/m_elstack.f90"
  "${dir}/m_entities.f90"
  "${dir}/m_fsm.f90"
  "${dir}/m_io_xmlparser.f90"
  "${dir}/m_reader.f90"
  "${dir}/m_xml_error.f90"
  "${dir}/m_xml_parser.f90"

)

set(srcs "${srcs}" PARENT_SCOPE)
