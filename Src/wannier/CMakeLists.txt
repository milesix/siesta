set(dir "${CMAKE_CURRENT_SOURCE_DIR}")

list(
  APPEND libsiesta_sources
  
  "${dir}/amn.F90"
  "${dir}/broadcast_projections.F"
  "${dir}/delk.F90"
  "${dir}/m_chempotwann.F90"
  "${dir}/m_wannier_in_nao.F90"
  "${dir}/m_w90_in_siesta.F90"
  "${dir}/siesta2wannier90.F90"
  "${dir}/wannier90_types.f"
  "${dir}/w90_in_types.f"
  "${dir}/write_inp_wannier.F90"
  "${dir}/compute_pw_matrix.F90"
  "${dir}/diagonalizeHk.F90"
  "${dir}/m_digest_nnkp.F90"
  "${dir}/m_noccbands.f"
  "${dir}/m_orderbands.F"
  "${dir}/m_overkkneig.F90"
  "${dir}/m_planewavematrix.F90"
  "${dir}/m_planewavematrixvar.F90"
  "${dir}/m_switch_local_proj.F90"
  "${dir}/m_trialorbitalclass.f90"
  "${dir}/m_writedelk.F"
  "${dir}/mmn.F90"
  "${dir}/reordpsi.F"

)

set(libsiesta_sources "${libsiesta_sources}" PARENT_SCOPE)
