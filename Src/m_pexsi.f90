module m_pexsi
use precision, only: dp
use iso_c_binding

integer(c_intptr_t), public :: plan

public :: pexsi_initialize_scfloop
public :: pexsi_finalize_scfloop

private

CONTAINS

subroutine pexsi_initialize_scfloop(World_Comm,npPerPole,mpirank)
  use f_ppexsi_interface
  integer, intent(in) :: npPerPole, mpirank
  integer, intent(in) :: World_Comm

  integer :: numProcRow, numProcCol
  integer :: outputFileIndex, info

  call get_row_col(npPerPole,numProcRow,numProcCol)

  ! Set the outputFileIndex to be the pole index.
  ! Starting from PEXSI v0.8.0, the first processor for each pole outputs
  ! information

  if( mod( mpirank, npPerPole ) .eq. 0 ) then
    outputFileIndex = mpirank / npPerPole;
  else
    outputFileIndex = -1;
  endif

  plan = f_ppexsi_plan_initialize(&
    World_Comm,&
    numProcRow,&
    numProcCol,&
    outputFileIndex,&
    info) 

  if (mpirank == 0) then
    print *, "Info in plan_initialize: ", info
  endif
end subroutine pexsi_initialize_scfloop

subroutine pexsi_finalize_scfloop() ! mpirank)
  use f_ppexsi_interface
  !  integer, intent(in) :: mpirank

  integer :: info

  call f_ppexsi_plan_finalize( plan, info )

  !  if (mpirank == 0) then
  !     print *, "Info in plan_finalize: ", info
  !  endif
end subroutine pexsi_finalize_scfloop

subroutine get_row_col(np,nrow,ncol)
  integer, intent(in)  :: np
  integer, intent(out) :: nrow, ncol
  !
  ! Finds the factors nrow and ncol such that nrow*ncol=np,
  ! are as similar as possible, and nrow>=ncol.
  ! For prime np, ncol=1, nrow=np.

  ncol  = floor(sqrt(dble(np)))
  do
    nrow = np/ncol
    if (nrow*ncol == np) exit
    ncol = ncol - 1
  enddo
end subroutine get_row_col

end module m_pexsi
