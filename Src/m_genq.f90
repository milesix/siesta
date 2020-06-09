!
module m_genq

  use precision, only: dp
  use sys,       only: die

  implicit none

  type, public :: genq_t
     integer                        :: nvars = 0
     real(dp), allocatable          :: q(:)
     character(len=32), allocatable :: qname(:)
  end type genq_t

  public :: geom_from_genq
  public :: cartF_to_genqF
  public :: read_genq

  type(genq_t), public, save  :: genq

  private

  CONTAINS

  subroutine geom_from_genq(M,q,xa,spec_no)
    ! Given the values of M generalized coordinates
    ! q(1), q(2), ... q(M)
    ! this routine specifies the values of the crystal coordinates xa

    ! It might be possible to generate this file automatically from
    ! some suitable symbolic representation of the coordinates
    ! This is just a proof of concept

    integer, parameter    :: p = selected_real_kind(10,100)

    integer, intent(in)                         :: M
    real(p), intent(in), target                 :: q(M)
    real(p), target, intent(inout)              :: xa(:,:)
    integer, intent(inout), optional            :: spec_no(:)


    integer :: natoms

    ! Make sure you include the right file
    include "geom.inc"
    
  end subroutine geom_from_genq

  subroutine cartF_to_genqF(M,q,ucell,na,fa,gF,ftol,epsgF)

    ! Given the values of M generalized coordinates
    ! q(1), q(2), ... q(M)
    ! and a routine geom_from_genq
    ! this routine computes \partial xa \partial q
    ! and the projection of the cartesian forces on the gen coords

    ! We assume linear dependencies
    ! In most (all?) cases, the derivatives will be simple integers

    Integer, intent(in)                    :: M
    real(dp), intent(in)                   :: q(M)
    real(dp), intent(in)                   :: ucell(3,3)
    integer, intent(in)                    :: na
    real(dp), intent(in)                   :: fa(3,na)
    real(dp), intent(out)                  :: gF(M)
    real(dp), intent(in)                   :: ftol
    real(dp), intent(out)                  :: epsgF(M)


    real(dp), allocatable           :: xa(:,:)
    real(dp), allocatable           :: xaf(:,:)
    real(dp), allocatable           :: dxdq(:,:)
    real(dp), allocatable           :: qq(:)

    real(dp)       :: h  = 0.01
    real(dp)       :: dxdq_cart(3), epsf(3)

    integer       :: ia, i, j

    allocate(xa(3,na))
    allocate(xaf(3,na))
    allocate(dxdq(3,na))
    allocate(qq(M))

    ! Compute partial derivatives by simple rule

    call geom_from_genq(M,q,xa)

    qq = q

    do i = 1, M
       gF(i) = 0.0_dp
       epsgF(i) = 0.0_dp
       !
       qq(i) = q(i) + h
       call geom_from_genq(M,qq,xaf)
!       dxdq(1:3,1:na) = nint((xaf-xa) / h)
       dxdq(1:3,1:na) = (xaf-xa) / h
       !
       !  Project the cartesian forces to the generalized coords
       !  to get the "generalized forces"
       !
       do ia = 1, na
          dxdq_cart(:) = matmul(ucell,dxdq(:,ia))
          gF(i) = gF(i) + dot_product(fa(:,ia),dxdq_cart(:))
          ! Estimate the force when almost converged
          do j=1,3
             epsf(j) = sign(ftol,fa(j,ia))
          enddo
          ! Estimate a tolerance for gF(i)
          epsgF(i) = epsgF(i) + dot_product(epsf(:),dxdq_cart(:))
       enddo
       !
       qq(i) = q(i)
    enddo

    deallocate(xa,xaf,dxdq,qq)

  end subroutine cartF_to_genqF

  subroutine read_genq(genq)
    ! For now, this routine reads from file GENQ

    use m_mpi_utils, only:  broadcast
    use parallel,    only:  Node

    type(genq_t) :: genq
    
    integer :: var_u, iostat, i

    if (Node .eq. 0) then
       open(unit=var_u,file="GENQ",form="formatted",status="old", &
                  position="rewind",action="read")

       read(var_u, fmt=*, iostat=iostat) genq%nvars
       if (iostat /= 0) then
          call die("Cannot read nvars from GENQ")
       endif
    endif

    call broadcast(genq%nvars)

    allocate(genq%q(genq%nvars),genq%qname(genq%nvars))
  
    if (Node .eq. 0) then
       do i = 1, genq%nvars
          read(var_u, fmt=*, iostat=iostat) genq%q(i), genq%qname(i)
          if (iostat /= 0) then
             call die("ERROR while reading q")
          endif
       enddo
       close(var_u)
    endif

    call broadcast(genq%q)
    !! No support yet.. call broadcast(genq%qname)

end subroutine read_genq

end module m_genq

    


