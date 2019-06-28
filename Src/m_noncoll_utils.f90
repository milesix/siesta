module m_noncoll_utils

  implicit none

  public :: flatten_noncoll
  public :: repack_noncoll

  
  integer, parameter, private :: dp = selected_real_kind(10,100)
  
CONTAINS

  subroutine flatten_noncoll(no_l, &
       numh,listhptr,maxnh,listh, &
       numh_u2,listhptr_u2,nnz_u2,listh_u2, &
       h_spin_dim, H, S, H_u2, S_u2)

    integer, intent(in) :: no_l
    integer, intent(in) :: numh(no_l), listhptr(no_l)
    integer, intent(in) :: maxnh
    integer, intent(in) :: listh(maxnh)
    integer, allocatable, intent(out) :: numh_u2(:), listhptr_u2(:)
    integer, allocatable, intent(out) :: listh_u2(:)
    integer, intent(in) :: h_spin_dim
    real(dp), intent(in) :: H(maxnh,h_spin_dim)
    real(dp), intent(in) :: S(maxnh)
    integer, intent(out) :: nnz_u2
    complex(dp), allocatable, intent(out) :: H_u2(:)
    complex(dp), allocatable, intent(out) :: S_u2(:)

    integer :: io, io2, jo, jo2, k, k2, ind, ind2
    
    allocate( numh_u2(2*no_l) )
    allocate( listhptr_u2(2*no_l) )

    do io = 1, no_l
       ! First rows of sub-blocks
       io2 = 2*io - 1
       numh_u2(io2) = 2 * numh(io)
       ! Second rows of sub-blocks
       io2 = 2*io
       numh_u2(io2) = 2 * numh(io)
    enddo
    ! Build listhptr_u2
    listhptr_u2(1) = 0
    do io2 = 2, 2*no_l
       listhptr_u2(io2) = listhptr_u2(io2-1) + numh_u2(io2-1)
    enddo

    nnz_u2 = sum(numh_u2)
    
    allocate( listh_u2(nnz_u2) )
    allocate( H_u2(nnz_u2), S_u2(nnz_u2) )

    ! Iterate over blocks
    do io = 1, no_l
       do k = 1, numh(io)
          ind = listhptr(io) + k
          jo = listh(ind)

          ! Place each element in its flattened place
          
          ! a_11
          io2 = 2*io - 1
          jo2 = 2*jo - 1
          k2 = 2*k-1
          ind2 = listhptr_u2(io2) + k2
          listh_u2(ind2) = jo2
          S_u2(ind2) = S(ind)
          H_u2(ind2) = cmplx( H(ind,1), H(ind,5), dp)
          
          ! a_12
          io2 = 2*io - 1
          jo2 = 2*jo
          k2 = 2*k
          ind2 = listhptr_u2(io2) + k2
          listh_u2(ind2) = jo2
          S_u2(ind2) = 0.0_dp
          H_u2(ind2) = cmplx( H(ind,3), -H(ind,4), dp)
          
          ! a_21
          io2 = 2*io 
          jo2 = 2*jo - 1
          k2 = 2*k-1
          ind2 = listhptr_u2(io2) + k2
          listh_u2(ind2) = jo2
          S_u2(ind2) = 0.0_dp
          H_u2(ind2) = cmplx( H(ind,7), H(ind,8), dp)

          ! a_22
          io2 = 2*io 
          jo2 = 2*jo 
          k2 = 2*k
          ind2 = listhptr_u2(io2) + k2
          listh_u2(ind2) = jo2
          S_u2(ind2) = S(ind)
          H_u2(ind2) = cmplx( H(ind,2), H(ind,6), dp)

       enddo
    enddo

  end subroutine flatten_noncoll
  
  subroutine repack_noncoll(no_l,  &
       numh,listhptr,maxnh, &
       numh_u2,listhptr_u2,nnz_u2, &
       h_spin_dim, DM_u2, EDM_u2, DM, EDM)

    integer, intent(in) :: no_l
    integer, intent(in) :: numh(no_l), listhptr(no_l)
    integer, intent(in) :: maxnh, nnz_u2
    integer, intent(in) :: numh_u2(:), listhptr_u2(:)
    integer, intent(in) :: h_spin_dim
    complex(dp), intent(in) :: DM_u2(:)
    complex(dp), intent(in) :: EDM_u2(:)
    real(dp), intent(out) :: DM(maxnh,h_spin_dim)
    real(dp), intent(out) :: EDM(maxnh,4)

    integer :: io, io2, k, k2, ind2, ind

    ! Iterate over blocks
    do io = 1, no_l
       do k = 1, numh(io)
          ind = listhptr(io) + k

          ! Get the pieces of each element from its flattened places
          
          ! a_11
          io2 = 2*io - 1
          k2 = 2*k-1
          ind2 = listhptr_u2(io2) + k2

          DM(ind,1) = real ( DM_u2(ind2) )
          DM(ind,5) = aimag ( DM_u2(ind2) )
          EDM(ind,1) = real ( EDM_u2(ind2) )
            !!! (ind2) = cmplx( H(ind,1), H(ind,5), dp)
          
          ! a_12
          io2 = 2*io - 1
          k2 = 2*k
          ind2 = listhptr_u2(io2) + k2

          DM(ind,3) = real ( DM_u2(ind2) )
          DM(ind,4) = -aimag ( DM_u2(ind2) )
          EDM(ind,3) = 0.5_dp * real ( EDM_u2(ind2) )
          EDM(ind,4) = -0.5_dp * aimag ( EDM_u2(ind2) )

           !!!   H_u2(ind2) = cmplx( H(ind,3), -H(ind,4), dp)
          
          ! a_21
          io2 = 2*io 
          k2 = 2*k-1
          ind2 = listhptr_u2(io2) + k2

          DM(ind,7) = real ( DM_u2(ind2) )
          DM(ind,8) = aimag ( DM_u2(ind2) )
          ! symmetrize just in case
          EDM(ind,3) = EDM(ind,3) + 0.5_dp * real ( EDM_u2(ind2) )
          EDM(ind,4) = EDM(ind,4) + 0.5_dp * aimag ( EDM_u2(ind2) )
     !!          H_u2(ind2) = cmplx( H(ind,7), H(ind,8), dp)

          ! a_22
          io2 = 2*io 
          k2 = 2*k
          ind2 = listhptr_u2(io2) + k2

          DM(ind,2) = real ( DM_u2(ind2) )
          DM(ind,6) = aimag ( DM_u2(ind2) )
          EDM(ind,2) =  real ( EDM_u2(ind2) )

          !!! H_u2(ind2) = cmplx( H(ind,2), H(ind,6), dp)

       enddo
    enddo

  end subroutine repack_noncoll
  
       
end module m_noncoll_utils
