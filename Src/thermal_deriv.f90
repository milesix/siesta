module derivation_routines
  !! The only way I could impement so far finite-difference derivation dispatch.
  !! Not the mostelegant way, probably, but works.
  !!
  !! Provides `apply_derivation' interface that works on 2, 3 and 4-dimentional
  !! arrays, with possible choice of 2, 3-midpoint and 5-midpoint derivative
  !! schemes. The last dimention of array `F' selects the time point, where
  !! index 1 stands for t-0 ( means F(..,1) ) where index 2 will contain the
  !! derivative after the function call ( means F(..,2) ).

  use precision, only: dp, grid_p
  implicit none

interface apply_derivation

   module procedure apply_derivation_dim2
   module procedure apply_derivation_dim3
   module procedure apply_derivation_dim4

end interface apply_derivation

contains

  subroutine apply_derivation_dim2(F, h, scheme)
    real(dp), intent(inout) :: F(:,:)
    real(dp), intent(in) :: h
    integer, intent(in)  :: scheme

    select case(scheme)
    case(2)
       call derivation_2pts_dim2(F, h)
    case(3)
       call derivation_3pts_mid_dim2(F, h)
    case(5)
       call derivation_5pts_mid_dim2(F, h)
    case default
       call die('Non-existent derivation scheme requested, aborting.')
    end select
  end subroutine apply_derivation_dim2

  subroutine apply_derivation_dim3(F, h, scheme)
    real(dp), intent(inout) :: F(:,:,:)
    real(dp), intent(in) :: h
    integer, intent(in)  :: scheme

    select case(scheme)
    case(2)
       call derivation_2pts_dim3(F, h)
    case(3)
       call derivation_3pts_mid_dim3(F, h)
    case(5)
       call derivation_5pts_mid_dim3(F, h)
    case default
       call die('Non-existent derivation scheme requested, aborting.')
    end select
  end subroutine apply_derivation_dim3

  subroutine apply_derivation_dim4(F, h, scheme)
    real(dp), intent(inout) :: F(:,:,:,:)
    real(dp), intent(in) :: h
    integer, intent(in)  :: scheme

    select case(scheme)
    case(2)
       call derivation_2pts_dim4(F, h)
    case(3)
       call derivation_3pts_mid_dim4(F, h)
    case(5)
       call derivation_5pts_mid_dim4(F, h)
    case default
       call die('Non-existent derivation scheme requested, aborting.')
    end select
  end subroutine apply_derivation_dim4


  subroutine derivation_2pts_dim2(F, h)
    real(dp), intent(inout) :: F(:,:)
    real(dp), intent(in) :: h

    F(:,2) = (F(:,2) - F(:,1)) / h
  end subroutine derivation_2pts_dim2

  subroutine derivation_3pts_mid_dim2(F, h)
    real(dp), intent(inout) :: F(:,:)
    real(dp), intent(in) :: h

    F(:,2) = (F(:,3) - F(:,2)) / (2.0_dp * h)
  end subroutine derivation_3pts_mid_dim2

  subroutine derivation_5pts_mid_dim2(F, h)
    real(dp), intent(inout) :: F(:,:)
    real(dp), intent(in) :: h

    F(:,2) = (F(:,2) - 8.0_dp * F(:,3) &
         & + 8.0_dp * F(:,4) - F(:,5)) / (12.0_dp * h)
  end subroutine derivation_5pts_mid_dim2

  subroutine derivation_2pts_dim3(F, h)
    real(dp), intent(inout) :: F(:,:,:)
    real(dp), intent(in) :: h

    F(:,:,2) = (F(:,:,2) - F(:,:,1)) / h
  end subroutine derivation_2pts_dim3

  subroutine derivation_3pts_mid_dim3(F, h)
    real(dp), intent(inout) :: F(:,:,:)
    real(dp), intent(in) :: h

    F(:,:,2) = (F(:,:,3) - F(:,:,2)) / (2.0_dp * h)
  end subroutine derivation_3pts_mid_dim3

  subroutine derivation_5pts_mid_dim3(F, h)
    real(dp), intent(inout) :: F(:,:,:)
    real(dp), intent(in) :: h

    F(:,:,2) = (F(:,:,2) - 8.0_dp * F(:,:,3) &
         & + 8.0_dp * F(:,:,4) - F(:,:,5)) / (12.0_dp * h)
  end subroutine derivation_5pts_mid_dim3

  subroutine derivation_2pts_dim4(F, h)
    real(dp), intent(inout) :: F(:,:,:,:)
    real(dp), intent(in) :: h

    F(:,:,:,2) = (F(:,:,:,2) - F(:,:,:,1)) / h
  end subroutine derivation_2pts_dim4

  subroutine derivation_3pts_mid_dim4(F, h)
    real(dp), intent(inout) :: F(:,:,:,:)
    real(dp), intent(in) :: h

    F(:,:,:,2) = (F(:,:,:,3) - F(:,:,:,2)) / (2.0_dp * h)
  end subroutine derivation_3pts_mid_dim4

  subroutine derivation_5pts_mid_dim4(F, h)
    real(dp), intent(inout) :: F(:,:,:,:)
    real(dp), intent(in) :: h

    F(:,:,:,2) = (F(:,:,:,2) - 8.0_dp * F(:,:,:,3) &
         & + 8.0_dp * F(:,:,:,4) - F(:,:,:,5)) / (12.0_dp * h)
  end subroutine derivation_5pts_mid_dim4

end module derivation_routines
