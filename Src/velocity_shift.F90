! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!

!< Module implementing an eigenspectrum shift according to the group velocity of an
!< eigenstate.
!<
!< This is developed in two parts:
!< 1. The diagkp_velocity code implements the diagonalization routines + the
!<    velocity calculations.
!<    The velocities are calculated using:
!<      < psi_i | dH_kx/ky/kz - epsilon_i dS_kx/ky/kz | psi_i >
!< 2. A shift of the eigenvalue epsilon_i depending on the applied bias and direction.
!<    The group velocity is projected onto a "field" direction and if the scalar projection
!<    is positive (> 0) the eigenspectrum is shifted V / 2 down, to fill it, while negative
!<    projections are shifted V / 2 up, to empty it.
module m_velocity_shift

  use precision, only: dp

  implicit none

  save
  private

  !< Logical to control whether the band-shifts are being performed.
  logical :: use_velocity_shift = .false.

  !< Logical to control whether the current is calculated
  logical :: calc_velocity_current = .true.

  !< The bias applied to the bands
  real(dp) :: velocity_h_bias = 0._dp

  !< The direction along which the bands will be displaced.
  real(dp) :: velocity_dir(3) = 0._dp

  !< The projection factor required before the velocity-shift is used
  real(dp) :: velocity_tolerance = 0._dp

  public :: use_velocity_shift
  public :: calc_velocity_current
  public :: read_velocity_shift
  public :: velocity_shift
  public :: velocity_current

contains

  subroutine read_velocity_shift()
    use fdf
    use parallel, only: IONode
    use units, only: eV, Ang
    use siesta_cml
    use intrinsic_missing, only: VNORM

    type(block_fdf) :: bfdf
    type(parsed_line), pointer :: pline

    real(dp) :: velocity_bias

    real(dp), parameter :: hbar_Rys = 4.8377647940592375e-17_dp
    character(len=*), parameter :: form_afa = '("redata: ",a,t53,"= ",f10.4,a)'
    character(len=*), parameter :: form_ava = '("redata: ",a,t53,"= ",e10.4,a)'
    character(len=*), parameter :: form_avvv = '("redata: ",a,t53,"=",3(tr1,e10.4))'
    character(len=*), parameter :: form_al = '("redata: ",a,t53,"= ",l1)'

    ! Read in the velocity related options
    velocity_bias = fdf_get('BulkBias.V', 0._dp, 'Ry')

    ! The velocity tolerance in atomic units
    velocity_tolerance = fdf_get('BulkBias.Tolerance', 1.e-15_dp)

    ! TODO libfdf
    ! Once libfdf 0.1 is merged with trunk we can use fdf_islist for reals as well
!    if ( fdf_islist('BulkBias.Direction') ) then
!      ! In this case the BulkBias.V is determining the size of the field
!      call fdf_list('BulkBias.Direction', 3, velocity_dir)
    if ( fdf_block('BulkBias.Direction', bfdf) ) then
      ! Read first line
      if ( .not. fdf_bline(bfdf, pline) ) then
        call die('velocity_shift: could not read BulkBias.Direction block')
      end if

      ! Read in the field
      velocity_dir(1) = fdf_bvalues(pline, 1)
      velocity_dir(2) = fdf_bvalues(pline, 2)
      velocity_dir(3) = fdf_bvalues(pline, 3)

      call fdf_bclose(bfdf)
    else
      velocity_dir = 0._dp
    end if

    ! Determine whether we should use the bulkbias
    ! Only for applied bias above 1e-6 eV
    use_velocity_shift = abs(velocity_bias) > 1.e-6_dp * eV

    calc_velocity_current = fdf_get('BulkBias.Current', .true.)
    
    ! Immediately return if not used
    if ( .not. use_velocity_shift ) return

    
    if ( VNORM(velocity_dir) == 0._dp ) then
      call die('velocity_shift: BulkBias.Direction direction has not been specified!')
    end if

    ! Normalize the direction
    velocity_dir = velocity_dir / VNORM(velocity_dir)

    if ( IONode ) then
      write(*,form_afa) 'Bulk velocity bias ', velocity_bias / eV, ' Volts'
      write(*,form_avvv) 'Bulk velocity bias direction', velocity_dir
      write(*,form_ava) 'Bulk velocity bias projection tolerance', velocity_tolerance * 1.e-12_dp / Ang / hbar_Rys, ' Ang/ps'
      write(*,form_al) 'Bulk velocity bias current calculation', calc_velocity_current
    end if

    if ( cml_p ) then
      call cmlStartPropertyList(xf=mainXML, title='BulkBias', &
          dictRef='siesta:BulkBias')
      call cmlAddProperty( xf=mainXML, value=velocity_bias/eV, &
          dictRef='siesta:BulkBias.V', units='siestaUnits:eV')
      call cmlAddProperty( xf=mainXML, value=velocity_tolerance * 1.e-12_dp / Ang / hbar_Rys, &
          dictRef='siesta:BulkBias.Tolerance', units='siestaUnits:Ang/ps')
      call cmlAddProperty( xf=mainXML, value=velocity_dir, &
          dictRef='siesta:BulkBias.Direction')
      call cmlEndPropertyList(mainXML)
    end if

    ! Only take half the bias (so we dont have to divide by 2 all the time)
    velocity_h_bias = velocity_bias * 0.5_dp

  end subroutine read_velocity_shift

  !< Shift the eigenspectrum according to the velocities
  subroutine velocity_shift(sign, ne, e, v)
    integer, intent(in) :: sign
    integer, intent(in) :: ne
    real(dp), intent(inout) :: e(ne)
    real(dp), intent(in) :: v(3, ne)

    integer :: ie
    real(dp) :: p

    if ( sign > 0 ) then

      do ie = 1, ne

        p = dot_product(velocity_dir, v(:,ie))
        if ( p > velocity_tolerance ) then
          ! The velocity is along the direction of the potential drop
          ! This means that the electron will be pushed in the same direction and filled
          ! up faster.
          e(ie) = e(ie) - velocity_h_bias
        else if ( p < - velocity_tolerance ) then
          e(ie) = e(ie) + velocity_h_bias
        end if
        
      end do

    else

      ! Shift back the eigenvalues
      
      do ie = 1, ne
        p = dot_product(velocity_dir, v(:,ie))
        if ( p > velocity_tolerance ) then
          e(ie) = e(ie) + velocity_h_bias
        else if ( p < - velocity_tolerance ) then
          e(ie) = e(ie) - velocity_h_bias
        end if
      end do
      
    end if

  end subroutine velocity_shift

  !< Calculate current from an eigenspectrum according to the velocities and calculate current
  subroutine velocity_current(ne, e, v, Ef, w, Temp, I)
    integer, intent(in) :: ne
    ! The eigenvalues `e` *must* be un-shifted
    real(dp), intent(in) :: e(ne), v(3,ne), Ef, w, Temp
    real(dp), intent(inout) :: I

    integer :: ie
    real(dp) :: p, lI

    lI = 0._dp
    do ie = 1, ne
      
      p = dot_product(velocity_dir, v(:,ie))
      ! this is velocity projected onto bias direction
      if ( p > velocity_tolerance ) then
        lI = lI + p * dnf(e(ie), Ef, velocity_h_bias, Temp)
      end if
      
    end do

    I = I + lI * w

  contains
    
    pure function dnf(E, Ef, Vh, T) result(nf)
      real(dp), intent(in) :: E, Ef, Vh, T
      real(dp) :: nf, x1, x2

      x1 = (E - Ef - Vh) / T
      x2 = (E - Ef + Vh) / T

      if ( x1 > 100._dp ) then
        x1 = 0._dp
      else if ( x1 < -100._dp ) then
        x1 = 1._dp
      else
        x1 = 1._dp / ( 1._dp + exp(x1) )
      end if

      if ( x2 > 100._dp ) then
        x2 = 0._dp
      else if ( x2 < -100._dp ) then
        x2 = 1._dp
      else
        x2 = 1._dp / ( 1._dp + exp(x2) )
      end if

      nf = x1 - x2
      !nf = 1._dp / (1._dp + exp((E - Ef - Vh) / T)) - &
      !    1._dp / (1._dp + exp((E - Ef + Vh) / T))
      
    end function dnf
    
  end subroutine velocity_current

end module m_velocity_shift
