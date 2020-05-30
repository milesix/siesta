module ion_flux_data
  use precision, only: dp, grid_p
  use m_erf, only: erf_func, erfc_func
  use cellsubs, only : volcel

  implicit none
  SAVE

  real(dp)  :: pi  = 4.0_dp * atan(1.0_dp)

  real(dp)  :: ion_flux_Jion(3)
  !! Ionic component of the heat flux.

  real(dp)  :: ion_flux_a(3)
  real(dp)  :: ion_flux_b(3)
  real(dp)  :: ion_flux_c(3)
  real(dp)  :: ion_flux_d(3)
  real(dp)  :: ion_flux_e(3)

  ! Data for gvecs:
  real(dp), dimension(3,3) :: ion_flux_cell
  integer, dimension(3)    :: ion_flux_ntml
  integer, dimension(3)    :: ion_flux_mesh

  ! FIXME: those parameters must be read from input
  real(dp) :: eta = 0.1_dp
  !! Ewald factor for convergence
  integer  :: n_max = 5
  !! ?? cutoff per somme in griglia reale - doesn't look like

  real(dp) :: I_prime, I_prime_rec
  real(dp), allocatable :: I_uno_g(:,:,:)
  real(dp), allocatable :: I_dos_g(:)

contains

  real(kind=dp) function h_(x)
    real(kind=dp) :: x
    h_=erfc_func(x)/x
  end function h_


  real(kind=dp) function hp_(x)
    real(kind=dp) :: x
    hp_=-(2.d0/sqrt(pi))*(1.d0/x)*exp(-x*x)-1.d0/(x*x)*erfc_func(x)
  end function hp_


  real(kind=dp) function modulus(vector)
    real(kind=dp) ::vector(3)
    modulus=sqrt(vector(1)**2+vector(2)**2+vector(3)**2)
  end function modulus


  real(dp) function alpha_0_lr(l_eta, G_sq, flag)
    real(dp)  :: l_eta, G_sq
    integer   :: flag
    if (flag==1) then
       alpha_0_lr = exp(-G_sq/(4.d0*l_eta))/G_sq
    else
       alpha_0_lr = 1.d0/G_sq
    end if
  end function alpha_0_lr


  subroutine pbc(vin,vout)
    !! Apply the minimum image distance in a cubic cell.
    use fdf, only: fdf_physical
    use utils, only: die

    real(DP), intent(in) :: vin(3)
    real(DP), intent(out) :: vout(3)

    real :: alat
    integer :: i,n

    alat = fdf_physical('LatticeConstant',0.0_dp,'Bohr')
    if (alat==0.0_dp) call die('ion_flux:pbc', 'VMD requires alat set')

    vout(:)=0.d0
    do i=1,3

       if (vin(i)>=0) then
          n=int(vin(i)/alat)
       else
          n=int(vin(i)/alat)-1
       end if

       vout(i)=vin(i)-dble(n)*alat
    end do

  end subroutine pbc


  subroutine add_local_uno(val,pos,a,b)
    integer,  intent(in) :: a, b
    real(dp), intent(in) :: pos(3)
    real(dp), intent(inout) :: val
    real(dp) :: u(3), u_mod, n(3)
    integer  :: n_x, n_y, n_z
    real(dp)  :: val1

    val1 = 0.d0

    do n_x=-n_max,n_max
       do n_y=-n_max,n_max
          do n_z=-n_max,n_max
             n(1:3) = n_x * ion_flux_cell(1:3,1) +&
                  & n_y * ion_flux_cell(1:3,2) +&
                  & n_z * ion_flux_cell(1:3,3)
             u(1:3) = pos(1:3)-n(1:3)
             u_mod=modulus(u)
             val1=val1+eta*hp_(sqrt(eta)*u_mod)*u(a)*u(b)/u_mod
             if (a==b) then
                val1=val1+sqrt(eta)*h_(sqrt(eta)*u_mod)
             end if
          end do
       end do
    end do

    val = val + val1

  end subroutine add_local_uno


  subroutine add_local_dos(val,pos)
    real(dp), intent(in) :: pos(3)
    real(dp), intent(inout) :: val
    real(dp) :: modul, n(3), erf_val
    integer  :: n_x, n_y, n_z
    real(DP) :: val1

    val1 = 0.d0

    do n_x=-n_max,n_max
       do n_y=-n_max,n_max
          do n_z=-n_max,n_max
             n(1:3) = n_x * ion_flux_cell(1:3,1) +&
                  & n_y * ion_flux_cell(1:3,2) +&
                  & n_z * ion_flux_cell(1:3,3)
             modul = modulus(pos(1:3)-n(1:3))
             erf_val = erfc_func(sqrt(eta)*modul)
             val1 = val1 + erf_val/modul
          end do
       end do
    end do

    val = val + val1

  end subroutine add_local_dos


  subroutine init_I_prime ()
    use gvecs, only: ngm, gstart, gg

    integer  :: igm, n_x, n_y, n_z
    real(dp) :: n(3), modul, erf_value, omega

    omega = volcel(ion_flux_cell)
    I_prime = 0.d0

    do n_x=-n_max,n_max
       do n_y=-n_max,n_max
          do n_z=-n_max,n_max
             if ((n_x/=0).or.(n_y/=0).or.(n_z/=0)) then
                n(1:3) = n_x * ion_flux_cell(1:3,1) +&
                     & n_y * ion_flux_cell(1:3,2) +&
                     & n_z * ion_flux_cell(1:3,3)
                modul=modulus(n(:))
                erf_value=erfc_func(sqrt(eta)*modul)
                I_prime = I_prime + erf_value / modul
             end if
          end do
       end do
    end do

    I_prime=I_prime-2.d0*sqrt(eta/pi)

    I_prime_rec=0.d0
    do igm=gstart,ngm
       I_prime_rec=I_prime_rec+2.d0*(4.d0 * pi)/omega *&
            & alpha_0_lr(eta,gg(igm),1) ! QE: scaled, gg(igm)*tpiba2
    end do

    I_prime_rec=I_prime_rec-pi/(eta*omega)
    I_prime=I_prime+I_prime_rec

  end subroutine init_I_prime


  subroutine init_I_uno_g ()
    use gvecs, only: ngm, gstart, g, gg

    integer  :: igm, a, b
    real(dp) :: omega

    omega = volcel(ion_flux_cell)
    allocate(I_uno_g(ngm,3,3))

    do a=1,3
       do b=1,3
          if ( a >= b ) then
             do igm=gstart,ngm
                I_uno_g(igm,a,b)=(4.d0*pi)/omega *&
                     & exp(-gg(igm)/(4.d0*eta))/(gg(igm)) *&
                     & (2 + gg(igm)/(2.d0*eta)) * &
                     & g(a,igm) * g(b,igm) / gg(igm)
             end do
             if (gstart==2) I_uno_g(1,a,b) = 0.d0
          end if
       end do
    end do

  end subroutine init_I_uno_g


  subroutine init_I_dos_g ()
    use gvecs, only: ngm, gstart, gg

    integer  :: igm
    real(dp) :: omega

    omega = volcel(ion_flux_cell)
    allocate(I_dos_g(ngm))

    do igm=gstart,ngm
       I_dos_g(igm)=(4.d0 * pi)/omega * alpha_0_lr(eta,gg(igm),1)
    end do

    if (gstart==2) I_dos_g(1) = -pi/(eta*omega)

  end subroutine init_I_dos_g


  subroutine I_dos_value(y,x)
    use gvecs, only: ngm, gstart, g, gg
    ! use wvfct,      only :npw

    real(DP), intent(in)  :: x(3)
    real(DP), intent(out) :: y
    integer  :: igm
    real(DP) :: scalar

    y = 0.d0
    ! do igm=gstart,npw
    do igm=gstart,ngm                          ! <- FIXME, or don't
       if (gg(igm)/(4.d0*eta) > 20.d0) exit    ! NOTE: npw? exit?
       y = y+2.d0*( I_dos_g(igm)*cos(dot_product(g(1:3,igm), x(1:3))) )
    end do

    if (gstart==2) y = y + I_dos_g(1)
    call add_local_dos(y,x)

  end subroutine I_dos_value


  subroutine I_uno_value(y,x,a,b)
    use gvecs, only: ngm, gstart, g, gg
    ! use wvfct, only :npw
    real(dp), intent(out) :: y
    real(dp),intent(in)   :: x(3)
    integer,intent(in)    :: a,b
    integer  :: igm
    real(dp) :: comp_iso

    y = 0.d0
    ! do igm=gstart,npw
    do igm=gstart,ngm                         ! <- FIXME, or don't
       if  (gg(igm)/(4.d0*eta) > 20.d0)	exit  ! NOTE: npw? gl? exit?
       y = y + 2.d0*(I_uno_g(igm,a,b)*cos(dot_product(g(1:3,igm), x(1:3))))
       ! if (gl(igm)>ng_max) exit
    end do

    if (gstart==2) y = y + I_uno_g(1,a,b)
    call add_local_uno(y,x,a,b)

    if (a==b) then
       call I_dos_value(comp_iso,x)
       y = y - comp_iso
    end if

  end subroutine I_uno_value


  subroutine init_ion_flux_data ()
    ion_flux_Jion(:) = 0.0_dp
    ion_flux_a(:) = 0.0_dp
    ion_flux_b(:) = 0.0_dp
    ion_flux_c(:) = 0.0_dp
    ion_flux_d(:) = 0.0_dp
    ion_flux_e(:) = 0.0_dp

    call init_I_prime()
    call init_I_uno_g()
    call init_I_dos_g()
  end subroutine init_ion_flux_data


  subroutine reset_ion_flux_data ()
    deallocate(I_uno_g)
    deallocate(I_dos_g)
  end subroutine reset_ion_flux_data

end module ion_flux_data
