!     This file is part of the SIESTA package.
!     
!     Copyright (c) Fundacion General Universidad Autonoma de Madrid:
!     E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal
!     and J.M.Soler, 1996-2006.
!     
!     Use of this software constitutes agreement with the full conditions
!     given in the SIESTA license, as signed by all legitimate users.
!     
! Auxiliary routines for pao genaration
module atom_pao_util
  use precision
  use pseudopotential, only:pseudopotential_t
  use atom_generation_types, only:basis_def_t,basis_parameters,shell_t,lshell_t,energies_t
  use radial
  use m_recipes, only: polint
  use atom_types, only: get_lmax_orbs 
  implicit none
  
  character(len=1)     :: sym(0:4) = (/ 's','p','d','f','g' /)

  public normalize_orbital,generate_vsoft,calculate_energies, sym, number_of_orbs
  public total_charge
  private

  real(dp) :: exponent, rcsan, rinn, vcte
  
  contains

  !-------------------------------------------------------------

  subroutine normalize_orbital(shell,izeta)
    !    Normalization of basis functions***

    type(shell_t),intent(inout)          :: shell

    integer, intent(in) :: izeta

    integer :: l

    l=shell%l

    call rad_normalize_r_l_1(shell%orb(izeta),l)

  end subroutine normalize_orbital

  !--------------------------------------------------------------------------

  subroutine calculate_energies(shell,iorb)
    !Kinetic, ionic and potential
    type(shell_t),intent(inout)        :: shell
    integer,intent(in)                 :: iorb

    !Internal vars
    
    real(dp)                          :: lambda,ekin,ionic,screened
    integer                           :: l
    type(rad_func_t)                  :: vsum 
 
    l = shell%l
    lambda = shell%lambda(iorb)

    !Init energies
    ionic    = 0.0_dp
    ekin     = 0.0_dp
    screened = 0.0_dp

    ekin = rad_kinetic_energy(shell%orb(iorb),l)
    
    ekin=ekin/lambda**2

    !The screend potential is the sum of ve_PAO and the pseudo
    vsum = rad_sum(shell%ve_PAO,shell%pseudo)
    screened = rad_potential_energy(shell%orb(iorb),vsum,l)

    ionic    = rad_potential_energy(shell%orb(iorb),shell%pseudo,l)
           
    shell%energies(iorb)%ekin     = ekin
    shell%energies(iorb)%screened = screened
    shell%energies(iorb)%ionic    = ionic
    shell%energies(iorb)%total    = ekin+screened
    call rad_dealloc(vsum)
    
  end subroutine calculate_energies

!----------------------------------------------------

  function generate_vsoft(shell,label,write_file) result(vtot)
    !This subroutine adds the soft confinemment potential
    !described in J. Junquera et. al (PRB) to the pseudopotential.
    !The resulting pseudopotential will be integrated to generate
    !the orbitals used by Siesta.
    type(shell_t), intent(inout)         :: shell
    character(len=*), intent(in)         :: label
    logical,intent(in)                   :: write_file
    type(rad_func_t) :: vtot

    integer   :: iu, ircsan !i/o unit
    character(len=80)   :: filename
  
    !rcsan = shell%rc(1)+1.0E-1
    ircsan = rad_get_ir_from_r(shell%ve_pao,shell%rc(1))
    rcsan = rad_get_r_from_ir(shell%ve_pao,ircsan+1)
    rcsan = rcsan + 1.0e-6_dp
    vcte = shell%vcte
    if (shell%rinn < 0.0_dp) then
       rinn = -shell%rinn*rcsan
    else
       rinn = shell%rinn
    endif
    
    vtot=rad_sum_function(shell%ve_pao,soft_confinement,rinn,rcsan)

    if (write_file) then
       write(filename,"(a,i1,a,i1,a,a)") trim(label) // ".L",shell%l, ".",shell%i_sm,".", "confpot"
    
       call io_assign(iu)
       open(unit=iu,file=filename,status='unknown')

       write(iu,'(2a)') '# Soft confinement potential for ', trim(label)
       write(iu,'(a,2i3)') '#    Soft confinement for shell l, nsm = ',shell%l,shell%i_sm
       write(iu,'(a,f10.4,a)') '#        Inner radius    (r_inn) = ', rinn,' bohrs'
       write(iu,'(a,f10.4,a)') '#        External radius (r_ext) = ', rcsan,' bohrs'
       write(iu,'(a,f10.4,a)') '#        Prefactor       (V_0)   = ', vcte,' Ry'

       call rad_dump_ascii(vtot,iu)
       call io_close(iu)
    endif

  end function generate_vsoft

  !------------------------------------------------------

  function soft_confinement(r) result(v)
    use precision
    real(dp), intent(in) :: r
    real(dp) :: v

    exponent = -(rcsan-rinn) / (r - rinn)
    
    if( abs(exponent) .gt. 500.d0) then
       v = 0.0_dp
    else
       if(r < rinn) then
          v=0.0_dp
       elseif ( r <= rcsan) then
          v = vcte/(rcSan-r)*exp(exponent)
       else
          v = 0.0_dp
       endif
    endif

  end function soft_confinement

  !---------------------------------------------------------------

  function number_of_orbs(isp) result(norbs)
    integer, intent(in) :: isp
    integer             :: norbs

    integer             :: nsm,l

    type(lshell_t),    pointer :: lshell
    type(shell_t),     pointer :: shell
    type(basis_def_t), pointer :: basp     !Parameters corresponding to this basis set.

    basp => basis_parameters(isp)

    !Find total number of orbs, including degeneracy
    norbs=0
    do l=0, get_lmax_orbs(isp)
       lshell => basp%lshell(l)
       do nsm=1,lshell%nn
          shell => lshell%shell(nsm)
          norbs = norbs+shell%nzeta 
       end do
    enddo
  end function number_of_orbs

   !---------------------------------------------------

   function total_charge(is,chval,rcocc) result(rho)
    integer, intent(in) :: is
    real(dp), intent(out) :: chval, rcocc
    
    
    type(rad_func_t)                :: rho, rphi2, rho_old, rho_shell
    type(shell_t),     pointer      :: shell    !Pointer to a shell.
    type(lshell_t),    pointer      :: lshell   !Pointer to a l-shell.
    real(dp), dimension(:), pointer :: values
    type(rad_grid_t)                :: grid
    type(pseudopotential_t),pointer :: vps
    integer :: l,nsm,izeta
    real(dp) :: ch_temp,pop
    type(basis_def_t),  pointer     :: basp
    
    basp=>basis_parameters(is)
    vps => basis_parameters(is)%pseudopotential

    grid = rad_get_grid(vps%vlocal)
    allocate(values(1:rad_grid_get_length(grid)))
    values=0.0_dp
    call rad_alloc(rho,values,grid)

    chval=0.0_dp
    rcocc=0.0_dp
    do l= 0,get_lmax_orbs(is)
       lshell => basp%lshell(l)

       do nsm=1,lshell%nn

          shell => lshell%shell(nsm)

          do izeta=1,shell%nzeta
             
             pop = shell%population(izeta)*dble((2*l+1))
             
             if(pop .gt. 0.0_dp) then

                call rad_copy(rho,rho_old)
                call rad_dealloc(rho)
                
                rphi2 = rad_multiply(shell%rphi(izeta),shell%rphi(izeta))
                rho_shell = rad_multiply_each_value(rphi2,pop)
                rcocc = max(rcocc,rad_cutoff(rphi2))
                
                rho = rad_sum(rho_old,rho_shell)
                call rad_dealloc(rho_old)
                ch_temp = rad_integral(rho_shell) 
                chval = chval + ch_temp
                
                call rad_dealloc(rphi2)
                call rad_dealloc(rho_shell)

             end if
          enddo
       enddo
    enddo

    call rad_set_origin(rho,0.0_dp)
    deallocate(values)
    call rad_grid_dealloc(grid)
  end function total_charge

end module atom_pao_util
