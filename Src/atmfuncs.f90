! 
! This file is part of the SIESTA package.
!
! Copyright (c) Fundacion General Universidad Autonoma de Madrid:
! E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal
! and J.M.Soler, 1996-2006.
! 
! Use of this software constitutes agreement with the full conditions
! given in the SIESTA license, as signed by all legitimate users.
!
module atmfuncs

  !     This file contains a set of routines which provide all the information
  !     about the basis set, pseudopotential, atomic mass, etc... of all the
  !     chemical species present in the calculation.

  !     The routines contained in this file can only be called after they
  !     are initialized by calling the subroutine 'atom' for all the 
  !     different chemical species in the calculation:

  use precision
  use sys, only: die
  use atom_types
  use atmfuncs_types
  use spher_harm, only: rlylm
  implicit none 

  private

  !
  type(species_info_t), pointer        :: spp
  !
  
  integer, parameter               :: max_l = 5
  integer, parameter               :: max_ilm = (max_l+1)*(max_l+1)

  
  private :: chk, max_l, max_ilm, message

  public  :: nofis, nkbfis, izofis, massfis
  public  :: rcore, rcut, chcore_sub, epskb, uion
  public  :: atmpopfio, psch, zvalfis, floating, psover
  public  :: lofio, symfio, cnfigfio, zetafio, mofio
  public  :: labelfis, lomaxfis, nztfl, rphiatm, lmxkbfis
  public  :: phiatm, all_phi, xphiatm, yphiatm, zphiatm
  public  :: check_atmfuncs
  public  :: func_t,orb_f,kbpj_f,vlocal_f,vna_f,chlocal_f,core_f,get_func_type,func_t_compare
 
  !
  
contains

  subroutine all_phi(is,func,r,nphi,phi, grphi)
    integer, intent(in) :: is     ! Species index
     type(func_t), intent(in) :: func  ! Function label
    !integer, intent(in) :: it     ! Orbital-type switch: basis, kbs, ldau
    real(dp), intent(in)  :: r(3)   ! Point vector, relative to atom
    integer, intent(out):: nphi   ! Number of phi's
    real(dp), intent(out) :: phi(:) ! Basis orbital, KB projector, or
    !  local pseudopotential
    real(dp), optional, intent(out) :: grphi(:,:) ! Gradient of phi

    !  Returns Kleynman-Bylander local pseudopotential, nonlocal projectors,
    !  and atomic basis orbitals (and their gradients).
    !  Same as phiatm but returns all orbitals or KB projectors of the atom
    !  Written by D.Sanchez-Portal and J.M.Soler. Jan. 2000 

    ! Distances in Bohr
    ! 1) Each projector and basis function has a well defined total
    !    angular momentum (quantum number l).
    ! 2) Basis functions are normalized and mutually orthogonal
    ! 3) Projection functions are normalized and mutually orthogonal
    ! 4) Normalization of KB projectors |Phi_lm> is such that 
    !     <Psi|V_KB|Psi'> = <Psi|V_local|Psi'> +
    !                   Sum_lm( epsKB_l * <Psi|Phi_lm> * <Phi_lm|Psi'> )
    !    where epsKB_l is returned by function EPSKB
    ! 5) Prints a message and stops when no data exits for IS
    ! 6) Returns exactly zero when |R| > RCUT(IS,IO)
    ! 8) If arrays phi or grphi are too small, returns with the required
    !    value of nphi

    integer i, jlm, l, lmax, m, maxlm, n
    double precision  rmod, phir
    real(dp) rly(max_ilm), grly(3,max_ilm)

    integer, parameter :: maxphi=100

    integer :: ilm(maxphi)
    double precision :: rmax(maxphi)
    logical :: within(maxphi)

    call chk('all_phi',is)
    spp => species(is)

    !     Find number of orbitals
    if (func==orb_f) then
       nphi=get_number_of_orbs(spp)
    elseif (func==kbpj_f) then
       nphi=get_number_of_kb_projs(spp)
    else
       call die("all_phi: Please use phiatm to get Vna...")
    endif

    if (nphi.gt.maxphi) call die('all_phi: maxphi too small')

    if (func==orb_f) then
       do i = 1, nphi
          l = get_orb_l(spp,i)
          m = get_orb_m(spp,i)
          ilm(i) = l*(l+1)+m+1
          rmax(i) = get_orb_cutoff(spp,i)
       enddo
    elseif(func==kbpj_f) then
       do i = 1, nphi
          rmax(i) =get_kb_proj_cutoff(spp,i)
          l = get_kb_proj_l(spp,i)
          m = get_kb_proj_m(spp,i)
          ilm(i) = l*(l+1)+m+1
       enddo
    else
       call die("all_phi: unknown function type")
    endif

    !     Check size of output arrays
    if (present(grphi)) then
       if (size(grphi,1).ne.3) call die('all_phi: incorrect first dimension of grphi')
       n = min( size(phi), size(grphi,2) )
    else
       n = size(phi)
    endif
    !     Return if the caller did not provide arrays large enough...
    if (n.lt.nphi) return

    !     Initialize orbital values
    phi(1:nphi) = 0._dp
    if (present(grphi)) grphi(:,1:nphi) = 0._dp

    if ((func == vna_f) .and. floating(is)) return
    !if ((it.lt.0) .and. floating(is)) return

    !     Find for which orbitals rmod < rmax and test for quick return
    rmod = sqrt(sum(r*r)) + tiny20
    within(1:nphi) = ( rmax(1:nphi) > rmod )
    if (.not.any(within(1:nphi))) return

    !     Find spherical harmonics
    maxlm = maxval( ilm(1:nphi), mask=within(1:nphi) )
    lmax=nint(sqrt(real(maxlm,dp)))-1
    call rlylm(lmax,r,rly,grly)

    !     Find values

    i_loop: do i=1,nphi

       !       Check if rmod > rmax
       if (.not.within(i)) cycle i_loop

       !       Find radial part

       if (func==orb_f) then
          call get_value_of_orb(spp,i,r,phir,grphi(:,i))          
       elseif(func==kbpj_f) then 
          call get_value_of_kb_proj(spp,i,r,phir,grphi(:,i))
       else
          call die("all_phi: unknown function type")
       endif

       !       Multiply radial and angular parts
       jlm = ilm(i)
       phi(i) = phir * rly(jlm)
       if (present(grphi))then
          grphi(:,i) = grphi(:,i) * rly(jlm) + phir * grly(:,jlm)
       endif
    enddo i_loop

  end subroutine all_phi

  !------------------------------------------------------------------------
  
   FUNCTION atmpopfio(IS,IO)
    real(dp) atmpopfio
    integer, intent(in) :: is    ! Species index
    integer, intent(in) :: io    ! Orbital index (within atom)

    ! Returns the population of the atomic basis orbitals in the atomic 
    ! ground state configuration.

    call chk('atmpopfio',is)
    if ( (io .gt. get_number_of_orbs(species(is))) .or. (io .lt. 1))   call die("atmpopfio: Wrong io")

    atmpopfio = get_orb_pop(species(is),io)
  end function atmpopfio

  !-------------------------------------------------------------------------

  subroutine chcore_sub(is,r,ch,grch)
    integer, intent(in) :: is      ! Species index
    real(dp), intent(in)  :: r(3)    ! Point vector, relative to atom
    real(dp), intent(out) :: ch      ! Value of pseudo-core charge dens.
    real(dp), intent(out) :: grch(3) ! Gradient of pseudo-core ch. dens.

    ! Returns returns pseudo-core charge density for non-linear core correction
    ! in the xc potential.
    ! Distances in Bohr, Energies in Rydbergs, Density in electrons/Bohr**3
    !  2) Returns exactly zero when |R| > Rcore

    call chk('chcore_sub',is)

    ch=0.0_dp
    grch(1:3)=0.0_dp

    if (floating(is)) return
    
    call get_value_of_core_charge(species(is), r,ch,grch)

  end subroutine chcore_sub
  
  !-------------------------------------------------------------------------

  subroutine chk(name,is)
    character(len=*), intent(in) :: name

    integer, intent(in) :: is

    character(len=79) message

    if ((is.lt.1).or.(is.gt.nspecies)) then 
       write(message,'(2a,i3,a,i3)') name, ": Wrong species", is, ". Have", nspecies
       call die(message)
    endif
  end subroutine chk

  !-------------------------------------------------------------------------
  
  subroutine check_atmfuncs
    !Sub which dumps all the stored info of all the species.
    !This is all the info which will be used by siesta

    real(dp):: r(3), val, grphi, gr(3),delta,rm,rc
    integer ::is1, is2, is,npoints=50,i,ir
    !type(rad_func_t), pointer :: func
    type(species_info_t), pointer :: spp
    character(len=40) :: filename
    real(dp) :: energy,dedr

    print ('(/,a)'),"check_atmfuncs begin"

    !------------elec_corr---------------------
    open(unit=88,file="elec_corr_check.dat",status="replace")
    do is1=1,nspecies
       do is2=1,nspecies
          if (.not. floating(is2) .and. .not. floating(is1))then
             rc=2.0_dp*max(get_vna_cutoff(species(is1)), &
                  get_vna_cutoff(species(is2)))
             delta = rc/dble(npoints)
             rm=0.0_dp
             write(88,*) "elec_corr: is1,is2",is1,is2
             do ir=1,npoints
                call psover(is1,is2,rm,energy,dedr)
                write(88,*) rm,energy,dedr
                rm=rm+delta
             enddo
          endif
       enddo
    enddo
    close(88)
    

    !---------------------------------------

    do is=1,nspecies
       print *, "Specie=",is,"label=",trim(get_label(spp))
       spp => species(is)
       write(filename,'(a,a)') trim(get_label(spp)),"-atomcheck.dat"
       open(unit=88,file=filename,status="replace")

       write(88,*) "Label=",labelfis(is)
       write(88,*) "Atomic number=",izofis(is)
       write(88,*) "Atomic mass=",massfis(is)
       write(88,*) "Valence charge=",zvalfis(is)
       write(88,*) "Self energy=",uion(is)
       write(88,*) "is floating?",floating(is)
       do i=1,4
          write(88,*) " l,n_orbs=",i,nztfl(is,i)
       enddo

       do i=1,4
          write(88,*) " l,n_kbs=",i,nkbl_func(is,i)
       enddo

       if(rcore(is) > 0.0_dp)then
          write(88,*) 'Chcore---------------'
          r=0.0_dp
          rc=rcore(is)
          delta=rc/dble(npoints)
          do i=1,npoints
             call chcore_sub(is,r,val,gr)
             write(88,'(7f9.3)') r,val,gr
             r=r+delta
          enddo
       endif

       !---------------------------------------
       write(88,*) 'Vna -----------------'
       if (has_neutral_atom_potential(spp)) then
          rc = rcut(is,vna_f,0)
          delta=rc/dble(npoints)
          write(88,*) 'Species_info ',is, ' rcut Vna: ', rc
          rm=0.0_dp
          do i=1,npoints
             call rphiatm(is,vna_f,0,rm,val,grphi)
             write(88,*) rm, val
             rm=rm+delta
          enddo
       endif
       
       !----------------------------------------
       write(88,*) 'Chloc-----------------'
       rc=5.0_dp
       delta=rc/dble(npoints)
       r=0.0_dp
       do i=1,npoints
          call psch(is,r,val,gr)
          write(88,'(7f9.3)') r,val,gr
          r=r+delta
       enddo
      

       !----------------------------------------
       !kbs
       write(88,*) 'kbs----------------------'
       do i=1,nkbfis(is)
          write(88,*) '*****************'
          write(88,*) 'l=',lofio(is,kbpj_f,i)
          write(88,*) 'm=',mofio(is,kbpj_f,i)
          write(88,*) 'rcut=',rcut(is,kbpj_f,i)
          write(88,*) 'Energy',epskb(is,i)
          write(88,*) 'Sym=',symfio(is,kbpj_f,i)

          r=0.0_dp
          rc = rcut(is,kbpj_f,i)
          delta=rc/(dble(npoints))
          do ir=1,npoints
             call phiatm(is,kbpj_f,i,r,val,gr)
             write(88,'(7f9.3)') r,val,gr
             r=r+0.1*delta
          enddo
       enddo

       !----------------------------------------
       !orb
       write(88,*) 'orbs----------------------'
       write(88,*) "norbs=",is,nofis(is)
       do i=1,nofis(is)
          write(88,*) '*****************'
          write(88,*) 'l=',lofio(is,orb_f,i)
          write(88,*) 'm=',mofio(is,orb_f,i)
          write(88,*) 'rcut=',rcut(is,orb_f,i)
          write(88,*) 'Sym=',symfio(is,orb_f,i)
          write(88,*) 'polarization?',pol(is,i)
          write(88,*) 'pop=',atmpopfio(is,i)
          r=0.0_dp
          rc = rcut(is,orb_f,i)
          delta=rc/(dble(npoints))
          do ir=1,npoints
             call phiatm(is,orb_f,i,r,val,gr)
             write(88,'(7f9.3)') r,val,gr
             r=r+0.1_dp*delta
          enddo
       enddo
       
    enddo

    print ('(/,a)'), "check_atmfuncs end"
  end subroutine check_atmfuncs

  !----------------------------------------------------------------------

  FUNCTION cnfigfio(IS,IO)
    integer cnfigfio
    integer, intent(in) :: is    ! Species index
    integer, intent(in) :: io    ! Orbital index (within atom)

    ! Returns the valence-shell configuration in the atomic ground state
    ! (i.e. the principal quatum number for orbitals of angular momentum l)

    !   INTEGER CNFIGFIO: Principal quantum number of the shell to what 
    !                     the orbital belongs ( for polarization orbitals
    !                     the quantum number corresponds to the shell which
    !                     is polarized by the orbital io) 


    call chk('cnfigfio',is)
    if (io .gt. get_number_of_orbs(species(is)) .or. (io .lt. 1))  &
         call die("cnfigfio: Wrong io")

    cnfigfio = get_orb_n(species(is),io)

  end function cnfigfio
  
  !------------------------------------------------------------------------

   FUNCTION epskb (IS,IO)
    real(dp) epskb
    integer, intent(in)   ::  is   ! Species index
    integer, intent(in)   ::  io   ! KB proyector index (within atom)
    ! May be positive or negative 
    ! (only ABS(IO) is used).

    !  Returns the energies epsKB_l of the Kleynman-Bylander projectors:
    !       <Psi|V_KB|Psi'> = <Psi|V_local|Psi'> +
    !                 Sum_lm( epsKB_l * <Psi|Phi_lm> * <Phi_lm|Psi'> )
    !  where Phi_lm is returned by subroutine PHIATM.
    !  Energy in Rydbergs.

    integer :: ik

    spp => species(is)
    ik = abs(io)
    if ((ik.gt. get_number_of_kb_projs(spp)) .or. (ik .lt. 1) )  call die("epskb: No such projector")
    epskb = get_kb_proj_energy(spp, ik) 
  end function epskb

  !--------------------------------------------------------------------

  function floating(is)
    !       Returns .true. if the species is really a "fake" one, intended
    !       to provide some floating orbitals.
    logical :: floating
    integer, intent(in) :: is

    floating = izofis(is) .lt. 0
  end function floating

  !--------------------------------------------------------------------

  FUNCTION izofis( IS )
    integer :: izofis ! Atomic number
    integer, intent(in) :: is ! Species index

    call chk('izofis',is)

    izofis = get_atomic_number(species(is))

  end function izofis

  !--------------------------------------------------------------------

  FUNCTION labelfis (IS)
    character(len=20) ::  labelfis  ! Atomic label
    integer, intent(in) :: is            ! Species index

    call chk('labelfis',is)
    labelfis= get_label(species(is))
  end function labelfis

  !-------------------------------------------------------------------

  FUNCTION lmxkbfis (IS)
    integer :: lmxkbfis    ! Maximum ang mom of the KB projectors
    integer, intent(in) :: is            ! Species index

    call chk('lmxkbfis',is)
    lmxkbfis= get_lmax_kb_proj(species(is))
  end function lmxkbfis
  
  !--------------------------------------------------------------------

  FUNCTION lofio (IS,func,io)
    integer lofio
    integer, intent(in)          :: is    ! Species index
    type(func_t), intent(in)     :: func  !Func label
    integer, intent(in)          :: io    ! Function index (within specie)

    !     Returns total angular momentum quantum number of a given atomic 
    !     basis orbital or Kleynman-Bylander projector.
    !
    !     INTEGER  IS   : Species index       
    !     INTEGER  IO   : Orbital index (within specie)
    !     CHARACTER(len=func_label_length) function_label  : Function label orb_f,kbpj_f,etc             
    !                   
    !     
    !***********************OUTPUT*****************************************
    !   INTEGER LOFIO  : Quantum number L of orbital or KB projector
    call chk('lofio',is)

    spp => species(is)

    if (func==orb_f) then
       if (io.gt. get_number_of_orbs(spp)) call die("lofio: No such orbital")
       lofio = get_orb_l(spp,io)
    else if (func==kbpj_f) then
       if (io.gt. get_number_of_kb_projs(spp)) call die("lofio: No such projector")
       lofio = get_kb_proj_l(spp,io)
    else
       lofio = 0
    endif

  end function lofio

  !------------------------------------------------------------------------

  FUNCTION lomaxfis (IS)
    integer :: lomaxfis  ! Maximum ang mom of the Basis Functions
    integer, intent(in) :: is            ! Species index

    call chk('lomaxfis',is)

    lomaxfis = get_lmax_orbs(species(is)) 
  end function lomaxfis

  !--------------------------------------------------------------------
  
  FUNCTION massfis(IS)
    real(dp) :: massfis            ! Mass
    integer, intent(in) :: is            ! Species index

    call chk('massfis',is)
    massfis=get_mass(species(is))
  end function massfis

  !--------------------------------------------------------------------

  FUNCTION mofio (IS,func,IO)
    integer mofio
    integer,    intent(in) :: is    ! Species index
    type(func_t), intent(in) :: func!Func label
    integer,    intent(in) :: io    ! Orbital index (within atom)

    !   Returns m quantum number of a given atomic basis
    !   basis orbital or Kleynman-Bylander projector.

    !    INTEGER  IO   : Orbital index (within atom)
    !
    !***********************OUTPUT*****************************************
    !   INTEGER MOFIO  : Quantum number m of orbital or KB projector

    call chk('mofio',is)

    spp => species(is)

    if (func == orb_f) then
       if (io.gt. get_number_of_orbs(spp))  call die("Mofio: No such orbital")
       mofio = get_orb_m(spp,io)
    else if (func == kbpj_f) then
       if (io.gt. get_number_of_kb_projs(spp))  call die("Mofio: No such projector")
       mofio = get_kb_proj_m(spp,io)
    else
       mofio = 0
    endif

  end function mofio

  !------------------------------------------------------------------------

  FUNCTION nkbfis(IS)
    integer :: nkbfis    ! Total number of KB projectors
    integer, intent(in) :: is            ! Species index

    call chk('nkbfis',is)
    nkbfis = get_number_of_kb_projs(species(is))
  end function nkbfis
  
  !--------------------------------------------------------------------

   FUNCTION nkbl_func (IS,L)
    integer nkbl_func
    integer, intent(in)  :: is   ! Species index
    integer, intent(in)  :: l    ! Angular momentum of the basis funcs

    ! Returns the number of different KB projectors
    ! with the same angular momentum and for a given species

    integer i

    call chk('nkbl_func',is)

    spp => species(is)

    nkbl_func = 0
    do i = 1, get_number_of_kb_projs(spp)
       if (get_kb_proj_l(spp,i).eq.l) nkbl_func = nkbl_func+1
    enddo

  end function nkbl_func

  !--------------------------------------------------------------------

  FUNCTION nofis(IS)
    integer :: nofis    ! Total number of Basis functions
    integer, intent(in) :: is            ! Species index

    call chk('nofis',is)
    nofis = get_number_of_orbs(species(is))
  end function nofis

  !--------------------------------------------------------------------

  !
  !     Deprecated
  !
  FUNCTION nztfl (IS,L)
    integer nztfl
    integer, intent(in)  :: is   ! Species index
    integer, intent(in)  :: l    ! Angular momentum of the basis funcs
    ! Returns the number of different basis functions
    ! with the same angular momentum and for a given species

    integer i

    call chk('nztfl',is)
    spp => species(is)

    nztfl = 0

    do i = 1, get_number_of_orbs(spp)
       if (get_orb_l(spp,i).eq.l) nztfl = nztfl+1
    enddo

  end function nztfl

  !-----------------------------------------------------------------

  subroutine phiatm(is,func,io,r,phi,grphi)
    integer, intent(in) :: is        ! Species index
    type(func_t), intent(in) :: func !Func label
    integer, intent(in) :: io        ! Orbital index (within atom)
    real(dp), intent(in)  :: r(3)    ! Point vector, relative to atom
    real(dp), intent(out) :: phi     ! Basis orbital, KB projector, or
    !  local pseudopotential
    real(dp), intent(out) :: grphi(3)! Gradient of BO, KB proj, or Loc ps

    !  Returns Kleynman-Bylander local pseudopotential, nonlocal projectors,
    !  and atomic basis orbitals (and their gradients).
    ! Distances in Bohr
    ! 1) Each projector and basis function has a well defined total
    !    angular momentum (quantum number l).
    ! 2) Basis functions are normalized and mutually orthogonal
    ! 3) Projection functions are normalized and mutually orthogonal
    ! 4) Normalization of KB projectors |Phi_lm> is such that 
    !     <Psi|V_KB|Psi'> = <Psi|V_local|Psi'> +
    !                   Sum_lm( epsKB_l * <Psi|Phi_lm> * <Phi_lm|Psi'> )
    !    where epsKB_l is returned by function EPSKB
    ! 5) Prints a message and stops when no data exits for IS and/or IO
    ! 6) Returns exactly zero when |R| > RCUT(IS,IO)
    ! 7) PHIATM with IO = 0 is strictly equivalent to VNA_SUB

    real(dp) rmod, phir,cutoff
    real(dp) rly(max_ilm), grly(3,max_ilm)
    integer i, l, m, ilm

    phi=0.0_dp
    grphi(1:3)=0.0_dp

    spp => species(is)
    if (func==orb_f) then
       if (io.gt. get_number_of_orbs(spp))  call die("phiatm: No such orbital")
       call get_value_of_orb(spp,io,r,phir,grphi)
       cutoff = get_orb_cutoff(spp,io)
       l = get_orb_l(spp,io)
       m = get_orb_m(spp,io)
    else if (func==kbpj_f) then
       if (floating(is)) return
       if (io.gt. get_number_of_kb_projs(spp))  call die("phiatm: No such projector")
       call get_value_of_kb_proj(spp,io,r,phir,grphi)
       cutoff = get_kb_proj_cutoff(spp,io)
       l = get_kb_proj_l(spp,io)
       m = get_kb_proj_m(spp,io)
    else if (func==vna_f)then
       if (floating(is)) return
       call get_value_vna(spp,r,phir,grphi)
       cutoff = get_vna_cutoff(spp)
       l = 0
       m = 0
    else
       call die("phiatm: unknown function")
    endif

    rmod = sqrt(sum(r*r)) + tiny20
    if(rmod .gt. cutoff-tiny12) return

    if (io.eq.0 .or. func == vna_f) then
       phi=phir
    else
       ilm = l*l + l + m + 1
       call rlylm( l, r, rly, grly )
       phi = phir * rly(ilm)
       do i = 1,3
          grphi(i)=grphi(i)*rly(ilm)+phir*grly(i,ilm)
       enddo

    endif

  end subroutine phiatm

  !---------------------------------------------------------------------

   FUNCTION pol (IS,IO)
    logical pol
    integer, intent(in) :: is    ! Species index
    integer, intent(in) :: io    ! Orbital index (within atom)
    ! io>0 => basis orbitals

    ! If true, the orbital IO is a perturbative polarization orbital

    spp => species(is)

    if ( (io .gt. get_number_of_orbs(spp)) .or. (io .le. 0))   call die("pol: Wrong io")
    spp => species(is)
    pol = get_orb_pol(spp,io)

  end function pol

  !---------------------------------------------------------------------

   subroutine psch(is,r,ch,grch)
    integer, intent(in) :: is      ! Species index
    real(dp), intent(in)  :: r(3)    ! Point vector, relative to atom
    real(dp), intent(out) :: ch      ! Local pseudopot. charge dens.
    real(dp), intent(out) :: grch(3) ! Gradient of local ps. ch. dens.

    ! Returns 'local-pseudotential charge density'.
    ! Distances in Bohr, Energies in Rydbergs
    ! Density in electrons/Bohr**3
    !  2) Returns exactly zero when |R| > Rchloc

    call chk('psch',is)

    ch=0.0_dp 
    grch(1:3)=0.0_dp 

    if (floating(is)) return

    call get_value_pseudo_local_charge(species(is), r,ch,grch)
  end subroutine psch

  !----------------------------------------------------------------------

   !
  !
  !     This routine takes two species arguments
  !
  subroutine psover(is1,is2,r,energ,dedr)
    integer, intent(in) :: is1, is2     ! Species indexes
    real(dp), intent(in)  :: r       ! Distance between atoms
    real(dp), intent(out) :: energ   ! Value of the correction
    !  interaction energy
    real(dp), intent(out) :: dedr    ! Radial derivative of the correction

    ! Returns electrostatic correction to the ions interaction energy
    ! due to the overlap of the two 'local pseudopotential charge densities'
    ! Distances in Bohr, Energies in Rydbergs
    !  2) Returns exactly zero when |R| > Rchloc

    integer ismx, ismn, indx
    

    call chk('psover',is1)
    call chk('psover',is2)

    energ=0.0_dp 
    dedr=0.0_dp 

    if (floating(is1) .or. floating(is2)) return

    ismx=max(is1,is2)
    ismn=min(is1,is2)
    indx=((ismx-1)*ismx)/2 + ismn
  
    call get_elec_corr(indx,r,energ,dedr)
    

  end subroutine psover

  !-----------------------------------------------------------------------

   FUNCTION rcore(is)
    real(dp) rcore
    integer, intent(in) :: is    ! Species index

    !  Returns cutoff radius of the pseudo-core charge density for the non-linear
    !   core corrections for xc potential.
    !  Distances in Bohr

    call chk('rcore',is)
    if(has_core_charge(species(is)))then
       rcore = get_rcore(species(is))
    else
       rcore = 0.0_dp
    endif
  end function rcore

  !-----------------------------------------------------------------------

  function rcut(is,func,io)
    real(dp) rcut
    type(func_t), intent(in) :: func  ! Function label
    integer, intent(in)      :: is    ! Species index
    integer, intent(in)      :: io    ! Function index (within atom)

    !  Returns cutoff radius of Kleynman-Bylander projectors and
    !  atomic basis orbitals.
    !  Distances in Bohr

    call chk('rcut',is)
    rcut = 0.0_dp
    spp => species(is)
    if (func==orb_f) then
       if (io.gt. get_number_of_orbs(spp))  call die("No such orbital")
       rcut = get_orb_cutoff(spp,io)
    else if (func==kbpj_f) then
       if (io.gt. get_number_of_kb_projs(spp))  call die("No such projector")
       rcut = get_kb_proj_cutoff(spp,io)
    else if (func==vna_f) then
       rcut = get_vna_cutoff(spp)
    else       
       call die("rcut: unknown function")
    endif

  end function rcut
  !
  !------------------------------------------------------------------------
  !

  subroutine rphiatm(is,func,io,r,phi,dphidr)
    integer, intent(in) :: is        ! Species index
    integer, intent(in) :: io        ! Orbital index (within atom)
    type(func_t), intent(in) :: func !Func label
    real(dp), intent(in)  :: r       ! Radial distance, relative to atom
    real(dp), intent(out) :: phi     ! Basis orbital, KB projector, or
    !  local pseudopotential
    real(dp), intent(out) :: dphidr  ! Radial derivative of BO, 
    !  KB proj, or Loc pseudopot.

    !  Returns the radial component of 
    !  Kleynman-Bylander local pseudopotential, nonlocal projectors,
    !  and atomic basis orbitals (and their radial drivatives)
    ! Distances in Bohr
    ! 1) Each projector and basis function has a well defined total
    !    angular momentum (quantum number l).
    ! 2) Basis functions are normalized and mutually orthogonal
    ! 3) Projection functions are normalized and mutually orthogonal
    ! 4) Normalization of KB projectors |Phi_lm> is such that 
    !     <Psi|V_KB|Psi'> = <Psi|V_local|Psi'> +
    !                   Sum_lm( epsKB_l * <Psi|Phi_lm> * <Phi_lm|Psi'> )
    !    where epsKB_l is returned by function EPSKB
    ! 6) Returns exactly zero when |R| > RCUT(IS,IO)
    ! 7) RPHIATM with ITYPE = 0 is strictly equivalent to VNA_SUB

    real(dp) :: phir,cutoff
    integer  :: l, m

    phi = 0.0_dp
    dphidr = 0._dp

    spp => species(is)
    if (func==orb_f) then
       if (io.gt. get_number_of_orbs(spp))  call die("rphiatm: No such orbital")
       call get_rvalue_of_orb(spp,io,r,phir,dphidr)
       cutoff = get_orb_cutoff(spp,io)
       l = get_orb_l(spp,io)
       m = get_orb_m(spp,io)
    else if (func==kbpj_f) then
       if (floating(is)) return
       if (io.gt. get_number_of_kb_projs(spp))  call die("rphiatm: No such projector")
       call get_rvalue_of_kb_proj(spp,io,r,phir,dphidr)
       l = get_kb_proj_l(spp,io)
       m = get_kb_proj_m(spp,io)
    else if(func==vna_f) then
       if (floating(is)) return
       call get_rvalue_na(spp,r,phir,dphidr)
       l = 0
       m = 0
    else 
       call die("rphiatm:unknown function")
    endif

    if (l.eq.0) then
       phi=phir
    elseif (l.eq.1) then
       phi=phir*r
       dphidr=dphidr*r
       dphidr=dphidr+phir 
    else
       phi=phir*r**l 
       dphidr=dphidr * r**l
       dphidr=dphidr + l * phir * r**(l-1)
    endif

  end subroutine rphiatm

  !------------------------------------------------------------------

   FUNCTION symfio (IS,func,IO)
    character(len=20) symfio
    integer, intent(in) :: is    ! Species index
    type(func_t), intent(in) :: func  ! Function label
    integer, intent(in) :: io    ! Orbital index (within atom)

    ! Returns a label describing the symmetry of the
    !   basis orbital or Kleynman-Bylander projector.
    !    INTEGER  IO   : Orbital index (within atom)
    !                    IO > 0 => Basis orbitals
    !                    IO < 0 => Kleynman-Bylander projectors

    !   INTEGER SYMFIO  : Symmetry of the orbital or KB projector
    !  2) Returns 's' for IO = 0


    integer ilm, i, lorb, morb
    integer, parameter  :: lmax_sym=3

    character(len=6)  sym_label((lmax_sym+1)*(lmax_sym+1)) 
    character(len=7)  paste

    external paste
    !
    data  sym_label(1)          / 's' /
    data (sym_label(i),i=2,4)   / 'py', 'pz', 'px' /
    data (sym_label(i),i=5,9)   / 'dxy', 'dyz', 'dz2', 'dxz', 'dx2-y2' / 
    data (sym_label(i),i=10,16) / 'f', 'f', 'f', 'f', 'f', 'f', 'f' /

    call chk('rcut',is)

    spp => species(is)
    if (func==orb_f) then
       if (io.gt. get_number_of_orbs(spp))  call die("No such orbital")
    else if (func==kbpj_f) then
       if (io.gt. get_number_of_kb_projs(spp))  call die("No such projector")
    else
       symfio = 's'
    endif

    lorb=lofio(is,func,io)
    morb=mofio(is,func,io)

    if(lorb.gt.lmax_sym ) then 
       symfio=' '
    else
       ilm=lorb*lorb+lorb+morb+1  
       if(func==orb_f) then
          if (pol(is,io)) then
             symfio=paste('P',sym_label(ilm))
          else
             symfio=sym_label(ilm)
          endif
       else
          symfio=sym_label(ilm) 
       endif
    endif

  end function symfio

   !-------------------------------------------------------------------------

  FUNCTION uion ( IS )
    real(dp) uion
    integer, intent(in) :: is    ! Species index
    call chk('uion',is)
    uion = get_self_energy(species(is))
  end function uion

  !-------------------------------------------------------------------------

  subroutine vna_sub(is,r,v,grv)
    integer, intent(in) :: is      ! Species index
    real(dp), intent(in)  :: r(3)    ! Point vector, relative to atom
    real(dp), intent(out) :: v       ! Value of local pseudopotential
    real(dp), intent(out) :: grv(3)  ! Gradient of local pseudopotential

    ! Returns local part of neutral-atom Kleynman-Bylander pseudopotential.
    ! Distances in Bohr,  Energies in Rydbergs
    !  2) Returns exactly zero when |R| > RCUT(IS,0)

    call chk('vna_sub',is)

    v=0.0_dp
    grv(1:3)=0.0_dp

    if (floating(is) .or. .not. has_neutral_atom_potential(species(is))) return

    call get_value_vna(species(is),r,v,grv)    
  end subroutine vna_sub

  !------------------------------------------------------------------------
  
  SUBROUTINE xphiatm(is,func,io,r,xphi,grxphi)
    !     Calculates x*phiatm and its gradient

    integer, intent(in)   :: is, io
    type(func_t), intent(in) :: func  ! Function label
    real(dp), intent(in)  :: r(3)
    real(dp), intent(out) :: xphi, grxphi(3)

    real(dp) phi, grphi(3), x

    call phiatm(is,func,io,r,phi,grphi)
    x = r(1)
    xphi = x * phi
    grxphi(1) = x * grphi(1) + phi
    grxphi(2) = x * grphi(2)
    grxphi(3) = x * grphi(3)
  END SUBROUTINE xphiatm

  !-------------------------------------------------------------------

  SUBROUTINE yphiatm(is,func,io,r,yphi,gryphi)
    !     Calculates y*phiatm and its gradient

    integer, intent(in)   :: is, io
    type(func_t), intent(in) :: func  ! Function label
    real(dp), intent(in)  :: r(3)
    real(dp), intent(out) :: yphi, gryphi(3)

    real(dp) phi, grphi(3), y

    call phiatm(is,func,io,r,phi,grphi)
    y = r(2)
    yphi = y * phi
    gryphi(1) = y * grphi(1)
    gryphi(2) = y * grphi(2) + phi
    gryphi(3) = y * grphi(3)
  END SUBROUTINE yphiatm
  
  !-------------------------------------------------------------------

  FUNCTION zetafio (IS,IO)
    integer zetafio
    integer, intent(in) :: is    ! Species index
    integer, intent(in) :: io    ! Orbital index (within atom)

    !   Returns zeta number of a
    !   basis orbital 

    !    INTEGER  IO   : Orbital index (within atom)
    !                    IO > 0 => Basis orbitals
    !***********************OUTPUT*****************************************
    !   INTEGER ZETAFIO  : Zeta number of orbital

   
    call chk('mofio',is)
    zetafio = 0
    spp => species(is)

    if (io.gt.0) then
       if (io.gt. get_number_of_orbs(spp))  call die("zetafio: No such orbital")
       zetafio = get_orb_zeta(spp,io)
    else 
       call die('zetafio only deals with orbitals')
    endif

  end function zetafio

  !------------------------------------------------------------------------

  SUBROUTINE zphiatm(is,func,io,r,zphi,grzphi)
    !     Calculates z*phiatm and its gradient

    integer, intent(in)   :: is, io
    type(func_t), intent(in) :: func  ! Function label
    real(dp), intent(in)  :: r(3)
    real(dp), intent(out) :: zphi, grzphi(3)

    real(dp) phi, grphi(3), z
    call phiatm(is,func,io,r,phi,grphi)
    z = r(3)
    zphi = z * phi
    grzphi(1) = z * grphi(1)
    grzphi(2) = z * grphi(2)
    grzphi(3) = z * grphi(3) + phi
  END SUBROUTINE zphiatm

  !-------------------------------------------------------------------------
  
  FUNCTION zvalfis( IS )
    real(dp) :: zvalfis          ! Valence charge
    integer, intent(in) :: is            ! Species index

    call chk('zvalfis',is)

    zvalfis= get_valence_charge(species(is))
  end function zvalfis


 
end module atmfuncs





