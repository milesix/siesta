!     This file is part of the SIESTA package.
!     
!     Copyright (c) Fundacion General Universidad Autonoma de Madrid:
!     E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal
!     and J.M.Soler, 1996-2006.
!     
!     Use of this software constitutes agreement with the full conditions
!     given in the SIESTA license, as signed by all legitimate users.
!     

module atom_generation_types
  !
  !=======================================================================
  !
  !     This module defines data structures to handle the specification
  !     of the basis set and KB projectors in SIESTA. This specification
  !     is read from the fdf file by routine 'basis_specs.basis_read', and
  !     then converted to the form expected by routine 'atom.atom_main'
  !
  !     At present, for historical reasons, these data structures are 
  !     different from those in module 'atm_types', but in principle could
  !     be made to contain the same information (except for the mapping
  !     of 'nl' to 'nlm' orbitals and projectors and the polarization
  !     orbitals) by the use of the 'rad_func' pointers indicated within 
  !     types 'shell_t' and 'kbshell_t'. (See module 'radial')
  !     Almost all the subroutines are public.

  use precision, only : dp
  use atomparams, only :
  use pseudopotential_new, only : pseudopotential_new_t
  use radial, only : rad_func_t, rad_cutoff

  Implicit None

  type ground_state_t
     integer                   ::  lmax_valence
     integer                   ::  n(0:3)
     real(dp)                  ::  occupation(0:3)
     logical                   ::  occupied(0:4)   ! note 0..4
     real(dp)                  ::  z_valence
  end type ground_state_t

  type energies_t
     real(dp) :: ekin,ionic,screened,total
  end type energies_t

  type shell_t
     integer                   ::  n          ! n quantum number
     integer                   ::  l          ! angular momentum
     integer                   ::  i_sm       ! This is the ith shell in the semicore list. 
     integer                   ::  nzeta      ! Number of PAOs
     logical                   ::  polarizes  !Does this shell polarize another?  
     logical                   ::  polarized = .false.  !Is this shell polarized?
     logical                   ::  auto_polarized = .false. !Is this shell polarized automagically?
     logical                   ::  auto_polarizes = .false. !Is this shell polarizing automagically?
     integer                   ::  nzeta_pol = 0
     real(dp)                  ::  split_norm = 0.0_dp ! Split norm value
     logical                   ::  split_norm_specified

     real(dp), pointer         :: spln(:) => Null()    ! Split norm used for each z

     real(dp)                  ::  rinn       ! Soft confinement
     real(dp)                  ::  vcte       ! Soft confinement
     real(dp), pointer         ::  rc(:) => Null()     ! rc's for PAOs
     real(dp), pointer         ::  lambda(:)  => Null() ! Contraction factors.
     real(dp), pointer         :: population(:)  => Null()!Population of each orb

     character(len=20)         :: multiple_z_kind !split, splitgauss, nodes, etc
     type(rad_func_t), pointer :: orb(:)  => Null()  ! Actual orbitals. 
     type(rad_func_t), pointer :: rphi(:)  => Null() ! Orb multiplied by r
     type(rad_func_t)          :: split_table !Norm for multiple z generation
     type(rad_func_t)          :: pseudo  ! Pseudopotential for this shell.    
    
     type(rad_func_t)          ::  ve   ! Potential used during the pseudo generation
     type(rad_func_t)          ::  ve_PAO ! Potential used for PAOs. May include a scaling factor

     integer                   ::  l_shell_polarized !l of the shell being polarized.
     integer                   ::  n_shell_polarized
     integer                   ::  l_shell_polarizes 
     type(energies_t), pointer :: energies(:)  !The different energies (kin, pot, screen.) of the orbs.
     real(dp)                  :: z_valence  !The z_valence is needed every time we solve schroed. eq.
  end type shell_t

  type lshell_t
     integer                   ::  l          ! angular momentum
     integer                   ::  nn         ! number of n's for this l
     type(shell_t), pointer    ::  shell(:)   ! One shell for each n
  end type lshell_t

  type kbshell_t
     integer                   ::  l          ! angular momentum
     integer                   ::  nkbl       ! No. of projs for this l
     real(dp), pointer         ::  erefkb(:)  ! Reference energies
  end type kbshell_t
  ! DSP_LDAU
  ! New data for the generation of ldaU projectors. This is similar
  ! to the generation of a SZ basis set (tipically with shorter radii)
  type ldaushell_t
     integer                   ::  n          ! n quantum number
     integer                   ::  l          ! angular momentum
     integer                   ::  i_sm       ! This is the ith shell in the semicore list. 
     integer                   ::  method     ! generation method (1 or 2)
     real(dp)                  ::  rinn       ! Soft confinement
     real(dp)                  ::  vcte       ! Soft confinement
     real(dp)                  ::  rc         ! rc's for PAOs
     real(dp)                  ::  lambda     ! Contraction factors.
     real(dp)                  ::  U          ! LDAU_U parameter
     real(dp)                 ::  J          ! LDAU_J parameter
     real(dp)                 ::  width      ! width for Fermi distribution 
                                              ! to cut the projector
     type(rad_func_t)          :: ldaupj ! Actual projector
     type(rad_func_t)          :: rldaupj ! Actual projector
     type(rad_func_t)          :: pseudo  ! Pseudopotential for this shell.    

     type(rad_func_t)          ::  ve   ! Potential used during the pseudo generation
     type(rad_func_t)          ::  ve_PAO ! Potential used for PAOs. May include a scaling factor
     real(dp)                  :: z_valence  !The z_valence is needed every time we solve schroed. eq.
  end type ldaushell_t

  type l_ldaushell_t
     integer                   ::  l          ! angular momentum
     integer                   ::  nn         ! number of n's for this l
     type(ldaushell_t), pointer    ::  ldaushell(:)   ! One shell for each n
  end type l_ldaushell_t
  !
  !     Main data structure
  !
  type basis_def_t
     character(len=20)         ::  label      ! Long label
     character(len=2)          ::  symbol     ! Atomic symbol
     integer                   ::  z          ! Atomic number
     type(ground_state_t)      ::  ground_state
     type(pseudopotential_new_t) ::  pseudopotential
     integer                   ::  lmxo       ! Max l for basis
     integer                   ::  lmxkb      ! Max l for KB projs
     integer                   ::  lmxldaupj  ! Max l for ldaU projs
     type(lshell_t), pointer   ::  lshell(:)  ! One shell per l 
     type(kbshell_t), pointer  ::  kbshell(:) ! One KB shell per l
     type(l_ldaushell_t), pointer::  l_ldaushell(:) ! One ldaU proj shell per l
     real(dp)                  ::  ionic_charge
     real(dp)                  ::  mass   
     !
     ! The rest of the components are auxiliary
     ! 
     logical                   ::  floating   
     logical                   ::  bessel
     logical                   ::  synthetic
     character(len=20)         ::  basis_type
     character(len=20)         ::  basis_size
     logical                   ::  semic      ! 
     integer                   ::  nshells_tmp
     integer                   ::  nkbshells
     integer                   ::  nldaushells
     integer                   ::  lmxkb_requested
     type(shell_t), pointer    ::  tmp_shell(:)
     type(ldaushell_t), pointer::  tmp_ldaushell(:)
  end type basis_def_t

  real(dp)                     :: filter_cutoff = 0.0_dp
  !integer, save                      :: nsp  ! Number of species
  type(basis_def_t),allocatable, save, target :: basis_parameters(:)

  !-----------------------------------------------------------------------

  interface destroy
     module procedure destroy_shell,destroy_lshell, destroy_basis_def, destroy_ldaushell, destroy_l_ldaushell

  end interface
  interface initialize
     module procedure init_shell, init_kbshell,  init_ldaushell, init_l_ldaushell, init_lshell, init_basis_def

  end interface

CONTAINS

  !--------------------------------------------------------------------------

  subroutine copy_shell(source,target)
    !
    !     This is a "deep-copy" of a structure, including the 
    !     *allocated memory* associated to the internal pointer components.
    !
    type(shell_t), intent(in)          :: source
    type(shell_t), intent(out)         :: target

    target%l = source%l
    target%n = source%n
    target%nzeta = source%nzeta
    target%polarizes = source%polarizes
    target%polarized    = source%polarized
    target%nzeta_pol = source%nzeta_pol
    target%rinn = source%rinn
    target%vcte = source%vcte
    target%split_norm = source%split_norm
    allocate(target%rc(1:size(source%rc)))
    allocate(target%lambda(1:size(source%lambda)))
    allocate(target%orb(1:size(source%orb)))
    allocate(target%energies(1:size(source%energies)))
    allocate(target%spln(1:size(source%spln)))
    allocate(target%population(1:size(source%population)))
    allocate(target%rphi(1:size(source%rphi)))
    target%spln(:)      = source%spln(:)
    target%rc(:)        = source%rc(:)
    target%lambda(:)    = source%lambda(:)
    target%orb          = source%orb
    target%energies(:)  = source%energies(:)
    target%population   = source%population
    target%rphi         = source%rphi
  end subroutine copy_shell
  !-----------------------------------------------------------------------

  subroutine copy_ldaushell(source,target)
    !
    !     This is a "deep-copy" of a structure, including the 
    !     *allocated memory* associated to the internal pointer components.
    !
    type(ldaushell_t), intent(in)          :: source
    type(ldaushell_t), intent(out)         :: target

    target%l = source%l
    target%n = source%n
    target%method=source%method
    target%i_sm = source%i_sm
    target%rinn = source%rinn
    target%vcte = source%vcte
    target%rc = source%rc
    target%width = source%width
    target%lambda = source%lambda
    target%U = source%U
    target%J = source%J
    target%ldaupj    = source%ldaupj
    target%rldaupj   = source%rldaupj
  end subroutine copy_ldaushell
  !-----------------------------------------------------------------------

  subroutine init_shell(p)
    type(shell_t)          :: p

    p%l = -1
    p%n = -1
    p%nzeta = 0
    p%polarizes = .false.
    p%polarized = .false.
    p%nzeta_pol = -1
    p%rinn = 0._dp
    p%vcte = 0._dp
    p%split_norm = 0.0_dp
    p%split_norm_specified = .false.
    nullify(p%orb,p%rphi,p%energies,p%rc,p%lambda)
  end subroutine init_shell

  !-----------------------------------------------------------------------

  subroutine init_ldaushell(p)
    type(ldaushell_t)          :: p

    p%l = -1
    p%n = -1
    p%method= 2
    p%rinn = 0._dp
    p%vcte = 0._dp
    p%rc   = 0._dp
    p%width = 0._dp
    p%lambda= 0._dp
    p%U    = 0._dp
    p%J    = 0._dp
!   nullify(p%ldaupj,p%rldaupj)
  end subroutine init_ldaushell


  !-----------------------------------------------------------------------
  subroutine init_kbshell(p)
    type(kbshell_t)          :: p
    p%l = -1
    p%nkbl = -1
    nullify(p%erefkb)
  end subroutine init_kbshell

  !-----------------------------------------------------------------------
  subroutine init_lshell(p)
    type(lshell_t)          :: p

    p%l = -1
    p%nn = 0
    nullify(p%shell)
  end subroutine init_lshell

  !-----------------------------------------------------------------------
  subroutine init_l_ldaushell(p)
      type(l_ldaushell_t)       :: p
      p%l = -1
      p%nn = 0
      nullify(p%ldaushell)
  end subroutine init_l_ldaushell
  !-----------------------------------------------------------------------

  subroutine init_basis_def(p)
    type(basis_def_t)          :: p

    p%lmxo = -1
    p%lmxkb = -1
    p%lmxldaupj = -1
    p%lmxkb_requested = -1
    p%nkbshells = -1
    p%nshells_tmp = -1
    p%nldaushells = -1
    p%label = 'Unknown'
    p%semic = .false.
    p%ionic_charge = 0.d0
    nullify(p%tmp_shell)
    nullify(p%lshell)
    nullify(p%kbshell)
    nullify(p%l_ldaushell)
    nullify(p%tmp_ldaushell)
  end subroutine init_basis_def

  !-----------------------------------------------------------------------
  subroutine destroy_shell(p)
    type(shell_t), pointer   :: p(:)

    integer                  :: i
    type(shell_t), pointer   :: q

    if (.not. associated(p)) return
    do i = 1, size(p)
       q=>p(i)
       deallocate(q%lambda,q%population)
       !if (rad_allocated(q%ve))     call rad_dealloc(q%ve)
       !if (rad_allocated(q%ve_pao)) call rad_dealloc(q%ve_pao)
       !print *, "ve, ve_pao"
       !do j=1,size(q%orb)
       !   if (rad_allocated(q%orb(j)))  call rad_dealloc(q%orb(j))
       !   if (rad_allocated(q%rphi(j))) call rad_dealloc(q%rphi(j))
       !enddo
       !print *, "orb, rphi"
       if (associated(q%orb)) deallocate(q%orb)
       if (associated(q%rphi)) deallocate(q%rphi)
       if (associated(q%energies)) deallocate(q%energies)
       if (associated(q%rc)) deallocate(q%rc)
       if (associated(q%spln)) deallocate(q%spln)
    enddo
    deallocate(p)
  end subroutine destroy_shell

  !-----------------------------------------------------------------------
  subroutine destroy_ldaushell(p)
    type(ldaushell_t), pointer   :: p(:)

    integer                  :: i, j
    type(ldaushell_t), pointer   :: q

    if (.not. associated(p)) return
    do i = 1, size(p)
       q=>p(i)
!       if (associated(q%ldaupj)) deallocate(q%ldaupj)
!       if (associated(q%rldaupj)) deallocate(q%rldaupj)
    enddo
    deallocate(p)
  end subroutine destroy_ldaushell

  !-----------------------------------------------------------------------
  subroutine destroy_lshell(p)
    type(lshell_t), pointer   :: p(:)

    integer i
    type(lshell_t), pointer   :: q

    if (.not. associated(p)) return
    do i = 1, size(p)
       q=>p(i)
       call destroy_shell(q%shell)
    enddo
    deallocate(p)
  end subroutine destroy_lshell

  !-----------------------------------------------------------------------

  subroutine destroy_l_ldaushell(p)
    type(l_ldaushell_t), pointer   :: p(:)

    integer i
    type(l_ldaushell_t), pointer   :: q

    if (.not. associated(p)) return
    do i = 1, size(p)
       q=>p(i)
       call destroy_ldaushell(q%ldaushell)
    enddo
    deallocate(p)
  end subroutine destroy_l_ldaushell

  !-----------------------------------------------------------------------
  subroutine destroy_basis_def(p)
    type(basis_def_t)          :: p

    call destroy_lshell(p%lshell)
    call destroy_shell(p%tmp_shell)
    call destroy_l_ldaushell(p%l_ldaushell)
    call destroy_ldaushell(p%tmp_ldaushell)
  end subroutine destroy_basis_def

  !-----------------------------------------------------------------------
  subroutine print_shell(p)
    type(shell_t)            :: p

    integer i

    write(6,*) 'SHELL-------------------------'
    write(6,'(5x,a20,i20)') 'Angular momentum',     p%l
    write(6,'(5x,a20,i20)') 'n quantum number',     p%n
    write(6,'(5x,a20,i20)') 'Nzeta'           ,     p%nzeta
    write(6,'(5x,a20,l20)') 'Polarized?       ',    p%polarized
    write(6,'(5x,a20,i20)') 'Nzeta pol'           , p%nzeta_pol
    write(6,'(5x,a20,g20.10)') 'rinn'           , p%rinn
    write(6,'(5x,a20,g20.10)') 'vcte'           , p%vcte
    write(6,'(5x,a)') 'rc and lambda for each nzeta:'
    do i = 1, p%nzeta
       write(6,'(5x,i2,2x,2g20.10)') i, rad_cutoff(p%orb(i)), p%lambda(i)
    enddo
    write(6,*) '--------------------SHELL'

  end subroutine print_shell

  !-----------------------------------------------------------------------
  subroutine print_kbshell(p)
    type(kbshell_t)            :: p

    integer i

    write(6,*) 'KBSHELL-------'
    write(6,'(5x,a20,i20)') 'Angular momentum', p%l
    write(6,'(5x,a20,i20)') 'number of projs',p%nkbl
    write(6,'(5x,a)') 'ref energy for each proj:'
    do i = 1, p%nkbl
       write(6,'(5x,i2,2x,g20.5)') i, p%erefkb(i)
    enddo
    write(6,*) '---------------------KBSHELL'

  end subroutine print_kbshell

  !-----------------------------------------------------------------------
  subroutine print_lshell(p)
    type(lshell_t) :: p

    integer i

    write(6,*) 'LSHELL --------------'
    write(6,'(5x,a20,i20)') 'Angular momentum',p%l
    write(6,'(5x,a20,i20)') 'Number of n shells',p%nn
    if (.not. associated(p%shell)) return
    do i=1, p%nn
       call print_shell(p%shell(i))
    enddo
    write(6,*) '--------------LSHELL'
  end subroutine print_lshell

  !-----------------------------------------------------------------------
  subroutine print_basis_def(p)
    type(basis_def_t) :: p

    integer i

    write(6,*) ' '
    write(6,*) 'SPECIES---'

    write(6,'(5x,a20,a20)') 'label',p%label
    write(6,'(5x,a20,i20)') 'atomic number',p%z
    write(6,'(5x,a20,a20)') 'basis type',p%basis_type
    write(6,'(5x,a20,a20)') 'basis size',p%basis_size
    write(6,'(5x,a20,g20.10)') 'ionic charge', p%ionic_charge
    write(6,'(5x,a20,i20)') 'lmax basis', p%lmxo

    if (.not. associated(p%lshell)) then
       write(6,*) 'No L SHELLS, lmxo=', p%lmxo
       return
    endif
    do i=0, p%lmxo
       call print_lshell(p%lshell(i))
    enddo
    if (.not. associated(p%kbshell)) then
       write(6,*) 'No KB SHELLS, lmxkb=', p%lmxkb
       return
    endif
    do i=0, p%lmxkb
       call print_kbshell(p%kbshell(i))
    enddo

    write(6,*) '------------SPECIES'
    write(6,*) 

  end subroutine print_basis_def

  !-----------------------------------------------------------------------
  
  subroutine write_basis_specs(lun,is)
    integer, intent(in):: lun
    integer, intent(in):: is
    type (basis_def_t), pointer :: basp
    type (shell_t)    , pointer :: nshell
    type (ldaushell_t), pointer :: nldaushell
    integer l,n,i

    basp => basis_parameters(is)

    write(lun,'(/a/79("="))') '<basis_specs>'
    write(lun,'(a70)') "="
    write(lun,'(a10,1x,a2,i4,4x,a5,g12.5,4x,a7,g12.5)') basp%label,'Z=',basp%z, &
         'Mass=', basp%mass, 'Charge=',basp%ionic_charge
    write(lun,'(a5,i1,1x,a6,i2,5x,a10,a10,1x,a6,l1)') &
         'Lmxo=',basp%lmxo, 'Lmxkb=', basp%lmxkb, &
         'BasisType=', basp%basis_type, 'Semic=', basp%semic
   
    do l=0,basp%lmxo
       write(lun,'(a2,i1,2x,a7,i1)') 'L=',l,'Nsemic=', basp%lshell(l)%nn-1
       do n=1,basp%lshell(l)%nn
          nshell => basp%lshell(l)%shell(n)
          if(nshell%nzeta == 0) exit
          write(lun,'(10x,a2,i1,2x,a6,i1,2x,a7,l2)') &
               'n=',nshell%n, 'nzeta=',nshell%nzeta, &
               'polorb=',nshell%polarizes
           write(lun,'(10x,a10,2x,g12.5)')  &
                'splnorm:', nshell%split_norm
          write(lun,'(10x,a10,2x,g12.5)') &
               'vcte:', nshell%vcte
          write(lun,'(10x,a10,2x,g12.5)') &
               'rinn:', nshell%rinn
          write(lun,'(10x,a10,2x,4g12.5)') 'rcs:',&
               (nshell%rc(i),i=1,nshell%nzeta)
          write(lun,'(10x,a10,2x,4g12.5)') 'lambdas:',&
               (nshell%lambda(i),i=1,nshell%nzeta)
       enddo
    enddo
    write(lun,"(a70)") "-"
    if (basp%lmxkb > 0)then
       do l=0,basp%lmxkb
          write(lun,'(a2,i1,2x,a5,i1,2x,a6,4g12.5)')&
               'L=', l, 'Nkbl=', basp%kbshell(l)%nkbl,&
               'erefs:  ', (basp%kbshell(l)%erefkb(i),i=1,&
               basp%kbshell(l)%nkbl)
       enddo
    endif
    write(lun,"(a70)") "-"
    if (basp%lmxldaupj > -1)then
      do l=0,basp%lmxldaupj
       if(basp%l_ldaushell(l)%nn == 0) cycle
       write(lun,'(a2,i1,2x,a12,i1)') 'L=',l,'Nldau_semic=', basp%l_ldaushell(l)%nn
       do n=1,basp%l_ldaushell(l)%nn
          nldaushell => basp%l_ldaushell(l)%ldaushell(n)
          write(lun,'(10x,a2,i1)') &
               'n=',nldaushell%n
          write(lun,'(10x,a10,2x,2g12.5)') &
               'U, J=:', nldaushell%U, nldaushell%J
          write(lun,'(10x,a10,2x,g12.5)') &
               'vcte:', nldaushell%vcte
          write(lun,'(10x,a10,2x,g12.5)') &
               'rinn:', nldaushell%rinn
          write(lun,'(10x,a10,2x,g12.5)') 'rcs:',&
               nldaushell%rc
          write(lun,'(10x,a10,2x,g12.5)') 'lambdas:',&
               nldaushell%lambda
       enddo
     enddo
    endif
    write(lun,'(a70)')"="
    write(lun,'(a)') '</basis_specs>'

  end subroutine write_basis_specs
  !-----------------------------------------------------------------------
end module atom_generation_types













