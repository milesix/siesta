!     This file is part of the SIESTA package.
!     
!     Copyright (c) Fundacion General Universidad Autonoma de Madrid:
!     E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal
!     and J.M.Soler, 1996-2006.
!     
!     Use of this software constitutes agreement with the full conditions
!     given in the SIESTA license, as signed by all legitimate users.
!     

module atm_types
  use precision, only : dp
  use radial
  use sys, only : die
  use hilbert_vector_collection
  
 
  implicit none

  !length of character strings.
  integer, parameter :: symbol_length = 2
  integer, parameter :: label_length = 20

  real(dp), parameter :: tiny20=1.e-20_dp
  real(dp), parameter :: tiny12=1.e-12_dp


  !
  !     Species_info: Consolidate all the pieces of information in one place
  !
  type species_info_t
     private
     logical                                    ::  fake = .false. !real species or not
     character(len=symbol_length)               ::  symbol
     character(len=label_length)                ::  label
     real(dp)                                   ::  atomic_number          
     real(dp)                                   ::  mass
     real(dp)                                   ::  valence_charge   
     real(dp)                                   ::  self_energy !Electrostatic

     type(rad_func_t), pointer                  ::  reduced_vlocal => Null()
     type(rad_func_t), pointer                  ::  neutral_atom_potential  => Null()
     type(rad_func_t), pointer                  ::  pseudo_local_charge  => Null()

     !Core charge for nonlinear core corrections
     type(rad_func_t), pointer                  ::  core_charge  => Null()
     logical                                    ::  has_core_charge
      
     type(hilbert_vector_collection_t)          ::  orbs
     type(hilbert_vector_collection_t), pointer ::  kb_proj, ldau_proj  => Null()
     logical                                    ::  read_from_file
  end type species_info_t

  !
  integer, save, public             :: nspecies
  integer, save, public             :: npairs

  real(dp)                          :: orbs_kc_max = 0.0_dp  !Maximum cutoff of all orbs in kspace

  type(species_info_t), allocatable, target, save   ::  species(:) 

!    Radial function with the difference between the electrostatic energy 
!    of two spherical charge-densities and two punctual charges with the 
!    same total charge as a function of the distance between the centers 
!    of these charge densities. 
  type(rad_func_t), allocatable, target, save   ::  elec_corr(:)
  !
  !-------------------------------------------------------------------------
 
contains

  subroutine broadcast_basis
!
!     Globalizes the basis and potentials data structure
!     Alberto Garcia, June 2000--

      use radial
      use parallel,   only : Node, Nodes

#ifdef MPI
      use mpi_siesta
#endif

      implicit none

#ifndef MPI
!
!     Do nothing...
!
    end subroutine broadcast_basis
#else
    integer MPIerror

    integer is, i
    type(species_info_t), pointer  :: spp => Null()
    logical :: kbs= .false., vna=.false., vlocal = .false., ldau = .false.
    logical :: pseudo_charge = .false.
    
    if (Nodes.eq.1) return
      
    call MPI_Bcast(nspecies,1,MPI_integer,0,MPI_Comm_World,MPIerror)
    call MPI_Bcast(npairs,1,MPI_integer,0,MPI_Comm_World,MPIerror)
    
    if (Node.ne.0) allocate (species(nspecies))

    print *, "species"
    do is=1,nspecies
       spp => species(is)
       call MPI_Bcast(spp%symbol,symbol_length,MPI_character,&
            0,MPI_Comm_World,MPIerror)
       call MPI_Bcast(spp%label,label_length,MPI_character, &
            0,MPI_Comm_World,MPIerror)
       call MPI_Bcast(spp%atomic_number,1,MPI_integer, &
            0,MPI_Comm_World,MPIerror)
       call MPI_Bcast(spp%mass,1,MPI_double_precision, &
            0,MPI_Comm_World,MPIerror)
       call MPI_Bcast(spp%valence_charge,1,MPI_double_precision, &
            0,MPI_Comm_World,MPIerror)
       call MPI_Bcast(spp%self_energy,1,MPI_double_precision, &
            0,MPI_Comm_World,MPIerror)
       call MPI_Bcast(spp%fake,1,MPI_logical,0,MPI_Comm_World,MPIerror)

       print *, "orbs begin"
       call broadcast_hilbert_vector_collection(spp%orbs)
       print *, "orbs done"
       if (Node .eq. 0)  kbs = has_kbs(spp)
       call MPI_Bcast(kbs,1,MPI_logical,0,MPI_Comm_World,MPIerror)
       if (Node .ne. 0 .and. kbs) allocate(spp%kb_proj)       
       if (kbs) call broadcast_hilbert_vector_collection(spp%kb_proj)
       
       if (Node .eq. 0)  ldau = associated(spp%ldau_proj)
       call MPI_Bcast(ldau,1,MPI_logical,0,MPI_Comm_World,MPIerror)
       if (Node .ne. 0 .and. ldau) allocate(spp%ldau_proj)
       if (ldau) call broadcast_hilbert_vector_collection(spp%ldau_proj)

       if (Node .eq. 0) vlocal = has_reduced_vlocal(spp) 
       call MPI_Bcast(vlocal,1,MPI_logical,0,MPI_Comm_World,MPIerror)
       if (Node .ne.0  .and. vlocal) allocate(spp%reduced_vlocal)
       if (vlocal) call rad_broadcast(spp%reduced_vlocal)

       if (Node .eq. 0) vna = has_neutral_atom_potential(spp) 
       call MPI_Bcast(vna,1,MPI_logical,0,MPI_Comm_World,MPIerror)
       if (Node .ne.0  .and. vna) allocate(spp%neutral_atom_potential)
       if (vna) call rad_broadcast(spp%neutral_atom_potential)

       if (Node .eq. 0) pseudo_charge = has_pseudo_local_charge(spp) 
       call MPI_Bcast(pseudo_charge,1,MPI_logical,0,MPI_Comm_World,MPIerror)
       if (Node .ne.0  .and. pseudo_charge) allocate(spp%pseudo_local_charge)
       if (pseudo_charge) call rad_broadcast(spp%pseudo_local_charge)
       
       call MPI_Bcast(spp%has_core_charge,1,MPI_logical,0,MPI_Comm_World,MPIerror)
       if (spp%has_core_charge) then
          if (Node .ne. 0) allocate(spp%core_charge)
          call rad_broadcast(spp%core_charge)
          if (Node .eq. 0) print *, "core charge"
       endif

       call MPI_Bcast(spp%read_from_file,1,MPI_logical,0,MPI_Comm_World,MPIerror)

  enddo

  if (Node .ne. 0) allocate(elec_corr(1:npairs))

  do i=1,npairs
     call rad_broadcast(elec_corr(i))
  end do
 
  
end subroutine broadcast_basis


#endif

!-----------------------------------------------------------------------

  subroutine get_elec_corr(index,r,energ,dedr)
    integer,   intent(in):: index
    real(dp), intent(in) :: r
    real(dp), intent(out):: energ
    real(dp), intent(out):: dedr

    real(dp) r_local

    energ = 0.0_dp
    dedr  = 0.0_dp

    if ( r .gt. rad_cutoff(elec_corr(index)) - tiny12 ) return

    call rad_get(elec_corr(index),r,energ,dedr)
    r_local = r+tiny20
    energ=2.0_dp*energ/r_local
    dedr=(-energ + 2.0_dp*dedr)/r_local

  end subroutine get_elec_corr

  !-----------------------------------------------------------------------
  
  function get_lmax_orbs(species) result(lmax)
    type(species_info_t), intent(in) :: species
    integer :: lmax
    lmax = get_lmax(species%orbs)
  end function get_lmax_orbs

  !-------------------------------------------------------------------------

  subroutine set_lmax_orbs(species,lmax)
    type(species_info_t), intent(inout) :: species
    integer, intent(in) :: lmax
    call set_lmax(species%orbs,lmax)
  end subroutine set_lmax_orbs

  !-------------------------------------------------------------------------
  
  function get_lmax_kb_proj(species) result(lmax)
    type(species_info_t), intent(in) :: species
    integer :: lmax
    if (associated(species%kb_proj)) then
       lmax = get_lmax(species%kb_proj)
    else
       lmax = -1
    endif
  end function get_lmax_kb_proj

  !-------------------------------------------------------------------------

  subroutine set_lmax_kb_proj(species, lmax)
    type(species_info_t), intent(inout) :: species
    integer,               intent(in) :: lmax
    call set_lmax(species%kb_proj, lmax)
  end subroutine set_lmax_kb_proj

  !-------------------------------------------------------------------------

  function get_lmax_ldau_proj(species) result(lmax)
    type(species_info_t), intent(in) :: species
    integer :: lmax
    if (associated(species%ldau_proj)) then
       lmax = get_lmax(species%ldau_proj)
    else
       lmax = -1
    end if
  end function get_lmax_ldau_proj

  !-------------------------------------------------------------------------

  function is_floating(species) result(floating)
    type(species_info_t), intent(in) :: species
    logical :: floating
    floating = species%fake
  end function is_floating

  !-------------------------------------------------------------------------

  subroutine set_floating(species, floating)
   type(species_info_t), intent(inout) :: species
    logical, intent(in) :: floating
    species%fake = floating
  end subroutine set_floating

  !-------------------------------------------------------------------------

  subroutine init_orbs(species,norbs)
    type(species_info_t), intent(inout) :: species
    integer, intent(in) :: norbs
    call allocate_collection(species%orbs,norbs)
  end subroutine init_orbs
 
  !-------------------------------------------------------------------------

  subroutine init_kb_proj(species,nkb)
    type(species_info_t), intent(inout) :: species
    integer, intent(in) :: nkb
    allocate(species%kb_proj)
    call allocate_collection(species%kb_proj,nkb)
  end subroutine init_kb_proj

  !-------------------------------------------------------------------------

  subroutine init_ldau_proj(species)
    type(species_info_t), intent(inout) :: species
    allocate(species%ldau_proj)
    stop "atm_types: init_ldau_proj"
  end subroutine init_ldau_proj

  !-------------------------------------------------------------------------

  function get_symbol(species)
    type(species_info_t), intent(inout) :: species
    character(len=symbol_length) :: get_symbol
    get_symbol = species%symbol
  end function get_symbol

  !-------------------------------------------------------------------------

  subroutine set_symbol(species,symbol)
    type(species_info_t), intent(inout) :: species
    character(len=symbol_length), intent(in) :: symbol
    species%symbol = symbol
  end subroutine set_symbol

  !--------------------------------------------------------------------------

  function get_label(species)
    type(species_info_t), intent(inout) :: species
    character(len=label_length) :: get_label
    get_label = species%label
  end function get_label

  !-------------------------------------------------------------------------

  subroutine set_label(species,label)
    type(species_info_t), intent(inout) :: species
    character(len=label_length), intent(in) :: label
    species%label = label
  end subroutine set_label

  !--------------------------------------------------------------------------

  function get_atomic_number(species)
    type(species_info_t), intent(inout) :: species
    integer:: get_atomic_number
    get_atomic_number = species%atomic_number
  end function get_atomic_number

  !-------------------------------------------------------------------------

  subroutine set_atomic_number(species,z)
    type(species_info_t), intent(inout) :: species
    integer, intent(in) :: z
    species%atomic_number = z
  end subroutine set_atomic_number
  
  !-------------------------------------------------------------------------

  function get_mass(species)
    type(species_info_t), intent(inout) :: species
    real(dp):: get_mass
    get_mass = species%mass
  end function get_mass

  !-------------------------------------------------------------------------

  subroutine set_mass(species,mass)
    type(species_info_t), intent(inout) :: species
    real(dp), intent(in) :: mass
    species%mass = mass
  end subroutine set_mass

  !--------------------------------------------------------------------------

  function get_valence_charge(species)
    type(species_info_t), intent(inout) :: species
    real(dp):: get_valence_charge
    get_valence_charge = species%valence_charge
  end function get_valence_charge

  !-------------------------------------------------------------------------

  subroutine set_valence_charge(species,charge)
    type(species_info_t), intent(inout) :: species
    real(dp), intent(in) :: charge
    species%valence_charge = charge
  end subroutine set_valence_charge

  !--------------------------------------------------------------------------
  
  function get_self_energy(species)
    type(species_info_t), intent(inout) :: species
    real(dp):: get_self_energy
    get_self_energy = species%self_energy
  end function get_self_energy

  !-------------------------------------------------------------------------
  subroutine set_no_kb(species) 
    type(species_info_t), intent(inout) :: species
    nullify(species%kb_proj)
  end subroutine set_no_kb

  !-------------------------------------------------------------------------
  subroutine set_no_neutral_atom_potential(species) 
    type(species_info_t), intent(inout) :: species
    nullify(species%neutral_atom_potential)
  end subroutine set_no_neutral_atom_potential

  !-------------------------------------------------------------------------

  subroutine set_self_energy(species,self_energy)
    type(species_info_t), intent(inout) :: species
    real(dp), intent(in) :: self_energy
    species%self_energy = self_energy 
  end subroutine set_self_energy

  !--------------------------------------------------------------------------
  
  function get_reduced_vlocal(species) result(vloc)
    type(species_info_t), intent(inout) :: species
    type(rad_func_t) :: vloc
    call rad_copy(species%reduced_vlocal,vloc)
  end function get_reduced_vlocal

  !-------------------------------------------------------------------------

  subroutine set_no_reduced_vlocal(species)
    type(species_info_t), intent(inout) :: species
    nullify(species%reduced_vlocal)
  end subroutine set_no_reduced_vlocal

  !--------------------------------------------------------------------------

  subroutine set_reduced_vlocal(species,vlocal)
    type(species_info_t), intent(inout) :: species
    type(rad_func_t), intent(in)  :: vlocal
    type(rad_func_t)              :: rad_tmp
    
    call rad_copy(vlocal,rad_tmp)
    if(rad_is_log(rad_tmp)) call rad_log_to_linear(rad_tmp)
    allocate(species%reduced_vlocal)
    call rad_copy(rad_tmp,species%reduced_vlocal)
    call rad_dealloc(rad_tmp)
  end subroutine set_reduced_vlocal

  !--------------------------------------------------------------------------
    
  function get_neutral_atom_potential(species) result(vna)
    type(species_info_t), intent(inout) :: species
    type(rad_func_t)                    :: vna
    if (associated(species%neutral_atom_potential)) then
      call rad_copy(species%neutral_atom_potential,vna)
   else
      call die("atm_types: get_neutral_atom_pot. species has no Vna!")
   endif
  end function get_neutral_atom_potential

  !-------------------------------------------------------------------------

  function has_neutral_atom_potential(species) result(has_vna)
    type(species_info_t), intent(inout) :: species
    logical                             :: has_vna 
    has_vna = .true.
    if (is_floating(species)) has_vna = .false.
  end function has_neutral_atom_potential

  !-------------------------------------------------------------------------
  subroutine set_neutral_atom_potential(species,vna)
    type(species_info_t), intent(inout) :: species
    type(rad_func_t), intent(in)  :: vna
    type(rad_func_t)              :: vna_tmp

    call rad_copy(vna,vna_tmp)
    if(rad_is_log(vna)) call rad_log_to_linear(vna_tmp)
    !call rad_fft(vna_tmp,0)
    allocate(species%neutral_atom_potential)
    call rad_copy(vna_tmp,species%neutral_atom_potential)
    call rad_dealloc(vna_tmp)
  end subroutine set_neutral_atom_potential

  !--------------------------------------------------------------------------
  
  function get_pseudo_local_charge(species) result(psloc)
    type(species_info_t), intent(inout) :: species
    type(rad_func_t)                    :: psloc
    call rad_copy(species%pseudo_local_charge,psloc)
  end function get_pseudo_local_charge

  !-------------------------------------------------------------------------

  subroutine set_pseudo_local_charge(species,local)
    type(species_info_t), intent(inout) :: species
    type(rad_func_t), intent(in)  :: local
    type(rad_func_t)              :: rad_tmp
    call rad_copy(local,rad_tmp)
    if(rad_is_log(rad_tmp)) call rad_log_to_linear(rad_tmp)
    allocate(species%pseudo_local_charge)
    call rad_copy(rad_tmp,species%pseudo_local_charge)
    call rad_dealloc(rad_tmp)
  end subroutine set_pseudo_local_charge

  !--------------------------------------------------------------------------
    
  function get_core_charge(species) result(chcore)
    type(species_info_t), intent(inout) :: species
    type(rad_func_t)                    :: chcore
    if (.not. has_core_charge(species)) call die("atm_types: species has no core charge!")
    call rad_copy(species%core_charge,chcore)
  end function get_core_charge

  !-------------------------------------------------------------------------

  subroutine set_core_charge(species,core)
    type(species_info_t), intent(inout) :: species
    type(rad_func_t), intent(in)  :: core

    type(rad_func_t)              :: rad_tmp
    call rad_copy(core,rad_tmp)
    if(rad_is_log(rad_tmp)) call rad_log_to_linear(rad_tmp)
    !call rad_fft(rad_tmp,0)
    allocate(species%core_charge)
    call rad_copy(rad_tmp,species%core_charge)        
    call rad_dealloc(rad_tmp)
  end subroutine set_core_charge

  !--------------------------------------------------------------------------
  
  function has_core_charge(species)
    type(species_info_t), intent(in) :: species
    logical                           :: has_core_charge
    has_core_charge = species%has_core_charge
  end function has_core_charge

  !-------------------------------------------------------------------------

   function has_kbs(species)
    type(species_info_t), intent(in) :: species
    logical                           :: has_kbs
    has_kbs = .true.
    if(is_floating(species)) has_kbs=.false.
  end function has_kbs

  !-------------------------------------------------------------------------

  function has_pseudo_local_charge(species)
    type(species_info_t), intent(in) :: species
    logical                           :: has_pseudo_local_charge
    has_pseudo_local_charge = .true.
    if (is_floating(species)) has_pseudo_local_charge=.false.
  end function has_pseudo_local_charge

  !-------------------------------------------------------------------------

  function has_reduced_vlocal(species)
    type(species_info_t), intent(inout) :: species
    logical                           :: has_reduced_vlocal
    has_reduced_vlocal = .true.
    if (is_floating(species)) has_reduced_vlocal=.false.
  end function has_reduced_vlocal

  !-------------------------------------------------------------------------

  subroutine set_has_core_charge(species,has_core)
    type(species_info_t), intent(inout) :: species
    logical, intent(in)                 :: has_core
    species%has_core_charge = has_core
  end subroutine set_has_core_charge

  !--------------------------------------------------------------------------
  
  function get_read_from_file(species)
    type(species_info_t), intent(inout) :: species
    logical                           :: get_read_from_file
    get_read_from_file = species%read_from_file
  end function get_read_from_file

  !-------------------------------------------------------------------------

  subroutine set_read_from_file(species,file_read)
    type(species_info_t), intent(inout) :: species
    logical                           :: file_read
    species%read_from_file = file_read
  end subroutine set_read_from_file

  !-------------------------------------------------------------------------

  function get_number_of_orbs(species) result(norbs)
    type(species_info_t), intent(in) :: species
    integer                          :: norbs

    norbs = get_n_funcs(species%orbs)

  end function get_number_of_orbs

  !-------------------------------------------------------------------------

  function get_number_of_orbs_non_deg(species) result(norbs)
    type(species_info_t),intent(in) :: species
    integer                         :: norbs

    norbs = get_number_of_vectors(species%orbs)

  end function get_number_of_orbs_non_deg

  !-------------------------------------------------------------------------

  function get_number_of_kb_projs(species) result(nkb)
    type(species_info_t), intent(in) :: species
    integer                        :: nkb

    if (is_floating(species)) then
       nkb = 0
    else
       nkb = get_n_funcs(species%kb_proj)
    endif

  end function get_number_of_kb_projs

  !-------------------------------------------------------------------------

  function get_number_of_kb_non_deg(species) result(nkbs)
    type(species_info_t),intent(in) :: species
    integer                         :: nkbs

    if (is_floating(species)) then
       nkbs = 0
    else
       nkbs = get_number_of_vectors(species%kb_proj)
    endif

  end function get_number_of_kb_non_deg

  !-------------------------------------------------------------------------

  function get_number_of_ldau_proj(species) result(nldau)
    type(species_info_t), intent(in) :: species
    integer                        :: nldau

    if (associated(species%ldau_proj)) then
       nldau = get_n_funcs(species%ldau_proj)
    else
       nldau = 0
    endif

  end function get_number_of_ldau_proj

  !------------------------------------------------------------------------

  function get_rcore(species) result(rcore)
    type(species_info_t), intent(in) :: species
    real(dp)                       :: rcore
    rcore = rad_cutoff(species%core_charge)
  end function get_rcore

  !------------------------------------------------------------------------

  function get_orb(species,io) result(orb)
    type(species_info_t), intent(in) :: species
    integer, intent(in)            :: io

    type(rad_func_t) :: orb
    call rad_copy(get_rad_func_p(species%orbs,io),orb)
  end function get_orb
  
  !------------------------------------------------------------------------

  function get_orb_v(species,io) result(orb)
    type(species_info_t), intent(in) :: species
    integer, intent(in)            :: io

    type(hilbert_vector_t) :: orb
    orb = get_vector(species%orbs,io)
  end function get_orb_v
  
  !------------------------------------------------------------------------

  function get_kb_v(species,io) result(kb)
    type(species_info_t), intent(in) :: species
    integer, intent(in)            :: io
    
    type(hilbert_vector_t) :: kb
    kb = get_vector(species%kb_proj,io)
  end function get_kb_v
  
  !------------------------------------------------------------------------

  function get_orb_pop(species,io) result(pop)
    type(species_info_t), intent(in) :: species
    integer, intent(in)            :: io
    real(dp)                       :: pop
    pop = get_pop(species%orbs,io)
  end function get_orb_pop

  !------------------------------------------------------------------------

  function get_orb_pol(species,io) result(pol)
    type(species_info_t), intent(in) :: species
    integer, intent(in)            :: io
    logical                        :: pol
    pol = get_pol(species%orbs,io)
  end function get_orb_pol

  !------------------------------------------------------------------------

  function get_orb_n(species,io) result(n)
    type(species_info_t), intent(in) :: species
    integer, intent(in)            :: io
    integer                        :: n
    n = get_n(species%orbs,io)
  end function get_orb_n

  !------------------------------------------------------------------------

  function get_orb_l(species,io) result(l)
    type(species_info_t), intent(in) :: species
    integer, intent(in)            :: io
    integer                        :: l
    l = get_l(species%orbs,io)
  end function get_orb_l

  !------------------------------------------------------------------------

  function get_kb_proj_l(species,io) result(l)
    type(species_info_t), intent(in) :: species
    integer, intent(in)            :: io
    integer                        :: l
    l = get_l(species%kb_proj,io)
  end function get_kb_proj_l

  !------------------------------------------------------------------------

  function get_orb_m(species,io) result(m)
    type(species_info_t), intent(in) :: species
    integer, intent(in)            :: io
    integer                        :: m
    m = get_m(species%orbs,io)
  end function get_orb_m

  !------------------------------------------------------------------------

  function get_kb_proj_m(species,io) result(m)
    type(species_info_t), intent(in) :: species
    integer, intent(in)            :: io
    integer                        :: m
    m = get_m(species%kb_proj,io)
  end function get_kb_proj_m
  
  !------------------------------------------------------------------------

  function get_orb_zeta(species,io) result(z)
    type(species_info_t), intent(in) :: species
    integer, intent(in)            :: io
    integer                        :: z
    z = get_zeta(species%orbs,io)
  end function get_orb_zeta

  !------------------------------------------------------------------------

  function get_orb_cutoff(species,io) result(cutoff)
    type(species_info_t), intent(in) :: species
    integer, intent(in)            :: io
    real(dp)                       :: cutoff
    cutoff = get_cutoff(species%orbs,io)
  end function get_orb_cutoff

  !------------------------------------------------------------------------

  function get_kb_proj_cutoff(species,io) result(cutoff)
    type(species_info_t), intent(in) :: species
    integer, intent(in)            :: io
    real(dp)                       :: cutoff
    if (has_kbs(species))then
       cutoff = get_cutoff(species%kb_proj,io)
    else
       cutoff = 0.0_dp
    endif
  end function get_kb_proj_cutoff

  !------------------------------------------------------------------------

  function get_vna_cutoff(species) result(cutoff)
    type(species_info_t), intent(in) :: species
    real(dp)                         :: cutoff
    if (associated(species%neutral_atom_potential)) then
       cutoff = rad_cutoff(species%neutral_atom_potential)
    else
       cutoff = 0.0_dp
    endif
  end function get_vna_cutoff

  !------------------------------------------------------------------------

  function get_local_charge_cutoff(species) result(cutoff)
    type(species_info_t), intent(in) :: species
    real(dp)                       :: cutoff
    cutoff = rad_cutoff(species%pseudo_local_charge)
  end function get_local_charge_cutoff

  !------------------------------------------------------------------------

  function get_kb_proj_energy(species,io) result(energy)
    type(species_info_t), intent(in) :: species
    integer, intent(in)            :: io
    real(dp)                       :: energy
    energy = get_energy(species%kb_proj,io)
  end function get_kb_proj_energy

  !------------------------------------------------------------------------
  
  subroutine get_value_vna(species,r,v,grv)
    type(species_info_t), intent(in) :: species
    real(dp), intent(in)  :: r(3)    ! Point vector, relative to atom
    real(dp), intent(out) :: v       ! Value of local pseudopotential
    real(dp), intent(out) :: grv(3)  ! Gradient of local pseudopotential

    !Internal vars
    type(rad_func_t)              :: vna
    real(dp)                      :: rmod, dvdr

    v=0.0_dp
    grv=0.0_dp
    
    vna  = species%neutral_atom_potential
    rmod = sqrt(sum(r*r))
    rmod = rmod+tiny20
    if (rmod .gt. rad_cutoff(vna)) return

    call rad_get(vna,rmod,v,dvdr)
    grv(1:3) = dvdr * r(1:3)/rmod

  end subroutine get_value_vna

  !------------------------------------------------------------------------

   subroutine get_rvalue_na(species,r,v,grv)
    type(species_info_t), intent(in) :: species
    real(dp), intent(in)  :: r       ! Point vector, relative to atom
    real(dp), intent(out) :: v       ! Value of local pseudopotential
    real(dp), intent(out) :: grv     ! Gradient of local pseudopotential

    !Internal vars
    type(rad_func_t)              :: vna
    real(dp)                      :: rmod

    v=0.0_dp
    grv=0.0_dp
    
    vna  = species%neutral_atom_potential
    rmod = r+tiny20
    if (rmod .gt. rad_cutoff(vna)) return

    call rad_get(vna,rmod,v,grv)

  end subroutine get_rvalue_na

  !------------------------------------------------------------------------

  subroutine get_value_pseudo_local_charge(species,r,v,grv)
    type(species_info_t), intent(in) :: species
    real(dp), intent(in)  :: r(3)    ! Point vector, relative to atom
    real(dp), intent(out) :: v       ! Value of local pseudopotential
    real(dp), intent(out) :: grv(3)  ! Gradient of local pseudopotential

    !Internal vars
    type(rad_func_t)      :: local_charge
    real(dp)              :: rmod, dvdr

    v=0.0_dp
    grv=0.0_dp
    
    local_charge = species%pseudo_local_charge
    rmod = sqrt(sum(r*r))
    rmod = rmod+tiny20

    if (rmod .gt. rad_cutoff(local_charge)) return

    call rad_get(local_charge,rmod,v,dvdr)
    
    grv(1:3) = dvdr * r(1:3)/rmod

  end subroutine get_value_pseudo_local_charge

  !------------------------------------------------------------------------

  subroutine get_value_of_core_charge(species,r,v,grv)
    type(species_info_t), intent(in) :: species
    real(dp), intent(in)  :: r(3)    ! Point vector, relative to atom
    real(dp), intent(out) :: v       ! Value of local pseudopotential
    real(dp), intent(out) :: grv(3)  ! Gradient of local pseudopotential

    !Internal vars
    type(rad_func_t)      :: core_charge
    real(dp)              :: rmod, dvdr

    v=0.0_dp
    grv=0.0_dp
    
    core_charge = species%core_charge
    rmod = sqrt(sum(r*r))
    rmod=rmod+tiny20

    if (rmod .gt. rad_cutoff(core_charge)) return

    call rad_get(core_charge,rmod,v,dvdr)
    
    grv(1:3) = dvdr * r(1:3)/rmod

  end subroutine get_value_of_core_charge

  !------------------------------------------------------------------------

  subroutine get_value_of_orb(species,i,r,v,grv)
    type(species_info_t), intent(in) :: species
    integer, intent(in)   :: i       ! orb index
    real(dp), intent(in)  :: r(3)    ! Point vector, relative to atom
    real(dp), intent(out) :: v       ! Value of local pseudopotential
    real(dp), intent(out) :: grv(3)  ! Gradient of local pseudopotential

    !Internal vars
    type(rad_func_t), pointer     :: orb
    real(dp)                      :: rmod, dvdr

    v=0.0_dp
    grv=0.0_dp
    
    orb => get_rad_func_p(species%orbs,i)
    rmod = sqrt(sum(r*r))
    rmod = rmod+tiny20
    if (rmod .gt. rad_cutoff(orb)) return
    
    call rad_get(orb,rmod,v,dvdr)
    grv(1:3) = dvdr * r(1:3)/rmod

  end subroutine get_value_of_orb

   !------------------------------------------------------------------------

  subroutine get_rvalue_of_orb(species,i,r,v,grv)
    type(species_info_t), intent(in) :: species
    integer, intent(in)   :: i       ! orb index
    real(dp), intent(in)  :: r       ! Point vector, relative to atom
    real(dp), intent(out) :: v       ! Value of local pseudopotential
    real(dp), intent(out) :: grv     ! Gradient of local pseudopotential

    !Internal vars
    type(rad_func_t),pointer      :: orb
    real(dp)                      :: rmod

    v=0.0_dp
    grv=0.0_dp
    
    orb => get_rad_func_p(species%orbs,i)
    rmod = r+tiny20
    if (rmod .gt. rad_cutoff(orb)) return

    call rad_get(orb,rmod,v,grv)

  end subroutine get_rvalue_of_orb
  
  !------------------------------------------------------------------------

  subroutine get_value_of_kb_proj(species,i,r,v,grv)
    type(species_info_t), intent(in) :: species
    integer, intent(in)   :: i       ! kb_proj index
    real(dp), intent(in)  :: r(3)    ! Point vector, relative to atom
    real(dp), intent(out) :: v       ! Value of local pseudopotential
    real(dp), intent(out) :: grv(3)  ! Gradient of local pseudopotential

    !Internal vars
    type(rad_func_t), pointer     :: kb_proj
    real(dp)                      :: rmod, dvdr

    v=0.0_dp
    grv=0.0_dp
    
    kb_proj => get_rad_func_p(species%kb_proj,i)
    rmod = sqrt(sum(r*r))
    rmod=rmod+tiny20
    if (rmod .gt. rad_cutoff(kb_proj)) return

    call rad_get(kb_proj,rmod,v,dvdr)
    
    grv(1:3) = dvdr * r(1:3)/rmod

  end subroutine get_value_of_kb_proj

  !------------------------------------------------------------------------

  subroutine get_rvalue_of_kb_proj(species,i,r,v,grv)
    type(species_info_t), intent(in) :: species
    integer, intent(in)   :: i       ! kb_proj index
    real(dp), intent(in)  :: r       ! Point vector, relative to atom
    real(dp), intent(out) :: v       ! Value of local pseudopotential
    real(dp), intent(out) :: grv     ! Gradient of local pseudopotential

    !Internal vars
    type(rad_func_t), pointer     :: kb_proj
    real(dp)                      :: rmod

    v=0.0_dp
    grv=0.0_dp
    
    kb_proj => get_rad_func_p(species%kb_proj,i)
    rmod = r+tiny20
    if (rmod .gt. rad_cutoff(kb_proj)) return

    call rad_get(kb_proj,rmod,v,grv)

  end subroutine get_rvalue_of_kb_proj

  !---------------------------------------------------------------------------

  subroutine set_orb(species,orb,iorb ) 
     type(species_info_t), intent(inout) :: species
     type(hilbert_vector_t),intent(in)    :: orb
     integer, intent(in)               :: iorb 

     call set_vector(species%orbs,orb,iorb)
  end subroutine set_orb
  
  !-------------------------------------------------------

  subroutine set_orbs_deg(species)
    type(species_info_t), intent(inout) :: species
    call set_deg(species%orbs)
  end subroutine set_orbs_deg
  
  !-------------------------------------------------------

   subroutine set_kb_projs_deg(species)
    type(species_info_t), intent(inout) :: species
    call set_deg(species%kb_proj)
  end subroutine set_kb_projs_deg
  
  !-------------------------------------------------------

  subroutine set_kb_proj(species,kb,ikb)
    type(species_info_t), intent(inout) :: species
    type(hilbert_vector_t), intent(in)  :: kb
    integer, intent(in)                 :: ikb
    call set_vector(species%kb_proj,kb,ikb)
  end subroutine set_kb_proj

  !---------------------------------------------------------

  subroutine filter_orbs(factor)
    !Filter orbitals of all the species.
    !The kc_max is found when the orbitals are generated.
    !(see module filter.f90, basis_gen.f90)
    real(dp), intent(in) :: factor !Instead of filtering the square of the orbs
                                   !we multiply them by this factor.

    integer :: is,norbs,io,l
    type(rad_func_t) :: filtered,non_filtered
    type(hilbert_vector_t) :: filtered_v,non_filtered_v

    write(6,'(a,f10.3,a)') "atm_types: Filtering orbitals. Kcutoff=",orbs_kc_max,' bohr^-1'
    do is=1,nspecies
       norbs = get_number_of_orbs_non_deg(species(is))
       do io=1,norbs
          non_filtered_v = get_vector(species(is)%orbs,io)
          non_filtered   = get_rad_func_v(non_filtered_v)
          l              = get_l_v(non_filtered_v)
          write(6,'(a,3i3)') "atm_types: species, orbital, l =", is,io,l

          filtered       = rad_filter(non_filtered,l,factor,2,orbs_kc_max)          
          call copy_vector(non_filtered_v,filtered_v)
          call set_rad_func_v(filtered_v,filtered)
          call set_orb(species(is),filtered_v,io)
          call destroy_vector(non_filtered_v)
          call destroy_vector(filtered_v)
          call rad_dealloc(non_filtered)
          call rad_dealloc(filtered)
       enddo
    enddo
  end subroutine filter_orbs

end module atm_types
