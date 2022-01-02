! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      module basis_specs
! 
! Alberto Garcia, August 2000, based on 'classic' redbasis.
! 
! Processes the basis information in an fdf file and populates
! the "basis specifications" data structures.
! 
! Here is a guide to the behavior of the main routine "read_basis_specs":
! 
! * Find out the number of species
! * Allocate storage for the data structures
! * Determine any "global" basis specification parameters:
!   - basis_size   (sz, dz, dzp, etc)
!   - basis_type   (split, nodes, etc) ("USER" is no longer valid)
! * Process the mandatory "Chemical-species-label" block to read the
!   species labels (not the same as chemical symbols) and atomic numbers.
!   (note that negative atomic numbers have a special meaning)
! * Assign default values to 'basis_size', 'basis_type'.
! * Assign default values to the 'mass' based on the atomic number.
! * Set up the information about the ground state of the atom.
! * Read information about the valence charge density in the
!   pseudopotential file and determine whether there are semicore
!   states. 
!   (A further semicore analysis is in 'autobasis')
! * Read the optional fdf blocks:
!   AtomicMass  (routine remass)
!   PAO.BasisSize - Overrides the default 'basis_size' on a per-species basis.
!   PAO.PolarizationScheme - Whether to use the standard perturbative approach
!                            for polarization orbitals, or to promote the shell             
!                            to a stand-alone status.
!   PAO.Basis - This is the most complex block, very flexible but in
!               need  of spelling-out the specific restrictions. Line by
!               line, these are:
! 
!   1st:   Species_label number_of_l_shells [basis_type] [ionic_charge] 
! 
!   For each l_shell:
!     [n= n] l nzeta [P [ nzeta_pol]] [E vcte rinn] [Q qcoe [qyuk [qwid]]] [F cutoff]
!   where 'n= n' is *mandatory* if the species has any semicore states,
!   and the 'P' (polarization), E (soft confinement potential), 
!   'Q' (charge confinement), and 'F' (filteret) sections are optional and can appear
!   in any order after nzeta. 
! 
!          rc_1  rc_2 .... rc_nzeta  
! 
!   are the cutoff radii in bohrs. This line is mandatory
!   If the number of rc's for a given shell is less than the number of
!   'zetas', the program will assign the last rc value to the remaining
!   zetas, rather than stopping with an error. This is particularly
!   useful for Bessel suites of orbitals.
!   A line containing contraction (scale) factors is optional, but if it
!   appears, the values *must be* real numbers. The same extension feature
!   as for the rc values works here.      
!   --------------------------------------------------------------------
! 
!   After processing PAO.Basis (if it exists), whatever PAO information
!   which is not already available is determined in routine 'autobasis', 
!   using the following defaults:
! 
!   (There is a first check to make sure that species
!   with Bessel floating orbitals  have been included in 
!   the PAO.Basis block)
!   
!   The maximum angular momentum for the basis (lmxo) (excluding any
!   polarization orbitals) is set to that of the ground state, as returned by
!   routine 'lmxofz' from Z. If there are any semicore states with a
!   higher l, lmxo is set to the maximum l needed.
! 
!   Each l-shell contains just one shell of PAOS, with n set
!   to the appropriate ground state value, except in the case
!   of semicore states, which add further lower-lying shells.      
! 
!   Nzeta is determined by the first two characters of the 'basis_size'
!   string: 1 for 'sz' and  2 for 'dz'. Lower-lying semicore states
!   get nzeta=1. Nzeta will be zero (for the top-most shell) if the
!   ground-state valence shell for this l is empty.      
!      
!   There are no 'per l-shell' polarization orbitals, except if the
!   third character of 'basis_size' is 'p' (as in 'dzp'), in which
!   case polarization orbitals are defined so they have the minimum      
!   angular momentum l such that there are no occupied orbitals 
!   with the same l in the valence shell of the ground-state 
!   atomic configuration. They polarize the corresponding l-1 shell.
!   (See above for generation scheme)
!
!   The Soft-Confinement parameters 'rinn' and 'vcte' are set to 0.0
!   The Charge-Confinement parameters 'qcoe', 'qyuk' and 'qwid' 
!   are set to 0.0, 0.0 and 0.01
! 
!   rc(1:nzeta) is set to 0.0
!   lambda(1:nzeta) is set to 1.0  (this is a change from old practice)
!
!  ----------------------------------
!  
!   Next come the blocks associated to the KB projectors:
! 
!   Block PS.lmax (routine relmxkb) assigns a per-species maximum l for
!   the KB projectors. Internally, it sets the 'lmxb_requested'
!   component of the data structure, which by default is -1.
!  
!   Block PS.KBprojectors (optional) sets the number of KB projectors
!   per angular momentum. The same routine which reads this block
!   (readkb) also sets any appropriate defaults:
! 
! 
!   If the species does not appear in the  PS.KBprojectors block, then
!       Setting of lmxkb:
!       If lmxkb is set in PS.lmax, 
!          set it to that
!       else
!          Use PAO information: (routine set_default_lmxkb)
!          (Set it to a meaningless value for floating orbitals)
!           Set it to lmxo+1, or to lpol+1, where lpol is the angular
!           momentum of the highest-l polarization orbital.
!       endif
!
!       There is a further check for lmxkb: If the pseudopotential
!       contains semilocal components up to lmax_pseudo, lmxkb is
!       set to this number.
!       (Note that this is appropriate for the new behavior of using
!       just Vlocal for any channels for which there is no V_l(r).)      
      
!       The  number of KB projectors per l is set to the number of
!       n-shells in the corresponding PAO shell with the same l. For l
!       greater than lmxo, it is set to 1. The reference energies are
!       in all cases set to huge(1.d0) to mark the default.
!
!       Archaeological note: The first implementation of the basis-set generation
!       module had only non-polarization orbitals. Polarization orbitals were added
!       later as "second-class" companions. This showed in details like "lmxo" (the
!       maximum l of the basis set) not taking into account polarization orbitals.
!       Polarization orbitals were tagged at the end, without maintaining l-shell
!       ordering.
!       Even later, support for "semicore" orbitals was added. A new ('nsm') index was
!       used to distinguish the different orbitals in a l-shell. Polarization orbitals
!       were not brought into this classification.
!       In this version, 'lmxo' takes into account any polarization orbitals.      
! 
!       Future work should probably remove the separate treatment of polarization orbitals.
!       (Note that if *all* orbitals are specified in a PAO.Basis block, in effect turning
!       perturbative polarization orbitals into 'normal' orbitals, this problem is not present.)
!     
! =======================================================================
!
      use precision
      use basis_types, only: basis_def_t, shell_t, lshell_t, kbshell_t
      use basis_types, only: nsp, basis_parameters, ground_state_t
      use basis_types, only: destroy, copy_shell, initialize
      use m_ncps,      only: pseudo_init_constant
      use periodic_table, only: qvlofz, lmxofz, cnfig, atmass
      use chemical
      use sys
      use fdf

      Implicit None

      type(basis_def_t), pointer :: basp => null()
      type(shell_t), pointer :: s => null()
      type(shell_t), pointer :: polarized_shell_ptr
      type(lshell_t), pointer :: ls => null()
      type(kbshell_t), pointer :: k => null()

      character(len=1), parameter   ::
     $                           sym(0:4) = (/ 's','p','d','f','g' /)

!     Default Soft-confinement parameters set by the user
!
      logical, save   :: lsoft
      real(dp), save  :: softRc, softPt

! Default norm-percentage for the automatic definition of
! multiple-zeta orbitals with the 'SPLIT' option

      real(dp), parameter       :: splnorm_default=0.15_dp
      real(dp), parameter       :: splnormH_default=-1.0_dp
      real(dp), save            :: global_splnorm, global_splnorm_H
      integer        isp  ! just an index dummy variable for the whole module

!
      logical, save, public     :: restricted_grid
      logical, parameter        :: restricted_grid_default = .true.

      real(dp), save, public    :: rmax_radial_grid
      
      public :: read_basis_specs
      public :: label2species

      private

      CONTAINS

!---
      subroutine read_basis_specs()

      use m_spin, only: SpOrb
      use m_ncps, only: pseudo_read
      use m_spin_orbit_potentials, only: valid_spin_orbit_potentials

      character(len=15), parameter  :: basis_size_default='standard'
      character(len=10), parameter  :: basistype_default='split'

      character(len=15) :: basis_size
      character(len=10) :: basistype_generic

      type(block_fdf)            :: bfdf
      type(parsed_line), pointer :: pline

      type(ground_state_t), pointer :: gs

      integer nns, noccs, i, ns_read, l
      logical synthetic_atoms, found, reparametrize_pseudos
      real(dp) :: new_a, new_b, new_rmax

!------------------------------------------------------------------------
      reparametrize_pseudos =
     $   fdf_boolean('ReparametrizePseudos',.false.)
!
!        r(i) = b * ( exp(a*(i-1)) - 1 )
!        spacing near zero: delta = a*b
!        If the desired spacing at r=R is delta*(1+beta), 
!        then:
!               a = beta*delta/R
!               b = R/beta
!               N at R = 1 + R*ln(1+beta)/(beta*delta)
!
      if (reparametrize_pseudos) then
         ! These default values will provide a grid spacing
         ! of 1.0e-5 near r=0, and 0.01 near r=10 (bohr units)
         ! with N at Rmax=100 on the order of 10000 points.
         new_a = fdf_double("NewAParameter",0.001_dp)
         new_b = fdf_double("NewBParameter",0.01_dp)
         new_rmax = fdf_double("NewGridRmax",0.0_dp)
      endif

!
!     Whether or not to restrict integration grids to an odd number
!     of points (and hence to shift rcs accordingly)     
!
      restricted_grid =
     $     fdf_boolean("Restricted.Radial.Grid",
     $                  restricted_grid_default)
!
!     If non-zero, the value will be the maximum value of the
!     radial coordinate
!
      if (reparametrize_pseudos) then
         rmax_radial_grid = fdf_double('Rmax.Radial.Grid',50.0_dp)
      else
         rmax_radial_grid = fdf_double('Rmax.Radial.Grid',0.0_dp)
      endif

!
      basis_size = fdf_string('PAO.BasisSize',basis_size_default)
      call size_name(basis_size)
      basistype_generic = fdf_string('PAO.BasisType',basistype_default)
      call type_name(basistype_generic)

C Read information about defaults for soft confinement

      lsoft  = fdf_boolean('PAO.SoftDefault',.false.)
      softRc = fdf_double('PAO.SoftInnerRadius',0.9d0)
      softPt = fdf_double('PAO.SoftPotential',40.0d0)

C Sanity checks on values

      softRc = max(softRc,0.00d0)
      softRc = min(softRc,0.99d0)
      softPt = abs(softPt)
!
!     Read defaults for split_norm parameter

      global_splnorm = fdf_double('PAO.SplitNorm',splnorm_default)
      global_splnorm_H = fdf_double('PAO.SplitNormH',splnormH_default)
      if (global_splnorm_H < 0.0_dp) global_splnorm_H = global_splnorm

!------------------------------------------------------------------
!
!     Use standard routine in chemical module to process the
!     chemical species
!
      nsp = size(basis_parameters)

      synthetic_atoms = .false.

      do isp=1,nsp
        basp=>basis_parameters(isp)

        basp%label = species_label(isp)
        basp%z = atomic_number(isp)
        basp%floating = is_floating(isp)
        basp%bessel = is_bessel(isp)
        basp%synthetic = is_synthetic(isp)

        basp%basis_size = basis_size
        basp%basis_type = basistype_generic
        if (basp%floating) then
          basp%mass = 1.d40   ! big but not too big, as it is used
                              ! later in computations
        else if (basp%synthetic) then
          basp%mass = -1.0_dp      ! Signal -- Set later
        else
          basp%mass = atmass(abs(int(basp%z)))
        endif

        write(6,"(/,a)") " ---- Processing specs for species: " //
     $                    trim(basp%label)
        if (basp%bessel) then
          ! Initialize a constant pseudo
          call pseudo_init_constant(basp%pseudopotential)
        else if (basp%synthetic) then
          synthetic_atoms = .true.
          ! Will set gs later
          call pseudo_read(basp%label,basp%pseudopotential,
     $         basp%psml_handle,basp%has_psml_ps,
     $         new_grid=reparametrize_pseudos,a=new_a,b=new_b,
     $         rmax=new_rmax)
        else
          call ground_state(abs(int(basp%z)),basp%ground_state)
          call pseudo_read(basp%label,basp%pseudopotential,
     $         basp%psml_handle,basp%has_psml_ps,
     $         new_grid=reparametrize_pseudos,a=new_a,b=new_b,
     $         rmax=new_rmax)
        endif
!        if (reparametrize_pseudos.and. .not. basp%bessel)
!     .    call pseudo_reparametrize(p=basp%pseudopotential,
!     .                             a=new_a, b=new_b,label=basp%label)

             if (SpOrb .and. basp%has_psml_ps) then
                if (.not. valid_spin_orbit_potentials(basp%psml_handle))
     $                    then
                   call die(
     $            "Cannot do spin-orbit without proper semilocal pots")
                endif
             endif

      enddo

      ! Allow manual specification of valence configuration,
      ! even for non-synthetic atoms, with the same block

      found = fdf_block('SyntheticAtoms',bfdf)
      if (.not. found) then
         if (synthetic_atoms)
     .          call die("Block SyntheticAtoms does not exist.")
      else
        ns_read = 0
        do while(fdf_bline(bfdf, pline))

          ns_read = ns_read + 1
          if (.not. fdf_bmatch(pline,'i'))
     .      call die("Wrong format in SyntheticAtoms")
          isp = fdf_bintegers(pline,1)
          if (isp .gt. nsp .or. isp .lt. 1)
     .      call die("Wrong specnum in SyntheticAtoms")
          basp => basis_parameters(isp)
          gs   => basp%ground_state
          if (.not. fdf_bline(bfdf, pline)) call die("No n info")
          nns = fdf_bnintegers(pline)
          if (nns .lt. 4)
     .      call die("Please give all valence n's " //
     .               "(up to l=3) " //
     .               "in SyntheticAtoms block")
          gs%n = 0
          do i = 1, nns
            gs%n(i-1) = fdf_bintegers(pline,i)
          enddo
          if (.not. fdf_bline(bfdf, pline))
     .      call die("No occupation info in Synthetic block")
          noccs = fdf_bnvalues(pline)
          if (noccs .lt. nns) call die("Need more occupations")
          gs%occupation(:) = 0.0_dp
          do i = 1, noccs
            gs%occupation(i-1) = fdf_bvalues(pline,i)
          enddo
          ! Determine the last occupied state in the atom
          do i = nns, 1, -1
            if (gs%occupation(i-1) /=0) then
              gs%lmax_valence = i-1
              exit
            endif
          enddo
          gs%occupied(0:4) = (gs%occupation .gt. 0.0_dp)
          gs%occupied(4) = .false.
          gs%z_valence = sum(gs%occupation(0:noccs-1))
          write(6,'(a,i2)',advance='no')
     .         'Ground state valence configuration (synthetic): '
          do l=0,3
            if (gs%occupied(l))
     .        write(6,'(2x,i1,a1,f8.5)',advance='no')
     .          gs%n(l),sym(l), gs%occupation(l)
          enddo
          write(6,'(a)') ''

        enddo
        write(6,"(a,i2)") "Number of synthetic species: ", ns_read
      endif
!
!  Defer this here in case there are synthetic atoms
      do isp=1,nsp
        basp=>basis_parameters(isp)
        if (basp%synthetic) then
          gs => basp%ground_state
          if (gs%z_valence .lt. 0.001)
     .      call die("Synthetic species not detailed")
        endif
        call semicore_check(isp)
      enddo

      call remass()
      call resizes()
      call read_polarization_scheme()
      call repaobasis()
      call autobasis()
      call relmxkb()
      call readkb()

!      do isp=1,nsp
!        call print_basis_def(basis_parameters(isp))
!      enddo

      end subroutine read_basis_specs
!-----------------------------------------------------------------------

      function label2species(str) result(index)
      integer index
      character(len=*), intent(in) ::  str

      integer i

      index = 0
      do i=1,nsp
        if(leqi(basis_parameters(i)%label,str)) index = i
      enddo
      end function label2species
!-----------------------------------------------------------------------

      subroutine ground_state(z,gs)
      integer, intent(in)               ::  z
      type(ground_state_t), intent(out) :: gs
!
!     Determines ground state valence configuration from Z
!
      integer l, latm

      gs%z_valence = 0.d0
      do l=0,4
        gs%occupation(l)=0.0d0
      enddo

      call lmxofz(z,gs%lmax_valence,latm)
      call qvlofz(z,gs%occupation(:))
      do l=0,gs%lmax_valence
        gs%z_valence = gs%z_valence + gs%occupation(l)
      enddo
      call cnfig(z,gs%n(0:4))

      write(6,'(a,i2)',advance='no')
     .     'Ground state valence configuration: '
      gs%occupied(4) = .false.         !! always
      do l=0,3
        gs%occupied(l) =  (gs%occupation(l).gt.0.1d0)
        if (gs%occupied(l))
     .    write(6,'(2x,i1,a1,i2.2)',advance='no')
     .       gs%n(l),sym(l),nint(gs%occupation(l))
      enddo
      write(6,'(a)') ''

      end subroutine ground_state

!---------------------------------------------------------------
      subroutine readkb()
      integer lpol, isp, ish, i, l
      integer :: lmax_pseudo
      character(len=20) unitstr

      type(block_fdf)            :: bfdf
      type(parsed_line), pointer :: pline

      lpol = 0

      if (fdf_block('PS.KBprojectors',bfdf) ) then
 
! First pass to find out about lmxkb and set any defaults.

        do while(fdf_bline(bfdf, pline))    !! over species
          if (.not. fdf_bmatch(pline,'ni'))
     .      call die("Wrong format in PS.KBprojectors")
          isp = label2species(fdf_bnames(pline,1))
          if (isp .eq. 0) then
            write(6,'(a,1x,a)')
     .        "WRONG species symbol in PS.KBprojectors:",
     .        trim(fdf_bnames(pline,1))
            call die()
          endif
          basp => basis_parameters(isp)
          basp%nkbshells = fdf_bintegers(pline,1)
          do ish= 1, basp%nkbshells
            if (.not. fdf_bline(bfdf, pline)) call die("No l nkbl")
            if (.not. fdf_bmatch(pline,'ii'))
     .        call die("Wrong format l nkbl")
            l = fdf_bintegers(pline,1)
            ! Check that we have enough V_ls in the pseudo to
            ! generate this l...
            ! if (l > (basp%pseudopotential%npotd-1)) then
            ! ... will do it later, once lmxkb has been
            ! determined at the end of the cascade of possible
            ! inputs.
            if (l .gt. basp%lmxkb) basp%lmxkb = l
            if (.not. fdf_bline(bfdf, pline)) then
               if (ish .ne. basp%nkbshells)
     .          call die("Not enough kb shells for this species...")
              ! There is no line with ref energies
            else  if (fdf_bmatch(pline,'ni')) then
                ! We are seeing the next species' section
               if (ish .ne. basp%nkbshells)
     .            call die("Not enough shells for this species...")
               if (.not. fdf_bbackspace(bfdf))
     .               call die('readkb: ERROR in PS.KBprojectors block')
            else if (fdf_bmatch(pline,'ii')) then
                ! We are seeing the next shell's section
               if (ish .gt. basp%nkbshells)
     .            call die("Too many kb shells for this species...")
               if (.not. fdf_bbackspace(bfdf))
     .            call die('readkb: ERROR in PS.KBprojectors block')
            endif
          enddo       ! end of loop over shells for species isp

        enddo
      endif

      do isp=1,nsp
!
!      Fix defaults and allocate kbshells
!
        basp=>basis_parameters(isp)
        if (basp%lmxkb .eq. -1) then ! not set in KBprojectors 
          if (basp%lmxkb_requested.eq.-1) then ! not set in PS.lmax
            basp%lmxkb = set_default_lmxkb(isp) ! Use PAO info
          else
            basp%lmxkb = basp%lmxkb_requested
          endif
          allocate(basp%kbshell(0:basp%lmxkb))
          do l=0,basp%lmxkb
            call initialize(basp%kbshell(l))
            k=>basp%kbshell(l)
            k%l = l
            if (l.gt.basp%lmxo) then
              k%nkbl = 1
            else
              ! Set equal to the number of PAO shells with this l
              k%nkbl = basp%lshell(l)%nn     
              ! Should include polarization orbs (as in Ti case: 3p..4p*)
              ! (See 'archaeological note' in the header of this file)
              ! ... but if the element was not in the PAO.Basis block, any such polarization orbital
              !     has been included in 'nn' already.
              if (l>0 .and. basp%in_pao_basis_block) then
                 do i = 1, basp%lshell(l-1)%nn
                    if (basp%lshell(l-1)%shell(i)%polarized) then
                       k%nkbl = k%nkbl + 1
                       write(6,"(a,i1,a)") trim(basp%label) //
     $                  ': nkbl increased for l=',l,
     $                  ' due to the presence of a polarization orbital'
                    endif
                 enddo
              endif
              if (k%nkbl.eq.0) then
                write(6,*) 'Warning: Empty PAO shell. l =', l
                write(6,*) 'Will have a KB projector anyway...'
                k%nkbl = 1
              endif
            endif
            allocate(k%erefkb(1:k%nkbl))
            k%erefkb(1:k%nkbl) = huge(1.d0)
          enddo
        else          ! Set in KBprojectors
          if (basp%lmxkb_requested.ne.-1) then ! set in PS.lmax
            if (basp%lmxkb.ne.basp%lmxkb_requested) then
              call die("LmaxKB conflict between " //
     $                 "PS.Lmax and PS.KBprojectors blocks")
            endif
          endif
          !! OK, we have a genuine lmxkb
          allocate(basp%kbshell(0:basp%lmxkb))
          do l=0,basp%lmxkb
            call initialize(basp%kbshell(l))
          enddo
        endif
        if (basp%z .le. 0) then
          if (basp%lmxkb .ne. -1)
     .      call die("Floating orbs cannot have KB projectors...")
        endif
      enddo
!
!     Now re-scan the block (if it exists) and fill in as instructed
!            
      if (fdf_block('PS.KBprojectors',bfdf) ) then

        do while(fdf_bline(bfdf, pline))     !! over species
          if (.not. fdf_bmatch(pline,'ni'))
     .      call die("Wrong format in PS.KBprojectors")
          isp = label2species(fdf_bnames(pline,1))
          if (isp .eq. 0) then
            write(6,'(a,1x,a)')
     .        "WRONG species symbol in PS.KBprojectors:",
     .        trim(fdf_bnames(pline,1))
            call die()
          endif
          basp => basis_parameters(isp)
          basp%nkbshells = fdf_bintegers(pline,1)
          do ish=1, basp%nkbshells
            if (.not. fdf_bline(bfdf,pline)) call die("No l nkbl")
            if (.not. fdf_bmatch(pline,'ii'))
     .        call die("Wrong format l nkbl")
            l = fdf_bintegers(pline,1)
            k => basp%kbshell(l)
            k%l = l
            k%nkbl = fdf_bintegers(pline,2)
            if (k%nkbl < 0) then
               call die("nkbl < 0 in PS.KBprojectors")
            endif
            if (k%nkbl == 0) then
               call message("WARNING","nkbl=0 in PS.KBprojectors")
            endif
            allocate(k%erefkb(k%nkbl))
            if (.not. fdf_bline(bfdf,pline)) then
              if (ish .ne. basp%nkbshells)
     .          call die("Not enough shells for this species...")
              ! There is no line with ref energies
              ! Use default values
              k%erefKB(1:k%nkbl) = huge(1.d0)
            else if (fdf_bmatch(pline,'ni')) then
                ! We are seeing the next species' section
               if (ish .ne. basp%nkbshells)
     .            call die("Not enough kb shells for this species...")
                ! Use default values for ref energies
               k%erefKB(1:k%nkbl) = huge(1.d0)
               if (.not. fdf_bbackspace(bfdf))
     .            call die('readkb: ERROR in PS.KBprojectors block')
            else if (fdf_bmatch(pline,'ii')) then
                ! We are seeing the next shell's section
                if (ish .gt. basp%nkbshells)
     .            call die("Too many kb shells for this species...")
                ! Use default values for ref energies
                k%erefKB(1:k%nkbl) = huge(1.d0)
                if (.not. fdf_bbackspace(bfdf))
     .            call die('readkb: ERROR in PS.KBprojectors block')
             else
                if (fdf_bnreals(pline) .ne. k%nkbl)
     .            call die("Wrong number of energies")
                unitstr = 'Ry'
                if (fdf_bnnames(pline) .eq. 1)
     .            unitstr = fdf_bnames(pline,1)
                ! Insert ref energies in erefkb
                do i= 1, k%nkbl
                  k%erefKB(i) =
     .                 fdf_breals(pline,i)*fdf_convfac(unitstr,'Ry')
                enddo
             endif
          enddo            ! end of loop over shells for species isp

          ! For those l's not specified in block, use default values
          do l=0, basp%lmxkb
            k => basp%kbshell(l)
            if (k%l.eq.-1) then
              k%l = l
              k%nkbl = 1
              allocate(k%erefkb(1))
              k%erefkb(1) = huge(1.d0) / 4.d0   ! signal different phase
            endif
          enddo
        enddo   !! Over species
      endif

      do isp=1,nsp
!
!      Check that we have enough semilocal components...
!
         basp=>basis_parameters(isp)
         lmax_pseudo = basp%pseudopotential%npotd - 1 
         if (basp%lmxkb > lmax_pseudo) then
            write(6,'(a,i1,a)')
     .           trim(basp%label) //
     .           " pseudopotential only contains V_ls up to l=",
     .           lmax_pseudo, " -- lmxkb reset."
            basp%lmxkb = lmax_pseudo
         endif
      enddo

      do isp=1,nsp
!
!      Check that we have enough semilocal components...
!
         basp=>basis_parameters(isp)
         lmax_pseudo = basp%pseudopotential%npotd - 1 
         if (basp%lmxkb > lmax_pseudo) then
            write(6,'(a,i1,a)')
     .           "Pseudopotential only contains V_ls up to l=",
     .           lmax_pseudo, " -- lmxkb reset."
            basp%lmxkb = lmax_pseudo
         endif
      enddo

      end subroutine readkb
!---------------------------------------------------------------

      subroutine repaobasis()
      
      use m_semicore_info_froyen, only: get_n_semicore_shells

      integer isp, ish, nn, i, ind, l, indexp, index_splnorm
      integer nrcs_zetas
      integer :: nsemic_shells(0:3)
      logical :: will_have_polarization_orb 
      
      type(block_fdf)            :: bfdf
      type(parsed_line), pointer :: pline

      if (.not. fdf_block('PAO.Basis',bfdf)) RETURN

      do while(fdf_bline(bfdf, pline))     !! over species
        if (.not. fdf_bmatch(pline,'ni'))
     .    call die("Wrong format in PAO.Basis")
        isp = label2species(fdf_bnames(pline,1))
        if (isp .eq. 0) then
          write(6,'(a,1x,a)')
     .      "WRONG species symbol in PAO.Basis:",
     .      trim(fdf_bnames(pline,1))
          call die()
        endif

        basp => basis_parameters(isp)
        basp%in_pao_basis_block = .true.
        basp%label = fdf_bnames(pline,1)
        basp%nshells_tmp = fdf_bintegers(pline,1)
        basp%lmxo = 0

        if (.not. (basp%synthetic .or.
     $             basp%bessel)) then
          ! To report on pseudized shells and, in the future,
          ! check on the specified structure of the block
          call get_n_semicore_shells(basp%pseudopotential,nsemic_shells)
        endif
        
        !! Check whether there are optional type and ionic charge
        if (fdf_bnnames(pline) .eq. 2)
     .    basp%basis_type = fdf_bnames(pline,2)
        if (fdf_bnvalues(pline) .eq. 2)
     .    basp%ionic_charge = fdf_bvalues(pline,2)
        allocate(basp%tmp_shell(basp%nshells_tmp))

        ! These are (l,n) shells
        shells: do ish= 1, basp%nshells_tmp
          s => basp%tmp_shell(ish)
          call initialize(s)
          if (.not. fdf_bline(bfdf,pline)) call die("No l nzeta, etc")

          if (fdf_bmatch(pline,'niii')) then
            s%n = fdf_bintegers(pline,1)
            s%l = fdf_bintegers(pline,2)
            basp%lmxo = max(basp%lmxo,s%l)
            s%nzeta = fdf_bintegers(pline,3)
          elseif (fdf_bmatch(pline,'ii')) then
            !    l, nzeta

            if (basp%semic)
     .        call die("Please specify n if there are semicore states")

            s%l = fdf_bintegers(pline,1)
            s%n = basp%ground_state%n(s%l)
            s%nzeta = fdf_bintegers(pline,2)
            basp%lmxo = max(basp%lmxo,s%l)
          else
            call die("Bad format of (n), l, nzeta line in PAO.Basis")
          endif
!
! If this is a filteret basis then the number of zetas input must be one
!
          if (basp%basis_type.eq.'filteret') then
            s%nzeta = 1
          endif
!
! Optional stuff: Polarization, Soft-confinement Potential and Filteret Cutoff
!
! Split norm
!
          if (fdf_bsearch(pline,'S',index_splnorm)) then
            if (fdf_bmatch(pline,'v',after=index_splnorm)) then
              s%split_norm = fdf_bvalues(pline, ind=1,
     .                                   after=index_splnorm)
              if (s%split_norm .eq. 0.0_dp)
     .          write(6,"(a)")
     .            "WARNING: zero split_norm after S in PAO.Basis"
              s%split_norm_specified = .TRUE.
            else
              call die("Specify split_norm after S in PAO.Basis")
            endif
          else
            if (abs(basp%z) .eq. 1) then
              s%split_norm = global_splnorm_H
            else
              s%split_norm = global_splnorm
            endif
          endif
!
! Polarization functions
!
          if (fdf_bsearch(pline,'P',indexp)) then
            s%polarized = .TRUE.
            if (fdf_bmatch(pline,'i',after=indexp)) then
              s%nzeta_pol = fdf_bintegers(pline,ind=1,after=indexp)
            else
              s%nzeta_pol = 1
            endif
            basp%lmxo = max(basp%lmxo,s%l+1)  ! NOTE new behavior
          endif
!
! Soft-confinement
!
          if (fdf_bsearch(pline,'E',indexp)) then
            if (fdf_bmatch(pline,'vv',after=indexp)) then
              s%vcte = fdf_bvalues(pline,ind=1,after=indexp)
              s%rinn = fdf_bvalues(pline,ind=2,after=indexp)
            else
              call die("Need vcte and rinn after E in PAO.Basis")
            endif
          elseif (lsoft) then
            s%vcte = softPt 
            s%rinn = -softRc
          else
            s%vcte = 0.0_dp
            s%rinn = 0.0_dp
          endif
!
! Charge confinement
! 
          if (fdf_bsearch(pline,'Q',indexp)) then
            if (fdf_bmatch(pline,'vvv',after=indexp)) then
               s%qcoe = fdf_bvalues(pline,ind=1,after=indexp)
               s%qyuk = fdf_bvalues(pline,ind=2,after=indexp)
               s%qwid = fdf_bvalues(pline,ind=3,after=indexp)
            elseif (fdf_bmatch(pline,'vv',after=indexp)) then
               s%qcoe = fdf_bvalues(pline,ind=1,after=indexp)
               s%qyuk = fdf_bvalues(pline,ind=2,after=indexp)
               s%qwid = 0.01_dp
            elseif (fdf_bmatch(pline,'v',after=indexp)) then
               s%qcoe = fdf_bvalues(pline,ind=1,after=indexp)
               s%qyuk = 0.0_dp
               s%qwid = 0.01_dp
            else
               call die("Need one, two or three real numbers after Q in 
     .                   PAO.Basis")
            endif
          else
            s%qcoe = 0.0_dp
            s%qyuk = 0.0_dp
            s%qwid = 0.01_dp
          endif
!
! Filteret cutoff
!
          if (fdf_bsearch(pline,"F",indexp)) then
            if (fdf_bmatch(pline,"v",after=indexp)) then
              s%filtercut = fdf_bvalues(pline,ind=1,after=indexp)
            else
              call die("Need cut-off after F in PAO.Basis")
            endif
          else
            s%filtercut = 0.0_dp
          endif

          allocate(s%rc(s%nzeta),s%lambda(s%nzeta))
          s%rc(:) = 0.d0
          s%lambda(:) = 1.d0
          if (.not. fdf_bline(bfdf,pline)) call die("No rc's")
          
          ! Use the last rc entered for the successive zetas
          ! if there are not enough values (useful for Bessel)
          nrcs_zetas = fdf_bnvalues(pline)
          if (nrcs_zetas < 1) then
           call die("Need at least one rc per shell in PAO.Basis block")
          endif
          do i= 1, s%nzeta
             if (i <= nrcs_zetas) then
                s%rc(i) = fdf_bvalues(pline,i)
             else
                s%rc(i) = s%rc(nrcs_zetas)
             endif
          enddo

          if (s%split_norm_specified) then
            do i = 2,s%nzeta
              if (s%rc(i) /= 0.0_dp) then
                write(6,"(/,a,i1,a,f8.4,/)")
     .            "*Warning: Per-shell split_norm parameter " //
     .            "will not apply to zeta-", i, ". rc=", s%rc(i)
              endif
            enddo
          endif

          ! Optional scale factors. They MUST be reals, or else...
          if (.not. fdf_bline(bfdf,pline)) then
            if (ish .ne. basp%nshells_tmp)
     .        call die("Not enough shells")
              ! Default values for scale factors
          else
            if (.not. fdf_bmatch(pline,'r')) then
              ! New shell or species
              ! Default values for the scale factors
              if (.not. fdf_bbackspace(bfdf)) 
     .          call die('repaobasis: ERROR in PAO.Basis block')
              cycle shells
            else
              ! Read scale factors
              ! Use the last scale factor entered for the successive zetas
              ! if there are not enough values 
               nrcs_zetas = fdf_bnreals(pline)
               if (nrcs_zetas < 1) then
                 call die("Need at least one scale factor in PAO.Basis")
               endif
               do i= 1, s%nzeta
                  if (i <= nrcs_zetas) then
                     s%lambda(i) = fdf_breals(pline,i)
                  else
                     s%lambda(i) = s%lambda(nrcs_zetas)
                  endif
               enddo
            endif
          endif

        enddo shells
        ! Clean up for this species
      enddo
!
!        OK, now classify the (l,n) states by l-shells
!
      do isp = 1, nsp
        basp => basis_parameters(isp)
        if (.not. basp%in_pao_basis_block) CYCLE

        polarized_shell_ptr => null()
        allocate (basp%lshell(0:basp%lmxo))

        loop_l: do l= 0, basp%lmxo
          ls => basp%lshell(l)
          call initialize(ls)
          ls%l = l
          will_have_polarization_orb = .false.
          
          ! Search for tmp_shells with given l
          nn = 0
          do ish= 1, basp%nshells_tmp
            s => basp%tmp_shell(ish)
            if (s%l .eq. l) nn=nn+1

            if (s%polarized .and. (s%l == (l-1))) then
               will_have_polarization_orb = .true.
            endif
            
          enddo

          if (will_have_polarization_orb) then
             ! We still give the option of treating the polarization
             ! orbital requested through a "P" option in the PAO.Basis block
             ! with a non-perturbative scheme. The main choices are handled
             ! by code in 'read_polarization_scheme', and result in the 
             ! setting of the variable 'non_perturbative_polorbs' in the
             ! specification derived type.

             ! Plus, we check whether there are already other (semicore) shells
             ! with the same angular momentum as the polarization orbital.
             ! In that case, the pol orb will have a node. The algorithm
             ! to generate it is not robust enough, so, if the fallback has
             ! been enabled (it is not by default), we switch to generating
             ! the orbital explicitly (unless forced by the user in the
             ! PAO.Polarization.Scheme block or with a global setting)
             if (nn >= 1) then
                if (basp%non_pert_polorbs_fallback) then
                   if (.not. basp%force_perturbative_polorbs) then
                      basp%polorb_with_semicore = .true.
                      basp%non_perturbative_polorbs = .true.
                      write(6,'(a,1x,a)')
     .                     "Fallback to non-perturbative " //
     $                     "polarization scheme for",
     $                     trim(basp%label)
                   endif
                endif
             endif

             ! If we have decided on the non-perturbative approach,
             ! make room for a new shell
             if (basp%non_perturbative_polorbs) then
                nn=nn+1
             endif
          endif
                
          ls%nn = nn
          if (nn.eq.0) then
            !! What else do we do here?
            cycle loop_l  
          endif
                                     !! 
          allocate(ls%shell(1:nn))
          ! Collect previously allocated shells
          ind = 0
          do ish=1, basp%nshells_tmp
            s => basp%tmp_shell(ish)

            if (s%l .eq. l) then
              ind = ind+1
              call copy_shell(source=s,target=ls%shell(ind))
              if (basp%non_perturbative_polorbs) then
                 if (s%polarized) then
                    ! Save for later. There should be just one 'polarized' shell
                    polarized_shell_ptr => ls%shell(ind)
                    ! Remove markers (in new objects)
                    ls%shell(ind)%polarized = .false.
                    ls%shell(ind)%was_polarized = .true.
                    ls%shell(ind)%nzeta_pol = 0
                 endif
              endif
            endif

            if (basp%non_perturbative_polorbs) then
               if (s%polarized .and. (s%l == (l-1))) then
                  ! Note that we have already seen this (parent) shell
                  ! in the previous iteration of loop_l
                  ind = ind + 1
                  ! Copy again to inherit data
                  call copy_shell(source=s,target=ls%shell(ind))
                  ls%shell(ind)%n = basp%ground_state%n(l)     !! earlier: -1   ! reset to find later
                  ls%shell(ind)%l = l
                  ls%shell(ind)%nzeta = s%nzeta_pol
                  ls%shell(ind)%nzeta_pol = 0
                  ls%shell(ind)%polarized = .false.
                  ls%shell(ind)%polarization_shell = .true.
                  ls%shell(ind)%shell_being_polarized =>
     $                 polarized_shell_ptr
                  polarized_shell_ptr => null()
                  ! rc's, lambdas, etc, will remain as in s.

                  ! ... except that we honor any explicit charge-confinement settings
                  ! for this species in the PAO.Polarization.Scheme block or the
                  ! global PAO.Polarization.Charge.Confinement option.
                  ! ... and maybe also the rc expansion factor
                  if (basp%qconf_options_polorbs_set) then
                     ls%shell(ind)%qcoe =
     $                    basp%qconf_options_polorbs%qcoe
                     ls%shell(ind)%qyuk =
     $                    basp%qconf_options_polorbs%qyuk
                     ls%shell(ind)%qwid =
     $                    basp%qconf_options_polorbs%qwid
                  else
                     ! (would apply the same to all species)
                     ls%shell(ind)%qcoe = fdf_get(
     $                    "PAO.Polarization.Charge.Confinement",
     $                    0.0_dp)
                  endif

               endif
            endif

         enddo

!!##         if (nn.eq.1) then
            ! If n was not specified, set it to ground state n
!!            if (ls%shell(1)%n.eq.-1)
!!     .        ls%shell(1)%n=basp%ground_state%n(l)
!!##          endif
          !! Do we have to sort by n value????
          !!
        enddo loop_l
        !! Destroy temporary shells in basp
        !! Warning: This does seem to destroy information!!
        call destroy(basp%tmp_shell)
      enddo

      end subroutine repaobasis
!_______________________________________________________________________


         function set_default_lmxkb(is) result(lmxkb)
         integer lmxkb
         integer, intent(in) :: is

         lmxkb = -1
         basp=>basis_parameters(is)
         if (basp%z .le. 0) return     ! Leave it at -1 for floating orbs.

         lmxkb = basp%lmxo + 1    ! This now includes any polarization orbitals

         write(6,'(3a,i1,/,2a,/,a)') 'For ', trim(basp%label),
     .              ', standard SIESTA heuristics set lmxkb to ',
     .              lmxkb,
     .              ' (one more than the basis l,',
     .              ' including polarization orbitals).',
     .  'Use PS.lmax or PS.KBprojectors blocks to override.'
!
!        But there is an upper limit for sanity: f is the highest
!
         if (lmxkb.gt.3) then
            write(6,'(3a,i1)') 'Warning: For ', trim(basp%label),
     .           ' lmxkb would have been set to ', lmxkb
            write(6,'(a)')
     .           'Setting it to maximum value of 3 (f projector)'
            lmxkb = 3
         endif

         end function set_default_lmxkb

!-----------------------------------------------------------------------

      subroutine resizes()

c Reading atomic basis sizes for different species.
c
c Reads fdf block. Not necessarily all species have to be given. The
c ones not given at input will be assumed to have the basis sizes
c given by the general input PAO.BasisSize, or its default value.

      type(block_fdf)            :: bfdf
      type(parsed_line), pointer :: pline

      integer isp

      if (fdf_block('PAO.BasisSizes',bfdf)) then
        do while(fdf_bline(bfdf,pline))
          if (.not. fdf_bmatch(pline,'nn'))
     .      call die("Wrong format in PAO.BasisSizes")
          isp = label2species(fdf_bnames(pline,1))
          if (isp .eq. 0) then
            write(6,'(a,1x,a)')
     .        "WRONG species symbol in PAO.BasisSizes:",
     .        trim(fdf_bnames(pline,1))
            call die()
          else
            basp => basis_parameters(isp)
            basp%basis_size = fdf_bnames(pline,2)
            call size_name(basp%basis_size)   !!! DEPRECATED
            write(6,'(4a)')
     .           'resizes: Read basis size for species ',
     .            trim(fdf_bnames(pline,1)),' = ',basp%basis_size
          endif
        enddo
      endif

      end subroutine resizes

!-----------------------------------------------------------------------
      subroutine read_polarization_scheme()

      ! What are considered as "polarization orbitals" (i.e., the first
      ! orbitals that are not occupied in the atomic ground-state) are
      ! commonly generated by formally applying an electric field to the
      ! orbital being polarized. This is the "perturbative" approach.
      ! The alternative is to treat the orbitals explicitly as part of
      ! another bona-fide shell. This "non-perturbative" approach can be
      ! made explicit in the PAO.Basis block.

      ! This routine contains code so that the scheme can be selected on
      ! a species-by-species basis through the use of the lighter
      ! 'PAO.Polarization.Scheme' block, or through a global fdf
      ! directive.

      use basis_types, only: qconf_options_t


      type(block_fdf)            :: bfdf
      type(parsed_line), pointer :: pline

      character(len=40) :: string
      logical           :: non_perturbative_pols
      logical           :: non_perturbative_pols_fallback
      integer isp, indexp

      type(qconf_options_t), pointer :: s
      
      ! Global options that might apply to all species
      
      non_perturbative_pols = 
     $     fdf_boolean('PAO.Polarization.NonPerturbative',
     $                 .false.)
      ! Allow legacy option
      non_perturbative_pols = 
     $     fdf_boolean('PAO.NonPerturbative.Polarization.Orbitals',
     $                 non_perturbative_pols)
      
      ! Enable fallback to non-perturbative polarization orbitals
      ! in cases where the algorithm might fail
      non_perturbative_pols_fallback = 
     $     fdf_boolean(
     $     'PAO.Polarization.NonPerturbative.Fallback',
     $                 .false.)

      loop: do isp=1, nsp
         basp=>basis_parameters(isp)

         basp%non_perturbative_polorbs = non_perturbative_pols
         basp%non_pert_polorbs_fallback =
     $                          non_perturbative_pols_fallback 

         if (non_perturbative_pols) then
            basp%non_pert_polorbs_req = .true.
         else
            basp%non_pert_polorbs_req = .false.
         endif
      end do loop

      ! Now process a block that allows per-species settings
      ! Note that this has precedence over the global setting

      if (fdf_block('PAO.PolarizationScheme',bfdf)) then
        do while(fdf_bline(bfdf,pline))

          if (.not. fdf_bmatch(pline,'nn'))  ! note: at least 2 strings
     .      call die("Wrong format in PAO.PolarizationScheme")

          isp = label2species(fdf_bnames(pline,1))
          if (isp .eq. 0) then
             write(6,'(a,1x,a)')
     .            "WRONG species symbol in PAO.PolarizationScheme:",
     .            trim(fdf_bnames(pline,1))
             call die("Wrong species in PAO.PolarizationScheme")
          endif

          basp => basis_parameters(isp)
             
          string = fdf_bnames(pline,2)
             
          select case (trim(string))
          case ("non-perturbative")
             
             basp%non_perturbative_polorbs = .true.
             basp%non_pert_polorbs_req = .true.
             
             ! Optional charge confinement options

             if (fdf_bsearch(pline,'Q',indexp)) then
                basp%qconf_options_polorbs_set = .true.
                s => basp%qconf_options_polorbs
                if (fdf_bmatch(pline,'vvv',after=indexp)) then
                   s%qcoe = fdf_bvalues(pline,ind=1,after=indexp)
                   s%qyuk = fdf_bvalues(pline,ind=2,after=indexp)
                   s%qwid = fdf_bvalues(pline,ind=3,after=indexp)
                elseif (fdf_bmatch(pline,'vv',after=indexp)) then
                   s%qcoe = fdf_bvalues(pline,ind=1,after=indexp)
                   s%qyuk = fdf_bvalues(pline,ind=2,after=indexp)
                   s%qwid = 0.01_dp
                elseif (fdf_bmatch(pline,'v',after=indexp)) then
                   s%qcoe = fdf_bvalues(pline,ind=1,after=indexp)
                   s%qyuk = 0.0_dp
                   s%qwid = 0.01_dp
                else
                   call die("Need one, two or three "//
     $                  "real numbers after 'Q' " //
     $                  "in PAO.Polarization.Scheme")
                endif
             endif
             
          case ("perturbative")
             
             if (basp%non_pert_polorbs_req) then
                write(6,'(a,1x,a)')
     .            "'Perturbative' Setting in PAO.PolarizationScheme" //
     $            " overrides PAO.Polarization.NonPerturbative" //
     $            " for:", trim(fdf_bnames(pline,1))
             endif
             basp%non_perturbative_polorbs = .false.
             basp%force_perturbative_polorbs = .true.
             
          case default
             call die("Bad keyword in PAO.PolarizationScheme")
          end select

          if (basp%non_perturbative_polorbs) then
             write(6,'(a,1x,a)')
     .            "Using non-perturbative " //
     $            "polarization scheme for",
     $            trim(basp%label)
          endif

        enddo
      endif


      end subroutine read_polarization_scheme

!-----------------------------------------------------------------------

      subroutine relmxkb()

c Reads the maximum angular momentum of the Kleinman-Bylander
c projectors for the different species.
c
c Reads fdf block. Not necessarily all species have to be given.

      type(block_fdf)            :: bfdf
      type(parsed_line), pointer :: pline

      integer isp

      if (fdf_block('PS.lmax',bfdf)) then
        do while(fdf_bline(bfdf,pline))
          if (.not. fdf_bmatch(pline,'ni'))
     .      call die("Wrong format in PS.lmax")
          isp = label2species(fdf_bnames(pline,1))
          if (isp .eq. 0) then
            write(6,'(a,1x,a)') "WRONG species symbol in PS.lmax:",
     .                           trim(fdf_bnames(pline,1))
            call die()
          else
            basp => basis_parameters(isp)
            basp%lmxkb_requested = fdf_bintegers(pline,1)
            write(6,"(a, i4, 2a)")
     .            'relmxkb: Read Max KB Ang. Momentum= ',
     .             basp%lmxkb_requested,
     .            ' for species ', trim(fdf_bnames(pline,1))
          endif
        enddo
      endif

      end subroutine relmxkb
!-----------------------------------------------------------------------

      subroutine remass()


      type(block_fdf)            :: bfdf
      type(parsed_line), pointer :: pline

      integer isp

c Read atomic masses of different species.
c
c Reads fdf block. Not necessarily all species have to be given. The
c ones not given at input will be assumed to have their natural mass
c (according to atmass subroutine).

      if (fdf_block('AtomicMass',bfdf)) then
        do while(fdf_bline(bfdf,pline))
          if (.not. fdf_bmatch(pline,'iv'))
     .      call die("Wrong format in AtomicMass")
          isp = fdf_bintegers(pline,1)                                 
          if (isp .gt. nsp .or. isp .lt. 1)                   
     .      call die("Wrong specnum in AtomicMass")        
          basp => basis_parameters(isp)                         
          basp%mass = fdf_bvalues(pline,2)                             
          write(6,"(a, i4, a, f12.5)")                        
     .         'remass: Read atomic mass for species ', isp,  
     .         ' as ', basp%mass                              
        enddo
      endif

      end subroutine remass
!-----------------------------------------------------------------------

      subroutine size_name(str)
      character(len=*), intent(inout)  ::  str


      if(leqi(str,'MINIMAL')) str='sz'
      if(leqi(str,'SZ'))  str='sz'
      if(leqi(str,'SZP')) str='szp'
      if(leqi(str,'SZP1')) str='szp'
      if(leqi(str,'SZSP')) str='szp'
      if(leqi(str,'SZ1P')) str='szp'
!
      if(leqi(str,'DZ')) str='dz'
      if(leqi(str,'STANDARD')) str='dzp'
      if(leqi(str,'DZP'))  str='dzp'
      if(leqi(str,'DZP1'))  str='dzp'
      if(leqi(str,'DZ1P'))  str='dzp'
      if(leqi(str,'DZSP'))  str='dzp'
      if(leqi(str,'DZP2'))  str='dzp2'
      if(leqi(str,'DZDP'))  str='dzp2'
      if(leqi(str,'DZ2P'))  str='dzp2'
!
      if(leqi(str,'TZ')) str='tz'
      if(leqi(str,'TZP')) str='tzp'
      if(leqi(str,'TZ1P')) str='tzp'
      if(leqi(str,'TZP1')) str='tzp'
      if(leqi(str,'TZSP')) str='tzp'
      if(leqi(str,'TZP2')) str='tzp2'
      if(leqi(str,'TZ2P')) str='tzp2'
      if(leqi(str,'TZDP')) str='tzp2'
      if(leqi(str,'TZP3')) str='tzp3'
      if(leqi(str,'TZ3P')) str='tzp3'
      if(leqi(str,'TZTP')) str='tzp3'
 
      if ( (str.ne.'szp').and.(str.ne.'sz').and.
     .     (str.ne.'dz') .and.(str.ne.'dzp') .and.
     .     (str.ne.'tz') .and.(str.ne.'tzp') .and.
     .     (str.ne.'dzp2') .and.
     .     (str.ne.'tzp2') .and. (str.ne.'tzp3') ) then

        write(6,'(/,2a,/,9(a,/))')
     .    'size_name: Incorrect basis-size option specified,',
     .    ' active options are:',
     .    '  SZ or MINIMAL', 
     .    '  SZP, SZSP, SZ1P, SZP1',
     .    '  DZ ',
     .    '  DZP, DZSP, DZP1, DZ1P or STANDARD',
     .    '  DZDP, DZP2, DZ2P ',
     .    '  TZ ',
     .    '  TZP, TZSP, TZP1, TZ1P',
     .    '  TZDP, TZP2, TZ2P',
     .    '  TZTP, TZP3, TZ3P'

        call die()
      endif

      end subroutine size_name
!-----------------------------------------------------------------------

      subroutine type_name(basistype)

      character basistype*(*)

      if(leqi(basistype,'NODES')) then
        basistype='nodes'
      elseif(leqi(basistype,'NONODES')) then
        basistype='nonodes'
      elseif(leqi(basistype,'SPLIT')) then
        basistype='split'
      elseif(leqi(basistype,'SPLITGAUSS')) then
        basistype='splitgauss'
      elseif (leqi(basistype,'FILTERET')) then
        basistype='filteret'
      else
        write(6,'(/,2a,(/,5(3x,a)),(/,2(3x,a)))')
     .    'type_name: Incorrect basis-type option specified,',
     .    ' active options are:',
     .    'NODES','SPLIT','SPLITGAUSS','NONODES','FILTERET'
         call die
      endif

      end subroutine type_name
!-----------------------------------------------------------------------
!
!     Find out whether semicore states are implied by the valence
!     charge density in the pseudopotential file.
!
      subroutine semicore_check(is)
      integer, intent(in)  :: is

      real(dp), parameter :: tiny = 1.d-5
      integer ndiff
      real(dp) zval, zval_vps, charge_loc

      basp => basis_parameters(is)

      basp%semic = .false.
      if (basp%bessel) return
      
      zval_vps = basp%pseudopotential%zval
      zval = basp%ground_state%z_valence
      
      if (abs(Zval-zval_vps).lt.tiny) return

      ndiff = nint(abs(Zval-zval_vps))
      if (abs(ndiff-abs(Zval-zval_vps)).gt.tiny) then
        write(6,'(2a)')
     .    'ERROR Fractional semicore charge for species ',
     .    basp%label
        call die()
      endif

      charge_loc = Zval_vps-Zval
      if (charge_loc > 0.0_dp) then
         basp%semic = .true.
         write(6,'(a,i2,a)')
     .        'Semicore shell(s) with ', nint(charge_loc),
     .        ' electrons included in the valence for ' //
     $        trim(basp%label)
      else
         write(6,'(a,i2,a)') 'Nominally valence shell(s) with ',
     $          nint(abs(charge_loc)),
     .        ' electrons are kept in the core for ' // trim(basp%label)
      endif

      end subroutine semicore_check
!----------------------------------------------------------------------
      subroutine autobasis()

      use m_semicore_info_froyen, only: get_n_semicore_shells
!
!     It sets the defaults if a species has not been included
!     in the PAO.Basis block
!
      integer l, nzeta, nzeta_pol
      integer :: nsemic_shells(0:3), nsh, i
      type(lshell_t), pointer :: ls_parent

      loop: do isp=1, nsp
         basp=>basis_parameters(isp)
         if (basp%in_pao_basis_block) cycle loop ! Species already set
                                                 ! in PAO.Basis block

         if (basp%bessel) then
            write(6,'(2a)') basp%label,
     .      ' must be in PAO.Basis (it is a floating Bessel function)'
            call die()
         endif
         !
         ! Set the default max l 
         !
         basp%lmxo = basp%ground_state%lmax_valence

         ! Check whether we need to consider larger l's due to
         ! semicore states.
         
         if (basp%bessel) then
            ! There are no semicore states 
            nsemic_shells(:) = 0

         else if (basp%synthetic) then

            write(6,'(a,1x,a)')
     .           "WARNING: Assuming absence of semicore " //
     .           "states for", trim(basp%label)
            write(6,'(a,1x,a)')
     .           "WARNING: If there are any, " //
     .           "use the PAO.Basis block"
            nsemic_shells(:) = 0

         else

            call get_n_semicore_shells(basp%pseudopotential,
     $                                 nsemic_shells)
         endif
         
         do l=0,ubound(nsemic_shells,dim=1)
            if (nsemic_shells(l) > 0) basp%lmxo = max(l,basp%lmxo)
         enddo
         
         ! Check whether we need to consider larger l's due to
         ! polarization orbitals (this is a new behavior -- in the
         ! case of perturbative pol orbs, the highest-l shell will remain empty)
         ! But lmxo in the 'specs' output will be correct, and not confusing

         if (basp%basis_size(3:3) .eq. 'p') then
            loop_pol: do l=1,4  ! Note that l starts at 1
            if ( (.not. basp%ground_state%occupied(l)) ) then
               basp%lmxo = max(l,basp%lmxo)
               EXIT loop_pol
            endif
            enddo loop_pol
         endif

         if (basp%basis_size(3:3) .eq. 'p') then
            loop_pol2: do l=1,4  ! Note that l starts at 1
            if ( (.not. basp%ground_state%occupied(l)) ) then
               ! Check whether there are already other (semicore) shells
               ! In that case, the pol orb will have a node. The algorithm
               ! to generate it is not robust enough, so, if the fallback has
               ! been enabled (it is not by default), we switch to generating
               ! the orbital explicitly (unless forced by the user in the block)
               if (nsemic_shells(l) >= 1) then
                  if (basp%non_pert_polorbs_fallback) then
                     if (.not. basp%force_perturbative_polorbs) then
                        basp%polorb_with_semicore = .true.
                        basp%non_perturbative_polorbs = .true.
                        write(6,'(a,1x,a)')
     .                     "Fallback to non-perturbative " //
     $                     "polarization scheme for",
     $                     trim(basp%label)
                     endif
                  endif
               endif
               EXIT loop_pol2
            endif
            enddo loop_pol2
         endif
        
         allocate (basp%lshell(0:basp%lmxo))

         if (basp%basis_size(1:2) .eq. 'sz') nzeta = 1
         if (basp%basis_size(1:2) .eq. 'dz') nzeta = 2
         if (basp%basis_size(1:2) .eq. 'tz') nzeta = 3

         loop_l: do l=0, basp%lmxo
            ls=>basp%lshell(l)
            call initialize(ls)
            ls%l = l
            !  Include semicore shells
            nsh = 1+nsemic_shells(l)
            ls%nn = nsh
            !
            allocate(ls%shell(1:nsh))
            !     
            ! Inner shells get a single-z orbital
            !
            do i = 1, nsh -1
               s => ls%shell(i)
               call initialize(s)
               s%l = l
               ! e.g.: nsh=3; i=1:  n-2; i=2: n-1
               s%n = basp%ground_state%n(l) - (nsh-i)
               s%nzeta = 1
               s%polarized = .false.
               s%polarization_shell = .false.
               s%split_norm = global_splnorm
               s%nzeta_pol = 0

               if (lsoft) then
                  s%vcte = softPt 
                  s%rinn = -softRc
               else
                  s%rinn = 0.d0
                  s%vcte = 0.d0
               endif

               ! Default filteret cutoff for shell
               s%filtercut = 0.0d0

               if (s%nzeta .ne.0) then
                  allocate(s%rc(1:s%nzeta))
                  allocate(s%lambda(1:s%nzeta))
                  s%rc(1:s%nzeta) = 0.0d0
                  s%lambda(1:s%nzeta) = 1.0d0
               endif
            enddo
            ! Do the final (canonical valence)  shell
            !  
            s => ls%shell(nsh)
            call initialize(s)
            s%l = l
            s%n = basp%ground_state%n(l)
            if (basp%ground_state%occupied(l)) then   
               s%nzeta = nzeta
            else
               s%nzeta = 0  ! This might be a polarization orbital... **, and still be counted as a "shell"
            endif
            s%polarized = .false.
            s%polarization_shell = .false.
            if (abs(basp%z).eq.1) then
               s%split_norm = global_splnorm_H
            else
               s%split_norm = global_splnorm
            endif
            s%nzeta_pol = 0

            if (lsoft) then
               s%vcte = softPt 
               s%rinn = -softRc
            else
               s%rinn = 0.d0
               s%vcte = 0.d0
            endif

            ! Default filteret cutoff for shell
            s%filtercut = 0.0d0

            if (s%nzeta .ne.0) then
               allocate(s%rc(1:s%nzeta))
               allocate(s%lambda(1:s%nzeta))
               s%rc(1:s%nzeta) = 0.0d0
               s%lambda(1:s%nzeta) = 1.0d0
               !
               ! This option needs scale factors for multiple-z
               if (basp%basis_type.eq.'nonodes') then
                  s%lambda(1) = 1.0d0
                  if (s%nzeta > 1)  then
                     write(6,"(a)") "WARNING: NONODES: " //
     $                 "Default scale factors set for shell (x0.8)"
                  endif
                  do i=2,s%nzeta
                     s%lambda(i) = 0.8d0 * s%lambda(i-1)
                  enddo
               endif
            endif
         enddo loop_l

         if (basp%basis_size(3:3) .eq. 'p') then

         ! Polarization orbitals are defined so they have the minimum      
         ! angular momentum l such that there are no occupied orbitals 
         ! with the same l in the valence shell of the ground-state 
         ! atomic configuration. They polarize the corresponding l-1    
         ! shell.

         ! note that we go up to l=4 to make the loop simpler.
         ! (that is the reason why 'occupied' is dimensioned to 0:4)

            select case (basp%basis_size(4:4))
              case (' ', '1')
                 nzeta_pol = 1
              case ('2')
                 nzeta_pol = 2
              case ('3')
                 nzeta_pol = 3
            end select

            loop_angmom: do l=1,4
               ! For example, in Ti, p is not occupied in the GS
               ! The following code will mark the top-most s shell
               ! as polarizable

               ! Note that l goes from 1 to 4
               if ( (.not. basp%ground_state%occupied(l)) ) then
    
                  if (basp%non_perturbative_polorbs) then
                     ls=>basp%lshell(l)
                     nsh = size(ls%shell)
                     ! This can happen, for example, for Ge (3d is semicore,
                     ! but the polarization shell is 4d (from 4p)
                     !if (nsh /= 1) call die("Empty l-shell with nsh>1?")
                     s => ls%shell(nsh)
                     s%nzeta = nzeta_pol
                     allocate(s%rc(1:s%nzeta))
                     allocate(s%lambda(1:s%nzeta))
                     s%rc(1:s%nzeta) = 0.0d0
                     s%lambda(1:s%nzeta) = 1.0d0
                     
                     s%polarization_shell = .true.
                     
                     ! Optional charge confinement
                     ! for non-perturbative polarization orbitals
                     ! This does not seem to be really necessary for
                     ! the "fallback" cases in which there are semicore
                     ! states with the same l as the polarization orbital.
                     if (basp%qconf_options_polorbs_set) then
                        s%qcoe = basp%qconf_options_polorbs%qcoe
                        s%qyuk = basp%qconf_options_polorbs%qyuk
                        s%qwid = basp%qconf_options_polorbs%qwid
                     else
                        ! (would apply the same to all species)
                        s%qcoe = fdf_get(
     $                       "PAO.Polarization.Charge.Confinement",
     $                       0.0_dp)
                     endif
                     
                     ls_parent=> basp%lshell(l-1)
                     nsh = size(ls_parent%shell) ! Parent is the last (l-1) shell
                     s%shell_being_polarized => ls_parent%shell(nsh)
                     ! Now, reset the polarization flags of the parent
                     ls_parent%shell(nsh)%polarized = .false.
                     ls_parent%shell(nsh)%was_polarized = .true.
                     ls_parent%shell(nsh)%nzeta_pol = 0
                     
                     EXIT loop_angmom       ! Polarize only one shell!
                  endif

                  ! Normal treatment 
                  ls=>basp%lshell(l-1)
                  nsh = size(ls%shell) ! Use only the last shell
                  s => ls%shell(nsh)

              ! Check whether shell to be polarized is occupied in the gs.
              ! (i.e., whether PAOs are going to be generated for it)
              ! If it is not, mark it for PAO generation anyway.
              ! This will happen for confs of the type s0 p0 dn

                  ! Default filter cutoff for shell
                  s%filtercut = 0.0d0

                  if (s%nzeta == 0) then
                     write(6,"(a,i2,a)") "Marking shell with l=",
     .                l-1, " for PAO generation and polarization."
                     s%nzeta = nzeta
                     allocate(s%rc(1:s%nzeta))
                     allocate(s%lambda(1:s%nzeta))
                     s%rc(1:s%nzeta) = 0.0d0
                     s%lambda(1:s%nzeta) = 1.0d0
                  endif
                  ! Put polarization flags
                  s%polarized = .true.
                  s%nzeta_pol = nzeta_pol

                  EXIT loop_angmom  ! Polarize only one shell!
               endif
            enddo loop_angmom

         endif

      enddo loop
            
      end subroutine autobasis

      End module basis_specs
