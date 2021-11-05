
! Copyright (C) 1996-2016       The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!

!> \brief General purpose of the m_ldau_so module:
!! Compute the LSDA+U potential in the case of non-collinear magnetism
!! or spin orbit coupling.
!!
!! Starting from Eq. (5) of the paper by Liechtenstein {\it et al.} 
!! \cite Liechtenstein-95,
!!
!! \f{eqnarray*}{
!!   V_{m m^\prime}^{\sigma} = \sum_{m^{\prime \prime} m^{\prime\prime\prime}} &
!!   \left[ \langle m, m^{\prime \prime} \vert
!!          V_{ee} \vert m^{\prime},m^{\prime\prime\prime} \rangle
!!          n_{m^{\prime \prime}m^{\prime\prime\prime}}^{-\sigma} \right.  +
!!   \nonumber \\
!!        & \left. \left(
!!        \langle m, m^{\prime \prime} \vert
!!          V_{ee} \vert m^{\prime},m^{\prime\prime\prime} \rangle -
!!        \langle m, m^{\prime \prime} \vert
!!          V_{ee} \vert m^{\prime\prime\prime}, m^{\prime} \rangle
!!        \right) n_{m^{\prime \prime}m^{\prime\prime\prime}}^{\sigma}
!!   \right]
!!   \nonumber \\
!!   & - U \left( n - \frac{1}{2} \right) + 
!!   J \left( n^{\sigma} - \frac{1}{2} \right),
!! \f}
!! where
!! \f{eqnarray*}{
!!   n^{\sigma} & = \mathrm{Tr} \left( n_{m m^{\prime}}^{\sigma} \right)
!!   \nonumber \\
!!   & = \sum_{m m^{\prime}} n_{m m^{\prime}}^{\sigma} \delta_{m m^{\prime}} ,
!! \f}
!! and
!! \f{eqnarray*}{
!!  n = n^{\uparrow} + n^{\downarrow}
!! \f}
!! Replacing the two last Equations in the double counting terms 
!! of the first Equation,
!! \f{eqnarray*}{
!!   V_{m m^\prime}^{\sigma} = \sum_{m^{\prime \prime} m^{\prime\prime\prime}} &
!!   \left[ \langle m, m^{\prime \prime} \vert
!!          V_{ee} \vert m^{\prime},m^{\prime\prime\prime} \rangle
!!          n_{m^{\prime \prime}m^{\prime\prime\prime}}^{-\sigma} \right.  +
!!   \nonumber \\
!!        & \left. \left(
!!        \langle m, m^{\prime \prime} \vert
!!          V_{ee} \vert m^{\prime},m^{\prime\prime\prime} \rangle -
!!        \langle m, m^{\prime \prime} \vert
!!          V_{ee} \vert m^{\prime\prime\prime}, m^{\prime} \rangle
!!        \right) n_{m^{\prime \prime}m^{\prime\prime\prime}}^{\sigma}
!!   \right]
!!   \nonumber \\
!!   & - U \left[ \sum_{m^{\prime\prime} m^{\prime \prime \prime}} \left(
!!       n_{m^{\prime\prime} m^{\prime \prime \prime}}^{\sigma} \delta_{m^{\prime\prime} m^{\prime \prime \prime}} +
!!       n_{m^{\prime\prime} m^{\prime \prime \prime}}^{-\sigma} \delta_{m^{\prime\prime} m^{\prime \prime \prime}}\right) \right] + U \frac{1}{2}
!!   \nonumber \\
!!   & + J \left[ \sum_{m^{\prime\prime} m^{\prime \prime \prime}} \left(  n_{m^{\prime\prime} m^{\prime \prime \prime}}^{\sigma} \delta_{m^{\prime\prime} m^{\prime \prime \prime}} \right) \right] - J \frac{1}{2}
!!   \nonumber \\
!!   = \sum_{m^{\prime \prime} m^{\prime\prime\prime}} &
!!   \left[ \left( \langle m, m^{\prime \prime} \vert
!!          V_{ee} \vert m^{\prime},m^{\prime\prime\prime} \rangle - U \delta_{m^{\prime\prime} m^{\prime \prime \prime}}   \right)
!!          n_{m^{\prime \prime}m^{\prime\prime\prime}}^{-\sigma} \right.  +
!!   \nonumber \\
!!        & \left. \left(
!!        \langle m, m^{\prime \prime} \vert
!!          V_{ee} \vert m^{\prime},m^{\prime\prime\prime} \rangle -
!!        \langle m, m^{\prime \prime} \vert
!!          V_{ee} \vert m^{\prime\prime\prime}, m^{\prime} \rangle
!!        - (U - J) \delta_{m^{\prime\prime} m^{\prime \prime \prime}}
!!        \right) n_{m^{\prime \prime}m^{\prime\prime\prime}}^{\sigma}
!!   \right]
!!   \nonumber \\
!!   & + \frac{1}{2} \left( U - J \right)
!! \f}
!! If we replace
!! \f{eqnarray*}{
!!   m & \rightarrow 1,
!!   \nonumber \\
!!   m^{\prime} & \rightarrow 2,
!!   \nonumber \\
!!   m^{\prime \prime} & \rightarrow 3,
!!   \nonumber \\
!!   m^{\prime \prime \prime} & \rightarrow 4,
!! \f}
!! in the former Equation, we arrive to Eq.(3) in the paper
!! by Bousquet and Spaldin \cite Bousquet-10,
!! \f{eqnarray*}{
!!    V_{1, 2}^{\sigma} = \sum_{3,4} &
!!   \left[ \left( \langle 1, 3 \vert
!!          V_{ee} \vert 2, 4 \rangle - U \delta_{3,4}   \right)
!!          n_{3,4}^{-\sigma} \right.  +
!!   \nonumber \\
!!        & \left. \left(
!!        \langle 1, 3 \vert
!!          V_{ee} \vert 2, 4 \rangle -
!!        \langle 1, 3 \vert
!!          V_{ee} \vert 4, 2 \rangle
!!        - (U - J) \delta_{3,4}
!!        \right) n_{3,4}^{\sigma}
!!   \right]
!!   \nonumber \\
!!   & + \frac{1}{2} \left( U - J \right)
!! \f}
!!
!! For the off-diagonal term, and following Ref.\cite Bousquet-10,
!! \f{eqnarray*}{
!!   V_{m m^\prime}^{\uparrow\downarrow} = 
!!        \sum_{m^{\prime \prime} m^{\prime\prime\prime}} &
!!        \left(
!!        - \langle m, m^{\prime \prime} \vert
!!          V_{ee} \vert m^{\prime\prime\prime}, m^{\prime} \rangle
!!          + J \delta_{m^{\prime\prime}m^{\prime\prime\prime}} \right)
!!         n_{m^{\prime \prime}m^{\prime\prime\prime}}^{\uparrow\downarrow}
!! \f}
!! and the corresponding for \f$V_{m m^\prime}^{\downarrow \uparrow}\f$



module m_ldau_so

  use precision,       only : dp            ! Double precision
  use parallel,        only : Node          ! Local processor number
  use parallel,        only : Nodes         ! Total number of processors in a
                                            !   parallel run
  use parallelsubs,    only : LocalToGlobalOrb
                                            ! Converts an orbital index in the 
                                            !   local frame to the global frame
  use atm_types,       only : nspecies      ! Number of species
  use atm_types,       only : species_info  ! Derived type with all the info
                                            !   about the radial functions
                                            !   (PAOs, KB projectors,
                                            !   LDA+U proj,
                                            !   VNA potentials, etc)
                                            !   for a given atomic specie
  use atm_types,       only : species       ! Actual array where the
                                            !   previous information is
                                            !   stored
  use atm_types,       only : ldau_so_integrals_type
                                            ! Derived type for the
                                            !   definition of the on-site
                                            !   four-center-integrals
                                            !   required for LDA+U+Spinorbit
  use atomlist,        only : no_l          ! Number of orbitals (local)
  use atomlist,        only : iaorb         ! Atomic index of each orbital
  use atomlist,        only : iphorb        ! Orbital index of each  orbital 
                                            !   in its atom
  use atomlist,        only : indxuo        ! Index of equivalent orbital in 
                                            !   the unit cell
  use atmfuncs,        only : cnfigfio      ! Returns principal quantum number
  use atmfuncs,        only : lofio         ! Returns angular mumentum number
  use atmfuncs,        only : mofio         ! Returns magnetic quantum number
  use atmfuncs,        only : pol           ! Returns whether an orbital is 
                                            !   polarized
  use atmfuncs,        only : symfio        ! Returns angular symmetry of orb
  use atmfuncs,        only : zetafio       ! Returns zeta number of orbital
  use siesta_geom,     only : na_u          ! Number of atoms in the unit cell
  use siesta_geom,     only : isa           ! Species index of each atom
  use sparse_matrices, only : maxnh         ! Number of non-zero elements in
                                            !   the sparse matrix
                                            !   (local to MPI Node)
  use sparse_matrices, only : numh          ! Number of non-zero elements per 
                                            !   row in the sparse matrices
                                            !   (local to MPI Node)
  use sparse_matrices, only : listh         ! Column indices in the sparse 
                                            !   matrices (local to MPI Node)
  use sparse_matrices, only : listhptr      ! Index pointer to `listh` for 
                                            !   each row in the CSR matrix, 
                                            !   0-based 
                                            !   (local to MPI Node)
  use sparse_matrices, only : Dscf          ! Sparse density matrix
  use m_spin,          only : spin          ! Variable containing the 
                                            !   information related with the 
                                            !   spin configuration
  use m_energies,      only : E_correc_dc   ! Correction energy required for 
                                            !   the LDA+U+SO calculations
  use m_energies,      only : E_ldau_so     ! 
!
! Allocation/Deallocation routines
!
  use alloc,           only: re_alloc       ! Reallocation routines
  use alloc,           only: de_alloc       ! Deallocation routines

#ifdef MPI
  use mpi_siesta
#endif


  implicit none

  integer, allocatable, save, dimension(:,:) :: nzeta
                                            ! Number of zeta radial functions 
                                            !   for the atom indicated in 
                                            !   the first index,
                                            !   for the correlated shell 
                                            !   indicated in the second index
  integer, allocatable, save, dimension(:)   :: index_at_correlated

  public ldau_so_hamil

  private


  CONTAINS 

  subroutine ldau_so_hamil( S, H_ldau_so )

     real(dp), intent(in) :: S(maxnh)
     complex(dp), intent(inout) :: H_ldau_so(maxnh,4)

     type(species_info),           pointer :: spp
     type(ldau_so_integrals_type), pointer :: ldauintegrals

     logical, save :: firstime = .true.          ! First time that this
                                      !   subroutine is called?
     integer, save :: maxldau = 0     ! Maximum number of 
                                      !   LDA+U projectors on a 
                                      !   given atom
                                      !   without including "m" copies
     integer, save :: maxl    = 0     ! 
     integer, save :: nat_correlated  ! Number of atoms in the unit cell with
                                      !   at least one correlated shell 
     integer, save :: maxnzeta        ! Maximum number of zetas
     integer :: ka                    ! Counter for loop on atoms
     integer :: iat_correlated        ! 
     integer :: ish                   ! Counter for loop on species
     integer :: iproj                 ! Counter for loop on projector
     integer :: ishell                ! Counter for loop on shells
     integer :: ildauso               ! Counter for loop on 4 center integrals
     integer :: l                     ! Angular quantum number
     integer :: m                     ! Magnetic quantum number of the 
                                      !   first  orbital
     integer :: mprime                ! Magnetic quantum number of the 
                                      !   second  orbital
     integer :: m2prime               ! Magnetic quantum number of the 
                                      !   third  orbital
     integer :: m3prime               ! Magnetic quantum number of the
                                      !   fourth  orbital
     integer :: io_local              ! Counter for the local atomic orbitals
                                      !   in a given node
     integer :: io_global             ! Global index of an atomic orbital
     integer :: j_neig                ! Counter for the loop on neighbours
     integer :: neig_orb              ! Index of the neighbour orbital
     integer :: ind                   ! Index of the neighbour orbital in listh
                                      !   between m and mprime
     integer :: ia                    ! Atom to which orbital belongs
     integer :: ia_neig               ! Atom to which orbital belongs
     integer :: is                    ! Atomic species index
     integer :: is_neig               ! Atomic species index
     integer :: iuo                   ! Equivalent orbital in first unit cell
     integer :: iua                   ! Equivalent atom in first unit cell
     integer :: iao                   ! Orbital index within atom
     integer :: iao_neig              ! Orbital index within atom
     integer :: n_qn                  ! Orbital's principal quantum number
     integer :: n_qn_neig             ! Orbital's principal quantum number
     integer :: l_qn                  ! Orbital's angular momentum number
     integer :: l_qn_neig             ! Orbital's angular momentum number
     integer :: zeta_qn               ! 'Zeta' index of orbital
     integer :: zeta_qn_neig          ! 'Zeta' index of orbital
     integer :: iz2prime              ! Counter for loop on zeta
     integer :: iz3prime              ! Counter for loop on zeta
     integer :: m_qn                  ! (Real) orbital's magnetic quantum number
     integer :: m_qn_neig             ! (Real) orbital's magnetic quantum number
     integer :: m_qn_2neig            ! (Real) orbital's magnetic quantum number
     integer :: m_qn_3neig            ! (Real) orbital's magnetic quantum number
     integer :: icounter              ! 
     logical :: polarized             ! Is this a polarization orbital?
     character(len=32) :: orb_sym     ! Name of orbital's angular symmetry
     character(len=32) :: orb_sym_neig! Name of orbital's angular symmetry
     integer :: n_correlated          ! Principal quantum number
                                      !   of the correlated shell
     integer :: l_correlated          ! Angular momentum quantum number
                                      !   of the correlated shell
     real(dp) :: U                    ! Value of U
     real(dp) :: J                    ! Value of J
     real(dp) :: direct_int           ! Direct four center integrals
     real(dp) :: fock_int             ! Fock (exchange) four center integrals
     real(dp) :: number_corr_elec_loc ! Number of correlated electrons
                                      !   in the local node
     real(dp) :: number_corr_elec     ! Number of correlated electrons
                                      !   sum over all nodes

     real(dp) :: E_correc_dc 
     real(dp) :: magnetic_field(3)
     complex(dp), parameter :: imag = cmplx( 0.0_dp, 1.0_dp, kind=dp )
     complex(dp), parameter :: one  = cmplx( 1.0_dp, 0.0_dp, kind=dp )
     complex(dp) ::  Dscf_cmplx_1
     complex(dp) ::  Dscf_cmplx_2
     complex(dp) ::  Dscf_cmplx_3
     complex(dp) ::  Dscf_cmplx_4
     complex(dp) ::  H_ldau_so_Hubbard(maxnh,4)
     complex(dp) ::  H_ldau_so_dc(maxnh,4)
     complex(dp) ::  H_Zeeman(maxnh,4)
     complex(dp) ::  magnetic_loc(3)
     complex(dp) ::  magnetic(3)
     real(dp), dimension(:,:,:,:), pointer :: Dscf_correlated => null()

#ifdef MPI
     integer     :: MPIerror
     real(dp), dimension(:,:,:,:), pointer :: auxloc => null()
                     ! Temporal array for the
                     !   the global reduction of Dscf_correlated
#endif


!    Initialization and allocation of matrices
     if( firstime ) then

!      Find maximum number of LDA+U projectors on a given atom
!      not including the "m copies"
       if( allocated(index_at_correlated) ) deallocate(index_at_correlated)
       allocate( index_at_correlated(na_u) )
       nat_correlated = 0
       index_at_correlated(:) = 0
       do ka = 1, na_u
         is  = isa(ka)
         spp  => species(is)
         if( spp%n_pjldaunl .gt. 0 ) then
           nat_correlated = nat_correlated + 1
           maxldau = max(maxldau,spp%n_pjldaunl)
           index_at_correlated(ka) = nat_correlated
         endif
       enddo 

!      Allocate the number of zeta functions in the correlated shell
       if( allocated(nzeta) ) deallocate(nzeta)
       allocate( nzeta(na_u,maxldau) )
       nzeta(:,:) = 0
       maxnzeta   = 0
       maxl = 0

!      Find maximum number of zetas on the correlated shell
       do ka = 1, na_u
         is  = isa(ka)
         spp => species(is)
         do iproj = 1, spp%nprojsldau   ! Number of LDA+U projectors
                                        !   counting the "m copies"
                                        !   This loop must include the 
                                        !   "m copies" since spp%pjldau_n
                                        !   is the only array containing
                                        !   the actual principal quantum numbe
                                        !   of the correlated shell,
                                        !   that will be compared with that
                                        !   of the atomic orbital shell

!!          For debugging
!           write(6,'(a,9i5)')                                                           &   
! &           'ldau_so_hamil: Node, Nodes, atom, species, iproj, index, n, l, m = ',     &
! &           Node, Nodes, ka, is, iproj,                                                &
! &           spp%pjldau_index(iproj), spp%pjldau_n(iproj),                              &
! &           spp%pjldau_l(iproj), spp%pjldau_m(iproj)
!!          End debugging

!          Here we know:
!          the principal quantum number
!          the angular quantum number 
!          of the correlated shell
!          We are going to detect how many zetas of the atomic orbitals
!          do we have for this shell
           do ishell = 1, spp%n_orbnl    ! Loop on atomic orbitals
             if( (spp%orbnl_n(ishell) .eq. spp%pjldau_n(iproj) ) .and.      &
 &               (spp%orbnl_l(ishell) .eq. spp%pjldau_l(iproj) ) ) then

!!              For debugging
!               write(6,'(a,8i5,l5)')                                                       &
! &              '  ldau_so_hamil: Node, Nodes, atom, species, ishell, n, l, zeta, pol = ', &
! &              Node, Nodes, ka, is, ishell,                                               &
! &              spp%orbnl_n(ishell), spp%orbnl_l(ishell), spp%orbnl_z(ishell),             &
! &              spp%orbnl_ispol(ishell)
!!              End debugging
               maxl = max( maxl, spp%pjldau_l(iproj) )
               nzeta(ka,spp%pjldau_index(iproj)) =                             &
 &               max( nzeta(ka,spp%pjldau_index(iproj)),spp%orbnl_z(ishell) )
             endif
           enddo   ! end loop in the number of atomic orbitals
           maxnzeta = max( maxnzeta, nzeta(ka,spp%pjldau_index(iproj)) )
         enddo     ! end loop in the number of projectors including "m" copies
       enddo     ! End loop on atoms in the unit cell

!!      For debugging
!       do ka = 1, na_u   
!         is  = isa(ka)
!         spp  => species(is)
!         write(6,'(/a,5i5)')                                          &
! &         'ldau_so_hamil: Node, Nodes, atom, maxldau, maxnzeta = ',  &
! &           Node, Nodes, ka, maxldau, maxnzeta
!         do iproj = 1, spp%n_pjldaunl   
!           write(6,'(a,5i5)')                                         &
! &           'ldau_so_hamil: Node, Nodes, atom, iproj, nzeta   = ',   &
! &            Node, Nodes, ka, iproj, nzeta(ka,iproj)  
!         enddo 
!       enddo   
!       do ka = 1, na_u
!         is  = isa(ka)
!         spp  => species(is)
!         do iproj = 1, spp%n_pjldaunl   ! Number of LDA+U projectors
!           write(6,'(a,7i5,2f12.5)')                                                    & 
! &           'ldau_so_hamil: Node, Nodes, atom, species, iproj, n, l, U, J = ',         &
! &           Node, Nodes, ka, is, iproj, spp%pjldaunl_n(iproj), spp%pjldaunl_l(iproj),  &
! &           spp%pjldaunl_U(iproj), spp%pjldaunl_J(iproj)
!           do ildauso = 0, 2*l
!             write(6,'(a,3i5,f12.5)')                                                   & 
! &             'ldau_so_hamil: Node, Nodes, index, Slater = ',                          &
! &             Node, Nodes, ildauso, ldauintegrals%Slater_F(ildauso)         
!           enddo
!           do m = 1, 2*l+1
!             do mprime = 1, 2*l+1
!               do m2prime = 1, 2*l+1
!                 do m3prime = 1, 2*l+1
!            write(6,'(a,6i5,f12.5)')                                                    &
! &  'ldau_so_hamil: Node, Nodes, m, mprime, m2prime, m3prime, vee_integral_cmplx = ',   &
! &            Node, Nodes, m, mprime, m2prime, m3prime,                                 &
! &            ldauintegrals%vee_4center_integrals(m, mprime, m2prime, m3prime)
!                 enddo 
!               enddo 
!             enddo 
!           enddo 
!         enddo 
!       enddo 
!
!       write(6,'(a,3i5)')                                                     &
!  &      'ldau_so_hamil: Node, Nodes, number correlated atoms             = ',&
!  &      Node, Nodes, nat_correlated
!       write(6,'(a,3i5)')                                                     &
!  &      'ldau_so_hamil: Node, Nodes, maximum number of LDA+U projectors  = ',&
!  &      Node, Nodes, maxldau
!       write(6,'(a,3i5)')                                                     &
!  &      'ldau_so_hamil: Node, Nodes, maximum value of l for the proj.    = ',&
!  &      Node, Nodes, maxl   
!       write(6,'(a,3i5)')                                                     &
!  &      'ldau_so_hamil: Node, Nodes, maximum number of zetas             = ',&
!  &      Node, Nodes, maxnzeta      
!       do ka = 1, na_u
!         write(6,'(a,4i5)')                                                     &
!  &        'ldau_so_hamil: Node, Nodes, iat, index_at_correlated            = ',&
!  &        Node, Nodes, ka, index_at_correlated(ka)
!       enddo
!       call die()
!!      End debugging

       nullify(  Dscf_correlated )
       allocate( Dscf_correlated(nat_correlated, maxldau, (2*maxl+1)**2 * maxnzeta**2, 8) )
       Dscf_correlated(:,:,:,:) = 0.0_dp
     endif   ! End if first-time
     firstime = .false.

!    Find the elements of the density matrix between orbitals
!    belonging to the correlated shell of an atom.
!    Remember that there might be more than one atomic shell with 
!    the same correlated l
       
!    Loop on all the atomic orbitals stored in the local node
     Dscf_correlated(:,:,:,:) = 0.0_dp
     do io_local = 1, no_l

!      Converts the orbital index from local frame to global frame
       call LocalToGlobalOrb(io_local,Node,Nodes,io_global)

!      Find the attributes of that particular orbital
       ia = iaorb(io_global)        ! Atom to which orbital belongs
       is = isa(ia)                 ! Atomic species index
       iuo = indxuo(io_global)      ! Equivalent orbital in first unit cell
       iua = iaorb(iuo)             ! Equivalent atom in first unit cell
       iao = iphorb(io_global)      ! Orbital index within atom
       n_qn = cnfigfio( is, iao)    ! Orbital's principal quantum number
       l_qn = lofio( is, iao )      ! Orbital's angular mumentum number
       m_qn = mofio( is, iao )      ! (Real) orbital's magnetic quantum number
       zeta_qn = zetafio( is, iao ) ! 'Zeta' index of orbital
       polarized = pol( is, iao )   ! Is this a polarization orbital?
       orb_sym = symfio( is, iao )  ! Name of orbital's angular symmetry
  
       spp  => species(is)

       do iproj = 1, spp%n_pjldaunl  ! Number of LDA+U projectors for this
                                     !   atomic species
                                     !   not including the "m copies"
         ldauintegrals => spp%ldau_so_integrals(iproj)
         l_correlated = spp%pjldaunl_l(iproj)

         if( l_qn .eq. l_correlated ) then ! We need to compute the new
                                           !    Hamiltonian matrix elements
                                           !    between atomic orbitals
                                           !    that belongs to the correlated
                                           !    shell

!!          For debugging
!           write(6,'(a,12i5,l5,1x,a10)')                                    &
! &           'ldau_so_hamil: Nodes, Node, io_l, io_g, ia, is, iuo, iua = ', &
! &                 Nodes, Node, io_local, io_global, ia, is, iuo, iua,      &
! &                 n_qn, l_qn, m_qn, zeta_qn, polarized, orb_sym 
!!          End debugging

!          Loop on all the neighbours of the atomic orbital
           do j_neig = 1, numh(io_local) 
             ind = listhptr(io_local) + j_neig

!            Identify index of the neighbour orbital, global frame
             neig_orb = listh(ind)

!            Identify relevant information for the neighbour orbital
!            Meaning of the different variables as before
             ia_neig      = iaorb(neig_orb)  
             is_neig      = isa(ia_neig) 
             iao_neig     = iphorb(neig_orb)
             n_qn_neig    = cnfigfio( is_neig, iao_neig)  
             l_qn_neig    = lofio( is_neig, iao_neig )
             m_qn_neig    = mofio( is_neig, iao_neig )  
             zeta_qn_neig = zetafio( is_neig, iao_neig ) 
             orb_sym_neig = symfio( is_neig, iao_neig ) 

             if( l_qn_neig .eq. l_correlated) then
               if( ia .eq. ia_neig ) then 
                 iat_correlated = index_at_correlated(ia)
                 icounter = (zeta_qn-1) * (2*l_correlated+1)**2 * maxnzeta                 + &
 &                          (m_qn + 1 + l_correlated - 1)  * (2*l_correlated+1) * maxnzeta + &
 &                          (zeta_qn_neig-1) * (2*l_correlated+1)                          + &
 &                          (m_qn_neig + 1 + l_correlated)
                 Dscf_correlated(iat_correlated,iproj,icounter,:) = Dscf(ind,:)

!!                For debugging
!                 write(6,'(a,11(i5,2x))')                                                         &
! &               '  ldau_so_hamil: Nodes, Node, iol, iog, j_neig, ind, listh, n, l, m =',         &
! &               Nodes, Node, io_local, io_global, j_neig, ind, neig_orb, n_qn_neig,              &
! &               l_qn_neig, m_qn_neig, zeta_qn_neig
!                 write(6,'(a,13(i7,2x))')                                                             &
! &               '  ldau_so_hamil: Nodes, Node, ia, ind, m_qn, zeta_qn, m_qn_neig, zeta_qn_neig, icounter =',  &
! &                                 Nodes, Node, ia, ind, m_qn, zeta_qn, m_qn_neig, zeta_qn_neig, icounter
!!                End debugging

               endif
             endif ! Close condition on the angular momentum of the neighbour
           enddo   ! Close loop on neighbours (j_neig)
         endif     ! Close condition on the angular momentum of the uc orbital

!!        For debugging
!         do zeta_qn = 1, nzeta(ia,iproj)
!           do zeta_qn_neig = 1, nzeta(ia,iproj)
!             do m_qn = -l, l
!               do m_qn_neig = -l, l
!                 write(6,'(a,7i5,4(i7,2x))')                                         &
! &                 'ldau_so_hamil: Node, Nodes, io_l, io_g, ia, iproj, l, index = ', &
! &                 Node, Nodes, io_local, io_global, ia, iproj, l_correlated,        &
! &                 m_qn, zeta_qn, m_qn_neig, zeta_qn_neig                           
!               enddo
!             enddo
!           enddo 
!         enddo
!!        End debugging

       enddo       ! Close loop on LDA+U projectors on the atom within the uc
     enddo         ! Close loop on the orbitals in the uc

#ifdef MPI
!    Allocate workspace array for global reduction
     call re_alloc( auxloc,                         & 
 &       1, nat_correlated,                         &
 &       1, maxldau,                                &
 &       1, (2*maxl+1)**2 * maxnzeta**2,            &
 &       1, 8,                                      &
 &       name='auxloc', routine='ldau_so_hamil' )
!    Global reduction of auxloc matrix
     auxloc(:,:,:,:) = 0.0_dp
     call MPI_AllReduce( Dscf_correlated(1,1,1,1),                            &
 &                       auxloc(1,1,1,1),                                     &
 &                       nat_correlated * maxldau *                           & 
 &                       (2*maxl+1)**2 * maxnzeta**2 * 8,             &
 &                       MPI_double_precision,MPI_sum,MPI_Comm_World,MPIerror )
     Dscf_correlated(:,:,:,:) = auxloc(:,:,:,:)
#endif

!!    For debugging
!     do ka = 1, na_u
!       if( index_at_correlated(ka) .ne. 0 ) then
!         is = isa(ka) 
!         spp  => species(is)
!         do iproj = 1, spp%n_pjldaunl  
!          write(6,*)ka, iproj
!           do j_neig = 1, (2*maxl+1)**2 * maxnzeta**2
!             do m = 1, 8
!               write(6,'(a,6i5,f12.5)')                                            &
! &               'Node, Nodes, ka, iproj, j_neig, m, Dscf_correlated = ',          &
! &                Node, Nodes, ka, iproj, j_neig, m, Dscf_correlated(ka,iproj,j_neig,m) 
!             enddo 
!           enddo 
!         enddo 
!       endif
!     enddo 
!     call die()
!!    End debugging

!     magnetic_field(1) = 0.025_dp
!     magnetic_field(2) = 0.025_dp
!     magnetic_field(3) = 0.03535533905932738_dp
     magnetic_field(1) = 0.0_dp
     magnetic_field(2) = 0.0_dp
     magnetic_field(3) = 0.0_dp
!!    For debugging
!     write(6,'(a,2i5,3f12.5)')' Node, Nodes, Magnetic field = ',    &
! &      Node, Nodes, magnetic_field(:)
!!    End debugging

!    Initialize the LSDA+U potential at every self-consistent step
     H_ldau_so         = cmplx( 0.0_dp, 0.0_dp, kind=dp )
     H_ldau_so_Hubbard = cmplx( 0.0_dp, 0.0_dp, kind=dp )
     H_ldau_so_dc      = cmplx( 0.0_dp, 0.0_dp, kind=dp )
     H_Zeeman          = cmplx( 0.0_dp, 0.0_dp, kind=dp )

!    Initialize the number of correlated electrons within this shell
     number_corr_elec_loc = 0.0_dp
     number_corr_elec     = 0.0_dp
     E_correc_dc          = 0.0_dp
     magnetic_loc(:) = cmplx(0.0_dp, 0.0_dp, kind=dp)
     magnetic(:)     = cmplx(0.0_dp, 0.0_dp, kind=dp)

!    Loop on all the atomic orbitals stored in the local node
     do io_local = 1, no_l

!      Converts the orbital index from local frame to global frame
       call LocalToGlobalOrb(io_local,Node,Nodes,io_global)

!      Find the attributes of that particular orbital
       ia = iaorb(io_global)         ! Atom to which orbital belongs
       is = isa(ia)                  ! Atomic species index
       iuo = indxuo(io_global)       ! Equivalent orbital in first unit cell
       iua = iaorb(iuo)              ! Equivalent atom in first unit cell
       iao = iphorb(io_global)       ! Orbital index within atom
       n_qn = cnfigfio( is, iao)     ! Orbital's principal quantum number
       l_qn = lofio( is, iao )       ! Orbital's angular mumentum number
       m_qn = mofio( is, iao )       ! (Real) orbital's magnetic quantum number
       zeta_qn = zetafio( is, iao )  ! 'Zeta' index of orbital
       polarized = pol( is, iao )    ! Is this a polarization orbital?
       orb_sym = symfio( is, iao )   ! Name of orbital's angular symmetry

       spp  => species(is)
       do iproj = 1, spp%n_pjldaunl  ! Number of LDA+U projectors for this
                                     !   atomic species
                                     !   not counting the "m copies"
         l_correlated = spp%pjldaunl_l(iproj)
         ldauintegrals => spp%ldau_so_integrals(iproj)
         U = spp%pjldaunl_U(iproj)
         J = spp%pjldaunl_J(iproj)

         if( l_qn .eq. l_correlated ) then   ! We need to compute the new
                                             !    Hamiltonian matrix elements
                                             !    between atomic orbitals
                                             !    that belongs to the correlated
                                             !    shell
           m = m_qn + 1 + l_correlated

!!          For debugging
!           write(6,'(a,11i5,l5,1x,a)')                                      &
! &           'ldau_so_hamil: Nodes, Node, io_l, io_g, ia, is, iuo, iua = ', &
! &                   Nodes, Node, io_local, io_global, ia, is, iuo, iua,    &
! &                   n_qn, l_qn, m_qn, polarized, orb_sym 
!!          End debugging

!          Loop on all the neighbours of the atomic orbital
           do j_neig = 1, numh(io_local) 
             ind = listhptr(io_local) + j_neig

!            Identify index of the neighbour orbital, global frame
             neig_orb = listh(ind)

!            Identify relevant information for the neighbour orbital
!            Meaning of the different variables as before
             ia_neig      = iaorb(neig_orb)  
             is_neig      = isa(ia_neig) 
             iao_neig     = iphorb(neig_orb)
             n_qn_neig    = cnfigfio( is_neig, iao_neig)  
             l_qn_neig    = lofio( is_neig, iao_neig )
             m_qn_neig    = mofio( is_neig, iao_neig )  
             zeta_qn_neig = zetafio( is_neig, iao_neig ) 
             orb_sym_neig = symfio( is_neig, iao_neig ) 

             H_Zeeman(ind,1) = -1.0_dp * magnetic_field(3) * S(ind)
             H_Zeeman(ind,2) = -1.0_dp * (-1.0_dp * magnetic_field(3)) * S(ind)
             H_Zeeman(ind,3) = -1.0_dp * (magnetic_field(1) - cmplx(0.0_dp,1.0_dp,kind=dp) * magnetic_field(2)) * S(ind)
             H_Zeeman(ind,4) = -1.0_dp * (magnetic_field(1) + cmplx(0.0_dp,1.0_dp,kind=dp) * magnetic_field(2)) * S(ind)

             if( l_qn_neig .eq. l_correlated) then
               if( ia .eq. ia_neig ) then
                 iat_correlated = index_at_correlated(ia)
                 mprime = m_qn_neig + 1 + l_correlated


!                Compute the contribution to the number of electrons in 
!                the correlated shell
!                 if( ( m .eq. mprime ) .and. (zeta_qn .eq. zeta_qn_neig) ) then
                 if( ( m .eq. mprime ) ) then
                   Dscf_cmplx_1 = cmplx(Dscf(ind,1),Dscf(ind,5), dp)
                   Dscf_cmplx_2 = cmplx(Dscf(ind,2),Dscf(ind,6), dp)
                   Dscf_cmplx_3 = cmplx(Dscf(ind,3),-Dscf(ind,4), dp)
                   Dscf_cmplx_4 = cmplx(Dscf(ind,7),Dscf(ind,8), dp)
!                   number_corr_elec = number_corr_elec +    &
! &                    real( Dscf_cmplx_1 + Dscf_cmplx_2 ) * S(ind)
                   number_corr_elec_loc = number_corr_elec_loc +    &
 &                    real( Dscf_cmplx_1 + Dscf_cmplx_2, dp ) 
                   magnetic_loc(1) = magnetic_loc(1) + (Dscf_cmplx_3 + Dscf_cmplx_4) * S(ind)
                   magnetic_loc(2) = magnetic_loc(2) + ((Dscf_cmplx_3 - Dscf_cmplx_4) * cmplx(0.0_dp,1.0_dp,kind=dp)) * S(ind)
                   magnetic_loc(3) = magnetic_loc(3) + (Dscf_cmplx_1 - Dscf_cmplx_2) * S(ind)
                 endif  ! End if m = mprime
                 


!                Compute the LSDA+U potential
!                For a given m and mprime, 
!                sum over all possible values of m2prime and m3prime
!                (i.e. loop over all possible zetas and
!                over all possible m values)
                 do iz2prime = 1, nzeta(ia,iproj)
                   do iz3prime = 1, nzeta(ia,iproj)
                     do m_qn_2neig = -l_correlated, l_correlated 
                       m2prime = m_qn_2neig + 1 + l_correlated
                       do m_qn_3neig = -l_correlated, l_correlated
                         m3prime = m_qn_3neig + 1 + l_correlated

                         icounter = (iz2prime-1)   * (2*l_correlated+1)**2 * maxnzeta + &
 &                                  (m2prime - 1)  * (2*l_correlated+1)    * maxnzeta + &
 &                                  (iz3prime - 1) * (2*l_correlated+1)               + &
 &                                  m3prime
                         direct_int =                                         &
 &    ldauintegrals%vee_4center_integrals(m,mprime,m2prime,m3prime)
                         fock_int   =                                         &
 &    ldauintegrals%vee_4center_integrals(m,m3prime,m2prime,mprime)

!!                        For debugging                     
!                         write(6,'(a,12i8)')                                   &
! &                         'm, mprime, m2prime, m3prime, ... = ',&
! &                        m_qn, m_qn_neig, m_qn_2neig, m_qn_3neig,             &
! &                        zeta_qn, zeta_qn_neig, iz2prime, iz3prime,           &
! &                        m, mprime, m2prime, m3prime                          
!!                        End debugging                     

!                        Define the real part for the (up,up) term
!                        First entry in the spin-component of the matrix element
                         H_ldau_so_Hubbard(ind,1) = H_ldau_so_Hubbard(ind,1) + one          *    &
 &                          ( ( direct_int *  Dscf_correlated(iat_correlated,iproj,icounter,2) )   +    &
 &                          ( direct_int - fock_int ) *  Dscf_correlated(iat_correlated,iproj,icounter,1) )  
                         if( (m .eq. mprime) .and. (m2prime .eq. m3prime) ) then
                           H_ldau_so_dc(ind,1) = H_ldau_so_dc(ind,1)  + one * 1.0_dp          *    &
 &                            (  - U * Dscf_correlated(iat_correlated,iproj,icounter,2) +    &
 &                               - (U - J) * Dscf_correlated(iat_correlated,iproj,icounter,1) )                         
                         endif

!                        Define the imaginary part for the (up,up) term
!                        Fifth entry in the spin-component of the matrix element
                         H_ldau_so_Hubbard(ind,1) = H_ldau_so_Hubbard(ind,1) + imag *             &
 &                          ( ( direct_int * Dscf_correlated(iat_correlated,iproj,icounter,6) ) +    &
 &                            ( direct_int - fock_int ) * Dscf_correlated(iat_correlated,iproj,icounter,5))   
                         if( (m .eq. mprime) .and. (m2prime .eq. m3prime) ) then
                           H_ldau_so_dc(ind,1) = H_ldau_so_dc(ind,1)  + imag * 1.0_dp          *    &
 &                            ( - U * Dscf_correlated(iat_correlated,iproj,icounter,6) +    &
 &                              - (U - J) * Dscf_correlated(iat_correlated,iproj,icounter,5) )                         
                         endif

!                        Define the real part for the (down,down) term
!                        Second entry in the spin-component of the matrix element
                         H_ldau_so_Hubbard(ind,2) = H_ldau_so_Hubbard(ind,2) + one           *    &
 &                          ( ( direct_int * Dscf_correlated(iat_correlated,iproj,icounter,1))+    &
 &                            ( direct_int - fock_int ) * Dscf_correlated(iat_correlated,iproj,icounter,2))                        
                         if( (m .eq. mprime) .and. (m2prime .eq. m3prime) ) then
                           H_ldau_so_dc(ind,2) = H_ldau_so_dc(ind,2)  + one *   1.0_dp         *    &
 &                          (  - U * Dscf_correlated(iat_correlated,iproj,icounter,1)+ &
 &                             -(U - J) * Dscf_correlated(iat_correlated,iproj,icounter,2))                        
                         endif

!                        Define the imaginary part for the (down,down) term
!                        Sixth entry in the spin-component of the matrix element
                         H_ldau_so_Hubbard(ind,2) = H_ldau_so_Hubbard(ind,2) + imag          *    &
 &                          ( ( direct_int * Dscf_correlated(iat_correlated,iproj,icounter,5))+    &
 &                            ( direct_int - fock_int ) * Dscf_correlated(iat_correlated,iproj,icounter,6))                         
                         if( (m .eq. mprime) .and. (m2prime .eq. m3prime) ) then
                           H_ldau_so_dc(ind,2) = H_ldau_so_dc(ind,2) + imag * 1.0_dp           *    &
 &                            ( - U * Dscf_correlated(iat_correlated,iproj,icounter,5)+    &
 &                              -(U - J) * Dscf_correlated(iat_correlated,iproj,icounter,6))                         
                         endif

!                        Define the real part for the (up,down) term
!                        Third entry in the spin-component of the matrix element
                         H_ldau_so_Hubbard(ind,3) = H_ldau_so_Hubbard(ind,3) + one           *    &
 &                          ( -fock_int * Dscf_correlated(iat_correlated,iproj,icounter,3) ) 
                         if( (m .eq. mprime) .and. (m2prime .eq. m3prime) ) then
                           H_ldau_so_dc(ind,3) = H_ldau_so_dc(ind,3) + one * 1.0_dp            *    &
 &                           (  J * Dscf_correlated(iat_correlated,iproj,icounter,3) ) 
                         endif

!                        Define the imaginary part for the (up,down) term
!                        Fourth entry in the spin-component of the matrix elemen
                         H_ldau_so_Hubbard(ind,3) = H_ldau_so_Hubbard(ind,3) + imag          *    &
 &                          ( -fock_int * (-1.0_dp)*Dscf_correlated(iat_correlated,iproj,icounter,4) ) 
                         if( (m .eq. mprime) .and. (m2prime .eq. m3prime) ) then
                           H_ldau_so_dc(ind,3) = H_ldau_so_dc(ind,3) + imag * 1.0_dp           *    &
 &                          (  J * (-1.0_dp)*Dscf_correlated(iat_correlated,iproj,icounter,4)) 
                         endif

!                        Define the real part for the (down,up) term
!                        Seventh entry in the spin-component of the matrix element
                         H_ldau_so_Hubbard(ind,4) = H_ldau_so_Hubbard(ind,4) + one           *    &
 &                          ( -fock_int * Dscf_correlated(iat_correlated,iproj,icounter,7)  ) 
                         if( (m .eq. mprime) .and. (m2prime .eq. m3prime) ) then
                            H_ldau_so_dc(ind,4) = H_ldau_so_dc(ind,4) + one * 1.0_dp            *    &
 &                          ( J * Dscf_correlated(iat_correlated,iproj,icounter,7) ) 
                         endif

!                        Define the imaginary part for the (down,up) term
!                        Eigth entry in the spin-component of the matrix element
                         H_ldau_so_Hubbard(ind,4) = H_ldau_so_Hubbard(ind,4) + imag          *    &
 &                          ( -fock_int * Dscf_correlated(iat_correlated,iproj,icounter,8) ) 
                         if( (m .eq. mprime) .and. (m2prime .eq. m3prime) ) then
                           H_ldau_so_dc(ind,4) = H_ldau_so_dc(ind,4)  + imag * 1.0_dp          *    &
 &                          ( J * Dscf_correlated(iat_correlated,iproj,icounter,8)  ) 
                         endif

!!                        For debugging
!                         write(6,'(a,10i6,16f12.5)')                                                      &
! &                         'Node, Nodes, io_local, io_global, ind, neig, m, m, m1prime, m1prime  = ',     &
! &                          Node, Nodes, io_local, io_global, ind, neig_orb, m_qn, m, m_qn_neig, mprime,  &
! &                          H_ldau_so(ind,1:4), Dscf(icounter,:)
!
!                         write(6,'(a,7(i6,2x),18f12.5)')                      &
! &                         'io, ind, neig, m1prime, m2prime, m3prime, m4prime,= ',    &
! &                         io_local, ind, neig_orb, m, mprime, m2prime, m3prime, & 
! &                         direct_int, fock_int, &
! &                         H_ldau_so(ind,1:4), Dscf(indm2m3,:)
!
!                write(6,'(a,8i5)')                                            &
! &                'ldau_so_hamil: Nodes, Node, j_neig, ind, listh, n, l, m =',&
! &                Nodes, Node, j_neig, ind, neig_orb, n_qn_neig,              &
! &                l_qn_neig, m_qn_neig
!!               End debugging
                       enddo ! Close loop on the fourth orbital m_qn_3neig 
                     enddo   ! Close loop on the third  orbital m_qn_2neig 
                   enddo     ! Close the loop on the zeta of the fourth orbital 
                             ! (iz3prime)
                 enddo       ! Close the loop on the zeta of the third orbital 
                             ! (iz2prime)
                 if( m .eq. mprime ) then
                   H_ldau_so_dc(ind,1) = H_ldau_so_dc(ind,1)  + one      *     &
 &                                0.5_dp * (U - J) 
                   H_ldau_so_dc(ind,2) = H_ldau_so_dc(ind,2)  + one      *     &
 &                                0.5_dp * (U - J) 
                 endif
               endif         ! Close condition on the same atom
             endif ! Close condition on the angular momentum of the neighbour
           enddo   ! Close loop on neighbours (j_neig)
         endif     ! Close condition on the angular momentum of the uc orbital
       enddo       ! Close loop on LDA+U projectors on the atom within the uc
     enddo         ! Close loop on the orbitals in the uc

     H_ldau_so = H_ldau_so_Hubbard + H_ldau_so_dc + H_Zeeman

!!    For debugging
!     do io_local = 1, no_l
!       call LocalToGlobalOrb(io_local,Node,Nodes,io_global)
!       do j_neig = 1, numh(io_local) 
!         ind = listhptr(io_local) + j_neig
!         neig_orb = listh(ind)
!         do m = 1, 4
!           if( abs(H_ldau_so(ind,m)) .gt. 1.d-8 ) then
!             write(6,'(a,7(i5,2x),2f12.5)')                         &
! &             'ldau_so_hamil: Node, Nodes, io_local, io_global, neig, ind, m, H_ldau_so(ind,m) ', & 
! &             Node, Nodes, io_local, io_global, neig_orb, ind, m, H_ldau_so(ind,m)
!! &             io_global, neig_orb, ind, m, H_ldau_so(ind,m), Dscf(ind,:), H_ldau_so_Hubbard(ind,m)+H_ldau_so_dc(ind,m)
!! &             io_global, neig_orb, ind, m, H_ldau_so_Hubbard(ind,m), H_ldau_so_dc(ind,m), Dscf(ind,:)
!           endif
!         enddo
!       enddo 
!     enddo  
!     call die()
!!    End debugging

      E_ldau_so   = 0.0_dp
      do ind = 1, maxnh
        Dscf_cmplx_1 = cmplx(Dscf(ind,1),Dscf(ind,5), dp)
        Dscf_cmplx_2 = cmplx(Dscf(ind,2),Dscf(ind,6), dp)
        Dscf_cmplx_3 = cmplx(Dscf(ind,3),-Dscf(ind,4), dp)
        Dscf_cmplx_4 = cmplx(Dscf(ind,7),Dscf(ind,8), dp)
!!       For debugging
!        write(6,'(a,3i5,6f12.5)')' Node, Nodes, ind, m_ldau_so: Dscf_cmplx_1: ', Node, Nodes, ind, & 
! &           Dscf_cmplx_1, H_ldau_so_Hubbard(ind,1), H_ldau_so_dc(ind,1)
!        write(6,'(a,3i5,6f12.5)')' Node, Nodes, ind, m_ldau_so: Dscf_cmplx_2: ', Node, Nodes, ind, & 
! &           Dscf_cmplx_2, H_ldau_so_Hubbard(ind,2), H_ldau_so_dc(ind,2)
!        write(6,'(a,3i5,6f12.5)')' Node, Nodes, ind, m_ldau_so: Dscf_cmplx_3: ', Node, Nodes, ind, & 
! &           Dscf_cmplx_3, H_ldau_so_Hubbard(ind,3), H_ldau_so_dc(ind,3)
!        write(6,'(a,3i5,6f12.5)')' Node, Nodes, ind, m_ldau_so: Dscf_cmplx_4: ', Node, Nodes, ind, & 
! &           Dscf_cmplx_4, H_ldau_so_Hubbard(ind,4), H_ldau_so_dc(ind,4)
!!       End debugging
        E_ldau_so = E_ldau_so +                                                &
 &        0.5_dp * ( real( H_ldau_so_Hubbard(ind,1)*conjg(Dscf_cmplx_1), dp)   +     &
 &                   real( H_ldau_so_Hubbard(ind,2)*conjg(Dscf_cmplx_2), dp)   +     &
 &                   real( H_ldau_so_Hubbard(ind,4)*conjg(Dscf_cmplx_4), dp)   +     &
 &                   real( H_ldau_so_Hubbard(ind,3)*conjg(Dscf_cmplx_3), dp) ) +     &
 &        1.0_dp * S(ind) * ( real( H_ldau_so_dc(ind,1)*Dscf_cmplx_1 , dp)     +     &
 &                   real( H_ldau_so_dc(ind,2)*Dscf_cmplx_2 , dp)        +     &
 &                   real( H_ldau_so_dc(ind,4)*Dscf_cmplx_3 , dp)        +     &
 &                   real( H_ldau_so_dc(ind,3)*Dscf_cmplx_4 , dp) ) 
      enddo

#ifdef MPI
      call MPI_AllReduce(number_corr_elec_loc, number_corr_elec, 1, &
 &                       MPI_double_precision,MPI_sum,MPI_Comm_World,MPIerror )
      call MPI_AllReduce(magnetic_loc, magnetic, 3, &
 &                       MPI_double_complex,MPI_sum,MPI_Comm_World,MPIerror )
#else
      number_corr_elec = number_corr_elec_loc
      magnetic(:)      = magnetic_loc(:)
#endif

      if (Node .eq. 0 ) then
        E_correc_dc = 0.5_dp * (U-J) * number_corr_elec**2 +  &
 &         J / 4.0_dp * ( number_corr_elec**2 - real(magnetic(1)**2 + magnetic(2)**2 + magnetic(3)**2 ) )   
      else
        E_correc_dc = 0.0_dp
      endif

!!     For debugging
!      write(6,'(a,2i5,f12.5)') ' Energy for the LDA+U+SO = ', Node, Nodes, E_ldau_so * 13.6058_dp
!      write(6,'(a,2i5,f12.5)') ' Correction for the LDA+U+SO = ', Node, Nodes, E_correc_dc * 13.6058_dp
!      write(6,'(a,2i5,f12.5)') ' Total energy for LDA+U+SO = ', Node, Nodes, (E_ldau_so + E_correc_dc) * 13.6058_dp
!      write(6,'(a,2i5,2f12.5)')' Total spin along x = ' , Node, Nodes, magnetic(1)
!      write(6,'(a,2i5,2f12.5)')' Total spin along y = ' , Node, Nodes, magnetic(2)
!      write(6,'(a,2i5,2f12.5)')' Total spin along z = ' , Node, Nodes, magnetic(3)
!      write(6,'(a,2i5,f12.5)') ' Total number of correlated electrons      = ', Node, Nodes, number_corr_elec
!!      call die()
!     End debugging

      E_ldau_so            = E_ldau_so + E_correc_dc

#ifdef MPI
   call de_alloc( auxloc,  'auxloc',  'ldau_so_hamil' )
#endif


  end subroutine ldau_so_hamil

end module m_ldau_so
