
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


  implicit none

  public ldau_so_hamil

  private


  CONTAINS 

  subroutine ldau_so_hamil( H_ldau_so )

     complex(dp), intent(inout) :: H_ldau_so(maxnh,spin%H)

     type(species_info),           pointer :: spp
     type(ldau_so_integrals_type), pointer :: ldauintegrals

     logical, save :: firstime = .true.          ! First time that this
                                      !   subroutine is called?
     integer, save :: maxldau = 0     ! Maximum number of 
                                      !   LDA+U projectors on a 
                                      !   given atom
     integer :: ka                    ! Counter for loop on atoms
     integer :: ish                   ! Counter for loop on species
     integer :: iproj                 ! Counter for loop on projector
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
     integer :: m_qn                  ! (Real) orbital's magnetic quantum number
     integer :: m_qn_neig             ! (Real) orbital's magnetic quantum number
     logical :: polarized             ! Is this a polarization orbital?
     character(len=32) :: orb_sym     ! Name of orbital's angular symmetry
     character(len=32) :: orb_sym_neig! Name of orbital's angular symmetry
     integer :: l_correlated          ! Angular momentum quantum number
                                      !   of the correlated shell

!    Initialization and allocation of matrices
     if( firstime ) then
!      Find maximum number of LDA+U projectors on a given atom
       do ka = 1, na_u
         maxldau = 0
         is  = isa(ka)
         spp  => species(is)
         maxldau = max(maxldau,spp%nprojsldau)
!        For debugging
         do iproj = 1, spp%n_pjldaunl   ! Number of LDA+U projectors
                                        !   not counting the "m copies"
           ldauintegrals => spp%ldau_so_integrals(iproj)
           l = spp%pjldaunl_l(iproj)
           nullify( ldauintegrals%index_neig_orb_corr )
           allocate( ldauintegrals%index_neig_orb_corr(-l:l,-l:l) )
           ldauintegrals%index_neig_orb_corr = 0

           write(6,'(a,5i5,2f12.5)')                                       & 
 &           'ldau_so_hamil: atom, species, iproj, n, l, U, J = ',         &
 &           ka, is, iproj, spp%pjldaunl_n(iproj), spp%pjldaunl_l(iproj),  &
 &           spp%pjldaunl_U(iproj), spp%pjldaunl_J(iproj)
           do ildauso = 0, 2*l
             write(6,'(a,i5,f12.5)')                                       & 
 &             'ldau_so_hamil: index, Slater = ',                          &
 &             ildauso, ldauintegrals%Slater_F(ildauso)         
           enddo
           do m = 1, 2*l+1
             do mprime = 1, 2*l+1
               do m2prime = 1, 2*l+1
                 do m3prime = 1, 2*l+1
            write(6,'(a,4i5,f12.5)')                                      &
 &  'ldau_so_hamil: m, mprime, m2prime, m3prime, vee_integral_cmplx = ',  &
 &            m, mprime, m2prime, m3prime,                                &
 &            ldauintegrals%vee_4center_integrals(m, mprime, m2prime, m3prime)

                 enddo 
               enddo 
             enddo 
           enddo 
         enddo 
         do iproj = 1, spp%nprojsldau   ! Number of LDA+U projectors
                                        !   counting the "m copies"
           write(6,'(a,7i5)')                                              & 
 &           'ldau_so_hamil: atom, species, iproj, index, n, l, m = ',     &
 &           ka, is, iproj,                                                &
 &           spp%pjldau_index(iproj), spp%pjldau_n(iproj),                 &
 &           spp%pjldau_l(iproj), spp%pjldau_m(iproj)
         enddo
!        End debugging
       enddo

!      Loop on all the atomic orbitals stored in the local node
       do io_local = 1, no_l

!        Converts the orbital index from local frame to global frame
         call LocalToGlobalOrb(io_local,Node,Nodes,io_global)

!        Find the attributes of that particular orbital
         ia = iaorb(io_global)        ! Atom to which orbital belongs
         is = isa(ia)                 ! Atomic species index
         iuo = indxuo(io_global)      ! Equivalent orbital in first unit cell
         iua = iaorb(iuo)             ! Equivalent atom in first unit cell
         iao = iphorb(io_global)      ! Orbital index within atom
         n_qn = cnfigfio( is, iao)    ! Orbital's principal quantum number
         l_qn = lofio( is, iao )      ! Orbital's angular mumentum number
         m_qn = mofio( is, iao )      ! (Real) orbital's magnetic quantum number
         polarized = pol( is, iao )   ! Is this a polarization orbital?
         orb_sym = symfio( is, iao )  ! Name of orbital's angular symmetry

         spp  => species(is)
         do iproj = 1, spp%n_pjldaunl  ! Number of LDA+U projectors for this
                                       !   atomic species
           ldauintegrals => spp%ldau_so_integrals(iproj)
           l_correlated = spp%pjldaunl_l(iproj)

           if( l_qn .eq. l_correlated ) then ! We need to compute the new
                                             !    Hamiltonian matrix elements
                                             !    between atomic orbitals
                                             !    that belongs to the correlated
                                             !    shell

!            For debugging
             write(6,'(a,11i5,l5,1x,a10)')                                    &
 &             'ldau_so_hamil: Nodes, Node, io_l, io_g, ia, is, iuo, iua = ', &
 &                   Nodes, Node, io_local, io_global, ia, is, iuo, iua,      &
 &                   n_qn, l_qn, m_qn, polarized, orb_sym 
!            End debugging

!            Loop on all the neighbours of the atomic orbital
             do j_neig = 1, numh(io_local) 
               ind = listhptr(io_local) + j_neig

!              Identify index of the neighbour orbital, global frame
               neig_orb = listh(ind)

!              Identify relevant information for the neighbour orbital
!              Meaning of the different variables as before
               ia_neig      = iaorb(neig_orb)  
               is_neig      = isa(ia_neig) 
               iao_neig     = iphorb(neig_orb)
               n_qn_neig    = cnfigfio( is_neig, iao_neig)  
               l_qn_neig    = lofio( is_neig, iao_neig )
               m_qn_neig    = mofio( is_neig, iao_neig )  
               orb_sym_neig = symfio( is_neig, iao_neig ) 

               if( l_qn_neig .eq. l_correlated) then
                 if( ia .eq. ia_neig ) then 
                   ldauintegrals%index_neig_orb_corr(m_qn,m_qn_neig) = ind
!                  For debugging
                   write(6,'(a,8i5)')                                          &
 &                 'ldau_so_hamil: Nodes, Node, j_neig, ind, listh, n, l, m =',&
 &                 Nodes, Node, j_neig, ind, neig_orb, n_qn_neig,              &
 &                 l_qn_neig, m_qn_neig
!                  End debugging
                 endif
               endif ! Close condition on the angular momentum of the neighbour
             enddo   ! Close loop on neighbours (j_neig)
           endif     ! Close condition on the angular momentum of the uc orbital
         enddo       ! Close loop on LDA+U projectors on the atom within the uc
       enddo         ! Close loop on the orbitals in the uc

     endif

     firstime = .false.

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
       polarized = pol( is, iao )    ! Is this a polarization orbital?
       orb_sym = symfio( is, iao )   ! Name of orbital's angular symmetry

       spp  => species(is)
       do iproj = 1, spp%n_pjldaunl  ! Number of LDA+U projectors for this
                                     !   atomic species
         l_correlated = spp%pjldaunl_l(iproj)

         if( l_qn .eq. l_correlated ) then   ! We need to compute the new
                                             !    Hamiltonian matrix elements
                                             !    between atomic orbitals
                                             !    that belongs to the correlated
                                             !    shell
           m = m_qn + 1 + l_correlated

!!          For debugging
!           write(6,'(a,11i5,l5,1x,a)')                                     &
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
             orb_sym_neig = symfio( is_neig, iao_neig ) 

             if( l_qn_neig .eq. l_correlated) then
               if( ia .eq. ia_neig ) then
                mprime = m_qn_neig + 1 + l_correlated
!!              For debugging
                write(6,'(a,5i5,8f12.5)')                                      &
 &                'm, m, m1prime, m1prime, ind = ',                            &
 &                m_qn, m, m_qn_neig, mprime,                                  &
 &                ldauintegrals%index_neig_orb_corr(m_qn,m_qn_neig),           &
 &                Dscf(ind,:)
!               write(6,'(a,8i5)')                                             &
! &               'ldau_so_hamil: Nodes, Node, j_neig, ind, listh, n, l, m =', &
! &               Nodes, Node, j_neig, ind, neig_orb, n_qn_neig,               &
! &               l_qn_neig, m_qn_neig
!!              End debugging
               endif
             endif ! Close condition on the angular momentum of the neighbour
           enddo   ! Close loop on neighbours (j_neig)
         endif     ! Close condition on the angular momentum of the uc orbital
       enddo       ! Close loop on LDA+U projectors on the atom within the uc
     enddo         ! Close loop on the orbitals in the uc

  end subroutine ldau_so_hamil

end module m_ldau_so
