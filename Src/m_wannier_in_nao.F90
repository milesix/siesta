!
! Copyright (C) 1996-2016       The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!

!> \brief General purpose of the m_wannier_in_nao module:
!! Here we write the expansion of a Wannier function, 
!! output of the WANNIER-90 code, in the numerical atomic orbital (NAO)
!! basis set of SIESTA.
!! 
!! The analytical expressions can be found in
!!
!! <https://personales.unican.es/junqueraj/JavierJunquera_files/Notes/Wannier/wannier_in_nao.pdf>
!!

module m_wannier_in_nao

  use precision,      only: dp                ! Real double precision type

  implicit none

!
! Variables related with the coefficients of the wavefunctions and
! eigenvalues at the Wannier90 k-point mesh
!
  complex(dp), pointer :: coeffs_wan_nao(:,:) => null()
                                             ! Coefficients of the 
                                             !   wannier functions in a basis
                                             !   of NAO
                                             !   First  index: wannier function
                                             !   Second index: NAO in the 
                                             !       supercell

  CONTAINS

  subroutine wannier_in_nao
    use siesta_geom,    only: na_s           ! Number of atoms in the supercell
    use siesta_geom,    only: xa             ! Atomic positions for all the
                                             !   atoms in the supercell
                                             !   In Bohrs
    use siesta_geom,    only: isa            ! Species index of each atom
    use atomlist,       only: no_s           ! Number of atomic orbitals in the
                                             !   supercell
    use atomlist,       only: no_u           ! Number of atomic orbitals in the
                                             !   unit cell
    use atomlist,       only: indxuo         ! Index of equivalent orbital in 
                                             !   the unit cell
    use atomlist,       only: iaorb          ! Atomic index of each orbital
    use atomlist,       only: iphorb         ! Orbital index of each orbital 
                                             !   in its atom
    use atmfuncs,       only: labelfis       ! Atomic label
    use atmfuncs,       only: cnfigfio       ! Principal quantum number of the
                                             !   shell to what the orbital 
                                             !   belongs (for polarization 
                                             !   orbitals the quantum number 
                                             !   corresponds to the shell which
                                             !   is polarized by the orbital io)
    use atmfuncs,       only: symfio         ! Symmetry of the orbital
    use w90_parameters, only: num_kpts       ! Total number of k-points that
                                             !   will be used within
                                             !   Wannier90
    use w90_parameters, only: kpt_latt       ! Coordinates of the k-points
                                             !   that will be used within
                                             !   Wannier90
                                             !   In fractional units, i. e.
                                             !   in units of the reciprocal
                                             !   lattice vectors
                                             !   First  index: component
                                             !   Second index: k-vector
    use w90_parameters, only: recip_lattice  ! Reciprocal lattice vectors
                                             !   as used inside
                                             !   Wannier90
                                             !   The factor 2.0 * pi
                                             !   is included
                                             !   BE CAREFUL: In WANNIER90,the
                                             !   order of the indices is
                                             !   inverted with respect to
                                             !   SIESTA
                                             !   First  index: vector
                                             !   Second index: component
                                             !   In Angstroms^-1
    use w90_parameters, only: num_proj       ! Number of projections
    use w90_parameters, only: u_matrix       ! Unitary rotations from the 
                                             !   optimal subspace to the
                                             !   optimally smooth states.
    use w90_parameters, only: gamma_only     ! Only the gamma point will be
                                             !   used within Wannier90?
    use w90_constants,  only: bohr_angstrom_internal  
                                             ! Conversion factor from 
                                             !   Bohr to Ang
    use m_switch_local_projection, only: coeffs
                                             ! Coefficients of the wavefunctions
                                             !   First  index: orbital
                                             !   Second index: band
                                             !   Third  index: k-point
    use m_switch_local_projection, only: seedname
                                             ! Name of the file where the 
                                             !   WANNIER90 code
                                             !   reads or dumps the
                                             !   information.
    use writewave,      only: writew         ! Subroutine to dump the 
                                             !   coefficients of the 
                                             !   wavefunction
                                             !   in a .WFSX file
    use writewave,      only: setup_wfs_list ! Subroutine to setup the
                                             !   list of wave functions 
                                             !   (in this case, Wannier func.)
                                             !   to be plotted
    use writewave,      only: wfs_filename   ! Name of the file where the 
                                             !   coefficients of the Wanniers
                                             !   will be dumped
    use writewave,      only: wwf            ! Write wave functions?

!
! Allocation/Deallocation routines
!
    use m_switch_local_projection, only: coeffs
                                             ! Coefficients of the wavefunctions
    use alloc,                     only: re_alloc  
                                             ! Reallocation routines
    use alloc,                     only: de_alloc  
                                             ! Deallocation routines

! 
! Internal variables 
! 
    integer  :: iu                           ! Logical unit
    integer  :: ik                           ! Counter for loop on k-points
    integer  :: iorb                         ! Counter for loop on atomic 
                                             !   orbitals
    integer  :: iproj                        ! Counter for loop on projectors
    integer  :: iprojn                       ! Counter for loop on projectors
    integer  :: iprojm                       ! Counter for loop on projectors
    integer  :: index_orbital_unit           ! Equivalent index in the unit cell
                                             !   of an orbital in the supercell
    integer  :: index_atom                   ! Index of the atom to which an
                                             !   atomic orbital belongs to
    integer  :: nk                           ! Number of k-points written in
                                             !   the .WFSX file
                                             !   Since the number of Wanniers 
                                             !   functions to be plotted do not 
                                             !   depend in k-points,
                                             !   it is set to 1
    integer  :: maxspn                       ! Maximum number of spin components
    logical  :: gamma                        ! Whether only the Gamma-point is 
                                             !    sampled.
                                             !    Since the dependency of the 
                                             !    phase is already included in 
                                             !    the coefficients, we set it
                                             !    up to .true.
    real(dp) :: kpoint(3)                    ! Coordinates of the k-point
    real(dp) :: kxmu                         ! Dot product of the k-point and
                                             !   the position of the atom in the
                                             !   supercell
    real(dp) :: ckxmu                        ! Cosine of kxmu
    real(dp) :: skxmu                        ! Sine of kxmu
    real(dp) :: kdummy(3)                    ! Dummy variable for the k-points
    real(dp), pointer :: psi(:,:,:)          ! Dummy variable to store the 
                                             !   coefficients of the Wanniers
                                             !   in a basis of atomic orbitals
                                             !   to call writew
    real(dp), target :: aux(2,no_s*5)        ! Dummy variable that will play
                                             !   the role of the eigenvalues
                                             !   in the call to writew
             

!   Allocate the array where the coefficients of the Wannier functions
!   in a basis of Numerical Atomic Orbitals will be stored
    nullify( coeffs_wan_nao )
    call re_alloc( coeffs_wan_nao,                                  &
 &                 1, num_proj,                                     &
 &                 1, no_s,                                         &
 &                 name='coeffs_wan_nao', routine='wannier_in_nao')
    coeffs_wan_nao = cmplx(0.0_dp,0.0_dp,kind=dp)

!   Allocate the array where the coefficients of the Wannier functions
!   in a basis of Numerical Atomic Orbitals will be stored
    nullify( psi )
    call re_alloc( psi,                                             &
 &                 1, 2,                                            &
 &                 1, no_s,                                         &
 &                 1, num_proj,                                     &
 &                 name='psi', routine='wannier_in_nao')
    psi = 0.0_dp


!    For debugging
!    do ik = 1, num_kpts
!      do iproj = 1, num_proj
!        do iorb = 1, no_u
!          write(6,'(a,3i5,4f12.5)')' ik, iproj, iorb, coeff = ', &
! &          ik, iproj, iorb, coeffs(iorb,iproj,ik),              &
! &          u_matrix(iproj,iproj,ik)
!        enddo 
!      enddo 
!    enddo 
!    End debugging

    do iprojn = 1, num_proj
      do ik = 1, num_kpts
!       Compute the coordinates of the k-point (in Ang^-1)
        kpoint(:) = kpt_latt(1,ik) * recip_lattice(1,:) +           &
 &                  kpt_latt(2,ik) * recip_lattice(2,:) +           &
 &                  kpt_latt(3,ik) * recip_lattice(3,:) 
!       Transform the coordinates of the k-point to Bohr^-1
        kpoint(:) = kpoint(:) * bohr_angstrom_internal
!!       For debugging        
!        write(6,'(a,2i5,3f12.5)')' iprojn, ik, kpoint = ',          &
! &                                 iprojn, ik, kpoint(:)
!!       End debugging        
        do iprojm = 1, num_proj
          do iorb = 1, no_s
            index_orbital_unit = indxuo(iorb)
            index_atom         = iaorb(iorb)
!!          For debugging        
!           write(6,'(a,5i5,6f12.5)')' iprojn, ik, iorb, kpoint = ',          &
! &                                    iprojn, ik, iorb, index_orbital_unit,  &
! &                                    index_atom, xa(:,index_atom), kpoint(:)
!!         End debugging        
   
           kxmu = kpoint(1) * xa(1,index_atom) +                        &
 &                kpoint(2) * xa(2,index_atom) +                        &
 &                kpoint(3) * xa(3,index_atom) 
           ckxmu = dcos(kxmu)
           skxmu = dsin(kxmu)

           coeffs_wan_nao(iprojn,iorb) = coeffs_wan_nao(iprojn,iorb) +  &
 &            u_matrix(iprojm,iprojn,ik)            *                   &
 &            coeffs(index_orbital_unit,iprojm,ik)  *                   &
 &            cmplx(ckxmu,skxmu,kind=dp) 
         
          enddo ! End loop on orbitals in the supercell
        enddo   ! End loop on projections m
      enddo     ! End loop on k-points
    enddo       ! End loop on projections n

    coeffs_wan_nao = coeffs_wan_nao / num_kpts

!!   For debugging
!    do iorb = 1, no_s
!      write(6,'(a,3i5,3f12.5)')'iorb, indxuo, iaorb = ',  &
! &      iorb, indxuo(iorb), iaorb(iorb), xa(:,iaorb(iorb))
!    enddo
     do iproj = 1, num_proj
       do iorb = 1, no_s
!         if( indxuo(iorb) .eq. 15) then
         if( aimag(coeffs_wan_nao(iproj,iorb)) .gt. 1.d-5 ) then
         write(6,'(a,2i5,5f12.5)') ' iproj, iorb, coeffs_wan_nao = ', &
 &         iproj, iorb, coeffs_wan_nao(iproj,iorb), xa(:,iaorb(iorb))
         endif
       enddo 
     enddo 
!!   End debugging

!   Set up the variables to call writew
    do iproj = 1, num_proj
      do iorb = 1, no_s
        psi(1,iorb,iproj) = real(coeffs_wan_nao(iproj,iorb))
        psi(2,iorb,iproj) = aimag(coeffs_wan_nao(iproj,iorb))
!       For debugging
!        if( psi(2,iorb,iproj) .gt. 1.d-5) then
        write(6,'(a,2i5,2f12.5)')' iproj, iorb, psi = ', &
 &        iproj, iorb, psi(1,iorb,iproj), psi(2,iorb,iproj)
!        endif
!       End debugging
      enddo 
    enddo
    aux          = 0.0_dp
    kdummy       = 0.0_dp
    wfs_filename = trim(seedname)//".WFSX"
    wwf          = .true. 
    nk           = 1
    gamma        = .false.
    maxspn       = 1

! 
!   Open the WFSX file and print the header of the file
!   with information of the atomic orbitals
!   This was done in the subroutine wwave when the coefficients of the 
!   wave functions at particular k-points are required
!
    call io_assign( iu )
    open(iu, file=wfs_filename,form="unformatted",status='unknown')

    write(iu) nk, gamma
    write(iu) maxspn
    write(iu) no_s
    write(iu) (iaorb(iorb),labelfis(isa(iaorb(iorb))),                    &
 &            iphorb(iorb), cnfigfio(isa(iaorb(iorb)),iphorb(iorb)),      &
 &            symfio(isa(iaorb(iorb)),iphorb(iorb)), iorb=1,no_s)
    call io_close( iu )

!!   For debugging
!    write(6,'(a,i5,l5)')' nk, gamma = ', nk, gamma
!    write(6,'(a,i5)')   ' maxspn = ', maxspn
!    write(6,'(a,i5)')   ' no_s   = ', no_s
!    do iorb = 1, no_s
!      write(6,*) iorb, iaorb(iorb),labelfis(isa(iaorb(iorb))),            &
! &               iphorb(iorb), cnfigfio(isa(iaorb(iorb)),iphorb(iorb)),   &
! &               symfio(isa(iaorb(iorb)),iphorb(iorb))
!    enddo 
!!   End debugging

    call setup_wfs_list(nk,no_s,1,num_proj,.false.,.false.)

    call writew( no_s, num_proj, 1, kdummy, 1, &
 &               aux, psi, gamma )

    
  end subroutine wannier_in_nao

end module m_wannier_in_nao
