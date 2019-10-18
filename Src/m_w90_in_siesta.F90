 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!

!> \brief General purpose of the m_w90_in_siesta module:
!! interface between SIESTA and WANNIER90 code.
!! 
!! In this module we perform a Wannier transformation of the Bloch functions
!! corresponding to a given manifold of bands, 
!! following the recipe given in \cite Marzari-97.
!! In order to compute the transformation, SIESTA uses the low-level modules
!! of the WANNIER90 code(\cite Wannier90), version 3.0.0, as building blocks
!! A complete description of the theory behind the maximally localized
!! Wannier functions and its applications can be found in Ref. \cite Marzari-12
!!
!! The tutorial with a full explanation of the new labels in the input file
!! can be found in
!!
!! <https://personales.unican.es/junqueraj/JavierJunquera_files/Metodos/Wannier/Exercise-Wannier90-within-siesta.pdf>
!!
!! Example: one of the Wannier functions obtained after a transformation of the
!! top of the valence band manifold of SrTiO3 (O-2p in character)
!! \image html SrTiO3.manifold.1_00001.png
module m_w90_in_siesta
#ifdef HAVE_WANNIER90


  use precision,      only: dp            ! Real double precision type
  use siesta_options, only: n_wannier_manifolds     
                                          ! Number of bands manifolds 
                                          !  that will be considered
                                          !  for Wannier transformation
  use siesta_options, only: r_between_manifolds     
                                          ! Will the position operator matrix
                                          !  elements be computed between bands
                                          !  of different manifolds
  use siesta_options, only: w90_in_siesta_compute_unk
                                          ! Shall SIESTA compute the files 
                                          !   which contains the periodic part 
                                          !   of the wave function in the 
                                          !   unit cell on a grid
  use parallel,       only: Node          ! Local processor number
  use parallel,       only: Nodes         ! Total number of processors in a 
                                          !  parallel run
  use parallel,       only: IONode        ! Node for Input/output
  use w90_in_siesta_types,   only: w90_in_manifold_t    
  use w90_in_siesta_types,   only: manifold_bands_w90_in
                                       ! Variable where the initial
                                       !   and final band of each
                                       !   manifold are stored
  use w90_in_siesta_types,   only: compute_chempotwann
                                       ! Compute the Hamiltonian matrix
                                       !   elements between NAO
                                       !   if a chemical potential in
                                       !   a Wannier is applied
  use w90_in_siesta_types,   only: first_chempotwann 
                                       !  First time the calculation of the
                                       !   matrix elements of the Hamiltonian
                                       !   with the chemical potential of the 
                                       !   Wanniers is called?
  use w90_in_siesta_types,   only: chempotwann_val
                                       ! Chemical potential
                                       !   applied to shift the energy of a
                                       !   the matrix elements in real space
                                       !   associated with a given
                                       !   Wannier function
  use w90_in_siesta_types,   only: numh_man_proj
                                       ! Number of projections that will be
                                       !   handled in the local node
                                       !   for a given manifolds
                                       !   Dimension: number of manifolds
  use w90_in_siesta_types,   only: listhptr_man_proj
                                       ! Index pointer to listh_man_proj such 
                                       ! listh_man_proj(listhptr_man_proj(1)+1)
                                       !   is the first projector of the first
                                       !   manifold handled by the local node 
                                       ! listh_man_proj(listhptr_man_proj(io)+1)
                                       !   is thus the first projector of 
                                       !   of manifold 'io' while 
                                       ! listh_man_proj(listhptr_man_proj(io) +
                                       !                numh_man_proj(io)) 
                                       !   is the last projectors of manifold 
                                       !   'io'.
                                       ! Dimension: number of manifolds
  use w90_in_siesta_types,   only: listh_man_proj
                                       ! The column indices for the projectors
                                       !   of all the manifolds handled by 
                                       !   the local node
  use w90_in_siesta_types,       only: coeffs_wan_nao
                                       ! Coefficients of the
                                       !   Wannier functions in a basis
                                       !   of NAO
                                       !   First  index: Index of the
                                       !       manifold and projector,
                                       !       handled by numh_man_proj,
                                       !       listhptr_man_proj, and
                                       !       listh_man_proj lists
                                       !   Second index: NAO in the
                                       !       supercell
                                       !   Third index: Spin component
  use atomlist,           only: no_u   ! Number of orbitals in unit cell
  use atomlist,           only: no_u   ! Number of orbitals in unit cell
                                       ! NOTE: When running in parallel,
                                       !   this is core independent
  use atomlist,           only: no_s   ! Number of atomic orbitals in the
                                       !   supercell
  use sys,                only : die   ! Termination routine
  use fdf

!
! Allocation/Deallocation routines
!
  use alloc,              only: re_alloc  ! Reallocation routines
  use alloc,              only: de_alloc  ! Deallocation routines

#ifdef MPI
  use mpi_siesta
  use parallelsubs,       only: LocalToGlobalOrb
#endif


  implicit none


! Routines
  public :: setup_w90_in_siesta
  public :: compute_matrices
  public :: compute_wannier
  public :: deallocate_wannier

  private

  CONTAINS

!> \brief General purpose of the subroutine setup_w90_in_siesta:
!! set up the interface between SIESTA and WANNIER90
!!
!! Within this subroutine:
!!
!! 1. Include in the SystemName.bib file a citation to the papers where
!!    the implementation is based. Here:
!!    the original paper by Marzari and Vanderbilt,
!!    the last paper on the implementation of WANNIER90,
!!    and the pre-print by Íñiguez and Yildirim in the arxiv where the 
!!    Lowdin orthogonalization is explained.
!!
!! 2. Read the WannierProjections block in the fdf 
!!    and set up the general information for every manifold,
!!    including initial (initial_band) and final (final_band) bands,
!!    number of bands that will be orthogonalized (numbands_w90_in),
!!    the orbitals that will be used as initial guess in the projections
!!    during the wannierization (orbital_indices),
!!    number of iterations that will be performed in the Wannierization 
!!    procedure (num_iter),
!!    the energy windows for the disentanglement,
!!    and other output options (if the wannier functions will be plotted,
!!    or the Fermi surface computed).
!!    This is done in the subroutine read_w90_in_siesta_specs. 
!!
!! 3. Set up the indexing of the bands that will be excluded from the 
!!    wannierization (isexcluded),
!!    the orbitals that will be used as the localized projector functions
!!    (orbexcluded),
!!    and the new indexing of the orbitals within a manifold (orb_in_manifold).
!!    This is done in the subroutine set_excluded_bands_w90_in.
!!
!! 4. Generate the trial localized functions from the atomic orbitals in the
!!    basis set of SIESTA.
!!    This is done in the subroutine define_trial_orbitals.
!!
!! 5. The derived variable manifold_bands_w90_in (where all the previous
!!    information is stored) is broadcast
!!    to the rest of the Nodes in the subroutine broadcast_w90_in_siesta.
!!
!! 6. Determine the number of bands that will be treated per node in a 
!!    parallel run.
!!
!! 7. Populate the lists required to handle the matrix of coefficients of
!!    a Wannier functions in terms of the numerical atomic orbitals of the
!!    supercells. These matrices are very similar to numh, listhptr, and 
!!    listh to handle the Hamiltonian and Overlap sparse matrices
!!
!! 8. Read all the information regarding the 
!!    k-point sampling (block kMeshforWannier), 
!!    and the neighbours for all the k-points in the BZ,
!!    required for the wannierization.
!!    This is done in the subroutine read_kpoints_wannier
!!    by all the nodes simultaneously, so no need for broadcasting
!!    these variables.
!!
  subroutine setup_w90_in_siesta()

    use m_cite,        only: add_citation          ! Subroutine used to cite 
                                                   !   the proper papers where 
                                                   !   the implementation is
                                                   !   based on
    use m_spin,        only: spin                  ! Spin configuration
                                                   !   for SIESTA
#ifdef MPI
    use parallelsubs,  only: set_blocksizedefault  ! Subroutine to find
                                                   !   a sensible default value
                                                   !   for the blocksize default
    use parallel,      only: BlockSize             ! Blocking factor used to
                                                   !   divide up the arrays over
                                                   !   the processes for the
                                                   !   Scalapack routines.
#endif

!
! Internal parameters
!

    integer :: index_manifold
    integer :: iproj_local
    integer :: iproj_global
    integer :: index
    integer :: blocksizeincbands_tmp 
    integer :: numincbands_tmp
    integer :: nincbands_loc_tmp
    integer :: numproj
    integer :: maxnh_man_proj

#ifdef MPI
    integer, external :: numroc                    ! Scalapack routine for 
                                                   !  block-cyclic distributions
    integer :: MPIError
    integer, allocatable :: blocksizeprojectors(:)
    integer :: blocksize_tmp 
#endif


    if( IONode ) then

!     Add a citation ...
!     ... to the Marzari and Vanderbilt paper on 
!     Maximally Localized Wannier Functions
      call add_citation("10.1103/PhysRevB.56.12847")
!     ... to the paper describing the WANNIER90 code
      call add_citation("10.1016/j.cpc.2014.05.003")
!     ... to the paper describing the Löwdin orthogonalization
      call add_citation("arXiv/cond-mat/0407677")

!
!     Read and set up the general information for every manifold,
!     including initial (initial_band) and final (final_band) bands,
!     number of bands that will be orthogonalized (numbands_w90_in),
!     number of iterations that will be performed in the Wannierization 
!     procedure (num_iter),
!     the orbitals that will be used in the wannierization (orbital_indices),
!     the energy windows for the disentanglement,
!     and other output options (if the wannier functions will be plotted,
!     or the Fermi surface computed).
!

      call read_w90_in_siesta_specs

!
!     Set up the indexing of the bands that will be excluded from the 
!     wannierization (isexcluded),
!     the orbitals that will be used as the localized projector functions
!     (orbexcluded) ,
!     and the new indexing of the orbitals within a manifold (orb_in_manifold)
!
      do index_manifold = 1, n_wannier_manifolds
        call set_excluded_bands_w90_in( index_manifold )
      enddo 

!     Generate the trial localized functions from the atomic orbitals in the
!     basis set of SIESTA

      do index_manifold = 1, n_wannier_manifolds
        call define_trial_orbitals( index_manifold )
      enddo

    endif   ! close if( IONode )

!
!   Up to know, all the information is known from the IONode,
!   Here we broadcast the information to all the other nodes.
!
    call broadcast_w90_in_siesta

!
!   Determine the number of bands that will be treated per node in a 
!   parallel run
!
    do index_manifold = 1, n_wannier_manifolds
       numincbands_tmp = manifold_bands_w90_in(index_manifold)%number_of_bands

!      Allocate memory related with the coefficients of the wavefunctions
#ifdef MPI
!      Find the number of included bands for Wannierization that will be stored
!      per node. Use a block-cyclic distribution of nincbands over Nodes.

       call set_blocksizedefault( Nodes, numincbands_tmp,                  &
 &                                blocksizeincbands_tmp )

       nincbands_loc_tmp = numroc( numincbands_tmp,                        &
 &                           blocksizeincbands_tmp, Node, 0, Nodes )
!!      For debugging
!       write(6,*) 
!       write(6,'(a,l5)')                                                   & 
! &       'setup_w90_in_siesta: r_between_manifolds = ',                    &
! &       r_between_manifolds
!       write(6,'(a,6i5)')                                                  & 
! &       'setup_w90_in_siesta: index_manifold, Node, Nodes, nincbands_loc_tmp, numincbands_tmp, Blocksize = ',  &
! &       index_manifold, Node, Nodes, nincbands_loc_tmp, numincbands_tmp,                  &
! &       blocksizeincbands_tmp
!!      End debugging
#else
       nincbands_loc_tmp = numincbands_tmp
#endif

       manifold_bands_w90_in(index_manifold)%nincbands_loc_w90_in =        &
 &        nincbands_loc_tmp
       manifold_bands_w90_in(index_manifold)%blocksizeincbands_w90_in =    &
 &        blocksizeincbands_tmp
    enddo  

!   Compute the lists that will be required to handle the coefficients of 
!   the Wannier functions as an expansion of Numerical Atomic Orbitals

    nullify( numh_man_proj )
    call re_alloc( numh_man_proj,  1, n_wannier_manifolds, &
 &                 'numh_man_proj', 'setup_w90_in_siesta' )
    numh_man_proj(:) = 0

    nullify( listhptr_man_proj )
    call re_alloc( listhptr_man_proj,  1, n_wannier_manifolds, &
 &                'listhptr_man_proj', 'setup_w90_in_siesta' )
    listhptr_man_proj(:) = 0

#ifdef MPI
    allocate(blocksizeprojectors(n_wannier_manifolds))
#endif

    maxnh_man_proj = 0
    do index_manifold = 1, n_wannier_manifolds
      numproj = manifold_bands_w90_in(index_manifold)%numbands_w90_in
#ifdef MPI
!     Find the number of projectors that will be stored
!     per node. Use a block-cyclic distribution of numproj over Nodes.
      call set_blocksizedefault( Nodes, numproj,                  &
 &                               blocksizeprojectors(index_manifold) )
      numh_man_proj(index_manifold) = numroc( numproj,           &
 &               blocksizeprojectors(index_manifold), Node, 0, Nodes )
#else
      numh_man_proj(index_manifold) = numproj
#endif
      maxnh_man_proj = maxnh_man_proj + numh_man_proj(index_manifold)
    enddo  

    listhptr_man_proj(1) = 0
    do index_manifold = 2, n_wannier_manifolds
      listhptr_man_proj(index_manifold)=listhptr_man_proj(index_manifold-1) +&
 &       numh_man_proj(index_manifold-1)
    enddo  

    call re_alloc( listh_man_proj,  1, maxnh_man_proj, &
 &                'listh_man_proj', 'setup_w90_in_siesta' )

    index = 0
#ifdef MPI
    blocksize_tmp = BlockSize
#endif
    do index_manifold = 1, n_wannier_manifolds
      do iproj_local = 1, numh_man_proj( index_manifold )
        index = index + 1
#ifdef MPI
        BlockSize = blocksizeprojectors(index_manifold)
        call LocalToGlobalOrb(iproj_local, Node, Nodes, iproj_global)
#else
        iproj_global = iproj_local
#endif
        listh_man_proj(index) = iproj_global
      enddo
    enddo 
#ifdef MPI
    BlockSize = blocksize_tmp
#endif

!!   For debugging
!    do index_manifold = 1, n_wannier_manifolds
!       write(6,'(a,4i5)') &
! &       'setup_w90_in_siesta: index_manifold, Node,numproj,numh_man_proj=',&
! &       index_manifold, Node, numproj, numh_man_proj(index_manifold) 
!       write(6,'(a,4i5)') &
! &       'setup_w90_in_siesta: index_manifold, Node, numproj, listhptr_man_proj =',&
! &       index_manifold, Node, numproj, listhptr_man_proj(index_manifold) 
!    enddo  
!    write(6,'(a,3i5)') &
! &    'setup_w90_in_siesta: Node, maxnh_man_proj =',               &
! &                          Node, maxnh_man_proj 
!    do index = 1, maxnh_man_proj
!      write(6,'(a,3i5)')                                             &
!        'setup_w90_in_siesta: Node, index, listh_man_proj  =',       &
! &                            Node, index, listh_man_proj(index)
!    enddo 
!!   End debugging

!   Allocate the array where the coefficients of the Wannier functions
!   in a basis of Numerical Atomic Orbitals will be stored
    nullify( coeffs_wan_nao )
    call re_alloc( coeffs_wan_nao,                                  &
 &                 1, maxnh_man_proj,                               &
 &                 1, no_s,                                         &
 &                 1, spin%H,                                       &
 &                 name='coeffs_wan_nao', routine='wannier_in_nao')
    coeffs_wan_nao = cmplx(0.0_dp,0.0_dp,kind=dp)

#ifdef MPI
    deallocate(blocksizeprojectors)
#endif

!! For debugging
!#ifdef MPI
!    call MPI_barrier(MPI_Comm_world,MPIError)
!#endif
!    call die()
!! End debugging

!
!   Read all the information regarding the k-point sampling, 
!   and the neighbours for all the k-points in the BZ
!   required for the wannierization.
!   This is done by all the nodes simultaneously, so no need for broadcasting
!   these variables
!
    if(IOnode) call read_kpoints_wannier


  end subroutine setup_w90_in_siesta

!
! ------------------------------------------------------------------------------
!


!> \brief  General purpouse of the subroutine read_w90_in_siesta_specs:
!!         read all the info in the fdf file related with the
!!         band manifolds that will be treated in the 
!!         Wannier transformation.
!!
!!         This information is contained in the 
!!         block WannierProjections.
!!         The derived variable manifold_bands_w90_in
!!         is populated. 
!!         The type of this variable is defined in the 
!!         module w90_in_siesta_types,
!!         in the file w90_in_types.f.
!!
!!         This subroutine is used only by the Node responsible for the
!!         input/output.

   subroutine read_w90_in_siesta_specs
!
!   Processes the information in an fdf file
!   regarding the bands that will enter into the orthonormalization procedure
!   and the atomic orbitals that will be orthonormalized.
!
    use fdf
    use files,    only: slabel          ! Short system label,
                                        !   used to generate file names
    use m_spin,   only: spin            ! Spin configuration
                                        !   for SIESTA

!   Internal variables
    integer :: index_manifold           ! Counter for the number of manifolds
    integer :: iorb                     ! Counter for the number of atomic orb.
    integer :: index_wannier            ! Variable to identify a given Wannier
    real(dp):: factor                   ! Conversion factor for energy units
    character(len=30) :: string_num     ! Check whether the number of 
                                        !    iterations is in the input      
    character(len=30) :: string_plot    ! 
    character(len=30) :: string_plot_f  !
    character(len=30) :: string_write_h !
    character(len=30) :: string_write_tb!

    type(block_fdf)            :: bfdf
    type(parsed_line), pointer :: pline

!   Allocate the pointer where all the data required for every manifold
!   will be stored
    allocate(manifold_bands_w90_in(n_wannier_manifolds))
    
!   Read the WannierProjections block
!   First check whether the block is present in the fdf file.
!   If it is not present, do nothing
    if (.not. fdf_block('WannierProjections',bfdf)) RETURN

!   Read the content of the block, line by line
    do while(fdf_bline(bfdf, pline))       

!     Read the index of the manifold to be considered
      if (.not. fdf_bmatch(pline,'i')) &   ! We expect that the first line
                                           !   of the block contains only one
                                           !   integer that is the sequential 
                                           !   number of the manifold. 
                                           !   That is the meaning of 'i'
 &      call die('Wrong format in the manifold index of the WannierProjections')
        
!     Assign the index of the manifold to the first integer found 
!     in the digested line
      index_manifold = fdf_bintegers(pline,1)

!     Go to the next line to 
!     read the initial and the final band of each manifold
      if (.not. fdf_bline(bfdf,pline)) &
 &      call die("No initial and final band specified for a manifold")
      if (.not. fdf_bmatch(pline,'ii')) &  ! We expect that the second line 
                                           !   of the block contains two
                                           !   integers that are, respectively
                                           !   the initial and the final bands
                                           !   of the manifold.
                                           !   That is the meaning of 'ii'
 &      call die('Wrong format in the initial/final band in WannierProjections')

!     Assign the indices of the initial and final band of the manifold 
!     as the first and second integer in the digested line, respectively
      manifold_bands_w90_in(index_manifold)%initial_band =fdf_bintegers(pline,1)
      manifold_bands_w90_in(index_manifold)%final_band   =fdf_bintegers(pline,2)

      if( manifold_bands_w90_in(index_manifold)%initial_band .gt.      &
 &        manifold_bands_w90_in(index_manifold)%final_band )           &
 &      call die('The final band should be larger than the initial one')

!     Compute the total number of bands in the manifold
      manifold_bands_w90_in(index_manifold)%number_of_bands=      &
 &        ( manifold_bands_w90_in(index_manifold)%final_band   -  &
 &          manifold_bands_w90_in(index_manifold)%initial_band )  + 1

!     Go to the next line to 
!     read the total number of bands that will enter the 
!     Wannier transformation
      if (.not. fdf_bline(bfdf,pline)) &
 &      call die("Number of bands to be orthogonalized not specified")
      if (.not. fdf_bmatch(pline,'i')) &   !  We expect that the third line
                                           !    of the block contains only one
                                           !    integer: 
                                           !    the number of "Wannier" 
                                           !    functions that will be produced.
                                           !    That is the meaning of 'i'
 &      call die('Wrong format in the number of Wannier functions')

!     Assign the number of Wannier functions required to the 
!     first integer in the digested line
      manifold_bands_w90_in(index_manifold)%numbands_w90_in = &
 &      fdf_bintegers(pline,1)  

      if( manifold_bands_w90_in(index_manifold)%numbands_w90_in .gt.   &
 &        manifold_bands_w90_in(index_manifold)%number_of_bands )      &
 &      call die('More Wannier functions requested than bands in the manifold')

!     Check whether the disentanglement procedure will be required:
      if( manifold_bands_w90_in(index_manifold)%number_of_bands /=      &
 &        manifold_bands_w90_in(index_manifold)%numbands_w90_in ) then
        manifold_bands_w90_in(index_manifold)%disentanglement = .true. 
      else
        manifold_bands_w90_in(index_manifold)%disentanglement = .false.
      endif

!     Go to the next line to 
!     read the indices of the atomic orbitals chosen as localized trial orbitals
      if (.not. fdf_bline(bfdf,pline)) &
 &      call die("No localized trial orbitals in WannierProjections")

!     We need as many localized trial orbitals as Wannier functions requested
      if( fdf_bnintegers(pline) .ne.                               &
 &        manifold_bands_w90_in(index_manifold)%numbands_w90_in )  &
 &      call die('We need as many localized trial orbitals as Wannier functions requested')

!     Allocate the array where the indices of the orbitals used as 
!     trial localized functions will be stored
      nullify( manifold_bands_w90_in(index_manifold)%orbital_indices )
      call re_alloc( manifold_bands_w90_in(index_manifold)%orbital_indices,  &
 &        1, manifold_bands_w90_in(index_manifold)%numbands_w90_in )

!     Read the indices of the orbitals used as 
!     trial localized functions will be stored
      do iorb = 1, manifold_bands_w90_in(index_manifold)%numbands_w90_in
        manifold_bands_w90_in(index_manifold)%orbital_indices(iorb) =        &
 &         fdf_bintegers(pline,iorb)
      enddo

!     Read the Number of iterations for the minimization of \Omega
      if (.not. fdf_bline(bfdf,pline)) &
 &         call die("No number of iterations for the minimization of \Omega")
      string_num    = fdf_bnames(pline,1)
      manifold_bands_w90_in(index_manifold)%num_iter = fdf_bintegers(pline,1)

!     Read whether the Wannier functions are plotted
      if (.not. fdf_bline(bfdf,pline)) &
 &         call die("No plotting of the Wannier functions")
      string_plot    = fdf_bnames(pline,1)
      manifold_bands_w90_in(index_manifold)%wannier_plot = .true.
      manifold_bands_w90_in(index_manifold)%wannier_plot_supercell(:) = &
 &       fdf_bintegers(pline,1)
!     If the UNK files are not written, then the Wannier functions will not
!     be plotted using the algorithms included in WANNIER90
!     They could be plotted using the coefficients of the wave functions
!     and DENCHAR
      if( .not. w90_in_siesta_compute_unk )    &
 &      manifold_bands_w90_in(index_manifold)%wannier_plot = .false.

!     Read whether the Fermi surface is plotted
      if (.not. fdf_bline(bfdf,pline)) &
 &         call die("No plotting of the Fermi surface")
      string_plot_f  = fdf_bnames(pline,1)
      manifold_bands_w90_in(index_manifold)%fermi_surface_plot = .true.

!     Read whether the Hamiltonian is written in a WF basis
      if (.not. fdf_bline(bfdf,pline)) &
 &         call die("No writting of the Hamiltonian in a WF basis")
      string_write_h  = fdf_bnames(pline,1)
      manifold_bands_w90_in(index_manifold)%write_hr = .true.

!     Read whether the lattice vectors, Hamiltonian, and position operator 
!     are written in a WF basis
      if (.not. fdf_bline(bfdf,pline)) &
 &         call die("No writting of the Hamiltonian and position operator in a WF basis")
      string_write_tb  = fdf_bnames(pline,1)
      manifold_bands_w90_in(index_manifold)%write_tb = .true.

!     If disentanglement is required, then read the outer and inner (frozen)
!     energy windows
      if( manifold_bands_w90_in(index_manifold)%disentanglement ) then
        if (.not. fdf_bline(bfdf,pline)) &
 &         call die("No outer energy window")
        manifold_bands_w90_in(index_manifold)%dis_win_min = fdf_bvalues(pline,1)
        manifold_bands_w90_in(index_manifold)%dis_win_max = fdf_bvalues(pline,2)

        if (.not. fdf_bline(bfdf,pline)) & 
 &         call die("No inner (frozen) energy window")
        manifold_bands_w90_in(index_manifold)%dis_froz_min =fdf_bvalues(pline,1)
        manifold_bands_w90_in(index_manifold)%dis_froz_max =fdf_bvalues(pline,2)
      else 
!       Set up the same default values as in WANNIER90 
!       (in subroutine parameters.F90)
        manifold_bands_w90_in(index_manifold)%dis_win_min =  -1.0_dp
        manifold_bands_w90_in(index_manifold)%dis_win_max =   0.0_dp
        manifold_bands_w90_in(index_manifold)%dis_froz_min = -1.0_dp
        manifold_bands_w90_in(index_manifold)%dis_froz_max =  0.0_dp
      endif

      write(manifold_bands_w90_in(index_manifold)%seedname_w90_in,  &
 &          "(a,'.manifold.',i1.1)")                                &
 &       trim(slabel),index_manifold

!    For debugging
!      write(6,'(/,a,i2)')                                                      &
! &      'read_w90_in_siesta_specs: Dumping the information for the manifold number ', &
! &       index_manifold
!
!!     Write the seed of the name of the file to dump the output
!      write(6,'(a,i1,1x,a)')                                                   &
! &      'read_w90_in_siesta_specs: Seed of the file to dump the output for manifold:',&
! &       index_manifold,                                                       &
! &       manifold_bands_w90_in(index_manifold)%seedname_w90_in
!
!!     Write the initial and the final band of the manifold
!      write(6,'(a,2i5)')                                                      &
! &      'read_w90_in_siesta_specs: Initial band, Final band              = ', &
! &      manifold_bands_w90_in(index_manifold)%initial_band,                   &
! &      manifold_bands_w90_in(index_manifold)%final_band                    
!
!!     Write the total number of bands in the manifold
!      write(6,'(a,i5)')                                                       &
! &      'read_w90_in_siesta_specs: Total number of bands in the manifold = ', &
! &      manifold_bands_w90_in(index_manifold)%number_of_bands
!
!!     Write the total number of bands to be orthogonalized
!      write(6,'(a,i5)')                                                        &
! &'read_w90_in_siesta_specs: Total number of bands in the Wannier transformation =', &
! &      manifold_bands_w90_in(index_manifold)%numbands_w90_in
!
!      if(manifold_bands_w90_in(index_manifold)%disentanglement) then
!         write(6,'(a)')                                                  &
! &        'read_w90_in_siesta_specs: Number of bands in the manifold is different of'
!         write(6,'(a)')                                                  &
! &        'read_w90_in_siesta_specs: the number of bands to be orthogonalized.'
!         write(6,'(a)')                                                  &
! &        'read_w90_in_siesta_specs: Disentanglement procedure required '
!      endif
!
!!     Write the indices of the orbitals used as initial localized guess
!      do iorb = 1, manifold_bands_w90_in(index_manifold)%numbands_w90_in
!        write(6,'(a,i5)')                                                      &
! & 'read_w90_in_siesta_specs: Index of the orbital used as initial localized guess =',&
! &      manifold_bands_w90_in(index_manifold)%orbital_indices(iorb)
!      enddo
!
!!     Write the values of the outer energy window for band disentanglement
!      write(6,'(a,l5)')                                             &
! &      'read_w90_in_siesta_specs: Disentanglement                    = ', &
! &      manifold_bands_w90_in(index_manifold)%disentanglement
!      write(6,'(a,i5,2f12.5)')                                      &
! &      'read_w90_in_siesta_specs: Outer energy window  (eV)          = ', &
! &      index_manifold,                                             &
! &      manifold_bands_w90_in(index_manifold)%dis_win_min,          &
! &      manifold_bands_w90_in(index_manifold)%dis_win_max 
!
!!     Write the values of the inner (frozen) energy window for 
!!     band disentanglement
!      write(6,'(a,i5,2f12.5)')                                      &
! &      'read_w90_in_siesta_specs: Inner (frozen) energy window  (eV) = ', &
! &      index_manifold,                                             &
! &      manifold_bands_w90_in(index_manifold)%dis_froz_min,         &
! &      manifold_bands_w90_in(index_manifold)%dis_froz_max 
!
!!     Write the values of the inner (frozen) energy window for 
!!     band disentanglement
!      write(6,'(a,i5)')                                                    &
! &      'read_w90_in_siesta_specs: Number of iterations for the minimization = ', &
! &      manifold_bands_w90_in(index_manifold)%num_iter
!
!!     Write the values of the inner (frozen) energy window for 
!!     band disentanglement
!      write(6,'(a,l5)')                                                    &
! &      'read_w90_in_siesta_specs: Plot the Wannier functions?               = ', &
! &      manifold_bands_w90_in(index_manifold)%wannier_plot
!      write(6,'(a,i5)')                                                    &
! &      'read_w90_in_siesta_specs: Plot the Wannier functions?               = ', &
! &      manifold_bands_w90_in(index_manifold)%wannier_plot_supercell(1)
!      write(6,'(a,l5)')                                                    &
! &      'read_w90_in_siesta_specs: Plot the Fermi surface?              ', &
! &      manifold_bands_w90_in(index_manifold)%fermi_surface_plot
!      write(6,'(a,l5)')                                                    &
! &      'read_w90_in_siesta_specs: Write the Hamiltonian in the basis of WF? = ', &
! &      manifold_bands_w90_in(index_manifold)%write_hr
!      write(6,'(a,l5)')                                                    &
! &      'read_w90_in_siesta_specs: Write the tight-binding elements in the basis of WF? = ', &
! &      manifold_bands_w90_in(index_manifold)%write_tb
!!     End debugging

    enddo ! end loop over all the lines in the block WannierProjections

!   Allocate the pointer where the chemical potential associated with 
!   a given Wannier function will be stored.
!   The correction will be applied only on the bands associated with 
!   the first manifold, that is why we are using this number of bands
!   to allocate the variable.
    nullify( chempotwann_val )
    call re_alloc( chempotwann_val,                              &
 &                 1, manifold_bands_w90_in(1)%numbands_w90_in,  &
 &                 1, spin%H,                                    & 
 &                 name='chempotwann_val',                       &
 &                 routine='read_w90_in_siesta_specs')
    chempotwann_val = 0.0_dp
    
!   Read the chemical potential associated with a given Wannier function
!   First check whether the block is present in the fdf file.
!   If it is not present, do nothing
    if (.not. fdf_block('ChemicalPotentialWannier',bfdf)) RETURN

!   If the block with the chemical potentials is present,
!   set to true the flag to compute the shifts of the matrix elements
    compute_chempotwann = .true. 

!   Read the content of the block, line by line
    do while(fdf_bline(bfdf, pline))       
      if( spin%H .eq. 1) then
        if (.not. fdf_bmatch(pline,'ivn')) &   ! We expect that the first line
 &        call die('Wrong format in ChemicalPotentialWannier')
      elseif( spin%H .eq. 2) then
        if (.not. fdf_bmatch(pline,'ivvn')) &  ! We expect that the first line
 &        call die('Wrong format in ChemicalPotentialWannier')
      endif
      factor = fdf_convfac( fdf_bnames(pline,1), 'Ry' )
      index_wannier = fdf_bintegers(pline,1)
      if( spin%H .eq. 1) then
        chempotwann_val(index_wannier,1) = fdf_breals(pline,1) * factor
      else if( spin%H .eq. 2) then
        chempotwann_val(index_wannier,1) = fdf_breals(pline,1) * factor
        chempotwann_val(index_wannier,2) = fdf_breals(pline,2) * factor
      endif
    enddo ! end loop over all the lines in the block ChemicalPotentialWannier

!!   For debugging
!    do index_wannier = 1, manifold_bands_w90_in(1)%numbands_w90_in
!      write(6,'(a,i5,a,2f12.5)')                                              &
! &      'read_w90_in_siesta_specs: Chemical potential for Wannier number ',  &
! &      index_wannier, ' = ' , chempotwann_val(index_wannier,:) 
!    enddo 
!!   End debugging
        
  end subroutine read_w90_in_siesta_specs

!> \brief General purpose of the subroutine read_kpoints_wannier:
!! process the information in the fdf file
!! required to generate the k-point sampling for Wannierization.
!! 
!! In this subroutine, we process the information in the fdf file 
!! required to generate the k-point sampling that will be used inside 
!! WANNIER90 for the Wannier transformation.
!! The block that is readed and digested is kMeshforWannier.
!! Example: 
!!
!! %block kMeshforWannier
!!
!!   20  20  1
!!
!! %endblock kMeshforWannier
!!
!! where the three integers are the number of k-points along the corresponding 
!! reciprocal lattice vector
!! The algorithm to generate the k-points in reciprocal units is borrowed from 
!! the utility kmesh.pl of WANNIER90.
!!
!! Then, and following the recipes given in the Appendix of \cite Marzari-97,
!! we compute the distance to nearest-neighbour shells of k-points
!! the vectors connecting neighbour k-points in the Monkhorst-Pack mesh,
!! and check whether the completeness relation is fully satisfied 
!! [Eq. (B1) of Ref. \cite Marzari-97]  
!! This is done within the subroutine kmesh_get, borrowed from 
!! WANNIER90 (version 3.0.0) \cite Wannier90.
!!
!! Finally, we populate the variables related with the neighbour k-points in 
!! the Monkhorst-Pack mesh in the module where all the Wannier90 parameters are
!! stored.

  subroutine read_kpoints_wannier

    use fdf
    use siesta_geom,    only: ucell            ! Unit cell lattice vectors 
                                               !   in real space as used inside
                                               !   SIESTA
                                               !   First  index: component
                                               !   Second index: vector
                                               !   In Bohrs
    use w90_in_siesta_types,   only: reclatvec_w90_in 
                                               ! Reciprocal lattice vectors 
                                               !   computed from real_lattice
                                               !   The factor 2.0 * pi 
                                               !   is included
                                               !   First  index: component
                                               !   Second index: vector
                                               !   In Angstroms^-1
    use w90_in_siesta_types,   only: kmesh_w90_in  
                                               ! Number of divisions along the
                                               !   reciprocal lattice vectors
    use w90_in_siesta_types,   only: numkpoints_w90_in
                                               ! Total number of k-points
    use w90_in_siesta_types,   only: kpointsfrac_w90_in   
                                               ! Coordinates of the k-points
                                               !   In fractional units, i. e.
                                               !   in units of the reciprocal
                                               !   lattice vectors
                                               !   First  index: component
                                               !   Second index: vector
    use w90_io,         only: stdout           ! Unit on which stdout is written
    use w90_parameters, only: real_lattice     ! Unit cell lattice vectors
                                               !   in real space as used inside
                                               !   WANNIER90
                                               !   BE CAREFUL: In WANNIER90,the 
                                               !   order of the indices is 
                                               !   inverted with respect to 
                                               !   SIESTA
                                               !   First  index: vector   
                                               !   Second index: component
                                               !   The same as the transpose 
                                               !   of ucell,
                                               !   but with different units
                                               !   In Angstroms
    use w90_parameters, only: recip_lattice    ! Reciprocal lattice vectors 
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
                                               !   The same as the transpose of
                                               !   reclatvec_w90_in
                                               !   In Angstroms^-1
    use w90_parameters, only: mp_grid          ! Number of divisions along the
                                               !   reciprocal lattice vectors.
                                               !   It is the same as 
                                               !   kmesh_w90_in
                                               !   But this is the variable
                                               !   that will be transferred to 
                                               !   the Wannier90 subroutines
    use w90_parameters, only: num_kpts         ! Total number of k-points that
                                               !   will be used within 
                                               !   Wannier90
                                               !   It is the same as 
                                               !   numkpoints_w90_in
                                               !   But this is the variable
                                               !   that will be transferred to 
                                               !   the Wannier90 subroutines
    use w90_parameters, only: gamma_only       ! Only the gamma point will be
                                               !   used within Wannier90?
    use w90_parameters, only: kpt_latt         ! Coordinates of the k-points
                                               !   that will be used within
                                               !   Wannier90
                                               !   In fractional units, i. e.
                                               !   in units of the reciprocal
                                               !   lattice vectors
                                               !   First  index: component
                                               !   Second index: vector
                                               !   It is the same as 
                                               !   kpointsfrac_w90_in
                                               !   but this is the variable that
                                               !   will be transferred to 
                                               !   Wannier90
    use w90_parameters, only : nntot           ! Number of nearest neighbours 
                                               !   belonging to each k-point of
                                               !   the Monkhorst-Pack mesh 
    use w90_parameters, only : nnlist          ! The list of nearest neighbours
                                               !   for each k-point
    use w90_parameters, only : nncell          ! The vector, in 
                                               !   fractional reciprocal lattice
                                               !   coordinates, that brings the
                                               !   nnth nearest neighbour
                                               !   of k-point nkp to its 
                                               !   periodic image that is needed
                                               !   for computing the overlap 
                                               !   M(k,b).
    use w90_parameters, only: param_read       ! Subroutine to read the 
                                               !    parameters required by
                                               !    WANNIER90 and populate 
                                               !    derived values
    use w90_kmesh,      only: kmesh_get        ! Main routine to calculate the 
                                               !    nearest neighbour vectors
                                               !    in reciprocal space 
                                               !    (the b-vectors),
                                               ! It also checks the completeness
                                               !    relation, Eq. (B1) in 
                                               !    N. Marzari et al.
                                               !    PRB 56, 12847 (1997)
    use w90_in_siesta_types,   only: nncount_w90_in   
                                               ! Same as nntot
                                               !    but in the module where the 
                                               !    variables that controls the
                                               !    Wannier90 in SIESTA 
                                               !    are stored
    use w90_in_siesta_types,   only: nnlist_w90_in    
                                               ! Same as nnlist
                                               !    but in the module where the 
                                               !    variables that controls the
                                               !    Wannier90 in SIESTA 
                                               !    are stored
    use w90_in_siesta_types,   only: nnfolding_w90_in 
                                               ! Same as nncell
                                               !    but in the module where the 
                                               !    variables that controls the
                                               !    Wannier90 in SIESTA 
                                               !    are stored
    use w90_in_siesta_types,   only: bvectorsfrac_w90_in
                                               ! Vectors that connect each mesh
                                               !    k-point
                                               !    to its nearest neighbours
    use m_digest_nnkp,  only: chosing_b_vectors! Subroutine that computes the b
                                               !    vectors that connect each
                                               !    mesh k-point to its 
                                               !    nearest neighbours.
    use w90_comms,      only: on_root          ! Set up the communications

!   Internal variables
    integer :: ik                   
    integer :: ikx, iky, ikz       
    integer :: i
    integer :: nkp
    integer :: nn


    type(block_fdf)            :: bfdf
    type(parsed_line), pointer :: pline
    
!   Setup the output unit for WANNIER90
    stdout = 6
    if( Node .eq. 0 ) then
!     Set on_root = .true. only if you want to print the information of the
!     k-point sampling
!      on_root = .true. 
      on_root = .false.
    else 
      on_root = .false.
    endif

!   Read the data to generate the grid in reciprocal space that will be used
!   for the Wannier Projections
    if (.not. fdf_block('kMeshforWannier',bfdf)) RETURN

    do while(fdf_bline(bfdf, pline))     
      if (.not. fdf_bmatch(pline,'iii'))        &   ! We expect that each line
                                                    !   contains three integers
                                                    !   That is the meaning of
                                                    !   iii
                                                    ! The first integer is the
                                                    !   number of divisions 
                                                    !   the first reciprocal 
                                                    !   lattice vector and so on
 &      call die('Wrong format in kMeshforWannier')
      kmesh_w90_in(1) = fdf_bintegers(pline,1)
      kmesh_w90_in(2) = fdf_bintegers(pline,2)
      kmesh_w90_in(3) = fdf_bintegers(pline,3)
    enddo 

!   Define the total number of k-points used in the Wannier projection
    numkpoints_w90_in = kmesh_w90_in(1) * kmesh_w90_in(2) * kmesh_w90_in(3)

!   Compute and store the components of the k-points in fractional units
    nullify( kpointsfrac_w90_in )
    call re_alloc( kpointsfrac_w90_in, 1, 3, 1, numkpoints_w90_in,  &
 &                 name='kpointsfrac_w90_in', routine='read_kpoints_wannier')

!   The algorithm to generate the k-points in fractional units 
!   is borrowed from the utility kmesh.pl of Wannier90
    ik = 0
    do ikx = 0, kmesh_w90_in(1) - 1
      do iky = 0, kmesh_w90_in(2) - 1
        do ikz = 0, kmesh_w90_in(3) - 1
          ik = ik + 1
          kpointsfrac_w90_in(1,ik) = (ikx*1.0_dp)/kmesh_w90_in(1)
          kpointsfrac_w90_in(2,ik) = (iky*1.0_dp)/kmesh_w90_in(2)
          kpointsfrac_w90_in(3,ik) = (ikz*1.0_dp)/kmesh_w90_in(3)
        enddo 
      enddo 
    enddo 

!!   For debugging
!    write(6,'(a,a,3i5)')'read_kpoints_wannier: Number of subdivisions ',  &  
! &    'of the reciprocal vectors for Wannier = ', kmesh_w90_in(:)
!    write(6,'(a,a,3i5)')'read_kpoints_wannier: Number of k-points used ', & 
! &    'in the Wannier90 projection = ', numkpoints_w90_in
!    do ik = 1, numkpoints_w90_in
!      write(6,'(a,a,i5,3f12.5)')'read_kpoints_wannier: k-points in ',     & 
! &      'fractional units:', ik, kpointsfrac_w90_in(:,ik)
!    enddo
!!   End debugging

!   Transform the units of the unit cell lattice vector from Bohrs to Ang,
!   as required by Wannier90
!   In WANNIER90 the real space lattice vector matrix is the 
!   transpose one used in SIESTA
!   In WANNIER90 the first component is the vector and the second the componen
!   In SIESTA the first component is the component and the second the vector
!   We make the transposition here to transfer the info to WANNIER90
    real_lattice  = transpose(ucell) * 0.529177_dp
!!   For debugging
!    do ik = 1, 3
!      write(6,'(a,3f12.5)')'read_kpoints_wannier: real_lattice = ',   &
! &      real_lattice(ik,:)
!    enddo
!!   End debugging

!   Compute the reciprocal lattice vectors as required by WANNIER90
!   The factor 2.0 * pi is included
    call reclat( real_lattice, reclatvec_w90_in, 1 )
!   Save the reciprocal lattice vectors in the variable that will be
!   transferred to WANNIER90
!   Since the reciprocal lattice vectors are already computed with the
!   transpose matrix of real space lattice vectors, there is no need
!   to transpose reclatvec_w90_in again
    recip_lattice = reclatvec_w90_in
!!   For debugging
!    do ik = 1, 3
!      write(6,'(a,3f12.5)')'read_kpoints_wannier: recip_lattice = ',   &
! &      recip_lattice(ik,:)
!    enddo
!!   End debugging

!   Save the number of subdivisions along the three reciprocal lattice vectors
!   that will be transferred to Wannier90
    mp_grid       = kmesh_w90_in
!!   For debugging
!    write(6,'(a,3i5)')'read_kpoints_wannier: mp_grid =',  mp_grid
!!   End debugging

!   Save the total number of k-points that will be transferred to Wannier90
    num_kpts      = numkpoints_w90_in
!!   For debugging
!    write(6,'(a,2i5)')'read_kpoints_wannier: Node, num_kpts =', Node, num_kpts 
!!   End debugging

!   Set if only the gamma points will be used
    if (num_kpts .ne. 1) then
      gamma_only    = .false.
    else
      gamma_only    = .true.
    endif 
!!   For debugging
!    write(6,'(a,l5)')'read_kpoints_wannier: gamma_only =', gamma_only
!!   End debugging

!   Allocate the variable where the k-points will be transferred to WANNIER90
!   in fractional units (i.e. in units of the reciprocal lattice vectors
    allocate ( kpt_latt(3,num_kpts) )
    kpt_latt = kpointsfrac_w90_in
!!   For debugging
!    do ik = 1, num_kpts
!      write(6,'(a,i5,3f12.5)')'read_kpoints_wannier: kpt_latt = ',   &
! &      ik, kpt_latt(:,ik)
!   enddo
!!   End debugging

!   Now, we generate the vectors that connect nearest-neighbour shells
!   in reciprocal space, together with the completeness relation,
!   Eq. (B1) in N. Marzari et al. Physical Review B 56, 12847 (1997)
!   For that, we use the subroutine kmesh_get, directly borrowed
!   from WANNIER90 (at this moment, from version 3.0.0)
!   But, before using this subroutine, we have to populate some parameters
!   required in WANNIER90, also using the subroutine param_read
!   borrowed from the version 3.0.0 of WANNIER90 with small modifications
!    
    call param_read
    call kmesh_get

!   Store the number of nearest neighbours belonging to each k-point of the 
!   Monkhorst-Pack mesh
    nncount_w90_in = nntot

!   Initialize the list of neighbour k-points
    nullify( nnlist_w90_in    )
    nullify( nnfolding_w90_in )

    call re_alloc( nnlist_w90_in, 1, numkpoints_w90_in, 1, nncount_w90_in,   &
 &                 name = "nnlist_w90_in", routine = "read_kpoints_wannier" )
    call re_alloc( nnfolding_w90_in, 1, 3, 1, numkpoints_w90_in,             &
 &                 1, nncount_w90_in, name = "nnfolding_w90_in",             &
 &                 routine = "read_kpoints_wannier" )

!   Store the list of nearest neighoburs for each k-point
    nnlist_w90_in     = nnlist

!   Store the vector, in fractional reciprocal lattice coordinates,
!   that brings the nnth nearest neighbour of k-point nkp to its periodic image
!   that is needed for computing the overlap matrices M_mn(k,b)
    nnfolding_w90_in  = nncell

!   Compute the vectors that connect each mesh k-point 
!   to its nearest neighbours
    call chosing_b_vectors( kpointsfrac_w90_in, nncount_w90_in,  &
 &                          nnlist_w90_in, nnfolding_w90_in,     &
 &                          bvectorsfrac_w90_in )

!!     For debugging
!      write(6,'(a)') 'begin nnkpts'
!      write(6,'(i4)') nntot
!      do nkp=1,num_kpts
!         do nn=1,nntot
!            write(6,'(2i6,3x,3i4)') &
!               nkp,nnlist(nkp,nn),(nncell(i,nkp,nn),i=1,3)
!         end do
!      end do
!      write(6,'(a/)') 'end nnkpts'
!!     End debugging

  end subroutine read_kpoints_wannier

!> \brief General purpose of the subroutine set_excluded_bands_w90_in:
!! identify the bands that will be considered for wannierization
!! and the atomic orbitals that will be used as initial guess for the
!! projections.
!!
!! Within this subroutine we:
!! 1. Identify the bands that will be considered for the wannierization of
!! a given manifold. These can be found in the variable:
!! manifold_bands_w90_in(index_manifold)%isexcluded.
!! 2. Identify the numerical atomic orbitals of the basis set that will
!! be considered as local functions for the initial guess.
!! These can be found in the variable
!! manifold_bands_w90_in(index_manifold)%orbexcluded.
!! 3. Set a sequential index for the orbitals that will be used
!! as initial guess, stored in the variable
!! manifold_bands_w90_in(index_manifold)%orb_in_manifold.
!! This routine is called only by the IOnode.

  subroutine set_excluded_bands_w90_in( index_manifold )
    
    integer, intent(in) :: index_manifold

    integer :: iband
    integer :: iorb
    integer :: iterex
    integer :: index
    integer :: index_orb_included

    nullify( manifold_bands_w90_in(index_manifold)%isexcluded )
    call re_alloc( manifold_bands_w90_in(index_manifold)%isexcluded,      &
 &                 1, no_u,                                               &
 &                 name='isexcluded',                                     &
 &                 routine='set_excluded_bands_w90_in' )

    nullify( manifold_bands_w90_in(index_manifold)%orbexcluded )
    call re_alloc( manifold_bands_w90_in(index_manifold)%orbexcluded,     &
 &                 1, no_u,                                               &
 &                 name='orbexcluded',                                    &
 &                 routine='set_excluded_bands_w90_in' )

    nullify( manifold_bands_w90_in(index_manifold)%orb_in_manifold )
    call re_alloc( manifold_bands_w90_in(index_manifold)%orb_in_manifold, &
 &                 1, no_u,                                               &
 &                 name='orb_in_manifold',                                &
 &                 routine='set_excluded_bands_w90_in' )


!   By default, all the bands are excluded from the calculation
    manifold_bands_w90_in(index_manifold)%isexcluded(:) = .true.

    do iband = 1, no_u
!     Exclude the corresponding bands from the computation
      if( (iband .ge. manifold_bands_w90_in(index_manifold)%initial_band) &
 &        .and.                                                           &
 &        (iband .le. manifold_bands_w90_in(index_manifold)%final_band) ) then
        manifold_bands_w90_in(index_manifold)%isexcluded( iband ) = .false.
      endif
    enddo

!!   For debugging
!    write(6,*) 
!    do iterex = 1, no_u
!     write(6,'(a,2i5,l5)')                                                     &
! &  'set_excluded_bands_w90_in: Node, excluded_bands_w90_in (iband, excluded)',&
! &     Node, iterex, manifold_bands_w90_in(index_manifold)%isexcluded(iterex)
!    enddo
!!   End debugging

!   By default, all the orbitals are excluded from the calculation
    manifold_bands_w90_in(index_manifold)%orbexcluded(:) = .true.
    do iorb = 1, manifold_bands_w90_in(index_manifold)%numbands_w90_in
      index_orb_included =                                          & 
 &      manifold_bands_w90_in(index_manifold)%orbital_indices(iorb)
      manifold_bands_w90_in(index_manifold)%orbexcluded(index_orb_included) = &
 &      .false. 
    enddo

    index = 0
    manifold_bands_w90_in(index_manifold)%orb_in_manifold(:) = 0
    do iorb = 1, manifold_bands_w90_in(index_manifold)%numbands_w90_in
      index = index + 1
      index_orb_included = manifold_bands_w90_in(index_manifold)%orbital_indices(iorb)
      manifold_bands_w90_in(index_manifold)%orb_in_manifold(index_orb_included) = index
    enddo

!!   For debugging
!    write(6,*) 
!    do iterex = 1, no_u
!     write(6,'(a,2i5,l5,i5)')                                            &
! &     'set_excluded_bands_w90_in: orbitals excluded: manifold, orb_unit cell, excluded, orb_mani', &
! &     index_manifold, iterex,                                           &
! &     manifold_bands_w90_in(index_manifold)%orbexcluded(iterex),        &
! &     manifold_bands_w90_in(index_manifold)%orb_in_manifold(iterex)    
!    enddo
!!   End debugging

  end subroutine set_excluded_bands_w90_in

!
! ------------------------------------------------------------------------------
!
!> \brief General purpose of the subroutine compute_matrices:
!! Compute the overlap between periodic parts of wave functions at neighbour
!! k-points, and the projections between the wave functions and the initial
!! guesses of the localized functions.
!!
!! In this subroutine we compute and dump in the corresponding files:
!! 1. The matrix elements of the plane waves connecting neighbour
!! points in the first-Brillouin zone for Wannierization in a basis
!! of atomic orbitals
!! (done in the subroutine compute_pw_matrix)
!! 2. The overlap matrices between the periodic part of wave functions
!! at neighbour k-points,
!! \f{eqnarray*}{
!! M_{m n}^{(\vec{k}, \vec{b})} =
!!        \langle u_{m \vec{k}} \vert u_{n \vec{k} + \vec{b}} \rangle, 
!! \f}
!! that is the Eq. (27) of the Ref. \cite Marzari-12, or
!! Eq. (25) of Ref. \cite Marzari-97
!! (done in the subroutine mmn)
!! 3. The overlap matrices between the eigenstates of the one-particle
!! Hamiltonian and the localized trial orbitals, taken as initial guess,
!! \f{eqnarray*}{
!! A_{mn} (\vec{k}) = \langle \psi_{m\vec{k}} \vert g_{n} \rangle,
!! \f}
!! as presented right before Eq. (17) of Ref. \cite Marzari-12
!! (done in the subroutine amn)
!! 4. The eigenvalues for the bands to be wannierized
!! (done in the subroutine writeeig)
!! 5. The values of the periodic part of the wavefunctions at the points
!! of a real space grid. This must be done if we want to plot the 
!! Wanniers
!! (done in the subroutine writeunk)

  subroutine compute_matrices( ispin, index_manifold )

    use m_switch_local_projection, only: seedname
    use m_switch_local_projection, only: nncount
    use m_switch_local_projection, only: numkpoints
    use m_switch_local_projection, only: bvectorsfrac
    use m_switch_local_projection, only: Amnmat
    use files,                     only: slabel        ! Short system label,
                                                       !   used to generate file
                                                       !   names
    use m_spin,                    only: spin          ! Spin configuration 
                                                       !   for SIESTA
    use sys,                       only: die

    integer, intent(in) :: ispin           ! Spin component
    integer, intent(in) :: index_manifold  ! Index of the manifold

!   Internal variables
    integer :: ik                          ! Counter for the loops on k-points
    integer :: n                           ! Counter for the loops on bands
    integer :: m                           ! Counter for the loops on bands

    integer :: number_of_bands_in_manifold_local
    integer :: number_of_bands_to_project

!   Variables related with the reading of the amn files
    integer            :: iman          ! Counter for the manifolds
    integer            :: amn_in        ! Logical unit for the file amn file
    character(len=100) :: filename      ! Name of the files where the amn
                                        !   matrices of the manifolds to be
                                        !   mixed are stored
    character(len=50)  :: dummy         ! Dummy variable
    integer            :: nb_tmp        ! Number of bands (read from the file)
    integer            :: nkp_tmp       ! Number of k-points (read from the file
    integer            :: nw_tmp        ! Number of Wanniers (read from the file
!    integer            :: m             ! Band index 
!                                        !    (read from the loop in the file)  
!    integer            :: n             ! Projector index 
!                                        !    (read from the loop in the file)  
    integer            :: nkp           ! k-point index 
                                        !    (read from the loop in the file)  
    integer            :: num_amn       ! Total number of entries in the Amn
                                        !   matrix
    integer            :: ncount        ! Counter for the Amn matrix elements
    integer            :: iproj         ! Counter for the loop on projectors
    integer            :: mband         ! Counter for the loop on bands
    real(dp)           :: a_real        ! Real part of the Amn matrix element
    real(dp)           :: a_imag        ! Imaginary part of the 
                                        !    Amn matrix element

    external     :: io_assign              ! Assign a logical unit
    external     :: io_close               ! Close a logical unit

#ifdef MPI
    integer            :: MPIError
#endif


    if( spin%H .eq. 1) then
      write(seedname,"(a,'.manifold.',i1.1)")                    &
 &         trim(slabel),index_manifold
    else if( spin%H .gt. 1) then
      write(seedname,                                            &
 &          "(a,'.manifold.',i1.1,'.spin.',i1.1)")               &
 &       trim(slabel),index_manifold,ispin
    endif

    number_of_bands_in_manifold_local =                                  &
 &        manifold_bands_w90_in(index_manifold)%nincbands_loc_w90_in
    number_of_bands_to_project =                                         &
 &        manifold_bands_w90_in(index_manifold)%numbands_w90_in

!   Compute the matrix elements of the plane wave,
!   for all the wave vectors that connect a given k-point to its nearest
!   neighbours
!    if( first_chempotwann ) call compute_pw_matrix( nncount, bvectorsfrac )
    call compute_pw_matrix( nncount, bvectorsfrac )

!   Compute the overlaps between Bloch orbitals at neighboring k points
!    if( first_chempotwann ) call mmn( ispin )
    call mmn( ispin )

!   Compute the overlaps between Bloch states onto trial localized orbitals
!    if( first_chempotwann ) call amn( ispin, index_manifold )
    call amn( ispin, index_manifold )

!   Check if the computations of the position operator matrix elements
!   between manifolds has to be carried out
    if( r_between_manifolds ) then
      if( index_manifold .eq. 3) then
#ifdef MPI
        call MPI_barrier(MPI_Comm_world,MPIError)
#endif
        if( IONode ) then
          do nkp = 1, numkpoints
            do iproj = 1, manifold_bands_w90_in(1)%number_of_bands 
              do mband = manifold_bands_w90_in(1)%number_of_bands+1,  &
                         manifold_bands_w90_in(1)%number_of_bands +   & 
                         manifold_bands_w90_in(2)%number_of_bands 
                Amnmat(mband, iproj, nkp) = cmplx(0.0_dp,0.0_dp,kind=dp)
              enddo
            enddo
            do iproj = manifold_bands_w90_in(1)%number_of_bands+1,    &
                       manifold_bands_w90_in(3)%number_of_bands
              do mband = 1, manifold_bands_w90_in(1)%number_of_bands 
                Amnmat(mband, iproj, nkp) = cmplx(0.0_dp,0.0_dp,kind=dp)
              enddo
            enddo
          enddo

!         Write the Amn overlap matrices in a file, in the format required
!         by Wannier90
          call writeamn( ispin )

        endif
#ifdef MPI
        call MPI_barrier(MPI_Comm_world,MPIError)
#endif
      endif
    endif

!   Write the selected eigenvalues within the manifold in the file
!   SystemLabel.eigW
    if( IOnode ) call writeeig( ispin )

    if( w90_in_siesta_compute_unk ) call writeunk( ispin )

    if (IONode) then
      write(6,'(/,a)')  &
 &     'compute_matrices: All the information dumped in the corresponding files'
      write(6,'(a)')  &
 &     'compute_matrices: End of the interface between Siesta and Wannier90'
    endif

!! For debugging
!#ifdef MPI
!    call MPI_barrier(MPI_Comm_world,MPIError)
!#endif
!    call die()
!! End debugging

     return 

102  call die('Error: Problem opening input file '//trim(filename))



  end subroutine compute_matrices


!> \brief General purpose of the subroutine compute_wannier:
!! populate all the variables required by WANNIER90 and call the
!! corresponding routines for Wannierization.
!!
!! Within this subroutine:
!! 1. We populate all the variables within the WANNIER90 modules
!!    required to run WANNIER90.
!!    Those variables are transferred from different modules in SIESTA,
!!    mostly m_switch_local_projection
!! 2. We call the different routines of the WANNIER90 code to perform
!!    the Wannierization.
!!    This part is a copy, almost verbatim, of the
!!    subroutine wannier_prog.F90 in WANNIER90, version 3.0.0

  subroutine compute_wannier( ispin, index_manifold )

!
!   General WANNIER90 variables 
!
    use w90_parameters, only : num_bands       ! number of bands
    use w90_parameters, only : num_wann        ! number of wannier 
                                               !   functions
    use w90_parameters, only : num_proj        ! number of projections
    use w90_parameters, only : num_iter        ! number of iterations for
                                               !   the minimization of
                                               !   \Omega
!
!   Variables related with the atomic structure coming from SIESTA
!
    use siesta_geom,    only : na_u            ! number of atoms in the
                                               !   unit cell
    use siesta_geom,    only : xa              ! atomic positions 
                                               !   in cartesian coordinates
                                               !   units in Bohrs
    use siesta_geom,    only : xa_last         ! atomic positions 
                                               !   in cartesian coordinates
                                               !   units in Bohrs
    use siesta_geom,    only : isa             ! species index of each atom
    use units,          only : Ang             ! conversion factor from 
                                               !   Ang to Bohrs
    use atm_types,      only : nspecies        ! number of different 
                                               !   chemical species
    use atm_types,      only : species         ! information about the different
                                               !   chemical species
!
!   Variables related with the atomic structure coming from WANNIER90
!
    use w90_parameters, only : lenconfac       ! conversion factor for
                                               !   unit cell
    use w90_parameters, only : num_atoms       ! number of atoms in the 
                                               !   unit cell
    use w90_parameters, only : num_species     ! number of atomic species
    use w90_parameters, only : atoms_symbol    ! atomic symbols
    use w90_parameters, only : atoms_label     ! atomic labels
    use w90_parameters, only : atoms_species_num    
                                               ! number of atoms of each
                                               !   species
    use w90_parameters, only : atoms_pos_cart  ! atomic positions
                                               !   in cartesian coordinates
                                               !   units in Angstroms
    use w90_parameters, only : lenconfac       ! conversion factor for
                                               !   length units
                                               !   lenconfac = 1.0 means
                                               !      that the lengths
                                               !      are in Angstroms

!
!   Variables related with the k-points sampling coming from SIESTA
!
    use w90_in_siesta_types,   only: kmesh_w90_in  
                                               ! Number of divisions along the
                                               !   reciprocal lattice vectors
    use m_switch_local_projection, only: numkpoints
                                               ! Number of k-points in the
                                               !    Monkhorst-Pack grid that
                                               !    will be used in the 
                                               !    Wannierization
    use m_switch_local_projection, only: nncount 
                                               ! The number of nearest
                                               !   neighbours belonging to
                                               !   each k-point of the 
                                               !   Monkhorst-Pack mesh
    use m_switch_local_projection, only: nnlist_neig
                                               ! Index of the
                                               !   inn-neighbour of ikp-point
                                               !   in the Monkhorst-Pack grid
    use m_switch_local_projection, only: nnfolding
                                               ! nnfolding(i,ikp,inn) is the
                                               !   i-component of the reciprocal
                                               !   lattice vector,
                                               !   in reduced units, that brings
                                               !   the inn-neighbour specified 
                                               !   in nnlist_neig (which is in
                                               !   the first BZ)
                                               !   to the actual 
                                               !   \vec{k} + \vec{b}
                                               !   that we need.
                                               !   In reciprocal lattice units.
!   Variables related with the k-points sampling coming from WANNIER90
!
    use w90_parameters, only: mp_grid          ! Number of divisions along the
                                               !   reciprocal lattice vectors.
                                               !   It is the same as 
                                               !   kmesh_w90_in
                                               !   But this is the variable
                                               !   that will be transferred to 
                                               !   the Wannier90 subroutines
    use w90_parameters, only : gamma_only      ! Only the Gamma point?
    use w90_parameters, only : num_kpts        ! number of k-points
    use w90_parameters, only : nntot           ! total number of neighbours
                                               !   for each k-point
    use w90_parameters, only : nnlist          ! list of neighbours for 
                                               !   each k-point
    use w90_parameters, only : nncell          ! The vector, in 
                                               !   fractional reciprocal 
                                               !   lattice coordinates,
                                               !   that brings the
                                               !   nnth nearest neighbour
                                               !   of k-point nkp to its
                                               !   periodic image that is 
                                               !   needed for computing 
                                               !   the overlap M(k,b).
    use w90_parameters, only : kpt_latt        ! kpoints in lattice vectors
    use w90_parameters, only : bk              ! the b-vectors that go 
                                               !   from each k-point to 
                                               !   its neighbours
    use w90_parameters, only : wb              ! weights associated with 
                                               !   neighbours of each 
                                               !   k-point
!
!   Variables related with the disentanglement procedure
!
    use w90_parameters, only : have_disentangled
    use w90_parameters, only : disentanglement ! logical value
                                               !   .true.:
                                               !   disentanglement active
                                               !   .false.:
                                               !   disentanglement inactive
    use w90_parameters, only : dis_win_min     ! lower bound of the 
                                               !   disentanglement outer 
                                               !   window
    use w90_parameters, only : dis_win_max     ! upper bound of the 
                                               !   disentanglement outer 
                                               !   window
    use w90_parameters, only : dis_froz_min    ! lower bound of the 
                                               !   disentanglement inner
                                               !   (frozen) window
    use w90_parameters, only : dis_froz_max    ! upper bound of the 
                                               !   disentanglement inner
                                               !   (frozen) window
!
!   Variables related with the input/output coming from SIESTA   
!
    use files,          only: slabel           ! Short system label,
                                               !   used to generate file names
    use m_spin,         only: spin             ! Spin configuration for SIESTA
!
!   Variables related with the input/output coming from WANNIER90
!
    use w90_io,         only: stdout           ! unit on which stdout is
                                               !    written
    use w90_io,         only : seedname        ! Seed for the name of the  
                                               !    file where the matrix 
                                               !    elements of the
                                               !    position and hamiltonian
                                               !    operator will be written
    use w90_io,         only : maxlen          ! max column width of 
                                               !    input file
    use w90_parameters, only : timing_level    ! Verbosity of timing output 
                                               !   info
!
!   Variables related with the post-processing coming from SIESTA
!
    use m_switch_local_projection, only: eo    ! Eigenvalues of the Hamiltonian
                                               !    at the numkpoints
    use m_energies,     only: ef               ! Fermi energy
    use units,          only: eV               ! Conversion factor from Ry to eV
    use siesta_options, only: n_wannier_manifolds
                                               ! Number of manifolds to be
                                               !   wannierized
!
!   Variables related with the post-processing coming from WANNIER90
!
    use w90_parameters, only : wannier_plot    ! are we going to plot the
                                               !   Wannier functions?
    use w90_parameters, only : wannier_plot_supercell 
                                               ! size of the supercell to 
                                               !   plot the WF
    use w90_parameters, only : fermi_surface_plot  
                                               ! are we going to plot the
                                               !   Fermi surface
    use w90_parameters, only : fermi_energy    ! Fermi energy (in eV)
    use w90_parameters, only : write_hr        ! are we going to dump in
                                               !   a file the matrix 
                                               !   elements of the 
                                               !   Hamiltonian? 
    use w90_parameters, only : write_tb        ! are we going to dump in
                                               !   a file the matrix 
                                               !   elements of the 
                                               !   Hamiltonian and position
                                               !   operator?
    use w90_parameters,  only : eigval         ! Eigenvalues of the 
                                               !   Hamiltonian (in eV)
    use w90_parameters,  only : lsitesymmetry  ! Symmetry-adapted 
                                               !   Wannier functions
    use w90_parameters,  only : transport      ! Transport calculation? 
    use w90_parameters,  only : tran_read_ht   ! Read the Hamiltonian for 
                                               !   transport calculation?

!
!   Variables related with the parallelization, coming from SIESTA
!
    use parallel,        only: Node            ! ID of this node
    use parallel,        only: Nodes           ! number of nodes
    use parallel,        only: IOnode          ! input/output node?
!
!   Variables related with the parallelization, coming from WANNIER90
!
    use w90_comms,       only: on_root         ! are we the root node?
    use w90_comms,       only: num_nodes       ! number of nodes
    use w90_comms,       only: my_node_id      ! ID of this node 
!
!   Subroutines coming from WANNIER90 (version 3.0.0) that will be called
!   from SIESTA 
!
    use w90_parameters,  only: param_read      ! Subroutine to read the 
                                               !    parameters required by
                                               !    WANNIER90 and populate 
                                               !    derived values
    use w90_parameters,  only: param_write_header  
                                               ! write a suitable header for the
                                               !    calculation 
                                               !    (version authors etc)
    use w90_parameters,  only: param_write     ! write wannier90 parameters
                                               !    to stdout 
    use w90_parameters,  only: param_dist      ! distribute the parameters 
                                               !    across processors 
    use w90_parameters,  only: param_write_chkpt
                                               ! write checkpoint file
    use w90_io,          only: io_time         ! subroutine to control the
                                               !    timing, borrowed from 
                                               !    WANNIER90
    use w90_io,          only: io_error        ! abort the code giving an 
                                               !    error message
    use w90_io,          only: io_print_timings! output timing information 
                                               !    to stdout
    use w90_io,          only: io_date         ! returns two strings containing
                                               !    the date and the time
                                               !    in human-readable format. 
                                               !    Uses a standard f90 call.
    use w90_overlap,     only: overlap_allocate  
                                               ! allocate memory to read Mmn
                                               !    and Amn from files
    use w90_overlap,     only: overlap_read    ! read the Mmn and Amn 
                                               !    from files
                                               !    and Amn from files
    use w90_wannierise,  only: wann_main       ! subroutine that calculates
                                               !    the Unitary Rotations
                                               !    to give Maximally 
                                               !    Localized Wannier 
                                               !    Functions.
    use w90_wannierise,  only: wann_main_gamma ! subroutine that calculates
                                               !    the Unitary Rotations
                                               !    to give Maximally 
                                               !    Localized Wannier 
                                               !    Functions.
                                               !    Gamma version.
    use w90_disentangle, only: dis_main        ! main disentanglement 
                                               !    routine
    use w90_plot,        only: plot_main       ! main plotting routine 
                                               !    of quantities related
                                               !    with WANNIER90
                                               !    (bands, Fermi surface)
    use w90_transport,   only: tran_main       ! main tranport routine in
                                               !    WANNIER90.
    use w90_sitesym,     only: sitesym_read    ! read the variables to impose
                                               !    the site symmetry during
                                               !    minimization of the spread
    use w90_comms,       only: comms_bcast     ! send integar array from 
                                               !    root node to all nodes 

!
! Input variables
!
    integer, intent(in) :: ispin            ! Spin index
    integer, intent(in) :: index_manifold   ! Index of the manifold 
                                            !   that is wannierized
!
! Internal variables
!
    integer :: nn
    integer :: n
    integer :: nkp              ! Counter for loop on k-points
    integer :: ik               ! Counter for loop on k-points
    integer :: nb               ! Counter for loop on bands
    integer :: m                ! Counter for loop on bands
    integer :: nsp              ! Counter for loop on species
    integer :: nat              ! Counter for loop on atoms
    integer :: counter 
    integer :: i, j
    integer :: max_sites        ! Maximum number of atomic species
    character(len=50) :: prog   ! Name of the program
    integer :: len_seedname
    character(len=9) :: cdate, ctime

! 
!   Variables to control the timing 
!
    real(kind=dp) time0
    real(kind=dp) time1
    real(kind=dp) time2

!   Set up the variables related with the writing of the Hamiltonian
    seedname=repeat(" ",len(seedname))
    if( spin%H .eq. 1) then
      write(seedname,"(a,'.manifold.',i1.1)")                    &
 &         trim(slabel),index_manifold
    else if( spin%H .gt. 1) then
      write(seedname,                                            &
 &          "(a,'.manifold.',i1.1,'.spin.',i1.1)")               &
 &       trim(slabel),index_manifold,ispin
    endif

!   Set up whether the Hamiltonian and tight-binding matrix elements will
!   be written in files
    write_hr = manifold_bands_w90_in(index_manifold)%write_hr
    write_tb = manifold_bands_w90_in(index_manifold)%write_tb

!   Set up general variables
    num_bands = manifold_bands_w90_in(index_manifold)%number_of_bands
    num_wann  = manifold_bands_w90_in(index_manifold)%numbands_w90_in
    num_proj  = num_wann

!   Set up the variables related with the structure
    num_atoms   = na_u
    num_species = nspecies

!   Define the atomic symbols and the atomic_labels
    if ( allocated(atoms_symbol)     ) deallocate(atoms_symbol)
    if ( allocated(atoms_label)      ) deallocate(atoms_label)
    allocate(atoms_symbol(num_species))
    allocate(atoms_label(num_species))
    do nsp = 1, num_species
      atoms_label(nsp)     = trim(species(nsp)%label)
      atoms_symbol(nsp)    = trim(species(nsp)%symbol)
    end do

    if ( allocated(atoms_species_num) ) deallocate(atoms_species_num)
    allocate(atoms_species_num(num_species))
    atoms_species_num(:)=0

    do nsp = 1, num_species
       do nat = 1, num_atoms
          if( trim(atoms_label(nsp))==trim(species(isa(nat))%label)) then
!!         For debugging
!           write(6,'(2i5,a20)') nsp, nat, atoms_label(nsp)
!!         End debugging
             atoms_species_num(nsp)=atoms_species_num(nsp)+1
          end if
       end do
    end do

!!   For debugging
!    write(6,*) 'atoms_label       = ', Node, atoms_label
!    write(6,*) 'atoms_species_num = ', Node, atoms_species_num
!!   End debugging

    max_sites=maxval(atoms_species_num)
    if ( allocated(atoms_pos_cart) )   deallocate(atoms_pos_cart)
    allocate(atoms_pos_cart(3,max_sites,num_species))

    do nsp = 1, num_species
       counter=0
       do nat = 1, num_atoms
          if( trim(atoms_label(nsp))==trim( species(isa(nat))%label)) then
             counter=counter+1
!             atoms_pos_cart(:,counter,nsp) = xa_last(:,nat)/Ang
             atoms_pos_cart(:,counter,nsp) = xa(:,nat)/Ang
          end if
       end do
    end do

!!   For debugging
!    do nsp = 1, nspecies
!      do nat = 1, max_sites
!        write(6,*)'Node, atoms_pos_cart = ', Node, atoms_pos_cart(:,nat,nsp)
!      enddo
!    enddo 
!    write(6,*)'Node, num_iter   = ', Node, num_iter
!    write(6,*)'Node, numkpoints = ', Node, numkpoints
!    write(6,*)'Node, nncount    = ', Node, nncount
!!   End debugging
  
!   Set up number of iterations for the minimization of \Omega
    num_iter  = manifold_bands_w90_in(index_manifold)%num_iter
    timing_level      = 1

!   Set up the variables related with the k-point sampling
    mp_grid   = kmesh_w90_in
    num_kpts  = numkpoints
    nntot     = nncount
    if ( allocated(nnlist) ) deallocate(nnlist)
    allocate(nnlist(num_kpts,nntot))
    nnlist    = nnlist_neig
    if ( allocated(nncell) ) deallocate(nncell)
    allocate(nncell(3,num_kpts,nntot))
    nncell = nnfolding

!!   For debugging
!    write(6,*)'Node, num_kpts = ', Node, num_kpts
!    write(6,*)'Node, nntot    = ', Node, nntot
!!   End debugging
  
!   Set up the variables related with the disentanglement
    disentanglement = manifold_bands_w90_in(index_manifold)%disentanglement
    dis_win_min     = manifold_bands_w90_in(index_manifold)%dis_win_min
    dis_win_max     = manifold_bands_w90_in(index_manifold)%dis_win_max
    dis_froz_min    = manifold_bands_w90_in(index_manifold)%dis_froz_min
    dis_froz_max    = manifold_bands_w90_in(index_manifold)%dis_froz_max

!   Set up the variables for post-processing
    wannier_plot    = manifold_bands_w90_in(index_manifold)%wannier_plot
    wannier_plot_supercell(1) =                                 &
 &     manifold_bands_w90_in(index_manifold)%wannier_plot_supercell(1)
    wannier_plot_supercell(2) = wannier_plot_supercell(1)     
    wannier_plot_supercell(3) = wannier_plot_supercell(1)     

    fermi_surface_plot =manifold_bands_w90_in(index_manifold)%fermi_surface_plot
    fermi_energy       = ef / eV

    write_hr           = manifold_bands_w90_in(index_manifold)%write_hr
    write_tb           = manifold_bands_w90_in(index_manifold)%write_tb

!   Store the eigenvalues of the Hamiltonian in the array that will be passed
!   to WANNIER90
    if ( allocated(eigval) ) deallocate(eigval)
    allocate( eigval(num_bands,num_kpts) )
    eigval = eo

!
!   Populate variables related with the parallelization
!
    num_nodes  = Nodes
    my_node_id = Node
    on_root    = IOnode

!!   For debugging
!    write(6,*)'my_node_id, num_nodes, on_root = ',  &
! &             my_node_id, num_nodes, on_root 
!!   End debugging
  

!   From this line till the end of the subroutine, it is a copy 
!   verbatim of the WANNIER90 main program

  time0 = io_time()

  if (on_root) then
    prog = 'wannier90'
    len_seedname = len(seedname)
  end if
  call comms_bcast(len_seedname, 1)
  call comms_bcast(seedname, len_seedname)

  if (on_root) then
    call io_date(cdate, ctime)
    if (.not. first_chempotwann) &
 &     write (stdout, *) 'Wannier90: Execution started on ', cdate, ' at ', ctime

    call param_read
    if( index_manifold .eq. 1 .and. (.not. first_chempotwann) ) &
 &      call param_write_header()
    if (num_nodes == 1) then
#ifdef MPI
    if (.not. first_chempotwann) &
 &    write (stdout, '(/,1x,a)') 'Running in serial (with parallel executable)'
#else
    if (.not. first_chempotwann) &
 &    write (stdout, '(/,1x,a)') 'Running in serial (with serial executable)'
#endif
    else
    if (.not. first_chempotwann) &
 &    write (stdout, '(/,1x,a,i3,a/)') &
        'Running in parallel on ', num_nodes, ' CPUs'
    endif
    if( index_manifold .eq. 1 .and. (.not. first_chempotwann) ) &
 &       call param_write()

    time1 = io_time()
    if (.not. first_chempotwann) &
 &    write (6, '(1x,a25,f11.3,a)') 'Time to read parameters  ', time1 - time0, ' (sec)'
  endif

  ! We now distribute the parameters to the other nodes
  call param_dist

  if (gamma_only .and. num_nodes > 1) &
    call io_error('Gamma point branch is serial only at the moment')

  if (transport .and. tran_read_ht) goto 3003

  if (lsitesymmetry) call sitesym_read()   ! update this to read on root and bcast - JRY

  call overlap_allocate()

  call overlap_read()

!! For debugging
!  call MPI_barrier(MPI_Comm_world,i)
!  call die()
!! End debugging

  time1 = io_time()
  if (on_root) then
    if (.not. first_chempotwann) &
 &    write (stdout, '(/1x,a25,f11.3,a)') 'Time to read overlaps    ', time1 - time2, ' (sec)'
  endif 

  have_disentangled = .false.

  if (disentanglement) then
    call dis_main()
    have_disentangled = .true.
    time2 = io_time()
    if (on_root) then
      if (.not. first_chempotwann) &
 &      write (stdout, '(1x,a25,f11.3,a)') 'Time to disentangle bands', time2 - time1, ' (sec)'
    endif 
  endif

  if (on_root) call param_write_chkpt('postdis')

1001 time2 = io_time()

  if (.not. gamma_only) then
    call wann_main()
  else
    call wann_main_gamma()
  end if

  time1 = io_time()
  if (on_root) then
    if (.not. first_chempotwann) &
 &    write (stdout, '(1x,a25,f11.3,a)') 'Time for wannierise      ', time1 - time2, ' (sec)'
  endif 

  if (on_root) call param_write_chkpt('postwann')

2002 continue
  if (on_root) then
    ! I call the routine always; the if statements to decide if/what
    ! to plot are inside the function
    time2 = io_time()
    call plot_main()
    time1 = io_time()
    ! Now time is always printed, even if no plotting is done/required, but
    ! it shouldn't be a problem.
    if (.not. first_chempotwann) &
 &    write (stdout, '(1x,a25,f11.3,a)') 'Time for plotting        ', time1 - time2, ' (sec)'
  endif

3003 continue
  if (on_root) then
    time2 = io_time()
    if (transport) then
      call tran_main()
      time1 = io_time()
      if (.not. first_chempotwann) &
 &      write (stdout, '(1x,a25,f11.3,a)') 'Time for transport       ', time1 - time2, ' (sec)'
      if (tran_read_ht) goto 4004
    end if
  endif

4004 continue

  if (on_root) then
    if (.not. first_chempotwann) &
 &    write (stdout, '(1x,a25,f11.3,a)') 'Total Execution Time     ', io_time(), ' (sec)'

    if (timing_level > 0) then
      if (.not. first_chempotwann) &
  &      call io_print_timings()
    endif 

    if (.not. first_chempotwann) then
      write (stdout, *)
      write (stdout, '(1x,a)') 'All done: wannier90 exiting'
    endif 

  endif

!!   For debugging
!    write(6,'(a,i5)')                            &
! &    'compute_wannier: num_bands       = ', num_bands
!    write(6,'(a,i5)')                            &
! &    'compute_wannier: num_wann        = ', num_wann
!    write(6,'(a,f12.5)')                         &
! &    'compute_wannier: lenconfac       = ', lenconfac
!    write(6,'(a,i5)')                            &
! &    'compute_wannier: num_kpts        = ', num_kpts
!    write(6,'(a,i5)')                            &
! &    'compute_wannier: nntot           = ', nntot
!    write(6,'(1x,a)') '|            No.         b_k(x)      b_k(y)      b_k(z)        w_b           |'
!    write(6,'(1x,a)') '|            ---        --------------------------------     --------        |'
!    do i = 1, nntot
!       write (6,'(1x,"|",11x,i3,5x,3f12.6,3x,f10.6,8x,"|")') &
!            i,(bk(j,i,1)/lenconfac,j=1,3),wb(i)*lenconfac**2
!    enddo
!    ! Nearest neighbour k-points
!    write(6,'(a)') 'begin nnkpts'
!    write(6,'(i4)') nntot
!    do nkp=1,num_kpts
!       do nn=1,nntot
!          write(6,'(2i6,3x,3i4)') &
!               nkp,nnlist(nkp,nn),(nncell(i,nkp,nn),i=1,3)
!       end do
!    end do
!    write(6,'(a/)') 'end nnkpts'
!    do nkp=1,num_kpts
!!       counter=0
!       do nb=1,num_bands
!!          if (lwindow(nb,nkp)) then
!!             counter=counter+1
!!             summ=0.0_dp
!!             do nw=1,num_wann
!!                summ=summ+abs(u_matrix_opt(counter,nw,nkp))**2
!!             enddo
!!             write(6,'(1x,16x,i5,1x,i5,1x,f14.6,2x,f14.8)') &
!!                  nkp,nb,eigval(nb,nkp),summ
!!          endif
!             write(6,'(1x,16x,i5,1x,i5,1x,f14.6)') &
!                  nkp,nb,eigval(nb,nkp)
!       enddo
!    enddo
!    write(6,'(1x,a78/)') repeat('-',78)
!
!    write(6,'(a,l5)')                            &
! &    'compute_wannier: disentanglement = ', disentanglement
!    write(6,'(a,f12.5)')                         &
! &    'compute_wannier: dis_win_min     = ', dis_win_min
!    write(6,'(a,f12.5)')                         &
! &    'compute_wannier: dis_win_max     = ', dis_win_max
!    write(6,'(a,f12.5)')                         &
! &    'compute_wannier: dis_froz_min    = ', dis_froz_min
!    write(6,'(a,f12.5)')                         &
! &    'compute_wannier: dis_froz_max    = ', dis_froz_max
!
!    do nsp=1,nspecies
!      do nat=1,atoms_species_num(nsp)
!        write(6,'(a,a2,3x,3f12.7)')                                 &
! &        'compute_wannier: atoms_symbol, atoms_pos_cart = ',       &
! &        atoms_symbol(nsp),(atoms_pos_cart(:,nat,nsp))
!      end do
!    end do
!!   End debugging

    return
  end subroutine compute_wannier

!> \brief General purpose of the define_trial_orbitals subroutine:
!! populate the derived variable where all the
!! information regarding the trial localized functions will be stored.
!! 
!! In this subroutine we populate the derived variable where all the 
!! information regarding the trial localized functions will be stored
!! This variable, called proj_w90_in, is defined within the derived variable
!! manifold_bands_w90_in.
!! The dimension of proj_w90_in should be the same as the number
!! of bands to be wannierized
!! Some default values (x-axis, z-axis, zovera, r) to define the 
!! hydrogenoid functions are taken directly from WANNIER90.
!! Indeed, this is not a major issue since they will not be used
!! in SIESTA, where we will take directly the shape of the numerical 
!! atomic orbitals in the basis set.
!! The center of the localized projection functions are directly given in 
!! cartesian Bohrs, and not in fractional units as it is done in WANNIER90.
!! The m-angular quantum number index is different in SIESTA and WANNIER90.
!! In WANNIER90 they are defined as in Tables 3.1 and 3.2 of the User Guide.
!! In SIESTA, they are defined as in the SystemName.ORB_INDX file.

  subroutine define_trial_orbitals( i_manifold )
    use trialorbitalclass,  only: trialorbital  
                                      ! Derived type to define the
                                      !    localized trial
                                      !    orbitals
    use atmfuncs,    only: lofio      ! Returns angular momentum number
    use atmfuncs,    only: mofio      ! Returns magnetic quantum number
    use atmfuncs,    only: rcut       ! Returns orbital cutoff radius
    use atmfuncs,    only: orb_gindex ! Returns the global index of a 
                                      !    basis orbital
    use atomlist,    only: iaorb      ! Atomic index of each orbital
    use atomlist,    only: iphorb     ! Orbital index of each orbital 
                                      !    in its atom
    use siesta_geom, only: xa         ! Atomic coordinates
    use siesta_geom, only: isa        ! Species index of each atom

    integer,  intent(in)  :: i_manifold
    integer  :: number_projections
    integer  :: iproj        ! Counter for the number of projections
    integer  :: iorb         ! Index of the atomic orbital 
    integer  :: ia           ! Pointer to atom to which orbital belongs
    integer  :: iao          ! Orbital index within atom
    integer  :: is           ! Species index
    real(dp) :: rc           ! Orbital's cutoff radius
    integer  :: l            ! Angular momentum of the localized trial function
    integer  :: m            ! Magnetic quantum number of the localized trial
                             !   function. 
                             ! The corresponding to the real spherical harmonics
                             !   can be found in Table 3.1 and 3.2 of the
                             !   WANNIER90 Users Guide
    integer  :: r            ! radial (optional): 
                             ! If r=2 – use a radial function with one node 
                             !   (ie second highest pseudostate with that 
                             !   angular momentum). 
                             !   Default is r=1. 
                             !   Radial functions associated with different 
                             !   values of r should be orthogonal to each other.
    real(dp) :: xaxis(3)     ! Sets the x-axis direction. Default is x=(1,0,0)
    real(dp) :: yaxis(3)     ! Sets the y-axis direction. Default is y=(0,1,0)
    real(dp) :: zaxis(3)     ! Sets the z-axis direction. Default is z=(0,0,1)
    real(dp) :: zovera       ! zovera (optional):
                             !   the value of Z for the radial part 
                             !   of the atomic orbital 
                             !   (controls the diffusivity of the radial 
                             !   function). 
                             !   Units always in reciprocal Angstrom. 
                             !   Default is zovera = 1.0.

!   Sets the default for the x-axis
    xaxis(1) = 1.0_dp
    xaxis(2) = 0.0_dp
    xaxis(3) = 0.0_dp

!   Sets the default for the y-axis
    yaxis(1) = 0.0_dp
    yaxis(2) = 1.0_dp
    yaxis(3) = 0.0_dp

!   Sets the default for the z-axis
    zaxis(1) = 0.0_dp
    zaxis(2) = 0.0_dp
    zaxis(3) = 1.0_dp

!   Sets the default for zovera
    zovera = 1.0_dp
!    zovera = 1.0_dp/0.529177_dp

!   Set up the default value for r (use radial functions without node)
    r = 1

!   There should be as many localized trial orbitals as the number of bands
!   to be projected
    number_projections = manifold_bands_w90_in(i_manifold)%numbands_w90_in

!   Allocate the derived variable where all the data for the 
!   localized trial orbitals will be stored
    if (allocated(manifold_bands_w90_in(i_manifold)%proj_w90_in)) &
 &     deallocate( manifold_bands_w90_in(i_manifold)%proj_w90_in )
    allocate(manifold_bands_w90_in(i_manifold)%proj_w90_in(number_projections))

    do iproj = 1, number_projections
!     Identify the atomic orbitals on which we are going to project
      iorb = manifold_bands_w90_in(i_manifold)%orbital_indices(iproj)
      ia = iaorb(iorb)               ! Atom to which orbital belongs
      is = isa(ia)                   ! Atomic species of the atom where the
                                     !   orbital is centered
      iao = iphorb( iorb )           ! Orbital index within atom
      l  = lofio( is, iao )          ! Orbital's angular mumentum number
      m  = mofio( is, iao )          ! (Real) orbital's magnetic quantum number
      rc = rcut(  is, iao )          ! Orbital's cutoff radius
      manifold_bands_w90_in(i_manifold)%proj_w90_in(iproj)%center = xa(:,ia)
      manifold_bands_w90_in(i_manifold)%proj_w90_in(iproj)%zaxis  = zaxis
      manifold_bands_w90_in(i_manifold)%proj_w90_in(iproj)%xaxis  = xaxis
      manifold_bands_w90_in(i_manifold)%proj_w90_in(iproj)%yaxis  = yaxis
      manifold_bands_w90_in(i_manifold)%proj_w90_in(iproj)%zovera = zovera
      manifold_bands_w90_in(i_manifold)%proj_w90_in(iproj)%r      = r
      manifold_bands_w90_in(i_manifold)%proj_w90_in(iproj)%l      = l
      manifold_bands_w90_in(i_manifold)%proj_w90_in(iproj)%mr     = m
      manifold_bands_w90_in(i_manifold)%proj_w90_in(iproj)%rcut   = rc
      manifold_bands_w90_in(i_manifold)%proj_w90_in(iproj)%lmax   = l
      manifold_bands_w90_in(i_manifold)%proj_w90_in(iproj)%from_basis_orbital =&
 &                                                         .true.
      manifold_bands_w90_in(i_manifold)%proj_w90_in(iproj)%iorb   = iorb
      manifold_bands_w90_in(i_manifold)%proj_w90_in(iproj)%iorb_gindex = &
 &                                     orb_gindex(is,iao) 
    enddo

!!   For debugging
!    write(6,'(a)') 'begin projections'
!    write(6,'(i6)') number_projections
!    do iproj = 1, number_projections
!      write(6,'(3(f10.5,1x),2x,3i3)') &
! &      manifold_bands_w90_in(i_manifold)%proj_w90_in(iproj)%center(1), &
! &      manifold_bands_w90_in(i_manifold)%proj_w90_in(iproj)%center(2), &
! &      manifold_bands_w90_in(i_manifold)%proj_w90_in(iproj)%center(3), &
! &      manifold_bands_w90_in(i_manifold)%proj_w90_in(iproj)%l,         &
! &      manifold_bands_w90_in(i_manifold)%proj_w90_in(iproj)%mr,        &
! &      manifold_bands_w90_in(i_manifold)%proj_w90_in(iproj)%r   
!      write(6,'(2x,3f11.7,1x,3f11.7,1x,f7.2)') &
! &      manifold_bands_w90_in(i_manifold)%proj_w90_in(iproj)%zaxis(1),  &
! &      manifold_bands_w90_in(i_manifold)%proj_w90_in(iproj)%zaxis(2),  &
! &      manifold_bands_w90_in(i_manifold)%proj_w90_in(iproj)%zaxis(3),  &
! &      manifold_bands_w90_in(i_manifold)%proj_w90_in(iproj)%xaxis(1),  &
! &      manifold_bands_w90_in(i_manifold)%proj_w90_in(iproj)%xaxis(2),  &
! &      manifold_bands_w90_in(i_manifold)%proj_w90_in(iproj)%xaxis(3),  &
! &      manifold_bands_w90_in(i_manifold)%proj_w90_in(iproj)%zovera
!    enddo
!    write(6,'(a)') 'end projections'
!!   End debugging

  end subroutine define_trial_orbitals

!> \brief General purpouse of the subroutine deallocate_wannier:
!! deallocate some of the arrays that were allocated during the WANNIER90 run
!!

  subroutine deallocate_wannier

    use w90_parameters,  only : lsitesymmetry       ! Symmetry-adapted 
                                                    !   Wannier functions
! 
! Deallocation routines
!
    use w90_parameters,  only: param_dealloc
    use w90_overlap,     only: overlap_dealloc
    use w90_sitesym,     only: sitesym_dealloc
    use w90_kmesh,       only: kmesh_dealloc
    use w90_transport,   only: tran_dealloc
    use w90_hamiltonian, only: hamiltonian_dealloc

    call tran_dealloc()
    call hamiltonian_dealloc()
    call overlap_dealloc()
    call kmesh_dealloc()
    call param_dealloc()
    if (lsitesymmetry) call sitesym_dealloc() !YN:

  end subroutine deallocate_wannier

#endif
endmodule m_w90_in_siesta
