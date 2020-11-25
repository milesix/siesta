! 
! Copyyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!

!  Interface module to call BigDFT solver instead of standar Siesta 
!  poison subroutine. It requires compilation of Siesta with additional
!  libraries.
!  More info at:
!  http://bigdft.org/Wiki/index.php?title=BigDFT_website 
! 
!  Created by Pablo Lopez-Tarifa and Daniel Sanchez Portal @CFM(2018)
!  Optimizations and made parallel by Nick Papior (2020)
module psolver_m

#ifdef SIESTA__PSOLVER
  use precision, only : dp, grid_p

  implicit none

  private
  save

  ! Options related to the PSolver algorithms
  logical, public :: use_psolver = .false. ! Use Poisson solver from BigDFT package instead of Siesta solver
  !< Use implicit solvent effects
  logical, public :: psolver_cavity = .false.
  !< Order of the interpolating scaling function family
  integer, public :: psolver_isf_order = 16
  !< Order of the finite-difference derivatives for the GPS solver
  integer, public :: psolver_fd_order = 16
  !< Ask PSolver to be verbose
  logical, public :: psolver_verbose = .false.
  !< Amplitude of the transition region in the rigid cavity [a.u.]
  real(dp), public :: psolver_delta = 2.0_dp
  !< Multiply factor for the whole rigid cavity
  real(dp), public :: psolver_fact_rigid = 1.12_dp
  !< Di-electric constant of the exterior region
  real(dp), public :: psolver_epsilon = 78.36_dp
  !< Cavitation term, surface tension of the solvent [dyn/cm]
  real(dp), public :: psolver_gammaS = 72.0_dp
  !< Proportionality of repulsion free energy in term of the surface integral [dyn/cm]
  real(dp), public :: psolver_alphaS = -22.0_dp
  !< Proportionality of dispersion free energy in term of volume integral [GPa].
  real(dp), public :: psolver_betaV = -0.35_dp
  !< Mapping of the radii that have to be used for each atomic species.
  integer, public :: psolver_atomic_radii = 0
  !< Define type of cavity: [none, soft-sphere, SCCS]
  character(len=16), public :: psolver_cavity_type = 'none'
  !< Type of radii type: [UFF, Bondi, Pauling]
  character(len=16), public :: psolver_radii_type = 'UFF'
  !< Algorithm used in the generalized Poisson equation: [PCG, PI]
  character(len=16), public :: psolver_gps_algorithm = 'PCG'
  !< Which kind of accelerator to use: [CUDA]
  character(len=16), public :: psolver_accel = 'none'
  !< Force PBC regardless of unit-cell
  logical, public :: psolver_force_pbc = .false.

  !< Offsets in grids for lattices
  !!
  !! The individual offsets in the grids to ensure PSolver
  !! method is obeyed. I.e. for molecules it is important that
  !! the atomic center of positions is at the center of the grid.
  !! Using offsets we can rotate values.
  integer, private :: lat_offset(3) = 0

  ! NOTE
  ! Currently we do not allow any permutation or offset
  ! in the third lattice vector. This is because this direction
  ! is the distributed direction and is currently complicated.
  ! It could be done, relatively easily, but currently
  ! we don't do it.

  public :: poisson_psolver_options
  public :: poisson_psolver_options_print
  public :: poisson_psolver_init
  public :: poisson_psolver

contains

  !< Initialize PSolver input variables from FDF
  subroutine poisson_psolver_options()

    use fdf

    character(len=80) :: ctmp
    type(block_fdf) :: bfdf
    type(parsed_line), pointer :: pline => null()

    ! Whether we should use the PSolver method
    ctmp = fdf_get('Poisson.Method', 'fft')
    use_psolver = leqi(ctmp, 'psolver')

    ! Quick return
    if ( .not. use_psolver ) return

    if ( dp /= grid_p ) call die('Error in compilation, please ensure -DGRID_DP')

    ! Parse the PSolver block
    ! If it does not exist, then use default values
    if ( .not. fdf_block('PSolver', bfdf) ) return

    do while ( fdf_bline(bfdf, pline) )

      if ( fdf_bnnames(pline) == 0 ) cycle

      ! Retrieve option
      ctmp = fdf_bnames(pline, 1)

      if ( leqi('isf.order', ctmp) ) then
        psolver_isf_order = fdf_bintegers(pline, 1)

      else if ( leqi('fd.order', ctmp) ) then
        psolver_fd_order = fdf_bintegers(pline, 1)

      else if ( leqi('cavity', ctmp) ) then
        psolver_cavity_type = fdf_bnames(pline, 2)
        psolver_cavity = .not. leqi('none', psolver_cavity_type)

      else if ( leqi('gps.algorithm', ctmp) ) then
        psolver_gps_algorithm = fdf_bnames(pline, 2)

      else if ( leqi('radii.type', ctmp) ) then
        psolver_radii_type = fdf_bnames(pline, 2)

      else if ( leqi('delta', ctmp) ) then
        psolver_delta = fdf_bvalues(pline, 1)

      else if ( leqi('epsilon', ctmp) ) then
        psolver_epsilon = fdf_bvalues(pline, 1)

      else if ( leqi('fact.rigid', ctmp) ) then
        psolver_fact_rigid = fdf_bvalues(pline, 1)

      else if ( leqi('gammaS', ctmp) ) then
        psolver_gammaS = fdf_bvalues(pline, 1)

      else if ( leqi('alphaS', ctmp) ) then
        psolver_alphaS = fdf_bvalues(pline, 1)

      else if ( leqi('betaV', ctmp) ) then
        psolver_betaV = fdf_bphysical(pline, -0.35_dp, 'GPa')

      else if ( leqi('atomic.radii', ctmp) ) then
        psolver_atomic_radii = fdf_bintegers(pline, 1)

      else if ( leqi('verbose', ctmp) ) then
        if ( fdf_bnnames(pline) > 1 ) then
          psolver_verbose = fdf_bboolean(pline, 2)
        else
          psolver_verbose = .true.
        end if

      else if ( leqi('accelerator', ctmp) .or. leqi('accel', ctmp) ) then
        psolver_accel = fdf_bnames(pline, 2)

      else if ( leqi('force.PBC', ctmp) ) then
        if ( fdf_bnnames(pline) > 1 ) then
          psolver_force_pbc = fdf_bboolean(pline, 2)
        else
          psolver_force_pbc = .true.
        end if

      end if

    end do

  end subroutine poisson_psolver_options

  !< Print out the options used in PSolver
  subroutine poisson_psolver_options_print()

    use parallel, only: IONode

    ! Quick return
    if ( .not. IONode ) return

    if ( .not. use_psolver ) then
      write(6,'(a,t53,"= ",a)') 'Poisson solver', 'FFT'
      return
    end if

    write(6,'(a,t53,"= ",a)') 'Poisson solver', 'PSolver'
    write(6,'(a,t53,"= ",i3)') 'PSolver interpolating scaling order', psolver_isf_order
    write(6,'(a,t53,"= ",i3)') 'PSolver finite difference order', psolver_fd_order

    write(6,'(a,t53,"= ",a)') 'PSolver cavity type', trim(psolver_cavity_type)
    if ( psolver_cavity ) then
      write(6,'(a,t53,"= ",a)') 'PSolver generalized Poisson equation algorithm', trim(psolver_gps_algorithm)
      write(6,'(a,t53,"= ",f10.5)') 'PSolver transition distance', psolver_delta
      write(6,'(a,t53,"= ",f10.5)') 'PSolver rigid cavity factor', psolver_fact_rigid
      write(6,'(a,t53,"= ",f10.5)') 'PSolver exterior dielectric constant', psolver_epsilon
      write(6,'(a,t53,"= ",f10.5)') 'PSolver surface tension of cavity', psolver_gammaS
      write(6,'(a,t53,"= ",f10.5)') 'PSolver free energy repulsion of surface', psolver_alphaS
      write(6,'(a,t53,"= ",f10.5)') 'PSolver free energy dispersion of volume', psolver_betaV
      write(6,'(a,t53,"= ",a)') 'PSolver atomic radii type', trim(psolver_radii_type)
      write(6,'(a,t53,"= ",i0)') 'PSolver atomic radii', psolver_atomic_radii
    end if
    write(6,'(a,t53,"=   ",L1)') 'PSolver verbose (for debugging)', psolver_verbose
    write(6,'(a,t53,"=   ",L1)') 'PSolver force PBC', psolver_force_pbc

  end subroutine poisson_psolver_options_print

  !< Initialize grids related to the Poisson solution using PSolver
  !!
  !! PSolver requires a different potential distribution than
  !! Siesta and we need to create routines for re-arranging them.
  subroutine poisson_psolver_init(cell, na_u, xa, nbcell, bcell, nm, ntm)

    use parallel, only: Node, Nodes
#ifdef MPI
    use mpi_siesta
#endif

    use intrinsic_missing, only: VEC_PROJ_SCA, VNORM
    use moremeshsubs, only: initMeshDistr, UNIFORM, UNIFORM_PSOLVER
    use poisson_solver, only: PS_dim4allocation

    integer, intent(in) :: nbcell, na_u
    real(dp), intent(in) :: xa(3,na_u)
    real(dp), intent(in) :: cell(3,3), bcell(3,nbcell)
    !< Total *big* points
    integer, intent(in) :: nm(3)
    !< Total mesh points (including sub-divisions) nm == ntm / nsm
    integer, intent(in) :: ntm(3)

    ! Local variables
    integer :: nz_d, nz_d_cum, nz_d_vion, nz_d_xc_offset, nz_d_start
    integer :: start_ntm(3), end_ntm(3)
    real(dp) :: sca_proj(3), cop(3), l_cell(3)
    real(dp), parameter :: ortho_tol = 0.0000001_dp
    character(len=1) :: boundary_type
#ifdef MPI
    logical :: check_subl, check_sub
    integer :: ierr
#endif

    ! Initialize permutations and offsets
    lat_offset(:) = 0

    ! Our first task is to figure out which type of boundary we are using.
    select case ( nbcell )
    case ( 0 ) ! molecule
      boundary_type = 'F'
    case ( 1 ) ! wire
      boundary_type = 'W'
    case ( 2 ) ! slab
      boundary_type = 'S'
    case ( 3 ) ! bulk
      boundary_type = 'P'
    end select
    if ( psolver_force_pbc ) boundary_type = 'P'

    ! Calculate center of positions
    cop(:) = sum(xa, dim=2) / na_u
    l_cell(1) = VNORM(cell(:,1))
    l_cell(2) = VNORM(cell(:,2))
    l_cell(3) = VNORM(cell(:,3))

    ! Then we have figure out how to correct the grid to ensure our
    ! boundaries are correct.
    ! We will allow a couple of cases
    select case ( boundary_type )
    case ( 'P' ) ! bulk
      ! We do not have to do anything. PSolver is compatible with this configuration
    case ( 'S' ) ! slab

      ! The open BC *has* to be along the 2nd lattice vector.
      ! So we have to check whether the bulk directions are along AC
      sca_proj(1) = abs(VEC_PROJ_SCA(bcell(:,1), cell(:,3)))
      sca_proj(2) = abs(VEC_PROJ_SCA(bcell(:,2), cell(:,3)))
      if ( all(sca_proj < ortho_tol) ) then
        write(6,'(/,a,/,a/,a,/,a,/,a,/,a)')                 &
            '#==================== WARNING PSolver ==================#',&
            ' Empty direction of the slab system has to be set        ',&
            ' along the 1st or 2nd lattice vector. NOT along 3rd.     ',&
            '#=======================================================#'

        call die('poisson_psolver_init: last lattice vector has to be periodic &
            &in slab calculations, please correct your geometry!')
      end if

      ! Figure out whether we have to swap 1 and 2
      sca_proj(1) = abs(VEC_PROJ_SCA(cell(:,2), bcell(:,1)))
      sca_proj(2) = abs(VEC_PROJ_SCA(cell(:,2), bcell(:,2)))
      if ( any(sca_proj > ortho_tol) ) then
        write(*,'(a)') 'PSolver lattice vector permutations'
        write(*,'(a)') '  A <-> B'
        call die('psolver: Siesta cannot convert the mesh-grid to be compatible with &
            &PSolver requirements. Please manually swap A-B lattice vectors &
            &to fulfil the above lattice permutation table.')
      end if

      ! Calculate integer offset
      ! First we project COP onto the vacuum region lattice
      ! vector.
      ! Then we want the projection point of the COP to
      ! be in the middle (0.5)
      sca_proj(1) = VEC_PROJ_SCA(cell(:,2), cop) / l_cell(2) - 0.5_dp

      ! Retrieve nearest integer offset
      lat_offset(2) = mod(nint(sca_proj(1) * ntm(2)), ntm(2))

    case ( 'W' )

      ! The open BC *has* to be along the 1st and 2nd lattice vector.
      ! So we have to check whether the bulk direction is along C.
      sca_proj(1) = abs(VEC_PROJ_SCA(bcell(:,1), cell(:,3)))
      if ( sca_proj(1) < ortho_tol ) then
        write(6,'(/,a,/,a/,a,/,a,/,a,/,a)')                 &
            '#==================== WARNING PSolver ==================#',&
            ' Empty direction of the wire system has to be set        ',&
            ' along the 1st and 2nd lattice vectors. NOT along 3rd.   ',&
            '#=======================================================#'

        call die('poisson_psolver_init: last lattice vector has to be periodic &
            &in wire calculations, please correct your geometry!')
      end if

      ! Calculate integer offset
      ! We do not need to offset grids
      ! First we project the COP onto the vacuum region lattice
      ! vector.
      ! Then we want the projection point of the COP to
      ! be in the middle (0.5)
      sca_proj(1) = VEC_PROJ_SCA(cell(:,1), cop) / l_cell(1) - 0.5_dp
      ! Retrieve nearest integer offset
      lat_offset(1) = mod(nint(sca_proj(1) * ntm(1)), ntm(1))
      sca_proj(2) = VEC_PROJ_SCA(cell(:,2), cop) / l_cell(2) - 0.5_dp
      lat_offset(2) = mod(nint(sca_proj(2) * ntm(2)), ntm(2))

    case ( 'F' )

      sca_proj(1) = VEC_PROJ_SCA(cell(:,1), cop) / l_cell(1) - 0.5_dp
      lat_offset(1) = mod(nint(sca_proj(1) * ntm(1)), ntm(1))
      sca_proj(2) = VEC_PROJ_SCA(cell(:,2), cop) / l_cell(2) - 0.5_dp
      lat_offset(2) = mod(nint(sca_proj(2) * ntm(2)), ntm(2))
      sca_proj(3) = VEC_PROJ_SCA(cell(:,3), cop) / l_cell(3) - 0.5_dp
      lat_offset(3) = mod(nint(sca_proj(3) * ntm(3)), ntm(3))

    end select

#ifdef MPI
    if ( Nodes > 1 ) then

      ! If moremeshsubs could be made *periodic* in its mesh
      ! interpretation, then we could add lat_offset to start_nm/end_nm
      ! and by that acheive the correct grid. However, moremeshsubs is
      ! not periodic.

      ! The current distributed version requires the mesh
      ! to be commensurate with PSolver's internal representation.
      ! In Siesta we have mesh sub-divisions.
      ! This is not something that is used in PSolver and thus
      ! there will be problems related to a distribution where
      ! PSolver expects ONE thing, but Siesta cannot provide routines
      ! for it.

      ! Figure out how the distribution in PSolver is required
      ! In our case we don't use PSolver for calculating XC, hence:
      !   grid3_d == grid3_d_v
      !   grid3_d == grid3_d_vion
      ! For *only* doing Poisson equation we *only* need:
      !   grid3_d, not any of the others.
      ! This is in the full grid (with mesh-sub divisions)
      call PS_dim4allocation(boundary_type, 'D', Node, Nodes, ntm(1), ntm(2), ntm(3), &
          .false., .false., 0, & !use_gradient, use_wb_corr, igpu
          nz_d_cum, nz_d, nz_d_vion, nz_d_xc_offset, nz_d_start)

      ! Initialize sub-divisioned data
      start_ntm(1:2) = 1
      start_ntm(3) = nz_d_start
      end_ntm(1:2) = ntm(1:2)
      end_ntm(3) = nz_d_start + nz_d - 1

      ! Create a new explicit mesh distribution for the PSolver library
      call initMeshDistr(UNIFORM_PSOLVER, start_ntm, end_ntm, nsm=1)

    end if
#endif

  end subroutine poisson_psolver_init

  !> \brief General purpose of the subroutine poisson_psolver
  !!
  !! This subroutine handles the full hartree potential calculation
  !! given by Psolver package of BigDFT. 
  !!
  !! The potential is calculated similarly in serial and parallel.
  !! To enable distributed grids we have to obey PSolver mesh distribution.
  !! This is not so easy so we currently advice users to do MeshSubdivisions 1
  !! to ensure PSolver mesh divisions is *always* commensurate with Siesta
  !! mesh.
  !! First the algorithm transfers data from the UNIFORM grid to
  !! the UNIFORM_PSOLVER which is another mesh distribution.
  !! Then we can compute the electrostatic potential using PSolver
  !! and finally we back-distribute the data from UNIFORM_PSOLVER
  !! to UNIFORM.
  !! It is important that the input density (rho) is given in SERIAL form.
  !!
  !> \brief Further improvements
  !! We are currently working on:
  !!
  !! Currently the stress-tensor is not calculated correctly.
  !! I do not know how to fix this but it is *not* a matter of units since
  !! I get a sign change.
  subroutine poisson_psolver(cell, na_u, xa, fal, ntm, rho, V, eh, stressl, calc_stress)

    use parallel, only: Node, Nodes
    use units,          only: eV, Ang
    use atmfuncs, only: floating
    use siesta_geom,    only: isa, shape
    use periodic_table, only: symbol
    use chemical, only: atomic_number
    use alloc,          only: re_alloc, de_alloc
    use intrinsic_missing, only: VNORM
    use moremeshsubs, only: setMeshDistr, distMeshData
    use moremeshsubs, only: UNIFORM, UNIFORM_PSOLVER, KEEP

    ! PSolver modules:
    use Poisson_Solver, only: coulomb_operator,  Electrostatic_Solver, &
        pkernel_set, pkernel_set_epsilon
    use PStypes, only: PSolver_energies, ps_soft_pcm_forces
    use PStypes, only: ps_set_options
    use PStypes, only: pkernel_init, pkernel_free, pkernel_get_radius
    use yaml_output
    use f_utils
    use f_enums, only: operator(.hasattr.)
    use dictionaries, dict_set => set

    ! PSolver uses atlab to handle geometry etc.
    use at_domain

    real(dp), intent(in) :: cell(3,3) ! Unit-cell vectors
    integer, intent(in) :: na_u ! number of atoms in unit cell
    real(dp), intent(in) :: xa(3,na_u) ! atomic coordinates
    real(dp), intent(inout) :: fal(3,na_u) !< Local atomic forces
    integer, intent(in) :: ntm(3) ! Total mesh divisions
    real(grid_p), intent(in) :: rho(:) ! Input density
    real(grid_p), intent(inout), target :: V(:) ! Output potential.
    real(dp), intent(out) :: eh ! Electrostatic energy.
    real(dp), intent(inout) :: stressl(3,3) ! stress-tensor contribution (local)
    logical, intent(in) :: calc_stress

    ! PSolver types
    type(dictionary), pointer :: dict => null() ! Input parameters.
    type(coulomb_operator) :: pkernel
    type(domain) :: dom
    type(PSolver_energies) :: energies

    real(grid_p), pointer :: Vaux(:) => null() ! 1D auxiliary array.

    real(dp), pointer :: xa_cavity(:,:) => null()
    real(dp), pointer :: radii_cavity(:) => null()
    real(dp) :: hgrid(3) ! Uniform mesh spacings in the three directions.
    real(dp) :: lat_doffset(3)
    integer :: ia
    real(dp) :: factor, stress_correction
    real(dp), external :: volcel

    integer :: p_nml(3), p_nnml, p_ntml(3), p_ntpl

    call timer('poisson_psolver', 1)

    ! f_util and dictionary initialisations:
    call f_lib_initialize()
    call dict_init(dict)

    ! Calculation of cell divisions
    hgrid(1) = VNORM(cell(:,1)) / ntm(1)
    hgrid(2) = VNORM(cell(:,2)) / ntm(2)
    hgrid(3) = VNORM(cell(:,3)) / ntm(3)

    ! Setting of dictionary variables, done at every calling:
    call set_variables(dom, cell, calc_stress, dict)

    if ( Nodes > 1 ) then

      call setMeshDistr(UNIFORM_PSOLVER, p_nml, p_nnml, p_ntml, p_ntpl)

      ! Allocate auxiliary array for correct size of grid
      call re_alloc(Vaux, 1, p_ntpl, 'Vaux', 'poisson_psolver')

      call distMeshData(UNIFORM, rho, UNIFORM_PSOLVER, Vaux, KEEP)

    else
      ! Copy Rho to V and make pointer
      ! Solution will come directly
      V(:) = Rho(:)
      Vaux => V(:)

      p_ntml(:) = ntm(:)

    end if

    ! Shift data to center molecule/slab in mesh
    call mesh_cshift(p_ntml, Vaux, -lat_offset)

    ! Initialize poisson-kernel
    pkernel = pkernel_init(Node, Nodes, dict, dom, ntm, hgrid)

    ! Free calling dictionary
    call dict_free(dict)

    ! Specify the verbosity
    call pkernel_set(pkernel, verbose=psolver_verbose)

    ! Implicit solvent?
    if ( psolver_cavity ) then
      call re_alloc(xa_cavity, 1, 3, 1, na_u, 'xa_cavity', 'poisson_psolver')
    end if

    if ( pkernel%method .hasattr. "rigid" ) then

      ! Calculate offset due to mesh shift
      lat_doffset(:) = 0._dp
      do ia = 1, 3
        lat_doffset(:) = lat_doffset(:) - (cell(:,ia) / ntm(ia)) * lat_offset(:)
      end do

      call re_alloc(radii_cavity, 1, na_u, 'radii_cavity', 'poisson_psolver')

      ! Set a radius for each atom in Angstroem (UFF, Bondi, Pauling, etc ...)
      ! Multiply for a constant prefactor (see the Soft-sphere paper JCTC 2017)
      ! The pkernel_get_radius returns in Ang, so we have to convert to Bohr
      factor = pkernel%cavity%fact_rigid * Ang
      do ia = 1 , na_u
        if ( floating(isa(ia)) ) then
          radii_cavity(ia) = 0._dp
        else
          radii_cavity(ia) = factor * &
              pkernel_get_radius(pkernel, atname=symbol(atomic_number(isa(ia))))
        end if
        xa_cavity(:,ia) = xa(:,ia) + lat_doffset(:)
      end do

      call pkernel_set_epsilon(pkernel, nat=na_u, &
          rxyz=xa_cavity, radii=radii_cavity)

      call de_alloc(radii_cavity, 'radii_cavity', 'poisson_psolver')

    end if

    ! Solve the Poisson equation and retrieve also energies and stress-tensor
    call Electrostatic_Solver(pkernel, Vaux, energies)

    if ( psolver_cavity ) then
      if ( calc_stress ) then
        xa_cavity(:,:) = 0._dp
        if ( pkernel%method .hasattr. "rigid" ) then
          ! First calculate the electrostatic forces (only)
          call PS_set_options(pkernel, cavitation_terms = .false.)
          ! ps_soft_PCM_forces adds the electrostatic forces
          ! We are calculating these outside
          call ps_soft_PCM_forces(pkernel, xa_cavity)
          ! Reverse forces to subtract the electrostatic forces
          ! in subsequent call (see psolver/src/environment[rigid_cavity_forces]
          do ia = 1 , na_u
            xa_cavity(:,ia) = - xa_cavity(:,ia)
          end do
          ! Now ask to also calculate the cavity forces
          ! The internal routine also calculates the electrostatic
          ! force, but we are now removing the electrostatic force
          ! and retain only cavity forces.
          call PS_set_options(pkernel, cavitation_terms = .true.)
          ! Calculate forces due to cavities
          call ps_soft_PCM_forces(pkernel, xa_cavity)
        end if

        ! Hartree => Rydberg
        do ia = 1, na_u
          fal(:,ia) = fal(:,ia) + 2._dp * xa_cavity(:,ia)
        end do
      end if

      ! Clean up memory
      call de_alloc(xa_cavity, 'xa_cavity', 'poisson_psolver')
    end if

    ! Shift data back
    call mesh_cshift(p_ntml, Vaux, lat_offset)

    ! In siesta poison.F the energies and stress tensors are calculated
    ! distributed. I.e. Siesta expects the Hartree energy and the
    ! stress tensor to be distributed.
    ! However, PSolver returns summed values.
    ! The factor two comes from Hartree => Rydberg

    Vaux(:) = 2._dp * Vaux(:)
    eh = 2._dp * energies%hartree / Nodes

    if ( shape == "bulk" ) then
      if ( eh < 0._dp ) then
        if ( Node == 0 ) then
          write(6,*) "WARNING: Psolver: U(rho-rhoatm) < 0 (changing to absolute value)"
        end if
        eh = abs(eh)
      end if
    end if

    ! The stress-tensor is definitely wrong. However, it seems that the Siesta
    ! convention is different than PSolver
    if ( calc_stress ) then
      ! Correct stress taking into account an extra 2*U/V term in the diagonal
      ! (See comment in routine poison in Siesta)
      ! factor of 2 here is intrinsic, not units!

!!$      ! This does NOT work.
!!$      if ( shape == "bulk" ) then
!!$        stress_correction = 2.0_dp * eh / volcel(cell)
!!$      else if ( shape == "slab" ) then
!!$        stress_correction = 2.0_dp * eh / area(cell(:, 1), cell(:,3))
!!$      else if ( shape == "wire" ) then
!!$        stress_correction = 2.0_dp * eh / VNORM(cell(:,3))
!!$      else
!!$        stress_correction = 0._dp
!!$      end if
      stress_correction = 2.0_dp * eh / volcel(cell)

      ! Here factor 2 is units
      stressl(1,1) = 2._dp * energies%strten(1) / Nodes + stress_correction
      stressl(2,2) = 2._dp * energies%strten(2) / Nodes + stress_correction
      stressl(3,3) = 2._dp * energies%strten(3) / Nodes + stress_correction
      stressl(3,2) = 2._dp * energies%strten(4) / Nodes
      stressl(2,3) = 2._dp * energies%strten(4) / Nodes
      stressl(3,1) = 2._dp * energies%strten(5) / Nodes
      stressl(1,3) = 2._dp * energies%strten(5) / Nodes
      stressl(1,2) = 2._dp * energies%strten(6) / Nodes
      stressl(2,1) = 2._dp * energies%strten(6) / Nodes
    end if

    if ( Nodes > 1 ) then

      ! Redistribute data back
      call setMeshDistr(UNIFORM, p_nml, p_nnml, p_ntml, p_ntpl)
      call distMeshData(UNIFORM_PSOLVER, Vaux, UNIFORM, V, KEEP)

      call de_alloc(Vaux, 'Vaux', 'poisson_psolver')

    end if

    ! Ending calculation
    call pkernel_free(pkernel)
    call f_lib_finalize_noreport()
    
    call timer('poisson_psolver', 2)

  contains

    function area(v, u) result(A)
      real(dp), intent(in) :: v(3), u(3)
      real(dp) :: A
      real(dp) :: vu_x(3)

      call cross(v, u, vu_x)
      A = VNORM(vu_x)
    end function area

  end subroutine poisson_psolver


  !< Grid mesh center-of-mass shifting
  !!
  !! Shift entire mesh depending on the center-of positions for each lattice
  !! vector.
  !! This should enable users to not rely on PSolver's requirements by
  !! shifting the data.
  !! Needless to say that performance will be greater if one ensures
  !! PSolvers requirements are fulfilled.
  !!
  !! Note: One cannot create the offsets in the initMeshDistr
  !! since that routine does not consider the mesh as periodic.
  !! If moremeshsubs would be altered for that it would make things
  !! much easier since it could *all* be handled in the distMeshData
  !! routine.
  subroutine mesh_cshift(ntml, V, shift)

    use parallel, only: Nodes
    use alloc, only: re_alloc, de_alloc
#ifdef MPI
    use mpi_siesta
#endif

    integer, intent(in) :: ntml(3)
    real(grid_p), intent(inout) :: V(ntml(1),ntml(2),ntml(3))
    integer, intent(in) :: shift(3)

    integer :: i1, i2, i3, j1, j2
    real(grid_p), pointer :: Vaux(:,:) => null()

    ! Luckily PSolver requires only distribution in C-axis.
    ! So A and B axis are easy to shift

    if ( shift(1) /= 0 .and. shift(2) /= 0 ) then
      call re_alloc(Vaux, 1, ntml(1), 1, ntml(2), 'Vaux', 'poisson_psolver_shift')

      ! Loop outer regions
      do i3 = 1, ntml(3)

        if ( shift(2) < 0 ) then
          j2 = ntml(2) + shift(2) + 1
        else
          j2 = shift(2)
        end if
        do i2 = 1, ntml(2)

          if ( shift(1) < 0 ) then
            j1 = ntml(1) + shift(1) + 1
          else
            j1 = shift(1)
          end if
          do i1 = 1 , ntml(1)
            Vaux(j1,j2) = V(i1,i2,i3)
            j1 = j1 + 1
            if ( j1 > ntml(1) ) j1 = 1
          end do

          j2 = j2 + 1
          if ( j2 > ntml(2) ) j2 = 1
        end do

        !call dcopy(ntml(1)*ntml(2), Vaux(1,1), 1, V(1,1,i3), 1)
        do i2 = 1, ntml(2)
          do i1 = 1, ntml(1)
            V(i1,i2,i3) = Vaux(i1,i2)
          end do
        end do

      end do

    else if ( shift(1) /= 0 ) then
      call re_alloc(Vaux, 1, ntml(1), 1, 1, 'Vaux', 'poisson_psolver_shift')

      ! Loop outer regions
      do i3 = 1, ntml(3)
        do i2 = 1, ntml(2)

          if ( shift(1) < 0 ) then
            j1 = ntml(1) + shift(1) + 1
          else
            j1 = shift(1)
          end if
          do i1 = 1 , ntml(1)
            Vaux(j1,1) = V(i1,i2,i3)
            j1 = j1 + 1
            if ( j1 > ntml(1) ) j1 = 1
          end do

          !call dcopy(ntml(1), Vaux(1,1), 1, V(1,i2,i3), 1)
          do i1 = 1, ntml(1)
            V(i1,i2,i3) = Vaux(i1,1)
          end do

        end do
      end do

    else if ( shift(2) /= 0 ) then
      call re_alloc(Vaux, 1, ntml(1), 1, ntml(2), 'Vaux', 'poisson_psolver_shift')

      ! Loop outer regions
      do i3 = 1, ntml(3)

        if ( shift(2) < 0 ) then
          j2 = ntml(2) + shift(2) + 1
        else
          j2 = shift(2)
        end if
        do i2 = 1, ntml(2)
          do i1 = 1, ntml(1)
            Vaux(i1,j2) = V(i1,i2,i3)
          end do
          j2 = j2 + 1
          if ( j2 > ntml(2) ) j2 = 1
        end do

        !call dcopy(ntml(1)*ntml(2), Vaux(1,1), 1, V(1,1,i3), 1)
        do i2 = 1, ntml(2)
          do i1 = 1, ntml(1)
            V(i1,i2,i3) = Vaux(i1,i2)
          end do
        end do

      end do

    end if

    call de_alloc(Vaux, 'Vaux', 'poisson_psolver_shift')

    ! Quick return if we don't need to shift last lattice vector
    if ( shift(3) == 0 ) return

    write(*,'(a,i4)') 'Number of grid-points required for shift: ', shift(3)
    call die('psolver: Please manually shift geometry such that &
        &the center of position along Z direction is at the center of &
        &the cell. &
        &AtomicCoordinatesOrigin COP!')

#ifdef MPI
    if ( Nodes > 1 ) then

    else
#endif

#ifdef MPI
    end if
#endif

  end subroutine mesh_cshift


  !< Setup PSolver variables to the dictionary
  !!
  !! This uses the PSolver dictionary to setup variables and passes
  !! from Siesta-options into the PSolver arguments.
  subroutine set_variables(dom, cell, calc_stress, dict)

    use parallel, only: Nodes, IONode
    use dictionaries, dict_set => set
    use at_domain
    use futile
    use siesta_geom, only : shape

    type(domain), intent(inout) :: dom
    real(dp), intent(in) :: cell(3,3)
    logical, intent(in) :: calc_stress
    type(dictionary), pointer :: dict ! Input parameters

    ! Initialize domain
    if ( psolver_force_pbc ) then
      dom = domain_new(ATOMIC_UNITS, [PERIODIC_BC, PERIODIC_BC, PERIODIC_BC], abc=cell)
    else
      select case ( shape )
      case ( 'molecule' )
        dom = domain_new(ATOMIC_UNITS, [FREE_BC, FREE_BC, FREE_BC], abc=cell)
      case ( 'slab' )
        dom = domain_new(ATOMIC_UNITS, [PERIODIC_BC, FREE_BC, PERIODIC_BC], abc=cell)
      case ( 'bulk' )
        dom = domain_new(ATOMIC_UNITS, [PERIODIC_BC, PERIODIC_BC, PERIODIC_BC], abc=cell)
      case ( 'wire' )
        dom = domain_new(ATOMIC_UNITS, [FREE_BC, FREE_BC, PERIODIC_BC], abc=cell)
      case default
        ! This should be enough to break PSolver ;)
        dom = domain_null()
      end select
    end if

    ! Data distribution will be always global:
    call set( dict//'setup'//'global_data', Nodes == 1) ! Hardwired, it cannot be otherwise.

    select case ( trim(psolver_accel) )
    case ( 'CUDA', 'cuda' )
      call set( dict//'setup'//'accel', "CUDA")
    case ( 'NONE', 'none' )
      ! Ok, just no accelerator
    case default
      call die("psolver_m: Unknown accelerator!")
    end select

    ! Set kernel and cavity variables:
    call set( dict//'kernel'//'isf_order', psolver_isf_order)
    call set( dict//'kernel'//'stress_tensor', calc_stress)
    call set( dict//'environment'//'fd_order', psolver_fd_order)
    call set( dict//'setup'//'verbose', psolver_verbose)
    if ( psolver_cavity ) then
      call set( dict//'environment'//'delta', psolver_delta)
      call set( dict//'environment'//'fact_rigid', psolver_fact_rigid)
      call set( dict//'environment'//'cavity', psolver_cavity_type)
      call set( dict//'environment'//'radii_set', psolver_radii_type)
      call set( dict//'environment'//'gps_algorithm', psolver_gps_algorithm)
      call set( dict//'environment'//'epsilon', psolver_epsilon)
      call set( dict//'environment'//'gammaS', psolver_gammaS)
      call set( dict//'environment'//'alphaS', psolver_alphaS)
      call set( dict//'environment'//'betaV', psolver_betaV)
      call set( dict//'environment'//'atomic_radii', psolver_atomic_radii)
    end if

  end subroutine set_variables

#else
contains

  subroutine psolver_dummy()
  end subroutine psolver_dummy

#endif
end module psolver_m
