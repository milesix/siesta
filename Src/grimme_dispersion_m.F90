!
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
module grimme_dispersion_m
  !! This module acts as an interface for Grimme's DFT-D3 library, which
  !! itself comprises two main sections: the D3 model itself (dftd3 module)
  !! and the MCTC library, which handles geometry and coordinates setup.
  !!
  !! We only need to call initialization once, and then energies and
  !! forces are calculated after SCF convergence.
  !!
  !! Even though the interfaced library holds different datatypes
  !! for parameters, we also keep our own dftd3_data_t so as to ease
  !! I/O from the fdf and avoid some problems in case the interface is
  !! changed in the future. This structure is initialized during module
  !! initialization in dftd3_initialize.
  !!
  !! Initialization can be called at any point before DFT-D3 forces
  !! contribution calculation, but it only needs to be called once. It
  !! also requires XC data to be already initialized in the code.
  !!
  !! dftd3_energy_forces, the subroutine which calculates forces, energy,
  !! and stress contributions, can also be called at any point after atom
  !! data is set up, since it only requires coordinates and atom identities.
  !! This is a post-SCF correction, so it does not enter the SCF loop.
  !!
  !! dftd3_get_periodic is an internal-only subroutine that attempts to guess
  !! system periodicity shape in order to send this information to the D3
  !! library. It relies on the shaper subroutine.
  !!
  !!
  !! THEORY BACKGROUND
  !! -----------------
  !!
  !! For ease of reading, we will present here three of the equations
  !! present in Grimme's work (DOI: 10.1063/1.3382344). Equation numbers
  !! correspond to those in the paper.
  !!
  !! (2) E(D3) = E(2-body) + E(3-body)
  !!
  !! (3) E(2-body) = sum_{A,B} ( s6 * C6AB / (rAB)^6 ) * f6(rAB)
  !!               + sum_{A,B} ( s8 * C8AB / (rAB)^8 ) * f8(rAB)
  !!
  !! Where A,B are atoms, sn and CnAB are parameters, rAB the distance
  !! between atoms, and fn a damping function:
  !!
  !! (4) fn(rAB) = 1 / ( 1 + 6 * ( rAB / (Srn * R0AB))^(-alp_n) )
  !!
  !! Where Srn, R0AB and alp_n are more parameters. Only CnAB and R0AB
  !! depend on the atoms, while Srn, Sn, and alp_n depend only on the DFT
  !! functional chosen.
  !!
  !! The 3-body term is a bit more complicated but at the same time
  !! requires less parametrization:
  !!
  !! (14) E(3-body) = Sum_{A,B,C} f3(rABC) * E(ABC)
  !!
  !! (11) E(ABC) = C9ABC * (3* cos(Ta) * cos(Tb) * cos(Tc) + 1) /
  !!                       (rAB * rBC * rCA)^3
  !!
  !! (13) C9ABC = - sqrt( C6AB * C6BC * C6AC )
  !!
  !! Where the value of C9ABC is, in truth, an approximation. Ta, Tb
  !! and Tc are the angles of the triangle formed by A, B and C, and
  !! the damping function f3 uses a value of 16 for alp and 4/3 for Sr.
#ifdef SIESTA__DFTD3
  use dftd3    , only : d3_param, rational_damping_param, zero_damping_param

  use precision, only : dp
  implicit none

  public :: dftd3_initialize
  public :: dftd3_energy_forces

  private

  type dftd3_data_t
    !! Datatype to contain all of the SIESTA
    !! input data for D3. s8 and rs8 are kept with those names
    !! (instead of s8 and rs8) for consistency with the D3 library.
    real(dp) :: s6       = 1.0_dp
      !! S6 coefficient that pre-multiplies all
      !! of the two-body C6 terms. (eq 3)
    real(dp) :: rs6      = 1.0_dp
      !! S_{r6} prefactor in the damping factor quotient. (eq 4)
    real(dp) :: s8       = 1.0_dp
      !! S8 coefficient that pre-multiplies all
      !! of the two-body C8 terms. (eq 3)
    real(dp) :: rs8      = 1.0_dp
      !! S_{r8} prefactor in the damping factor quotient.
      !! Usually set to 1.0. (eq 4)
    real(dp) :: s9       = 1.0_dp
      !! Weight for 3-body interaction terms. When zero, no 3-body
      !! terms are calculated.
    real(dp) :: a1       = 0.4_dp
      !! First parameter in Becke-Johnson damping
    real(dp) :: a2       = 5.0_dp
      !! First parameter in Becke-Johnson damping
    real(dp) :: alp      = 14.0_dp
      !! Exponent for the C6 damping factor. The C8 factor is alp + 2.0.
    logical  :: BJ_damp  = .true.
      !! Whether we use Becke-Johnson damping or zero damping variants.
    real(dp) :: cutoff_cn = 10.0_dp
      !! Cut-off for coordination.
    real(dp) :: cutoff_2b = 60.0_dp
      !! Cut-off for 2-body interactions.
    real(dp) :: cutoff_3b = 40.0_dp
      !! Cut-off for 3-body interactions.
  end type dftd3_data_t

  type(d3_param)               :: grimme_d3_dat
    !! A datatype within the external library that contains all relevant D3
    !! parameters.
  type(rational_damping_param) :: grimme_d3_dparam
    !! Contains damping information in the case of rational BJ damping.
  type(zero_damping_param)     :: grimme_d3_zparam
    !! Contains damping information in the case of zero damping.
  type(dftd3_data_t)           :: siesta_d3_dat
    !! All of the SIESTA input data that is passed to the D3 library.

contains

  subroutine dftd3_initialize( )
    !! Reads input data and initializes the D3 library.
    use fdf      , only : fdf_get
    use sys      , only : die

    use dftd3    , only : get_rational_damping, get_zero_damping, &
                          new_rational_damping, new_zero_damping
    use mctc_env , only : error_type
    use gridXC   , only : getXC=>gridxc_getXC
    use sys      , only : message

    implicit none
    integer           :: n_xc
    type(error_type) , allocatable :: d3_error
    character(len=20), allocatable :: fun_type(:), fun_author(:)

    siesta_d3_dat%BJ_damp = fdf_get( 'DFTD3.BJdamping', siesta_d3_dat%BJ_damp )

    call getXC( n = n_xc )
    if ( n_xc > 1 ) &
      call die( 'DFT-D3 functional initialization is not supported'&
              &' for mixed/cocktail functionals.'  )

    allocate( fun_type(n_xc), fun_author(n_xc) )
    call getXC( n = n_xc, func = fun_type, auth = fun_author )

    if ( fdf_get("DFTD3.UseXCDefaults", .true.) ) then
      if ( fun_type(1) == 'LDA' ) &
        call die( 'DFT-D3 functional initialization is not supported'&
                &' for LDA.' )
      if ( fun_type(1) == 'VDW' ) &
        call die( 'DFT-D3 functional initialization is not supported'&
                &' for VDW functionals.' )

      select case ( fun_author(1) )
      case ( 'PBE', 'pbe', 'PBESOL', 'pbesol', 'PBEsol', 'REVPBE', &
             'revpbe', 'revPBE', 'BLYP', 'LYP', 'blyp', 'lyp',    &
             'rpbe', 'RPBE', 'hse6', 'HSE6', 'pbe0', 'PBE0' )
        if ( (fun_author(1) == 'hse6') .or. (fun_author(1) == 'HSE6') ) &
          fun_author(1) = 'hse06'
          call message( 'INFO',  'DFT-D3: loading default parameters for' //&
                        ' functional ' // trim(fun_author(1)) // '.' )
      case default
        if ( allocated(fun_type)   ) deallocate( fun_type   )
        if ( allocated(fun_author) ) deallocate( fun_author )
        call die( 'This functional is not available with DFT-D3 '&
                &'functional initialization. Set DFTD3.UseXCDefaults '&
                &'to false and add custom parameters for s6, rs6 and s8.')
      end select

      siesta_d3_dat%s9 = fdf_get( 'DFTD3.s9', siesta_d3_dat%s9 )
      if ( siesta_d3_dat%BJ_damp ) then
        call get_rational_damping( grimme_d3_dat, fun_author(1), d3_error, &
                                   s9 = siesta_d3_dat%s9 )
      else
        call get_zero_damping( grimme_d3_dat, fun_author(1), d3_error, &
                               s9 = siesta_d3_dat%s9 )
      endif
      if ( allocated( d3_error ) ) then
         call message('WARNING',"d3_error in dftd3_initialize")
         return
      endif
    endif


    if ( fun_type(1) == 'VDW' ) &
      call message( 'WARNING', 'D3 corrections should ideally be used in '//&
                    'combination with vdW functionals.' )

    if ( allocated(fun_type)   ) deallocate( fun_type   )
    if ( allocated(fun_author) ) deallocate( fun_author )

    ! We overwrite parameters with custom values if present in the fdf.
    siesta_d3_dat%s6  = fdf_get( 'DFTD3.s6'   , siesta_d3_dat%s6  )
    siesta_d3_dat%rs6 = fdf_get( 'DFTD3.rs6'  , siesta_d3_dat%rs6 )
    siesta_d3_dat%s8  = fdf_get( 'DFTD3.s8'   , siesta_d3_dat%s8  )
    siesta_d3_dat%rs8 = fdf_get( 'DFTD3.rs8'  , siesta_d3_dat%rs8 )
    siesta_d3_dat%a1  = fdf_get( 'DFTD3.a1'   , siesta_d3_dat%a1  )
    siesta_d3_dat%a2  = fdf_get( 'DFTD3.a2'   , siesta_d3_dat%a2  )
    siesta_d3_dat%alp = fdf_get( 'DFTD3.alpha', siesta_d3_dat%alp )
    siesta_d3_dat%s9  = fdf_get( 'DFTD3.s9'   , siesta_d3_dat%s9  )

    grimme_d3_dat%s6  = siesta_d3_dat%s6
    grimme_d3_dat%rs6 = siesta_d3_dat%rs6
    grimme_d3_dat%s8  = siesta_d3_dat%s8
    grimme_d3_dat%s9  = siesta_d3_dat%s9
    grimme_d3_dat%rs8 = siesta_d3_dat%rs8
    grimme_d3_dat%a1  = siesta_d3_dat%a1
    grimme_d3_dat%a2  = siesta_d3_dat%a2
    grimme_d3_dat%alp = siesta_d3_dat%alp
    if ( siesta_d3_dat%BJ_damp ) then
      call new_rational_damping( grimme_d3_dparam, grimme_d3_dat )
    else
      call new_zero_damping( grimme_d3_zparam, grimme_d3_dat )
    endif

    ! These should be in atomic units.
    siesta_d3_dat%cutoff_cn = fdf_get( 'DFTD3.CoordinationCutoff', &
                                       siesta_d3_dat%cutoff_cn , 'Bohr' )
    siesta_d3_dat%cutoff_2b = fdf_get( 'DFTD3.2BodyCutoff', &
                                       siesta_d3_dat%cutoff_2b , 'Bohr' )
    siesta_d3_dat%cutoff_3b = fdf_get( 'DFTD3.3BodyCutoff', &
                                       siesta_d3_dat%cutoff_3b , 'Bohr' )

  end subroutine dftd3_initialize

  subroutine dftd3_energy_forces( natoms, atm_crd, iza, isa, lattice_vecs, &
                                  Edisp, Frc, stress )
    !! Calculates D3 corrections to energies, forces and stress.
    use alloc         , only : re_alloc, de_alloc
    use dftd3         , only : new_d3_model, get_dispersion, d3_model, &
                               realspace_cutoff
    use mctc_io       , only : structure_type, new_structure
    use periodic_table, only : symbol

    implicit none
    integer , intent(in)    :: natoms
      !! Number of atoms in unit cell.
    real(dp), intent(in)    :: atm_crd(3, natoms)
      !! Atom coordinates.
    integer , intent(in)    :: iza(natoms)
      !! Atomic numbers.
    integer , intent(in)    :: isa(natoms)
      !! Atomic species index.
    real(dp), intent(in)    :: lattice_vecs(3,3)
      !! Cell vectors.
    real(dp), intent(out)   :: Edisp
      !! Dispersion correction to energies.
    real(dp), intent(inout) :: Frc(3, natoms)
      !! Dispersion correction to forces.
    real(dp), intent(inout) :: stress(3,3)
      !! Dispersion correction to stress.

    external                      :: timer
    real(dp), external            :: volcel

    integer                       :: ii, jj
    logical                       :: isperiodic(3)
    real(dp)                      :: cell_vol
    real(dp), pointer, contiguous :: Gdisp(:,:), Gstress(:,:)
    character(len=4), pointer     :: atsyms(:)

    type(realspace_cutoff)        :: cuts
      ! Stores cut-off data.
    type(structure_type)          :: molec
      ! Stores molecule information and PBC.
    type(d3_model)                :: grimme_d3_model
      ! Stores all coefficients and pre-calculations needed for the D3 model.

    call timer( 'dftd3', 1 )

    nullify( Gdisp, Gstress )
    call re_alloc( Gdisp  , 1, 3, 1, natoms, 'Gdisp' , &
                   'dftd3_energy_forces' )
    call re_alloc( Gstress, 1, 3, 1,      3, 'Gstress', &
                   'dftd3_energy_forces' )

    nullify( atsyms )
    call re_alloc( atsyms, 1, natoms, 'atsyms' , 'dftd3_energy_forces' )

    do jj = 1, natoms
      atsyms(jj) = adjustl( symbol( iza(jj) ) )
    enddo

    call dftd3_get_periodic( natoms, isa, atm_crd, lattice_vecs, isperiodic )
    call new_structure( molec, iza, atsyms, atm_crd, lattice = lattice_vecs, &
                        periodic = isperiodic )

    call new_d3_model( grimme_d3_model, molec )

    cuts%cn    = siesta_d3_dat%cutoff_cn
    cuts%disp2 = siesta_d3_dat%cutoff_2b
    cuts%disp3 = siesta_d3_dat%cutoff_3b

    if ( siesta_d3_dat%BJ_damp ) then
      call get_dispersion( molec, grimme_d3_model, grimme_d3_dparam,        &
                           cuts, Edisp, gradient = Gdisp, sigma = Gstress )
    else
      call get_dispersion( molec, grimme_d3_model, grimme_d3_zparam,        &
                           cuts, Edisp, gradient = Gdisp, sigma = Gstress )
    endif

    call de_alloc( atsyms, 'atsyms' , 'dftd3_energy_forces' )

    ! The *2 factor is due to the conversion from Hartree to Rydberg
    Edisp = Edisp * 2.0_dp

    cell_vol = volcel( lattice_vecs )
    do jj = 1, 3
    do ii = 1, 3
      stress(ii,jj) = stress(ii,jj) - Gstress(ii,jj) * 2.0_dp / cell_vol
    enddo
    enddo

    do jj = 1, natoms
    do ii = 1, 3
      Frc(ii,jj) = Frc(ii,jj) - Gdisp(ii,jj) * 2.0_dp
    enddo
    enddo

    call de_alloc( Gdisp  , 'Gdisp'  , 'dftd3_energy_forces' )
    call de_alloc( Gstress, 'Gstress', 'dftd3_energy_forces' )

    call timer( 'dftd3', 2 )
  end subroutine dftd3_energy_forces

  subroutine dftd3_get_periodic( natoms, isa, at_crd, lat_vec, is_periodic )
    !! Attempts to guess system periodicity using the shaper routine.
    !! Possible results are molecule (i.e. non-periodic), chain (1D),
    !! slab (2D) and bulk (3D). This can also be manually set up via the
    !! DFTD3.Periodic input variable; in that case, periodicity guess is
    !! completely skipped.
    use fdf      , only : fdf_list, fdf_islist
    use precision, only : dp
    use sys      , only : message, die

    implicit none
    integer , intent(in)  :: natoms
      !! Total number of atoms.
    integer , intent(in)  :: isa(*)
      !! Atomic species index.
    real(dp), intent(in)  :: at_crd(3,natoms)
      !! Atomic coordinates.
    real(dp), intent(in)  :: lat_vec(3,3)
      !! Cell lattice vectors.
    logical , intent(out) :: is_periodic(3)

    external :: shaper
    character(len=8) :: cell_shape
    integer          :: nvec, ivec, jvec, per_list(3)
    real(dp)         :: per_vec(3,3), tmp_cellv(3,3), norm, vdiff, currv(3)

    if ( fdf_islist('DFTD3.Periodic') ) then
      nvec = -1
      call fdf_list( 'DFTD3.Periodic', nvec, per_list )

      if ( nvec > 3 ) then
        call message( 'WARNING', 'DFTD3.Periodic has more than 3 elements.'//&
                      ' Only the first 3 elements will be read.')
        nvec = 3
      endif

      if ( nvec > 0 ) then
        is_periodic(1:3) = .false.
        call fdf_list( 'DFTD3.Periodic', nvec, per_list )

        do ivec = 1, nvec
          if ( (per_list(ivec) > 3) .or. (per_list(ivec) < 0) ) &
            call die( 'DFTD3.Periodic indices must be 1, 2 or 3.' )

          is_periodic( per_list(ivec) ) = .true.
        enddo

      else
        call die( 'DFTD3.Periodic declared but no elements found in list.' )

      endif

    else
      nvec = 0
      per_vec(1:3,1:3) = 0.0_dp
      call shaper( lat_vec, natoms, isa, at_crd, cell_shape, nvec, per_vec )

      if ( nvec < 3 ) then
        call message( 'WARNING', 'System symmetry is of type '//&
                      trim(cell_shape)//&
                      ' but no DFTD3.Periodic list was found.'//&
                      ' Attempting to guess symmetry for D3.')

        is_periodic(1:3) = .false.
        if ( nvec == 0 ) then
          call message( 'INFO', 'System symmetry is of type '//&
                                trim(cell_shape)//'.')
          call message( 'INFO', 'Turning off periodicity for D3.')
          return
        endif

        ! Here the real guesswork starts: We'll try to see if any of the
        ! vectors provided by shaper() coincides with the cell vector,
        ! at least in direction.

        ! Normalize cell vectors.
        do ivec = 1, 3
          norm = lat_vec(1,ivec) ** 2 + lat_vec(2,ivec) ** 2 + &
                 lat_vec(3,ivec) ** 2
          norm = sqrt( norm )

          tmp_cellv(:,ivec) = lat_vec(:,ivec) / norm
        enddo

        do ivec = 1, nvec
          norm = per_vec(1,ivec) ** 2 + per_vec(2,ivec) ** 2 + &
                 per_vec(3,ivec) ** 2
          norm = sqrt( norm )
          per_vec(:,ivec) = per_vec(:,ivec) / norm

          ! We see the difference between the (normalized) periodicity
          ! vector and
          do jvec = 1, 3
            currv(:) = tmp_cellv(:,jvec)

            ! Here we check the sign of the scalar product. If
            ! currv == per_vec, then the scalar product is always positive.
            ! It might happen that currv == (- per_vec), in which case it
            ! is always negative, so we need to correct the sign.
            ! Lastly, if currv /= per_vec or -per_vec, this check changes
            ! nothing.
            norm = per_vec(1,ivec) * currv(1) + per_vec(2,ivec) * currv(2) +&
                   per_vec(3,ivec) * currv(3)

            if ( abs(norm) < 1e-6_dp ) cycle ! They are orthogonal...
            if ( norm < 0.0_dp ) currv(:) = - currv(:)

            vdiff = ( ( per_vec(1,ivec) - currv(1) ) ** 2 +&
                      ( per_vec(2,ivec) - currv(2) ) ** 2 +&
                      ( per_vec(3,ivec) - currv(3) ) ** 2 )
            vdiff = sqrt( vdiff )

            if ( vdiff < 1e-6_dp ) then
              ! Vectors are the same.
              is_periodic(jvec) = .true.
            endif

          enddo
        enddo

        ! We verify if we guessed the correct amount of
        ! periodicity vectors.
        if ( count( is_periodic(:) .eqv. .true. ) /= nvec ) then
          call message( 'INFO', 'Could not guess the right amount of'//&
                        ' periodicity directions for D3.')
          call message( 'INFO', &
                        'Enabling periodicity on all directions.')
          is_periodic(1:3) = .true.
        else
          call message( 'INFO', 'Periodicity direction guess succeeded.')
        endif

      else ! Gremlins were here
        is_periodic(1:3) = .true.

      endif
    endif

  end subroutine dftd3_get_periodic
#endif
end module grimme_dispersion_m
