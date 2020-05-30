module hartree_flux_procs

  implicit none

  private
  public :: compute_Jhart

contains

  pure function g_in_Gplus (i1, i2, i3)
    integer, intent(in) :: i1, i2, i3
    logical :: g_in_Gplus

    g_in_Gplus = ( i1 > 0 .or. &
         & ( i1 == 0 .and. i2 > 0 ) .or. &
         & ( i1 == 0 .and. i2 == 0 .and. i3 > 0 ))
  end function g_in_Gplus


  subroutine compute_Jhart(CELL, N1, N2, N3, Mesh, RHO, NSM)
    use precision,   only : dp, grid_p
    use parallel,    only : Node, Nodes, ProcessorY
    use sys,         only : die
    use alloc,       only : re_alloc, de_alloc
    use m_fft,       only : fft     ! 3-D fast Fourier transform
    use cellsubs,    only : reclat  ! Finds reciprocal lattice vectors
    use cellsubs,    only : volcel  ! Finds unit cell volume
    use m_chkgmx,    only : chkgmx  ! Checks planewave cutoff of a mesh

    use iso_c_binding, only: c_loc, c_f_pointer

    use m_virtual_step_data, only: BASE_STEP, VIRTUAL_STEP
    use m_virtual_step_data, only: substep, istep_vmd
    use siesta_options, only: virtual_dt, virtual_md_Jhart
    use hartree_flux_data, only: Vhart_deriv, h_flux_Jhart

!     Input/output variables
    integer, intent(in)    :: N1, N2, N3
    !! Number of mesh divisions in each cell vector
    integer, intent(in)    :: Mesh(3)
    !! Number of global mesh divisions

    real(grid_p), intent(in) :: RHO(N1*N2*N3)
    !! Densitiy at mesh points
    real(dp), intent(in)     :: CELL(3,3)
    !! Unit cell vectors

    integer, intent(in)  :: NSM
    !! Number of sub-mesh points per mesh point along each axis

!     Local variables
    integer               :: I, I1, I2, I3, IX, J, J1, J2, J3, JX
    integer               :: NP, NG, NG2, NG3
    integer               :: ProcessorZ, Py, Pz, J2min, J2max
    integer               :: J3min, J3max, J2L, J3L, NRemY, NRemZ
    integer               :: BlockSizeY, BlockSizeZ
    real(dp)              :: C, B(3,3), DU, G(3), G2, G2MAX
    real(dp)              :: PI, VG, VOLUME, PI8
    real(dp),   parameter :: K0(3)= (/ 0.0, 0.0, 0.0 /), TINY= 1.0e-15

    real(grid_p), pointer :: charge_g(:,:)
    !! Charge density in the reciprocal space.
    !! Turns into Hartree potential during code execution.
    !! It is not returned - instead, stored in `m_hartree_flux_data`
    !! for computation of the derivative.
    complex(grid_p), dimension (:), pointer :: tc_v => null()
    complex(grid_p), dimension (:), pointer :: tc_dot => null()
    real(dp), allocatable                   :: fac(:)


    ! Start time counter
    call timer('COMP_JHART', 1)

    NP = N1 * N2 * N3
    NG = Mesh(1)*Mesh(2)*Mesh(3)

    ! Allocate local memory
    nullify( charge_g )
    call re_alloc( charge_g, 1, 2, 1, NP, 'charge_g', 'comp_Jh' )

    if (substep .eq. BASE_STEP) then
       call re_alloc(Vhart_deriv,1,2,1,NP,'Vhart_deriv','comp_Jh')
    end if

    allocate(fac(NP)); fac(:) = 0.0_dp ! allocate factors array
    ! for elements in the required hemisphere they will be re-assigned

    ! Find unit cell volume
    volume = volcel(cell)

    ! Find reciprocal lattice vectors
    call reclat(cell, B, 1 )

    ! Find maximun planewave cutoff
    G2MAX = 1.0e30_dp
    call CHKGMX( K0, B, Mesh, G2MAX )

    ! Copy density to complex array
!$OMP parallel do default(shared), private(I)
    do i = 1, NP
       charge_g(1,i) = rho(i)
       charge_g(2,i) = 0.0_grid_p
    enddo
!$OMP end parallel do

    ! Forward Fourier transform of density
    call fft(charge_g, Mesh, -1)  ! -1 means Forward
    charge_g(:,:) = charge_g(:,:) * volume / dble(NG)

    ! Work out processor grid dimensions
    ProcessorZ = Nodes/ProcessorY
    if (ProcessorY*ProcessorZ.ne.Nodes) then
       call die('ERROR: ProcessorY must be a factor of the &
            &number of processors!')
    end if
    Py = (Node/ProcessorZ) + 1
    Pz = Node - (Py - 1)*ProcessorZ + 1

    ! Multiply by 8*PI/G2 to get the potential
    PI  = 4.0_dp * atan(1.0_dp)
    PI8 = PI * 8._dp
    NG2 = Mesh(2)
    NG3 = Mesh(3)
    BlockSizeY = ((NG2/NSM)/ProcessorY)*NSM
    BlockSizeZ = ((NG3/NSM)/ProcessorZ)*NSM
    NRemY = (NG2 - BlockSizeY*ProcessorY)/NSM
    NRemZ = (NG3 - BlockSizeZ*ProcessorZ)/NSM
    J2min = (Py-1)*BlockSizeY + NSM*min(Py-1,NRemY)
    J2max = J2min + BlockSizeY - 1
    if (Py-1.lt.NRemY) J2max = J2max + NSM
    J2max = min(J2max,NG2-1)
    J3min = (Pz-1)*BlockSizeZ + NSM*min(Pz-1,NRemZ)
    J3max = J3min + BlockSizeZ - 1
    if (Pz-1.lt.NRemZ) J3max = J3max + NSM
    J3max = min(J3max,NG3-1)

    ! First loop over J-indeces: to get the Hartree potential.

    do J3 = J3min,J3max
       J3L = J3 - J3min
       if (J3.gt.NG3/2) then
          I3 = J3 - NG3
       else
          I3 = J3
       endif
       do J2 = J2min,J2max
          J2L = J2 - J2min
          if (J2.gt.NG2/2) then
             I2 = J2 - NG2
          else
             I2 = J2
          endif
          do J1 = 0,N1-1
             if (J1.gt.N1/2) then
                I1 = J1 - N1
             else
                I1 = J1
             endif

             G(1)= B(1,1) * I1 + B(1,2) * I2 + B(1,3) * I3
             G(2)= B(2,1) * I1 + B(2,2) * I2 + B(2,3) * I3
             G(3)= B(3,1) * I1 + B(3,2) * I2 + B(3,3) * I3

             G2 = G(1)**2 + G(2)**2 + G(3)**2
             J = 1 + J1 + N1 * J2L + N1 * N2 * J3L

             ! Hemisphere selection with indeces order priority:
             ! I1 - I3 - I2, to keep similarity with QE.
             if ((G2.LT.G2MAX .AND. G2.GT.TINY) &
                  &.and. g_in_Gplus(I1, I3, I2)) then
                fac(J) = PI8 / G2 / VOLUME
             end if

          enddo
       enddo
    enddo

    ! This means that V_hart(G)=0 if G is in the "negative" hemisphere
    ! I assume that this is done to do further operations in the same way
    ! as QE. In particular using a half-space avoids having to find -G
    ! as G's are generated.
    ! Just DO NOT use this as the "real" V_H(G) later on

    charge_g(1,:) = charge_g(1,:)*fac(:)
    charge_g(2,:) = charge_g(2,:)*fac(:)

    ! If in the `BASE` step: store Hartree potential in `Vhart_deriv`.
    ! If in the `VIRTUAL` step:
    ! 1. get the derivative of Hartree potential;
    ! 2. grab C-to-Fortran pointers to the potential and its derivative
    ! (as complex numbers);
    ! 3. second loop over J-indeces: to compute the Hartree flux component
    ! according to the formula 1.17 in Aris Marcolongo's manual.

    ! NOTE: The sign of the terms for the sum in [1.17] is opposite to ours:
    ! that is due to opposite phases of FFT used in SIESTA and QE.

    if (substep .eq. BASE_STEP) then
       Vhart_deriv(:,:) = charge_g(:,:)
    else
       Vhart_deriv(:,:) = (charge_g(:,:) - Vhart_deriv(:,:))/virtual_dt
       h_flux_Jhart(:) = 0.d0

       ! In A. Marcolongo's notes, they take the semi-sum of base+virtual
       ! for tc_v...
       call c_f_pointer(c_loc(charge_g), tc_v, [NP])
       call c_f_pointer(c_loc(Vhart_deriv), tc_dot, [NP])

       do J3 = J3min,J3max
          J3L = J3 - J3min
          if (J3.gt.NG3/2) then
             I3 = J3 - NG3
          else
             I3 = J3
          endif
          do J2 = J2min,J2max
             J2L = J2 - J2min
             if (J2.gt.NG2/2) then
                I2 = J2 - NG2
             else
                I2 = J2
             endif
             do J1 = 0,N1-1
                if (J1.gt.N1/2) then
                   I1 = J1 - N1
                else
                   I1 = J1
                endif

                G(1)= B(1,1) * I1 + B(1,2) * I2 + B(1,3) * I3
                G(2)= B(2,1) * I1 + B(2,2) * I2 + B(2,3) * I3
                G(3)= B(3,1) * I1 + B(3,2) * I2 + B(3,3) * I3

                G2 = G(1)**2 + G(2)**2 + G(3)**2
                J = 1 + J1 + N1 * J2L + N1 * N2 * J3L

                ! Hemisphere selection with indeces order priority:
                ! I1 - I3 - I2, to keep similarity with QE.
                if (G2.LT.G2MAX .AND. G2.GT.TINY) then
                   if (g_in_Gplus(I1, I3, I2)) then
                      h_flux_Jhart(1:3) = h_flux_Jhart(1:3) + &
                           & dimag(conjg(-tc_dot(J))*tc_v(J)) * &
                           & G(1:3) * volume / 4.d0 / PI
                   end if
                end if
             enddo
          enddo
       enddo

       print*, "[Jhart]: ", h_flux_Jhart

       nullify(tc_v)
       nullify(tc_dot)
    end if

    call de_alloc( charge_g, 'charge_g', 'gradVh' )

    ! Stop time counter
    call timer('COMP_JHART', 2)

  end subroutine compute_Jhart

end module hartree_flux_procs
