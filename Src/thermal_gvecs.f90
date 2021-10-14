module thermal_gvecs

  use precision,   only : dp, grid_p
  use thermal_flux_data

  implicit none

contains

  subroutine gvecs_init()
    !! Allocate g-vectors arrays with size based on
    !! stored number of global mesh divisions `mesh_vmd`.

    ngm_max_vmd = product(mesh_vmd)

    if(.not.(allocated(gg_vmd))) allocate(gg_vmd(ngm_max_vmd))
    if(.not.(allocated(g_vmd))) allocate(g_vmd(3, ngm_max_vmd))
    if(.not.(allocated(igsrt_vmd))) allocate(igsrt_vmd(ngm_max_vmd))
    if(.not.(allocated(igtongl_vmd))) allocate(igtongl_vmd(ngm_max_vmd))
  end subroutine gvecs_init


  subroutine gvecs_teardown()
    if(allocated(g_vmd)) deallocate(g_vmd)
    if(allocated(gg_vmd)) deallocate(gg_vmd)
    if(allocated(gl_vmd)) deallocate(gl_vmd)
    if(allocated(igsrt_vmd)) deallocate(igsrt_vmd)
    if(allocated(igplus_vmd)) deallocate(igplus_vmd)
    if(allocated(igtongl_vmd)) deallocate(igtongl_vmd)
  end subroutine gvecs_teardown


  pure function g_in_Gplus (i1, i2, i3)
    integer, intent(in) :: i1, i2, i3
    logical :: g_in_Gplus

    g_in_Gplus = ( i1 > 0 .or. &
         & ( i1 == 0 .and. i2 > 0 ) .or. &
         & ( i1 == 0 .and. i2 == 0 .and. i3 >= 0 ))
  end function g_in_Gplus


  subroutine gshells()
    !! Replicates `gshells` from `recvec_subs` module of Quantum Espresso.
    !! Calculate number of G shells: `ngl`, and the index ng = `igtongl`(ig)
    !! that gives the shell index `ng` for (local) G-vector of index `ig`.
    !!@note Non-variable cell @endnote
    !!@note Sorted with respect to `igplus_vmd` in G>. @endnote

    real(dp), parameter :: eps8 = 1.0E-8_DP  !! small constant
    integer :: ng, igl

    ngl_vmd = 1
    igtongl_vmd(1) = 1
    do ng = 2, ngm_plus_vmd ! ngm_max_vmd
       if (gg_vmd(igplus_vmd(ng)) > gg_vmd(igplus_vmd(ng - 1)) + eps8) then
          ngl_vmd = ngl_vmd + 1
       endif
       igtongl_vmd(ng) = ngl_vmd
    enddo

    if(.not.(allocated(gl_vmd))) allocate(gl_vmd(ngl_vmd))

    gl_vmd(1) = gg_vmd(1)
    igl = 1
    do ng = 2, ngm_plus_vmd ! ngm_max_vmd
       if (gg_vmd(igplus_vmd(ng)) > gg_vmd(igplus_vmd(ng - 1)) + eps8) then
          igl = igl + 1
          gl_vmd(igl) = gg_vmd(igplus_vmd(ng))
       endif
    enddo

    if (igl /= ngl_vmd) call die("gshells: igl <> ngl_vmd")

  end subroutine gshells


  subroutine hpsort_eps (n, ra, ind, eps)
    !---------------------------------------------------------------------
    ! sort an array ra(1:n) into ascending order using heapsort algorithm,
    ! and considering two elements being equal if their values differ
    ! for less than "eps".
    ! n is input, ra is replaced on output by its sorted rearrangement.
    ! create an index table (ind) by making an exchange in the index array
    ! whenever an exchange is made on the sorted data array (ra).
    ! in case of equal values in the data array (ra) the values in the
    ! index array (ind) are used to order the entries.
    ! if on input ind(1)  = 0 then indices are initialized in the routine,
    ! if on input ind(1) != 0 then indices are assumed to have been
    !                initialized before entering the routine and these
    !                indices are carried around during the sorting process
    !
    ! no work space needed !
    ! free us from machine-dependent sorting-routines !
    !
    ! adapted from Numerical Recipes pg. 329 (new edition)
    !
    ! use kinds, only : DP
    implicit none
    !-input/output variables
    integer, intent(in) :: n
    integer, intent(inout) :: ind (*)
    real(dp), intent(inout) :: ra (*)
    real(dp), intent(in) :: eps
    !-local variables
    integer :: i, ir, j, l, iind
    real(dp) :: rra
    ! initialize index array
    if (ind (1) .eq.0) then
       do i = 1, n
          ind (i) = i
       enddo
    endif
    ! nothing to order
    if (n.lt.2) return
    ! initialize indices for hiring and retirement-promotion phase
    l = n / 2 + 1

    ir = n

    sorting: do

       ! still in hiring phase
       if ( l .gt. 1 ) then
          l    = l - 1
          rra  = ra (l)
          iind = ind (l)
          ! in retirement-promotion phase.
       else
          ! clear a space at the end of the array
          rra  = ra (ir)
          !
          iind = ind (ir)
          ! retire the top of the heap into it
          ra (ir) = ra (1)
          !
          ind (ir) = ind (1)
          ! decrease the size of the corporation
          ir = ir - 1
          ! done with the last promotion
          if ( ir .eq. 1 ) then
             ! the least competent worker at all !
             ra (1)  = rra
             !
             ind (1) = iind
             exit sorting
          endif
       endif
       ! wheter in hiring or promotion phase, we
       i = l
       ! set up to place rra in its proper level
       j = l + l
       !
       do while ( j .le. ir )
          if ( j .lt. ir ) then
             ! compare to better underling
             if ( abs(ra(j)-ra(j+1)).ge.eps ) then
                if (ra(j).lt.ra(j+1)) j = j + 1
             else
                ! this means ra(j) == ra(j+1) within tolerance
                if (ind (j) .lt.ind (j + 1) ) j = j + 1
             endif
          endif
          ! demote rra
          if ( abs(rra - ra(j)).ge.eps ) then
             if (rra.lt.ra(j)) then
                ra (i) = ra (j)
                ind (i) = ind (j)
                i = j
                j = j + j
             else
                ! set j to terminate do-while loop
                j = ir + 1
             end if
          else
             !this means rra == ra(j) within tolerance
             ! demote rra
             if (iind.lt.ind (j) ) then
                ra (i) = ra (j)
                ind (i) = ind (j)
                i = j
                j = j + j
             else
                ! set j to terminate do-while loop
                j = ir + 1
             endif
          end if
       enddo
       ra (i) = rra
       ind (i) = iind

    end do sorting
    !
  end subroutine hpsort_eps


  subroutine setup_thermal_gvecs()
    !! Sets up entities for VMD-schematics.
    !! Called from `dhscf` on each substep (`Base` and `Virtual`):
    !! `siesta()` -> `siesta_forces()` -> `dhscf()` -> `setup_thermal_gvecs()`
    use fdf, only: fdf_physical
    use parallel,    only : Node, Nodes, ProcessorY
    use sys,         only : die
    use alloc,       only : re_alloc, de_alloc
    use m_fft,       only : fft     ! 3-D fast Fourier transform
    use cellsubs,    only : reclat  ! Finds reciprocal lattice vectors
    use cellsubs,    only : volcel  ! Finds unit cell volume
    use m_chkgmx,    only : chkgmx  ! Checks planewave cutoff of a mesh

    integer               :: N1, N2, N3, iplus
    integer               :: I, I1, I2, I3, IX, J, J1, J2, J3, JX
    integer               :: NP, NG, NG2, NG3
    integer               :: ProcessorZ, Py, Pz, J2min, J2max
    integer               :: J3min, J3max, J2L, J3L, NRemY, NRemZ
    integer               :: BlockSizeY, BlockSizeZ
    real(dp)              :: C, B(3,3), DU, G(3), G2, G2MAX
    real(dp)              :: PI, VG, VOLUME, PI8
    real(dp),   parameter :: K0(3)= (/ 0.0, 0.0, 0.0 /), TINY= 1.0e-15

    ! Temporary arrays for sorting of g-vectors
    real(dp), allocatable :: g2sort_g(:)
    integer,  allocatable :: gplus_flag(:)

    ! Start time counter
    call timer('SETUP_THERMAL_FLUX_GVECS', 1)

    alat_vmd = fdf_physical('LatticeConstant',0.0_dp,'Bohr')
    if (alat_vmd==0.0_dp) call die("VMD requires alat set")

    at_vmd(:,:) = cell_vmd(:,:)/alat_vmd

    call gvecs_init()           ! <- ngm_max_vmd init here also

    N1 = ntml_vmd(1)
    N2 = ntml_vmd(2)
    N3 = ntml_vmd(3)

    NP = product(ntml_vmd)
    NG = product(mesh_vmd)      ! <- =ngm_max_vmd

    ! Allocate memory for the Fourier image of charge:
    if (substep .eq. 1) then
       call re_alloc(charge_g,1,2,1,NP,'charge_g','setup_vmd')
    end if

    allocate(g2sort_g(NG))
    g2sort_g(:) = 1.0e20_dp
    allocate(gplus_flag(NG))
    gplus_flag(:) = 0

    ! Find unit cell volume
    volume = volcel(cell_vmd)

    ! Find reciprocal lattice vectors
    call reclat(cell_vmd, B, 1 )

    ! Find maximun planewave cutoff
    G2MAX = 1.0e30_dp
    call CHKGMX( K0, B, mesh_vmd, G2MAX )

    ! Copy density to complex array
!$OMP parallel do default(shared), private(I)
    do i = 1, NP
       ! charge_g(1,i) = thtr_Rho(i)
       charge_g(1,i) = Rho_save(i,substep)
       charge_g(2,i) = 0.0_grid_p
    enddo
!$OMP end parallel do

    ! Forward Fourier transform of density
    call fft(charge_g, mesh_vmd, -1)  ! -1 means Forward
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
    NG2 = mesh_vmd(2)
    NG3 = mesh_vmd(3)
    BlockSizeY = ((NG2/nsm_vmd)/ProcessorY)*nsm_vmd
    BlockSizeZ = ((NG3/nsm_vmd)/ProcessorZ)*nsm_vmd
    NRemY = (NG2 - BlockSizeY*ProcessorY)/nsm_vmd
    NRemZ = (NG3 - BlockSizeZ*ProcessorZ)/nsm_vmd
    J2min = (Py-1)*BlockSizeY + nsm_vmd*min(Py-1,NRemY)
    J2max = J2min + BlockSizeY - 1
    if (Py-1.lt.NRemY) J2max = J2max + nsm_vmd
    J2max = min(J2max,NG2-1)
    J3min = (Pz-1)*BlockSizeZ + nsm_vmd*min(Pz-1,NRemZ)
    J3max = J3min + BlockSizeZ - 1
    if (Pz-1.lt.NRemZ) J3max = J3max + nsm_vmd
    J3max = min(J3max,NG3-1)

    ! First loop over J-indeces: to get the Hartree potential.

    ngm_plus_vmd = 0 ! init the counter

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

             !NOTE: (/I1, I2, I3/) corresponds to (/i, j, k/)
             ! but the order is different from that in QE.
             ! Look at following invocation in `g_in_Gplus`.
             G(:)= B(:,1) * I1 + B(:,2) * I2 + B(:,3) * I3
             G2 = G(1)**2 + G(2)**2 + G(3)**2
             J = 1 + J1 + N1 * J2L + N1 * N2 * J3L

             g_vmd(:,J) = G(:)
             gg_vmd(J)  = G2

             if (G2.LT.G2MAX) then
                ! mill_unsorted(1:3, J) = (/I3, I2, I1/) !DEBUG: also strange order here
                if(G2.GT.TINY) then
                   g2sort_g(J) = G2
                else
                   g2sort_g(J) = 0.0_dp
                end if
                !NOTE: The order of index selection I1-I3-I2
                ! to correspond with the G> hemisphere in QE.
                if(g_in_Gplus(I1, I3, I2)) then
                   ngm_plus_vmd = ngm_plus_vmd + 1
                   gplus_flag(J) = 1
                end if
             end if
          enddo
       enddo
    enddo

    igsrt_vmd(1) = 0
    call hpsort_eps(ngm_max_vmd, g2sort_g, igsrt_vmd, TINY) !DEBUG

    ! Reduce indeces of vectors in G> to separate `igplus_vmd` array:
    if(.not.(allocated(igplus_vmd))) allocate(igplus_vmd(ngm_plus_vmd))

    iplus = 0
    do i=1,NG
       if(gplus_flag(igsrt_vmd(i)).eq.1) then
          iplus = iplus + 1
          igplus_vmd(iplus) = igsrt_vmd(i)
          ! write(995, *) "ig:", i, "ngtmp", ngtmp, "igsrt:", igsrt_vmd(i), "igplus:", igplus_vmd(ngtmp),&
          !      & "gp:", g_vmd(:,igplus_vmd(ngtmp)) !DEBUG
       end if
    end do
    if (iplus /= ngm_plus_vmd) call die("setup_thermal_gvecs: iplus <> ngm_plus_vmd")

    ! determine first nonzero g vector
    if (gg_vmd(1).le.TINY) then
       gstart_vmd = 2
    else
       gstart_vmd = 1
    endif

    ! determine g-vector shells:
    call gshells()              !NOTE: Sorted with respect to `igplus_vmd`.
    ! gl_debug_loop: do i=1,ngl_vmd
    !    write(994, *) "igl:", i, "gl_vmd(igl)", gl_vmd(i)
    ! end do gl_debug_loop
    ! igtongl_debug_loop: do i=1,ngm_plus_vmd
    !    write(993, *) "ng:", i, "igtongl_vmd(ng)", igtongl_vmd(i)
    ! end do igtongl_debug_loop

    deallocate(g2sort_g)
    deallocate(gplus_flag)

    ! Stop time counter
    call timer('SETUP_THERMAL_FLUX_GVECS', 2)

  end subroutine setup_thermal_gvecs

end module thermal_gvecs
