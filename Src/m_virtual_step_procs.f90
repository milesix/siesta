module m_virtual_step_procs
  !! Core virtual_step_procs are responsible for saving/restoring
  !! the state of the system between `BASE` and `VIRTUAL` steps
  !! in VMD scheme, and for calling subroutines that implement
  !! VMD logic for corresponding substep.
  !! For the setup phase, they also aim to replicate some functionality
  !! of reciprocal lattice-related routines similar to what is realized
  !! for QE in `gvect`, `recvec` etc.
  !! Since some of these entities are common for separate components
  !! calculated in VMD-scheme, but not used in general in Siesta,
  !! and since they are bound to calculation of `charge_g`,
  !! I place them here.

  use precision, only: dp, grid_p

  use m_virtual_step_data
  use siesta_options

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


  subroutine setup_step_vmd()
    !! Sets up entities for VMD-schematics.
    !! Called from `dhscf` on each substep (`Base` and `Virtual`):
    !! `siesta()` -> `siesta_forces()` -> `dhscf()` -> `setup_step_vmd()`
    use fdf, only: fdf_physical
    use precision,   only : dp, grid_p
    use parallel,    only : Node, Nodes, ProcessorY
    use sys,         only : die
    use alloc,       only : re_alloc, de_alloc
    use m_fft,       only : fft     ! 3-D fast Fourier transform
    use cellsubs,    only : reclat  ! Finds reciprocal lattice vectors
    use cellsubs,    only : volcel  ! Finds unit cell volume
    use m_chkgmx,    only : chkgmx  ! Checks planewave cutoff of a mesh

    use xc_flux_data, only: thtr_Rho

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
    call timer('SETUP_VMD', 1)

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
    if (substep .eq. BASE_STEP) then
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
       charge_g(1,i) = thtr_Rho(i)
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
    ! ALWAYS init the counters, dammit
    ! I SPENT AN HOUR ON EXPLODING STACK COS I DID NOT SET IT HERE TO ZERO

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
    if (iplus /= ngm_plus_vmd) call die("setup_step_vmd: iplus <> ngm_plus_vmd")

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
    call timer('SETUP_VMD', 2)

  end subroutine setup_step_vmd


  subroutine base_md_step_logic ()
    !! Wrapper over logic performed in the `Base` MD step,
    !! e.g. output integrated energy components etc.

    use sparse_matrices, only: Dscf

    use ks_flux_data,  only: Dscf_deriv
    use ks_flux_data,  only: init_ks_flux_data
    use xc_flux_data,  only: init_xc_flux_data
    use ion_flux_data, only: init_ion_flux_data
    use zero_flux_data,only: init_zero_flux_data
    use hartree_flux_procs, only: init_Jhart_BASE

    if ( virtual_md_Jks ) then
       call init_ks_flux_data()

       Dscf_deriv(:,:) = Dscf(:,:)
    end if

    if ( virtual_md_Jhart ) call init_Jhart_BASE()

    if ( virtual_md_Jxc )   call init_xc_flux_data()

    if ( virtual_md_Jion )  call init_ion_flux_data()

    if ( virtual_md_Jzero ) call init_zero_flux_data()

  end subroutine base_md_step_logic


  subroutine virtual_md_step_logic ()
    !! Wrapper over logic performed in the `Virtual` MD step,
    !! e.g. computation of derivatives of heat flux components etc.

    use sparse_matrices, only: Dscf, gradS

    use hartree_flux_data
    use hartree_flux_procs

    use ks_flux_data, only: Dscf_deriv
    use ks_flux_procs, only: compute_Jks

    use xc_flux_data, only: thtr_Rho, thtr_Rho_deriv, thtr_dexcdGD
    use xc_flux_data, only: xc_flux_Jxc, xc_flux_dvol

    use ion_flux_data, only: ion_flux_a, ion_flux_b
    use ion_flux_data, only: ion_flux_c, ion_flux_d, ion_flux_e
    use ion_flux_data, only: ion_flux_Jion
    use ion_flux_procs, only: compute_Jion_a, compute_Jion_b
    use ion_flux_procs, only: compute_Jion_cde

    use zero_flux_data,  only: zero_flux_Jzero
    use zero_flux_procs, only: compute_Jzero

    use ks_flux_data,      only: reset_ks_flux_data
    use xc_flux_data,      only: reset_xc_flux_data
    use ion_flux_data,     only: reset_ion_flux_data

    ! use ion_flux_data, only: ion_flux_mesh, ion_flux_cell

    integer :: i

    if ( virtual_md_Jhart ) then
       call compute_Jhart_VIRTUAL()
       print*, "[Jhart] ", h_flux_Jhart(:)
    end if

    if ( virtual_md_Jks ) then
       ! Dscf_deriv(:,:) = (Dscf_deriv(:,:)-Dscf(:,:))/virtual_dt
       Dscf_deriv(:,:) = (Dscf(:,:)-Dscf_deriv(:,:))/virtual_dt

       call compute_Jks()

       ! call reset_ks_flux_data()
    end if

    if ( virtual_md_Jxc ) then
       ! compute Rho derivative on a virtual step:
       thtr_Rho_deriv(:) = (thtr_Rho_deriv(:)-thtr_Rho(:))/virtual_dt

       ! xc_flux_Jxc(1:3) = 0.0_dp <- CHECK: should be reset to 0 already from `init_xc_flux_data`

       do i=1,size(thtr_Rho)
          xc_flux_Jxc(1:3) = xc_flux_Jxc(1:3) &
               & - thtr_Rho_deriv(i) * thtr_dexcdGD(i,1:3,1) &
               & * xc_flux_dvol
       end do

       print*, "[Jxc] ", xc_flux_Jxc(:)

    end if

    if ( virtual_md_Jion ) then

       call compute_Jion_a()

       print*, "[Jion] flux A: ", ion_flux_a(:)

       !NOTE:
       ! The following contribution behaves like a non-diffusive mass flux.
       ! According to Davide Tisi, It does not contribute to the resulting
       ! thermal conductivity and is therefore excluded from computation.
       !
       ! call compute_Jion_b()
       ! print*, "[Jion] flux B: ", ion_flux_b(:)

       call compute_Jion_cde()

       print*, "[Jion] flux C: ", ion_flux_c(:)
       print*, "[Jion] flux D: ", ion_flux_d(:)
       print*, "[Jion] flux E: ", ion_flux_e(:)

       ! call reset_ion_flux_data()
    end if

    if ( virtual_md_Jzero ) then
       call compute_Jzero()

       print*, "[Jzero] ", zero_flux_Jzero(:)
    end if

    !FIXME: regroup resetters in corresponding subroutine
    if ( virtual_md_Jion) call reset_ion_flux_data()
    if ( want_virtual_step_md ) call reset_xc_flux_data()
    if ( virtual_md_Jks ) call reset_ks_flux_data()

    ! if ( want_vmd_in_dhscf ) call reset_xc_flux_data()

  end subroutine virtual_md_step_logic


  subroutine reset_virtual_step_md ()
    !! Deallocates auxiliary arrays for before- and after-Base
    !! MD move system state.

    use class_Sparsity
    use alloc, only : de_alloc
    use ks_flux_data, only: base_step_sparse_pattern

    use ks_flux_data,      only: reset_ks_flux_data
    use xc_flux_data,      only: reset_xc_flux_data
    use ion_flux_data,     only: reset_ion_flux_data
    use hartree_flux_data, only: reset_hartree_flux_data

    if ( want_virtual_step_md ) then !TODO: make this check in `siesta_end`

       if (virtual_md_Jhart) call reset_hartree_flux_data()
       ! if ( virtual_md_Jion) call reset_ion_flux_data()
       ! if ( want_vmd_in_dhscf ) call reset_xc_flux_data()
       ! if ( virtual_md_Jks ) call reset_ks_flux_data()

       ! clean the aux `Base`-step sparse pattern
       call delete( base_step_sparse_pattern )

       call de_alloc(charge_g, "charge_g", "virtual_step_md")
       nullify(charge_g)

       call gvecs_teardown()    ! deallocate g-vectors data

       ! deallocate full density matrices
       call de_alloc(Dfull,  "Dfull",  "virtual_step_md")
       call de_alloc(Dderiv, "Dderiv", "virtual_step_md")

       nullify(Dfull)
       nullify(Dderiv)

       ! deallocate generalized coordinates and forces arrays
       call de_alloc(fa_before_move, "fa_before_move", "virtual_step_md")
       call de_alloc(xa_before_move, "xa_before_move", "virtual_step_md")
       call de_alloc(va_before_move, "va_before_move", "virtual_step_md")

       call de_alloc(fa_after_move, "fa_after_move", "virtual_step_md")
       call de_alloc(xa_after_move, "xa_after_move", "virtual_step_md")
       call de_alloc(va_after_move, "va_after_move", "virtual_step_md")

       nullify(fa_before_move)
       nullify(xa_before_move)
       nullify(va_before_move)

       nullify(fa_after_move)
       nullify(xa_after_move)
       nullify(va_after_move)
    end if

  end subroutine reset_virtual_step_md


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

end module m_virtual_step_procs
