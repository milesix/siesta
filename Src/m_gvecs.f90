module gvecs
  !! This module aims to replicate some funcrionality
  !! of reciprocal lattice-related routines similar to
  !! what is realized for QE in `gvect`, `recvec` etc.

  use precision, only: dp, grid_p

  implicit none
  save

  integer :: ngm  = 0 !! local  number of G vectors (on this processor)
  integer :: ngl = 0  !! number of G-vector shells
  integer, allocatable :: nl(:), nlm(:)
  !! nl  = fft index for G-vectors (with gamma tricks, only for G>)
  !! nlm = as above, for G< (used only with gamma tricks)
  integer :: gstart = 2 !! index of the first G vector whose module is > 0

  real(dp), allocatable, target :: gg(:)
  !! G^2 in increasing order
  !! (for QE it's in units of tpiba2=(2pi/a)^2)

  real(dp), pointer :: gl(:)
  !! gl(i) = i-th shell of G^2 (in units of tpiba2)

  integer, allocatable, target :: igtongl(:)
  !! shell index for n-th G-vector

  real(dp), allocatable, target :: g(:,:)
  !! G-vectors cartesian components
  !! (for QE it's in units tpiba =(2pi/a)  )

  integer, allocatable, target :: mill(:,:)
  !! mill = miller index of G vectors

  integer, allocatable, target :: ig_l2g(:)

  real(dp) :: gcutm = 0.0_DP   !! ecutrho/(2 pi/a)^2, cut-off for |G|^2

contains

  subroutine gvecs_init(mesh)
    integer, intent(in) :: mesh(3)
    !! Number of global mesh divisions

    ngm = mesh(1)*mesh(2)*mesh(3)

    allocate( gg(ngm) )
    allocate( g(3, ngm) )
    allocate( mill(3, ngm) )
    allocate( nl (ngm) )
    allocate( nlm(ngm) )
    allocate( ig_l2g(ngm) )
    allocate( igtongl(ngm) )

  end subroutine gvecs_init


  subroutine gvecs_teardown()
    if ( allocated(nlm) ) deallocate(nlm)
    if ( allocated(nl) ) deallocate(nl)
    if ( allocated(mill) ) deallocate(mill)
    if ( allocated(g) ) deallocate(g)
    if ( allocated(gg) ) deallocate(gg)
    if ( allocated(ig_l2g) ) deallocate(ig_l2g)
    if ( allocated(igtongl) ) deallocate(igtongl)
  end subroutine gvecs_teardown


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
    real(DP), intent(inout) :: ra (*)
    real(DP), intent(in) :: eps
    !-local variables
    integer :: i, ir, j, l, iind
    real(DP) :: rra
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


  pure function g_in_Gplus (i1, i2, i3)
    integer, intent(in) :: i1, i2, i3
    logical :: g_in_Gplus

    g_in_Gplus = ( i1 > 0 .or. &
         & ( i1 == 0 .and. i2 > 0 ) .or. &
         & ( i1 == 0 .and. i2 == 0 .and. i3 > 0 ))
  end function g_in_Gplus


  subroutine gvecs_gen(CELL, Mesh)
    use precision,   only : dp, grid_p
    use parallel,    only : Node, Nodes, ProcessorY
    use sys,         only : die
    use alloc,       only : re_alloc, de_alloc
    use m_fft,       only : fft     ! 3-D fast Fourier transform
    use cellsubs,    only : reclat  ! Finds reciprocal lattice vectors
    use cellsubs,    only : volcel  ! Finds unit cell volume
    use m_chkgmx,    only : chkgmx  ! Checks planewave cutoff of a mesh

    use iso_c_binding, only: c_loc, c_f_pointer

    ! use m_virtual_step_data, only: BASE_STEP, VIRTUAL_STEP
    ! use m_virtual_step_data, only: substep, istep_vmd
    ! use siesta_options, only: virtual_dt, virtual_md_Jhart
    ! use hartree_flux_data, only: Vhart_deriv, h_flux_Jhart

!     Input/output variables
    ! integer, intent(in)    :: N1, N2, N3
    !! Number of mesh divisions in each cell vector
    integer, intent(in)    :: Mesh(3)
    !! Number of global mesh divisions

    ! real(grid_p), intent(in) :: RHO(N1*N2*N3)
    !! Densitiy at mesh points
    real(dp), intent(in)     :: CELL(3,3)
    !! Unit cell vectors

    ! integer, intent(in)  :: NSM
    !! Number of sub-mesh points per mesh point along each axis

!     Local variables
    integer               :: NP, NG, NG1, NG2, NG3
    real(dp)              :: B(3,3), G2, G2MAX
    real(dp)              :: PI, VG, VOLUME, PI8
    real(dp),   parameter :: K0(3)= (/ 0.0, 0.0, 0.0 /), TINY= 1.0e-15

    integer :: i, j, k, ngm_max, ngm_save, ig
    logical :: gamma_only = .true. ! since so far I'm interested in this case
    real(dp):: t(3), tt
    real(dp), allocatable :: g2sort_g(:)
    integer, allocatable :: igsrt(:)
    integer, allocatable :: mill_g(:,:), mill_unsorted(:,:)
    !! array containing all g vectors generators, on all processors: replicated data

    ! complex(grid_p), dimension (:), pointer :: tc_v => null()
    ! complex(grid_p), dimension (:), pointer :: tc_dot => null()
    ! real(dp), allocatable                   :: fac(:)


    ! Start time counter
    call timer('COMP_GVECS_GEN', 1)

    ! Init internal (and not only) variable values

    ! NP = N1 * N2 * N3
    NG = Mesh(1)*Mesh(2)*Mesh(3)

    PI  = 4.0_dp * atan(1.0_dp)
    PI8 = PI * 8._dp

    ngm_max = ngm
    ! ngm_save = ngm
    ngm = 0

    ! gg(:) = gcutm + 1.d0

    ! Allocate local memory

    ! allocate(fac(NP)); fac(:) = 0.0_dp ! allocate factors array
    ! for elements in the required hemisphere they will be re-assigned
    allocate( mill_g( 3, ngm_max ),mill_unsorted( 3, ngm_max ) )
    allocate( igsrt( ngm_max ) )
    allocate( g2sort_g( ngm_max ) )
    !
    g2sort_g(:) = 1.0d20


    ! Find unit cell volume
    volume = volcel(cell)

    ! Find reciprocal lattice vectors
    call reclat(cell, B, 1 )

    ! Find maximun planewave cutoff
    G2MAX = 1.0e30_dp
    call CHKGMX( K0, B, Mesh, G2MAX )

    ! NOTE: I'm not estimating block-grid-dimensions for now
    ! since this is for now a one-processor non-parallel version.
    ! An exaple how to do it should exist in similar routines like `rhofft`
    ! Instead I assign NG(1-3) from the Mesh(1-3):
    NG1 = (Mesh(1)-1)/2     ! ni = (dfftp%nr1-1)/2
    NG2 = (Mesh(2)-1)/2     ! nj = (dfftp%nr2-1)/2
    NG3 = (Mesh(3)-1)/2     ! nk = (dfftp%nr3-1)/2
    !
    ! write (666,*) ' ni,nj,nk ', NG1, NG2, NG3

    iloop: do i = -NG1, NG1
       !
       ! gamma-only: exclude space with x < 0
       !
       if ( gamma_only .and. i < 0) cycle iloop
       jloop: do j = -NG2, NG2
          !
          ! gamma-only: exclude plane with x = 0, y < 0
          !
          if ( gamma_only .and. i == 0 .and. j < 0) cycle jloop
          kloop: do k = -NG3, NG3
             !
             ! gamma-only: exclude line with x = 0, y = 0, z < 0
             !
             if ( gamma_only .and. i == 0 .and. j == 0 .and. k < 0) cycle kloop
             t(:) = i * B(:,1) + j * B(:,2) + k * B(:,3) ! G; rename it to e.g. GV
             tt = t(1)**2 + t(2)**2 + t(3)**2            ! G2
             if(tt.LE.G2MAX) then
                ngm = ngm + 1
                ! if (ngm > ngm_max) call die("Too many g-vectors")
                mill_unsorted( :, ngm ) = (/ i,j,k /)
                if ( tt > TINY ) then
                   g2sort_g(ngm) = tt
                else
                   g2sort_g(ngm) = 0.d0
                endif
             end if
          end do kloop
       end do jloop
    end do iloop

    ! write (667,*) ' ngm, ngm_max', ngm,ngm_max
    ! if (ngm  /= ngm_max) &      ! which is true in my case - why?
    !      call die('ggen', 'g-vectors missing !', abs(ngm - ngm_max))
    ! FIXME: hack for ngm, explicitly setting ngm_max limit:
    ngm_max = ngm
    ! !!!!!
    ! write (667,*) ' ngm, ngm_max', ngm,ngm_max

    igsrt(1) = 0
    call hpsort_eps( ngm, g2sort_g, igsrt, TINY )

    ! mill_g(1,:) = mill_unsorted(1,igsrt(:))
    ! mill_g(2,:) = mill_unsorted(2,igsrt(:))
    ! mill_g(3,:) = mill_unsorted(3,igsrt(:))
    mill_g(1,1:ngm_max) = mill_unsorted(1,igsrt(1:ngm_max))
    mill_g(2,1:ngm_max) = mill_unsorted(2,igsrt(1:ngm_max))
    mill_g(3,1:ngm_max) = mill_unsorted(3,igsrt(1:ngm_max))

    deallocate( g2sort_g, igsrt, mill_unsorted )

    ngm = 0
    ngloop: do ig = 1, ngm_max
       i = mill_g(1, ig)
       j = mill_g(2, ig)
       k = mill_g(3, ig)

       ngm = ngm + 1
       ig_l2g( ngm ) = ig

       g(1:3, ngm) = i * B(:,1) + j * B(:, 2) + k * B(:, 3)
       gg(ngm) = sum(g (1:3, ngm)**2)
       ! if (ngm > ngm_max) call die("Too many g-vectors")

       ! write(999,*) g(:, ngm), gg(ngm)

    end do ngloop
    ! write (668,*) ' ngm, ngm_save', ngm,ngm_save

    deallocate( mill_g )

    ! Stop time counter
    call timer('COMP_GVECS_GEN', 2)

  end subroutine gvecs_gen


  subroutine gshells(vc)
    !! Replicates `gshells` from `recvec_subs` module of Quantum Espresso.
    !! Calculate number of G shells: `ngl`, and the index ng = `igtongl`(ig)
    !! that gives the shell index `ng` for (local) G-vector of index `ig`.

    ! USE gvect,              ONLY : gg, ngm, gl, ngl, igtongl

    logical, intent(in) :: vc
    real(dp), parameter :: eps8 = 1.0E-8_DP  !! small constant

    integer :: ng, igl

    if ( vc ) then
       ! in case of a variable cell run each G vector has its shell
       ngl = ngm
       gl => gg
       do ng = 1, ngm
          igtongl(ng) = ng
       enddo
    else
       ! G vectors are grouped in shells with the same norm
       ngl = 1
       igtongl(1) = 1
       do ng = 2, ngm
          if (gg(ng) > gg(ng - 1) + eps8) then
             ngl = ngl + 1
          endif
          igtongl(ng) = ngl
       enddo

       allocate (gl(ngl))

       gl(1) = gg(1)
       igl = 1
       do ng = 2, ngm
          if (gg(ng) > gg(ng - 1) + eps8) then
             igl = igl + 1
             gl(igl) = gg(ng)
          endif
       enddo

       if (igl /= ngl) call die("gshells: igl <> ngl")

    endif

  end subroutine gshells

end module gvecs
