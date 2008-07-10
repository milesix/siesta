! 
! This file is part of the SIESTA package.
!
! Copyright (c) Fundacion General Universidad Autonoma de Madrid:
! E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal
! and J.M.Soler, 1996-2008.
! 
! Use of this software constitutes agreement with the full conditions
! given in the SIESTA license, as signed by all legitimate users.
!
!******************************************************************************
! MODULE m_filter
! Handles k-filtering of strictly confined radial functions
! Written by J.M.Soler. April 2008
!******************************************************************************
!
!   PUBLIC procedures available from this module:
! function   kcPhi  : Returns planewave cutoff of a radial function
! subroutine filter : Filters out high-k components of a radial function
!
!   PUBLIC parameters, types, and variables available from this module:
! none
!
!******************************************************************************
!
!   USED module procedures:
! use sys,       only: die    ! Termination routine
! use m_radfft,  only: radfft ! Radial Fourier transform
!
!   USED module parameters:
! use precision, only: dp     ! Real double precision type
!
!   USED MPI procedures, interfaces and types:
! none
!
!   EXTERNAL procedures used:
! function   besshp : Returns spherical Bessel functions
! subroutine rdiag  : Diagonalization of real symmetric matrices
!
!******************************************************************************
!
! subroutine filter( l, nr, r, f, kc, norm_opt )
!
! Filters out high-k components of a radial function
!
! Input:
!  integer  l     : Angular momentum of the function to filter
!  integer  nr    : Number of radial points
!  real(dp) r(nr) : Radial mesh, in ascending order (in bohr)
!  real(dp) kc    : Cutoff in reciprocal space (in bohr^-1)
!  integer  norm_opt : Option to renormalize (to initial norm) after filtering:
!                        0 => do not renormalize
!                        1 => renormalize f
!                        2 => renormalize f**2
! Input/output:
!  real(DP) f(nr) : Radial function to be filtered
!
!******************************************************************************
!
! function kcPhi( l, nr, r, phi, etol, kexp )
!
! Returns planewave cutoff of a radial function
!
! Input:
!  integer  l        : Angular momentum of the radial function
!  integer  nr       : Number of radial points
!  real(dp) r(nr)    : Radial mesh, in ascending order (bohr)
!  real(dp) phi(nr)  : Radial function, assumed zero beyond r(nr)
!  real(dp) etol     : Kinetic energy allowed above kc (Ry)
!
! Output:
!  real(dp) kcPhi : Planewave cutoff kc (in bohr^-1) such that
!                     etol = Integral_kc^infinity(dk*k^4*phi(k)^2) 
!                          / Integral_0^infinity(dk*k^2*phi(k)^2)
!
!******************************************************************************

MODULE m_filter

! Used module procedures
  use sys,       only: die    ! Termination routine
  use m_radfft,  only: radfft ! Radial Fourier transform

! Used module parameters
  use precision, only: dp     ! Real double precision type
  use fdf, only: fdf_physical

! All public procedures (there are no public types, parameters, or variables):
PUBLIC:: &
  kcPhi, &! Returns planewave cutoff of a radial function
  filter  ! Filters out high-k components of a radial function

PRIVATE ! Nothing is declared public beyond this point

! Internal parameters for filter subroutine
  integer, parameter:: nmesh = 128        ! Number of radial integr. points
  integer, parameter:: minj  = 10         ! Min. num. of Bessel functions
  real(dp),parameter:: njkr  = 0.65_dp    ! Num. Bess. funcs. / (kc*rc)
  real(dp),parameter:: emax  = 1.e-2_dp   ! Min. accepted eigval
  real(dp),parameter:: krtol = 1.e-12_dp  ! Tol. for roots of Bess. funcs.
!  real(dp),parameter:: k2mix = 1.e-6_dp  ! Mix weight of k2 in 'Hamiltonian'
  real(dp),parameter:: k2mix = 0.10_dp  ! Mix weight of k2 in 'Hamiltonian'

! Interfaces to external routines
  interface
    DOUBLE PRECISION FUNCTION BESSPH (L,X)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    END
    subroutine rdiag(H,S,n,nm,nml,w,Z,neigvec,iscf,ierror)
      use precision, only: dp
      implicit none
      integer :: ierror, iscf, n, neigvec, nm, nml
      real(dp):: H(nml,nm), S(nml,nm), w(nml), Z(nml,nm)
    end subroutine rdiag
  end interface

  real(dp), save :: kc_needed = 0.0_dp

CONTAINS

!------------------------------------------------------------------------------

  subroutine filter(l,nr,r,func,factor,norm_opt)
    ! Filters out high-k components of a radial function
    ! Interface of J.M.Soler filter subroutine. E. Anglada 2007-2008
    !---------------------------------------------------------------------------
    ! Input:
    !  integer l     : Angular momentum of the function to filter
    !  integer nr    : Number of radial points
    !  real*8  r(nr) : Radial mesh, in ascending order
    !  real*8  factor: Scaling factor of the function (func*factor)
    !                  It's 1 except for orbitals as we're interested 
    !                  in their square.  
    !  integer norm_opt : Option to renormalize after filtering:
    !                       0 => do not renormalize
    !                       1 => renormalize f
    !                       2 => renormalize f**2
    ! Input/output:
    !  real*8  func(nr) : Radial function to be filtered
    !-------------------------------------------------------------------------
    ! Logic:
    !    The user imposes a kc (cutoff in reciprocal space) or we call kcPhi 
    !in order to determine the kc from the value of the kinetic energy above kc
    !-------------------------------------------------------------------------

    integer,  intent(in)                    :: l 
    integer,  intent(in)                    :: nr
    real(dp), dimension(1:nr), intent(in)   :: r
    real(dp), dimension(1:nr), intent(inout):: func
    real(dp), intent(in)                    :: factor 
    integer,  intent(in)                    :: norm_opt

    !Internal vars.
    real(dp) :: kmax, meshcutoff,etol,kmax_tol
    Integer  :: ir
    real(dp), parameter :: EkinTolDefault = 0.003_dp

    kmax = fdf_physical("FilterCutoff",0.0_dp,"Bohr^-1")
    
     !Siesta stores the orbs divided by r**l
    if (l /= 0) then
       do ir=2,nr
          func(ir)=func(ir)*r(ir)**(l)
       enddo
       func(1)=func(2)
    endif

    if (kmax .eq. 0.0)then
       etol = fdf_physical('Filter.Tol',EkinTolDefault,'Ry')
       !etol = fdf_double('Fil.Tol',0.003)
       !print *, "Tol=",etol
       kmax_tol = kcPhi(l,nr,r,func,etol)
       print *,"kmax for tol=",kmax_tol, etol      
       !if kmax_tol < 7.07 Bohr^-1 (50Ry) we fix it.
       !!if(kmax_tol < 7.07 ) then
       !!   kmax = factor*7.07
       !!else
       !!   kmax = factor *kmax_tol
       !!endif
       kmax = kmax_tol
    else
       kmax = factor*kmax
    endif

    write(6,'(a,f8.2,a)') "Equivalent plane wave cutoff:",kmax,' Bohr^-1'
    write(6,'(a,f8.2,a)') "Equivalent plane wave cutoff:",kmax**2,' Ry'

    call filter_kernel(l,nr,r(1:nr),func(1:nr),kmax,norm_opt)
 
    kmax_tol = (kcPhi(l,nr,r,func,etol))

    write(6,'(a,f8.2,a)') "Filter: The appropiate mesh cutoff is (at least):", &
         kmax_tol**2*factor,' Ry'

    !Siesta spects the orbs divided by r**l
    if (l /= 0) then
       do ir=2,nr
          func(ir)=func(ir)/(r(ir)**(l))
       enddo
       func(1)=func(2)
    endif

  end subroutine filter

!------------------------------------------------------------------------------

function kcPhi( l, nr, r, phi, etol )

! Returns planewave cutoff of a radial function

  implicit none

! Arguments
  integer,  intent(in) :: l        ! Angular momentum of function to filter
  integer,  intent(in) :: nr       ! Number of radial points
  real(dp), intent(in) :: r(nr)    ! Radial mesh, in ascending order (bohr)
  real(dp), intent(in) :: phi(nr)  ! Radial function to be filtered
  real(dp), intent(in) :: etol     ! Kinetic energy allowed above kc (Ry)
  real(dp)             :: kcPhi    ! Plane wavevector cutoff (bohr^-1)

! Internal parameters
  integer, parameter:: nmesh    = 1024    ! Number of radial integr. points
  integer, parameter:: kexp     = 2       ! k exponent in integral
  real(dp),parameter:: rmaxByrc = 3._dp   ! rmax/rc

! Internal variables and arrays
  integer :: ik, ik1, ik2, ir, nrc
  real(dp):: de, dk, dr, e, f(0:nmesh), kmax, kmesh(0:nmesh), k1, k2, &
             norm, pi, rc, rmax, rmesh(0:nmesh)

! Trap bad arguments  
  if (etol<=0._dp .or. etol>=1._dp) call die('kcPhi: ERROR: bad tol')
  if (kexp/=0 .and. kexp/=2)        call die('kcPhi: ERROR: bad kexp')

! Find cutoffs for integration mesh
  pi = acos(-1._dp)
  rc = r(nr)
  rmax = rc * rmaxByrc
  dr = rmax / nmesh
  nrc = nint( rc / dr )
  dr = rc / nrc   ! Make rc one of the mesh points
  rmax = nmesh * dr
  kmax = pi / dr
  dk = kmax / nmesh

! Find integration meshes
  forall(ir=0:nmesh) rmesh(ir) = ir*dr
  forall(ik=0:nmesh) kmesh(ik) = ik*dk

! Find function at integration points
  f(0:nrc) = interpolation( phi, r, rmesh(0:nrc) )
  f(nrc+1:nmesh) = 0

! Find Fourier transform of function
  call radfft( l, nmesh, rmax, f(0:nmesh), f(0:nmesh) )

! Normalize function
  norm = sum( kmesh(:)**2 * f(:)**2 ) * dk
  if (norm==0._dp) call die('kcPhi: ERROR: phi=0')
  f(:) = f(:) / sqrt(norm)

! Find cutoff kc such that tol = Integral_kc^infinity dk * k^(2*kexp) * f(k)^2
  e = 0
  do ik2 = nmesh,1,-1
    ik1 = ik2 - 1
    k1 = kmesh(ik1)
    k2 = kmesh(ik2)
    dk = k2 - k1
    de = (k1**(2+kexp) * f(ik1)**2 + k2**(2+kexp) * f(ik2)**2) * dk/2
    if (e+de > etol) then
      kcPhi = k2 - dk * (etol-e)/de
      exit
    end if
    e = e + de
  end do

end function kcPhi

!------------------------------------------------------------------------------

subroutine filter_kernel( l, nr, r, f, kc, norm_opt )

! Filters out high-k components of a radial function
! Ref: J.M.Soler notes of 17/10/2006.
! Written by J.M.Soler. Apr.2007

  implicit none

! Arguments
  integer,  intent(in)    :: l        ! Angular momentum of function to filter
  integer,  intent(in)    :: nr       ! Number of radial points
  real(dp), intent(in)    :: kc       ! Cutoff in reciprocal space
  real(dp), intent(in)    :: r(nr)    ! Radial mesh, in ascending order
  real(dp), intent(inout) :: f(nr)    ! Radial function to be filtered
  integer,  intent(in)    :: norm_opt ! Option to renormalize after filtering

! Internal variables and arrays
  integer :: i, ierror, ij, ij1, ij2, ik, ir, m, n, nj
  real(dp):: cij, dk, dkr, dr, f0norm, fg, fnorm, jl0, jl1, jl2, &
             k, k4j1(0:nmesh), kmesh(0:nmesh), kr, kr0, kr1, kr2, &
             pi, rc, rmesh(0:nmesh)
!             k, k4j1(0:nmesh), kmesh(0:2*nmesh), kr, kr0, kr1, kr2, &
!             pi, rc, rmesh(0:2*nmesh)
  integer, allocatable:: indx(:)
  real(dp),allocatable:: aux(:), c(:,:), e(:), fm(:), g(:,:), gr(:,:), &
                         h(:,:), jl(:,:), jlk(:,:), jlnorm(:), jlr(:), &
                         kj(:), s(:,:)

! Fix number of basis Bessel functions by empirical rule
  pi = acos(-1._dp)
  rc = r(nr)
  n = minj + nint(njkr*kc*rc)

! Find integration meshes
  dr = rc / nmesh
  dk = kc / nmesh
  do i = 0,nmesh
!  do i = 0,2*nmesh
    rmesh(i) = i*dr
    kmesh(i) = i*dk
  end do

! Allocate arrays
  allocate( aux(n), c(n,n), e(n), fm(0:nmesh), &
            g(0:nmesh,n), gr(nr,n), h(n,n), indx(n), &
            jl(0:nmesh,n), jlk(0:nmesh,n), jlnorm(n), jlr(n), kj(n), s(n,n) )
!            jl(0:nmesh,n), jlk(0:2*nmesh,n), jlnorm(n), jlr(n), kj(n), s(n,n) )

! Find roots of spherical Bessel functions 
  dkr = pi/2         ! Safe enough not to miss any root
  kr1 = dkr/2
  jl1 = bessph(l,kr1)
  nj = 0
  sampling: do
    kr2 = kr1 + dkr
    jl2 = bessph(l,kr2)
    if (jl1*jl2 < 0._dp) then  ! Root bracketted
      bisection: do
        kr0 = (kr1+kr2)/2
        if (kr2-kr1 < krtol) exit bisection
        jl0 = bessph(l,kr0)
        if (jl0*jl1 > 0._dp) then
          kr1 = kr0
          jl1 = jl0
        else
          kr2 = kr0
          jl2 = jl0
        end if
      end do bisection
      nj = nj + 1
      kj(nj) = kr0 / rc
      if (nj==n) exit sampling
    end if
    kr1 = kr2
    jl1 = jl2
  end do sampling

! Find normalized spher. Bessel funcs. at integration mesh points
  jl(:,:) = 0
  do ij = 1,n
    jlnorm(ij) = 0
    do ir = 0,nmesh
      kr = kj(ij)*rmesh(ir)
      jl(ir,ij) = bessph( l, kr )
      jlnorm(ij) = jlnorm(ij) + rmesh(ir)**2 * jl(ir,ij)**2 * dr
    end do ! ir
    jlnorm(ij) = sqrt(jlnorm(ij))
    jl(0:nmesh,ij) = jl(0:nmesh,ij) / jlnorm(ij)
  end do ! ij

! Find Fourier transform of cutted Bessel functions
! Ref: J.M.Soler notes of 16/04/2008
  do ij = 1,n
    cij = (-1._dp)**ij * 2*sqrt(rc/pi) * kj(ij)
    do ik = 0,nmesh
!    do ik = 0,2*nmesh
      k = kmesh(ik)
      if (abs(k-kj(ij)) > 100*krtol) then
        jlk(ik,ij) = cij * bessph(l,k*rc) / (k**2 - kj(ij)**2)
      else ! denominator too small
        ! Use j'(l,x)=-j(l+1,x)+l*jl(l,x)/x  and  j(l,kj(ij)*rc)=0
        k = (k + kj(ij)) / 2
        jlk(ik,ij) = cij * rc/(2*k) * &
                   ( -bessph(l+1,k*rc) + l*bessph(l,k*rc)/(k*rc) )
      end if
    end do
  end do ! ij

! Find 'Hamiltonian' (integral between kc and infinity of k**2*jl1*jl2)
  h(:,:) = 0
  do ij1 = 1,n
    h(ij1,ij1) = kj(ij1)**2 ! Integral from zero to infinity of diagonal terms
    k4j1(0:nmesh) = kmesh(0:nmesh)**4 * jlk(0:nmesh,ij1)
    do ij2 = 1,ij1
      ! Subtract integral up to kc using trapezoidal rule
!      h(ij1,ij2) = h(ij1,ij2) - sum( k4j1(:)*jlk(:,ij2) )*dk + &
!                   k4j1(nmesh)*jlk(nmesh,ij2)*dk/2
      ! Subtract integral up to kc using Simpson's rule
      h(ij1,ij2) = h(ij1,ij2) - &
        ( 4*sum(k4j1(1:nmesh-1:2)*jlk(1:nmesh-1:2,ij2)) + &
          2*sum(k4j1(2:nmesh-2:2)*jlk(2:nmesh-2:2,ij2)) + &
                k4j1(nmesh)      *jlk(nmesh,ij2)          ) * dk/3
      ! Subtract integral up to kc using higher-order routine
!      h(ij1,ij2) = h(ij1,ij2) - &
!        integral( nmesh+1, k4j1(:)*jlk(:,ij2), x=kmesh(:) )
      ! Symmetrize
      h(ij2,ij1) = h(ij1,ij2)
      write(33,*) ij1,ij2,h(ij1,ij2)
    end do ! ij2
  end do ! ij1

! Add some mixing weight of total kinetic energy to 'Hamiltonian'
  do ij = 1,n
    h(ij,ij) = (1-k2mix) * h(ij,ij) + k2mix * kj(ij)**2
    
  end do

! Initialize unity overlap matrix to call rdiag
  s = 0
  forall(i=1:n) s(i,i) = 1

! Diagonalize 'Hamiltonian'
  call filter_rdiag( h, s, n, n, n, e, c, n, 0, ierror )

! Subtract mixing weight of total kinetic energy from 'Hamiltonian'
!  do ij = 1,n
!    h(ij,ij) = ( h(ij,ij) - k2mix * kj(ij)**2 ) / (1-k2mix)
!  end do

! Print eigenvalues for debugging
!  print'(a,i4,e15.6)', ('filter: i,eigval/k2 =',i,e(i)/kj(i)**2-k2mix,i=1,n)

! Write leaked kinetic energy
!  open( unit=1, file='eigval.out' )
!  do i = 1,n
!    write(1,'(4f12.6)') kj(i)/kc, h(i,i)/kj(i)**2, &
!      sum(c(:,i)*matmul(h(:,:),c(:,i))) / kj(i)**2, e(i)/kj(i)**2
!  end do
!  close( unit=1 )

! Find how many eigenvalues are within filter tolerance
  m = 0
  do i = 1,n
    if (e(i)/kj(i)**2-k2mix > emax) exit
    m = i
  end do
  print*, 'filter: n, m =', n, m

! Find eigenvectors at integration mesh points
  g(0:nmesh,1:n) = matmul( jl(0:nmesh,1:n), c(1:n,1:n) )

! Write eigenvectors and Bessel functions
!  open( unit=1, file='bessel_k.out' )
!  do ik = 0,2*nmesh
!    write(1,'(20f12.6)') kmesh(ik)/kc, (kmesh(ik)**2*jlk(ik,i)**2,i=1,15)
!  end do
!  close( unit=1 )
!  open( unit=1, file='eigvec_k.out' )
!  do ik = 0,2*nmesh
!    write(1,'(20f12.6)') kmesh(ik)/kc, &
!                         (kmesh(ik)**2*sum(jlk(ik,:)*c(:,i))**2,i=1,15)
!  end do
!  close( unit=1 )

! Find function to be filtered at integration mesh points
  fm(0:nmesh) = interpolation( f, r, rmesh(0:nmesh) )

! Find eigenvectors at input mesh points
  do ir = 1,nr
    do ij = 1,n
      kr = kj(ij) * r(ir)
      jlr(ij) = bessph( l, kr ) / jlnorm(ij)
    end do
    gr(ir,:) = matmul( jlr(:), c(:,:) )
  end do

! Filter radial function (times r)
  f = 0
  do i = 1,n
    fg = sum( rmesh(0:nmesh)**2 * fm(:) * g(:,i) ) * dr
    if (i<=m) f(1:nr) = f(1:nr) + fg * gr(1:nr,i)
  end do

! Renormalize filtered function
  if (norm_opt==1) then
    f0norm = sum( matmul(rmesh(0:nmesh)**2*fm,g(:,1:n))*dr )
    fnorm  = sum( matmul(rmesh(0:nmesh)**2*fm,g(:,1:m))*dr )
    f = f * f0norm/fnorm
  else if (norm_opt==2) then
    f0norm = sum( (matmul(rmesh(0:nmesh)**2*fm,g(:,1:n))*dr)**2 )
    fnorm  = sum( (matmul(rmesh(0:nmesh)**2*fm,g(:,1:m))*dr)**2 )
    f = f * sqrt(f0norm/fnorm)
  else if (norm_opt/=0) then
    call die('filter: ERROR: invalid value of norm_opt')
  end if

! Deallocate arrays
  deallocate( aux, c, e, fm, g, gr, h, indx, jl, jlk, jlnorm, jlr, kj, s )

end subroutine filter_kernel

!------------------------------------------------------------------------------

  function interpolation( yold, xold, xnew ) result(ynew)

    ! Interpolates function yold(xold) at new points xnew
    ! Assumes that xold and xnew are in increasing order

    implicit none
    real(dp), intent(in) :: yold(:), xold(:), xnew(:)
    real(dp) :: ynew(size(xnew))

    integer :: i1, i2, inew, iold, nnew, nold

    nold = size(xold)
    nnew = size(xnew)

    iold = 1
    do inew = 1,nnew
      ! Find iold such that xold(iold) < xnew(inew) < xold(iold+1)
      do
        if (iold+1==nold .or. xold(iold+1) > xnew(inew)) exit
        iold = iold + 1
      end do
      ! Four-point Lagrange interpolation (3-point at extreme intervals)
      i1 = max(iold-1,1)
      i2 = min(iold+2,nold)
      ynew(inew) = lagrange( yold(i1:i2), xold(i1:i2), xnew(inew) )
    end do

  end function interpolation

!------------------------------------------------------------------------------

  function lagrange( yold, xold, xnew ) result(ynew)

    ! Lagrange interpolation of function yold(xold) at point xnew
    ! Based on routine polint of Numerical Recipes

    implicit none
    real(dp), intent(in) :: yold(:), xold(:), xnew
    real(dp) :: ynew

    integer :: i, i0, m, n
    real(dp):: c(size(xold)), d(size(xold)), den, dy, ho, hp, w

    n = size(xold)
    c = yold
    d = yold
    i0 = n/2
    ynew = yold(i0)
    i0 = i0 - 1
    do m = 1,n-1
      do i=1,n-m
        ho = xold(i) - xnew
        hp = xold(i+m) - xnew
        w = c(i+1) - d(i)
        den = ho - hp
        if (den==0._dp) call die('filter:lagrange: ERROR: den=0')
        d(i) = hp * w / den
        c(i) = ho * w / den
      end do
      if (2*i0 < n-m) then
        dy = c(i0+1)
      else
        dy = d(i0)
        i0 = i0 - 1
      endif
      ynew = ynew + dy
    end do

  end function lagrange


  ! ***************************************************************************
! subroutine rdiag(H,S,n,nm,nml,w,Z,neigvec,iscf,ierror)
!
! Simple replacement to subroutine rdiag of siesta
! J.M.Soler, April 2008
! ************************** INPUT ******************************************
! real*8 H(nml,nm)                 : Symmetric H matrix
! real*8 S(nml,nm)                 : Symmetric S matrix, ignored in this version
! integer n                        : Order of the generalized  system
! integer nm                       : Right hand dimension of H and S matrices
! integer nml                      : Left hand dimension of H and S matrices
!                                    which is greater than or equal to nm
! integer neigvec                  : No. of eigenvectors to calculate
! integer iscf                     : SCF cycle, ignored in this version
! ************************** OUTPUT *****************************************
! real*8 w(nml)                    : Eigenvalues
! real*8 Z(nml,nm)                 : Eigenvectors
! integer ierror                   : Flag indicating success code for routine
!                                  :  0 = success
!                                  : -1 = repeat call as memory is increased
!                                  :  1 = fatal error
! ***************************************************************************

subroutine filter_rdiag(H,S,n,nm,nml,w,Z,neigvec,iscf,ierror)

  use precision, only: dp
  implicit none

! Passed variables
  integer  :: ierror
  integer  :: iscf
  integer  :: n
  integer  :: neigvec
  integer  :: nm
  integer  :: nml
  real(dp) :: H(nml,nm)
  real(dp) :: S(nml,nm)
  real(dp) :: w(nml)
  real(dp) :: Z(nml,nm)

! Internal variables and arrays
  integer :: indx(n)
  real(dp):: aux(n), c(n,n), e(n)

! Diagonalize Hamiltonian
  c(1:n,1:n) = H(1:n,1:n)
  call tred2( c, n, n, e, aux )
  call tqli( e, aux, n, n, c )

! Order eigenvalues and eigenvectors by increasing eigval
  call ordix( e, 1, n, indx )
  call order( e, 1, n, indx )
  call order( c, n, n, indx )

! Copy eigenvectors and eigenvalues to output arrays
  w = 0
  Z = 0
  w(1:neigvec) = e(1:neigvec)
  Z(1:n,1:neigvec) = c(1:n,1:neigvec)
  ierror = 0

end subroutine filter_rdiag


END MODULE m_filter
