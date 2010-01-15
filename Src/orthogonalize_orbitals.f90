module m_orthogonalize_orbitals

use precision, only: dp
use sys,       only: die
use atom_options, only: swap_zetas_only, just_refill_zetas
use atom_options, only: keep_zeta_ordering, keep_vna

implicit none

integer, parameter :: MAXNZ = 10     ! Maximum number of zetas per shell
integer, parameter :: NPTS  = 1000   ! Number of points for grid
real(dp)           :: xgrid(NPTS)    ! radial values in grid
real(dp)           :: delta          ! common grid spacing
real(dp)           :: max_cutoff     ! common maximum radius

public :: orthogonalize_orbitals
private

CONTAINS

subroutine orthogonalize_orbitals()
use atm_types, only: species_info, species, nspecies

type(species_info), pointer        :: spp
integer  :: i, is, n, zeta
logical  :: in_group 

integer  :: orb_list(MAXNZ)

do is = 1, nspecies
   spp => species(is)
   in_group = .false.
!   write(0,*) "Number of nl orbitals: ", spp%n_orbnl
   i = 0
   orbitals: do 
      i = i + 1
      if (i > spp%n_orbnl) exit orbitals
      zeta = spp%orbnl_z(i)
      write(0,*) "processing orbital : ", i, " with zeta: ", zeta
      if (zeta == 1) then
         if (in_group) then
            !
            ! Will deal with group, saving this orbital for
            ! the next round
            i = i - 1

            ! if our group was single-zeta, continue
            if (n == 1) then
               in_group = .false.
               cycle orbitals
            endif

            ! close_group and dispatch
            call orthogonalize_subspace(n,orb_list(1:n),spp)

            in_group = .false.

         else   ! not in_group

            !  create group
            n = 1
            ! add this orbital to group
            orb_list(n) = i
            in_group = .true.

         endif

      else    ! z /= 1
         ! add this orbital to group
         n = n + 1
         if (n > MAXNZ) call die("Increase MAXNZ")
         orb_list(n) = i

      endif

   enddo  orbitals
   !
   if (just_refill_zetas .or. keep_zeta_ordering) then
      ! do nothing
   else
      ! We have changed the first zeta, so we need to recompute
      ! the neutral-atom potential
      if (keep_vna) then
         ! do nothing, the user says
         else
            call recompute_vna(spp)
         endif
   endif
   !
enddo    ! species

end subroutine orthogonalize_orbitals
!----------------------------------------------------------

subroutine orthogonalize_subspace(n,orb_list,spp)
  use atm_types, only: species_info
  use radial,    only: rad_func, evaluate
  use m_recipes, only: sort

  integer, intent(in)   :: n
  integer, intent(in)   :: orb_list(n)
  type(species_info), pointer    :: spp

  integer  :: io, ii, i, j, iold, ip, l
  real(dp) :: norm
  type(rad_func)        :: tmp_func

  real(dp), allocatable   :: v(:,:)
  integer                 :: indx(MAXNZ)
  real(dp)                :: cutoff(MAXNZ)
  logical                 :: already_done(MAXNZ)

  allocate(v(NPTS,n))

  do io = 1, n
     ii = orb_list(io)
     cutoff(io)   = spp%orbnl(ii)%cutoff
  end do
!  write(0,*) "Processing orbitals: ", orb_list(1:n)
!  write(0,*) cutoff(1:n)


  ! sort, with orb with shortest rc first
  ! the sort routine is programmed for sorting
  ! in ascending order
  !
  if (just_refill_zetas .or. keep_zeta_ordering) then
     do i = 1, n
        indx(i) = i
     enddo
     write(0,*) "Not swapping zetas"
  else
     call sort(n,cutoff,indx)
  endif
  ! debug
  !
  ! Use a common grid for all orbitals
  !
  max_cutoff = maxval(cutoff(1:n))
  delta = max_cutoff / (NPTS-1)   ! from 0 to max_cutoff
  ! Auxiliary array for interpolation
  do i = 1, NPTS
     xgrid(i) = (i-1) * delta
  enddo
  !
  do io=1,n
     ii = orb_list(indx(io))
     l = spp%orbnl_l(ii)
     v(1:NPTS,io) = evaluate(spp%orbnl(ii),xgrid(1:NPTS))
     if ( l /= 0) v(:,io) = v(:,io) * xgrid(:)**l
  enddo
  !
  if (swap_zetas_only .or. just_refill_zetas) then
     !
     write(0,*) "Not actually orthogonalizing zetas"
  else
     call gram_schmidt(v)
  endif
  ! orthonormalize
  do io=1,n
     norm = sqrt(scalar_product(v(:,io),v(:,io)))
     v(:,io) = v(:,io) / norm
  enddo

  ! put back in tables
  do io=1,n
     ii = orb_list(indx(io))
     if (keep_zeta_ordering) then
        ! All the orbitals extend to the maximum cutoff
        call insert_in_table(NPTS,v(:,io),xgrid,delta,     &
             spp%orbnl(ii),spp%orbnl_l(ii),max_cutoff)
     else
        ! Each orbital keeps its old cutoff
        call insert_in_table(NPTS,v(:,io),xgrid,delta,     &
             spp%orbnl(ii),spp%orbnl_l(ii))
     endif
  enddo

  ! re-sort orbitals so that the shortest zeta is first
  already_done(1:n) = .false.
  do io=1,n
     if (already_done(io)) cycle
     ip = indx(io)
     ii = orb_list(ip)
     iold = orb_list(io)
     tmp_func = spp%orbnl(ii)
     spp%orbnl(ii) = spp%orbnl(iold)
     spp%orbnl(iold) = tmp_func
     spp%orbnl_z(iold) = io
     spp%orbnl_z(ii) = indx(io)
     already_done(ip) = .true.
  enddo
  !
  deallocate(v)

end subroutine orthogonalize_subspace
!--------------------------------------------------------

subroutine gram_schmidt(v)
! Performs an orthogonalization of the columns of v

real(dp), intent(inout) :: v(:,:)

integer :: i, j, k

k = size(v,dim=2)

do j = 1, k
   do i = 1, j-1
      v(:,j) = v(:,j) - proj(v(:,i),v(:,j))
   enddo
enddo
end subroutine gram_schmidt
!----------------------------------------------------------

function proj(v1,v2) result (projection)
! Projects v2 onto v1
!
! The underlying scalar product is a functional
! scalar product

real(dp), intent(in) :: v1(:), v2(:)
real(dp)             :: projection(size(v1))

if (size(v1) /= size(v2)) call die("not same length")
projection = v1 * scalar_product(v1,v2) / scalar_product(v1,v1)

end function proj
!---------------------------------------------------------
function scalar_product(v1,v2) result(scalar)
real(dp), intent(in) :: v1(:), v2(:)
real(dp)             :: scalar

integer :: i

scalar = 0.0_dp
do i = 1, NPTS
   scalar = scalar + v1(i) * v2(i) * xgrid(i)*xgrid(i)
enddo
scalar = scalar * delta

end function scalar_product
!---------------------------------------------------------

subroutine insert_in_table(np,phi,rgrid,dr,func,l,req_cutoff)

use radial,    only: rad_func, rad_setup_d2
use m_recipes, only: polint

integer, intent(in)           :: np
real(dp), intent(in)          :: phi(1:np)
real(dp), intent(in)          :: rgrid(1:np)
real(dp), intent(in)          :: dr
type(rad_func), intent(inout) :: func
integer, intent(in)           :: l
real(dp), intent(in), optional :: req_cutoff ! Force new cutoff

integer, parameter  :: npoint = 1    ! half-degree of interpolation

integer   :: ntbmax, nr, nn, nmin, nmax, i, j
real(dp)  :: r, dy

!
! Polynomial interpolation into points of radial function grid
!
ntbmax = func%n             ! total number of points (n=1 ==> r=0)
!
if (present(req_cutoff)) then
   func%cutoff = req_cutoff
   func%delta = func%cutoff / (ntbmax - 1)
else
   ! leave default cutoff in record
end if

!
do j=2,ntbmax-1
   r= func%delta*(j-1) ! value of radius
   nr=nint(r/dr)+1  ! approximate point in input grid
   ! Set up range of interpolation
   nmin=max(1,nr-npoint)
   nmax=min(np,nr+npoint)
   nn=nmax-nmin+1
   call polint(rgrid(nmin:nmax),phi(nmin:nmax),nn,r,func%f(j),dy)
   if (l /= 0 ) func%f(j) = func%f(j) / (r**l)
enddo
func%f(ntbmax) = 0.0_dp
! quadratic extrapolation to zero in uniform grid
func%f(1) = (4.0_dp * func%f(2) - func%f(3) ) / 3.0_dp

! Setup the spline information part
call rad_setup_d2(func,yp1=0.0_dp,ypn=huge(1.0_dp))

end subroutine insert_in_table
!--------------------------------------------------------

subroutine recompute_vna(spp)
use atm_types, only: species_info
use radial,    only: evaluate

type(species_info), pointer    :: spp

integer, parameter  :: NR   = 1000
integer      :: i, l, ir, nVna
real(dp)     :: rmax , dr, pop, eps, qtot
real(dp), allocatable :: charge(:), vna(:), vlocal(:), ve(:), rgrid(:)
real(dp), allocatable :: aux(:)
real(dp), parameter :: pi = 3.14159265358979_dp

allocate(charge(1:NR), vna(1:NR), ve(1:NR), vlocal(1:NR), rgrid(1:NR))
allocate(aux(1:NR))

! New grid
! Compute rmax using max orbital cutoff and chlocal cutoff
rmax = 0.0_dp
do i = 1, spp%n_orbnl
  rmax = max(rmax,spp%orbnl(i)%cutoff)
enddo
rmax = max(rmax,spp%chlocal%cutoff) 
dr = (rmax + 4.0_dp) / (NR-1)

do i=1,NR
   rgrid(i) = dr * (i-1)
enddo


! Compute electronic charge   (divided by 4pi r^2)
charge(1:NR) = 0.0_dp
qtot = 0.0_dp
do i = 1, spp%n_orbnl
   pop = spp%orbnl_pop(i)
   l = spp%orbnl_l(i)
   if (pop /= 0.0_dp) then
      qtot = qtot + pop
      aux(1:nr) = evaluate(spp%orbnl(i),rgrid(1:nr))
      if ( l /= 0) aux = aux * rgrid**l
      charge(1:nr) = charge(1:nr) + pop*aux(1:nr)*aux(1:nr) / ( 4*pi)
   endif
enddo
write(0,*) "qtot: ", qtot
!
! Compute Ve
! 
call normalize_charge(NR,rgrid(1:NR),charge(1:NR),qtot)
!
call radial_poisson(NR-1,dr,charge(1:NR),ve(1:NR))
open(unit=1,file="ELEC",form="formatted",status="unknown", &
     position="rewind", action="write")
do ir=1,NR
  write(1,*) (ir-1)*dr, charge(ir), ve(ir)
enddo
close(1)
!
! chlocal charg (divided by 4pi r^2)
!
charge(1:nr) = evaluate(spp%chlocal,rgrid(1:nr))
call normalize_charge(NR,rgrid(1:NR),charge(1:NR),-qtot)
!
! Compute Vlocal
! 
call radial_poisson(NR-1,dr,charge(1:NR),vlocal(1:NR))
open(unit=1,file="LOCAL",form="formatted",status="unknown", &
     position="rewind", action="write")
do ir=1,NR
  write(1,*) (ir-1)*dr, charge(ir), vlocal(ir)
enddo
close(1)
!
open(unit=1,file="VNA",form="formatted",status="unknown", &
     position="rewind", action="write")
do ir=1,NR
  ! Factor of 2 to convert from hartree to rydberg
  vna(ir) = 2.0_dp * (ve(ir) + vlocal(ir))
  write(1,*) (ir-1)*dr, vna(ir)
enddo
close(1)
!
! Cutoff for Vna
eps=1.0e-5_dp
do ir=NR,2,-1
  if (abs(vna(ir)).gt.eps) then
     nVna=ir+1
     exit
  endif
enddo
spp%vna%cutoff = max(dr*(nVna-1),rmax)  ! Make it as big as the max orbitals' cutoff
spp%vna%delta = spp%vna%cutoff / (spp%vna%n - 1)
do i=1,NR
   rgrid(i) = dr * (i-1)
enddo
call insert_in_table(NR,vna,rgrid,dr,spp%vna,l=0)
!
deallocate(charge, vna, ve, aux, vlocal, rgrid)

end subroutine recompute_vna
!-----------------------------------------------------------------

subroutine normalize_charge(nr,r,q,qtot)
integer, intent(in)      :: nr
real(dp), intent(in)     :: r(nr)
real(dp), intent(inout)  :: q(nr)
real(dp), intent(in)     :: qtot

real(dp), parameter :: pi = 3.14159265358979_dp

integer i
real(dp) :: qsum

qsum = 0.0_dp
do i = 1, nr
 qsum = qsum + 4*pi*r(i)*r(i)*q(i)
enddo
qsum = qsum * (r(2)-r(1))
write(0,*) "qsum: ", qsum
do i = 1, nr
 q(i) = q(i) * qtot / qsum
enddo
end subroutine normalize_charge
!--------------------------------------------------------------------

SUBROUTINE radial_poisson(NR,DELT,Q,V)
  integer, intent(in)  :: nr
  REAL(DP), intent(in)   :: delt
  REAL(DP), intent(in)   :: q(0:nr)
  REAL(DP), intent(out)  :: v(0:nr)

!   Being Q(r) a spherical charge density in a homogeneus radial mesh
!   with distance DELT between consecutive points, this routine returns
!   the electrostatic potential generated by this charge distribution.
!   Written by D. Sanchez-Portal, March 1997.

!   INTEGER NR      :    Number of radial points.
!   REAL(DP)  DELT    :    Distance between consecutive points.
!   REAL(DP)  Q(0:NR) :    Spherical charge density.
!   REAL(DP)  V(0:NR) :    Electrostatic potential at mesh points.

!   Qtot/r asimptotic behaviour is imposed.
!
!   Note that then the V units are hartree !

  integer ir
  REAL(DP) pi, fourpi, qtot, r, cons

  PI=4.0_DP*ATAN(1.0_DP)
  FOURPI=4.0_DP*PI

!     NUMEROV ALGORITHM* 
!
  V(0)=0.0_DP
  V(1)=1.0_DP

  DO IR=2,NR
     V(IR)=2.0_DP*V(IR-1)-V(IR-2) - FOURPI*DELT**3*            &
           ( Q(IR)*IR+10.0_DP*Q(IR-1)*(IR-1)+Q(IR-2)*(IR-2) )/12.0_DP
  ENDDO

! CALCULATE TOTAL CHARGE
   
  QTOT=0.0_DP
  DO IR=1,NR
     R=IR*DELT
     QTOT=QTOT+R*R*Q(IR)
  ENDDO
  QTOT=4.0_DP*PI*QTOT*DELT
  write(0,*) "Qtot: ", Qtot
  
  ! FIXING QTOT/R ASYMPTOTIC BEHAVIOUR*

  CONS=(QTOT-V(NR))/(NR*DELT)
             
  DO IR=1,NR
     R=IR*DELT
     V(IR)=V(IR)/(IR*DELT)+CONS
  ENDDO
  V(0)=(4.0_DP*V(1)-V(2))/3.0_DP

END subroutine radial_poisson

end module m_orthogonalize_orbitals
