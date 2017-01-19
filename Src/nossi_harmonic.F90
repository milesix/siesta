module nossi_harmonic
! Copyright  (c) 2010,2011,2012 INRIA    
!
! author O. Coulaud - Olivier.Coulaud@inria.fr 
!
! This software is a computer program whose purpose is to [descr1ibe
! functionalities and technical features of your software].
!
! This software is governed by the CeCILL C license under French law and
! abiding by the rules of distribution of free software.  You can  use, 
! modify and/ or redistribute the software under the terms of the CeCILL
! license as circulated by CEA, CNRS and INRIA at the following URL
! "http://www.cecill.info". 
!
! As a counterpart to the access to the source code and  rights to copy,
! modify and redistribute granted by the license, users are provided only
! with a limited warranty  and the software's author,  the holder of the
! economic rights,  and the successive licensors  have only  limited
! liability. 
! 
! In this respect, the user's attention is drawn to the risks associated
! with loading,  using,  modifying and/or developing or reproducing the
! software by the user in light of its specific status of free software,
! that may mean  that it is complicated to manipulate,  and  that  also
! therefore means  that it is reserved for developers  and  experienced
! professionals having in-depth computer knowledge. Users are therefore
! encouraged to load and test the software's suitability as regards their
! requirements in conditions enabling the security of their systems and/or 
! data to be ensured and,  more generally, to use and operate it in the 
! same conditions as regards security. 
!------------------------------------------------------------------------------------
!
! The fact that you are presently reading this means that you have had
! knowledge of the CeCILL license and that you accept its terms.                               
!
  integer, private             :: P_multipole
  real(8), private,allocatable :: Anm(:),values(:)
!
  public :: init_harmonic, evalylm
contains
  subroutine init_harmonic(num_coeffMultipole,P)
    !
    integer, intent(in) :: num_coeffMultipole,P

    integer :: ind, n, m
    !
    P_multipole = P
    if(.NOT. allocated(Anm) ) then 
       allocate ( Anm( ( P+1)* ( P+2)/2), values(num_coeffMultipole) )
    endif
    !
    ! Set up the coefficients Anm = sqrt( (n-m)!/(n+m)! ) m = 0,...,n
    !
    Anm(1)  = 1 ;  Anm(2)  = 1 ; Anm(3)  = 2 
    ! start for n = 2
    ind = 3
    do n = 2 , P_multipole
       ind = ind + 1
       Anm(ind) = 1 
       do m = 1,  n
          Anm(ind+1) = Anm(ind)*((n+m)*(n-m+1))
          ind = ind + 1
       end do
    end do
    do m = 1, size(anm)
       Anm(m) = 1.0_8/sqrt(Anm(m))
    end do
  end subroutine init_harmonic
  pure subroutine cart2sph(x,y,z,r,costheta,phi)
    real(8), intent(in)  :: x,y,z
    real(8), intent(out) :: r,costheta,phi
    !
    real(8) :: r2
    !
    r2 =  x**2 + y**2
    r  = sqrt(r2 + z**2 )
    if( r > 1.0e-10_8) then
       costheta  =  z/r 
       if(r2  > 1.0e-10_8 ) then 
          phi = atan2(y,x) 
       else
          phi = 0.0_8
       endif
    else
       costheta  = 1.0_8  ; phi = 0.0_8
    end if
    !
  end subroutine cart2sph
#if DEBUG
  subroutine polyLegendre(lmax,x,values)
#else
    pure subroutine polyLegendre(lmax,x,values)
#endif
    !
    ! Compute associated  legendre polynom for l >= 0
    !
    implicit none
    integer, intent(in) :: lmax
    real(8), intent(in) :: x
    real(8), intent(out) :: values(:)
    !
    real(8)  :: pmm, somx2,fact,plgndr,pmmp1,pll,sign
    integer  :: ind,m,ll,ind1,ind2,l  ,i,fact0,fact1,fact2
!    logical :: bug
    !
    ! p_l^m is stored i, (l+1)*(l+2)/2 + m+1 
    !
    values(:) = 0.0_8
#if DEBUG
    if ( size(values) < (lmax+1)**2 ) then 
       write(*,*) "Wrong size of values ", size(values)," <", (lmax+1)**2 
       stop "in  polyLegendre subroutine"
    end if
    if (abs(x) > 1.0_8) then 
       write(*,*) "Wrong value of x - must be less than 1 ",x,x/abs(x)
       write(120,'(A,3E14.5)') 'polyLegendre::Error greater than 1 ',x,x-1,x+1
!       x = x/abs(x)
       write(120,'(A,3E14.5)') 'polyLegendre::Error greater than 1 ',x,x-1,x+1
!       stop 'in subroutine polyLegendre'
    end if
#endif
    !
    !  Step 1 : compute P_m^m = (-1)^m (2m-1)!! (1-x^2)^(m/2)
    !    They are stored in values  at (m+1)(m+2)/2
    !    
    values(1) = 1.0_8
    !
    somx2 = sqrt(1.0_8 -x*x) 
    fact  = 1.0_8 
    ind = 1   
    do i=1,P_multipole
       ind1 = ind + i + 1
       values(ind1)  = - values(ind)*fact*somx2
       fact = fact + 2.0_8
       ind = ind1
    enddo
    !
    ! step 2 Compute P_(m+1)^m = x (2m+1) P_(m)^m
    !   stored at (m+1)*(m+4)/2
    ind = 0 ; ind1 = 0
    do i=0,P_multipole-1
       ind  = ind  + (i+1) 
       ind1 = ind1 + (i+2) 
       values(ind1) = x*real(2*i+1,8)*values(ind)
    end do
    !
    ! step 3 Compute P_ (l)^m <-  P_(l-1)^m ,  P_ (l-2)^m
    !    (l-m) P_l^m(x) = (2*l -1) x  P_(l-1)^m(x)  -(l-1+m) P_(l-2)^m(x)
    !
    do m = 0, P_multipole-2
       !  start with P_(m)^m  ind2 = (m+1)(m+2+1)
       ind =   (m+1)*(m+2)/2 ; ind1 = ind + (m+1)
       do l = m+2, P_multipole  ! start with l-2 = m ! l - 1 stored l*(l-1) + m+1 
           ind2 = ind1 + (l-1)  + 1
          values(ind2)   = (x*(2*l-1)*values(ind1)-(l+m-1)*values(ind))/(l-m) 
          ind   = ind1 ; ind1 = ind2
       enddo
    enddo
#if DEBUG
    do i = 1, size(values) 
       if (values(i) /= values(i) ) then 
          print*, 'Not a Number', i, values(i)
          stop 'in subroutine polyLegendre'
       end if
    end do
#endif
  end subroutine polyLegendre
  subroutine evalylm(x,y,z,ylm) 
    !
    ! Our definition of spherical harmonics coincides with that of Epton and Dembart 
    !
    implicit none
    real(8), intent(in)  :: x,y,z
    real(8),intent(out) :: ylm(:)
    !
    integer :: l,m,ind,ind1
    real(8), parameter :: sqrt2 =sqrt(2.0_8) , pi=4.0_8*atan(1.0) 
    real(8) :: r,cosinus,sinus,theta,phi,coeff,sign
    !
    call cart2sph(x,y,z,r,cosinus,phi)
    !    write(260,'(A,4E14.5)')  '  eval --       (r,CosTheta,SinTheta,phi): ', r,cosinus,sinus, phi
    !
    call  polyLegendre(P_multipole,cosinus,values)
    !
    ind1 = 1
    do l = 0, P_multipole
       ind =  l*(l+1) + 1 ! position pour l et m = 0
       ylm(ind) = values(ind1)*Anm(ind1)
       ind1 = ind1 + 1
       do m = 1, l
          coeff      = sqrt2*values(ind1)*Anm(ind1)
          ylm(ind+m) = coeff*cos(m*phi)
          ylm(ind-m) = coeff*sin(m*phi)
          ind1       = ind1 + 1 
      end do
    end do
  end subroutine evalylm
end module nossi_harmonic
