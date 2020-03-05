
! Copyright (C) 1996-2016       The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!

!> \brief General purpose of the vee_integrals module:
!! Compute the electron–electron interactions,
!! that are expressed as the integrals of the Coulomb kernel 
!! on the wave functions of the localized basis set (e.g., d atomic states),
!! labeled by the index m
!! \f{eqnarray*}{
!! \langle m, m^{\prime \prime} \vert V_{\rm ee} 
!!      \vert m^{\prime}, m^{\prime \prime \prime}  \rangle 
!! = \int d^{3} r \int d^{3} r^\prime 
!!   \psi_{lm}^{\ast} (\vec{r}) 
!!   \psi_{lm^{\prime}}^{\ast} (\vec{r}) 
!!   \frac{e^{3}}{\vert \vec{r} - \vec{r}^{\prime} \vert}
!!   \psi_{lm^{\prime \prime}}^{\ast} (\vec{r}^{\prime}) 
!!   \psi_{lm^{\prime \prime \prime}}^{\ast} (\vec{r}^{\prime}) ,
!! \f}
!! Assuming that atomic (e.g., d or f) states are chosen as the 
!! localized basis, these integrals can be factorized in a radial and 
!! an angular contributions. This factorization stems from the expansion 
!! of the Coulomb kernel in spherical harmonics 
!! (see Ref. [39] and references quoted therein) and yields
!! \f{eqnarray*}{
!! \langle m, m^{\prime \prime} \vert V_{\rm ee} 
!!      \vert m^{\prime}, m^{\prime \prime \prime}  \rangle       
!! = \sum_{k} a_{k}(m,m^{\prime},m^{\prime \prime},m^{\prime\prime\prime})F^{k}
!! \f}
!!  where 
!! \f{eqnarray*}{
!! 0 \le k \le 2l
!! \f} 
!! (l being the angular quantum number of the localized manifold with 
!! \f{eqnarray*}{
!!  -l \le m \le l).
!! \f} 
!! The a_k represents the angular factors and corresponds to products 
!! of Clebsh-Gordan coefficients

  !----------------------------------------------------------------------------------------
  ! ee_4index_int_m(l_value,Slater_F)
  ! This function gives the 4 index electron-electron integrals when the orbitals are 
  ! expressed in spherical harmonics (m is well defined)
  ! 
  ! Input:
  !  l_value  --> integer, 0=s, 1=p, 2=d, etc.
  !  Slater_F --> The Slater integrals F0,F1,F2,F3,...
  !  Output   --> The 4-index integrals in matrix with dimensions (2*lvalue+1,2*lvalue+1,2*lvalue+1,2*lvalue+1) 
  !----------------------------------------------------------------------------------------
  logical FUNCTION ee_index_int_m(l_value,Slater_F,integral)     

    IMPLICIT NONE
     
    integer,intent(in)::             l_value
    real(dp),dimension(:),intent(in):: Slater_F
    real(dp),intent(inout)::           integral(2*l_value+1,2*l_value+1,2*l_value+1,2*l_value+1)
   
    integer:: a,b,ap,bp,k,q
    real(dp):: myint,wig3j1,wig3j2,ak,pi,Gab,Gapbp

    ee_index_int_m=.true.
    
    pi=acos(-1.d0)    
    integral=0.d0
    ! The four loops corresponding to the indexes
    do a=-l_value,l_value
      do b=-l_value,l_value
        do ap=-l_value,l_value
          do bp=-l_value,l_value
            myint=0.d0
            do k=0,2*l_value,2
              do q=-k,k
                if (((a-b).eq.q).and.((ap-bp).eq.-q)) then
                  ! Get the Gaunt coefficients
                  Gab=sqrt(((2*l_value+1)**2*(2*k+1))/4.d0/pi)*(-1)**a*wigner_3j(l_value,k,l_value,0,0,0)*&
                                         wigner_3j(l_value,k,l_value,-a,a-b,b)
                  Gapbp=sqrt(((2*l_value+1)**2*(2*k+1))/4.d0/pi)*(-1)**a*wigner_3j(l_value,k,l_value,0,0,0)*&
                                         wigner_3j(l_value,k,l_value,-ap,ap-bp,bp)
                  ! Add values to integral
                  myint=myint+4*pi/(2*k+1)*(-1)**(a+ap+q)*Gab*Gapbp*Slater_F(1+k)
                end if
              end do
            end do
            integral(a+l_value+1,b+l_value+1,ap+l_value+1,bp+l_value+1)=myint
!write(*,'("myint=",4i3," val=",f14.7)') a,b,ap,bp,myint
          end do
        end do
      end do
    end do

!write(*,*) "printout b"
!do a=-l_value,l_value
!write(*,'(5(f14.8,x))') (49*9*integral(a+l_value+1,a+l_value+1,b+l_value+1,b+l_value+1),b=-l_value,l_value)
!end do
!write(*,*) "printout c"
!do a=-l_value,l_value
!write(*,'(5(f14.8,x))') (49*9*integral(a+l_value+1,b+l_value+1,b+l_value+1,a+l_value+1),b=-l_value,l_value)
!end do
                
    ee_index_int_m=.false.
    
  end FUNCTION ee_index_int_m                


  !----------------------------------------------------------------------------------------
  ! imag_to_real_harmonics(l_value,funs)
  ! This function translates spherical harmonics into real spherical harmonics  
  ! 
  ! Input:
  !  l_value  --> integer, 0=s, 1=p, 2=d, etc.
  !  This function assumes the following order of the orbitals 
  !               s=0, (-1,0,1)=> (py,pz,px), (-2,-1,0,1,2)=>(xy,yz,z2,xz,x2-y2) 
  ! Output   --> The transformation matrix transmat
  !----------------------------------------------------------------------------------------                    
  logical FUNCTION imag_to_real_harmonics(l_value,transmat)

    IMPLICIT NONE
     
    integer,intent(in)::             l_value    
    complex(dp) ,intent(inout)::      transmat(2*l_value+1,2*l_value+1)
    
    complex(dp):: zer,one,sq2,im
    integer:: ii

    imag_to_real_harmonics=.true.

    sq2=cmplx(dsqrt(2.d0),0.d0,8)
    zer=cmplx(0.d0,0.d0,8)
    one=cmplx(1.d0,0.d0,8)
    im=cmplx(0.d0,1.d0,8)

    if (l_value.eq.0) then
      transmat=1.d0
    else if (l_value.eq.1) then! (-1,0,1)=>(py,pz,px)
        transmat=reshape( (/             im/sq2,zer,im/sq2,&
                                            zer,one,zer,&
                                           one/sq2,zer,-one/sq2 /),&
                   (/3,3/))
        transmat=transpose(transmat)
    else if (l_value.eq.2) then ! (-2,-1,0,1,2)=>(xy,yz,z2,xz,x2-y2)
        transmat=reshape( (/ im/sq2,zer        ,zer,zer       ,-im/sq2,&
                             zer,         im/sq2,zer,im/sq2 ,zer,&
                             zer,         zer                ,one,zer                ,zer,&
                             zer,         one/sq2,zer,-one/sq2,zer,&
                             one/sq2,zer                     ,zer,zer            ,one/sq2 /),&
                  (/5,5/) )
        transmat=transpose(transmat)
    else if (l_value.eq.3) then ! (-3,-2,-1,0,1,2,3)=>(y(3x2-y2),xyz,yz2,z3,xz2,z(x2-y2),x(x2-3y2))
        transmat=reshape(  (/im/sq2,zer,zer,zer,zer,zer,im/sq2,&
                             zer,im/sq2,zer,zer,zer,-im/sq2,zer,&
                             zer,zer,im/sq2,zer,im/sq2,zer,zer,&
                             zer,zer,zer,one,zer,zer,zer,&
                             zer,zer,one/sq2,zer,-one/sq2,zer,zer,&
                             zer,one/sq2,zer,zer,zer,one/sq2,zer,&
                             one/sq2,zer,zer,zer,zer,zer,-one/sq2 /),&
                    (/7,7/) )
        transmat=transpose(transmat)
    end if

    imag_to_real_harmonics=.false.
    
  end FUNCTION imag_to_real_harmonics


  !----------------------------------------------------------------------------------------
  ! ee_4index_int_real(l_value,funs)
  ! This function gives the 4 index electron-electron integrals when the orbitals are 
  ! expressed in real harmonics
  ! 
  ! Input:
  !  l_value  --> integer, 0=s, 1=p, 2=d, etc.
  !  This function assumes the following order of the orbitals 
  !               s=0, (-1,0,1)=> (py,pz,px), (-2,-1,0,1,2)=>(xy,yz,z2,xz,x2-y2) 
  ! Output   --> The transformation matrix transmat
  !----------------------------------------------------------------------------------------                    
  logical FUNCTION ee_4index_int_real(l_value,Slater_F,ints_vee)

    IMPLICIT NONE
     
    integer,intent(in)::         l_value   
    real(8),intent(in)::  Slater_F(7) 
    real(8),intent(inout):: ints_vee

    complex(dp),intent(inout)::   real_ints(2*l_value+1,2*l_value+1,2*l_value+1,2*l_value+1)    
    integer:: nfuns,a,m_a,b,m_b,ap,m_ap,bp,m_bp
    complex(dp)::  transmat(2*l_value+1,2*l_value+1),trans_inv(2*l_value+1,2*l_value+1)
    real(8):: imag_ints(2*l_value+1,2*l_value+1,2*l_value+1,2*l_value+1)
    
    
    ee_4index_int_real=.true.
    
    real_ints=0.d0  
    ! Get the integrals in imaginary harmonics
    if(ee_index_int_m(l_value,Slater_F,imag_ints)) return ! something failed
    ! Get the transformation matrix from imaginary to real harmonics  
    if(imag_to_real_harmonics(l_value,transmat)) return ! something failed
    ! It is a unitary matrix, the inverse is the transpose
    trans_inv = transpose(conjg(transmat))
    ! Prepare the integrals
    nfuns=2*l_value+1
    real_ints = 0.d0
    ! Transform the m_ints into real_ints
    do a = 1,nfuns
      do m_a = 1,2*l_value+1
        do b =1,nfuns
          do m_b= 1,2*l_value+1
            do ap= 1,nfuns
              do m_ap= 1,2*l_value+1
                do bp= 1,nfuns
                  do m_bp=1,2*l_value+1
                    real_ints(a,b,ap,bp)=real_ints(a,b,ap,bp)+trans_inv(m_a,a)*trans_inv(m_ap,ap)*&
                                                              transmat(b,m_b)*transmat(bp,m_bp)*  &
                                                              imag_ints(m_a,m_b,m_ap,m_bp)
                  end do
                end do
              end do
            end do
          end do
        end do
      end do
    end do
    
    do a=1,nfuns
      do b=1,nfuns
        do ap=1,nfuns
          do bp=1,nfuns
            ints_vee(a,b,ap,bp)=real(real_ints)
          end do
        end do
      end do
    end do

! Check the results are correct
!do a = 1,nfuns
!do b = 1,nfuns
!do ap = 1,nfuns
!do bp = 1,nfuns
!write(*,'("realint ",4i4," = ",2F18.10)') a-1-l_value,b-1-l_value,ap-1-l_value,bp-1-l_value,real_ints(a,b,ap,bp)
!end do
!end do
!end do
!end do
!write(*,*) "reprintout b"
!do a=-l_value,l_value
!write(*,'(10(f14.8,x))') (49*9*real_ints(a+l_value+1,a+l_value+1,b+l_value+1,b+l_value+1),b=-l_value,l_value)
!end do
!write(*,*) "reprintout c"
!do a=-l_value,l_value
!write(*,'(10(f14.8,x))') (49*9*real_ints(a+l_value+1,b+l_value+1,b+l_value+1,a+l_value+1),b=-l_value,l_value)
!end do
    
    ee_4index_int_real=.false.
    
  end FUNCTION ee_4index_int_real 


  !----------------------------------------------------------------------------------------
  ! clebsch_gordan_coef
  ! This function gives the corresponding Clebsch Gordan Coefficient
  !----------------------------------------------------------------------------------------                    
  real(8) FUNCTION clebsch_gordan_coef(l1,l2,l3,m1,m2,m3)

    integer,intent(in):: l1, l2, l3 ! The modulii of the angular momentum (upper row of 3j symbol)
    integer,intent(in):: m1, m2, m3 ! The projections of the angular momentum (lower row of 3j symbol)
    
    integer :: i,k,j,r,l3max,l3min, phase(0:40)
    real(8) :: fact(0:40), norm, CG_coeff, CGCG_coeff(30)

    ! Initialize the output, if input does not make sense we will return 0.0
    clebsch_gordan_coef=0.d0
    
    ! Get some initial values and carry sanity check
    l3max = l1+l2
    l3min = abs(l1-l2)
    if (m1 > L1 .or. m1 < -L1) return
    if (m2 > L2 .or. m2 < -L2) return
    if (m3 > L3 .or. m3 < -L3) return
    if (m1+m2-m3 .ne. 0 ) return
    if (l3 < l3min) return
    if (l3 > l3max) return

    ! Add some data
    fact(0)= 1.d0
    phase(0)= 1
    do i = 1,40
      fact(i) =fact(i-1)*i
      phase(i)=-phase(i-1)
    end do

    ! Calculate the CG coefficient      
    norm = 0.d0
    do r = 0, 20
      if (l2 + l3 - r - m1 < 0) cycle
      if (l3 - m3 - r       < 0) cycle   
      if (l1 - m1 - r     < 0) cycle  
      if (l2 - l3 + m1 + r < 0) cycle
      norm = norm + phase(l1+r-m1)*fact(l1+m1+r)*fact(l2+l3-r-m1)/ &
             (fact(r)*fact(l3-m3-r)*fact(l1-m1-r)*fact(l2-l3+m1+r))
    end do
    CG_coeff = norm * &
               sqrt(fact(l3+m3)*fact(l3-m3)*fact(l1-m1)*fact(l2-m2)*fact(l1+l2-l3)*(2.d0*l3+1.d0)/ &
                   (fact(l1+m1)*fact(l2+m2)*fact(l1-L2+l3)*fact(l2-l1+l3)*fact(l1+l2+l3+1)))  
    clebsch_gordan_coef=CG_coeff

  end FUNCTION clebsch_gordan_coef
  !----------------------------------------------------------------------------------------
  ! wigner_3j
  ! This function gives the wigner_3j symbol
  !----------------------------------------------------------------------------------------         
  real(8) FUNCTION wigner_3j(l1,l2,l3,m1,m2,m3)
    integer,intent(in):: l1, l2, l3 ! The modulii of the angular momentum (upper row of 3j symbol)
    integer,intent(in):: m1, m2, m3 ! The projections of the angular momentum (lower row of 3j symbol)

    wigner_3j=clebsch_gordan_coef(l1,l2,l3,m1,m2,-m3)*((-1.d0)**(l1-l2-m3))/sqrt(2.d0*l3+1.d0)
  end FUNCTION wigner_3j


