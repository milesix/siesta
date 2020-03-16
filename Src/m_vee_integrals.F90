
! Copyright (C) 1996-2016       The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!

!> \brief General purpose of the m_vee_integrals module:
!! Compute the electron–electron interactions,
!! that are expressed as the integrals of the Coulomb kernel 
!! on the wave functions of the localized basis set (e.g. \f$d\f$ atomic states),
!! labeled by the index \f$m\f$.
!!
!! Following Ref. \cite Himmetoglu:2014:RHC,  
!! \f{eqnarray*}{
!! \langle m, m^{\prime \prime} \vert V_{\rm ee} 
!!      \vert m^{\prime}, m^{\prime \prime \prime}  \rangle 
!! = \int d^{3} r \int d^{3} r^\prime 
!!   \psi_{lm}^{\ast} ({\bf r}) 
!!   \psi_{lm^{\prime\prime}}^{\ast} ({\bf r}^\prime) 
!!   \frac{e^{2}}{\vert {\bf r} - {\bf r}^{\prime} \vert}
!!   \psi_{lm^{\prime}} ({\bf r}) 
!!   \psi_{lm^{\prime \prime \prime}} ({\bf r}^{\prime}). 
!! \f}
!! Assuming that atomic (e.g. \f$d\f$ or \f$f\f$) states are chosen as the 
!! localized basis, these integrals can be factorized in a radial and 
!! an angular contributions. This factorization stems from the expansion 
!! of the Coulomb kernel in spherical harmonics and yields
!! \f{eqnarray*}{
!! \langle m, m^{\prime \prime} \vert V_{\rm ee} 
!!      \vert m^{\prime}, m^{\prime \prime \prime}  \rangle       
!! = \sum_{k} a_{k}(m,m^{\prime},m^{\prime \prime},m^{\prime\prime\prime})F^{k},
!! \f}
!!  where \f$ 0 \le k \le 2l \f$
!! (\f$l\f$ being the angular quantum number of the localized manifold with 
!! \f$-l \le m \le l\f$).
!!
!! The \f$a_{k}\f$ represents the angular factors and corresponds to products 
!! of Clebsh-Gordan coefficients
!! \f{eqnarray*}{
!!   a_{k}(m,m^\prime,m^{\prime\prime},m^{\prime\prime\prime}) = 
!!   \frac{4\pi}{2k+1} \sum_{q=-k}^{k} \langle lm \vert Y_{kq} 
!!                                                \vert lm^{\prime} \rangle
!!   \langle lm^{\prime \prime} \vert Y^\ast_{kq} \vert 
!!                                             lm^{\prime\prime\prime} \rangle.
!! \f}      
!! The quantities \f$ F^{k} \f$ are the Slater integrals involving 
!! the radial part of the
!! atomic wave functions \f$ R_{nl} \f$ (\f$ n \f$ indicating the atomic shell
!! they belong to). They have the following expression
!! \f{eqnarray*}{
!!   F^{k} = e^{2} \int d^{3} r \int d^{3} r^\prime r^{2} r^{\prime 2}
!!           R^{2}_{nl} ({\bf r})
!!           \frac{r^{k}_{<}}{r^{k+1}_{>}} R^{2}_{nl}({\bf r}^\prime),
!! \f}
!! where \f$ r_{<} \f$  and \f$ r_{>} \f$ indicate, respectively, 
!! the shorter and the
!! larger radial distances between \f$ r \f$  and \f$ r^{\prime} \f$.
!! For \f$ d \f$ electrons, only \f$ F^{0} \f$, \f$ F^{2} \f$, and \f$ F^{4} \f$
!! are needed to compute the \f$ V_{\rm ee} \f$ matrix elements
!! (for higher \f$ k \f$ values the corresponding \f$ a_{k} \f$ vanish)
!! while the \f$ f \f$ electrons also require \f$ F^{6} \f$.
!! 


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
module m_vee_integrals

  use precision, only: dp    ! Double precision

  implicit none

  private

  CONTAINS 

!> \brief General purpose of the function ee_index_int_m
!!
!! The integration over the product of three spherical harmonics can be 
!! simplified as a simple expression involving only a normalization 
!! constant and two Wigner 3\f$ j \f$ symbols.
!! 
!! The following expressions are borrowed from [Eq.(1.23)] of
!!
!! https://sahussaintu.files.wordpress.com/2014/03/spherical_harmonics.pdf
!!
!! So, the two products that appear in the expression needed to 
!! compute the four center matrix elements,
!! will be referred to in the function as \f$ {\tt Gab} \f$ and
!! \f$ {\tt Gapbp} \f$ respectively, and can be expressed as
!!
!! \f{eqnarray*}{
!!   {\tt Gab } \equiv 
!!   \langle lm \vert Y_{kq} \vert lm^\prime \rangle & = 
!!        \int_{0}^{2 \pi} \int_{0}^{\pi} Y^{\ast}_{l,m} (\theta, \phi) 
!!                                     Y_{k,q} (\theta, \phi) 
!!                                      Y_{l,m^{\prime}} (\theta, \phi) 
!!                                      \sin \theta d \theta d \phi
!!   \nonumber \\
!!       & = \int_{0}^{2 \pi} \int_{0}^{\pi} (-1)^{m} Y_{l,-m} (\theta, \phi) 
!!                                      Y_{k,m-m^{\prime}} (\theta, \phi) 
!!                                      Y_{l,m^{\prime}} (\theta, \phi) 
!!                                      \sin \theta d \theta d \phi
!!   \nonumber \\
!!       & = (-1)^{m} \sqrt{\frac{(2l+1)\times(2l+1)\times(2k+1)}{4 \pi}} 
!!            \left( \begin{array}{ccc}
!!                     l & k & l \\\
!!                     0 & 0 & 0 \end{array} \right)
!!            \left( \begin{array}{ccc}
!!                     l & k & l \\\
!!                     -m & m-m^{\prime} & m^{\prime} \end{array} \right)
!! \f}
!! \f{eqnarray*}{
!!   {\tt Gapbp } \equiv 
!!   \langle lm^{\prime\prime} \vert Y^{\ast}_{kq} \vert 
!!           lm^{\prime\prime\prime} \rangle & = 
!!       \int_{0}^{2 \pi} \int_{0}^{\pi} 
!!                                 Y^{\ast}_{l,m^{\prime\prime}} (\theta, \phi) 
!!                                 Y^{\ast}_{k,q} (\theta, \phi) 
!!                                 Y_{l,m^{\prime\prime\prime}} (\theta, \phi) 
!!                                 \sin \theta d \theta d \phi
!!   \nonumber \\
!!       & = \int_{0}^{2 \pi} \int_{0}^{\pi} 
!!                                 Y^{\ast}_{l,m^{\prime\prime}} (\theta, \phi) 
!!                          (-1)^q Y_{k,-q} (\theta, \phi) 
!!                                 Y_{l,m^{\prime\prime\prime}} (\theta, \phi) 
!!                                 \sin \theta d \theta d \phi
!!   \nonumber \\
!!       & = \int_{0}^{2 \pi} \int_{0}^{\pi} 
!!              (-1)^{m^{\prime\prime}} Y_{l,-m^{\prime\prime}} (\theta, \phi) 
!!                          (-1)^q Y_{k,-q} (\theta, \phi) 
!!                                 Y_{l,m^{\prime\prime\prime}} (\theta, \phi) 
!!                                 \sin \theta d \theta d \phi
!!   \nonumber \\
!!       & = \int_{0}^{2 \pi} \int_{0}^{\pi} 
!!              (-1)^{m^{\prime\prime}} Y_{l,-m^{\prime\prime}} (\theta, \phi) 
!!          (-1)^q Y_{k,m^{\prime\prime}-m^{\prime\prime\prime}} (\theta, \phi) 
!!                                 Y_{l,m^{\prime\prime\prime}} (\theta, \phi) 
!!                                 \sin \theta d \theta d \phi
!!   \nonumber \\
!!       & = (-1)^{m^{\prime\prime} + q} 
!!           \sqrt{\frac{(2l+1)\times(2l+1)\times(2k+1)}{4 \pi}} 
!!            \left( \begin{array}{ccc}
!!                     l & k & l \\\
!!                     0 & 0 & 0 \end{array} \right)
!!            \left( \begin{array}{ccc}
!!                     l & k & l \\\
!!                     -m^{\prime\prime} & 
!!                      m^{\prime\prime}-m^{\prime\prime\prime} & 
!!                      m^{\prime\prime\prime}
!!                     \end{array} \right)
!! \f}
!! 
!! In these expressions we have made use of the fact that:
!! 
!! \f$ {\tt Gab} \ne 0 \f$ if and only if 
!! \f$ q + m^{\prime} = m \Rightarrow q = m - m^{\prime}, \f$
!!
!! and
!! 
!! \f$ {\tt Gapbp} \ne 0 \f$ if and only if 
!! \f$ -q + m^{\prime\prime\prime} = m^{\prime\prime} 
!!     \Rightarrow q = m^{\prime\prime\prime} - m^{\prime\prime}. \f$
!!



  logical function ee_index_int_m(l_value,Slater_F,vee_integral)     

    implicit none
     
    integer,  intent(in)    :: l_value     ! Angular momentum of the 
                                           !    localized basis of atomic states
                                           !    (l = 2 => d-shell)
                                           !    (l = 3 => f-shell)
    real(dp), intent(in)    :: Slater_F(:) ! Value of the Slater integrals
    real(dp), intent(inout) :: vee_integral(2*l_value+1,2*l_value+1,2*l_value+1,2*l_value+1)
                                           ! Value of the four center integral
                                           !    <m m'' | V_ee | m' m'''>
                                           !    Since -l <= m <= +l,
                                           !    we need to compute and store
                                           !    (2l+1)^4 integrals.

! 
!   Internal variables
!  

    integer  :: m       ! Magnetic quantum number of the first  orbital
    integer  :: mprime  ! Magnetic quantum number of the second orbital
    integer  :: m2prime ! Magnetic quantum number of the third  orbital
    integer  :: m3prime ! Magnetic quantum number of the fourth orbital
   
    integer  :: k       ! Counter for the loop on the Clebsh-Gordan coefficients
                        !    (0 <= k <= 2l_value)
    integer  :: q       ! Counter for the loop on spherical harmonic matrix
                        !    elements
                        !    (-k <= q <= k)
    real(dp) :: pi      ! Value of pi
    real(dp) :: myint,wig3j1,wig3j2,ak,Gab,Gapbp

    ee_index_int_m=.true.
    
!
!   Initialize pi
!
    pi=acos(-1.0_dp)    

!
!   Initialize the value of the integrals
!
    vee_integral = 0.0_dp

!   The four loops corresponding to the indexes
    do m = -l_value, l_value
      do mprime = -l_value, l_value
        do m2prime = -l_value, l_value
          do m3prime = -l_value, l_value
            myint = 0.0_dp

!           Sum on k to compute the four center matrix elements
!           < m, m'' | V_ee | m' m''' > = \sum_{k} a_{k} (m, m', m'', m''') F^k
            do k = 0, 2*l_value, 2

!             For each k, we have to compute the corresponding Clebsh-Gordan
!             coefficient
!             a_k(m, m', m'', m''') = (4pi/(2k+1)) * 
!                           sum_{q=-k}^{k}  <lm|Y_{kq}|lm'> <lm''|Y*_{kq}|lm''>
              do q = -k, k

!               The integral <lm|Y_{kq}|lm'> does not vanish if and only if
!               (m-m') = q
!               The integral <lm''|Y*_{kq}|lm'''> does not vanish if and only if
!               (m''-m''') = -q
!               Only in these two cases we need to compute the product
                if ( ( (m - mprime) .eq. q ) .and.               &
 &                   ( (m2prime - m3prime) .eq. -q) ) then

!                 Get the Gaunt coefficients 
!                 (see the expressions in the documentation above)

!
!                 Gab = <lm|Y_{kq}|lm'>
!
                  Gab = dsqrt(( (2*l_value+1)**2*(2*k+1) )/(4.0_dp*pi)) *    &
 &                      (-1)**m * wigner_3j(l_value,k,l_value,0,0,0)    *    &
 &                      wigner_3j(l_value,k,l_value,-m,m-mprime,mprime)

!
!                 Gapbp = <lm''|Y*_{kq}|lm'''>
!

                  Gapbp = dsqrt(( (2*l_value+1)**2*(2*k+1) )/(4.d0*pi)) *    &
 &                       (-1)**(m2prime+q)                              *    &
 &                 wigner_3j(l_value,k,l_value,0,0,0)                   *    &
 &                 wigner_3j(l_value,k,l_value,-m2prime,m2prime-m3prime,m3prime)

!                 Add values to integral

                  myint = myint + (4*pi/(2*k+1))        *        &
! &                                (-1)**(m+m2prime+q)   *        &
 &                                Gab * Gapbp * Slater_F(1+k)
                endif
              enddo     ! Close the loop on the sum over q
            enddo       ! Close the loop on the sum over k
            vee_integral( m+l_value+1,                &
 &                        mprime+l_value+1,           &
 &                        m2prime+l_value+1,          &
 &                        m3prime+l_value+1 ) = myint

! For debugging
!          write(*,'("myint=",4i3," val=",f14.7)')  &
! &          m, mprime, m2prime, m3prime, myint
! End debugging

          enddo         ! Close the loop on m3prime
        enddo           ! Close the loop on m2prime
      enddo             ! Close the loop on mprime
    enddo               ! Close the loop on m

!write(*,*) "printout b"
!do a=-l_value,l_value
!write(*,'(5(f14.8,x))') (49*9*integral(a+l_value+1,a+l_value+1,b+l_value+1,b+l_value+1),b=-l_value,l_value)
!end do
!write(*,*) "printout c"
!do a=-l_value,l_value
!write(*,'(5(f14.8,x))') (49*9*integral(a+l_value+1,b+l_value+1,b+l_value+1,a+l_value+1),b=-l_value,l_value)
!end do
                
    ee_index_int_m=.false.
    
  end function ee_index_int_m                


!-----------------------------------------------------------------------------
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

    use precision, only : dp      ! Double precision

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

    use precision, only : dp      ! Double precision

    IMPLICIT NONE
     
    integer,intent(in)::         l_value   
    real(8),intent(in)::         Slater_F(7) 
    real(8),intent(inout)::      ints_vee

    complex(dp)  real_ints(2*l_value+1,2*l_value+1,2*l_value+1,2*l_value+1)    
!    complex(dp),intent(inout)::  real_ints(2*l_value+1,2*l_value+1,2*l_value+1,2*l_value+1)    
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
!            ints_vee(a,b,ap,bp)=real(real_ints)
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


!> \brief General purpose of the function clebsch_gordan_coef: It computes the
!!        Clebsch–Gordan coefficient
!!
!! The explicit expression for the Clebsch–Gordan coefficients, taken from
!!
!! https://en.wikipedia.org/wiki/Table_of_Clebsch–Gordan_coefficients
!!
!! is
!! \f{eqnarray*}{
!!   {\tt clebsch\_gordan\_coef } & \equiv 
!!      \langle j_{1}, j_{2} ; m_{1}, m_{2} \vert j_{1}, j_{2}, J, M \rangle =
!!      \delta_{M, m_{1} + m_{2}} 
!!      \sqrt{ \frac{(2J+1) (J + j_{1} -j_{2}) !
!!             (J - j_{1} +j_{2})! (j_{1} +j_{2} -J)! }
!!           {( j_{1} + j_{2} +J +1)!}} 
!!      \nonumber \\
!!      & \times \sqrt{ (J+M)! (J-M)! (j_{1}-m_{1})! (j_{1} + m_{1})! 
!!                      (j_{2}-m_{2})! (j_{2}+m_{2})! }
!!      \nonumber \\
!!      & \times \left( \sum_{k} \frac{(-1)^k}{k! (j_{1}+j_{2}-J-k)! 
!!                      (j_{1}-m_{1}-k)! (j_{2}+m_{2}-k)! 
!!                      (J-j_{2} + m_{1} +k)! (J - j_{1} -m_{2} +k)!}  \right)
!! \f}

  real(dp) function clebsch_gordan_coef(l1,l2,l3,m1,m2,m3)

    integer, intent(in):: l1, l2, l3 ! Angular momentum quantum numbers 
                                     !   (upper row of 3j symbol)
    integer,intent(in):: m1, m2, m3  ! Magnetic quantum numbers 
                                     !   (lower row of 3j symbol)
    
    integer :: i,k,j,r,l3max,l3min
    integer  :: phase(0:40)          ! (-1)^n, where n is the index within the
                                     !    array
    real(dp) :: pi                   ! Value of pi
    real(dp) :: fact(0:40)           ! Factorial of the first fourty integers
    real(dp) :: norm                 ! Normalization factor of the 
                                     !    Clebsch-Gordan coefficient
    real(dp) :: CG_coeff             ! Clebsch-Gordan coefficient

!   Initialize the output, if input does not make sense we will return 0.0
    clebsch_gordan_coef = 0.0_dp

    pi = acos(-1.0_dp)
    
!   Get some initial values and carry sanity check
    l3max = l1 + l2
    l3min = abs(l1-l2)
    if (m1 > l1 .or. m1 < -l1) return    ! m1 has to take a value between 
                                         !    -l1 < m1 < l1
                                         ! Same for m2 and m3 below   
    if (m2 > l2 .or. m2 < -l2) return
    if (m3 > l3 .or. m3 < -l3) return
    if (m1 + m2 - m3 .ne. 0  ) return    ! There is a Kronecker delta 
                                         !   delta_{M,m1+m2}
                                         !   in the definition of the 
                                         !   Clebsch-Gordan coefficient
    if (l3 < l3min) return
    if (l3 > l3max) return

!   Initialize some data
!   Between them, the factorial of the first fourty integers
!   and (-1) up to the index of the array
    fact(0)= 1.0_dp
    phase(0)= 1
    do i = 1, 40
      fact(i)  =  fact(i-1) * i
      phase(i) = -phase(i-1)
    enddo

!   Calculate the Clebsch-Gordan coefficient      
    norm = 0.0_dp
    do r = 0, 20
!     jjunquer
!     Version by Pablo Garcia
!      if (l2 + l3 - r - m1 < 0) cycle
!      if (l3 - m3 - r      < 0) cycle   
!      if (l1 - m1 - r      < 0) cycle  
!      if (l2 - l3 + m1 + r < 0) cycle
!      norm = norm + phase(l1+r-m1)* fact(l1+m1+r)*fact(l2+l3-r-m1)/ &
!             (fact(r)*fact(l3-m3-r)*fact(l1-m1-r)*fact(l2-l3+m1+r))
!     New version by Javier
      if ( l3 - l2 + r + m1 < 0 ) cycle
      if ( l3 - l1 + r - m2 < 0 ) cycle   
      if ( l1 + l2 - l3 - r < 0 ) cycle  
      if ( l1 - r - m1      < 0 ) cycle
      if ( l2 - r + m2      < 0 ) cycle
      norm = norm + phase(r) * 1.0_dp /             &
 &           ( fact(r) * fact(l3 - l2 + r + m1) *   &
 &             fact( l3 - l1 + r - m2 )         *   &
 &             fact( l1 + l2 - l3 - r )         *   &
 &             fact( l1 - r - m1 )              *   &
 &             fact( l2 - r + m2 ) )   
!   end jjunquer
    end do
!    jjunquer
!    Version by Pablo Garcia
!    CG_coeff = norm * &
!               sqrt( fact(l3+m3) * fact(l3-m3) * fact(l1-m1) *      &
! &               fact(l2-m2)* fact(l1+l2-l3)* (2.dp*l3+1.dp) /      &
! &             ( fact(l1+m1) * fact(l2+m2) * fact(l1-l2+l3)  *      &
! &               fact(l2-l1+l3) * fact(l1+l2+l3+1) ) )  
!   New version by Javier
    CG_coeff = norm * &
               dsqrt( fact(l1+m1) * fact(l1-m1)                      * &
 &                    fact(l2+m2) * fact(l2-m2)                      * &
 &                    fact(l3+m3) * fact(l3-m3) )                    * & 
 &             dsqrt( (fact(l1 + l2 - l3) * fact(l2 + l3 - l1)       * &
 &                     fact(l3 + l1 - l2 ) ) / (l1 + l2 + l3 + 1) )  * &
 &             dsqrt( 2.0_dp*l3 + 1.0_dp ) 
    clebsch_gordan_coef = CG_coeff

  end function clebsch_gordan_coef

!> \brief General purpose of the function wigner_3j: It computes the
!!        Wigner_3j symbols
!!
!! The 3-j symbols are given in terms of the Clebsch–Gordan coefficients by
!! 
!! The following expressions are borrowed from Wikipedia
!!
!! https://en.wikipedia.org/wiki/3-j_symbol
!!
!! \f{eqnarray*}{
!!   {\tt wigner\_3j } & \equiv 
!!            \left( \begin{array}{ccc}
!!                     j_{1} & j_{2} & j_{3} \\\
!!                     m_{1} & m_{2} & m_{3} \end{array} \right)
!!            \equiv \frac{(-1)^{j_{1}-j_{2}-m_{3}}}{\sqrt{2j_{3}+1}}
!!            \langle j_{1} m_{1} j_{2} m_{2} \vert j_{3} (-m_{3}) \rangle,
!! \f}
!! where \f$ \langle j_{1} m_{1} j_{2} m_{2} \vert j_{3} (-m_{3}) \rangle\f$
!! is a Clebsch–Gordan coefficient

  real(dp) function wigner_3j(l1,l2,l3,m1,m2,m3)
    integer, intent(in):: l1, l2, l3 ! Angular momentum quantum numbers 
                                     !   (upper row of 3j symbol)
    integer,intent(in):: m1, m2, m3  ! Magnetic quantum numbers 
                                     !   (lower row of 3j symbol)

    wigner_3j = clebsch_gordan_coef(l1,l2,l3,m1,m2,-m3) *   &
 &     ((-1.0_dp)**(l1-l2-m3))/dsqrt(2.0_dp*l3+1.0_dp)

  end function wigner_3j


end module m_vee_integrals
