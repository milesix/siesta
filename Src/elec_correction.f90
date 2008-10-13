!     This file is part of the SIESTA package.
!     
!     Copyright (c) Fundacion General Universidad Autonoma de Madrid:
!     E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal
!     and J.M.Soler, 1996-2006.
!     
!     Use of this software constitutes agreement with the full conditions
!     given in the SIESTA license, as signed by all legitimate users.
!     

      module elec_correction

      use precision
      use atom_types
      use radial
      use m_radfft
      use m_recipes, only:ratint
      implicit none

      private
      public elec_corr_setup

      contains

      subroutine elec_corr_setup

      integer is, is2, i
      real(dp) :: rchloc, rchloc2, cutoff,values(1:ntbmax)

      type(rad_func_t)  :: chlocal1, chlocal2
      type(rad_grid_t)  :: grid

      npairs = ((nspecies+1)*nspecies)/2
      allocate(elec_corr(npairs))
      
      rchloc = 0.0_dp
      rchloc2 = 0.0_dp

      do is=1,nspecies
         if (.not. is_floating(species(is))) then
            chlocal1 = get_pseudo_local_charge(species(is))
            rchloc = rad_cutoff(chlocal1)            
         endif

         do is2=is,1,-1
            
            i = ((is-1)*is)/2+is2
            
            if (is_floating(species(is)) .or. is_floating(species(is2))) then
               grid = rad_grid_alloc(ntbmax,delta=0.0001_dp)
               values=0.0_dp
               call rad_alloc(elec_corr(i),values,grid)
               !call rad_zero(elec_corr(i))
               call rad_grid_dealloc(grid)
            else
               
               chlocal2 = get_pseudo_local_charge(species(is2))
               rchloc2 = rad_cutoff(chlocal2)
               cutoff = rchloc + rchloc2 + 0.2_dp
               grid = rad_get_grid(chlocal1)
               elec_corr(i) = ch_overlap(is,is2,cutoff,grid)
               call rad_grid_dealloc(grid)
               call rad_dealloc(chlocal2)
            endif
            
         enddo
         
         if (.not. is_floating(species(is))) call rad_dealloc(chlocal1)

      enddo

      end subroutine elec_corr_setup
!
!======================================================================
!
!     This "atomic" routines use information **in the new structures**,
!     so they cannot be replaced by those in atom.f...
!
      function CH_OVERLAP(IS1,IS2,RMX,grid) result(elec_corr)
      !use atmfuncs, only: zvalfis, psch
      use parallel, only: IOnode

      integer, intent(in)   :: is1, is2
      real(dp), intent(in )    :: rmx
      type(rad_grid_t), intent(in)    :: grid
      type(rad_func_t) :: elec_corr
     

!    Returns a table with the difference between the electrostatic energy 
!    of two spherical charge-densities and two punctual charges with the 
!    same total charge as a function of the distance between the centers 
!    of these charge densities. 
!    Written by D.Sanchez-Portal. March, 1997.(from routine MATEL, written 
!    by Jose M. Soler)

!    INTEGER IS1,IS2             :  Species indexes.
!    RMX                         :  Maximum range of the correction.
!    CORR(NTBMAX)                :  Electrostatic correction energy.

!    Distances in Bohr. Energy in Rydbergs.

!    Internal precision parameters  ------------------------------------
!    NQ is the number of radial points in reciprocal space.
!    Npoint , 2npoint+1 is the number of points used by RATINT in the 
!    interpolation.
!    Q2CUT is the required planewave cutoff for the expansion of
!    the 'local-pseudopotential atomic charge density'
!    (in Ry if lengths are in Bohr).
!    CHERR is a small number to check the precision of the charge density
!    integration.

      real(dp), pointer        :: aux(:)   
      real(dp), pointer        :: corr(:)

      integer nq, npoint
      real(dp)            :: q2cut, cherr
      PARAMETER ( NQ     =  512  )
      PARAMETER ( NPOINT =  4     ) 
      PARAMETER ( Q2CUT  =  2.5e3 )
      PARAMETER ( CHERR   =  5.e-2 )

      real(dp)::  CH(0:NQ,2),VTB(NTBMAX,2), V(0:NQ,2),GRCH(3),RX(3),RAUX(2*NPOINT+1)
      

      REAL(DP) cons, qmax, rmax, delt, c, dlt, z1, z2, ch1, ch2, pi
      REAL(DP) r, vd, vv1, vv2, energ1, energ2, bessph,iz1,iz2
      integer  itb, nr, nmin, nmax, nn, iq, ir

      REAL(DP) QTMP             

      allocate(aux(1:lin_rad_default_length()),corr(1:lin_rad_default_length()))
      

      PI= 4._DP * ATAN(1._DP)       
      CONS= 1.0_dp/(2.0_dp*PI)**1.5_DP
!    
!C***  CUT-OFF IN REAL AND RECIPROCAL SPACE**
!    
      QMAX =  SQRT( Q2CUT )
      RMAX = PI * NQ / QMAX
      IF(RMX.GT.RMAX) THEN  
         if (IOnode) then
            WRITE(6,*) 'CH_OVERLAP: THE NUMBER OF INTEGRATION', &
                 ' POINTS MUST BE INCREASED'
            write(6,'(a,2f15.6)') 'ch_overlap: rmx,rmax =', rmx, rmax
         endif
         call die
      ENDIF 

      DELT=PI/QMAX
      C=4.0_DP*PI*DELT
      DLT=RMX/(NTBMAX-1)

      !IZ1=ZVALFIS(IS1)
      IZ1=get_valence_charge(species(is1))
      
      !IZ2=ZVALFIS(IS2)
      IZ2=get_valence_charge(species(is2))

      Z1=0.0_DP
      Z2=0.0_DP

      RX(2)=0.0_DP
      RX(3)=0.0_DP 

      DO IR=0,NQ
         R=IR*DELT
         RX(1)=R
             
         if (.not. is_floating(species(is1))) &
              call get_value_pseudo_local_charge(species(is1),rx,ch1,grch)
          !CALL PSCH(IS1,RX,CH1,GRCH)

         if (.not. is_floating(species(is2))) &
              call get_value_pseudo_local_charge(species(is2),rx,ch2,grch)
         !CALL PSCH(IS2,RX,CH2,GRCH)

         CH(IR,1)=-CH1
         CH(IR,2)=-CH2

         Z1=Z1-C*CH1*R*R    
         Z2=Z2-C*CH2*R*R

      ENDDO

      
      IF((ABS(Z1-IZ1).GT.CHERR).OR.(ABS(Z2-IZ2).GT.CHERR)) THEN 
         if (IOnode) then
            WRITE(6,*) 'CH_OVERLAP: THE NUMBER OF INTEGRATION', &
                 ' POINTS MUST BE INCREASED'
            WRITE(6,*) 'CH_OVERLAP: Z1=',Z1,' IZ1=',IZ1
            WRITE(6,*) 'CH_OVERLAP: Z2=',Z2,' IZ2=',IZ2
         endif
         call die
      ENDIF

      DO IR=0,NQ
         CH(IR,1)=real(IZ1,dp)*CH(IR,1)/Z1
         CH(IR,2)=real(IZ2,dp)*CH(IR,2)/Z2
      ENDDO 

!    REAL SPACE INTEGRATION OF POISSON EQUATION
!         
          
      CALL NUMEROV(NQ,DELT,CH(0,1),V(0,1))
      CALL NUMEROV(NQ,DELT,CH(0,2),V(0,2))
           
      DO ITB=1,NTBMAX
         R=DLT*(ITB-1)
         NR=NINT(R/DELT)
         NMIN=MAX(0,NR-NPOINT)
         NMAX=MIN(NQ,NR+NPOINT)
         NN=NMAX-NMIN+1
         DO IR=1,NN
            RAUX(IR)=DELT*(NMIN+IR-1) 
         ENDDO 
         CALL RATINT(RAUX,V(NMIN,1),NN,R,VV1,VD)
         CALL RATINT(RAUX,V(NMIN,2),NN,R,VV2,VD)
 
         VTB(ITB,1)=VV1
         VTB(ITB,2)=VV2
      ENDDO 
         
!      C****FOURIER-TRANSFORM OF RADIAL CHARGE DENSITY****
!      C
      CALL RADFFT( 0, NQ, RMAX, CH(0:NQ,1), CH(0:NQ,1) )
      CALL RADFFT( 0, NQ, RMAX, CH(0:NQ,2), CH(0:NQ,2) )
!C

!CNEUTRALIZE CHARGE DENSITY FOR FOURIER-SPACE CALCULATION
!C
      DO IQ=0,NQ
         R=IQ*QMAX/NQ
         CH1 = (CH(IQ,1)-IZ1*CONS)*CH(IQ,2)
         CH2=  (CH(IQ,2)-IZ2*CONS)*CH(IQ,1)
         CH(IQ,1) = CH1
         CH(IQ,2) = CH2
      ENDDO
!C
!    THE ELECTROSTATI!ENERGY CORRECTION IS STORED IN 'CORR'
! 
      DO IR=1,NTBMAX

         R=DLT*(IR-1)
         ENERG1=0.0_dp
         ENERG2=0.0_dp


         DO IQ=0,NQ
            QTMP=IQ*QMAX/NQ
            QTMP=QTMP*R 
            ENERG1=ENERG1+BESSPH(0,QTMP)*CH(IQ,1)
            ENERG2=ENERG2+BESSPH(0,QTMP)*CH(IQ,2)
         ENDDO 

         ENERG1=ENERG1*QMAX/NQ
         ENERG2=ENERG2*QMAX/NQ
   
         ENERG2=ENERG2*4.0_DP*(2.0_dp*PI)**2
         ENERG1=ENERG1*4.0_DP*(2.0_dp*PI)**2
              
         ENERG1=-(ENERG1*R)-(IZ2*(VTB(IR,1)*R-IZ1))
         ENERG2=-(ENERG2*R)-(IZ1*(VTB(IR,2)*R-IZ2))
  
         CORR(IR)=0.5_DP*(ENERG1+ENERG2)

      ENDDO 

      call rad_alloc(elec_corr,corr,grid)
      
      deallocate(corr,aux)
      END function ch_overlap


      SUBROUTINE NUMEROV(NR,DELT,Q,V)
      integer, intent(in)  :: nr
      REAL(DP), intent(in)   :: delt
      REAL(DP), intent(in)   :: q(0:nr)
      REAL(DP), intent(out)  :: v(0:nr)

!  Being Q(r) a spherical charge density in a homogeneus radial mesh
!  with distance DELT between consecutive points, this routine returns
!  the electrostati!potential generated by this charge distribution.
!  Written by D. Sanchez-Portal, March 1997.

!  INTEGER NR      :    Number of radial points.
!  REAL(DP)  DELT    :    Distance between consecutive points.
!  REAL(DP)  Q(0:NR) :    Spherical charge density.
!  REAL(DP)  V(0:NR) :    Electrostati!potential at mesh points.

!  Qtot/r asimptoti!behaviour is imposed.


      integer ir
      REAL(DP) pi, fourpi, qtot, r, cons

      PI=4.0_DP*DATAN(1.0_DP)
      FOURPI=4.0_DP*PI

!    NUMEROV ALGORITHM* 
!C
      V(0)=0.0_DP
      V(1)=1.0_DP

      DO IR=2,NR
         V(IR)=2.0_DP*V(IR-1)-V(IR-2) - FOURPI*DELT**3*( Q(IR)*IR+10.0_DP*Q(IR-1)*(IR-1)+Q(IR-2)*(IR-2) )/12.0_DP
      ENDDO 

!C***CALCULATE TOTAL CHARGE***
   
      QTOT=0.0_DP
      DO IR=1,NR
         R=IR*DELT
         QTOT=QTOT+R*R*Q(IR)
      ENDDO
      QTOT=4.0_DP*PI*QTOT*DELT

!C** FIXING QTOT/R ASIMPTOTI!BEHAVIOUR*

      CONS=(QTOT-V(NR))/(NR*DELT)
             
      DO IR=1,NR
         R=IR*DELT
         V(IR)=V(IR)/(IR*DELT)+CONS
      ENDDO 
      V(0)=(4.0_DP*V(1)-V(2))/3.0_DP

      END subroutine numerov

!-----------------------
      

      end module elec_correction










