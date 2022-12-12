! 
! This file is part of the SIESTA package.
!
! Copyright (c) Fundacion General Universidad Autonoma de Madrid:
! E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal
! and J.M.Soler, 1996- .
! 
! Use of this software constitutes agreement with the full conditions
! given in the SIESTA license, as signed by all legitimate users.
!
      subroutine vmb2(nat,ttemp,mass,xa,isa,va,natoms_1)
C *************************************************************************
C This routine creates a velocity distribution according to the
C Maxwell-Boltzman distribution with a given temperature.
C It also imposes the constraint of zero total velocity 
C (to avoid center of mass drifts).
C
C Temperatures in Kelvin
C Mass in atomic mass units
C Velocities in bohr/fs
C
C Writen by J. Junquera  Nov. 96,  P. Ordejon  Nov 97.
C Modified by P. Ordejon to take into account constraints, June 2003
C ************************** INPUT ****************************************
C integer nat              : Number of atoms
C real*8  ttemp            : Temperature desired for the distribution (K)
C real*8  mass(nat)        : Atomic masses (in amu)
C real*8  xa(3,nat)        : Atomic positions
C integer isa(nat)         : Species indexes
C ************************* OUTPUT ****************************************
C real*8  va(3,nat)        : Atomic Velocities
C *************************************************************************
C
C  Modules
C
      use precision,   only : dp
      use parallel,    only : Node
      use m_fdf_global,   only: fdf_global_get

#ifdef MPI
      use mpi_siesta
#endif

      implicit none

#ifdef MPI
      integer  MPIerror
#endif

      integer  nat, isa(nat), natoms_1

      real(dp)  mass(nat), ttemp(1,2), va(3,nat), xa(3,nat)

C Internal variables ..................

      integer  iseed, i, ix, iy, ntcon

      real(dp)
     .  cell(3,3), stress(3,3), fa(3,nat), 
     .  massi, tempe(1,2), velo2, pcm(3,2), mtot(1,2),
     .  tmpstress(3,3), tmpfa(3,nat)

      external velo2
C .....................

C Only generate the velocities on one Node in case they are truely random
        call fdf_global_get(iseed,'MD.Vseed',-17)

      if (Node.eq.0) then
       write(*,*) 'siesta: Generating initial',
     &  'velocities with seed',iseed

C Loop over atoms to assign velocities .................
        mtot(1,:) = 0.0_dp
        pcm(:,:) = 0.0_dp
        !Nose 1
        do i = 1,natoms_1
          massi = mass(i)
          mtot(1,1) = mtot(1,1) + massi  
          va(1,i) = velo2(iseed,ttemp(1,1),massi)
          va(2,i) = velo2(iseed,ttemp(1,1),massi)
          va(3,i) = velo2(iseed,ttemp(1,1),massi)
          pcm(:,1) = pcm(:,1) + massi * va(:,i)        
        enddo
        !Nose 2
        do i = natoms_1+1,nat
          massi = mass(i)
          mtot(1,2) = mtot(1,2) + massi  
          va(1,i) = velo2(iseed,ttemp(1,2),massi)
          va(2,i) = velo2(iseed,ttemp(1,2),massi)
          va(3,i) = velo2(iseed,ttemp(1,2),massi)
          pcm(:,2) = pcm(:,2) + massi * va(:,i)
        enddo
C Impose constraint on zero center of mass velocity ....................
!Nose 1
        do i = 1,natoms_1
          do ix = 1,3
            va(ix,i) = va(ix,i) - pcm(ix,1)/mtot(1,1)
          enddo
        enddo
!Nose 2
        do i = natoms_1+1,nat
          do ix = 1,3
            va(ix,i) = va(ix,i) - pcm(ix,2)/mtot(1,2)
          enddo
        enddo
      endif
C Distribute the velocities over the Nodes

#ifdef MPI
      call MPI_Bcast(va(1,1),3*nat,MPI_double_precision,0,
     .  MPI_Comm_World,MPIerror)
#endif


C Impose other constraints for MD, as specified in fixed.F ..................
C Note that fixed.F imposes constraints on forces; here we
C need those on velocities, so must transform multiplying by atomic mass
C

C Zero dummy matrices for stress and cell (needed in
C call to constrains routine)

      cell = 0.0_dp
      stress = 0.0_dp
      do i = 1,nat
        do i x =1,3
          fa(ix,i) = va(ix,i)*mass(i)
        enddo
      enddo

      ! NAG's F95 compiler does not like reuse of an
      ! argument with intent(in) and intent(out), such
      ! as stress and fa in the call to fixed. Quick 
      ! workaround is to make the use of a temporary
      ! variable thus:
      call fixed( cell, stress, nat, isa, mass, xa, fa,
     .                  tmpstress, tmpfa, ntcon )
      stress = tmpstress
      fa = tmpfa
      ! But is it better to use intent(inout) in fixed? 
      ! If this is done, remember to remove the two tmp
      ! variable definitions.

      do i = 1,nat
        do ix = 1,3
          va(ix,i) = fa(ix,i)/mass(i)
        enddo
      enddo

      if (nat .le. 1) return

C Correct velocity to exact temperature .....................
      call temp2(2,nat,mass,va,ntcon,tempe,natoms_1)
!Nose 1
              if ((abs(tempe(1,1)-ttemp(1,1)) .ge. 1.e-4) .and. 
     .  (tempe(1,1) .ge. 1.e-4)) then
                do i = 1,natoms_1
                  do ix = 1,3
                    va(ix,i) = va(ix,i) * dsqrt(ttemp(1,1)/tempe(1,1))
                  enddo
                enddo
              endif
!Nose 2
             if ((abs(tempe(1,2)-ttemp(1,2)) .ge. 1.e-4) .and. 
     .  (tempe(1,2) .ge. 1.e-4)) then
        do i = natoms_1+1, Nat
          do ix = 1,3
            va(ix,i) = va(ix,i) * dsqrt(ttemp(1,2)/tempe(1,2))
          enddo
        enddo
      endif

      return
      end



      real*8 function velo2(iseed,temp,mass)
C *************************************************************************
C This function assigns velocities to atoms according to the 
C Maxwell-Boltzman distribution.
C It generates random numbers drawn from the normal distribution,
C using as input random numbers distributed uniformly from 0 to 1,  
C which are provided by routine ran3.
C
C Temperatures in Kelvin
C Mass in atomic mass units
C Velocities in bohr/fs
C
C Writen by J. Junquera  Nov. 96,  P. Ordejon  Nov 97.
C ************************** INPUT ****************************************
C integer iseed            : Seed for random number generator
C real*8  temp             : Temperature desired for the distribution (K)
C real*8  mass             : Atomic mass (in amu)
C ************************* OUTPUT ****************************************
C real*8  velo             : Velocity component 
C *************************************************************************

      implicit none

      integer 
     . iseed

      real*8
     . mass,temp

C Internal variables .................

      real*8
     .  arg1, arg2, gauss, med, pi, ran3, var
      
      external
     .  ran3

C ...........

C  For other distributions med may be different from cero.
      med = 0.0
      pi = acos(-1.0d0)
      var = sqrt(temp/mass)
C  conversion factor to bohr/fs
      var = var * 0.00172309d0
 
      arg1 = sqrt((-2.) * log(ran3(iseed)))

      arg2 = 2.0d0 * pi * ran3(iseed)
      gauss = arg1 * cos(arg2)

      velo2 = med + var * gauss

      return
      end



      subroutine temp2(iunit,natoms,ma,va,ntcon,tempe,natoms_1)
C *************************************************************************
C Subroutine to calculate instantaneous temperature
C
C Written by P.Ordejon, November'97
C ************************* UNITS ******************************************
C Temperature in Kelvin
C Atomic masses in atomic mass units
C
C Space units depend on input option:
C
C   Option iunit = 1:
C     Distances are in Angstrom
C   Option iunit = 2:
C     Distances are in Bohr
C ************************* INPUT *********************************************
C integer iunit         : Units option: 1 or 2 (see UNITS above)
C integer natoms        : Number of atoms in the simulation cell
C real*8 ma(natoms)     : Atomic masses 
C real*8 va(3,natoms)   : Atomic velocities 
C integer               : Total number of position constraints imposed
C ************************* OUTOPUT *******************************************
C real*8 tempe          : Instantaneous system temperature 
C *****************************************************************************

      use precision, only : dp
      use sys,       only : die

      implicit none

      integer 
     .   natoms,iunit,ntcon,natoms_1

      real(dp)
     .  ma(natoms),tempe(1,2),va(3,natoms)

C Internal variables ..........................................................
 
      integer
     .  i,ia

      real(dp)
     .  Ang,eV,fovermp,kin(1,2)
        tempe(1,:)=0.0_dp

C ........................

      if (iunit .ne. 1 .and. iunit .ne. 2)
     $     call die('temp: Wrong iunit option;  must be 1 or 2')

C Define constants and conversion factors .....................................

      Ang = 1.0d0 / 0.529177d0
      eV  = 1.0d0 / 13.60580d0

      if (iunit .eq. 1) then
C  convert F/m in (eV/Amstrong)/amu  to  Amstrong/fs**2
        fovermp = 0.009579038d0
      else
C  convert F/m in (Ry/Bohr)/amu  to  Bohr/fs**2
        fovermp = 0.009579038d0 *Ang**2 / eV
      endif

C ........................

C Calculate kinetic energy and temperature ...................
C Kinetic energy of atoms
      kin(1,:) = 0.0d0
!Nose 1
      do ia = 1,natoms_1
        do i = 1,3
          kin(1,1) = kin(1,1) + 0.5d0 * ma(ia) * va(i,ia)**2 / fovermp
        enddo
      enddo
!Nose 2
      do ia = natoms_1+1,natoms
        do i = 1,3
          kin(1,2) = kin(1,2) + 0.5d0 * ma(ia) * va(i,ia)**2 / fovermp
        enddo
      enddo


C Instantaneous temperature (Kelvin)
      if (natoms_1 .eq. 1) then 
        ntcon=-3 !for the single atom case
      endif

      if (iunit .eq. 1) then
        tempe(1,1) = 2.0d0*kin(1,1)/(3.0d0*natoms_1-3.0d0-ntcon)
     .  /8.617d-5
        tempe(1,2) = 2.0d0*kin(1,2)/(3.0d0*(natoms-natoms_1)-3.0d0
     .  -ntcon)/8.617d-5
      else
        tempe(1,1) = 2.0d0*kin(1,1)/(3.0d0*natoms_1-3.0d0-ntcon)
     .  /8.617d-5/eV
        tempe(1,2) = 2.0d0*kin(1,2)/(3.0d0*(natoms-natoms_1)-3.0d0
     .  -ntcon)/8.617d-5/eV
      endif

C .....................

      return
      end
    
