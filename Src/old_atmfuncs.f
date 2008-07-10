! 
! This file is part of the SIESTA package.
!
! Copyright (c) Fundacion General Universidad Autonoma de Madrid:
! E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal
! and J.M.Soler, 1996-2006.
! 
! Use of this software constitutes agreement with the full conditions
! given in the SIESTA license, as signed by all legitimate users.
!
      module old_atmfuncs

C     This module contains a set of routines which provide all the information
C     about the basis set, pseudopotential, atomic mass, etc... of all the
C     chemical species present in the calculation.

!     AG: It is now a "legacy" module to interface the old and new data
!     structures.

      use precision, only: dp
      use sys
      use atmparams, only: nzetmx, lmaxd, nsemx
      use atmparams, only: maxos, nkbmx
      use radial, only : ntbmax

      implicit none 

      integer, save, public                   ::  nsmax

      integer, save, public, allocatable    ::  izsave(:)
      integer, save, public, allocatable    ::  nomax(:)
      integer, save, public, allocatable    ::  nkbmax(:)

      integer, save, public, allocatable    ::  lmxosave(:)
      integer, save, public, allocatable    ::  nkblsave(:,:)
      integer, save, public, allocatable    ::  npolorbsave(:,:,:)
      integer, save, public, allocatable    ::  nsemicsave(:,:)
      integer, save, public, allocatable    ::  nzetasave(:,:,:)
      integer, save, public, allocatable    ::  cnfigtb(:,:,:)

      logical, save, public, allocatable    ::  semicsave(:)
     
      real(dp), save, public, allocatable        :: zvaltb(:)
      integer, save, public, allocatable         :: lmxkbsave(:)
      real(dp), save, public, allocatable        :: smasstb(:)
      real(dp), save, public, allocatable        :: chargesave(:)
      real(dp), save, public, allocatable        :: slfe(:)
      real(dp), save, public, allocatable        :: lambdatb(:,:,:,:)

      real(dp), save, public, allocatable  :: qtb(:,:)

      real(dp), save, public, allocatable  :: table(:,:,:)
      real(dp), save, public, allocatable  :: tabpol(:,:,:)
      real(dp), save, public, allocatable  :: tab2(:,:,:)
      real(dp), save, public, allocatable  :: tab2pol(:,:,:)


      real(dp), save, public, allocatable  ::  coretab(:,:,:)
      real(dp), save, public, allocatable  ::  chloctab(:,:,:)
      real(dp), save, public, allocatable  ::  vlocaltab(:,:,:)
      real(dp), save, public, allocatable  ::  rctb(:,:,:)
      real(dp), save, public, allocatable  ::  rcotb(:,:,:,:)
      real(dp), save, public, allocatable  ::  rcpoltb(:,:,:,:)
!
!     PGI compiler cannot allocate these!!!!
!
      character(len=20), save, public    :: label_save(40)
      character(len=10), save, public    :: basistype_save(40)  
!--------------------------------------------------------------
!
!     Public routines
!
      public :: labelfis, izofis, zvalfis, nkbfis
      public :: massfis, lomaxfis, nofis, lmxkbfis
      public :: cnfigfio, lofio, mofio
      public :: atmpopfio , epskb, rcut
      public :: clear_tables, allocate_old_arrays
      public :: deallocate_old_arrays

!--------------------------------------------------------------

      PRIVATE

      CONTAINS !================================================
!
      subroutine allocate_old_arrays()

        allocate(rcotb(nzetmx,0:lmaxd,nsemx,nsmax))
        allocate(rcpoltb(nzetmx,0:lmaxd,nsemx,nsmax))
        allocate(lambdatb(nzetmx,0:lmaxd,nsemx,nsmax))
        allocate(qtb(maxos,nsmax))
        allocate(slfe(nsmax))
        allocate(rctb(nkbmx,0:lmaxd,nsmax))
        allocate(smasstb(nsmax))
        allocate(chargesave(nsmax))
!
!       Table: This is a hybrid
!            negative values of the second index correspond to KB projectors
!            positive values of the second index correspond to orbitals
!            A second index of zero (0) corresponds to the local potential
! 
!            The first index has ntbmax "real points" and two extra
!            entries for bookeeping
!            The total number of angular momentum entries is lmaxd+1 (since
!            s=0, p=1, etc)
!
!
        allocate(table((ntbmax+2),-nkbmx*(lmaxd+1):nzetmx*nsemx*
     .    (lmaxd+1),nsmax))
        allocate(tab2(ntbmax,-nkbmx*(lmaxd+1):nzetmx*nsemx*
     .    (lmaxd+1),nsmax))
!
        allocate(tabpol((ntbmax+2),nzetmx*nsemx*(lmaxd+1),nsmax))
        allocate(tab2pol(ntbmax,nzetmx*nsemx*(lmaxd+1),nsmax))
        allocate(coretab(ntbmax+1,2,nsmax))
        allocate(chloctab((ntbmax+1),2,nsmax))
        allocate(vlocaltab((ntbmax+1),2,nsmax))
        allocate(izsave(nsmax))
        allocate(lmxkbsave(nsmax))
        allocate(lmxosave(nsmax))
        allocate(npolorbsave(0:lmaxd,nsemx,nsmax))
        allocate(nsemicsave(0:lmaxd,nsmax))
        allocate(nzetasave(0:lmaxd,nsemx,nsmax))
        allocate(nomax(nsmax))
        allocate(nkbmax(nsmax))
        allocate(zvaltb(nsmax))
        allocate(cnfigtb(0:lmaxd,nsemx,nsmax))
        allocate(nkblsave(0:lmaxd,nsmax))
!
!       PGI compiler cannot allocate these...
!        allocate(label_save(nsmax))
!        allocate(basistype_save(nsmax))
!
        allocate(semicsave(nsmax))

        end subroutine allocate_old_arrays
!---------------------------------------------------------------
      subroutine deallocate_old_arrays()

        deallocate(rcotb)
        deallocate(rcpoltb)
        deallocate(lambdatb)
        deallocate(qtb)
        deallocate(slfe)
        deallocate(rctb)
        deallocate(smasstb)
        deallocate(chargesave)
!
        deallocate(table)
        deallocate(tab2)
!
        deallocate(tabpol)
        deallocate(tab2pol)
        deallocate(coretab)
        deallocate(chloctab)
        deallocate(vlocaltab)
        deallocate(izsave)
        deallocate(lmxkbsave)
        deallocate(lmxosave)
        deallocate(npolorbsave)
        deallocate(nsemicsave)
        deallocate(nzetasave)
        deallocate(nomax)
        deallocate(nkbmax)
        deallocate(zvaltb)
        deallocate(cnfigtb)
        deallocate(nkblsave)
        deallocate(semicsave)

        end subroutine deallocate_old_arrays

!-------------------------------------------------------
      subroutine clear_tables()

      integer is

           do is=1,nsmax
              izsave(is)=0
              lmxosave(is)=0
              lmxkbsave(is)=0
              label_save(is)='  '
              nkbmax(is)=0
              nomax(is)=0  
              semicsave(is)=.false.
              
              nsemicsave(:,is) = 0
              nzetasave(:,:,is) = 0
              rcotb(:,:,:,is) = 0.0_dp
              lambdatb(:,:,:,is) = 0.0_dp
              rcpoltb(:,:,:,is) = 0.0_dp

              table(:,:,is) = 0.0_dp
              tab2(:,:,is) = 0.0_dp
             tabpol(:,:,is) = 0.0_dp
             tab2pol(:,:,is) = 0.0_dp

             qtb(1:maxos,is)=0.00_dp

           enddo 
       end subroutine clear_tables
!----------------------------------------------------------------

      subroutine check_is(name,is)
      character(len=*), intent(in) :: name
      integer, intent(in) :: is

      if ((is.lt.1).or.(is.gt.nsmax)) then 
            write(6,*) trim(name),': THERE ARE NO DATA FOR IS=',IS
            write(6,*) trim(name),': ISMIN= 1, NSMAX= ',nsmax
         call die
      endif
      end subroutine check_is
!
!
      FUNCTION IZOFIS( IS )
      integer :: izofis ! Atomic number
      integer, intent(in) :: is ! Species index

      call check_is('izofis',is)
      izofis=izsave(is)

      end function izofis
!
      FUNCTION ZVALFIS( IS )
      real(dp) :: zvalfis          ! Valence charge
      integer, intent(in) :: is            ! Species index

      call check_is('zvalfis',is)
 
      zvalfis=zvaltb(is)
      end function zvalfis
!
      FUNCTION LABELFIS (IS)
      character(len=20) ::  labelfis  ! Atomic label
      integer, intent(in) :: is            ! Species index

      call check_is('labelfis',is)
      labelfis=label_save(is)
      end function labelfis
!
      FUNCTION LMXKBFIS (IS)
      integer :: lmxkbfis    ! Maximum ang mom of the KB projectors
      integer, intent(in) :: is            ! Species index

      call check_is('lmxkbfis',is)
      lmxkbfis=lmxkbsave(is)
      end function lmxkbfis
!
      FUNCTION LOMAXFIS (IS)
      integer :: lomaxfis  ! Maximum ang mom of the Basis Functions
      integer, intent(in) :: is            ! Species index

      integer lmx, nsm

      call check_is('lomaxfis',is)

      lomaxfis=0           
      lmx=lmxosave(is)
      do nsm=1,nsemicsave(lmx,is)+1
         if(npolorbsave(lmx,nsm,is).gt.0)   lomaxfis=lmx+1
      enddo     
      
      lomaxfis=max(lomaxfis,lmx)
      end function lomaxfis
!

      FUNCTION MASSFIS(IS)
      real(dp) :: massfis            ! Mass
      integer, intent(in) :: is            ! Species index

      call check_is('massfis',is)
      massfis=smasstb(is)
      end function massfis
!
      FUNCTION NKBFIS(IS)
      integer :: nkbfis    ! Total number of KB projectors
      integer, intent(in) :: is            ! Species index

      call check_is('nkbfis',is)
      nkbfis=nkbmax(is)
      end function nkbfis
!

      FUNCTION NOFIS(IS)
      integer :: nofis    ! Total number of Basis functions
      integer, intent(in) :: is            ! Species index

      call check_is('nofis',is)
      nofis=nomax(is)
      end function nofis

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! AMENoFIS

      FUNCTION ATMPOPFIO (IS,IO)
      real(dp) atmpopfio
      integer, intent(in) :: is    ! Species index
      integer, intent(in) :: io    ! Orbital index (within atom)

C Returns the population of the atomic basis orbitals in the atomic 
C ground state configuration.

      call check_is('atmpopfio',is)
      if((io.gt.nomax(is)).or.(io.lt.1)) then
            write(6,*) 'ATMPOPFIO: THERE ARE NO DATA FOR IO=',IO
            write(6,*) 'ATMPOPFIO: IOMIN= 1', ' IOMAX= ',nomax(is)
         call die
      endif
 
      atmpopfio=qtb(io,is) 
      end function atmpopfio

!
!
      FUNCTION CNFIGFIO(IS,IO)
      integer cnfigfio
      integer, intent(in) :: is    ! Species index
      integer, intent(in) :: io    ! Orbital index (within atom)

C Returns the valence-shell configuration in the atomic ground state
C (i.e. the principal quatum number for orbitals of angular momentum l)

C   INTEGER CNFIGFIO: Principal quantum number of the shell to what 
C                     the orbital belongs ( for polarization orbitals
C                     the quantum number corresponds to the shell which
C                     is polarized by the orbital io) 

      integer l, norb, lorb, izeta, ipol,nsm
      integer  indx, nsmorb

C
      call check_is('cnfigfio',is)
      if ((io.gt.nomax(is)).or.(io.lt.1)) then
            write(6,*) 'CNFIGFIO: THERE ARE NO DATA FOR IO=',IO
            write(6,*) 'CNFIGFIO: IOMIN= 1',
     .           ' IOMAX= ',nomax(is)
         call die
      endif

        norb=0
        indx=0
        do 10 l=0,lmxosave(is)
         do 8 nsm=1,nsemicsave(l,is)+1
          do 5 izeta=1,nzetasave(l,nsm,is)
            norb=norb+(2*l+1)
            indx=indx+1
            if(norb.ge.io) goto 30
 5        continue
 8       continue
10      continue

        indx=0
        do  20 l=0,lmxosave(is)
          do 18 nsm=1,nsemicsave(l,is)+1
            do 15 ipol=1, npolorbsave(l,nsm,is)
              norb=norb+(2*(l+1)+1)
              indx=indx+1
              if(norb.ge.io) goto 40
15          continue
18        continue
20      continue
           write(6,*) 'CNFIGFIO: ERROR: ORBITAL INDEX IO=',IO
           write(6,*) 'CNFIGFIO: ERROR: NOT FOUND'
        call die

30      lorb=l
        nsmorb=nsm
        cnfigfio=cnfigtb(lorb,nsmorb,is)
        return

40      lorb=l
        nsmorb=nsm
        cnfigfio=cnfigtb(lorb,nsmorb,is)  
        return

      end function cnfigfio
!
!
      FUNCTION LOFIO (IS,IO)
      integer lofio
      integer, intent(in) :: is    ! Species index
      integer, intent(in) :: io    ! Orbital index (within atom)

C Returns total angular momentum quantum number of a given atomic basis
C   basis orbital or Kleynman-Bylander projector.

C    INTEGER  IO   : Orbital index (within atom)
C                    IO > 0 => Basis orbitals
C                    IO < 0 => Kleynman-Bylander projectors
C************************OUTPUT*****************************************
C   INTEGER LOFIO  : Quantum number L of orbital or KB projector

      integer l, norb, izeta, ipol, nkb, nsm

      call check_is('lofio',is)
      if ((io.gt.nomax(is)).or.(io.lt.-nkbmax(is))) then 
            write(6,*) 'LOFIO: THERE ARE NO DATA FOR IO=',IO
            write(6,*) 'LOFIO: IOMIN= ',-nkbmax(is),
     .           ' IOMAX= ',nomax(is)
         CALL DIE
      endif
 
       if (io.gt.0) then

        norb=0
        do 10 l=0,lmxosave(is)
          do 8 nsm=1,nsemicsave(l,is)+1
            do 5 izeta=1,nzetasave(l,nsm,is)
              norb=norb+(2*l+1)
              if(norb.ge.io) goto 30
 5          continue
 8        continue
10      continue

        do  20 l=0,lmxosave(is)
          do 18 nsm=1,nsemicsave(l,is)+1
            do 15 ipol=1, npolorbsave(l,nsm,is)
              norb=norb+(2*(l+1)+1)
              if(norb.ge.io) goto 40
15          continue
18        continue
20      continue
          write(6,*) 'LOFIO: ERROR: ORBITAL INDEX IO=',IO
          write(6,*) 'LOFIO: ERROR: NOT FOUND'
        call die

30      lofio=l
        return

40      lofio=l+1
        return

       elseif(io.lt.0) then

        nkb=0
        do 50 l=0,lmxkbsave(is)
          do 45 izeta=1,nkblsave(l,is)
             nkb=nkb-(2*l+1)
             if(nkb.le.io) goto 60
45        continue
50      continue 

60      lofio=l       

c       elseif (io.eq.0) then
        else

        lofio=0

        endif
      end  function lofio
!
!
      FUNCTION MOFIO (IS,IO)
      integer mofio
      integer, intent(in) :: is    ! Species index
      integer, intent(in) :: io    ! Orbital index (within atom)

C Returns magnetic quantum number of a given atomic basis
C   basis orbital or Kleynman-Bylander projector.

C    INTEGER  IO   : Orbital index (within atom)
C                    IO > 0 => Basis orbitals
C                    IO < 0 => Kleynman-Bylander projectors
C************************OUTPUT*****************************************
C   INTEGER MOFIO  : Quantum number M of orbital or KB projector

      integer l, norb, izeta, ipol, nkb, lorb, lkb, nsm

      call check_is('mofio',is)
      if((io.gt.nomax(is)).or.(io.lt.-nkbmax(is))) then
            write(6,*) 'MOFIO: THERE ARE NO DATA FOR IO=',IO
            write(6,*) 'MOFIO: IOMIN= ',-nkbmax(is),
     .           ' IOMAX= ',nomax(is)
         CALL DIE
      endif

      if (io.gt.0) then

        norb=0
        do 10 l=0,lmxosave(is)
          do 8 nsm=1,nsemicsave(l,is)+1
            do 5 izeta=1,nzetasave(l,nsm,is)
              norb=norb+(2*l+1)
              if(norb.ge.io) goto 30
 5          continue
 8        continue
10      continue 

        do  20 l=0,lmxosave(is)
          do 18 nsm=1,nsemicsave(l,is)+1
            do 15 ipol=1, npolorbsave(l,nsm,is)
              norb=norb+(2*(l+1)+1)
              if(norb.ge.io) goto 40
15          continue
18        continue
20      continue
        write(6,*) 'MOFIO: ERROR: ORBITAL INDEX IO=',IO
        write(6,*) 'MOFIO: ERROR: NOT FOUND'
        call die

30      lorb=l 
        mofio=io-norb+lorb
        return

40      lorb=l+1 
        mofio=io-norb+lorb
        return


       elseif(io.lt.0) then


        nkb=0
        do 50 l=0,lmxkbsave(is)
          do 45 izeta=1,nkblsave(l,is)
             nkb=nkb-(2*l+1)
             if(nkb.le.io) goto 60
45        continue
50      continue

60      lkb=l
        mofio=-io+nkb+lkb
c       elseif (io.eq.0) then
        else

        mofio=0

        endif
        
      end function mofio
!

!  End of FIOs ----------------------------------------------------
!

      FUNCTION EPSKB (IS,IO)
      real(dp) epskb
      integer, intent(in)   ::  is   ! Species index
      integer, intent(in)   ::  io   ! KB proyector index (within atom)
                                     ! May be positive or negative 
                                     ! (only ABS(IO) is used).

C  Returns the energies epsKB_l of the Kleynman-Bylander projectors:
C       <Psi|V_KB|Psi'> = <Psi|V_local|Psi'> +
C                 Sum_lm( epsKB_l * <Psi|Phi_lm> * <Phi_lm|Psi'> )
C  where Phi_lm is returned by subroutine PHIATM.


C  REAL(DP) EPSKB  : Kleynman-Bylander projector energy
C  Energy in Rydbergs.

      integer  ik, nkb, indx, l, ikb
C
C******************************************************************
C*****************Variables in the common blocks*******************
C
      call check_is('epskb',is)
 
         ik=-abs(io)
         if (ik.eq.0) then 
             write(6,*) 'EPSKB: FUNCTION CANNOT BE CALLED WITH'
     .         ,' ARGUMENT EQUAL TO ZERO' 
           CALL DIE 
         endif

         if(ik.lt.-nkbmax(is)) then
             write(6,*) 'EPSKB: THERE ARE NO DATA FOR IO=',IK
             write(6,*) 'EPSKB: IOMIN= ',-nkbmax(is)
           CALL DIE
         endif

         nkb=0
         indx=0
         do 10  l=0,lmxkbsave(is)
             do 5 ikb=1,nkblsave(l,is)
                indx=indx+1
                nkb=nkb-(2*l+1)
                if(nkb.le.ik) goto 20
5            continue
10       continue 

20        indx=-indx
 
        epskb=table(2,indx,is)

      end function epskb
!
!
      function rcut(is,io)
      real(dp) rcut
      integer, intent(in) :: is    ! Species index
      integer, intent(in) :: io    ! Orbital index (within atom)
                                   ! io> => basis orbitals
                                   ! io<0  => KB projectors
                                   ! io=0 : Local screened pseudopotential

C  Returns cutoff radius of Kleynman-Bylander projectors and
C  atomic basis orbitals.
C  Distances in Bohr

      integer l, norb, lorb, nzetorb, izeta, ipol, nkb,lkb,nsm
      integer  indx, nsmorb        
C
      call check_is('rcut',is)

      if ((io.gt.nomax(is)).or.(io.lt.-nkbmax(is))) then
          write(6,*) 'RCUT: THERE ARE NO DATA FOR IO=',IO
          write(6,*) 'RCUT: IOMIN= ',-nkbmax(is),
     .      ' IOMAX= ',nomax(is)
        CALL DIE
      endif


       if (io.gt.0) then

        norb=0 
        indx=0
        do 10 l=0,lmxosave(is)
         do 8 nsm=1,nsemicsave(l,is)+1
          do 5 izeta=1,nzetasave(l,nsm,is)
            norb=norb+(2*l+1)
            indx=indx+1 
            if(norb.ge.io) goto 30 
 5        continue      
 8       continue 
10      continue   

        indx=0
        do  20 l=0,lmxosave(is)
          do 18 nsm=1,nsemicsave(l,is)+1
            do 15 ipol=1, npolorbsave(l,nsm,is)
              norb=norb+(2*(l+1)+1)  
              indx=indx+1  
              if(norb.ge.io) goto 40
15          continue
18        continue
20      continue     
          write(6,*) 'RCUT: ERROR: ORBITAL INDEX IO=',IO
          write(6,*) 'RCUT: ERROR: NOT FOUND'
        call die

30      lorb=l
        nzetorb=izeta
        nsmorb=nsm
        rcut=rcotb(nzetorb,lorb,nsmorb,is)  
c       rcut=table(1,indx,is)*(ntbmax-1)
        return  

40      lorb=l
        nzetorb=ipol
        nsmorb=nsm
        rcut=rcpoltb(nzetorb,lorb,nsmorb,is) 
c       rcut=tabpol(1,indx,is)*(ntbmax-1)
        return 


       elseif(io.lt.0) then 


        nkb=0
        do 50 l=0,lmxkbsave(is)
          do 45 izeta=1,nkblsave(l,is)
            nkb=nkb-(2*l+1)
            if(nkb.le.io) goto 60
45        continue
50      continue 

60      lkb=l 
        indx=-(lkb+1)
        rcut=rctb(izeta,lkb,is) 
c       rcut=table(1,indx,is)*(ntbmax-1)
        return 

c       elseif (io.eq.0) then
        else

        rcut=table(2,0,is)

        endif

      end function rcut
!
!
      end module old_atmfuncs
