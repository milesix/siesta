      Program SIESTA_QMMM 

C****************************************************************************
C Program for QM-MM calculations. 
C It uses the SIESTA code to treat at DFT level the QM subsystem. 
C It uses an own implementation of Amber99 force fiel parametrization
C to treat the MM subsystem. 
C Original idea by D. Scherlis, D. Estrin and P. Ordejon. 2001.
C Original QM-MM interfase with Siesta by D. Scherlis and P. Ordejon. 2001.
C Original subroutines by D. Scherlis. 2001/2002.
C Solvent implementation using the Amber99 force field parametrization
C by A. Crespo and M. Marti. subrouitnes by M. Marti. 2003.
C Further modifications of QM-MM Siesta by A. Crespo and P. Ordejon. 2004.
C Developing of SIESTA-QMMM: programf or QM-MM calculations using the SIESTA
C as SUBROUTINE by A. Crespo and P. Ordejon. 2004.
C This code uses several subroutines of SIESTA  as well as the 
C SiestaAsSubroutine and FDF packages.
C Crespo, May 2005. 
C C. F. Sanz-Navarro (2008-2012)
C*****************************************************************************

C Modules
      use m_siesta_init
      use m_steps
      use atomlist, only: rmaxo, no_s, iaorb, iphorb, amass, datm, 
     .                    indxua
      use atmfuncs,      only : nofis, rcut
      use mesh,          only : xdsp, nsm, nsp
      use siesta_geom
      use m_siesta_move
      use alloc,         only : re_alloc

      use m_intramol_pressure
      use m_kinetic, only : tempion
      use m_energies,   only: Ekinion, Etot
      use m_stress
      use fsiesta
      use precision
      use m_timestamp, only : timestamp
      use sys
      use fdf
      use m_qmmm_fdf
      use parallel, only: Node, Nodes, IOnode    
      use files,       only : slabel

      use siesta_cml, only : cml_p
      USE siesta_options
      use m_forces
      use zmatrix,      only: lUseZmatrix, iofaZmat,
     .                        CartesianForce_to_ZmatForce
      use m_broyden_optim, only: broyden_optimizer
      use m_fire_optim, only: fire_optimizer
      use m_zm_broyden_optim,      only: zm_broyden_optimizer
      use m_zm_fire_optim,      only: zm_fire_optimizer
      use m_cell_broyden_optim,      only: cell_broyden_optimizer
      use m_cell_fire_optim,      only: cell_fire_optimizer
      use m_iorho,       only : check_rho, write_rho, read_rho
      use iogrid_netcdf, only: write_grid_netcdf
      use write_subs, only : siesta_write_forces, 
     .            siesta_write_stress_pressure
#ifdef QMMM_BSC
      use m_fixed, only : fixed
      use meshsubs, only : initMesh, sameMeshAndAtoms
#endif
      use m_state_init
      use m_state_analysis
      use m_siesta_analysis
      use m_siesta_end

      use siesta_qmmm_options
      use m_fill_auxrho
      use linkatoms
      use qmmm_neighbour

C General Variables
      implicit none

#ifndef QMMM_BSC
      include 'fdf/fdf.h'                                                                 
#endif
                                                                                      
      logical :: overflow

      integer ::
     .     i, ia, 
     .     istep, iunit, ix, j, k, jx, ntm(3), 
     .     nspec, ntpl, qmmm_nsm,  step, 
     .     iadispl, ixdispl, nspin, ntplext, nspinext

      integer :: restart_istep = 0

      integer :: nm(3), nml(3), ntml(3), nmpl, ntopl

      integer :: io, iphi, is, noa 

      real(dp) :: alat, g2max, G2mesh

      real(dp) :: volume, volcel

      logical :: harrisdens = .false.
      logical :: state_analaysis_flag=.false.

      logical :: null_mm_charges=.true.

      real(dp) dvol

C      real(dp) :: volcel

      logical :: reorder_rho = .false.
      real(grid_p), dimension(:), pointer :: Vqm, Rho, 
     .        resonance_Rho, auxrho, Chlocal
      integer :: resonance_id

      integer :: ns          ! Number of species
      integer, dimension(:), allocatable :: isa_qm,iza_qm
      real(dp), dimension(:), allocatable :: atm_rmaxo

      character, dimension(:), allocatable ::
     .  atsym*2

      logical
     .  found, foundxv, relaxd, qm , mm, lex

      logical :: IsDiag, samesh, samexa

      real(dp) :: cell(3,3)

      character
     .  siesta_qmmmslabel*30, siestaslabel*40,
     .  filerho*30, filevqm*30, paste*30, coulombtype*10,
     .  unix_command*80

      real(dp) ewald_alpha
      real(dp) kewald_cutoff

      parameter ( qmmm_nsm = 2 )

      external qmmm_ioxv, iorho, qmmm_fixed1,  
     .  assign, iofa, memory, paste, qmmm_reinit, read_qm, 
     .  read_md, qmmm_timer, write_lab, qmmm_fixed2

C Solvent General variables

      character*4,  dimension(:), allocatable, save ::
     .   atname,aaname,attype,qmattype

      character*5, dimension(:), allocatable, save::
     .  bondtype

      character*8, dimension(:), allocatable, save::
     .  angletype

      character*11, dimension(:), allocatable, save::
     .  dihetype,imptype

      real(dp), dimension(:), allocatable, save ::
     .  kbond,bondeq,kangle,angleeq,kdihe,diheeq,
     .  perdihe,kimp,impeq,perimp,desplaz

      integer, dimension(:,:,:), allocatable, save ::
     .  atange,atangm,atdihe,atdihm,atimp

      integer, dimension(:,:), allocatable, save :: ng1

      integer, dimension(:), allocatable, save ::
     .  aanum,multidihe,multiimp, bondxat,angexat,
     .  dihexat,dihmxat,angmxat,impxat,
     .  scalexat,nonbondedxat,izs, qmstep

      integer, dimension(:), allocatable :: graphite_layer_no

      integer na_qm, na_mm,natot,nfree,nparm, imm, 
     .  nroaa, nbond, nangle, ndihe, nimp , 
     .  atxres(20000), totcoor

      real(dp)  Etots, Etot_amber,sfc,timestep,kcal,rcorteqm,
     .  rcortemm,Elj,Eelec

      real(dp) xmax, xmin

      real(dp), dimension(:,:), allocatable, save :: fce_amber, rclas

      real(dp), dimension(3,3) :: stress_amber

      real(dp), dimension(:,:), allocatable, save :: rqm, rmm

      real(dp), dimension(:), allocatable, save :: Em, Rm,
     .  pc, masst

      integer, dimension(:,:,:), allocatable, save :: evaldihe,
     .  evaldihm
      integer, dimension(:,:), allocatable, save :: ng1type, angetype,
     .  angmtype,dihety,dihmty,impty,nonbonded,scale
      logical, dimension(:,:), allocatable, save :: evaldihelog,
     .  evaldihmlog

      logical firstpc,constropt,actualiz, free, deuter 

      integer :: ncount
      real(dp) :: fovermp, kin, tempe

C ConstrOpt variables
        integer nconstr,typeconstr(20),atmsconstr(20,20)
        integer nstepconstr,istepconstr, ndists(20)
        real(dp) kforce(20),ro(20),rt(20),rini,rfin,dr,
     .                   coef(20,10)
 
C Cut Off QM-MM variables
        real(dp) rcorteqmmm, radbloqmmm
        integer, dimension(:), allocatable, save :: blocklist,
     .  blockqmmm,listqmmm

C more variables
        logical water,foundcrd,foundvat,foundva,writeipl,debug 

C Solvent external variables
        external solv_assign, solv_ene_fce, qmmm_lst_blk, wrtcrd,
     .  centermol, centerdyn, pcpot, ljef,
     .  mmforce, readcrd,  subfree, prversion, ioxvconstr, elecfield, 
     .  wripdb, wriene, wrirtc, subconstr1, subconstr2, subconstr3
C--------------------------------------------------------------------

C Periodic boundary conditions
        character lattice_type, get_lattice_type

        character(len=13) :: string_path

        character(len=20) filein, fileout

C Initialise some variables
        relaxd=.false.
C        varcel=.false.
        ntcon=0
        tempion=0.0
        tp=0.0
        vcell=0.0
        strtol=0.0
        bulkm=0.0
        ntm=0
        ntpl=1
C     nml: Number of Mesh divisions of each cell vector (in serial: (1,1,1)
        nml(1:3)=1
        qm=.true.
        mm=.true.
        lattice_type='D'

C Initialize IOnode
      IOnode = .true. ! So far only serial version of siesta-qmmm

C Print version information
      call prversion
      write(6,'(/,a,i4,a)') 'siesta-qmmm: Running in serial mode' 

C Start time counter 
      call timestamp('Start of run')      
      call qmmm_timer( 'all', 0 )
      call qmmm_timer( 'Siesta-qmmm', 1 )
      call qmmm_timer( 'Setup', 1 )

C Factors of conversion to internal units 
      kcal   = 1.602177E-19_dp * 6.02214E23_dp / 4184.0_dp

      call qmmm_reinit(sname)
      siesta_qmmmslabel=slabel

C Read the number of QM atoms
      na_qm=fdf_integer('NumberOfAtoms',0)
      if (na_qm.eq.0) then
c     call die("You must specify number of QM atoms")
         write(6,'(/a)') 'siesta-qmmm: Running with no QM atoms'
         qm=.false.
      else
C     Initialize slabel, filerho and filevqm 
         siestaslabel = paste( siesta_qmmmslabel, '.siesta' )
      endif

C     Read the number of species
      nspec = fdf_integer('NumberOfSpecies',0)
      if(qm.and.nspec.eq.0) then
         call die("siesta-qmmm: You must specify the number of species")
      endif

      if(qm) then
C     Read QM coordinates and labels
         write(6,*)
         write(6,"('read: ',73(1h*))")
         allocate(rqm(3,na_qm))
         call memory('A','D',3*na_qm,'siesta-qmmm')  
         allocate(isa_qm(na_qm))
         call memory('A','D',na_qm,'siesta-qmmm')
         allocate(iza_qm(na_qm))
         call memory('A','D',na_qm,'siesta-qmmm')
         allocate(atsym(nspec))
         call memory('A','D',nspec,'siesta-qmmm')
         write(6,*)
         write(6,"('read: ',73(1h*))")
         call read_qm(na_qm,nspec,isa_qm,iza_qm,rqm,atsym)
         allocate(qmattype(na_qm))
         call memory('A','D',na_qm,'siesta-qmmm')
      endif                     !qm

C     Read the number of MM atoms
      na_mm = fdf_integer('NumberOfSolventAtoms',0)
      if (na_mm.eq.0) then
c     call die("You must specify number of solvent atoms")
         write(6,'(/a)') 'siesta-qmmm: Running with no MM atoms'
         mm=.false.
      endif

      call redcel( alat, ucell, scell, mscell )

      lattice_type=get_lattice_type(ucell)

C Read and assign Solvent variables 
      if(mm) then

         allocate(rmm(3,na_mm))
         call memory('A','D',3*na_mm,'siesta-qmmm')
         allocate(pc(na_mm))
         call memory('A','D',na_mm,'siesta-qmmm')
         allocate(attype(na_mm))
         call memory('A','D',na_mm,'siesta-qmmm')
         allocate(atname(na_mm))
         call memory('A','D',na_mm,'siesta-qmmm')
         allocate(aaname(na_mm))
         call memory('A','D',na_mm,'siesta-qmmm')
         allocate(aanum(na_mm))
         call memory('A','D',na_mm,'siesta-qmmm')
         allocate(ng1(na_mm,6))
         call memory('A','D',6*na_mm,'siesta-qmmm')
         allocate(graphite_layer_no(na_mm))
         call memory('A','D',na_mm,'siesta-qmmm')
         allocate(blockqmmm(na_mm))
         call memory('A','D',na_mm,'siesta-qmmm')
         blockqmmm = 0
         allocate(listqmmm(na_mm))
         call memory('A','D',na_mm,'siesta-qmmm')
         listqmmm = 0
         allocate(fce_amber(3,na_mm))
         call memory('A','D',3*na_mm,'siesta-qmmm')
         allocate(ng1type(na_mm,6))
         call memory('A','D',6*na_mm,'siesta-qmmm')
         allocate(angetype(na_mm,25))
         call memory('A','D',25*na_mm,'siesta-qmmm')
         allocate(angmtype(na_mm,25))
         call memory('A','D',25*na_mm,'siesta-qmmm')
         allocate(evaldihe(na_mm,100,5))
         call memory('A','D',500*na_mm,'siesta-qmmm')
         allocate(evaldihm(na_mm,100,5))
         call memory('A','D',500*na_mm,'siesta-qmmm')
         allocate(dihety(na_mm,100))
         call memory('A','D',100*na_mm,'siesta-qmmm')
         allocate(dihmty(na_mm,100))
         call memory('A','D',100*na_mm,'siesta-qmmm')
         allocate(impty(na_mm,25))
         call memory('A','D',25*na_mm,'siesta-qmmm')
         allocate(nonbonded(na_mm,100))
         call memory('A','D',100*na_mm,'siesta-qmmm')
         allocate(scale(na_mm,100))
         call memory('A','D',100*na_mm,'siesta-qmmm')
         allocate(evaldihelog(na_mm,100))
         call memory('A','D',100*na_mm,'siesta-qmmm')
         allocate(evaldihmlog(na_mm,100))
         call memory('A','D',100*na_mm,'siesta-qmmm')
         allocate(desplaz(na_mm))
         call memory('A','D',na_mm,'siesta-qmmm')
         allocate(scalexat(na_mm))
         call memory('A','D',na_mm,'siesta-qmmm')
         allocate(nonbondedxat(na_mm))
         call memory('A','D',na_mm,'siesta-qmmm')   
         nparm = 500
         allocate(kbond(nparm),bondeq(nparm),
     .        bondtype(nparm))
         call memory('A','D',3*nparm,'siesta-qmmm')
         allocate(kangle(nparm),angleeq(nparm),
     .        angletype(nparm))
         call memory('A','D',3*nparm,'siesta-qmmm')
         allocate(kdihe(nparm),diheeq(nparm),
     .        dihetype(nparm),multidihe(nparm),
     .        perdihe(nparm))
         call memory('A','D',5*nparm,'siesta-qmmm')
         allocate(kimp(nparm),impeq(nparm),
     .        imptype(nparm),multiimp(nparm),
     .        perimp(nparm))
         call memory('A','D',5*nparm,'siesta-qmmm')
         allocate(atange(na_mm,25,2))
         call memory('A','D',50*na_mm,'siesta-qmmm')
         allocate(atangm(na_mm,25,2))
         call memory('A','D',50*na_mm,'siesta-qmmm')
         allocate(atdihe(na_mm,100,3))
         call memory('A','D',100*na_mm,'siesta-qmmm')
         allocate(atdihm(na_mm,100,3))
         call memory('A','D',100*na_mm,'siesta-qmmm')
         allocate(bondxat(na_mm))
         call memory('A','D',na_mm,'siesta-qmmm') 
         allocate(angexat(na_mm))
         call memory('A','D',na_mm,'siesta-qmmm')
         allocate(dihexat(na_mm))
         call memory('A','D',na_mm,'siesta-qmmm')
         allocate(dihmxat(na_mm))
         call memory('A','D',na_mm,'siesta-qmmm')
         allocate(angmxat(na_mm))
         call memory('A','D',na_mm,'siesta-qmmm')
         allocate(impxat(na_mm))
         call memory('A','D',na_mm,'siesta-qmmm')
         allocate(atimp(na_mm,25,4))
         call memory('A','D',100*na_mm,'siesta-qmmm')
         allocate(Em(na_qm+na_mm))
         call memory('A','D',na_qm+na_mm,'siesta-qmmm')
         allocate(Rm(na_qm+na_mm))
         call memory('A','D',na_qm+na_mm,'siesta-qmmm')

         call solv_assign(na_qm,na_mm,nroaa,Em,Rm,attype,pc,
     .        ng1,bondxat,angexat,atange,angmxat,atangm,dihexat,atdihe,
     .        dihmxat,atdihm,impxat,atimp,
     .        nbond,kbond,bondeq,bondtype,
     .        nangle,kangle,angleeq,angletype,
     .        ndihe,kdihe,diheeq,dihetype,multidihe,perdihe,
     .        nimp,kimp,impeq,imptype,multiimp,perimp,
     .        nparm,aaname,atname,aanum,qmattype,rmm,
     .        graphite_layer_no,
     .        rcorteqmmm,rcorteqm,rcortemm,sfc,timestep,
     .        radbloqmmm,atxres,ucell,lattice_type)

C Substract the MM charge from the total charge in order to
C calculate the charge on the QM region.

         null_mm_charges=.true.
         do i=1,na_mm
            if (pc(i)/=0.0_dp) null_mm_charges=.false.
            charnet=charnet-pc(i)
         enddo

         if (abs(charnet)<0.0010) charnet=0.0d0  

      endif                     !mm
C     
      debug = fdf_boolean('SolventDebug',.false.)

      linkatom=.false.

      if (qm.and..not.mm) then

         allocate(resonance(num_resonances))   

         resonance(1)%path='./'
         resonance(1)%weight=1.0_dp

      else if (qm.and.mm) then

         call get_numlink(na_qm,na_mm,rqm,rmm,attype,qmattype,nparm,
     .           nbond,bondtype,bondeq,debug,ucell,lattice_type)

         allocate(resonance(num_resonances))   

         resonance(1)%path='./'
         resonance(1)%weight=1.0_dp

         if (numlink>0) then

            linkatom=.true.

            do i=1,num_resonances
               allocate(resonance(i)%atype(numlink))
               allocate(resonance(i)%isa(numlink))
               allocate(resonance(i)%iza(numlink))
            enddo

            if (num_resonances>1) then
               do i=2,num_resonances
                  write(string_path,'(a,i2.2,a)')'RESONANCE-',i,'/'
                  resonance(i)%path=trim(string_path)
                  call system('mkdir -p '//trim(resonance(i)%path))
                  call system('cp *.psf '//trim(resonance(i)%path))
               enddo
            endif

            do i=1,num_resonances
               do j=1,numlink
                  resonance(i)%iza(j)=0 
               enddo
            enddo

!     Allocate LinkAtom variables
            allocate(rlink(num_resonances,3,numlink))
            call memory('A','D',num_resonances*3*numlink,'siesta-qmmm')  
            rlink=0.0_dp
            allocate(fa_link(num_resonances,3,numlink))
            call memory('A','D',num_resonances*3*numlink,'siesta-qmmm')  
            fa_link=0.0_dp
            allocate(num_linkqm(numlink))
            num_linkqm=0
            call memory('A','I',numlink,'siesta-qmmm') 
            allocate(num_linkmm(numlink))
            num_linkmm=0
            call memory('A','I',numlink,'siesta-qmmm') 
            allocate(linkqm(numlink,4))
            call memory('A','I',4*numlink,'siesta-qmmm')
            linkqm=0
            allocate(linkmm(numlink,4))
            call memory('A','I',4*numlink,'siesta-qmmm')
            linkmm=0
            allocate(linkmm2(numlink,4,3))
            call memory('A','I',12*numlink,'siesta-qmmm')
            linkmm2=0
            allocate(parametro(num_resonances,numlink,25,4))
            call memory('A','I',num_resonances*100*numlink,
     .                  'siesta-qmmm')
            allocate(link_bondeq(num_resonances,numlink))
            call memory('A','D',num_resonances*numlink,'siesta-qmmm')
            allocate(link_kbond(num_resonances,numlink))
            call memory('A','D',num_resonances*numlink,'siesta-qmmm')
            allocate(link_kangle(num_resonances,numlink,3))
            call memory('A','D',num_resonances*numlink*3,'siesta-qmmm')
            allocate(pclinkmm(numlink,4))
            call memory('A','D',numlink**2,'siesta-qmmm')
            allocate(Emlink(numlink,4))
            call memory('A','D',4*numlink,'siesta-qmmm')
            allocate(distl(numlink))
            call memory('A','D',numlink,'siesta-qmmm')
            allocate(kch(numlink))
            call memory ('A','D',numlink,'siesta-qmmm')
            allocate(rch(numlink))
            call memory ('A','D',numlink,'siesta-qmmm')
            allocate(linkqmtype(numlink,4))
            call memory ('A','I',4*numlink,'siesta-qmmm')
            
C     Sets parameters for all LinkAtoms
            call link1(ng1,na_qm,na_mm,attype,qmattype,nparm,nbond,
     .           bondtype,bondeq,rqm,rmm,nspec,debug,ucell,
     .           lattice_type)
           call get_link_ff_parameters(na_mm,nparm,nbond,nangle,ndihe,
     .      nimp,bondtype,angletype,dihetype,
     .      linkqmtype,attype,imptype,perdihe,kbond,bondeq,kangle,kdihe)
c     Sets to zero Em for HL and CQM
            do i=1,numlink
               Em(linkqm(i,1:num_linkqm(i)))=0.0
            enddo   

         endif

      endif                     !qm & mm

C Allocation of solvent+solute variables
      natot = na_mm + na_qm

      allocate(izs(natot))
      call memory('A','D',natot,'siesta-qmmm')
      allocate(rclas(3,natot))
      call memory('A','D',3*natot,'siesta-qmmm')

      if (qm) then
         do i=1,na_qm
            rclas(1:3,i)=rqm(1:3,i)
         enddo
         deallocate(rqm)
      endif

      if (na_mm>0) then
         do i=1,na_mm
            rclas(1:3,na_qm+i)=rmm(1:3,i)
         enddo
         deallocate(rmm)
      endif

      allocate(masst(natot))
      call memory('A','D',natot,'siesta-qmmm')
      allocate(blocklist(natot))
      call memory('A','D',natot,'siesta-qmmm')

C Some definitions 
      Elj=0.0
      Elink=0.0
      Etots=0.0
      Etot=0.0
      Eelec=0.0
      nfree=natot
      blocklist = 0
      step = 0
      nstepconstr = 0
      constropt = .false.
      constropt = fdf_block_qmmm('ConstrainedOpt',iunit)
      free = .false.
      free = fdf_block_qmmm('FreeEnergy',iunit)
      foundvat = .false.
      writeipl = fdf_boolean('WriteIniParLas',.false.)

C     Read simulation data 
      call read_siesta_qmmm_options( na_qm, natot)

C     assignation of masses and species 
      call assign(na_qm,na_mm,atname,iza_qm,izs,masst)

      if(qm.and.center_qm_system) then
C     Center system  to make the QM system to coincide with SIESTA QM grid.
         call centermol(na_qm,rclas,ucell,natot)

      endif                     !qm & mm

      call reinserting_atoms_in_box(lattice_type,natot,na_qm,na_mm,
     .     aaname,rclas,ucell)

C     Sets LinkAtoms positions
      if(linkatom) then         
         call link3(rclas,na_qm,natot,nparm,bondeq,ucell,
     .        lattice_type)
C     Sets to zero pc and Em for MMLink 
         do i=1,numlink
            pclinkmm(i,1:num_linkmm(i))=pc(linkmm(i,1:num_linkmm(i)))
            pc(linkmm(i,1:1))=0.0
            pc(linkmm(i,2:num_linkmm(i)))=
     .           pc(linkmm(i,2:num_linkmm(i)))+
     .           pclinkmm(i,1)/real(num_linkmm(i)-1,kind=dp)
C     two options: Em(linkmm(i,1:1))=0.0 or Em(linkmm(i,1:4))=0.0
            Emlink(i,1:num_linkmm(i))=
     .           Em(na_qm+linkmm(i,1:num_linkmm(i)))
            Em(na_qm+linkmm(i,1:num_linkmm(i)))=0.0
         enddo
      endif                     !LA

      coulombtype = fdf_string('Coulombtype','cut-off')

C     Calculate Ewald's coefficients
!     if (coulombtype.eq.'ewald') 
      call ewald_coeff(rcortemm,ewald_alpha,kewald_cutoff)

      if (.not.qm.and.mm) then

         call write_input_debug_mm(siesta_qmmmslabel,na_mm,rclas,
     .                 ucell)

         call fdf_shutdown  ! Closes the 'head' file

      else

        call write_input_debug(siesta_qmmmslabel,na_mm,na_qm,nspec,
     .           isa_qm,rclas,ucell)

         call fdf_shutdown  ! Closes the 'head' file

      endif

      deallocate(rclas)

C     the following call assumes that the subroutine reinit is going to
C     read data from INPUT_DEBUG rather than from STANDARD INPUT.
      call siesta_init()

C     Now the file 'INPUT_DEBUG' must be removed to prevent the SIESTA program
C     to read from this file.
      call system('rm INPUT_DEBUG')

      if (natot.ne.na_u-numlink) then
         write(*,*)'natot= ',natot
         write(*,*)'numlink= ',numlink
         write(*,*)'na_u= ',na_u
         call die('natot/=na_u-numlink')
      endif

      natot=na_u
      na_u=na_qm+na_mm
c$$$      if (fc_idyn.and.ia2>na_u) then
c$$$         ia2=na_u
c$$$         fincoor=(ia2-ia1+1)*3*2
c$$$      endif 
      usesavexv=siesta_qmmm_usesavexv
      va=0.0d0
      amass(1:na_u)=masst(1:na_u)
      deallocate(masst)

      foundxv=.false.
      foundvat=.false.
C     Read cell shape and atomic positions from a former run
      if (siesta_qmmm_usesavexv) then
         call qmmm_ioxv( siesta_qmmmslabel, 'read', restart_istep, na_u, 
     .                    ucell, vcell, xa, va, foundxv, foundvat) 
         lattice_type=get_lattice_type(ucell)
      endif

C     Read a crd file from a former run
      call readcrd(na_qm,na_mm,amass,xa,va,foundcrd,foundvat)

      if(qm.and.mm) then
C     Calculate Rcut & block list QM-MM 
         call qmmm_lst_blk(na_qm,na_mm,nroaa,atxres,xa,
     .        rcorteqmmm,radbloqmmm,blockqmmm,listqmmm,rcorteqm,
     .        siesta_qmmmslabel)
      endif                     !qm & mm

C     Read fixed atom constraints
      call qmmm_fixed1(na_qm,na_mm,na_u,nroaa,xa,blocklist,
     .     atname,aaname,aanum,water)

      if ( (.not.opt_idyn ) .and. ( .not. fc_idyn ) .and.
     .                       (.not. foundvat)) then
         
! Build initial velocities according to Maxwell-Bolzmann distribution....
         call vmb(na_u,tempinit,amass,xa,isa,va)

C     Impose constraints to atomic movements by changing forces
         call qmmm_fixed2(na_qm,na_mm,na_u,nfree,blocklist,
     .     blockqmmm,fa,cfa,va)

         call add_siesta_qmmm_ntcon(ucell,na_qm,na_mm,isa,amass,xa,
     .     blocklist,blockqmmm,ntcon)

         if(qm.and.mm) then
C  Correct initial velocities for blocked MM atoms
C  convert F/m in (Ry/Bohr)/amu  to  Bohr/fs**2
           fovermp = 0.009579038d0 *Ang**2 / eV
C     Calculate kinetic energy and temperature ...................
C     Kinetic energy of atoms
            ncount=0
            kin = 0.0d0
            do ia = 1,na_u
               if (sqrt(va(1,ia)**2+va(2,ia)**2+va(3,ia)**2)==0.0_dp) 
     .                    cycle
               ncount=ncount+1
               do i = 1,3
                  kin = kin + 0.5d0 * amass(ia) * va(i,ia)**2 / fovermp
               enddo
            enddo
            tempe = 2.0d0*kin/(3.0d0*ncount-3.0d0)/8.617d-5/eV
            do i = 1,na_u
               if (sqrt(va(1,i)**2+va(2,i)**2+va(3,i)**2)==0.0_dp)
     .          cycle
               do ix = 1,3
                  va(ix,i) = va(ix,i) * dsqrt(tempinit/tempe)
               enddo
            enddo
         endif
C     Calculate kinetic energy and temperature ...................
C     Kinetic energy of atoms
            ncount=0
            kin = 0.0d0
            do ia = 1,na_u
               if (sqrt(va(1,ia)**2+va(2,ia)**2+va(3,ia)**2)==0.0_dp) 
     .                    cycle
               ncount=ncount+1
               do i = 1,3
                  kin = kin + 0.5d0 * amass(ia) * va(i,ia)**2 / fovermp
               enddo
            enddo
            tempe = 2.0d0*kin/(3.0d0*ncount-3.0d0)/8.617d-5/eV
            do i = 1,na_u
               if (sqrt(va(1,i)**2+va(2,i)**2+va(3,i)**2)==0.0_dp)
     .          cycle
               do ix = 1,3
                  va(ix,i) = va(ix,i) * dsqrt(tempinit/tempe)
               enddo
            enddo

         if (3*ncount/=3*natot-ntcon) 
     .              call die('3*ncount/=3*natot-ntcon')

      endif

C     Dump initial coordinates after qmmm_ioxv for debuging in PDB format
      if (debug) then       
         write(6,'(/a)') 'Dump of initial coordinates for debuging: '
         if(qm) then
            write(6,'(A4,I7,2x,A2,2x,A4,A,I4,4x,3f8.3)')
     .           ('ATOM',i,atsym(isa_qm(i)),'STO ','A',i,
     .           (xa(ia,i)/Ang,ia=1,3), i=1,na_qm)
         endif                  !qm
         if(mm) then
            write(6,'(A4,I7,2x,A4,A4,A,I4,4x,3f8.3)')
     .           ('ATOM',na_qm+i,atname(i),aaname(i),'A',aanum(i),
     .           (xa(ia,na_qm+i)/Ang,ia=1,3), i=1,na_mm)
         endif                  !mm
         write(6,'(A3)') 'END'
      endif                     !debug

C     Sets LinkAtoms positions for new positions
      if(linkatom.and.(foundxv.or.foundcrd)) then
            call link3(xa,na_qm,natot,nparm,bondeq,ucell,
     .                     lattice_type)
      endif                     !LA

C     Dump initial coordinates to output 
      if (writeipl) then
         write(6,*) 
         write(6,"('siesta-qmmm: ',71(1h=))")
         if(qm) then
            write(6,'(a)') 'siesta-qmmm: Atomic coordinates (Ang) and '
     .                                   //'species'
            write(6,"(i6,2x,3f10.5,2x,i3)")
     .          (ia, (xa(ix,ia)/Ang, ix=1,3), isa_qm(ia), ia=1,na_qm)
         endif                  !qm

C     Dump initial Solvent coordinates to output
         if(mm) then      
            write(6,*)
            write(6,"('siesta-qmmm: Number of solvent atoms:',i6)") 
     .                  na_mm
            write(6,"('siesta-qmmm: Number of residues:',i4)") nroaa
            write(6,'(/a)')'siesta-qmmm: Solvent coordinates (Ang)'
            write(6,"(i6,2x,A4,2x,A4,2x,i5,2x,3f12.6)")
     .           (i,atname(i),aaname(i),aanum(i),
     .           (xa(ix,i+na_qm)/Ang, ix=1,3),i=1,na_mm)
            write(6,'(/a)') 'siesta-qmmm: Atoms parameters'
         endif                  !mm
         if(qm.and.mm) then
            write(6,"(i6,2x,A4,2x,f12.7,2x,f12.7)")
     .           (i,qmattype(i),Em(i),Rm(i),i=1,na_qm)
         endif                  !qm & mm
         if(mm) then
            write(6,"(i6,2x,A4,2x,f9.6,2x,f9.6,2x,f9.6)")
     .           (i,attype(i),Em(i+na_qm),Rm(i+na_qm),pc(i),i=1,na_mm)
         endif                  !mm
         write(6,"('siesta-qmmm: ',71(1h=))")
      endif                     !writeipl
      call pxfflush(6)

      if (qm.and.mm.and..not.null_mm_charges) then

         filerho = paste( siestaslabel, '.TOCH' )
         filevqm = paste( siestaslabel, '.VEXT' ) 

         nsm    = fdf_integer('MeshSubDivisions',2)
         nsm    = max(nsm,1)

C     Set mesh sub-division variables & perform one off allocation
         nsp = nsm*nsm*nsm
         nullify( xdsp )
         call re_alloc( xdsp, 1, 3, 1, nsp, name='xdsp',
     &        routine='dhscf' )

         allocate(atm_rmaxo(na_qm + numlink))

C     ----------------------------------------------------------------------
C     Mesh initialization 
C     ----------------------------------------------------------------------
         G2max = g2cut
#ifdef QMMM_BSC
         call InitMesh( na_s, ucell, no_s, iaorb, iphorb,         
     .        isa, rmaxo, G2max, G2mesh, nsc, nmpl, nm,
     .        nml, ntm, ntml, ntpl, dvol)
#else
         call InitMesh( na_s, ucell, no_s, iaorb, iphorb,         
     .        isa, rmaxo, G2max, G2mesh, nsc, nmpl, nm,
     .        nml, ntm, ntml, ntpl, ntopl, dvol)
#endif

C     calculate cell volume
         volume = volcel( ucell )
         dvol=volume/ntpl

         nullify(Vqm,Rho,resonance_Rho,auxrho,Chlocal)
C Allocate memory related to ntpl
c     allocate Vqm: the external potential (Vext in siesta)
         call re_alloc( Vqm, 1, ntpl, name='Vqm', 
     .            routine='siesta-qmmm' ) 
c     allocate Rho: total (ionic+electronic) charge
         call re_alloc( Rho, 1, ntpl, name='Rho', 
     .             routine='siesta-qmmm' )
         call re_alloc( resonance_Rho, 1, ntpl, name='Vqm', 
     .             routine='siesta-qmmm' )
         call re_alloc( auxrho, 1, ntpl,name='auxrho',
     .                  routine='siesta-qmmm')
         call re_alloc( Chlocal, 1, ntpl,name='Chlocal',
     .          routine='siesta-qmmm')

         do ia = 1,na_qm
            atm_rmaxo(ia)=0.0_dp
            is = isa_qm(ia)
            noa  = nofis(is)
            do io = 1,noa
               atm_rmaxo(ia) = max( atm_rmaxo(ia), rcut(is,io) )
            enddo
         enddo

         if (numlink>0) then
            do resonance_id=1,num_resonances
               do ia = 1,numlink
                  atm_rmaxo(na_qm+ia)=0.0_dp
                  is = resonance(1)%isa(ia)
                  noa  = nofis(is)
                  do io = 1,noa
                     atm_rmaxo(ia) = max( atm_rmaxo(ia), rcut(is,io) )
                  enddo
               enddo
            enddo
         endif

         Rho=0.0
         auxrho=0.0
         reorder_rho=.true.
         if (usesavedm) then
            do resonance_id=1,num_resonances
               ntplext = ntpl
               nspinext = 1 
               call check_rho( trim(resonance(resonance_id)%path)//
     .              filerho, ntplext, nspinext, qmmm_nsm,  
     .              found, overflow )
               if (.not. found) then
                  if (harrisdens) then
C     Initialise quantities relating to the atom-mesh positioning
                    call InitAtomMesh( na_s, xa)
C     Find diffuse ionic charge on mesh
                    call LocalChargeOnMesh( na_s, isa, ntpl, Chlocal,  
     .                          indxua, .false., dvol, 0.0, 0.0_dp )
C     Calculate orbital values on mesh
                    call PhiOnMesh( nmpl, no_s, iaorb, iphorb, isa )
C     Find Harris electron density (sum of atomic densities)
C     for every point ,so that Vqm is only calculated 
C     in points where Rho /= 0.
                    call rhooda( no_s, nmpl, datm, resonance_Rho, iaorb,
     .                               iphorb, isa )
                    auxrho=auxrho+Chlocal+
     .                resonance(resonance_id)%weight*resonance_Rho
                  else
                    call fill_auxrho(na_qm + numlink,ucell,xa,ntm,
     .                    ntpl,atm_rmaxo,auxrho) 
                  endif
               else
                  if (overflow) then
                     if (Node .eq. 0) then
                        write(6,100)
                        write(6,101)
                        write(6,'(a)')
     .                       'siesta-qmmm:   Format or size of the '//
     .                       'external potential file'
                        write(6,'(a)')
     .                       'siesta-qmmm:   does not match with '//
     .                       'this system'
                        write(6,102)
                     endif
                     call die("check_rho: Unable to read  "//
     .                    "TotCharge from file")
                  endif
                  write(6,*)'Reading previous Rho from directory: ',
     .                            trim(resonance(resonance_id)%path)

c     allocate Rho: total (ionic+electronic) charge
                  call read_rho( trim(resonance(resonance_id)%path)//
     .                 filerho, ucell, ntm, qmmm_nsm, ntpl, 1, 
     .                 resonance_Rho )
                  reorder_rho=.false.
                  auxrho=auxrho+
     .              resonance(resonance_id)%weight*resonance_Rho 
               endif
            enddo
         else
            do resonance_id=1,num_resonances
               if (harrisdens) then
C     Initialise quantities relating to the atom-mesh positioning
                 call InitAtomMesh( na_s, xa)
C     Find diffuse ionic charge on mesh
                 call LocalChargeOnMesh( na_s, isa, ntpl, Chlocal,  
     .                          indxua, .false., dvol, 0.0, 0.0_dp )
C     Calculate orbital values on mesh
                 call PhiOnMesh( nmpl, no_s, iaorb, iphorb, isa )
C     Find Harris electron density (sum of atomic densities)
C     for every point ,so that Vqm is only calculated 
C     in points where Rho /= 0.
                 call rhooda( no_s, nmpl, datm, resonance_Rho, iaorb,  
     .                               iphorb, isa )
                 auxrho=auxrho+Chlocal+
     .              resonance(resonance_id)%weight*resonance_Rho
               else
                 call fill_auxrho(na_qm + numlink,ucell,xa,ntm,
     .                    ntpl,atm_rmaxo,auxrho) 
               endif
            enddo
         endif
         if (reorder_rho) call reord( auxrho, auxrho, nml, nsm, +1)
         Rho=auxrho

      endif  ! qm & mm & .not.null_mm_charges

      if (qm) then
!     Launch siesta process
         call system('rm -f directories')
         do i=1,num_resonances
            call system('echo '//trim(resonance(i)%path)//
     .                          ' >> directories') 
         enddo

C     Write the 'siesta_qmmmslabel.siesta.fdf' file to be read by siesta as subroutine
         call write_lab( siesta_qmmmslabel, siestaslabel, mm, na_qm,   
     .                         isa_qm, xa, null_mm_charges )

!     Set physical units for communication with siesta
         call siesta_units( 'Bohr', 'Ry' )

         call create_fifos( siestaslabel)

C     writes cell
         write(6,'(/,a,3(/,a,3f12.6))')
     .        'siesta-qmmm: unit cell vectors (Ang) from siesta run:',
     .        ('    ' , (ucell(ix,ia)/Ang,ix=1,3), ia =1,3)

      endif                     !qm.and.mm

      call qmmm_timer( 'Setup', 2 )
      call printmemory( 6, 0 )

C     Start loop over constrained optimization steps
      if(constropt) then
         call subconstr1(nconstr,typeconstr,kforce,nstepconstr,
     .        rini,rfin,atmsconstr,dr,ro,ndists,coef,constropt)
      endif
      allocate(qmstep(nstepconstr+1))
      call memory('A','D',nstepconstr+1,'siesta-qmmm')
      qmstep=0
      Etot=0.0_dp
      Etots=0.0_dp
      fa=0.0_dp
      stress=0.0_dp
C      cstress=0.0_dp
      lattice_type=get_lattice_type(ucell)

      do istepconstr=1,nstepconstr+1
         if(constropt) then
            write(6,*)
            write(6,'(A)')    '*******************************'
            write(6,'(A,i5)') '  Constrained Step : ', istepconstr
            write(6,'(A)')    '*******************************'
         endif

         totcoor = fincoor - inicoor + 1

C     Start loop over coordinate changes
         istp = 0
         if (restart_istep>0) inicoor=restart_istep+1

         if (inicoor>fincoor) then
            relaxd=.true.
            goto 10 
         endif

         do istep = inicoor,fincoor

            if (state_analaysis_flag) then
            ! Initialization tasks for a given geometry
               call state_init( istep )
            else
               cml_p=.false.
               call timer( 'IterMD', 1 )
            endif

            call qmmm_timer( 'IterMD', 1 )      
            istp = istp + 1

            write(6,'(2/,7x,A)')    '********************************'
            write(6,'(9x,a,2x,i6)')'siesta-qmmm:  Begin Step= ',istep
            write(6,'(7x,A,2/)')    '********************************'

c     actualize the Neb list every mneigh_freq steps
            actualiz=.false.
            if(mod(istep,mneigh_freq).eq.0.or.istep.eq.inicoor) 
     .                       actualiz=.true.

C reinsert atoms in the box to have better images unless the Nose thermostat is used since
C the self-consistency in the Nose algorithm is not compatible with the reinsertion of atoms 
            if (actualiz.and.(.not. therm_idyn)) then
               call reinserting_atoms_in_box(lattice_type,na_u,na_qm,
     .              na_mm,aaname,xa,ucell)
            endif

            if (qm.and.mm.and.actualiz)
     .             call actualize_qmmm_mneighb(na_qm,na_mm,na_u,atname,
     .                 aaname,xa,ucell,pc,rcorteqmmm,
     .                 rmaxo,actualiz,lattice_type,include_qmqm_vdw)

C     Calculate the external potential Vqm
            if(qm.and.mm.and..not.null_mm_charges) then

               if (harrisdens) then
                 auxrho=Rho
               else
                 auxrho=0.0
                 call fill_auxrho(na_qm + numlink,ucell,xa,ntm,
     .                    ntpl,atm_rmaxo,auxrho) 
               endif

C     Select Ewald method or do direct summation limited to rcorteqm
               if (coulombtype.eq.'ewald') then
                  call pcpot_ewald(na_qm,na_mm,na_u,ntpl,
     .                 ntm,listqmmm,xa,ucell,pc,
     .                 rcorteqm,ewald_alpha,
     .                 kewald_cutoff,auxrho,Vqm,lattice_type)
               else if (coulombtype.eq.'cut-off') then
                  call pcpot(na_qm,na_mm,na_u,ntpl,ntm,
     .                 listqmmm,xa,ucell,pc,rcorteqm,
     .                 auxrho,Vqm,lattice_type)
               else
                  write(*,*) 'coulombtype= ',coulombtype
                  stop 'Unrecognised coulombtype'
               endif

C     Write Vqm to file
C     call reord( Vqm, Vqm, nml, qmmm_nsm, +1)

               call write_rho( filevqm, ucell, ntm, qmmm_nsm, ntpl, 
     $              1, Vqm)
               if (num_resonances>1) then
                  do resonance_id=2,num_resonances
                     call system('cp '//filevqm//' '//
     $                    trim(resonance(resonance_id)%path)//filevqm)
                  enddo
               endif
               call write_grid_netcdf( ucell, ntm, 1, ntpl,
     $              Vqm, "ExternalPotential")
C     call reord( Vqm, Vqm, nml, qmmm_nsm, -1)

            endif                !qm & mm & 

!     Calculate Energy and Forces using Siesta as Subroutine
            if(qm) then
               call qmmm_timer( 'Siesta', 1 )
               call siesta_send_coordinates( siestaslabel, na_qm, xa, 
     $                    ucell )
               call siesta_receive_forces( siestaslabel, na_qm, Etots,
     $                                        fa, stress )
               call qmmm_timer( 'Siesta', 2 )   
               qmstep(istepconstr) = qmstep(istepconstr) + 1
            endif               !qm

C     Read Rho from file 
            if(qm.and.mm.and..not.null_mm_charges) then
               Rho=0.0_dp
               do resonance_id=1,num_resonances
                  ntplext = ntpl
                  nspinext = 1 
                  call check_rho( trim(resonance(resonance_id)%path)//
     .                 filerho, ntplext, nspinext, qmmm_nsm,  
     .                 found, overflow )
                  if (.not. found) then
                     if (Node .eq. 0) then
                        call die('siesta-qmmm: WARNING: total '//
     .                       'charge file not found!')
                     endif
                  else
                     if (overflow) then
                        if (Node .eq. 0) then
                           write(6,100)
                           write(6,101)
                           write(6,'(a)')
     .                          'siesta-qmmm:   Format or size of the '
     .                          //'external potential file'
                           write(6,'(a)')
     .                          'siesta-qmmm:   does not match with '//
     .                          'this system'
                           write(6,102)
                        endif
                        call die("check_rho: Unable to read  "//
     .                       "TotCharge from file")
                     endif
c     allocate Rho: total (ionic+electronic) charge
                     call read_rho( trim(resonance(resonance_id)%path)//
     .                    filerho, ucell, ntm, qmmm_nsm, ntpl, 1, 
     .                            resonance_Rho )
                     Rho=Rho+
     .                     resonance(resonance_id)%weight*resonance_Rho
                  endif
               enddo
C               call reord( Rho, Rho, nml, qmmm_nsm, -1)
            endif               !qm & mm & .not.null_mm_charges

C     Start MMxQM loop
            do imm=1,mmsteps

               step = step +1
               if(mmsteps.ne.1) then
                  write(6,*)
                  write(6,'(A)')    '*******************************'
                  write(6,'(A,i5)') '   MM x QM Step : ', imm 
                  write(6,'(A)')    '*******************************'
               endif

               if(qm.and.mm)call qmmm_timer( 'QM-MMcoupl', 1 )
C     Calculation of Forces over the classical atoms from QM-MM interaction
               if(qm.and.mm.and..not.null_mm_charges) then
                  if (coulombtype.eq.'ewald') then
                     call mmforce_ewald(na_qm,na_mm,na_u,nspec,ntpl,
     .                    ntm,xa,fa,stress,pc,ucell,rcorteqmmm,
     .                    rcorteqm,ewald_alpha,kewald_cutoff,dvol,Rho,
     .                    lattice_type)
                  else
                     call mmforce(na_qm,na_mm,na_u,nspec,ntpl,ntm,
     .                    xa,fa,stress,pc,ucell,rcorteqmmm,rcorteqm,
     .                    dvol,Rho,lattice_type)
                  endif
               endif            !qm & mm & .not.null_mm_charges

C     Calculation of last QM-MM interaction: LJ Energy and Forces only 
               if(qm.and.mm) then
                  call ljef(na_qm,na_mm,na_u,xa,Em,Rm,fa,stress,Elj,
     .                           rcorteqmmm,ucell,lattice_type)
                  call qmmm_timer( 'QM-MMcoupl', 2 )
               endif            !qm & mm

C     LinkAtom: set again linkmm atoms parameters
               if(qm.and.mm) then
                  if(linkatom) then
                     do i=1,numlink
                        pc(linkmm(i,1:num_linkmm(i)))=
     .                           pclinkmm(i,1:num_linkmm(i))
                        Em(na_qm+linkmm(i,1:num_linkmm(i)))=
     .                                    Emlink(i,1:num_linkmm(i))
                     enddo
                  endif         !LA
               endif            !qm & mm

C     Calculate pure Solvent energy and forces 
               if(mm) then
                  call qmmm_timer( 'MMenergy', 1 )
                  call solv_ene_fce(na_u,na_qm,na_mm,ng1,xa,
     .                 Em,Rm,pc,Etot_amber,fce_amber,stress_amber,
     .                 attype,nbond,nangle,ndihe,nimp,multidihe, 
     .                 multiimp,kbond,bondeq,
     .                 kangle,angleeq,kdihe,diheeq,kimp,impeq,
     .                 perdihe,perimp,bondtype,angletype,dihetype,
     .                 imptype,bondxat,angexat,angmxat,dihexat,
     .                 dihmxat,impxat,atange,atangm,atdihe,atdihm,
     .                 atimp,ng1type,angetype,angmtype,evaldihe,
     .                 evaldihm,dihety,dihmty,impty,nonbonded,
     .                 scale,nonbondedxat,scalexat,
     .                 evaldihelog,evaldihmlog,step,nparm,
     .                 graphite_layer_no,
     .                 actualiz,rcortemm,
     .                 atname,aaname,sfc,dt,
     .                 water,amass,ucell,lattice_type,coulombtype,
     .                 ewald_alpha,kewald_cutoff)
                  call qmmm_timer( 'MMenergy', 2 )
               endif            !mm

C     converts fa to Kcal/mol/Ang  
               fa(1:3,1:na_u)=fa(1:3,1:na_u)*Ang/eV*kcal
               if (numlink>0) then
                  fa_link(1:num_resonances,1:3,1:numlink)=
     .               fa_link(1:num_resonances,1:3,1:numlink)*Ang/eV*kcal
               endif
               stress(1:3,1:3)=stress(1:3,1:3)*Ang*Ang*Ang/eV*kcal

C     add famber to fa  
               if(mm) then
                  fa(1:3,na_qm+1:na_u)=fa(1:3,na_qm+1:na_u)
     .                 +fce_amber(1:3,1:na_mm)
                 stress(1:3,1:3)=stress(1:3,1:3)+stress_amber(1:3,1:3)
               endif            !mm

C     Calculation of LinkAtom Energy and Forces
               if(qm.and.mm) then
                  if(linkatom) then
                     call link2(xa,na_u,na_mm,fa,stress,ng1,nparm,
     .                    ndihe,multidihe,multiimp,kbond,bondeq,kangle,
     .                    angleeq,kdihe,diheeq,kimp,impeq,perdihe,
     .                    perimp,bondxat,Elink,step,
     .                    kch,rch,lattice_type,ucell)
C     Set again link atmos parameters to zero for next step  
                     do i=1,numlink
                        pclinkmm(i,1:num_linkmm(i))=
     .                               pc(linkmm(i,1:num_linkmm(i)))
                        pc(linkmm(i,1:num_linkmm(i)))=0.0
                        pc(linkmm(i,2:num_linkmm(i)))=
     .                        pc(linkmm(i,2:num_linkmm(i)))+
     .                      pclinkmm(i,1)/real(num_linkmm(i)-1,kind=dp)
                        Em(na_qm+linkmm(i,1:num_linkmm(i)))=0.0
                     enddo
                  endif         ! LA
               endif            !qm & mm

C     Calculation of Constrained Optimization Energy and Forces 
               if(imm.eq.1) then
                  if(constropt) then
                     call subconstr2(nconstr,typeconstr,kforce,rini,
     .                    rfin,ro,rt,nstepconstr,atmsconstr,na_u,
     .                    xa,fa,istp,istepconstr,
     .                    ndists,coef,ucell,lattice_type)
                  endif 
               endif            !imm

C     Free Energy calculation using Jarzynski formula 
               if(imm.eq.1) then
                  if(free) then
                     call subfree(na_u,xa,fa,dt,free,totcoor,
     .                    ucell,lattice_type)
                  endif
               endif            !imm

C     call to subroutine that calculates the Electric field energy
               if(mm) then
                  call elecfield(na_mm,xa(1:3,na_qm+1:na_u),
     .                 fa(1:3,na_qm+1:na_u),pc,Eelec)
               endif            !mm

C     Converts fa to Ry/Bohr since the SIESTA subroutines for siesta_move has got
C      targeted parameters in SIESTA internal units (Ry, Bohr, etc..)
               fa(1:3,1:na_u)=fa(1:3,1:na_u)/Ang*eV/kcal
               stress(1:3,1:3)=stress(1:3,1:3)/(Ang*Ang*Ang)*eV/kcal

C     Writes final energy decomposition
               Etot=Etots+Elj+((Etot_amber+Elink+Eelec)/kcal*eV)
               write(6,*)
               write(6,'(/,a)') 
     .                  'siesta-qmmm: Energy Decomposition (eV):'
               if(qm) write(6,'(a,2x,F16.6)')           'Esiesta:',
     .              Etots/eV
               if(qm.and.mm) write(6,'(a,2x,F16.6)')    'Elj:    ',
     .              Elj/eV
               if(mm) write(6,'(a,2x,F16.6)')      'Esolv:  ',
     .              Etot_amber/kcal
               if(Elink.ne.0.0) write(6,'(a,2x,F16.6)') 'Elink:  ',
     .              Elink/kcal
               if(Eelec.ne.0.0) write(6,'(a,2x,F16.6)') 'Eefield:',
     .              Eelec/kcal
               if(qm.or.mm) write(6,'(a,2x,F16.6)')    'Etot:  ',
     .              Etot/eV
               call pxfflush(6)

C     Sets fa to zero inside mmxqm step
               if(qm.and.mm) then
                  if(imm.ne.1) then
                     fa(1:3,1:na_qm) = 0.0
                     va(1:3,1:na_qm)=0.0
                     if(linkatom) then
                        do i=1,numlink
                           fa(1:3,linkmm(i,1)+na_qm)=0.0
                           va(1:3,linkmm(i,1)+na_qm)=0.0
                        enddo
                     endif      !LA
                  endif         !imm
               endif            !qm & mm

C
      if (state_analaysis_flag) then


C     Impose constraints to atomic movements by changing forces
        call qmmm_fixed2(na_qm,na_mm,na_u,nfree,blocklist,
     .     blockqmmm,fa,cfa,va)
        call state_analysis(istep)
        call add_siesta_qmmm_ntcon(ucell,na_qm,na_mm,isa,amass,xa,
     .     blocklist,blockqmmm,ntcon)

      else

      !  Compute stress without internal molecular pressure
      call remove_intramol_pressure(ucell,stress,na_u,xa,fa,mstress)

C Impose constraints to atomic movements by changing forces ...........
      if (RemoveIntraMolecularPressure) then
         ! Consider intramolecular pressure-removal as another
         ! kind of constraint
         call fixed(ucell,mstress,na_qm,isa, amass, xa, fa,
     $        cstress, cfa, ntcon )
      else
         call fixed(ucell,stress,na_qm,isa, amass, xa, fa,
     $        cstress, cfa, ntcon )
      endif

C     Impose constraints to atomic movements by changing forces
      call qmmm_fixed2(na_qm,na_mm,na_u,nfree,blocklist,
     .     blockqmmm,fa,cfa,va)
      call add_siesta_qmmm_ntcon(ucell,na_qm,na_mm,isa,amass,xa,
     .     blocklist,blockqmmm,ntcon)

! Calculate and output Zmatrix forces
      if (lUseZmatrix .and. opt_idyn) then
        call CartesianForce_to_ZmatForce(na_u,xa,fa)
        if (IOnode) call iofaZmat()
      endif

C     Compute kinetic contribution to stress
               kin_stress(1:3,1:3) = 0.0_dp
               volume = volcel(ucell)
               volume_of_some_cell=volume
               do ia = 1,na_u
                  do jx = 1,3
                     do ix = 1,3
                        kin_stress(ix,jx) = kin_stress(ix,jx) -
     .                       amu * amass(ia) * va(ix,ia) * va(jx,ia) 
     .                                     / volume
                     enddo
                  enddo
               enddo
C     Add kinetic term to stress tensor
               tstress = stress + kin_stress

! Force output 
      if (IOnode) then
        call siesta_write_forces()
        call siesta_write_stress_pressure()
      endif

      endif

C     Accumulate coordinates in PDB/CRD file for animation
      call wripdb(na_qm,siesta_qmmmslabel,xa,na_u,step,wricoord,
     .     na_mm,atname,aaname,aanum,nspec,atsym,isa_qm,
     .     listqmmm,blockqmmm,ucell)

C      if(nfce.ne.na_u) call iofa(na_u,cfa)

C     Move atoms 
               call siesta_move( istep, relaxd )

               if((.not.opt_idyn) .and. (.not.fc_idyn)) then
                  write(6,'(/,a,f12.6)')
     .                 'siesta-qmmm: Kinetic Energy (eV):',Ekinion/eV
                  write(6,'(/,a,f12.6)')
     .                 'siesta-qmmm: Total Energy + Kinetic (eV):',(Etot
     .                         +Ekinion)/eV
                  write(6,'(/,a,f12.3,a)')
     .                 'siesta-qmmm: System Temperature:', tempion, ' K'
C     if(qm) call centerdyn(na_qm,xa,ucell,na_u)
               endif

C     Write Energy in file
               call wriene(istep,dt,siesta_qmmmslabel,idyn,Etot,cfa,
     .                        Ekinion,tempion)

C     Save last atomic positions and velocities
               call qmmm_ioxv( siesta_qmmmslabel, 'write', istep, na_u,  
     .                     ucell, vcell, xa, va, foundxv, foundva )

C     write atomic constraints each step
               call wrtcrd(na_u,xa,ucell,lattice_type)              

c     sets variables for next siesta cycle
               fa = 0.0_dp
               stress = 0.0_dp
               lattice_type=get_lattice_type(ucell)

C     Calculation Hlink's New positions 
               if(qm.and.mm) then
                  if(linkatom) then
                     call link3(xa,na_qm,na_u,nparm,bondeq,ucell,
     .                                     lattice_type)
                  endif         !LA
               endif            !qm & mm

C     Exit MMxQM loop
            enddo               !imm
            if((mmsteps.ne.1).and.(imm.ne.1)) relaxd = .false.

C     Exit coordinate relaxation loop
            call qmmm_timer( 'IterMD', 2 )     

      if (qm.and..not.null_mm_charges) then

C     ----------------------------------------------------------------------
C     Check whether mesh initialisation is needed
C     ----------------------------------------------------------------------
         call SameMeshAndAtoms( na_s, xa, ucell, rmaxo, G2max, G2mesh,
     .        samesh, samexa )

C     ----------------------------------------------------------------------
C     Initialization due to change of mesh
C     ----------------------------------------------------------------------
         if ( .not.samesh ) then

C     ----------------------------------------------------------------------
C     Unit cell handling
C     ----------------------------------------------------------------------
C     Find diagonal unit cell and supercell
            call digcel( ucell, mscell, cell, scell, nsc, IsDiag )
            if (.not.IsDiag) then
               if (Node.eq.0) then
                  write(6,'(/,a,3(/,a,3f12.6,a,i6))')
     .                 'siesta-qmmm: WARNING: New shape of unit cell 
     .                  and supercell:',('siesta-qmmm:',  
     .                  (cell(i,j),i=1,3),'x',nsc(j),j=1,3)
               endif
            endif

C     Find system volume
            volume = volcel( ucell )

            G2max=g2cut
C     ----------------------------------------------------------------------
C     Mesh initialization 
C     ----------------------------------------------------------------------
#ifdef QMMM_BSC
         call InitMesh( na_s, ucell, no_s, iaorb, iphorb,         
     .        isa, rmaxo, G2max, G2mesh, nsc, nmpl, nm,
     .        nml, ntm, ntml, ntpl, dvol)
#else
         call InitMesh( na_s, ucell, no_s, iaorb, iphorb,         
     .        isa, rmaxo, G2max, G2mesh, nsc, nmpl, nm,
     .        nml, ntm, ntml, ntpl, ntopl, dvol)
#endif

C     Allocate memory related to ntpl
c     allocate Vqm: the external potential (Vext in siesta)
            call re_alloc( Vqm, 1, ntpl, name='Vqm',
     .            routine='siesta-qmmm' ) 
c     allocate Rho: total (ionic+electronic) charge
            call re_alloc( Rho, 1, ntpl, name='Rho', 
     .          routine='siesta-qmmm' )
            call re_alloc( resonance_Rho, 1, ntpl, name='Vqm', 
     .           routine='siesta-qmmm' )
            call re_alloc( auxrho, 1, ntpl,name='auxrho',
     .                   routine='siesta-qmmm')
            call re_alloc( Chlocal, 1, ntpl,name='Chlocal',
     .                   routine='siesta-qmmm')

         endif

      endif

         enddo
 10      continue
C     End of coordinate-relaxation loop ==================================

C     End of constrain optimization loop
         if(constropt) then
            call subconstr3(ro(1),rt(1),dr,Etot)
            call wrirtc(siesta_qmmmslabel,Etot,rt(1),istepconstr,na_qm,
     .           na_mm,na_u,xa,atname,aaname,aanum,nspec,atsym,isa_qm)
            call ioxvconstr( na_u, ucell, xa, va, istepconstr )
         endif
      enddo                     !istepconstr

!     Stop siesta process
        if(qm) call siesta_quit( 'all' )

C      call siesta_analysis( relaxd )

C      call siesta_end()

C     Mulliken population analysis
c$$$      if(qm) then
c$$$         if(mullipop.eq.1) then
c$$$            call qmmm_mulliken(na_qm,isa_qm,iza_qm,nspec,atsym,
c$$$     .           siestaslabel,qmstep,nstepconstr+1)
c$$$         endif
c$$$      endif                     !qm    



C     Dump last coordinates to output
      if (writeipl) then
         if(qm) then
            write(6,'(/a)')'siesta-qmmm: Last atomic coordinates '//
     .           '(Ang) and'//' species'
            write(6,"(i6,2x,3f10.5,2x,i3)")
     .          (ia, (xa(ix,ia)/Ang, ix=1,3), isa_qm(ia), ia=1,na_qm)
         endif                  !qm
C     Dump last Solvent coordinates to output
         if(mm) then
            write(6,'(/a)')'siesta-qmmm: Last solvent coordinates (Ang)'
            write(6,'(A4,I7,2x,A4,A4,A,I4,4x,3f8.3)')
     .           ('ATOM',i,atname(i),aaname(i),'A',aanum(i),
     .           (xa(ia,na_qm+i)/Ang,ia=1,3), i=1,na_mm)
         endif                  !mm
      endif                     !writeipl
      call pxfflush(6)

C     Stop time counter
      call qmmm_timer( 'Siesta-qmmm', 2 )
      call qmmm_timer( 'all', 3 )
      call printmemory( 6, 0 )

C     Print final date and time
      call timestamp('End of run')

100   format(/,'siesta-qmmm: ',71(1h*))
101   format('siesta-qmmm:                  INPUT ERROR')
102   format('siesta-qmmm: ',71(1h*))

      end


