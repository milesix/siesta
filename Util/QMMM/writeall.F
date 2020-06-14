c subroutine that writes the energy  

      subroutine wriene(istp,dt,slabel,idyn,Etots,cfa,Ekinion,tempion)  

      use precision, only : dp
      use siesta_geom, only : na_u
      use siesta_qmmm_options, only : opt_idyn, fc_idyn

      implicit          none
      character         slabel*30, paste*30
      integer           istp,idyn
      real(dp)  dt, Etots, cfmax, Ekinion, tempion 
      real(dp) cfa(3,na_u) 
      external          io_assign, io_close, paste 

      character         fname*30 
      logical           frstme
      integer           unit
      real(dp)  eV,Ang
      save              frstme,Ang,eV,unit,fname
      data              frstme /.true./

c -------------------------------------------------------------------

      if ( frstme ) then
        fname = paste(slabel,'.ene')
        Ang  = 1.d0 / 0.529177d0
        eV  = 1.d0 / 13.60580d0
        frstme = .false.
      endif

      cfmax = maxval(abs(cfa))

      call io_assign(unit)
      open( unit, file=fname, form = 'formatted', position='append',
     .      status='unknown')

	if(opt_idyn .or. fc_idyn ) then
        if (istp .eq. 0)
     . write(unit,*)'#step   Energy(eV)   Max. Force (eV/Ang)   '
      write(unit,'(i5,2x,F14.7,2x,F14.7)') istp,Etots/eV,cfmax*Ang/eV
	else
        if (istp .eq. 1)
     . write(unit,*)'# time(fs)      K.E.(eV)     P.E.(ev)      '//
     .         'T.E.(eV)    Temp(k)'
      write(unit,'(F9.2,2x,F12.4,2x,F12.4,2x,F12.4,2x,F7.2)')
     .istp*dt,Ekinion/eV,Etots/eV,(Etots+Ekinion)/eV,tempion
	endif

      call io_close(unit)
 
      return
      end

c*************************************************************************
c subroutine that writes the coordinates in PDB/CRD format

      subroutine wripdb(na,slabel,rclas,natot,istp,wricoord,nac,atname,
     .                  aaname,aanum,nesp,atsym,isa,list,block,mmcell)

      use precision, only : dp
      use linkatoms, only : numlink, rlink, resonance

      implicit          none
      character         slabel*30, paste*30
      integer           na,nac,natot,istp,wricoord
      integer           nesp, isa(na)
      integer           list(nac),block(nac)
      real(dp), dimension(3,3) :: mmcell
      character*2       atsym(nesp)       
      real(dp)  rclas(3,natot)
      external          paste

      character         fnamep*30,fnamec*30,fnamel*30,title*30,ch*4
      character*1       chain(natot)
      logical           frstme,foundc
      integer           unitp,unitc,unitl
      real(dp)  eV,Ang
      save              frstme,Ang,eV,unitc,fnamec,unitl,fnamel
      data              frstme /.true./

c Varibles added to write a pdb file
      integer     i,ia,aanum(nac)
      character*4 aaname(nac),atname(nac)

      real(dp), save        :: a, b, c, alpha, beta, gamma
      real(dp)              :: degtorad

c -------------------------------------------------------------------
      chain(1:natot)='A'
      do i=1,nac
        if(list(i).eq.0)  chain(i+na)='L'
        if(block(i).eq.0) chain(i+na)='B'
        if(list(i).eq.0.and.block(i).eq.0) chain(i+na)='C'
      enddo

      if ( frstme ) then
        Ang  = 1.d0 / 0.529177d0
        eV  = 1.d0 / 13.60580d0
      fnamep = paste(slabel,'.init.pdb')
      fnamec = paste(slabel,'.mdcrd')
      fnamel = paste(slabel,'.last.pdb')

      degtorad = 4.0d0*atan(1.0d0)/180.0d0
      call uncell(mmcell,a,b,c,alpha,beta,gamma,degtorad)

c writes initial PDB coords in .init.pdb file 
      call io_assign(unitp)
      open( unitp, file=fnamep, form = 'formatted', status='unknown')
       write(unitp,'(A,3(F9.3),3(F7.2),A)')'CRYST1',a/Ang,b/Ang,c/Ang,
     .       alpha,beta,gamma,' P 1           1'
       write(unitp,'(A4,I7,2x,A2,2x,A4,A,I4,4x,3f8.3)')
     . ('ATOM',i,atsym(isa(i)),'STO ',chain(i),i,
     . (rclas(ia,i)/Ang,ia=1,3), i=1,na)
       if (nac>9999) then
          write(unitp,'(A4,I7,2x,A4,A4,I5,4x,3f8.3)')
     .         ('ATOM',na+i,atname(i),aaname(i),aanum(i),
     .         (rclas(ia,na+i)/Ang,ia=1,3), i=1,nac)
       else
          write(unitp,'(A4,I7,2x,A4,A4,A,I4,4x,3f8.3)')
     .         ('ATOM',na+i,atname(i),aaname(i),chain(na+i),aanum(i),
     .         (rclas(ia,na+i)/Ang,ia=1,3), i=1,nac)
       endif
       write(unitp,'(A4,I7,2x,A2,2X,A4,A,I4,4x,3f8.3)')
     . ('ATOM',na+nac+i,atsym(resonance(1)%isa(i)),'LA  ','A',
     .  i,(rlink(1,ia,i)/Ang,ia=1,3), i=1,numlink)
       write(unitp,'(A3)') 'END'
       call io_close(unitp)

c writes in .crd file at begining 
      inquire( file=fnamec, exist=foundc )
	if (foundc) then
      call io_assign(unitc)
      open( unitc, file=fnamec, status='old' )
      read(unitc,'(a30)',err=1,end=1) title
      call io_close(unitc)
      ch=title(1:4)
      if(ch.ne.'File') then
      call io_assign(unitc)
      open( unitc, file=fnamec, form = 'formatted', status='unknown')
      write(unitc,'(20a)') 'File CRD'
      call io_close(unitc)
      endif
	else
      call io_assign(unitc)
      open( unitc, file=fnamec, form = 'formatted', status='unknown')  
      write(unitc,'(20a)') 'File CRD'
      call io_close(unitc)
	endif
      frstme = .false.
      endif !frstme
	call pxfflush(6)

c writes in .crd acumulately 
      if(MOD(istp,wricoord).eq.0) then         
      call io_assign(unitc)
      open( unitc, file=fnamec, form = 'formatted', position='append',
     .      status='unknown')
       write(unitc,'(10F8.3)') (rclas(1:3,i)/Ang,i=1,natot),
     .                                  (rlink(1,1:3,i)/Ang,i=1,numlink)
      call io_close(unitc)
      endif

c writes actual coords in PDB format in .last.pdb file 
      call io_assign(unitl)
      open( unitl, file=fnamel, form = 'formatted', status='unknown')
       write(unitl,'(A,3(F9.3),3(F7.2),A)')'CRYST1',a/Ang,b/Ang,c/Ang,
     .       alpha,beta,gamma,' P 1           1'
       write(unitl,'(A4,I7,2x,A2,2x,A4,A,I4,4x,3f8.3)')
     . ('ATOM',i,atsym(isa(i)),'STO ',chain(i),i,
     . (rclas(ia,i)/Ang,ia=1,3), i=1,na)
       if (nac>9999) then
          write(unitl,'(A4,I7,2x,A4,A4,I5,4x,3f8.3)')
     .         ('ATOM',na+i,atname(i),aaname(i),aanum(i),
     .         (rclas(ia,na+i)/Ang,ia=1,3), i=1,nac)
       else
          write(unitl,'(A4,I7,2x,A4,A4,A,I4,4x,3f8.3)')
     .         ('ATOM',na+i,atname(i),aaname(i),chain(na+i),aanum(i),
     .         (rclas(ia,na+i)/Ang,ia=1,3), i=1,nac)
       endif

       write(unitl,'(A4,I7,2x,A2,2X,A4,A,I4,4x,3f8.3)')
     . ('ATOM',na+nac+i,atsym(resonance(1)%isa(i)),'LA  ','A',
     .  i,(rlink(1,ia,i)/Ang,ia=1,3), i=1,numlink)
       write(unitl,'(A3,i5)') 'END',istp
       call io_close(unitl)

	return
 1      stop 'write: problem reading from CRD file'
	end

c******************************************************************************
c subrutine that writes the reaction coordinate and its PDB

       subroutine wrirtc(slabel,Etots,rtot,constrpaso,na,nac,natot,
     .                   rclas,atname,aaname,aanum,nesp,atsym,isa)

      use precision, only : dp
      implicit          none
      character         slabel*30, paste*30
      integer           na,nac,natot,istp,wricoord
      integer           nesp,isa(na)
      real(dp)  rclas(3,natot)
      character*2       atsym(nesp)
      real(dp)  Etots, rtot   
      external          paste
 
      character         fnamee*40,fnamec*40 
      logical           frstme
      integer           i,j,ia,unite,unitc
      real(dp)  eV, Ang
      save              frstme,Ang,eV,unite,fnamee,unitc,fnamec
      data              frstme /.true./
      integer           aanum(nac),constrpaso
      character*4       aaname(nac),atname(nac)

c -------------------------------------------------------------------

      if ( frstme ) then
        Ang  = 1.d0 / 0.529177d0
        eV  = 1.d0 / 13.60580d0
      fnamec = paste(slabel,'.ctr.pdb')
      fnamee = paste(slabel,'.ctr.ene')
      frstme = .false.
      endif

c writes .rce file
      call io_assign(unite)
      open( unite, file=fnamee, form = 'formatted', position='append',
     .      status='unknown')
       write(unite,'(F10.4,2x,F14.7)') rtot, Etots/eV
      call io_close(unite)

c wirtes .rcg file
      call io_assign(unitc)
      open( unitc, file=fnamec, form = 'formatted', position='append',
     .      status='unknown')
       write(unitc,'(A4,I7,2x,A2,2x,A4,A,I4,4x,3f8.3)')
     . ('ATOM',i,atsym(isa(i)),'STO ','A',i,
     . (rclas(ia,i)/Ang,ia=1,3), i=1,na)
       write(unitc,'(A4,I7,2x,A4,A4,A,I4,4x,3f8.3)')
     . ('ATOM',na+i,atname(i),aaname(i),'A',aanum(i),
     . (rclas(ia,na+i)/Ang,ia=1,3), i=1,nac)
      write(unitc,'(A3,i3)') 'END',constrpaso
      call io_close(unitc) 

      return
      end

c*********************************************************************************
c writes the slabel.siesta.fdf with its modifications
      subroutine write_input_debug(siesta_qmmmslabel,na_mm,na_qm,
     .          nspec,isa,rclas,mmcell)

      use precision, only : dp
      use fdf
      use m_qmmm_fdf, only : chrlen_qmmm
      use sys
      use siesta_qmmm_options, only : fc_idyn
      use linkatoms, only : numlink, rlink, num_resonances, resonance

      implicit none

      integer :: na_mm, na_qm, nspec   
      real(dp)  rclas(3,na_qm+na_mm)
      integer, dimension(na_qm) :: isa

      logical :: add_latvecs_flag
      real(dp), dimension(3,3) :: mmcell

      character siesta_qmmmslabel*30

C  Internal variables 
      character
     .  filein*30,fileout*43, line(100000)*150, paste*43,
     .  lineaux*150 

      integer i, j, nlines, length(100000), unit, iunit
      integer resonance_id

      external paste

      real(dp) x1, x2, x3
      integer fdf_isa
      integer :: local_ia1, local_ia2

      real(dp) :: fscale
      character :: acf*22, acf_default*22

      logical :: paobasis_read_flag
    
#ifndef QMMM_BSC
      logical :: leqi
#endif

      paobasis_read_flag=.false.

c Read slabel.fdf
      filein = paste(siesta_qmmmslabel,'.fdf')
      call io_assign(iunit)
      open( iunit, file=filein)
      i=1
 10   continue
      read(iunit, err=20,end=20,fmt='(a)') line(i)
      call chrlen_qmmm(line(i),0,length(i))
      if (length(i) .ne. 0) i=i+1
      goto 10
 20   continue
      nlines=i-1
      call io_close(iunit)

! Format of atomic coordinates
      acf_default = 'Ang'
      acf = fdf_string('AtomicCoordinatesFormat',acf_default)

      if (leqi(acf,'NotScaledCartesianBohr') .or. 
     . leqi(acf,'Bohr') ) then
      fscale = 1.0_dp
      write(6,'(a,a)') 
     . 'write_input_debug: Atomic-coordinates input format  = ',
     . '    Cartesian coordinates (in Bohr)'
      else if (leqi(acf,'NotScaledCartesianAng') .or.
     .    leqi(acf,'Ang') ) then
         fscale = 0.529177_dp
         write(6,'(a,a)')
     .    'write_input_debug: Atomic-coordinates input format  = ',
     .    '    Cartesian coordinates (in Ang)'
      else
         write(6,"(/,'write_input_debug: ',73(1h*))")
         write(6,"('write_input_debug:                  INPUT ERROR')")
         write(6,'(2a)') 
     .    'write_input_debug: You must use one of the following',
     .    'coordinate options:'
         write(6,'(a)') 
     .   'write_input_debug:     - NotScaledCartesianBohr (or '//
     .                          'Bohr)'
         write(6,'(a)') 
     .    'write_input_debug:     - NotScaledCartesianAng (or '//
     .                    'Ang) '
         write(6,"('read: ',73(1h*))")
         call die
      endif


      do resonance_id=1,num_resonances

C     Write slabel.siesta.fdf
         fileout = trim(resonance(resonance_id)%path)//
     .                        'INPUT_DEBUG'

         call io_assign(unit)
         open( unit, file=fileout)

         write(unit,FMT='(a)')'LatticeConstant   1.0   Bohr'
         write(unit,FMT='(a)')'%block LatticeVectors'
         do i=1,3
            write(unit,FMT='(3(2x,f12.6))')(mmcell(j,i), j=1,3)
         enddo
         write(unit,FMT='(a,/)')'%endblock LatticeVectors'
         i=0    
         do 
            i=i+1
            if (i>nlines) exit
            lineaux=line(i)
            if (lineaux(1:10).eq.'SystemName') then
               write(unit,'(a,2x,a)') 'SystemName',siesta_qmmmslabel 
            elseif(lineaux(1:12).eq.'MD.UseSaveXV') then
               ! Do nothing
            elseif(lineaux(1:12).eq.'DM.UseSaveDM') then
               ! Do nothing
            elseif(lineaux(1:11).eq.'WriteMDXmol') then
               write(unit,'(a)') 'WriteMDXmol .false.'
            elseif(lineaux(1:11).eq.'WriteMDXmol') then
               write(unit,'(a)') 'WriteMDXmol .false.'
            elseif(lineaux(1:15).eq.'NumberOfSpecies') then
              write(unit,'(a,2X,I5)')'NumberOfSpecies   ',nspec+1 
            elseif(lineaux(1:13).eq.'NumberOfAtoms') then
               write(unit,'(a,3x,I6)') 'NumberOfAtoms ',na_qm+na_mm+
     .                       numlink
            elseif(lineaux(1:15).eq.'LatticeConstant') then
               ! Do nothing
            elseif(lineaux(1:21).eq.'%block LatticeVectors') then
               ! Jump 4 lines
               i=i+4
            elseif(lineaux(1:30).eq.'%endblock ChemicalSpeciesLabel') 
     .                               then
               write(unit,'(I2,2X,A)') nspec+1,'-1  MM'
               write(unit,'(a)')'%endblock ChemicalSpeciesLabel' 
               write(unit,'(a)')
     .                    '%block AtomicCoordinatesAndAtomicSpecies'

               if (na_qm>0) then
                  do j = 1, na_qm
                     write(iunit,'(1X,3F12.6,I4)')fscale*rclas(1:3,j),
     .                           isa(j)
                  enddo
               endif

               if (na_mm>0) then
                  do j = 1, na_mm
                     write(iunit,'(1X,3F12.6,I4)')
     .                    fscale*rclas(1:3,na_qm+j),nspec+1
                  enddo
               endif

               if (numlink>0) then
                  do j=1,numlink
                     write(iunit,'(1X,3F12.6,I4)')fscale*
     .                    rlink(resonance_id,1:3,j),
     .                    resonance(resonance_id)%isa(j)
                  enddo
               endif

               write(unit,'(a)')
     .                    '%endblock AtomicCoordinatesAndAtomicSpecies'

            elseif(lineaux(1:40).eq.
     .            '%block AtomicCoordinatesAndAtomicSpecies' ) then
               do
                 i=i+1
                 lineaux=line(i)
                 if (lineaux(1:43).eq.
     .              '%endblock AtomicCoordinatesAndAtomicSpecies') exit
               enddo
            elseif(lineaux(1:19).eq.'%endblock PAO.Basis') then
               write(unit,FMT='(A)')'MM    0    0.00000'
               write(unit,FMT='(a)')'%endblock PAO.Basis' 
               paobasis_read_flag=.true.
            elseif(lineaux(1:22).eq.'%block SoluteAtomTypes') then
                do
                 i=i+1
                 lineaux=line(i)
                 if (lineaux(1:25).eq.
     .              '%endblock SoluteAtomTypes') exit
               enddo
            elseif(lineaux(1:19).eq.'%block SolventInput') then
                do
                 i=i+1
                 lineaux=line(i)
                 if (lineaux(1:22).eq.
     .              '%endblock SolventInput') exit
               enddo
            elseif(lineaux(1:19).eq.'%block CutOffRadius') then
                do
                 i=i+1
                 lineaux=line(i)
                 if (lineaux(1:22).eq.
     .              '%endblock CutOffRadius') exit
               enddo
            elseif(lineaux(1:20).eq.'NumberOfSolventAtoms') then
               ! Do nothing
            elseif(lineaux(1:12).eq.'MD.TypeOfRun') then
               if (fc_idyn) then
! First and last atoms to displace for calculation of force constants
                  local_ia1=fdf_integer('MD.FCfirst',1)
                  local_ia2=fdf_integer('MD.FClast',na_mm+na_qm)
                  if (local_ia1>na_mm+na_qm) local_ia1=na_mm
                  if (local_ia2>na_mm+na_qm) local_ia2=na_mm
                  write(iunit,'(a)')'MD.TypeOfRun    FC'
                  write(iunit,'(a,2x,I6)')'MD.FCfirst',local_ia1
                  write(iunit,'(a,2x,I6)')'MD.FClast',local_ia2
               else
                  write(unit,'(a)') trim(lineaux)
               endif
            elseif (lineaux(1:10).eq.'MD.FCfirst') then
               ! Do nothing
            elseif (lineaux(1:9).eq.'MD.FClast') then
               ! Do nothing
            elseif (lineaux(1:10).eq.'LongOutput') then
               ! Do nothing
            elseif (lineaux(1:12).eq.'WriteKpoints') then
               ! Do nothing
            elseif (lineaux(1:16).eq.'WriteEigenvalues') then
               ! Do nothing
            elseif (lineaux(1:11).eq.'WriteKbands') then
               ! Do nothing
            elseif (lineaux(1:10).eq.'WriteBands') then
               ! Do nothing
            elseif (lineaux(1:18).eq.'WriteWaveFunctions') then
               ! Do nothing
            elseif (lineaux(1:16).eq.'WriteMullikenPop') then
               write(iunit,'(a)')'WriteMullikenPop  0'
            elseif (lineaux(1:27).eq.'%block LocalDensityOfStates') then
!     Jump 2 more lines to end of block
               i=i+2
            elseif (lineaux(1:31).eq.'%block ProjectedDensityOfStates') 
     .                        then
!     Jump 2 more lines to end of block
               i=i+2
            elseif (lineaux(1:24).eq.'%block PolarizationGrids') then
!     Jump 4 more lines to end of block
               i=i+4
            elseif (lineaux(1:10).eq.'COOP.Write') then
               ! Do nothing
            elseif (lineaux(1:7).eq.'SaveRho') then
               ! Do nothing
            elseif (lineaux(1:12).eq.'SaveDeltaRho') then
               ! Do nothing
            elseif (lineaux(1:26).eq.'SaveElectrostaticPotential') then
               ! Do nothing
            elseif (lineaux(1:18).eq.'SaveTotalPotential') then
               ! Do nothing
            elseif (lineaux(1:15).eq.'SaveIonicCharge') then
               ! Do nothing
            elseif (lineaux(1:15).eq.'SaveTotalCharge') then
               ! Do nothing
            elseif (lineaux(1:18).eq.'OpticalCalculation') then
               ! Do nothing
            elseif (lineaux(1:18).eq.'BornCharge') then
               ! Do nothing
            else
               write(unit,'(a)') lineaux(1:length(i)) 
            endif
         enddo

         if (.not.paobasis_read_flag) then
            write(unit,FMT='(a)')'%block PAO.Basis'          
            write(unit,FMT='(A)')'MM    0    0.00000'
            write(unit,FMT='(a)')'%endblock PAO.Basis'             
         endif

         write(unit,'(a)') 'WriteXML .false.'
         write(unit,'(a)') 'MD.UseSaveXV .false.'
         write(unit,'(a)') 'DM.UseSaveDM .false.'

         call io_close(unit)

      enddo

        end

c writes the slabel.siesta.fdf with its modifications
      subroutine write_input_debug_mm(siesta_qmmmslabel,na_mm,rclas,
     .      mmcell)

      use precision, only : dp
      use fdf
      use m_qmmm_fdf, only : chrlen_qmmm
      use sys

      implicit none

      integer :: na_mm
      real(dp)  rclas(3,na_mm)
      real(dp) :: mmcell(3,3)

      character siesta_qmmmslabel*30

C  Internal variables 
      character
     .  filein*30,fileout*43, line(100000)*150, paste*43, lineaux*150 

      integer i, j, nlines, length(100000), unit, iunit
      integer resonance_id

      external paste

      real(dp) x1, x2, x3
      integer fdf_isa

      real(dp) :: fscale
      character :: acf*22, acf_default*22
#ifndef QMMM_BSC
      logical :: leqi
#endif

c Read slabel.fdf
      filein = paste(siesta_qmmmslabel,'.fdf')
      call io_assign(iunit)
      open( iunit, file=filein)
      i=1
 10   continue
      read(iunit, err=20,end=20,fmt='(a)') line(i)
      call chrlen_qmmm(line(i),0,length(i))
      if (length(i) .ne. 0) i=i+1
      goto 10
 20   continue
      nlines=i-1
      call io_close(iunit)

! Format of atomic coordinates
      acf_default = 'Ang'
      acf = fdf_string('AtomicCoordinatesFormat',acf_default)

      if (leqi(acf,'NotScaledCartesianBohr') .or. 
     . leqi(acf,'Bohr') ) then
      fscale = 1.0_dp
      write(6,'(a,a)') 
     . 'write_input_debug_mm: Atomic-coordinates input format  = ',
     . '    Cartesian coordinates (in Bohr)'
      else if (leqi(acf,'NotScaledCartesianAng') .or.
     .    leqi(acf,'Ang') ) then
         fscale = 0.529177_dp
         write(6,'(a,a)')
     .    'write_input_debug_mm: Atomic-coordinates input format  = ',
     .    '    Cartesian coordinates (in Ang)'
      else
         write(6,"(/,'write_input_debug_mm: ',73(1h*))")
         write(6,
     .         "('write_input_debug_mm:                  INPUT ERROR')")
         write(6,'(2a)') 
     .    'write_input_debug_mm: You must use one of the following',
     .    'coordinate options:'
         write(6,'(a)') 
     .   'write_input_debug_mm:     - NotScaledCartesianBohr (or '//
     .                          'Bohr)'
         write(6,'(a)') 
     .    'write_input_debug_mm:     - NotScaledCartesianAng (or '//
     .                    'Ang) '
         write(6,"('read: ',73(1h*))")
         call die
      endif

C     Write slabel.siesta.fdf
         fileout = 'INPUT_DEBUG'

         call io_assign(unit)
         open( unit, file=fileout)

         write(unit,FMT='(a)')'LatticeConstant   1.0   Bohr'
         write(unit,FMT='(a)')'%block LatticeVectors'
         do i=1,3
            write(unit,FMT='(3(2x,f12.6))')(mmcell(j,i), j=1,3)
         enddo
         write(unit,FMT='(a)')'%endblock LatticeVectors'
         i=0    
         do 
            i=i+1
            if (i>nlines) exit
            lineaux=line(i) 
            if (lineaux(1:10).eq.'SystemName') then
               write(unit,'(a,2x,a)') 'SystemName',siesta_qmmmslabel 
            elseif(lineaux(1:12).eq.'MD.UseSaveXV') then
               ! Do nothing
            elseif(lineaux(1:12).eq.'DM.UseSaveDM') then
               ! Do nothing
            elseif(lineaux(1:23).eq.'AtomicCoordinatesFormat') then
               ! Do nothing
            elseif(lineaux(1:11).eq.'SystemLabel') then
               write(unit,'(a,2x,a)') 'SystemLabel',siesta_qmmmslabel
            elseif(lineaux(1:11).eq.'WriteMDXmol') then
               write(unit,'(a)') 'WriteMDXmol .false.'
            elseif(lineaux(1:20).eq.'NumberOfSolventAtoms') then

               write(unit,'(a)')'NumberOfSpecies   1' 

               write(unit,'(a)')'%block ChemicalSpeciesLabel'
               write(unit,'(A)')' 1   -1  MM'
               write(unit,'(a)')'%endblock ChemicalSpeciesLabel' 
               
               write(unit,'(a,3x,I6)') 'NumberOfAtoms ',na_mm

               write(unit,'(a,a)') 'AtomicCoordinatesFormat  ',acf

               write(unit,'(a)')
     .                    '%block AtomicCoordinatesAndAtomicSpecies'
               if (na_mm>0) then
                  do j = 1, na_mm
                     write(iunit,'(1X,3F12.6,I4)')fscale*rclas(1:3,j),
     .                    1
                  enddo
               endif
               write(unit,'(a)')
     .                    '%endblock AtomicCoordinatesAndAtomicSpecies'

               write(unit,FMT='(a)')'%block PAO.Basis'
               write(unit,FMT='(A)')' MM    0    0.00000'
               write(unit,FMT='(a)')'%endblock PAO.Basis'
            elseif(lineaux(1:15).eq.'LatticeConstant') then
               ! Do nothing
            elseif(lineaux(1:21).eq.'%block LatticeVectors') then
               ! Jump 4 lines
               i=i+4
            elseif(lineaux(1:17).eq.'LatticeParameters') then
               ! Jump 2 lines
               i=i+2
            elseif(lineaux(1:15).eq.'NumberOfSpecies') then
C              Do nothing               
            elseif(lineaux(1:27).eq.'%block ChemicalSpeciesLabel') then
               do
                 i=i+1
                 lineaux=line(i)
                 if (lineaux(1:30).eq.
     .              '%endblock ChemicalSpeciesLabel') exit
               enddo
            elseif(lineaux(1:13).eq.'NumberOfAtoms') then
               ! Do nothing
            elseif(lineaux(1:40).eq.
     .             '%block AtomicCoordinatesAndAtomicSpecies') then
               do
                 i=i+1
                 lineaux=line(i)
                 if (lineaux(1:43).eq.
     .              '%endblock AtomicCoordinatesAndAtomicSpecies') exit
               enddo
            elseif(lineaux(1:16).eq.'%block PAO.Basis') then
                do
                 i=i+1
                 lineaux=line(i)
                 if (lineaux(1:19).eq.
     .              '%endblock PAO.Basis') exit
               enddo
            elseif(lineaux(1:22).eq.'%block SoluteAtomTypes') then
                do
                 i=i+1
                 lineaux=line(i)
                 if (lineaux(1:25).eq.
     .              '%endblock SoluteAtomTypes') exit
               enddo
            elseif(lineaux(1:19).eq.'%block SolventInput') then
                do
                 i=i+1
                 lineaux=line(i)
                 if (lineaux(1:22).eq.
     .              '%endblock SolventInput') exit
               enddo
            elseif(lineaux(1:19).eq.'%block CutOffRadius') then
                do
                 i=i+1
                 lineaux=line(i)
                 if (lineaux(1:22).eq.
     .              '%endblock CutOffRadius') exit
               enddo
            elseif(lineaux(1:20).eq.'NumberOfSolventAtoms') then
               ! Do nothing
            elseif (lineaux(1:10).eq.'LongOutput') then
               ! Do nothing
            elseif (lineaux(1:12).eq.'WriteKpoints') then
               ! Do nothing
            elseif (lineaux(1:16).eq.'WriteEigenvalues') then
               ! Do nothing
            elseif (lineaux(1:11).eq.'WriteKbands') then
               ! Do nothing
            elseif (lineaux(1:10).eq.'WriteBands') then
               ! Do nothing
            elseif (lineaux(1:18).eq.'WriteWaveFunctions') then
               ! Do nothing
            elseif (lineaux(1:16).eq.'WriteMullikenPop') then
               write(iunit,'(a)')'WriteMullikenPop  0'
            elseif (lineaux(1:27).eq.'%block LocalDensityOfStates') then
!     Jump 2 more lines to end of block
               i=i+2
            elseif (lineaux(1:31).eq.'%block ProjectedDensityOfStates') 
     .                        then
!     Jump 2 more lines to end of block
               i=i+2
            elseif (lineaux(1:24).eq.'%block PolarizationGrids') then
!     Jump 4 more lines to end of block
               i=i+4
            elseif (lineaux(1:10).eq.'COOP.Write') then
               ! Do nothing
            elseif (lineaux(1:7).eq.'SaveRho') then
               ! Do nothing
            elseif (lineaux(1:12).eq.'SaveDeltaRho') then
               ! Do nothing
            elseif (lineaux(1:26).eq.'SaveElectrostaticPotential') then
               ! Do nothing
            elseif (lineaux(1:18).eq.'SaveTotalPotential') then
               ! Do nothing
            elseif (lineaux(1:15).eq.'SaveIonicCharge') then
               ! Do nothing
            elseif (lineaux(1:15).eq.'SaveTotalCharge') then
               ! Do nothing
            elseif (lineaux(1:18).eq.'OpticalCalculation') then
               ! Do nothing
            elseif (lineaux(1:18).eq.'BornCharge') then
               ! Do nothing
            else
               write(unit,'(a)') lineaux(1:length(i)) 
            endif
         enddo

         write(unit,'(a)') 'WriteXML .false.'
         write(unit,'(a)') 'MD.UseSaveXV .false.'
         write(unit,'(a)') 'DM.UseSaveDM .false.'

         call io_close(unit)

        end


c*********************************************************************************
c writes the slabel.siesta.fdf with its modifications
      subroutine write_lab(slabel, siestaslabel, mm, na_qm, isa, rqm,
     .                               null_mm_charges)

      use precision, only : dp
      use fdf
      use m_qmmm_fdf, only : chrlen_qmmm
      use sys
      use linkatoms, only : numlink, rlink, num_resonances, resonance

      implicit none

      integer :: na_qm    
      real(dp)  rqm(3,na_qm)
      integer, dimension(na_qm) :: isa

      character slabel*30,siestaslabel*30

      logical mm, null_mm_charges
C  Internal variables 
      character
     .  filein*30,fileout*43, line(100000)*150, paste*43,
     .  chsn*10,chsl*11,chtr*12,chmdx*11,chna*13, 
     .  chlv*23, chcas*40,lineaux*150 

      integer i, j, nlines, length(100000), unit, iunit
      integer resonance_id

      external paste

      real(dp) x1, x2, x3
      integer fdf_isa

      real(dp) :: fscale
      character :: acf*22, acf_default*22      
#ifndef QMMM_BSC
      logical :: leqi
#endif

c Read slabel.fdf
      filein = paste(slabel,'.fdf')
      call io_assign(iunit)
      open( iunit, file=filein)
      i=1
 10   continue
      read(iunit, err=20,end=20,fmt='(a)') line(i)
      call chrlen_qmmm(line(i),0,length(i))
      if (length(i) .ne. 0) i=i+1
      goto 10
 20   continue
      nlines=i-1
      call io_close(iunit)

! Format of atomic coordinates
      acf_default = 'Ang'
      acf = fdf_string('AtomicCoordinatesFormat',acf_default)

      if (leqi(acf,'NotScaledCartesianBohr') .or. 
     . leqi(acf,'Bohr') ) then
      fscale = 1.0_dp
      write(6,'(a,a)') 
     . 'write_lab: Atomic-coordinates input format  = ',
     . '    Cartesian coordinates (in Bohr)'
      else if (leqi(acf,'NotScaledCartesianAng') .or.
     .    leqi(acf,'Ang') ) then
         fscale = 0.529177_dp
         write(6,'(a,a)')
     .    'write_lab: Atomic-coordinates input format  = ',
     .    '    Cartesian coordinates (in Ang)'
      else
         write(6,"(/,'write_lab: ',73(1h*))")
         write(6,"('write_lab:                  INPUT ERROR')")
         write(6,'(2a)') 'write_lab: You must use one of the following',
     .    'coordinate options:'
         write(6,'(a)') 'write_lab:     - NotScaledCartesianBohr (or '//
     .                          'Bohr)'
         write(6,'(a)') 'write_lab:     - NotScaledCartesianAng (or '//
     .                    'Ang) '
         write(6,"('read: ',73(1h*))")
         call die
      endif

      do resonance_id=1,num_resonances

C     Write slabel.siesta.fdf
         fileout = paste(trim(resonance(resonance_id)%path)//
     .                         siestaslabel,'.fdf')

         call io_assign(unit)
         open( unit, file=fileout)

         i=0    
         do 
            i=i+1
            if (i>nlines) exit
            lineaux=line(i)
            chsn=lineaux(1:10)
            chsl=lineaux(1:11)
            chtr=lineaux(1:12)
            chna=lineaux(1:13)
            chmdx=lineaux(1:11) 
            chcas=lineaux(1:40)
            chlv=lineaux(1:21)
            if (chsn.eq.'SystemName') then
               write(unit,'(a,2x,a)') 'SystemName',siestaslabel 
            elseif(chsl.eq.'SystemLabel') then
               write(unit,'(a,2x,a)') 'SystemLabel',siestaslabel
            elseif(chtr.eq.'MD.TypeOfRun') then
               write(unit,'(a)') 'MD.TypeOfRun forces' 
            elseif(chmdx.eq.'WriteMDXmol') then
               write(unit,'(a)') 'WriteMDXmol .false.'
            elseif(chna.eq.'NumberOfAtoms') then
               write(unit,'(a,3x,I6)') 'NumberOfAtoms ',na_qm+numlink
            elseif(chlv.eq.'%block LatticeVectors') then
               write(unit,'(a)')'%block LatticeVectors' 
               do j=1,3
                  i=i+1
                  lineaux=line(i)
                  write(unit,'(a)')trim(lineaux)
               enddo
               i=i+1
               write(unit,'(a,3x,I6)')'%endblock LatticeVectors' 
            elseif(chcas.eq.'%block AtomicCoordinatesAndAtomicSpecies'
     .                         ) then

               write(unit,'(a)')lineaux(1:length(i))

               do j = 1, na_qm
                  i=i+1
                  write(iunit,'(1X,3F8.3,I4)')fscale*rqm(1:3,j),isa(j)
               enddo

               if (numlink>0) then
                  do j=1,numlink
                     write(iunit,'(1X,3F8.3,I4)')fscale*
     .                    rlink(resonance_id,1:3,j),
     .                    resonance(resonance_id)%isa(j)
                  enddo
               endif

               i=i+1
               lineaux=line(i)
               write(unit,'(a)')lineaux(1:length(i))

            elseif(lineaux(1:22).eq.'%block SoluteAtomTypes') then
                do
                 i=i+1
                 lineaux=line(i)
                 if (lineaux(1:25).eq.
     .              '%endblock SoluteAtomTypes') exit
               enddo
            elseif(lineaux(1:19).eq.'%block SolventInput') then
                do
                 i=i+1
                 lineaux=line(i)
                 if (lineaux(1:22).eq.
     .              '%endblock SolventInput') exit
               enddo
            elseif(lineaux(1:19).eq.'%block CutOffRadius') then
                do
                 i=i+1
                 lineaux=line(i)
                 if (lineaux(1:22).eq.
     .              '%endblock CutOffRadius') exit
               enddo

            elseif(lineaux(1:20).eq.'NumberOfSolventAtoms') then
               ! Do nothing
            else
               write(unit,'(a)') lineaux(1:length(i)) 
            endif
         enddo

         if (.not.null_mm_charges) then
            if (mm) write(unit,'(a)') 'ReadVext .true.'
            write(unit,'(a)') 'SaveTotalCharge .true.'
         endif
         write(unit,'(a)') 'WriteXML .false.'

         call io_close(unit)

      enddo

        end

c***********************************************************************************
	subroutine wrtcrd(natot,rclas,cell,lattice_type)

        use precision, only : dp
	use fdf
        use m_qmmm_fdf, only : fdf_block_qmmm     
	use sys
	implicit none
	integer natot
	integer i,unit,iunit,nconstr,iconstr,typeconstr(20)
	integer atmsconstr(20,20),ndists(20),npi
	real(dp) rclas(3,natot),coef(20,10),rtot(20)
	character exp
        character slabel*30, paste*30,fname*30
        external  paste
        logical   frstme,constrlog
        data      frstme /.true./
        data      constrlog /.true./
        integer at1,at2,at3,at4,at5,at6,at7,at8
        real(dp) rp(3),r12,r34,r56,r78,dist_v2,angle_v2,dihedro2_v2
        real(dp) dx12, dx32, dy12, dy32, dz12, dz32
        real(dp) dx23, dy23, dz23, dx43, dy43, dz43
        real(dp) dx34, dy34, dz34, dx56, dy56, dz56
        real(dp) dx78, dy78, dz78
        real(dp)  cell(3,3),kcell(3,3)
        character lattice_type
        save frstme,fname,atmsconstr,ndists,coef,nconstr,typeconstr,
     .       constrlog 

        call reccel(3,cell,kcell,0)

c 8 ways of defining the reaction coordinate 
c 1 = r1 - r2 coupled
c 2 = distance 
c 3 = angle
c 4 = dihedral
c 5 = r1 + r2 coupled
c 6 = ( r1 + r2 ) - ( r3 + r4 ) coupled
c 7 = distance atom to 2-4 atoms average 
c 8 = c1*r1 + c2*r2 + c3*r3 + ....

	if(frstme) then
c read variables
	if ( fdf_block_qmmm('WrtCrd',iunit) ) then
	read(iunit,'(A)',advance='no',err=100,end=100) exp
	if(exp.eq.'%') then 
	constrlog=.false.
	goto 10
	endif
	read(iunit,*,err=100,end=100) exp,nconstr
	if(nconstr.gt.20) then
	call die('wrtcrd: nconstr must be lower than 20')
	endif
	do iconstr=1,nconstr
	read(iunit,*,err=100,end=100) exp,typeconstr(iconstr)

         if     (typeconstr(iconstr).eq.1) then          
         read(iunit,*,err=100,end=100) exp,(atmsconstr(iconstr,i),i=1,4)
         elseif (typeconstr(iconstr).eq.2) then
         read(iunit,*,err=100,end=100) exp,(atmsconstr(iconstr,i),i=1,2)
         elseif (typeconstr(iconstr).eq.3) then
         read(iunit,*,err=100,end=100) exp,(atmsconstr(iconstr,i),i=1,3)
         elseif (typeconstr(iconstr).eq.4) then
         read(iunit,*,err=100,end=100) exp,(atmsconstr(iconstr,i),i=1,4)
         elseif (typeconstr(iconstr).eq.5) then
         read(iunit,*,err=100,end=100) exp,(atmsconstr(iconstr,i),i=1,4)
         elseif (typeconstr(iconstr).eq.6) then
         read(iunit,*,err=100,end=100) exp,(atmsconstr(iconstr,i),i=1,8) 
         elseif (typeconstr(iconstr).eq.7) then
         read(iunit,*,err=100,end=100) exp,(atmsconstr(iconstr,i),i=1,5)
         elseif (typeconstr(iconstr).eq.8) then
         read(iunit,*,err=100,end=100) exp,ndists(iconstr)

       if(ndists(iconstr).gt.10) then
       call die('wrtcrd: ndists in typeconstr 8 must not exceed 10')
       endif
       read(iunit,*,err=100,end=100) 
     . exp,(coef(iconstr,i),i=1,ndists(iconstr))
       read(iunit,*,err=100,end=100) 
     . exp,(atmsconstr(iconstr,i),i=1,ndists(iconstr)*2)
         else
         call die('wrtcrd: typeconstr must be 1-8')
         endif

	if(i.gt.20) then
	call die('wrtcrd: atoms with constrain must be lower than 20')
	endif

	enddo !nconstr
	else
	constrlog=.false.
	goto 10
	endif !fdf

c name file
        slabel = fdf_string( 'SystemLabel', 'siesta' )
        fname = paste(slabel,'.wrt')

10	continue
	frstme=.false.
	endif !frstme

c write only if constrlog is true
        if(constrlog) then

c change units
        rclas(1:3,1:natot)=rclas(1:3,1:natot)*0.529177

c loop over nconstr
        do iconstr=1,nconstr

        if (typeconstr(iconstr).eq.1) then
        at1=atmsconstr(iconstr,1)
        at2=atmsconstr(iconstr,2)
        at3=atmsconstr(iconstr,3)
        at4=atmsconstr(iconstr,4)  
 
        dx12=rclas(1,at1)-rclas(1,at2)
        dy12=rclas(2,at1)-rclas(2,at2)
        dz12=rclas(3,at1)-rclas(3,at2)
        call pbc_displ_vector(lattice_type,cell,kcell,dx12,dy12,
     .                         dz12)
        r12=dist_v2(dx12,dy12,dz12)

        dx34=rclas(1,at3)-rclas(1,at4)
        dy34=rclas(2,at3)-rclas(2,at4)
        dz34=rclas(3,at3)-rclas(3,at4)
        call pbc_displ_vector(lattice_type,cell,kcell,dx34,dy34,
     .                         dz34)
        r34=dist_v2(dx34,dy34,dz34)
        rtot(iconstr)=r34-r12

        elseif (typeconstr(iconstr).eq.2) then

        at1=atmsconstr(iconstr,1)
        at2=atmsconstr(iconstr,2)
 
        dx12=rclas(1,at1)-rclas(1,at2)
        dy12=rclas(2,at1)-rclas(2,at2)
        dz12=rclas(3,at1)-rclas(3,at2) 
        call pbc_displ_vector(lattice_type,cell,kcell,dx12,dy12,
     .                              dz12)
        rtot(iconstr)=dist_v2(dx12,dy12,dz12)

        elseif(typeconstr(iconstr).eq.3) then
        at1=atmsconstr(iconstr,1)
        at2=atmsconstr(iconstr,2)
        at3=atmsconstr(iconstr,3)

	dx12=rclas(1,at1)-rclas(1,at2)
	dy12=rclas(2,at1)-rclas(2,at2)
	dz12=rclas(3,at1)-rclas(3,at2)
        call pbc_displ_vector(lattice_type,cell,kcell,dx12,
     .                              dy12,dz12)
	dx32=rclas(1,at3)-rclas(1,at2)
	dy32=rclas(2,at3)-rclas(2,at2)
	dz32=rclas(3,at3)-rclas(3,at2)
        call pbc_displ_vector(lattice_type,cell,kcell,dx32,
     .                              dy32,dz32)
        rtot =angle_v2(dx12,dy12,dz12,dx32,dy32,dz32)

        elseif(typeconstr(iconstr).eq.4) then
        at1=atmsconstr(iconstr,1)
        at2=atmsconstr(iconstr,2)
        at3=atmsconstr(iconstr,3)
        at4=atmsconstr(iconstr,4)

        rtot(iconstr)=
     .   dihedro2_v2(rclas(1,at1),rclas(2,at1),rclas(3,at1),
     .   rclas(1,at2),rclas(2,at2),rclas(3,at2),rclas(1,at3),
     .   rclas(2,at3),rclas(3,at3),rclas(1,at4),rclas(2,at4),
     .   rclas(3,at4),cell,lattice_type)

        elseif (typeconstr(iconstr).eq.5) then
        at1=atmsconstr(iconstr,1)
        at2=atmsconstr(iconstr,2)
        at3=atmsconstr(iconstr,3)
        at4=atmsconstr(iconstr,4)

	dx12=rclas(1,at1)-rclas(1,at2)
	dy12=rclas(2,at1)-rclas(2,at2)
	dz12=rclas(3,at1)-rclas(3,at2)
        call pbc_displ_vector(lattice_type,cell,kcell,dx12,
     .                              dy12,dz12)
        r12=dist_v2(dx12,dy12,dz12)

	dx34=rclas(1,at3)-rclas(1,at4)
	dy34=rclas(2,at3)-rclas(2,at4)
	dz34=rclas(3,at3)-rclas(3,at4)
        call pbc_displ_vector(lattice_type,cell,kcell,dx34,
     .                              dy34,dz34)
        r34=dist_v2(dx34,dy34,dz34)

        rtot(iconstr)=r34+r12

        elseif (typeconstr(iconstr).eq.6) then
        at1=atmsconstr(iconstr,1)
        at2=atmsconstr(iconstr,2)
        at3=atmsconstr(iconstr,3)
        at4=atmsconstr(iconstr,4)
        at5=atmsconstr(iconstr,5)
        at6=atmsconstr(iconstr,6)
        at7=atmsconstr(iconstr,7)
        at8=atmsconstr(iconstr,8)

 	dx12=rclas(1,at1)-rclas(1,at2)
	dy12=rclas(2,at1)-rclas(2,at2)
	dz12=rclas(3,at1)-rclas(3,at2)
        call pbc_displ_vector(lattice_type,cell,kcell,dx12,
     .                              dy12,dz12)
        r12=dist_v2(dx12,dy12,dz12)

	dx34=rclas(1,at3)-rclas(1,at4)
	dy34=rclas(2,at3)-rclas(2,at4)
	dz34=rclas(3,at4)-rclas(3,at4)
        call pbc_displ_vector(lattice_type,cell,kcell,dx34,
     .                              dy34,dz34)
        r34=dist_v2(dx34,dy34,dz34)

	dx56=rclas(1,at5)-rclas(1,at6)
	dy56=rclas(2,at5)-rclas(2,at6)
	dz56=rclas(3,at5)-rclas(3,at6)
        call pbc_displ_vector(lattice_type,cell,kcell,dx56,
     .                              dy56,dz56)
        r56=dist_v2(dx56,dy56,dz56)

	dx78=rclas(1,at7)-rclas(1,at8)
	dy78=rclas(2,at7)-rclas(2,at8)
	dz78=rclas(3,at7)-rclas(3,at8)
        call pbc_displ_vector(lattice_type,cell,kcell,dx78,
     .                              dy78,dz78)
        r78=dist_v2(dx78,dy78,dz78)

        rtot(iconstr)=(r34+r12)-(r56+r78)

        elseif (typeconstr(iconstr).eq.7) then
        npi=0
        do i=2,5
        if(atmsconstr(iconstr,i).ne.0) npi=npi+1
        enddo
        if(npi.eq.0) stop 'constr opt: atoms in average can not be zero'
        at1=atmsconstr(iconstr,1)
        at2=atmsconstr(iconstr,2)
        at3=atmsconstr(iconstr,3)
        at4=atmsconstr(iconstr,4)
        at5=atmsconstr(iconstr,5)

        rp(1)=(rclas(1,at2)+rclas(1,at3)+rclas(1,at4)+rclas(1,at5))/npi
        rp(2)=(rclas(2,at2)+rclas(2,at3)+rclas(2,at4)+rclas(2,at5))/npi
        rp(3)=(rclas(3,at2)+rclas(3,at3)+rclas(3,at4)+rclas(3,at5))/npi

        dx12=rclas(1,at1)-rp(1)
        dy12=rclas(2,at1)-rp(2)
        dz12=rclas(3,at1)-rp(3)
        call pbc_displ_vector(lattice_type,cell,kcell,dx12,
     .                              dy12,dz12)
        rtot(iconstr)=dist_v2(dx12,dy12,dz12)

        elseif (typeconstr(iconstr).eq.8) then
        rtot(iconstr)=0.0
        do i=1,ndists(iconstr)
        at1=atmsconstr(iconstr,i)
        at2=atmsconstr(iconstr,i+1)

        dx12=rclas(1,at1)-rclas(1,at2)
        dy12=rclas(2,at1)-rclas(2,at2)
        dz12=rclas(3,at1)-rclas(3,at2)
        call pbc_displ_vector(lattice_type,cell,kcell,dx12,
     .                              dy12,dz12)
        rtot(iconstr)=dist_v2(dx12,dy12,dz12)
        rtot(iconstr)=rtot(iconstr)+coef(iconstr,i)*12

        enddo !ndists
        endif !typeconstr
        enddo !nconstr

c write file
        call io_assign(unit)
        open( unit, file=fname, form = 'formatted', position='append',
     .        status='unknown')
        write(unit,'(20f14.8)') rtot(1:nconstr)
        call io_close(unit)

c change units 
        rclas(1:3,1:natot)=rclas(1:3,1:natot)/0.529177

        endif !constrlog
        return
 100    stop 'wrtcrd: problem reading WrtCrd block'
        end

c************************************************************************************************


