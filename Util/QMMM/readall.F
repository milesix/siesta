c Reads the required QM information

      subroutine read_qm(na_qm,nesp,isa_qm,iza_qm,rqm,atsym)

      use precision
      use sys
      use fdf
      use m_qmmm_fdf, only : fdf_block_qmmm 

      implicit          none
      character         acf*22, acf_default*22
      integer           i, j, iunit, na_qm, nesp, iscale 
      integer           isa_qm(na_qm),iza_qm(na_qm),atnum(nesp)
      real(dp)          rqm(3,na_qm), spin
      character*2       atsym(nesp)
#ifndef QMMM_BSC
      logical           leqi
#endif
      logical fix
      integer           nelec,charge
      real(dp)  charnet,charnet_def

C Format of atomic coordinates
      acf_default = 'Ang'
      acf = fdf_string('AtomicCoordinatesFormat',acf_default)

      if (leqi(acf,'NotScaledCartesianBohr') .or.
     .    leqi(acf,'Bohr') ) then
        iscale = 0
        write(6,'(a,a)')
     .   'read_qm: Atomic-coordinates input format  = ',
     .   '    Cartesian coordinates (in Bohr)'
      else if (leqi(acf,'NotScaledCartesianAng') .or.
     .         leqi(acf,'Ang') ) then
        iscale = 1
        write(6,'(a,a)')
     .   'qmmm_read: Atomic-coordinates input format  = ',
     .   '    Cartesian coordinates (in Ang)'
      else
        write(6,"(/,'qmmm_read: ',73(1h*))")
        write(6,"('qmmm_read:                  INPUT ERROR')")
        write(6,'(2a)') 'qmmm_read: You must use one of the following',
     .                            'coordinate options:'
        write(6,'(a)') 'qmmm_read:     - NotScaledCartesianBohr (or '//
     .                      'Bohr)'
        write(6,'(a)') 'qmmm_read:     - NotScaledCartesianAng (or '//
     .                    'Ang) '
        write(6,"('qmmm_read: ',73(1h*))")
       call die
      endif                                                   

c Read atomic coordinates and species
      if ( fdf_block_qmmm('AtomicCoordinatesAndAtomicSpecies',iunit) )
     .     then
         do i = 1, na_qm
            read(iunit,*,err=1,end=1) (rqm(j,i), j=1,3), isa_qm(i)
         enddo
      else
         call die("qmmm_read: You must specify the atomic coordinates")
      endif

c change coordinates to Siesta's format
c       Scale atomic coordinates
c       Coord. option = 0 => Do nothing
c       Coord. option = 1 => Multiply by 1./0.529177 (Ang --> Bohr)

      if (iscale .eq. 1) then
        rqm(1:3,1:na_qm) = rqm(1:3,1:na_qm) / 0.529177d0
      endif

c Read the atomic labels
      if ( fdf_block_qmmm('ChemicalSpeciesLabel',iunit) ) then
         do i=1,nesp
         read(iunit,*,err=2,end=2) j,atnum(i),atsym(i)
C change sign of atnum if there are ghost atoms
         if (atnum(i).lt.0) atnum(i)=-atnum(i)
         enddo
      else 
         call die("qmmm_read: You must specify the atomic labels")
      endif

c Assignates atomic number (iza_qm)
       iza_qm=0 
       do i=1,na_qm
          do j=1,nesp
             if(isa_qm(i).eq.j) iza_qm(i)=atnum(j)
          enddo
             if(iza_qm(i).eq.0) then
             call die("qmmm_read: There are atomos without "//
     .                       "atomic number")
             endif
       enddo

c Calculates the total number of electrons 
      charnet_def=0.0
      charnet=fdf_double('NetCharge',charnet_def)
      charge=int(charnet)
      nelec=0
      do i=1,na_qm
        if(iza_qm(i).eq.1)                   nelec=nelec+iza_qm(i)
        if(iza_qm(i).gt.2.and.iza_qm(i).lt.10)  nelec=nelec+iza_qm(i)-2
        if(iza_qm(i).gt.10.and.iza_qm(i).lt.18) nelec=nelec+iza_qm(i)-10
        if(iza_qm(i).gt.18.and.iza_qm(i).lt.36) nelec=nelec+iza_qm(i)-18
      enddo
      nelec=nelec-charge
        write(6,'(a,i5)')
     .  'qmmm_read: Total number of electrons        = ',nelec

C Check mumber of electrons with solute spin multiplicity
c$$$       fix = fdf_boolean('FixSpin',.false.)
c$$$       if (fix) then
c$$$       spin = fdf_double('TotalSpin',0.0d0)
c$$$       if(mod(spin,2.0_dp).eq.0.0.and.mod(nelec,2).ne.0) then
c$$$         call die("qmmm_read:
c$$$     . Problem between the number of electrons and the total spin")
c$$$       endif
c$$$       if(mod(spin,2.0_dp).ne.0.0.and.mod(nelec,2).eq.0) then
c$$$       call die("qmmm_read:
c$$$     . Problem between the number of electrons and the total spin")
c$$$       endif !spin
c$$$       endif !fix

       return
 1     stop 'qmmm_read: problem reading solute coordinates'
 2     stop 'qmmm_read: problem reading solute chemical specie label' 
        end

c**********************************************************************************
C Subroutine to read a crd file from an Amber run
	subroutine readcrd(na_qm,na_mm,masst,rclas,vat,foundcrd,foundvat)

        use precision
        use sys
        use fdf

	implicit none
	integer i,j,na_qm,na_mm,nat,natoms,iu
	real(dp) rclas(1:3,na_qm+na_mm),masst(na_qm+na_mm),
     .                   vat(1:3,na_qm+na_mm),time
	character sname*30,fname*30,paste*30,title*30,ch*4
	logical   foundcrd,foundvat
	external io_assign, io_close, paste

        foundcrd = .false.
        sname = fdf_string( 'SystemLabel', 'siesta' )
        fname = paste( sname, '.mdcrd' )

        inquire( file=fname, exist=foundcrd )

        if (foundcrd) then
          call io_assign( iu )
          open( iu, file=fname, status='old' )

	read(iu,err=1,end=1,fmt='(a30)') title
	ch=title(1:4)
	if(ch.eq.'File') then
	foundcrd=.false.
	else
	read(iu,*,err=1,end=1) natoms,time
	
	nat=na_qm+na_mm
	if(nat.ne.natoms)  stop 'Wrong number of atoms!!!' 
	
	read(iu,err=1,end=1,fmt='(6f12.7)') 
     .              ((rclas(i,j),i=1,3),j=1,natoms)
	read(iu,err=1,end=1,fmt='(6f12.7)') ((vat(i,j),i=1,3),j=1,natoms)

        rclas(1:3,1:natoms)=rclas(1:3,1:natoms)/0.5291772
        vat(1:3,1:natoms)=vat(1:3,1:natoms)/0.5291772*20.455/1000

        write(6,'(/a)') 
     .   'readcrd: Reading Coordinates and Velocities from CRD file'
	endif
	call io_close(iu)
        endif !foundcrd

c change velocities if using deuterium
	if(foundcrd) then
	do i=1,na_qm+na_mm
	if(masst(i).eq.2.0) vat(1:3,i)=vat(1:3,i)/sqrt(2.)
	enddo
	endif

c check if velocities are zero
	do j=1,3
	  do i=1,na_qm+na_mm
	    if(vat(j,i).ne.0.0) foundvat=.true.
	  enddo
	enddo

	return
 1      stop 'readcrd: problem reading CRD file'
	end

c**************************************************************************
