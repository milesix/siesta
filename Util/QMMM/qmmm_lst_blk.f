c subroutine that calculates the rcut and block QM-MM only in the first step

	subroutine qmmm_lst_blk(na_qm,na_mm,nroaa,atxres,rclas,
     .  rcorteqmmm,radiobloqmmm,blockqmmm,listqmmm,rcorteqm,slabel)

	use precision 
	use sys
	implicit none
        integer na_qm,na_mm,natot,iu,nroaa,atxres(20000)
        real(dp) cm(3,20000),dist,dist2,rcorteqm,
     .  rcorteqmmm,radiobloqmmm,rclas(3,na_qm+na_mm),Ang
        integer i,j,k,l,count,blockqmmm(na_mm),listqmmm(na_mm)
        logical bloqmmm,liqmmm,found
        character slabel*30,fname*30,paste*30
        external io_assign, io_close, paste
        Ang    = 1._dp / 0.529177_dp

c change units
        rclas=rclas/Ang

        natot=na_qm+na_mm

c calculate center of masses of all residues
	k=na_qm+1
	do i=1,nroaa
	 do j=1,atxres(i)
	cm(1:3,i)=cm(1:3,i)+rclas(1:3,k)
	k=k+1
	 enddo
	cm(1:3,i)=cm(1:3,i)/atxres(i)
	enddo

c rcut QM-MM
        found=.false.
        liqmmm=.false.
c find file name
        fname = paste(slabel,'.lst')
c check if the input file exists
        inquire( file=fname, exist=found )
       if (found) then
c Open file
        call io_assign( iu )
        open( iu, file=fname, status='old' )
c read listqmmm
        do i=1,na_mm
        read(iu,*,err=1,end=1) listqmmm(i)
        enddo
c Close file
        call io_close( iu )
        write(6,'(/a)') 'qm-mm: Reading QM-MM neighbours list from file'
       else
       if(rcorteqmmm.le.99.0) liqmmm=.true.
        if(liqmmm) then
c QM-MM neigh list
      write(6,'(/a,f12.6)') 
     .'qm-mm: cut off radius QM-MM (Ang):',rcorteqmmm
      if(rcorteqm.ge.rcorteqmmm) then
      call die("qm-mm: cut off QM-MM must be greater than cut off QM")
      endif
        dist=0.0
        dist2=rcorteqmmm**2
	listqmmm=1
	do l=1,na_qm
	k=1
	 do i=1,nroaa
	  dist=(rclas(1,l)-cm(1,i))**2+
     .         (rclas(2,l)-cm(2,i))**2+
     .         (rclas(3,l)-cm(3,i))**2
	  do j=1,atxres(i)
	   if(dist.le.dist2) listqmmm(k)=0
	   k=k+1
	  enddo
	 enddo
	enddo
	count=0
        do i=1,na_mm
          if(listqmmm(i).eq.0) count=count+1
        enddo
	write(6,'(/a,2x,i5)') 
     .  'qm-mm: Number of QM-MM interacting atoms:',count
c Open file
        call io_assign( iu )
        open( iu, file=fname, form='formatted', status='unknown' )
c write listqmmm
        do i=1,na_mm
        write(iu,*) listqmmm(i)
        enddo
c Close file
        call io_close( iu )
c        write(6,'(/a)') 'qm-mm: Writing QM-MM neighbours list in file'
        else ! liqmmm=.false.
        write(6,'(/a)') 'qm-mm: Warning QM-MM cutoff too large'
        endif
       endif

c block QM-MM
        found=.false.
        bloqmmm=.false.
c find file name
        fname = paste(slabel,'.blk')
c check if the input file exists
        inquire( file=fname, exist=found )
       if (found) then
c Open file
        call io_assign( iu )
        open( iu, file=fname, status='old' )
c read blockqmmm
        do i=1,na_mm
        read(iu,*,err=2,end=2) blockqmmm(i)
        enddo
c Close file
        call io_close( iu )
        write(6,'(/a)') 'qm-mm: Reading blocked QM-MM atoms from file'
       else
       if(radiobloqmmm.le.99.0) bloqmmm=.true.
        if(bloqmmm) then
c fixing MM atoms beyond block cut off  
        write(6,'(/a,f12.6)') 
     .  'qm-mm: cut off radius Block (Ang):',radiobloqmmm
        dist=0.0
        dist2=radiobloqmmm**2
        blockqmmm=1
        do l=1,na_qm
        k=1
         do i=1,nroaa
          dist=(rclas(1,l)-cm(1,i))**2+
     .         (rclas(2,l)-cm(2,i))**2+
     .         (rclas(3,l)-cm(3,i))**2
          do j=1,atxres(i)
           if(dist.le.dist2) blockqmmm(k)=0
           k=k+1
          enddo
         enddo
        enddo
        count=0
        do i=1,na_mm
        if(radiobloqmmm.eq.0.0) blockqmmm(i)=1
        if(blockqmmm(i).eq.0) count=count+1
        enddo
        write(6,'(/a,2x,i5)') 
     .  'qm-mm: Number of QM-MM free atoms:',count 
c Open file
        call io_assign( iu )
        open( iu, file=fname, form='formatted', status='unknown' )
c write blockqmmm
        do i=1,na_mm
        write(iu,*) blockqmmm(i)
        enddo
c Close file
        call io_close( iu )
c        write(6,'(/a)') 'qm-mm: Writing blocked QM-MM atoms in file'
        else ! bloqmmm=.false.
        write(6,'(/a)') 'qm-mm: Warning Block cutoff too large'
        endif
       endif

c change units
      rclas=rclas*Ang
      rcorteqm=rcorteqm*Ang       
      call pxfflush(6)

      return
 1    stop 'qm-mm: Problem reading from .lst file'
 2    stop 'qm-mm: Problem reading from .blk file'
      end

