c subroutine that read the constrained atom block 

	subroutine qmmm_fixed1(na_qm,na_mm,natot,nroaa,rclas,
     .                    blocklist,atname,aaname,aanum,wat)

        use precision
        use fdf
        use m_qmmm_fdf, only : fdf_block_qmmm
        use sys
	implicit none
	integer i,j,k,l,iunit,na_qm,na_mm,natot,nroaa,aanum(na_mm),
     .                 blocklist(natot)
	character*4 atname(na_mm),aaname(na_mm)
	real(dp) rclas(3,natot)
	logical wat

	integer nt,type(9),aa,con1(4,natot),con2(4,natot),ncon(4)
	integer nc,c1(na_qm),c2(na_qm)
	real(dp) r(3,na_mm),cut,cqm(3),mdist,dist,dist2
	character ch1*1,ch4*4,exp
	wat = .false.

c reads 'PositionConstraints' block
	if ( fdf_block_qmmm('PositionConstraints',iunit) ) then
        read(iunit,'(A)',advance='no',err=100,end=100) exp 
        if(exp.eq.'%') goto 50 
	read(iunit,*,err=100,end=100) exp,nt
	do i=1,nt
	  read(iunit,*,err=100,end=100) exp,type(i)

	if(na_mm.eq.0.and.type(i).gt.1) then
	call die('fixed: constrained type for only QM atoms must be only
     .     1')
	endif
	  if(type(i).le.2) then

       k=1
 10    continue
       if(k.gt.natot) then
       call die('fixed: sets of constraints must not exeed natot')
       endif
       read(iunit,'(A)',advance='no',err=100,end=100) exp
       if(exp.eq.'f'.or.exp.eq.'F'.or.exp.eq.'p'.or.exp.eq.'P') then
       read(iunit,*,err=100,end=100) exp,con1(i,k),exp,con2(i,k)
       k=k+1
       goto 10 
       else
       continue
       endif
       ncon(i)=k-1

	if(ncon(i).eq.0) then
	call die('fixed: sets of constraints must not be zero') 
	endif

	do k=1,ncon(i)
	if(con2(i,k).lt.con1(i,k)) then
	call die('fixed: sets of constraints must be defined in correct 
     .                 order')
	endif
	enddo

	  elseif(type(i).eq.6) then
            read(iunit,*,err=100,end=100) exp,aa
            read(iunit,*,err=100,end=100) exp,cut
	  elseif(type(i).lt.-2 .or. type(i).gt.7) then
	  call die('fixed: Wrong type of constraint')
	  endif !type
	enddo !nt
	write(6,'(/,a)') 'fixed: reading "PositionConstraints" block'
	else 
	goto 50
	endif !fdf

c blocklist assignation
      do i=1,nt
	if(type(i).eq.1) then
          do k=1,ncon(i)
        if(con2(i,k).gt.natot) then
        call die('fixed: atoms in constraint must not exeed natot')
        endif
            do j=con1(i,k),con2(i,k)
            blocklist(j)=1
            enddo
          enddo
	elseif(type(i).eq.-1) then
        blocklist=1
          do k=1,ncon(i)
        if(con2(i,k).gt.natot) then
        call die('fixed: atoms in constraint must not exeed natot')
        endif
            do j=con1(i,k),con2(i,k)
            blocklist(j)=0
            enddo
          enddo
	elseif(type(i).eq.2) then
          do k=1,ncon(i)
        if(con2(i,k).gt.nroaa) then
        call die('fixed: residues in constraint must not exeed nroaa')
        endif

            do j=con1(i,k),con2(i,k)
              do l=1,na_mm
                if(j.eq.aanum(l)) then
                blocklist(na_qm+l)=1
                endif
              enddo
            enddo
          enddo
	elseif(type(i).eq.-2) then
        blocklist=1
          do k=1,ncon(i)
        if(con2(i,k).gt.nroaa) then
        call die('fixed: residues in constraint must not exeed nroaa')
        endif
            do j=con1(i,k),con2(i,k)
              do l=1,na_mm
                if(j.eq.aanum(l)) then
                blocklist(na_qm+l)=0
                endif
              enddo
            enddo
          enddo
	elseif(type(i).eq.3) then
          do l=1,na_mm
            ch4=atname(l)
            ch1=ch4(1:1)
            if(ch1.ne.'H') then
            blocklist(na_qm+l)=1
            endif
          enddo
	elseif(type(i).eq.4) then
          do l=1,na_mm
            if(atname(l).eq.'CA'.or.
     .         atname(l).eq.'C' .or.
     .         atname(l).eq.'N') then
            blocklist(na_qm+l)=1
            endif
          enddo
	elseif(type(i).eq.5) then
          do l=1,na_mm
            if(aaname(l).ne.'HOH') then
            blocklist(na_qm+l)=1
            endif
          enddo
	elseif(type(i).eq.6) then
        r(1:3,1:na_mm)=rclas(1:3,na_qm+1:natot)*0.529177
        cqm=0.0
        k=0
        do j=1,na_mm
        if(aanum(j).eq.aa) then
        k=k+1
        cqm(1)=cqm(1)+r(1,j)
        cqm(2)=cqm(2)+r(2,j)
        cqm(3)=cqm(3)+r(3,j)
        endif
        enddo
        cqm(1:3)=cqm(1:3)/k

        mdist=0.0
        dist=0.0
        do j=1,na_mm
        if(aanum(j).eq.aa) then
        dist=(r(1,j)-cqm(1))**2+
     .       (r(2,j)-cqm(2))**2+
     .       (r(3,j)-cqm(3))**2
        if(dist.gt.mdist) mdist=dist
        endif
        enddo
        mdist=sqrt(mdist)
        cut=cut+mdist

        dist=0.0
        dist2=cut**2
        do j=1,na_mm
        if(aaname(j).eq.'HOH'.and.atname(j).eq.'O') then
                        dist=(r(1,j)-cqm(1))**2+
     .                       (r(2,j)-cqm(2))**2+
     .                       (r(3,j)-cqm(3))**2
                        if(dist.gt.dist2) then
                        blocklist(j+na_qm)=1
                        blocklist(j+1+na_qm)=1
                        blocklist(j+2+na_qm)=1
                        endif
        elseif(aaname(j).ne.'HOH') then
                        dist=(r(1,j)-cqm(1))**2+
     .                       (r(2,j)-cqm(2))**2+
     .                       (r(3,j)-cqm(3))**2
                        if(dist.gt.dist2) then
                        blocklist(j+na_qm)=1
                        endif
        endif
        enddo
	elseif(type(i).eq.7) then
	wat = .true.
	write(6,'(/,a)') 'fixed: Running restraining water cap'
	endif !type
      enddo !nt

 50    continue 
c read 'GeometryConstraints' block
c$$$	if ( fdf_block_qmmm('GeometryConstraints',bfdf) ) then
c$$$        k=1
c$$$ 60     continue
c$$$        if(k.gt.na_qm) then
c$$$        call die('fixed: sets of QM constraints must not exeed na_qm')
c$$$        endif
c$$$        read(iunit,'(A)',advance='no',err=200,end=200) exp
c$$$        if(exp.eq.'%'.and.k.eq.1) goto 70
c$$$        if(exp.eq.'%'.and.k.ne.1) goto 65
c$$$        if(exp.eq.'P'.or.exp.eq.'p') then
c$$$        read(iunit,*,err=200,end=200) exp,exp,c1(k),exp,c2(k)
c$$$        k=k+1
c$$$        goto 60
c$$$        else
c$$$        call die('fixed: wrong syntaxis in GeometryConstraints block')
c$$$        endif
c$$$ 65     continue
c$$$        nc=k-1
c$$$ 
c$$$        do i=1,nc
c$$$        if(c2(i).gt.na_qm) then
c$$$        call die('fixed: atoms in QM constraint must not exeed na_qm')
c$$$        endif
c$$$
c$$$        if(c2(i).lt.c1(i)) then
c$$$        call die('fixed: sets of QM constraints must be defined 
c$$$     .            in correct order')
c$$$        endif
c$$$
c$$$c       blocklist assignation
c$$$        do j=c1(i),c2(i)
c$$$	   blocklist(j)=1
c$$$        enddo
c$$$        enddo
c$$$
c$$$        write(6,'(/,a)') 'fixed: reading "GeometryConstraints" block'
c$$$	endif			!fdf 
c$$$
c$$$ 70    continue
       return
 100   stop 'fixed: Problem reading "PositionConstraints" block'
 200   stop 'fixed: Problem reading "GeometryConstraints" block'
       end

c*************************************************************
c subroutine that imposes fce and vel constraints 

       subroutine  qmmm_fixed2(na_qm,na_mm,natot,nfree,
     .             blocklist,blockqmmm,fdummy,cfdummy,vat)

        use precision, only: dp
	implicit none
	integer i,k,na_qm,na_mm,natot,nfree,
     .           blockqmmm(na_mm), blocklist(natot)
	real(dp) fdummy(3,natot),cfdummy(3,natot),vat(3,natot) 
        logical    frstme
        save       frstme
        data frstme /.true./

C copy fdummy array to cfdummy
	cfdummy = fdummy 

c nullify fce and vel 
        do i=1,na_qm
           if (blocklist(i).eq.1) then
           cfdummy(1:3,i)=0.0_dp
           vat(1:3,i)=0.0_dp
           endif
        enddo

        do i=1,na_mm
           if (blockqmmm(i).eq.1.or.blocklist(i+na_qm).eq.1) then
           cfdummy(1:3,i+na_qm)=0.0_dp
           vat(1:3,i+na_qm)=0.0_dp
           endif
        enddo

c set total number of free atoms
	if(frstme) then
	k=0
	do i=1,na_qm
	if((blocklist(i).eq.0)) k=k+1
	enddo
	do i=1,na_mm
	if((blocklist(i+na_qm).eq.0).and.(blockqmmm(i).eq.0)) k=k+1
	enddo
	nfree=k
	write(6,'(/a,2x,i5)') 'siesta-qmmm: Total Free Atoms:', nfree 
	frstme=.false.
	endif

	end

        subroutine add_siesta_qmmm_ntcon(cell,na_qm,na_mm,isa,
     .          amass,xa,blocklist,blockqmmm,ntcon)

        use precision, only: dp
#ifdef QMMM_BSC
        use m_fixed
#endif

        implicit none

        real(dp) cell(3,3)
        integer na_qm, na_mm, isa(na_qm+na_mm)
        real(dp) amass(na_qm+na_mm), xa(na_mm+na_qm)
        integer blockqmmm(na_qm+na_mm), blocklist(na_qm+na_mm)
        integer ntcon

        real(dp) stress(3,3), cstress(3,3)
        real(dp) fa(3,na_qm+na_mm), cfa(3,na_qm+na_mm)
        integer i

        stress=0.0_dp
        fa=1.0_dp

        call fixed( cell, stress, na_qm+na_mm, isa, amass, xa, fa,
     .                  cstress, cfa, ntcon )

c set total number of free atoms
	do i=1,na_qm
	   if(blocklist(i).eq.1) then
	      if (cfa(1,i)>0.0_dp) ntcon=ntcon+1
	      if (cfa(2,i)>0.0_dp) ntcon=ntcon+1
	      if (cfa(3,i)>0.0_dp) ntcon=ntcon+1
	   endif
	enddo
	do i=1,na_mm
	   if((blocklist(i+na_qm).eq.1).or.(blockqmmm(i).eq.1)) then
	      if (cfa(1,i)>0.0_dp) ntcon=ntcon+1
	      if (cfa(2,i)>0.0_dp) ntcon=ntcon+1
	      if (cfa(3,i)>0.0_dp) ntcon=ntcon+1
	   endif
	enddo

        end

c******************************************************************

