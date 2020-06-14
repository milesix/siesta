c subroutine that calculates the free energy using
c the multiple steering molecualr dynamics approach 
c developed by Jarzynski 
	subroutine subfree(natot,rclas,fdummy,dt,freelog,totcoor,cell,
     .                lattice_type)

        use precision, only: dp
	use fdf  
        use m_qmmm_fdf, only : fdf_block_qmmm   
	use sys
	implicit none
	integer natot,i,j,k,l,constrats(20),ndists,totcoor
	integer at1,at2,at3,at4,at5,at6,at7,at8
	integer iunit,freetype,istp,at9,at10,at11,at12
	real(dp) kf,rclas(3,natot),work,fold,
     .  fdummy(3,natot),fnew(3,8),zo,zot,zt,uzt,v,dt,
     .  r12,r34,r56,r78,rtot,dist,dx,dy,dz,fce,rp(3)
        real(dp) dx12, dx32, dy12, dy32, dz12, dz32
        real(dp) dx23, dy23, dz23, dx43, dy43, dz43
        real(dp) dx34, dy34, dz34, dx56, dy56, dz56
        real(dp) dx78, dy78, dz78
        real(dp) pi,fdihe(12),dihedro2_v2
        real(dp) angle,scal,scalar_v2,dscalar,
     .  r32,dr12r32,coef(10)
	character exp
        character slabel*30, paste*30
        external  io_assign, io_close, paste
        character fname*30,fnameumb*30
        logical   frstme,found,freelog
        integer   unit
        real(dp) cell(3,3)
        character lattice_type
        real(dp) amber_cell(3,3), amber_kcell(3,3)
        real(dp) dist_v2, angle_v2
        save      frstme,unit,fname,zo,kf,v,freetype,constrats,istp
        save      fnameumb,fold,work,ndists,coef
        data      frstme /.true./
	data      found /.false./ 

c 8 ways of defining the reaction coordinate
c 1 = r1 - r2 coupled
c 2 = distance
c 3 = angle
c 4 = dihedral
c 5 = r1 + r2 coupled
c 6 = ( r1 + r2 ) - ( r3 + r4 ) coupled
c 7 = plane to atom distance
c 8 = c1*r1 + c2*r2 + c3*r3 + .... 

c change units every step
        rclas(1:3,1:natot)=rclas(1:3,1:natot)*0.529177
        amber_cell(1:3,1:3)=cell(1:3,1:3)*0.529177
        call reccel(3,amber_cell,amber_kcell,0)
        pi = acos(-1.0d0)
        fce = 0.0

c reading varables for first time
	if ( frstme ) then
	istp=0
	work=0.0
	fold=0.0
	if ( fdf_block_qmmm('FreeEnergy',iunit) ) then
        read(iunit,'(A)',advance='no',err=100,end=100) exp
        if(exp.eq.'%') then
        freelog=.false.
        return
        endif
        read(iunit,*,err=100,end=100) exp,freetype
        read(iunit,*,err=100,end=100) exp,kf
        read(iunit,*,err=100,end=100) exp,v
         if     (freetype.eq.1) then
         read(iunit,*,err=100,end=100) exp,(constrats(i),i=1,4)
         elseif (freetype.eq.2) then
         read(iunit,*,err=100,end=100) exp,(constrats(i),i=1,2)
         elseif (freetype.eq.3) then
         read(iunit,*,err=100,end=100) exp,(constrats(i),i=1,3)
         elseif (freetype.eq.4) then
         read(iunit,*,err=100,end=100) exp,(constrats(i),i=1,4)
         elseif (freetype.eq.5) then
         read(iunit,*,err=100,end=100) exp,(constrats(i),i=1,4)
         elseif (freetype.eq.6) then
         read(iunit,*,err=100,end=100) exp,(constrats(i),i=1,8)
         elseif (freetype.eq.7) then
         read(iunit,*,err=100,end=100) exp,(constrats(i),i=1,5)
         elseif (freetype.eq.8) then
         read(iunit,*,err=100,end=100) exp,ndists
       if(ndists.gt.10) then
       call die('constr opt: ndists in typeconstr 8 must not exceed 10')
       endif
         read(iunit,*,err=100,end=100) exp,(coef(i),i=1,ndists)
         read(iunit,*,err=100,end=100) exp,(constrats(i),i=1,ndists*2)
         else
         call die('free energy: freetype must be 1-8')
         endif
        if(i.gt.20) then
        call die('free energy: atoms in block must be lower than 20')
        endif
	else
	call die('free energy: You must specify the free energy block')  
	endif

        write(*,'(/,a)')
     .  'free energy: Starting a free energy calculation'

c zo calculaton if first time
	if(freetype.eq.1) then
	at1=constrats(1)
	at2=constrats(2)
	at3=constrats(3)
	at4=constrats(4)  

        dx12=rclas(1,at1)-rclas(1,at2)
        dy12=rclas(2,at1)-rclas(2,at2)
        dz12=rclas(3,at1)-rclas(3,at2)
        call pbc_displ_vector(lattice_type,amber_cell,amber_kcell,
     .                         dx12,dy12,dz12)
        r12=dist_v2(dx12,dy12,dz12)
        dx34=rclas(1,at3)-rclas(1,at4)
        dy34=rclas(2,at3)-rclas(2,at4)
        dz34=rclas(3,at3)-rclas(3,at4)
        call pbc_displ_vector(lattice_type,amber_cell,amber_kcell,
     .                         dx34,dy34,dz34)
        r34=dist_v2(dx34,dy34,dz34)

	zo=r34-r12

	elseif(freetype.eq.2) then

        at1=constrats(1)
        at2=constrats(2)

        dx12=rclas(1,at1)-rclas(1,at2)
        dy12=rclas(2,at1)-rclas(2,at2)
        dz12=rclas(3,at1)-rclas(3,at2)
        call pbc_displ_vector(lattice_type,amber_cell,amber_kcell,
     .                         dx12,dy12,dz12)
        r12=dist_v2(dx12,dy12,dz12)

        zo=r12      

        elseif(freetype.eq.3) then
        at1=constrats(1)
        at2=constrats(2)
        at3=constrats(3)

	dx12=rclas(1,at1)-rclas(1,at2)
	dy12=rclas(2,at1)-rclas(2,at2)
	dz12=rclas(3,at1)-rclas(3,at2)
        call pbc_displ_vector(lattice_type,amber_cell,amber_kcell,dx12,
     .                              dy12,dz12)
	dx32=rclas(1,at3)-rclas(1,at2)
	dy32=rclas(2,at3)-rclas(2,at2)
	dz32=rclas(3,at3)-rclas(3,at2)
        call pbc_displ_vector(lattice_type,amber_cell,amber_kcell,dx32,
     .                              dy32,dz32)
        rtot =angle_v2(dx12,dy12,dz12,dx32,dy32,dz32)

        zo=rtot

        elseif(freetype.eq.4) then
        at1=constrats(1)
        at2=constrats(2)
        at3=constrats(3)
        at4=constrats(4)

        rtot=dihedro2_v2(rclas(1,at1),rclas(2,at1),rclas(3,at1),
     .   rclas(1,at2),rclas(2,at2),rclas(3,at2),rclas(1,at3),
     .   rclas(2,at3),rclas(3,at3),rclas(1,at4),rclas(2,at4),
     .   rclas(3,at4),cell,lattice_type)

        zo=rtot

        elseif (freetype.eq.5) then
        at1=constrats(1)
        at2=constrats(2)
        at3=constrats(3)
        at4=constrats(4)

	dx12=rclas(1,at1)-rclas(1,at2)
	dy12=rclas(2,at1)-rclas(2,at2)
	dz12=rclas(3,at1)-rclas(3,at2)
        call pbc_displ_vector(lattice_type,amber_cell,amber_kcell,dx12,
     .                              dy12,dz12)
        r12=dist_v2(dx12,dy12,dz12)

	dx34=rclas(1,at3)-rclas(1,at4)
	dy34=rclas(2,at3)-rclas(2,at4)
	dz34=rclas(3,at3)-rclas(3,at4)
        call pbc_displ_vector(lattice_type,amber_cell,amber_kcell,dx34,
     .                              dy34,dz34)
        r34=dist_v2(dx34,dy34,dz34)

        rtot=r34+r12

        zo=rtot

        elseif (freetype.eq.6) then
        at1=constrats(1)
        at2=constrats(2)
        at3=constrats(3)
        at4=constrats(4)
        at5=constrats(5)
        at6=constrats(6)
        at7=constrats(7)
        at8=constrats(8)

	dx12=rclas(1,at1)-rclas(1,at2)
	dy12=rclas(2,at1)-rclas(2,at2)
	dz12=rclas(3,at1)-rclas(3,at2)
        call pbc_displ_vector(lattice_type,amber_cell,amber_kcell,dx12,
     .                              dy12,dz12)
        r12=dist_v2(dx12,dy12,dz12)

	dx34=rclas(1,at3)-rclas(1,at4)
	dy34=rclas(2,at3)-rclas(2,at4)
	dz34=rclas(3,at4)-rclas(3,at4)
        call pbc_displ_vector(lattice_type,amber_cell,amber_kcell,dx34,
     .                              dy34,dz34)
        r34=dist_v2(dx34,dy34,dz34)

	dx56=rclas(1,at5)-rclas(1,at6)
	dy56=rclas(2,at5)-rclas(2,at6)
	dz56=rclas(3,at5)-rclas(3,at6)
        call pbc_displ_vector(lattice_type,amber_cell,amber_kcell,dx56,
     .                              dy56,dz56)
        r56=dist_v2(dx56,dy56,dz56)

	dx78=rclas(1,at7)-rclas(1,at8)
	dy78=rclas(2,at7)-rclas(2,at8)
	dz78=rclas(3,at7)-rclas(3,at8)
        call pbc_displ_vector(lattice_type,amber_cell,amber_kcell,dx78,
     .                              dy78,dz78)
        r78=dist_v2(dx78,dy78,dz78)

        rtot=(r34+r12)-(r56+r78)

        zo=rtot

        elseif (freetype.eq.7) then
        at1=constrats(1)
        at2=constrats(2)
        at3=constrats(3)
        at4=constrats(4)
        at5=constrats(5)

        rp(1)=(rclas(1,at2)+rclas(1,at3)+rclas(1,at4)+rclas(1,at5))/4.0
        rp(2)=(rclas(2,at2)+rclas(2,at3)+rclas(2,at4)+rclas(2,at5))/4.0
        rp(3)=(rclas(3,at2)+rclas(3,at3)+rclas(3,at4)+rclas(3,at5))/4.0

        dx12=rclas(1,at1)-rp(1)
        dy12=rclas(2,at1)-rp(2)
        dz12=rclas(3,at1)-rp(3)
        call pbc_displ_vector(lattice_type,amber_cell,amber_kcell,dx12,
     .                              dy12,dz12)
        rtot=dist_v2(dx12,dy12,dz12)

        zo=rtot

        elseif (freetype.eq.8) then
        rtot=0.0
        do i=1,ndists
        at1=constrats(i)
        at2=constrats(i+1)

        dx12=rclas(1,at1)-rclas(1,at2)
        dy12=rclas(2,at1)-rclas(2,at2)
        dz12=rclas(3,at1)-rclas(3,at2)
        call pbc_displ_vector(lattice_type,amber_cell,amber_kcell,dx12,
     .                              dy12,dz12)
        r12=dist_v2(dx12,dy12,dz12)

        rtot=rtot+coef(i)*r12
c ndists
        enddo
        zo=rtot

c freetype
        endif

c reads from .fce of a former run 
        slabel = fdf_string( 'SystemLabel', 'siesta' )
        fname = paste(slabel,'.fce')
        fnameumb = paste(slabel,'.umb')        

	inquire( file=fname, exist=found )   
	if(found) then

	call io_assign(unit)
        open( unit, file=fname ) 

	zo=0.0
	k=0
 10     continue
	k=k+1
        read(unit,*,err=200,end=20) zo,fold,work
        goto 10
 20     continue
        call io_close(unit)
        istp=k
        zo=zo-v*dt*istp
	istp=istp+1

        write(*,'(/,a)')
     .  'free energy: re-starting a free energy calculation from file'
	endif !found

C write in output file free energy values only the first step
        write(*,'(/,a)')     'free energy: calculation values'
        write(*,'(A,i4)')    'type       :', freetype 
        write(*,'(A,F8.3)')  'constant   :', kf
        write(*,'(A,F8.3)')  'velocity   :', v
        write(*,'(A,i7)')    'steps      :', totcoor
        write(*,'(A,F8.3)')  'initial    :', zo
        write(*,'(A,F8.3)')  'final      :', zo+v*dt*(totcoor-istp)

        if     (freetype.eq.1) then
        write(*,'(A,4i4)')   'atoms      :', (constrats(i),i=1,4)
        elseif (freetype.eq.2) then
        write(*,'(A,2i4)')   'atoms      :', (constrats(i),i=1,2)
        elseif (freetype.eq.3) then
        write(*,'(A,3i4)')   'atoms      :', (constrats(i),i=1,3)
        elseif (freetype.eq.4) then
        write(*,'(A,4i4)')   'atoms      :', (constrats(i),i=1,4)
        elseif (freetype.eq.5) then
        write(*,'(A,4i4)')   'atoms      :', (constrats(i),i=1,4)
        elseif (freetype.eq.6) then
        write(*,'(A,8i4)')   'atoms      :', (constrats(i),i=1,8)
        elseif (freetype.eq.7) then
        write(*,'(A,5i4)')   'atoms      :', (constrats(i),i=1,5)
        elseif (freetype.eq.8) then
        write(*,'(A,i4)')    'ndists     :', ndists
        do i=1,ndists
        write(*,'(A,F5.2)')  'coefs      :', coef(i)
        enddo
        do i=1,ndists*2
        write(*,'(A,i4)')    'atoms      :', constrats(i)
        enddo
        endif

        endif !frstme

c zt time calculation 
        if(freetype.eq.1) then    
        at1=constrats(1)
        at2=constrats(2)
        at3=constrats(3)
        at4=constrats(4)
 
        dx12=rclas(1,at1)-rclas(1,at2)
        dy12=rclas(2,at1)-rclas(2,at2)
        dz12=rclas(3,at1)-rclas(3,at2)
        call pbc_displ_vector(lattice_type,amber_cell,amber_kcell,
     .                         dx12,dy12,dz12)
        r12=dist_v2(dx12,dy12,dz12)

        dx34=rclas(1,at3)-rclas(1,at4)
        dy34=rclas(2,at3)-rclas(2,at4)
        dz34=rclas(3,at3)-rclas(3,at4)
        call pbc_displ_vector(lattice_type,amber_cell,amber_kcell,
     .                         dx34,dy34,dz34)
        r34=dist_v2(dx34,dy34,dz34)

        zt=r34-r12

c u(z,t)
        zot=zo+v*dt*istp
        uzt=(kf/2.)*(zt-zot)**2

c forces
          fce=-kf*(zt-zot)
          work=work+(fce+fold)*(v*dt/2.)
          fold=fce
 
c atom1: -dr12       
        dx=0.0
        dy=0.0
        dz=0.0
          dx=(1.0/r12)*dx12
          dy=(1.0/r12)*dy12
          dz=(1.0/r12)*dz12

        fnew(1,1)=-dx*fce
        fnew(2,1)=-dy*fce
        fnew(3,1)=-dz*fce

c atom3: dr34
        dx=0.0
        dy=0.0
        dz=0.0
          dx=(1.0/r34)*dx34
          dy=(1.0/r34)*dy34
          dz=(1.0/r34)*dz34

        fnew(1,3)=dx*fce
        fnew(2,3)=dy*fce
        fnew(3,3)=dz*fce

c atom2: dr12
        fnew(1,2)=-fnew(1,1)
        fnew(2,2)=-fnew(2,1)
        fnew(3,2)=-fnew(3,1)

c atom4: -dr34
        fnew(1,4)=-fnew(1,3)
        fnew(2,4)=-fnew(2,3)
        fnew(3,4)=-fnew(3,3)

c add fnew to fdummy 
        fdummy(1:3,at1)= fdummy(1:3,at1)+fnew(1:3,1)
        fdummy(1:3,at2)= fdummy(1:3,at2)+fnew(1:3,2)
        fdummy(1:3,at3)= fdummy(1:3,at3)+fnew(1:3,3)
        fdummy(1:3,at4)= fdummy(1:3,at4)+fnew(1:3,4)      

        elseif(freetype.eq.2) then    
        at1=constrats(1)
        at2=constrats(2)
 
        dx12=rclas(1,at1)-rclas(1,at2)
        dy12=rclas(2,at1)-rclas(2,at2)
        dz12=rclas(3,at1)-rclas(3,at2) 
        call pbc_displ_vector(lattice_type,amber_cell,amber_kcell,
     .                              dx12,dy12,dz12)
        r12=dist_v2(dx12,dy12,dz12)
 
        zt=r12
 
c u(z,t)
        zot=zo+v*dt*istp
        uzt=(kf/2.)*(zt-zot)**2

c forces 
          fce=-kf*(zt-zot)
          work=work+(fce+fold)*(v*dt/2.)
          fold=fce

c atom1: -dr12
        dx=0.0
        dy=0.0
        dz=0.0
          dx=(1.0/r12)*dx12
          dy=(1.0/r12)*dy12
          dz=(1.0/r12)*dz12
 
        fnew(1,1)=dx*fce
        fnew(2,1)=dy*fce
        fnew(3,1)=dz*fce

c atom2: dr12
        fnew(1,2)=-fnew(1,1)
        fnew(2,2)=-fnew(2,1)
        fnew(3,2)=-fnew(3,1)
 
c add fnew to fdummy 
        fdummy(1:3,at1)= fdummy(1:3,at1)+fnew(1:3,1)
        fdummy(1:3,at2)= fdummy(1:3,at2)+fnew(1:3,2)

        elseif(freetype.eq.3) then
        at1=constrats(1)
        at2=constrats(2)
        at3=constrats(3)

	dx12=rclas(1,at1)-rclas(1,at2)
	dy12=rclas(2,at1)-rclas(2,at2)
	dz12=rclas(3,at1)-rclas(3,at2)
        call pbc_displ_vector(lattice_type,amber_cell,amber_kcell,dx12,
     .                              dy12,dz12)
	dx32=rclas(1,at3)-rclas(1,at2)
	dy32=rclas(2,at3)-rclas(2,at2)
	dz32=rclas(3,at3)-rclas(3,at2)
        call pbc_displ_vector(lattice_type,amber_cell,amber_kcell,dx32,
     .                              dy32,dz32)

        rtot =angle_v2(dx12,dy12,dz12,dx32,dy32,dz32)

        zt=rtot

c u(z,t)
        zot=zo+v*dt*istp
        uzt=(pi/180.)**2*(kf/2.)*(zt-zot)**2

c forces
       fce=-(pi/180.)*kf*(zt-zot)

	scal=scalar_v2(dx12,dy12,dz12,dx32,dy32,dz32)
	r12=dist_v2(dx12,dy12,dz12)
	r32=dist_v2(dx32,dy32,dz32)

c atom1:
       dr12r32=r32*dx12/(r12)
       dx=(dx32*r12*r32-scal*dr12r32)/(r12*r32)**(2.0)
       dx=-1.0/(sqrt(1.0-(scal/(r12*r32))**2.0))*dx
       dx=fce*dx
       dr12r32=r32*dy12/(r12)
       dy=(dy32*r12*r32-scal*dr12r32)/(r12*r32)**(2.0)
       dy=-1.0/(sqrt(1.0-(scal/(r12*r32))**2.0))*dy
       dy=fce*dy
       dr12r32=r32*dz12/(r12)
       dz=(dz32*r12*r32-scal*dr12r32)/(r12*r32)**(2.0)
       dz=-1.0/(sqrt(1.0-(scal/(r12*r32))**2.0))*dz
       dz=fce*dz
       fnew(1,1)=-dx
       fnew(2,1)=-dy
       fnew(3,1)=-dz

c atom2:
       dscalar=-dx12-dx32
       dr12r32=(r32*(-dx12)/r12)+(r12*(-dx32)/r32)
       dx=(dscalar*r12*r32-scal*dr12r32)/(r12*r32)**(2.0)
       dx=-1.0/(sqrt(1.0-(scal/(r12*r32))**2.0))*dx
       dx=fce*dx
       dscalar=-dy12-dy32
       dr12r32=(r32*(-dy12)/r12)+(r12*(-dy32)/r32)
       dy=(dscalar*r12*r32-scal*dr12r32)/(r12*r32)**(2.0)
       dy=-1.0/(sqrt(1.0-(scal/(r12*r32))**2.0))*dy
       dy=fce*dy
       dscalar=-dz12-dz32
       dr12r32=(r32*(-dz12)/r12)+(r12*(-dz32)/r32)
       dz=(dscalar*r12*r32-scal*dr12r32)/(r12*r32)**(2.0)
       dz=-1.0/(sqrt(1.0-(scal/(r12*r32))**2.0))*dz
       dz=fce*dz
       fnew(1,2)=-dx
       fnew(2,2)=-dy
       fnew(3,2)=-dz

c atom3:
       fnew(1,3)=-fnew(1,1)
       fnew(2,3)=-fnew(2,1)
       fnew(3,3)=-fnew(3,1)

c adding fnew to fdummy
        fdummy(1:3,at1)= fdummy(1:3,at1)+fnew(1:3,1)
        fdummy(1:3,at2)= fdummy(1:3,at2)+fnew(1:3,2)
        fdummy(1:3,at3)= fdummy(1:3,at3)+fnew(1:3,3)

c accumulated work
       fce=fce*(pi/180.)
       work=work+(fce+fold)*(v*dt/2.)
       fold=fce

        elseif(freetype.eq.4) then
        at1=constrats(1)
        at2=constrats(2)
        at3=constrats(3)
        at4=constrats(4)

        rtot=dihedro2_v2(rclas(1,at1),rclas(2,at1),rclas(3,at1),
     .   rclas(1,at2),rclas(2,at2),
     .   rclas(3,at2),
     .   rclas(1,at3),rclas(2,at3),
     .   rclas(3,at3),
     .   rclas(1,at4),rclas(2,at4),
     .   rclas(3,at4),amber_cell,amber_kcell,
     .   lattice_type)

        zt=rtot

c u(z,t)
        zot=zo+v*dt*istp
        uzt=(pi/180.)**2*(kf/2.)*(zt-zot)**2

c forces
       if(zot.lt.90 .and.zt.gt.180) zt=zt-360.
       if(zot.gt.270.and.zt.lt.180) zt=zt+360.
       fce=-(pi/180.)*kf*(zt-zot)
        call diheforce2_v2(natot,rclas,at1,at2,at3,at4,
     .          1,fce,fdihe,amber_cell,amber_kcell,
     .   lattice_type)
        fnew(1:3,1)=fdihe(1:3)
        call diheforce2_v2(natot,rclas,at1,at2,at3,at4,
     .          2,fce,fdihe,amber_cell,amber_kcell,
     .   lattice_type)
        fnew(1:3,2)=fdihe(4:6)
        call diheforce2_v2(natot,rclas,at1,at2,at3,at4,
     .          3,fce,fdihe,amber_cell,amber_kcell,
     .   lattice_type)
        fnew(1:3,3)=fdihe(7:9)
        call diheforce2_v2(natot,rclas,at1,at2,at3,at4,
     .          4,fce,fdihe,amber_cell,amber_kcell,
     .   lattice_type)
        fnew(1:3,4)=fdihe(10:12)

c adding fnew to fdummy
        if((zt.ge.0..and.zt.le.180.).or.(zt.gt.360)) then
         fnew(1:3,1:4)=(-1.0)*fnew(1:3,1:4)
        elseif((zt.gt.180..and.zt.lt.360).or.(zt.lt.0)) then
         fnew(1:3,1:4)=fnew(1:3,1:4)
        else
         call die('free energy: Wrong dihedral angle value')
        endif
        fdummy(1:3,at1)= fdummy(1:3,at1)+fnew(1:3,1)
        fdummy(1:3,at2)= fdummy(1:3,at2)+fnew(1:3,2)
        fdummy(1:3,at3)= fdummy(1:3,at3)+fnew(1:3,3)
        fdummy(1:3,at4)= fdummy(1:3,at4)+fnew(1:3,4)

c accumulated work
       fce=fce*(pi/180.)
       work=work+(fce+fold)*(v*dt/2.)
       fold=fce

        elseif(freetype.eq.5) then
        at1=constrats(1)
        at2=constrats(2)
        at3=constrats(3)
        at4=constrats(4)

	dx12=rclas(1,at1)-rclas(1,at2)
	dy12=rclas(2,at1)-rclas(2,at2)
	dz12=rclas(3,at1)-rclas(3,at2)
        call pbc_displ_vector(lattice_type,amber_cell,amber_kcell,dx12,
     .                              dy12,dz12)
        r12=dist_v2(dx12,dy12,dz12)

	dx34=rclas(1,at3)-rclas(1,at4)
	dy34=rclas(2,at3)-rclas(2,at4)
	dz34=rclas(3,at3)-rclas(3,at4)
        call pbc_displ_vector(lattice_type,amber_cell,amber_kcell,dx34,
     .                              dy34,dz34)
        r34=dist_v2(dx34,dy34,dz34)

        zt=r34+r12

c u(z,t)
        zot=zo+v*dt*istp
        uzt=(kf/2.)*(zt-zot)**2

c forces
          fce=-kf*(zt-zot)
          work=work+(fce+fold)*(v*dt/2.)
          fold=fce

c atom1: dr12
        dx=0.0
        dy=0.0
        dz=0.0
          dx=(1.0/r12)*dx12
          dy=(1.0/r12)*dy12
          dz=(1.0/r12)*dz12

        fnew(1,1)=dx*fce
        fnew(2,1)=dy*fce
        fnew(3,1)=dz*fce

c atom3: dr34
        dx=0.0
        dy=0.0
        dz=0.0
          dx=(1.0/r34)*dx34
          dy=(1.0/r34)*dy34
          dz=(1.0/r34)*dz34

        fnew(1,3)=dx*fce
        fnew(2,3)=dy*fce
        fnew(3,3)=dz*fce

c atom2: -dr12
        fnew(1,2)=-fnew(1,1)
        fnew(2,2)=-fnew(2,1)
        fnew(3,2)=-fnew(3,1)

c atom4: -dr34
        fnew(1,4)=-fnew(1,3)
        fnew(2,4)=-fnew(2,3)
        fnew(3,4)=-fnew(3,3)

c add fnew to fdummy
        fdummy(1:3,at1)= fdummy(1:3,at1)+fnew(1:3,1)
        fdummy(1:3,at2)= fdummy(1:3,at2)+fnew(1:3,2)
        fdummy(1:3,at3)= fdummy(1:3,at3)+fnew(1:3,3)
        fdummy(1:3,at4)= fdummy(1:3,at4)+fnew(1:3,4)

        elseif(freetype.eq.6) then
        at1=constrats(1)
        at2=constrats(2)
        at3=constrats(3)
        at4=constrats(4)
        at5=constrats(5)
        at6=constrats(6)
        at7=constrats(7)
        at8=constrats(8)

	dx12=rclas(1,at1)-rclas(1,at2)
	dy12=rclas(2,at1)-rclas(2,at2)
	dz12=rclas(3,at1)-rclas(3,at2)
        call pbc_displ_vector(lattice_type,amber_cell,amber_kcell,dx12,
     .                              dy12,dz12)
        r12=dist_v2(dx12,dy12,dz12)

	dx34=rclas(1,at3)-rclas(1,at4)
	dy34=rclas(2,at3)-rclas(2,at4)
	dz34=rclas(3,at4)-rclas(3,at4)
        call pbc_displ_vector(lattice_type,amber_cell,amber_kcell,dx34,
     .                              dy34,dz34)
        r34=dist_v2(dx34,dy34,dz34)

	dx56=rclas(1,at5)-rclas(1,at6)
	dy56=rclas(2,at5)-rclas(2,at6)
	dz56=rclas(3,at5)-rclas(3,at6)
        call pbc_displ_vector(lattice_type,amber_cell,amber_kcell,dx56,
     .                              dy56,dz56)
        r56=dist_v2(dx56,dy56,dz56)

	dx78=rclas(1,at7)-rclas(1,at8)
	dy78=rclas(2,at7)-rclas(2,at8)
	dz78=rclas(3,at7)-rclas(3,at8)
        call pbc_displ_vector(lattice_type,amber_cell,amber_kcell,dx78,
     .                              dy78,dz78)
        r78=dist_v2(dx78,dy78,dz78)

        zt=(r34+r12)-(r56+r78)

c u(z,t)
        zot=zo+v*dt*istp
        uzt=(kf/2.)*(zt-zot)**2

c forces
          fce=-kf*(zt-zot)
          work=work+(fce+fold)*(v*dt/2.)
          fold=fce

c atom1: dr12
        dx=0.0
        dy=0.0
        dz=0.0
          dx=(1.0/r12)*dx12
          dy=(1.0/r12)*dy12
          dz=(1.0/r12)*dz12

        fnew(1,1)=dx*fce
        fnew(2,1)=dy*fce
        fnew(3,1)=dz*fce

c atom3: dr34
        dx=0.0
        dy=0.0
        dz=0.0
          dx=(1.0/r34)*dx34
          dy=(1.0/r34)*dy34
          dz=(1.0/r34)*dz34

        fnew(1,3)=dx*fce
        fnew(2,3)=dy*fce
        fnew(3,3)=dz*fce

c atom2: dr12
        fnew(1,2)=-fnew(1,1)
        fnew(2,2)=-fnew(2,1)
        fnew(3,2)=-fnew(3,1)

c atom4: -dr34
        fnew(1,4)=-fnew(1,3)
        fnew(2,4)=-fnew(2,3)
        fnew(3,4)=-fnew(3,3)

c atom5: dr56
        dx=0.0
        dy=0.0
        dz=0.0
          dx=(1.0/r56)*dx56
          dy=(1.0/r56)*dy56
          dz=(1.0/r56)*dz56

        fnew(1,5)=-dx*fce
        fnew(2,5)=-dy*fce
        fnew(3,5)=-dz*fce

c atom7: dr78
        dx=0.0
        dy=0.0
        dz=0.0
          dx=(1.0/r78)*dx78
          dy=(1.0/r78)*dy78
          dz=(1.0/r78)*dz78

        fnew(1,7)=-dx*fce
        fnew(2,7)=-dy*fce
        fnew(3,7)=-dz*fce

c atom6: -dr56
        fnew(1,6)=-fnew(1,5)
        fnew(2,6)=-fnew(2,5)
        fnew(3,6)=-fnew(3,5)

c atom8: -dr78
        fnew(1,8)=-fnew(1,7)
        fnew(2,8)=-fnew(2,7)
        fnew(3,8)=-fnew(3,7)

c adding fnew to fdummy
        fdummy(1:3,at1)= fdummy(1:3,at1)+fnew(1:3,1)
        fdummy(1:3,at2)= fdummy(1:3,at2)+fnew(1:3,2)
        fdummy(1:3,at3)= fdummy(1:3,at3)+fnew(1:3,3)
        fdummy(1:3,at4)= fdummy(1:3,at4)+fnew(1:3,4)
        fdummy(1:3,at5)= fdummy(1:3,at5)+fnew(1:3,5)
        fdummy(1:3,at6)= fdummy(1:3,at6)+fnew(1:3,6)
        fdummy(1:3,at7)= fdummy(1:3,at7)+fnew(1:3,7)
        fdummy(1:3,at8)= fdummy(1:3,at8)+fnew(1:3,8)

        elseif(freetype.eq.7) then
        at1=constrats(1)
        at2=constrats(2)
        at3=constrats(3)
        at4=constrats(4)
        at5=constrats(5)

        rp(1)=(rclas(1,at2)+rclas(1,at3)+rclas(1,at4)+rclas(1,at5))/4.0
        rp(2)=(rclas(2,at2)+rclas(2,at3)+rclas(2,at4)+rclas(2,at5))/4.0
        rp(3)=(rclas(3,at2)+rclas(3,at3)+rclas(3,at4)+rclas(3,at5))/4.0

        dx12=rclas(1,at1)-rp(1)
        dy12=rclas(2,at1)-rp(2)
        dz12=rclas(3,at1)-rp(3)
        call pbc_displ_vector(lattice_type,amber_cell,amber_kcell,dx12,
     .                              dy12,dz12)
        rtot=dist_v2(dx12,dy12,dz12)
        
        zt=rtot 

c u(z,t)
        zot=zo+v*dt*istp
        uzt=(kf/2.)*(zt-zot)**2

c forces
          fce=-kf*(zt-zot)
          work=work+(fce+fold)*(v*dt/2.)
          fold=fce

c atom1: dr12
          dx=0.0
          dy=0.0
          dz=0.0
          dx=(1.0/rtot)*dx12
          dy=(1.0/rtot)*dy12
          dz=(1.0/rtot)*dz12

        fnew(1,1)=dx*fce
        fnew(2,1)=dy*fce
        fnew(3,1)=dz*fce

c adding fnew to fdummy
        fdummy(1:3,at1)= fdummy(1:3,at1)+fnew(1:3,1)

        elseif (freetype.eq.8) then
        rtot=0.0
        do i=1,ndists
        at1=constrats(i)
        at2=constrats(i+1)

        dx12=rclas(1,at1)-rclas(1,at2)
        dy12=rclas(2,at1)-rclas(2,at2)
        dz12=rclas(3,at1)-rclas(3,at2)
        call pbc_displ_vector(lattice_type,amber_cell,amber_kcell,dx12,
     .                              dy12,dz12)
        r12=dist_v2(dx12,dy12,dz12)

        rtot=rtot+coef(i)*r12
c ndists
        enddo
        zt=rtot

c u(z,t)
        zot=zo+v*dt*istp
        uzt=(kf/2.)*(zt-zot)**2

c forces
          fce=-kf*(zt-zot)
          work=work+(fce+fold)*(v*dt/2.)
          fold=fce

        do i=1,ndists
        at1=constrats(i)
        at2=constrats(i+1)

        dx12=rclas(1,at1)-rclas(1,at2)
        dy12=rclas(2,at1)-rclas(2,at2)
        dz12=rclas(3,at1)-rclas(3,at2)
        call pbc_displ_vector(lattice_type,amber_cell,amber_kcell,dx12,
     .                              dy12,dz12)

c atom1: -dr12
        dx=0.0
        dy=0.0
        dz=0.0
          dx=(1.0/rtot)*dx12
          dy=(1.0/rtot)*dy12
          dz=(1.0/rtot)*dz12

        fnew(1,1)=dx*fce*coef(i)
        fnew(2,1)=dy*fce*coef(i)
        fnew(3,1)=dz*fce*coef(i)

c atom2: dr12
        fnew(1,2)=-fnew(1,1)
        fnew(2,2)=-fnew(2,1)
        fnew(3,2)=-fnew(3,1)

c add fnew to fdummy
        fdummy(1:3,at1)= fdummy(1:3,at1)+fnew(1:3,1)
        fdummy(1:3,at2)= fdummy(1:3,at2)+fnew(1:3,2)

c ndists
        enddo

c freetype
        endif

c writing .FCE or .UMB with RC and FCE
      call io_assign(unit)
      open( unit, file=fname, form = 'formatted', position='append',
     .      status='unknown')
      write(unit,'(F10.6,2x,F14.7,2x,F14.7)') zot,fce,work 
      call io_close(unit)

      if(v.eq.0.0) then
      call io_assign(unit)
      open( unit, file=fnameumb, form = 'formatted', position='append',
     .      status='unknown')
      write(unit,'(F7.4,2x,F10.6)')  dt*istp,zt
      call io_close(unit)
      endif

c change units 
        rclas(1:3,1:natot)=rclas(1:3,1:natot)/0.529177

	istp=istp+1
	frstme=.false.
	return
 100    stop 'read: problem reading FreeEnergy block'
 200    stop 'read: problem reading from .fce file' 
	end

