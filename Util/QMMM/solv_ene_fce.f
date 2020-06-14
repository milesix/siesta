c this subroutine calculates the solvent energy

      subroutine solv_ene_fce(natot,na_qm,na_mm,ng1,rclas,Em,Rm,pc,
     .    Etot_amber,fcetot_amber,stress_amber,attype,
     .    nbond,nangle,ndihe,nimp,multidihe, multiimp,kbond,bondeq,
     .    kangle,angleeq,kdihe,diheeq,kimp,impeq,perdihe,perimp,
     .    bondtype,angletype,dihetype,imptype,
     .    bondxat,angexat,angmxat,dihexat,dihmxat,impxat,atange,
     .    atangm,atdihe,atdihm,atimp,
     .    ng1type,angetype,angmtype,evaldihe,evaldihm,
     .    dihety,dihmty,impty,nonbonded,
     .    scale,nonbondedxat,scalexat,
     .    evaldihelog,evaldihmlog,paso,nparm,
     .    graphite_layer_no,
     .    actualiz,listcut,
     .    noat,noaa,sfc,timestep,
     .    water,masst,cell,lattice_type,coulombtype,
     .    ewald_alpha,kewald_cutoff)         

        use precision, only: dp

        implicit none
c      vbles grales
       integer   i,j,k,l,na_qm,natot,na_mm,ng1(na_mm,6),paso 
       real(dp)  pcA(na_mm),rclas(3,natot),ramber(3,na_mm),
     .    EmA(na_mm),RmA(na_mm),pc(na_mm),Em(natot),Rm(natot),
     .    Etot_amber,Ebond_amber,Eangle_amber,Edihe_amber,Eimp_amber,
     .    Elj_amber,Eelec_amber,Elj_amber14,Eelec_amber14
       real(dp) fcetot_amber(3,na_mm),
     .   fcebond_amber(3,na_mm),fceangle_amber(3,na_mm),
     .   fcedihe_amber(3,na_mm),fceimp_amber(3,na_mm),
     .   fcelj_amber(3,na_mm),fceelec_amber(3,na_mm),
     .   fcelj_amber14(3,na_mm),fceelec_amber14(3,na_mm),
     .   fcetotaxes_amber(3)
       character  attype(na_mm)*4,noat(na_mm)*4,noaa(na_mm)*4, 
     .                        coulombtype*10
       real(dp) listcut,sfc,timestep,masst(natot),ewat,fwat(3,na_mm)
       logical water
c      vbles de los params de union
       integer   nbond,nangle,ndihe,nimp,nparm,multidihe(nparm),
     .    multiimp(nparm)
       real(dp)   kbond(nparm),bondeq(nparm),kangle(nparm),
     .   angleeq(nparm),kdihe(nparm),diheeq(nparm),kimp(nparm),
     .   impeq(nparm),perdihe(nparm),perdihe2(nparm),perimp(nparm)
       character   bondtype(nparm)*5,angletype(nparm)*8,
     .   dihetype(nparm)*11,imptype(nparm)*11
c      vbles de bond, angle, dihe e imp
       integer   bondxat(na_mm),angexat(na_mm),angmxat(na_mm),
     .   dihexat(na_mm),dihmxat(na_mm),impxat(na_mm)
       integer   atange(na_mm,25,2),atangm(na_mm,25,2),
     .   atdihe(na_mm,100,3),atdihm(na_mm,100,3),atimp(na_mm,25,4)
c	vbles q faltaban
       integer  ng1type(na_mm,6),angetype(na_mm,25),angmtype(na_mm,25),
     .          evaldihe(na_mm,100,5),evaldihm(na_mm,100,5),
     .          dihety(na_mm,100),dihmty(na_mm,100),impty(na_mm,25),
     .          nonbonded(na_mm,100),scale(na_mm,100),
     .          nonbondedxat(na_mm),scalexat(na_mm)
       logical  evaldihelog(na_mm,100),evaldihmlog(na_mm,100),actualiz


       integer graphite_layer_no(na_mm)
       character lattice_type
       real(dp) cell(3,3), amber_cell(3,3), amber_kcell(3,3)

       real(dp) stress_amber(3,3)
       real(dp) ewald_alpha
       real(dp) kewald_cutoff

c inicializa las energias y fuerzas
      Etot_amber=0.0
      Ebond_amber=0.0
      Eangle_amber=0.0
      Edihe_amber=0.0
      Eimp_amber=0.0
      Elj_amber=0.0
      Eelec_amber=0.0
      Elj_amber14=0.0
      Eelec_amber14=0.0
      fcetot_amber=0.0
      fcebond_amber=0.0
      fceangle_amber=0.0
      fcedihe_amber=0.0
      fceimp_amber=0.0
      fcelj_amber=0.0
      fceelec_amber=0.0
      fcetotaxes_amber=0.0   
      ewat=0.0
      fwat=0.0
      stress_amber=0.0

c cambia variables
      k=1
      do j=1,na_mm
      pcA(k)=pc(j)
      k=k+1
      enddo
      k=1
      do i=na_qm+1,natot
      ramber(1:3,k)=rclas(1:3,i)
      k=k+1
      enddo
      k=1
      do i=1,na_mm
      EmA(k)=Em(na_qm+i)
      RmA(k)=Rm(na_qm+i)
      k=k+1
      enddo
 
c  pasa a las unidades de Amber
      do i=1,na_mm
      RmA(i) = RmA(i)*(2.0**(1.0/6.0))*0.529177/2.0
      EmA(i) = EmA(i)*627.5108
      ramber(1:3,i)=ramber(1:3,i)*0.529177
      enddo
      amber_cell(1:3,1:3)=cell(1:3,1:3)*0.529177
      call reccel(3,amber_cell,amber_kcell,0)

c  llama a subrutina q calcula la energia y fuerzas de bonds
        call amber_bonds(na_mm,ng1,ramber,Ebond_amber,attype,
     .       nbond,kbond,bondeq,bondtype,bondxat,fcebond_amber,
     .       stress_amber,ng1type,paso,nparm,amber_cell,amber_kcell,
     .       lattice_type)

c  llama a subrutina q calcula la energia y fuerzas de angles
        call amber_angles(na_mm,ng1,ramber,Eangle_amber,attype,
     .       nangle,kangle,angleeq,angletype,angexat,angmxat,atange,
     .       atangm,fceangle_amber,stress_amber,angetype,angmtype,paso,
     .       nparm,amber_cell,amber_kcell,lattice_type)

c  llama a subrutina q calcula la energia y fuerzas de dihedros     
	perdihe2=perdihe  
        call amber_dihes(na_mm,ng1,ramber,Edihe_amber,
     .            attype,ndihe,kdihe,diheeq,perdihe2,multidihe,
     .            dihetype,dihexat,dihmxat,atdihe,atdihm,
     .            fcedihe_amber,stress_amber,evaldihelog,evaldihe,
     .            dihety,evaldihmlog,evaldihm,dihmty,paso,nparm,
     .            amber_cell,amber_kcell,lattice_type) 

c  llama a subrutina q calcula la energia y fuerzas de impropers
        call amber_improper(na_mm,ng1,ramber,Eimp_amber,attype,
     .       nimp,kimp,impeq,perimp,multiimp,imptype,impxat,atimp,
     .       fceimp_amber,stress_amber,impty,paso,nparm,amber_cell,
     .       amber_kcell,lattice_type)

c  llama a subrutina q calcula la energia y fuerzas de terminos non-bonded
        call amber_nonbonded(na_mm,ng1,ramber,Elj_amber,
     .       Eelec_amber,Elj_amber14,Eelec_amber14,attype,
     .       EmA,RmA,pcA,bondxat,
     .       angexat,angmxat,atange,atangm,
     .       dihexat,dihmxat,atdihe,atdihm,
     .       fceelec_amber,fcelj_amber,stress_amber,nonbonded,
     .       coulombtype,scale,nonbondedxat,scalexat,paso,actualiz,
     .       listcut,
     .       noat,noaa,sfc,timestep,
     .       na_qm,natot,rclas,
     .       graphite_layer_no,water,masst,ewat,fwat,amber_cell,
     .       ewald_alpha,kewald_cutoff,amber_kcell,lattice_type)       

c  calcula la energia total de amber
       Etot_amber=Ebond_amber+Eangle_amber+Edihe_amber+Eimp_amber+
     .    Elj_amber+Eelec_amber+Elj_amber14+Eelec_amber14+ewat

c  calcula la fuerza total de amber
	do i=1,na_mm
	fcetot_amber(1:3,i)=fcebond_amber(1:3,i)+fceangle_amber(1:3,i)
     .  +fcedihe_amber(1:3,i)+fceimp_amber(1:3,i)+
     .  fcelj_amber(1:3,i)+fceelec_amber(1:3,i)+fwat(1:3,i)
       fcetot_amber(1:3,i)=(-1.0)*fcetot_amber(1:3,i)     
       fcetotaxes_amber(1:3)=fcetotaxes_amber(1:3)+fcetot_amber(1:3,i)   
	enddo

      return
      end
c
c****************************************************************
c subrutina q calcula la energia y fuerzas de bonds

        subroutine amber_bonds(na_mm,ng1,ramber,Ebond_amber,
     .             attype,nbond,kbond,bondeq,bondtype,bondxat,
     .             fcebond_amber,stress_amber,ng1type,paso,nparm,cell
     .             ,kcell,lattice_type)

       use precision, only: dp
       implicit none      
c      vbles grales 
       integer   na_mm,ng1(na_mm,6),i,j,k,l,m,n,paso
       real(dp)   ramber(3,na_mm),Ebond_amber,
     .                    fcebond_amber(3,na_mm)
       character   attype(na_mm)*4
c      vbles de los params de union
       integer   nbond,nparm
       real(dp)   kbond(nparm),bondeq(nparm)
       character   bondtype(nparm)*5
c      vbles de bond, angle, dihe e imp
       integer   bondxat(na_mm)
c      vbles de asignacion
       character*4 ty1,ty2
       character*5 tybond
       integer ng1type(na_mm,6)
       real(dp)   rij,dist_v2,dx1,dx2,dy1,dy2,dz1,dz2,
     .                    ftotbond(3)
       real(dp) volcel
       real(dp) cell(3,3), kcell(3,3)
       real(dp) stress_amber(3,3)
       real(dp) :: stress_fact
       real(dp), dimension(3) :: dr
       character lattice_type
       logical           first                                                                  
       save              first                                                                  
       data              first /.true./    

c     asignacion de tipos de union
      if(first) then
         do i=1,na_mm
            do j=1,bondxat(i)
               do k=1,nbond
                  tybond=bondtype(k)
                  ty1=tybond(1:2)
                  ty2=tybond(4:5)
                  if(attype(i).eq.ty1.and.attype(ng1(i,j)).eq.ty2) then
                     ng1type(i,j)=k
                  elseif(attype(i).eq.ty2.and.attype(ng1(i,j)).eq.ty1) 
     .                                     then
                     ng1type(i,j)=k
                  endif
               enddo
            enddo
         enddo
         first=.false.
      endif                     !asignacion

      stress_fact=1.0/volcel(cell)
c     calculo de la E y F de bond
      do i=1,na_mm 
         do j=1,bondxat(i) 
            dx1=ramber(1,i)-ramber(1,ng1(i,j))
            dy1=ramber(2,i)-ramber(2,ng1(i,j))
            dz1=ramber(3,i)-ramber(3,ng1(i,j))
            call pbc_displ_vector(lattice_type,cell,kcell,dx1,dy1,dz1)
            rij=dist_v2(dx1,dy1,dz1)
            Ebond_amber= Ebond_amber+kbond(ng1type(i,j))*
     .           (rij-bondeq(ng1type(i,j)))**2 
            dx1=(1.0/rij)*dx1
            dr(1)=2.0*kbond(ng1type(i,j))*(rij-bondeq(ng1type(i,j)))*dx1
            dy1=(1.0/rij)*dy1
            dr(2)=2.0*kbond(ng1type(i,j))*(rij-bondeq(ng1type(i,j)))*dy1  
            dz1=(1.0/rij)*dz1
            dr(3)=2.0*kbond(ng1type(i,j))*(rij-bondeq(ng1type(i,j)))*dz1
            do k=1,3
               fcebond_amber(k,i)=fcebond_amber(k,i)+dr(k)
               do l=1,3
                  stress_amber(l,k)=stress_amber(l,k)+
     .                       stress_fact*ramber(l,i)*dr(k)
               enddo
            enddo
         enddo
      enddo

      Ebond_amber= Ebond_amber/2
      ftotbond=0.0

      end
c******************************************************
c  subrutina q calcula la energia y fuerzas de angles
 
        subroutine  amber_angles(na_mm,ng1,ramber,
     .   Eangle_amber,attype,nangle,kangle,angleeq,angletype,
     .   angexat,angmxat,atange,atangm,fceangle_amber,stress_amber,
     .   angetype,angmtype,paso,nparm,cell,kcell,lattice_type)

       use precision, only: dp
        implicit none
c      vbles grales
       integer   na_mm,ng1(na_mm,6),i,j,k,l,m,n,paso
       real(dp)   ramber(3,na_mm),Eangle_amber,
     .                    fceangle_amber(3,na_mm)
       character   attype(na_mm)*4
c      vbles de los params de union
       integer   nangle,nparm
       real(dp) dx, dy, dz
       real(dp) kangle(nparm),angleeq(nparm)
       character angletype(nparm)*8
c      vbles de bond, angle, dihe e imp
       integer   angexat(na_mm),angmxat(na_mm)
       integer   atange(na_mm,25,2),atangm(na_mm,25,2)
c      vbles de asignacion
       character*4 ty1,ty2,ty3
       character*8 tyangle
       integer angetype(na_mm,25),angmtype(na_mm,25)
       real(dp)  angulo,angle,angle_v2,pi,dx12,dx32,dy12,dy32,
     .                   dz12,dz32,scal,r12,r32,
     .                   scalar,scalar_v2,ftotangle(3),dr12r32,dscalar,
     .                   dist,dist_v2,fesq(3,na_mm),fmedio(3,na_mm)
       real(dp) volcel
       real(dp)  cell(3,3), kcell(3,3)
       real(dp) stress_amber(3,3)
      character lattice_type
      logical           first                                                                  
      save              first                                                                  
      data              first /.true./ 
      real(dp), dimension(3) :: dr   
      real(dp) :: stress_fact
      
       pi=DACOS(-1.d0)
       fesq=0.0
       fmedio=0.0
       ftotangle=0.0

c asignacion de los tipos de angulos 
      if(first) then
       do i=1,na_mm
        do j=1,angexat(i)
         do k=1,nangle 
          tyangle=angletype(k)
          ty1=tyangle(1:2)
          ty2=tyangle(4:5)
          ty3=tyangle(7:8)
          if(attype(i).eq.ty1.and.attype(atange(i,j,1)).eq.ty2.and.
     .       attype(atange(i,j,2)).eq.ty3) then
          angetype(i,j)=k
          elseif(attype(i).eq.ty3.and.attype(atange(i,j,1)).eq.ty2.and.
     .       attype(atange(i,j,2)).eq.ty1) then
          angetype(i,j)=k
          endif
         enddo
        enddo
       enddo

       do i=1,na_mm
        do j=1,angmxat(i)
         do k=1,nangle
          tyangle=angletype(k)
          ty1=tyangle(1:2)
          ty2=tyangle(4:5)
          ty3=tyangle(7:8)
          if(attype(i).eq.ty2.and.attype(atangm(i,j,1)).eq.ty1.and.
     .       attype(atangm(i,j,2)).eq.ty3) then
          angmtype(i,j)=k
          elseif(attype(i).eq.ty2.and.attype(atangm(i,j,1)).eq.ty3.and.
     .       attype(atangm(i,j,2)).eq.ty1) then
          angmtype(i,j)=k
          endif
         enddo
        enddo
       enddo
      first=.false.
      endif ! asignacion

      stress_fact=1.0_dp/volcel(cell)

c calcula la E y F de angles
c para el angulo en la esquina
       do i=1,na_mm
        do j=1,angexat(i)
           dx12=ramber(1,i)-ramber(1,atange(i,j,1))
           dy12=ramber(2,i)-ramber(2,atange(i,j,1))
           dz12=ramber(3,i)-ramber(3,atange(i,j,1))
           dx32=ramber(1,atange(i,j,2))-ramber(1,atange(i,j,1))
           dy32=ramber(2,atange(i,j,2))-ramber(2,atange(i,j,1))
           dz32=ramber(3,atange(i,j,2))-ramber(3,atange(i,j,1))
           call pbc_displ_vector(lattice_type,cell,kcell,dx12,dy12,dz12)
           call pbc_displ_vector(lattice_type,cell,kcell,dx32,dy32,dz32)
           angulo=angle_v2(dx12,dy12,dz12,dx32,dy32,dz32)
           Eangle_amber = Eangle_amber + kangle(angetype(i,j))*
     .          ((angulo-angleeq(angetype(i,j)))*
     .          (pi/180))**2
           scal=scalar_v2(dx12,dy12,dz12,dx32,dy32,dz32)
           r12=dist_v2(dx12,dy12,dz12)
           r32=dist_v2(dx32,dy32,dz32)
           dr12r32=r32*dx12/(r12)
           dx=(dx32*r12*r32-scal*dr12r32)/(r12*r32)**(2.0)
           dx=-1.0/(sqrt(1.0-(scal/(r12*r32))**2.0))*dx
           dr(1)=2.0*kangle(angetype(i,j))*
     .          (angulo-angleeq(angetype(i,j)))*(pi/180)*dx
           dr12r32=r32*dy12/(r12)
           dy=(dy32*r12*r32-scal*dr12r32)/(r12*r32)**(2.0)
           dy=-1.0/(sqrt(1.0-(scal/(r12*r32))**2.0))*dy
           dr(2)=2.0*kangle(angetype(i,j))
     .          *(angulo-angleeq(angetype(i,j)))*(pi/180)*dy
           dr12r32=r32*dz12/(r12)
           dz=(dz32*r12*r32-scal*dr12r32)/(r12*r32)**(2.0)
           dz=-1.0/(sqrt(1.0-(scal/(r12*r32))**2.0))*dz
           dr(3)=2.0*kangle(angetype(i,j))
     .          *(angulo-angleeq(angetype(i,j)))*(pi/180)*dz
           do k=1,3
              fesq(k,i)=fesq(k,i)+dr(k)
              do l=1,3 
                 stress_amber(l,k)=stress_amber(l,k)+
     .                      stress_fact*ramber(l,i)*dr(k)
              enddo
           enddo
        enddo
       enddo  

c para el angulo en el medio   
       do i=1,na_mm
        do j=1,angmxat(i)
           dx12=ramber(1,atangm(i,j,1))-ramber(1,i)
           dy12=ramber(2,atangm(i,j,1))-ramber(2,i)
           dz12=ramber(3,atangm(i,j,1))-ramber(3,i)
           dx32=ramber(1,atangm(i,j,2))-ramber(1,i)
           dy32=ramber(2,atangm(i,j,2))-ramber(2,i)
           dz32=ramber(3,atangm(i,j,2))-ramber(3,i)
           call pbc_displ_vector(lattice_type,cell,kcell,dx12,dy12,dz12)
           call pbc_displ_vector(lattice_type,cell,kcell,dx32,dy32,dz32)
           angulo=angle_v2(dx12,dy12,dz12,dx32,dy32,dz32)
           Eangle_amber = Eangle_amber + kangle(angmtype(i,j))*
     .          ((angulo-angleeq(angmtype(i,j)))*
     .          (pi/180))**2
           scal=scalar_v2(dx12,dy12,dz12,dx32,dy32,dz32)
           r12=dist_v2(dx12,dy12,dz12)
           r32=dist_v2(dx32,dy32,dz32)
           dscalar=-dx12-dx32
           dr12r32=(r32*(-dx12)/r12)+(r12*(-dx32)/r32)
           dx=(dscalar*r12*r32-scal*dr12r32)/(r12*r32)**(2.0)
           dx=-1.0/(sqrt(1.0-(scal/(r12*r32))**2.0))*dx
           dr(1)=2.0*kangle(angmtype(i,j))
     .          *(angulo-angleeq(angmtype(i,j)))*(pi/180)*dx
           dscalar=-dy12-dy32
           dr12r32=(r32*(-dy12)/r12)+(r12*(-dy32)/r32)
           dy=(dscalar*r12*r32-scal*dr12r32)/(r12*r32)**(2.0)
           dy=-1.0/(sqrt(1.0-(scal/(r12*r32))**2.0))*dy
           dr(2)=2.0*kangle(angmtype(i,j))
     .          *(angulo-angleeq(angmtype(i,j)))*(pi/180)*dy
           dscalar=-dz12-dz32
           dr12r32=(r32*(-dz12)/r12)+(r12*(-dz32)/r32)
           dz=(dscalar*r12*r32-scal*dr12r32)/(r12*r32)**(2.0)
           dz=-1.0/(sqrt(1.0-(scal/(r12*r32))**2.0))*dz
           dr(3)=2.0*kangle(angmtype(i,j))
     .          *(angulo-angleeq(angmtype(i,j)))*(pi/180)*dz
           do k=1,3
              fmedio(k,i)=fmedio(k,i)+dr(k)
              do l=1,3 
                 stress_amber(l,k)=stress_amber(l,k)+
     .                    stress_fact*ramber(l,i)*dr(k)
              enddo
           enddo
        enddo
       enddo
    
       Eangle_amber= Eangle_amber/3

       do i=1,na_mm
          fceangle_amber(1,i)=fesq(1,i)+fmedio(1,i)
          fceangle_amber(2,i)=fesq(2,i)+fmedio(2,i)    
          fceangle_amber(3,i)=fesq(3,i)+fmedio(3,i)    
       enddo
        
	end
c****************************************************************
c  subrutina q calcula la energia y fuerzas de dihedros
 
      subroutine  amber_dihes(na_mm,ng1,ramber,Edihe_amber,
     .            attype,ndihe,kdihe,diheeq,perdihe,multidihe,
     .            dihetype,dihexat,dihmxat,atdihe,atdihm,
     .            fcedihe_amber,stress_amber,evaldihelog,evaldihe,
     .            dihety,
     .            evaldihmlog,evaldihm,dihmty,paso,nparm,cell,
     .            kcell,lattice_type)

       use precision, only: dp	
        implicit none
 
c      vbles grales
       integer   na_mm,ng1(na_mm,6),i,j,k,l,m,n,paso
       real(dp)   ramber(3,na_mm),Edihe_amber,
     .                    fcedihe_amber(3,na_mm)
       character   attype(na_mm)*4
c      vbles de los params de union
       integer   ndihe,nparm,multidihe(nparm)
       real(dp) kdihe(nparm),diheeq(nparm),perdihe(nparm)
       character dihetype(nparm)*11
c      vbles de bond, angle, dihe e imp
       integer   dihexat(na_mm),dihmxat(na_mm)
       integer   atdihe(na_mm,100,3),atdihm(na_mm,100,3)
c      vbles de asignacion
       character*4 ty1,ty2,ty3,ty4
       character*11 tydihe
       integer dihety(na_mm,100),dihmty(na_mm,100),evaldihe(na_mm,100,5)
     .         ,evaldihm(na_mm,100,5)
       real(dp)  pi,dihedro_v2,dihedral,E1,dist,
     .                   dx,dy,dz,ftotdihe(3),
     .                   fesq(3,na_mm),fmedio(3,na_mm),
     .			 fce(12)
       logical evaldihelog(na_mm,100),evaldihmlog(na_mm,100)
      logical           first                                                                  
      save              first                                                                  
      data              first /.true./  

      real(dp) volcel
      real(dp) cell(3,3), kcell(3,3)
      character lattice_type
      real(dp) stress_amber(3,3)
      real(dp) :: stress_fact

       pi=DACOS(-1.d0)
       ftotdihe=0.0
       fesq=0.0
       fmedio=0.0

c asignacion de los tipos de dihedros
      if(first) then
       evaldihelog=.false.
       evaldihmlog=.false.

       do i=1,na_mm
        do j=1,dihexat(i)
	dihety(i,j)=1
         m=0
         do k=1,ndihe
         tydihe=dihetype(k)
         ty1=tydihe(1:2)
         ty2=tydihe(4:5)   
         ty3=tydihe(7:8)    
         ty4=tydihe(10:11)   
         if(ty1.eq.'X ') then
	 if(attype(atdihe(i,j,1)).eq.ty2.and.
     .      attype(atdihe(i,j,2)).eq.ty3)  then
         dihety(i,j)=k
         elseif(attype(atdihe(i,j,1)).eq.ty3.and.
     .      attype(atdihe(i,j,2)).eq.ty2)  then
         dihety(i,j)=k
         endif
         elseif(ty1.ne.'X ') then
         if(attype(i).eq.ty1.and.attype(atdihe(i,j,1)).eq.
     .   ty2.and.attype(atdihe(i,j,2)).eq.ty3.and.
     .   attype(atdihe(i,j,3)).eq.ty4) then
         dihety(i,j)=k            
	  if(perdihe(k).lt.0) then
          evaldihelog(i,j)=.true.
          m=m+1
 	  evaldihe(i,j,m)=k
          evaldihe(i,j,m+1)=k+1
          endif
         elseif(attype(i).eq.ty4.and.attype(atdihe(i,j,1)).eq.
     .   ty3.and.attype(atdihe(i,j,2)).eq.ty2.and.
     .   attype(atdihe(i,j,3)).eq.ty1) then
         dihety(i,j)=k
          if(perdihe(k).lt.0) then
          evaldihelog(i,j)=.true.
          m=m+1
          evaldihe(i,j,m)=k
          evaldihe(i,j,m+1)=k+1    
          endif
         endif
         endif
         enddo

	if(dihety(i,j).eq.1) then
C	write(*,*) 'dihety: ',i,j,'sin parametro'
	endif

        enddo
       enddo

       do i=1,na_mm
        do j=1,dihmxat(i)
	dihmty(i,j)=1
         m=0
         do k=1,ndihe
         tydihe=dihetype(k)
         ty1=tydihe(1:2)
         ty2=tydihe(4:5)
         ty3=tydihe(7:8)
         ty4=tydihe(10:11)
         if(ty1.eq.'X ') then
         if(attype(i).eq.ty2.and.
     .      attype(atdihm(i,j,2)).eq.ty3)  then
         dihmty(i,j)=k
         elseif(attype(i).eq.ty3.and.
     .      attype(atdihm(i,j,2)).eq.ty2)  then
         dihmty(i,j)=k
         endif
         elseif(ty1.ne.'X ') then
         if(attype(atdihm(i,j,1)).eq.ty1.and.attype(i).eq.
     .   ty2.and.attype(atdihm(i,j,2)).eq.ty3.and.
     .   attype(atdihm(i,j,3)).eq.ty4) then
         dihmty(i,j)=k
          if(perdihe(k).lt.0.0) then
          evaldihmlog(i,j)=.true.
          m=m+1
          evaldihm(i,j,m)=k
          evaldihm(i,j,m+1)=k+1    
          endif
         elseif(attype(atdihm(i,j,1)).eq.ty4.and.attype(i).eq.
     .   ty3.and.attype(atdihm(i,j,2)).eq.ty2.and.
     .   attype(atdihm(i,j,3)).eq.ty1) then
         dihmty(i,j)=k
          if(perdihe(k).lt.0.0) then
          evaldihmlog(i,j)=.true.
          m=m+1
          evaldihm(i,j,m)=k
          evaldihm(i,j,m+1)=k+1    
          endif
         endif
         endif
         enddo
        if(dihmty(i,j).eq.1) then
C        write(*,*) 'dihmty: ',i,j,'sin parametro'
        endif
        enddo
       enddo
      first=.false.
      endif !asignacion

c calcula la E y F de dihedros 
c para los dihes en la esquina
        do i=1,ndihe
        if(perdihe(i).lt.0.0) then
        perdihe(i)=-perdihe(i)
        endif
        enddo

        stress_fact=1.0_dp/volcel(cell)

        do i=1,na_mm
         do j=1,dihexat(i)
         dihedral=dihedro_v2(ramber(1,i),ramber(2,i),ramber(3,i),
     .   ramber(1,atdihe(i,j,1)),ramber(2,atdihe(i,j,1)),
     .   ramber(3,atdihe(i,j,1)),
     .   ramber(1,atdihe(i,j,2)),ramber(2,atdihe(i,j,2)),
     .   ramber(3,atdihe(i,j,2)), 
     .   ramber(1,atdihe(i,j,3)),ramber(2,atdihe(i,j,3)),
     .   ramber(3,atdihe(i,j,3)),cell,kcell,lattice_type)


c si el dihedro es 500 es error y sale
        if(dihedral.lt.500.1.and.dihedral.gt.499.9) then
        write(*,*) 'ERROR EN EL DIHEDRO(esq)',i,j
        print*,'atom1 ',i,ramber(1,i),ramber(2,i),ramber(3,i)
        print*,'atom2 ',atdihe(i,j,1),ramber(1,atdihe(i,j,1)),
     .     ramber(2,atdihe(i,j,1)),ramber(3,atdihe(i,j,1))
        print*,'atom3 ',atdihe(i,j,2),ramber(1,atdihe(i,j,2)),
     .     ramber(2,atdihe(i,j,2)),ramber(3,atdihe(i,j,2))
        print*,'atom4 ',atdihe(i,j,3),ramber(1,atdihe(i,j,3)),
     .     ramber(2,atdihe(i,j,3)),ramber(3,atdihe(i,j,3))
        STOP
        endif

       if(evaldihelog(i,j)) then
	do m=1,5
	 if(evaldihe(i,j,m).ne.0) then

       k=evaldihe(i,j,m)
       E1=(kdihe(k)/multidihe(k))*
     .  (1+dCOS((pi/180)*(perdihe(k)*dihedral-diheeq(k))))
	Edihe_amber=Edihe_amber+E1
	call diheforce_v2(na_mm,ramber,
     .                 i,atdihe(i,j,1),atdihe(i,j,2),atdihe(i,j,3),1,
     .        kdihe(k),diheeq(k),perdihe(k),multidihe(k),fce,cell,kcell,
     .         lattice_type)    
        do n=1,3
           fesq(n,i)=fesq(n,i)+fce(n)
           do l=1,3
              stress_amber(l,n)=stress_amber(l,n)+
     .                            stress_fact*ramber(l,i)*fce(n)
           enddo
        enddo
      endif
      enddo
      else
         k=dihety(i,j)
         E1=(kdihe(k)/multidihe(k))*
     .        (1+dCOS((pi/180)*(perdihe(k)*dihedral-diheeq(k))))
         Edihe_amber=Edihe_amber+E1
         call diheforce_v2(na_mm,ramber,
     .        i,atdihe(i,j,1),atdihe(i,j,2),atdihe(i,j,3),1,
     .        kdihe(k),diheeq(k),perdihe(k),multidihe(k),fce,cell,kcell,
     .         lattice_type)
        do n=1,3
           fesq(n,i)=fesq(n,i)+fce(n)
           do l=1,3
              stress_amber(l,n)=stress_amber(l,n)+
     .                   stress_fact*ramber(l,i)*fce(n)
           enddo
        enddo
      endif
      enddo
      enddo

c para los dihes en el medio
       do i=1,na_mm
         do j=1,dihmxat(i)
         dihedral=dihedro_v2(
     .   ramber(1,atdihm(i,j,1)),ramber(2,atdihm(i,j,1)),
     .   ramber(3,atdihm(i,j,1)),
     .   ramber(1,i),ramber(2,i),ramber(3,i),
     .   ramber(1,atdihm(i,j,2)),ramber(2,atdihm(i,j,2)),
     .   ramber(3,atdihm(i,j,2)),
     .   ramber(1,atdihm(i,j,3)),ramber(2,atdihm(i,j,3)),
     .   ramber(3,atdihm(i,j,3)),cell,kcell,lattice_type)

c si el dihedro es 500 es error y sale
	if(dihedral.lt.500.1.and.dihedral.gt.499.9) then
	write(*,*) 'ERROR EN EL DIHEDRO(dihm)',i,j
	STOP
	endif

        if(evaldihmlog(i,j)) then
           do m=1,5
              if(evaldihm(i,j,m).ne.0) then
                 k=evaldihm(i,j,m)
                 E1=(kdihe(k)/multidihe(k))*
     .                (1+dCOS((pi/180)*(perdihe(k)*dihedral-diheeq(k))))
                 Edihe_amber=Edihe_amber+E1
                 call diheforce_v2(na_mm,ramber,
     .                atdihm(i,j,1),i,atdihm(i,j,2),atdihm(i,j,3),2,
     .                kdihe(k),diheeq(k),perdihe(k),multidihe(k),fce,
     .                    cell,kcell,lattice_type)
                 do n=1,3
                    fmedio(n,i)=fmedio(n,i)+fce(3+n)
                    do l=1,3
                       stress_amber(l,n)=stress_amber(l,n)+
     .                           stress_fact*ramber(l,i)*fce(3+n)
                    enddo
                 enddo
              endif
           enddo
        else
           k=dihmty(i,j)
           E1=(kdihe(k)/multidihe(k))*
     .          (1+dCOS((pi/180)*(perdihe(k)*dihedral-diheeq(k))))
           Edihe_amber=Edihe_amber+E1
           call diheforce_v2(na_mm,ramber,
     .          atdihm(i,j,1),i,atdihm(i,j,2),atdihm(i,j,3),2,
     .          kdihe(k),diheeq(k),perdihe(k),multidihe(k),fce,cell,
     .          kcell,lattice_type)
           do n=1,3
              fmedio(n,i)=fmedio(n,i)+fce(3+n)
              do l=1,3
                 stress_amber(l,n)=stress_amber(l,n)+
     .                stress_fact*ramber(l,i)*fce(3+n)
              enddo
           enddo
        endif
      enddo
      enddo

      Edihe_amber=Edihe_amber/4

      do i=1,na_mm
         fcedihe_amber(1,i)=fesq(1,i)+fmedio(1,i)
         fcedihe_amber(2,i)=fesq(2,i)+fmedio(2,i)
         fcedihe_amber(3,i)=fesq(3,i)+fmedio(3,i)
      enddo 

       end
c******************************************************************
c  subrutina q calcula la energia y fuerzas de impropers 

       subroutine amber_improper(na_mm,ng1,ramber,Eimp_amber,attype,
     .            nimp,kimp,impeq,perimp,multiimp,imptype,impxat,
     .            atimp,fimp,stress_amber,impty,paso,nparm,cell,kcell,
     .             lattice_type)

       use precision, only: dp
        implicit none
c      vbles grales
       integer   na_mm,ng1(na_mm,6),i,j,k,l,m,n,paso
       real(dp)   ramber(3,na_mm),Eimp_amber
       character   attype(na_mm)*4
c      vbles de los params de union
       integer   nimp,nparm,multiimp(nparm)
       real(dp) kimp(nparm),impeq(nparm),perimp(nparm)
       character imptype(nparm)*11 
c      vbles de bond, angle, dihe e imp
       integer   impxat(na_mm)
       integer   atimp(na_mm,25,4)
c      vbles de asignacion
       character*4 ty1,ty2,ty3,ty4
       character*11 tyimp
       integer impty(na_mm,25)
       real(dp)  pi,dihedro_v2,dihedral
c	varianles de fza
	real(dp) fimp(3,na_mm),fce(12),fimptot(3),dx,dy,dz
	integer atom
      logical           first                                                                  
      save              first                                                                  
      data              first /.true./
      real(dp) volcel
      real(dp) cell(3,3), kcell(3,3)
      character lattice_type
      real(dp) stress_amber(3,3)
      real(dp) stress_fact
    
       pi=DACOS(-1.d0)
       fimptot=0.0

c asignacion de los tipos de impropers
      if(first) then
      do i=1,na_mm
       do j=1,impxat(i)
	impty(i,j)=1
        do k=1,nimp
        tyimp=imptype(k)
        ty1=tyimp(1:2)
        ty2=tyimp(4:5)
        ty3=tyimp(7:8)
        ty4=tyimp(10:11)
	if(ty1.eq.'X '.and.ty2.eq.'X '.and.ty4.eq.'X ') then
	       if(attype(atimp(i,j,3)).eq.ty3) then
        impty(i,j)=k
		elseif(attype(atimp(i,j,2)).eq.ty3) then
        impty(i,j)=k
		endif
        elseif(ty1.eq.'X '.and.ty2.eq.'X ') then
        if(attype(atimp(i,j,3)).eq.ty3.and.
     .     attype(atimp(i,j,4)).eq.ty4) then
        impty(i,j)=k
        elseif(attype(atimp(i,j,1)).eq.ty4.and.
     .     attype(atimp(i,j,2)).eq.ty3) then
        impty(i,j)=k   
        endif
        elseif(ty1.eq.'X ') then
        if(attype(atimp(i,j,2)).eq.ty2.and.
     .     attype(atimp(i,j,3)).eq.ty3.and.
     .     attype(atimp(i,j,4)).eq.ty4) then
        impty(i,j)=k   
        elseif(attype(atimp(i,j,3)).eq.ty2.and.
     .     attype(atimp(i,j,2)).eq.ty3.and.
     .     attype(atimp(i,j,1)).eq.ty4) then
        impty(i,j)=k
        endif
        else
        if(attype(atimp(i,j,1)).eq.ty1.and.
     .     attype(atimp(i,j,2)).eq.ty2.and.
     .     attype(atimp(i,j,3)).eq.ty3.and.
     .     attype(atimp(i,j,4)).eq.ty4) then
        impty(i,j)=k
        elseif(attype(atimp(i,j,4)).eq.ty1.and.
     .     attype(atimp(i,j,3)).eq.ty2.and.
     .     attype(atimp(i,j,2)).eq.ty3.and.
     .     attype(atimp(i,j,1)).eq.ty4) then
        impty(i,j)=k
        endif
        endif
        enddo
       if(impty(i,j).eq.1) then
C        write(*,*) 'impty: ',i,j,'sin parametro'
        endif
       enddo
      enddo
      first=.false.
      endif !asignacion

      stress_fact=1.0_dp/volcel(cell)

c calcula la E y F de impropers
      do i=1,nimp
        if(perimp(i).lt.0.0) then
        perimp(i)=-perimp(i)
        endif
        enddo
        do i=1,na_mm
         do j=1,impxat(i)
         dihedral=dihedro_v2(ramber(1,atimp(i,j,1)),
     .   ramber(2,atimp(i,j,1)),ramber(3,atimp(i,j,1)),
     .   ramber(1,atimp(i,j,2)),ramber(2,atimp(i,j,2)),
     .   ramber(3,atimp(i,j,2)),
     .   ramber(1,atimp(i,j,3)),ramber(2,atimp(i,j,3)),
     .   ramber(3,atimp(i,j,3)),
     .   ramber(1,atimp(i,j,4)),ramber(2,atimp(i,j,4)),
     .   ramber(3,atimp(i,j,4)),cell,kcell,lattice_type)

c si el dihedro es 500 es error y sale
        if(dihedral.lt.500.1.and.dihedral.gt.499.9) then
        write(*,*) 'ERROR EN EL IMPROP',i,j
        STOP
        endif

       k=impty(i,j)
       Eimp_amber=Eimp_amber+(kimp(k)/multiimp(k))*
     .  (1+COS((pi/180)*(perimp(k)*dihedral-impeq(k))))

c     busca que atomo del impropio es el i	
       if (i.eq.atimp(i,j,1)) atom=1
       if (i.eq.atimp(i,j,2)) atom=2 
       if (i.eq.atimp(i,j,3)) atom=3 
       if (i.eq.atimp(i,j,4)) atom=4 
       call diheforce(na_mm,ramber,
     .      atimp(i,j,1),atimp(i,j,2),atimp(i,j,3),atimp(i,j,4),atom,
     .      kimp(k),impeq(k),perimp(k),multiimp(k),fce,cell,kcell,
     .      lattice_type)
       do n=1,3
          fimp(n,i)=fimp(n,i)+fce((atom-1)*3+n)
          do l=1,3
             stress_amber(l,n)=stress_amber(l,n)+stress_fact*ramber(l,i)
     .            *fce((atom-1)*3+n)
          enddo
       enddo
      enddo
      enddo

      Eimp_amber=Eimp_amber/4

      end
c *******************************************************
c subrutina que calcula la energia y fuerzas de terminos non-bonded

       subroutine amber_nonbonded(na_mm,ng1,ramber,Elj_amber,
     .       Eelec_amber,Elj_amber14,Eelec_amber14,attype,
     .       Em,Rm,pc,bondxat,
     .       angexat,angmxat,atange,atangm,
     .       dihexat,dihmxat,atdihe,atdihm,
     .       felec,flj,stress_amber,nonbonded,coulombtype,
     .       scale,nonbondedxat,scalexat,paso,
     .       actualiz,listcut,noat,noaa,
     .       sfc,timestep,
     .       na_qm,natot,rclas,
     .       graphite_layer_no,water,masst,ewat,fwat,cell,
     .       ewald_alpha,kewald_cutoff,kcell,lattice_type)       

       use precision, only: dp
      use qmmm_neighbour

        implicit none
c      vbles grales
       integer   na_mm,ng1(na_mm,6),i,j,k,l,m,n,paso,dimvec,
     . n_pointer,na_qm,natot
       real(dp)   ramber(3,na_mm),Em(na_mm),Rm(na_mm),pc(na_mm),
     .    Elj_amber,Eelec_amber,Elj_amber14,Eelec_amber14
       character*4 attype(na_mm),noat(na_mm),noaa(na_mm), coulombtype*10
       real(dp) ewat,fwat(3,na_mm),masst(natot),rclas(3,natot)
       logical water
c      vbles de bond, angle, dihe e imp
       integer   bondxat(na_mm),angexat(na_mm),angmxat(na_mm),
     .   dihexat(na_mm),dihmxat(na_mm)
       integer   atange(na_mm,25,2),atangm(na_mm,25,2),
     .   atdihe(na_mm,100,3),atdihm(na_mm,100,3)
c      vbles de asignacion      
       integer nonbonded(na_mm,100),scale(na_mm,100),scalexat(na_mm),
     . nonbondedxat(na_mm),acs
       real(dp) A,B,Rij,Eij,dist,distancia,unidades,
     .                  factorlj,factorelec,E1,E2,pi,epsilon,
     .                  dist2,distancia2,rcut,Rij6,distancia2_3,
     .                  fac,dfac,ca,cb,cc,cd,sfc,timestep,e1f,e2f,
     .                  x0,x1,rinn,rout
       logical scale14,nonbond,actualiz 
c variables asoc a las fzas
       real(dp) dr(3)
       real(dp) dx1,dy1,dz1,dx2,dy2,dz2,felec(3,na_mm),
     .  felectot(3),fel,flj(3,na_mm),fljtot(3),listcut
C variable for interaction between graphite layers
       integer graphite_layer_no(na_mm)
c variables de la lista de vecinos
       integer, dimension(:), allocatable ::  veclist, veclistemp,       
     .   veclistxat
       save veclist, veclistxat
       integer, dimension(:,:), allocatable ::  nr, nrtemp
       save nr
       integer, dimension(3) :: nr_indx
       real(dp), dimension(3) :: drij
       character lattice_type
      real(dp) dx, dy, dz, const1, const2, const3, const4, const5
      real(dp) const6, const7, kmod2, kronij
      integer n1, n2, n3, n1max, n2max, n3max
      real(dp) kewald_cutoff, kcut2
      logical           first                                        
      save              first                                    
      data              first /.true./   
 
      real(dp) cell(3,3), kcell(3,3)
      real(dp) stress_amber(3,3)
      real(dp) De, De2
      real(dp) ewald_alpha,sqrt_ewald_alpha
      integer ewald_nmax
      parameter (ewald_nmax=50)
      real(dp) S_real, S_imag
      real(dp) twopi, kr, lattice_volume
      real(dp) pc_cos_kr(na_mm), pc_sin_kr(na_mm), SS
      real(dp) volcel, scalar_v2
      real(dp) Erecip_amber
      real(dp) krecip(3)
      real(dp) :: stress_fact

      integer   nna
      save nna
      data nna /200/  
      integer in
      real(dp) lattice_vector_len, vector_len_max
      real(dp) rcoor
      save rcoor
      data rcoor /0.0d0/
      logical bigcell
      save bigcell
      data bigcell /.true./

      pi=DACOS(-1.d0)
      twopi=2.0d0*pi
      unidades=((1.602177E-19**2)*6.02E23)/(1.0E-10*4*
     .     pi*8.8541878E-12*4.184*1000)
      fel=0.0
      flj=0.0
      felec=0.0
      felectot=0.0
      fljtot=0.0
      epsilon=1.0
      factorlj=2.0
      factorelec=1.2
      acs=200	
      ewat=0.0
      fwat=0.0
      dimvec=na_mm*3000
      if(dimvec.gt.100000000) stop 
     .     'Solvent Energy and Forces: "dimvec" too large!'

      sqrt_ewald_alpha=sqrt(ewald_alpha)
      const1=sqrt_ewald_alpha/sqrt(acos(-1.0d0))

      lattice_volume=volcel(cell)

      stress_fact=1.0_dp/lattice_volume

      const2=4.0*acos(-1.0d0)/lattice_volume

      call reccel(3,cell,kcell,0)

      n1max=INT(kewald_cutoff/(twopi*
     .                  sqrt(scalar_v2(kcell(1,1),kcell(1,2),kcell(1,3),
     .                        kcell(1,1),kcell(1,2),kcell(1,3)))))
      n2max=INT(kewald_cutoff/(twopi*
     .                  sqrt(scalar_v2(kcell(2,1),kcell(2,2),kcell(2,3),
     .                        kcell(2,1),kcell(2,2),kcell(2,3)))))
      n3max=INT(kewald_cutoff/(twopi*
     .                  sqrt(scalar_v2(kcell(3,1),kcell(3,2),kcell(3,3),
     .                        kcell(3,1),kcell(3,2),kcell(3,3)))))
      kcut2=kewald_cutoff*kewald_cutoff

c     asigna los atomos nonbonded y su tipo
      if(first) then       

         do i=1,na_mm 
            m=1
            do j=1,bondxat(i)
               if(i.lt.ng1(i,j))then
                  nonbonded(i,m)=ng1(i,j)
                  m=m+1
               endif
            enddo
            do j=1,angexat(i)
               if(i.lt.atange(i,j,2)) then
                  nonbonded(i,m)=atange(i,j,2)
                  m=m+1
               endif
            enddo
            nonbondedxat(i)=m-1

c     se fija los que estan 1-4(scaled)
            m=1
            do j=1,dihexat(i)
               if(i.lt.atdihe(i,j,3)) then
c     se fija si ya lo puso (por si hay dihedros repetidos)
                  do n=1,m-1	  
                     if(atdihe(i,j,3).eq.scale(i,n)) goto 10  
                  enddo
c     se fija si no lo puso en nonbonded (se no es tambien 1-3)
                  do n=1,nonbondedxat(i)
                     if(atdihe(i,j,3).eq.nonbonded(i,n)) goto 10	
                  enddo	

                  scale(i,m)=atdihe(i,j,3)
                  m=m+1
               endif
 10         enddo
            scalexat(i)=m-1
         enddo

C         actualiz=.true.
         first=.false.

         vector_len_max=0.0d0
         do i=1,3
            lattice_vector_len=sqrt(cell(i,1)*cell(i,1)+
     .            cell(i,2)*cell(i,2)+cell(i,3)*cell(i,3))
            if (lattice_vector_len>vector_len_max) vector_len_max=
     .                                       lattice_vector_len
         enddo

!     sfc: smooth-function cut-off
!     skin in the neighbor list of 2 angstroms
         rcoor=listcut+sfc+2.

         if (2.0*rcoor<vector_len_max) then
             bigcell=.true.
         else
             bigcell=.false.
         endif

C     Initialize routine for neighbour search. This call creates the cells
         call mm_mneighb(cell,rcoor,na_mm,ramber,0,1,nna) 
         if(associated(mm_jan)) deallocate(mm_jan)
         allocate(mm_jan(maxnna))
         if(associated(mm_r2ij)) deallocate(mm_r2ij)
         allocate(mm_r2ij(maxnna))
         if(associated(mm_xij)) deallocate(mm_xij)
         allocate(mm_xij(3,maxnna))

      endif                     !asignacion

c     ahora crea la lista de vecinos si estan detro de listcut+sfc, 
c     se debe actualizar cada 10 fs. 
C      acs=int(10/timestep)

C      if(mod(paso,acs).eq.0.or.paso.eq.1) actualiz=.true.

c     llama a la sub q agrega el water restrain potential
      if(water.and.actualiz) then
         call waters(na_qm,na_mm,natot,rclas,masst,
     $        noaa,noat,ewat,fwat,cell,kcell)
      endif

c     si listcut > 100A, la lista se hace SOLO el 1er paso 
C      if(listcut.ge.100.and.paso.ne.1) actualiz=.false.

      if (actualiz) then

C     Initialize routine for neighbour search. This call creates the cells
         call mm_mneighb(cell,rcoor,na_mm,ramber,0,1,nna) 
         if(associated(mm_jan)) deallocate(mm_jan)
         allocate(mm_jan(maxnna))
         if(associated(mm_r2ij)) deallocate(mm_r2ij)
         allocate(mm_r2ij(maxnna))
         if(associated(mm_xij)) deallocate(mm_xij)
         allocate(mm_xij(3,maxnna))

         if(allocated(veclist)) deallocate(veclist)
         allocate(veclistemp(dimvec))
         if(allocated(nr)) deallocate(nr)
         allocate(nrtemp(3,dimvec))
         if(allocated(veclistxat)) deallocate(veclistxat)
         allocate(veclistxat(na_mm))

c     crea la lista veclist y el indice veclistxat
c     write(6,*) 'Neighbour list actualization in step:',paso

         veclistemp=0
         nrtemp=0
         veclistxat=0

         rcoor=listcut+sfc+2.

         rcut=rcoor**2
         n=1

         if (bigcell) then

            do i=1,na_mm
C     Look for neighbours of atom ia
               call mm_mneighb(cell,rcoor,na_mm,ramber,i,0,nna)
               IF (NNA .GT. MAXNNA) STOP 'MAXNNA too small'
               do in=1,nna
                  if (mm_r2ij(in).eq.0.0) cycle
                  j=mm_jan(in)
                  if (i.gt.j) cycle
                  if(noaa(j).eq.'HOH') then
                     if (noat(j).ne.'O') cycle
                  endif
!     exclude atoms connected by 1 or 2 covalent bonds
                  do k=1,nonbondedxat(i)
                     if(nonbonded(i,k).eq.j) goto 30
                  enddo
!     exclude atoms 1 and 4 in covalently-bonded chain 1-2-3-4
                  do m=1,scalexat(i)
                     if(scale(i,m).eq.j) goto 30
                  enddo
                  drij(1:3)=mm_xij(1:3,in)-ramber(1:3,j)+
     .                 ramber(1:3,i)
                  call get_pbc_vectors(lattice_type,cell,kcell,
     .                 drij,nr_indx)
                  if(noaa(j).eq.'HOH') then
                     veclistemp(n)=j
                     nrtemp(1:3,n)=nr_indx(1:3)
                     n=n+1
                     veclistemp(n)=j+1
                     nrtemp(1:3,n)=nr_indx(1:3)
                     n=n+1
                     veclistemp(n)=j+2
                     nrtemp(1:3,n)=nr_indx(1:3)
                     n=n+1
                  else if(noaa(j).ne.'HOH') then
                     veclistemp(n)=j
                     nrtemp(1:3,n)=nr_indx(1:3)
                     n=n+1
                  endif
c     fin loop de la lista para cada atomo
 30            enddo
               veclistxat(i)=n-1
c     se fija si la dimension de veclist es suficiente
               if((n-1).gt.dimvec) then
                  write(6,*)'Dimension Neighbour list (required, used)='
     .                 ,n-1,dimvec
                  write(6,*) 'Solvent Energy and Forces: '//
     .                 'Stopping Program'
                  STOP
               endif
            enddo

         else

            do i=1,na_mm

C     Look for neighbours of atom ia
               call mm_mneighb(cell,rcoor,na_mm,ramber,i,0,nna)
               IF (NNA .GT. MAXNNA) STOP 'MAXNNA too small'

               do in=1,nna
                  if (mm_r2ij(in).eq.0.0) cycle
                  j=mm_jan(in)
                  if (i.gt.j) cycle
                  if(noaa(j).eq.'HOH') then
                     if (noat(j).ne.'O') cycle
                  endif
                  drij(1)=ramber(1,i)-ramber(1,j)
                  drij(2)=ramber(2,i)-ramber(2,j)
                  drij(3)=ramber(3,i)-ramber(3,j)
                  call pbc_displ_vector(lattice_type,cell,kcell,drij(1)
     .                 ,drij(2),drij(3))
                  distancia2 = drij(1)*drij(1) + drij(2)*drij(2)  + 
     .                            drij(3)*drij(3)
                  if (abs(mm_r2ij(in)-distancia2)<0.1d0) then
!     exclude atoms connected by 1 or 2 covalent bonds
                     do k=1,nonbondedxat(i)
                        if(nonbonded(i,k).eq.j) goto 40
                     enddo
!     exclude atoms 1 and 4 in covalently-bonded chain 1-2-3-4
                     do m=1,scalexat(i)
                        if(scale(i,m).eq.j) goto 40
                     enddo
                  endif
                  drij(1:3)=mm_xij(1:3,in)-ramber(1:3,j)+
     .                 ramber(1:3,i)
                  call get_pbc_vectors(lattice_type,cell,kcell,
     .                 drij,nr_indx)
                  if(noaa(j).eq.'HOH') then
                     veclistemp(n)=j
                     nrtemp(1:3,n)=nr_indx(1:3)
                     n=n+1
                     veclistemp(n)=j+1
                     nrtemp(1:3,n)=nr_indx(1:3)
                     n=n+1
                     veclistemp(n)=j+2
                     nrtemp(1:3,n)=nr_indx(1:3)
                     n=n+1
                  else if(noaa(j).ne.'HOH') then
                     veclistemp(n)=j
                     nrtemp(1:3,n)=nr_indx(1:3)
                     n=n+1
                  endif
c     fin loop de la lista para cada atomo
 40            enddo
               veclistxat(i)=n-1
c     se fija si la dimension de veclist es suficiente
               if((n-1).gt.dimvec) then
                  write(6,*)'Dimension Neighbour list (required, used)='
     .                 ,n-1,dimvec
                  write(6,*) 'Solvent Energy and Forces: '//
     .                 'Stopping Program'
                  STOP
               endif
            enddo
c     fin loop de lista para todos atomos
         endif

c     alocate a veclist posta
         allocate(veclist(n-1))
         veclist(1:n-1)=veclistemp(1:n-1)
         deallocate(veclistemp)
         allocate(nr(3,n-1))
         nr(1:3,1:n-1)=nrtemp(1:3,1:n-1)
         deallocate(nrtemp)

         actualiz=.false.

      endif

c     calculo de la E y F de terminos non-bonded
c     Starts loop 1,4-scaled nonbonded interactions
      do i=1,na_mm
         do k=1,scalexat(i)
            j=scale(i,k)
            if (.not.(noaa(i).eq.'GRAP'.and.noaa(j).eq.'GRAP')) then	
               dr(1)=ramber(1,i)-ramber(1,j)
               dr(2)=ramber(2,i)-ramber(2,j)
               dr(3)=ramber(3,i)-ramber(3,j)
               call pbc_displ_vector(lattice_type,cell,kcell,dr(1),
     .                          dr(2),dr(3))
               distancia2 = dr(1)*dr(1) + dr(2)*dr(2) + dr(3)*dr(3)
               distancia = sqrt(distancia2)
               if (distancia.le.listcut) then
                  if (graphite_layer_no(i)/=graphite_layer_no(j)) then
                     Rij=Rm(i)+Rm(j)
                     Eij=sqrt(Em(i)*Em(j))
                     Rij6 = Rij**6
                     distancia2_3 = distancia2**3
                     E1 = Eij*Rij6/distancia2_3*((Rij6/distancia2_3)-2.)
                     E1=E1/factorlj
                     Elj_amber14=Elj_amber14+E1
                     fel = -12.0*Eij*Rij6/distancia2**4* 
     .                   (Rij6/distancia2_3 - 1.0)
                     fel = fel/factorlj
                     do m=1,3
                        flj(m,i)=flj(m,i)+dr(m)*fel
                        flj(m,j)=flj(m,j)-dr(m)*fel
                        do n=1,3
                           stress_amber(n,m)=stress_amber(n,m)+
     .                         stress_fact*dr(n)*dr(m)*fel
                        enddo
                     enddo
                  endif
                  if (coulombtype.eq.'ewald') then
                     call coulomb_ewald_real(pc(i),pc(j),distancia,
     .                    E2,fel,ewald_alpha,sqrt_ewald_alpha,const1)
                  else
                     call coulomb_cutoff(pc(i),pc(j),distancia,
     .                    distancia2,E2,fel)
                  endif
                  E2=E2/factorelec*unidades/epsilon
                  Eelec_amber14=Eelec_amber14+E2
                  fel=fel/factorelec*unidades/epsilon
                  do m=1,3
                     felec(m,i)=felec(m,i)+dr(m)*fel
                     felec(m,j)=felec(m,j)-dr(m)*fel
                     do n=1,3
                        stress_amber(n,m)=stress_amber(n,m)+
     .                            stress_fact*dr(n)*dr(m)*fel
                     enddo
                  enddo
               endif
            endif
         enddo
      enddo
c     end of scaled nonbonden

c     Starts loop nonscaled nonbonded
      n_pointer=1      
      x0=listcut
      x1=listcut+sfc
      if (coulombtype.eq.'ewald') then
         cb=-sfc/x0*(erfc(sqrt_ewald_alpha
     .      *x0)/x0+2.0*const1*exp(-ewald_alpha*x0**2))
         ca=-cb/2.
         cc=ca
         cd=cc-1./x0*erfc(sqrt_ewald_alpha*x0)
      else
         cb=-sfc/x0**2
         ca=-cb/2.
         cc=ca
         cd=cc-1./x0
      endif
      rinn=x0**2
      rout=x1**2

      do i=1,na_mm
C     n_pointer: points the first neb atom of i in the neb list.
         do k=n_pointer,veclistxat(i)
            j=veclist(k)
            if (graphite_layer_no(i)/=graphite_layer_no(j).or.
     .               graphite_layer_no(i)*graphite_layer_no(j)==0) then
               if (lattice_type=='D') then
                  drij(1)=ramber(1,i)-ramber(1,j)+nr(1,k)*cell(1,1)
                  drij(2)=ramber(2,i)-ramber(2,j)+nr(2,k)*cell(2,2)
                  drij(3)=ramber(3,i)-ramber(3,j)+nr(3,k)*cell(3,3)
               else
                  drij(1)=ramber(1,i)-ramber(1,j)
                  drij(2)=ramber(2,i)-ramber(2,j)
                  drij(3)=ramber(3,i)-ramber(3,j)
                  do l=1,3
                     do m=1,3
                        drij(l)=drij(l)+nr(m,k)*cell(l,m)
                     enddo
                  enddo
               endif
               distancia2=drij(1)*drij(1)+drij(2)*drij(2)
     .                             +drij(3)*drij(3)
C     Non-bonded interaction is smoothly decayed from sqrt(rinn) 
C     to sqrt(rout)
               if(distancia2.le.rinn) then
                  distancia = sqrt(distancia2)
                  Rij=Rm(i)+Rm(j)
                  Eij=sqrt(Em(i)*Em(j))
                  Rij6 = Rij**6
                  distancia2_3 = distancia2**3
                  E1 = Eij*Rij6/distancia2_3*((Rij6/distancia2_3)-2.)
                  Elj_amber=Elj_amber+E1
                  fel = -12.0*Eij*Rij6/distancia2**4*(Rij6/distancia2_3 
     .                 - 1.0)
                  do m=1,3
                     flj(m,i)=flj(m,i)+drij(m)*fel
                     flj(m,j)=flj(m,j)-drij(m)*fel
                     do n=1,3
                        stress_amber(n,m)=stress_amber(n,m)+
     .                             stress_fact*drij(n)*drij(m)*fel
                     enddo
                  enddo
                  if (coulombtype.eq.'ewald') then
                     call coulomb_ewald_real(pc(i),pc(j),distancia,
     .                    E2,fel,ewald_alpha,sqrt_ewald_alpha,const1)
                  else
                     call coulomb_cutoff(pc(i),pc(j),distancia,
     .                    distancia2,E2,fel)
                  endif
                  E2=E2*unidades/epsilon
C     Add a constant to E2 to make the energy go to 
C     E=pc(i)*pc(j)/listcut*cc=
C     pc(i)*pc(j)/listcut*(sfc/(2.0*listcut**2))
C     when distancia=listcut and where
C     sfc: smooth-function cut-off
                  E2F=E2+(pc(i)*pc(j)*unidades/epsilon*cd)
                  Eelec_amber=Eelec_amber+E2F
                  fel=fel*unidades/epsilon
                  do m=1,3
                     felec(m,i)=felec(m,i)+drij(m)*fel
                     felec(m,j)=felec(m,j)-drij(m)*fel
                     do l=1,3
                        stress_amber(l,m)=stress_amber(l,m)+
     .                            stress_fact*drij(l)*drij(m)*fel
                     enddo
                  enddo
               elseif(distancia2.gt.rinn.and.distancia2.lt.rout) then
                  distancia = sqrt(distancia2) 
                  Rij=Rm(i)+Rm(j)
                  Eij=sqrt(Em(i)*Em(j))
                  Rij6 = Rij**6
                  distancia2_3 = distancia2**3
                  E1 = Eij*Rij6/distancia2_3*((Rij6/distancia2_3)-2.)
                  Elj_amber=Elj_amber+E1
                  fel = -12.0*Eij*Rij6/distancia2**4*(Rij6/distancia2_3 
     .                 - 1.0)
                  do m=1,3
                     flj(m,i)=flj(m,i)+drij(m)*fel
                     flj(m,j)=flj(m,j)-drij(m)*fel
                     do n=1,3
                        stress_amber(n,m)=stress_amber(n,m)+
     .                                stress_fact*drij(n)*drij(m)*fel
                     enddo
                  enddo
                  E2=pc(i)*pc(j)*unidades/epsilon
C     Add cutt-off switching function
                  fac=(distancia-listcut)/sfc
                  E2F=E2*(ca*fac**2+cb*fac+cc)
                  Eelec_amber=Eelec_amber+E2F
                  fel=E2/distancia
                  fel=fel*(2*ca*fac+cb)/sfc
                  do m=1,3
                     felec(m,i)=felec(m,i)+drij(m)*fel
                     felec(m,j)=felec(m,j)-drij(m)*fel
                     do n=1,3
                        stress_amber(n,m)=stress_amber(n,m)+
     .                             stress_fact*drij(n)*drij(m)*fel
                     enddo
                  enddo
               endif
            endif
c     fin del loop sobre bondlist
         enddo
         n_pointer = veclistxat(i) + 1
c     fin del loop sobre todos los atomos
      enddo

      if (coulombtype.eq.'ewald') then

         Erecip_amber=0.0_dp
C     RECIPROCAL-SPACE PART OF EWALD SUM

C     Calculate structure factors for all MM atoms:
         const4=0.5_dp/lattice_volume
         const6=1.0_dp/(4.0_dp*ewald_alpha)
         if (lattice_type.eq.'D') then
!     Reciprocal-space sum
            do n1=-n1max,n1max
               do  n2=-n2max,n2max
                  do  n3=-n3max,n3max
                     if (.not.(n1.eq.0.and.n2.eq.0.and.n3.eq.0)) then
!     loop over MM atoms
                        krecip(1)=n1*twopi*kcell(1,1)
                        krecip(2)=n2*twopi*kcell(2,2)
                        krecip(3)=n3*twopi*kcell(3,3)
                        kmod2=krecip(1)*krecip(1)+krecip(2)*krecip(2)
     .                       +krecip(3)*krecip(3)
                        if (kmod2<=kcut2) then
                           S_real=0.0d0
                           S_imag=0.0d0
                           do i=1,na_mm
                              kr=krecip(1)*ramber(1,i)
     .                                  +krecip(2)*ramber(2,i)
     .                                       +krecip(3)*ramber(3,i)
                              pc_cos_kr(i)=pc(i)*cos(kr)
                              S_real=S_real+pc_cos_kr(i)
                              pc_sin_kr(i)=pc(i)*sin(kr)
                              S_imag=S_imag+pc_sin_kr(i)
                           enddo
                           const3=const2/kmod2*
     .                          exp(-kmod2/(4.0d0*ewald_alpha))*
     .                             unidades/epsilon
                           SS=0.5d0*const3*(S_real*S_real+S_imag*S_imag)
                           Eelec_amber=Eelec_amber + SS
                           Erecip_amber=Erecip_amber + SS
                           do i=1,na_mm
                              De=const3*(pc_cos_kr(i)*S_imag
     $                             -pc_sin_kr(i)*S_real)
                              do k=1,3
                                 felec(k,i) = felec(k,i) + De*krecip(k)
                              enddo
                           enddo
                           const5=const3*const4
                           De2=const5*(S_imag*S_imag+S_real*S_real)
                           const7=2.0_dp*(1.0_dp+kmod2*const6)/kmod2
                           do k=1,3
                              do l=1,3
                                 kronij=real(int(((l+k)-abs(l-k))/
     $                                ((l+k)+abs(l-k))),kind=dp)
                                 stress_amber(l,k)=stress_amber(l,k)+
     $                                De2*(kronij-const7*krecip(l)*
     $                                krecip(k))
                              enddo
                           enddo
                        endif
                     endif
                  enddo
               enddo
            enddo
         else
!     Reciprocal-space sum
            do n1=-n1max,n1max
               do  n2=-n2max,n2max
                  do  n3=-n3max,n3max
                     if (.not.(n1.eq.0.and.n2.eq.0.and.n3.eq.0)) then
!     loop over MM atoms
                        krecip(1)=twopi*(n1*kcell(1,1)+
     .                       n2*kcell(2,1)+n3*kcell(3,1))
                        krecip(2)=twopi*(n1*kcell(1,2)+
     .                       n2*kcell(2,2)+n3*kcell(3,2))
                        krecip(3)=twopi*(n1*kcell(1,3)+
     .                       n2*kcell(2,3)+n3*kcell(3,3))
                        kmod2=krecip(1)*krecip(1)+krecip(2)*krecip(2)
     .                       +krecip(3)*krecip(3)
                        if (kmod2<=kcut2) then
                           S_real=0.0d0
                           S_imag=0.0d0
                           do i=1,na_mm
                              kr=krecip(1)*ramber(1,i)
     .                                  +krecip(2)*ramber(2,i)
     .                                       +krecip(3)*ramber(3,i)
                              pc_cos_kr(i)=pc(i)*cos(kr)
                              S_real=S_real+pc_cos_kr(i)
                              pc_sin_kr(i)=pc(i)*sin(kr)
                              S_imag=S_imag+pc_sin_kr(i)
                           enddo
                           const3=const2/kmod2*
     .                          exp(-kmod2/(4.0*ewald_alpha))*
     .                             unidades/epsilon
                           SS=0.5d0*const3*
     .                          (S_real*S_real + S_imag*S_imag)
                           Eelec_amber=Eelec_amber + SS
                           Erecip_amber=Erecip_amber + SS
                           do i=1,na_mm
                              De=const3*(pc_cos_kr(i)*S_imag
     $                             -pc_sin_kr(i)*S_real)
                              do k=1,3
                                 felec(k,i) = felec(k,i) + De*krecip(k)
                              enddo
                           enddo
                           const5=const3*const4
                           De2=const5*(S_imag*S_imag+S_real*S_real)
                           const7=2.0_dp*(1.0_dp+kmod2*const6)/kmod2
                           do k=1,3
                              do l=1,3
                                 kronij=real(int(((l+k)-abs(l-k))/
     $                                ((l+k)+abs(l-k))),kind=dp)
                                 stress_amber(l,k)=stress_amber(l,k)+
     $                                De2*(kronij-const7*krecip(l)*
     $                                krecip(k))
                              enddo
                           enddo
                        endif
                     endif
                  enddo
               enddo
            enddo
         endif

C     From the reciprocal sum, substract contributions to the electrostatic 
C     interactions already contained in the many-body energy terms (1-,2-,
C     3-body and 1,4-scaled interactions) 
         do i=1,na_mm

C     consider atoms connected by 1 or 2 covalent bonds
            do k=1,nonbondedxat(i)
               j=nonbonded(i,k)
C     Check that the atom i and j are not connected by a dihedral
               do m=1,dihexat(i)
                  if (atdihe(i,m,3)==j) goto  50
               enddo
               dr(1)=ramber(1,i)-ramber(1,j)
               dr(2)=ramber(2,i)-ramber(2,j)
               dr(3)=ramber(3,i)-ramber(3,j)
               call pbc_displ_vector(lattice_type,cell,kcell,
     .                            dr(1),dr(2),dr(3))
               distancia2=dr(1)*dr(1)+dr(2)*dr(2)+dr(3)*dr(3)
               distancia=sqrt(distancia2)
!     Substract the real part of the energy added in the reciprocal-space sum
!     The bond and angle terms have already the contribution from electrostatics
               Eelec_amber=Eelec_amber-(pc(i)*pc(j)/distancia)*
     .              erf(sqrt_ewald_alpha*distancia)*
     .              unidades/epsilon
               Erecip_amber=Erecip_amber-(pc(i)*pc(j)/distancia)*
     .              erf(sqrt_ewald_alpha*distancia)*
     .              unidades/epsilon
               De=pc(i)*pc(j)/distancia2*(erf(sqrt_ewald_alpha
     .              *distancia)/distancia-2.0*const1*
     .              exp(-ewald_alpha*distancia2))*
     .              unidades/epsilon
               do m=1,3
                  felec(m,i)=felec(m,i)+dr(m)*De
                  felec(m,j)=felec(m,j)-dr(m)*De
                  do n=1,3
                     stress_amber(n,m)=stress_amber(n,m)+
     .                    stress_fact*dr(n)*dr(m)*De
                  enddo
               enddo
 50         enddo

C     consider atoms 1 and 4 in covalently-bonded chain 1-2-3-4
            do k=1,scalexat(i)
               j=scale(i,k)
               dr(1)=ramber(1,i)-ramber(1,j)
               dr(2)=ramber(2,i)-ramber(2,j)
               dr(3)=ramber(3,i)-ramber(3,j)
               call pbc_displ_vector(lattice_type,cell,kcell,dr(1),
     .                                  dr(2),dr(3))
               distancia2=dr(1)*dr(1)+dr(2)*dr(2)+dr(3)*dr(3)
               distancia=sqrt(distancia2)
!     Substract the real part of the energy added in the reciprocal-space sum
!     The bond and angle terms have already the contribution from electrostatics
               Eelec_amber14=Eelec_amber14-(pc(i)*pc(j)/distancia)*
     .              erf(sqrt_ewald_alpha*distancia)*
     .              (factorelec-1.0d0)/factorelec*unidades/epsilon
               Erecip_amber=Erecip_amber-(pc(i)*pc(j)/distancia)*
     .              erf(sqrt_ewald_alpha*distancia)*
     .              (factorelec-1.0d0)/factorelec*unidades/epsilon
               fel=pc(i)*pc(j)/distancia2*(erf(sqrt_ewald_alpha
     .              *distancia)/distancia-2.0*const1*
     .              exp(-ewald_alpha*distancia2))*(factorelec-1.0d0)
     .              /factorelec*unidades/epsilon
               do m=1,3
                  felec(m,i)=felec(m,i)+dr(m)*fel
                  felec(m,j)=felec(m,j)-dr(m)*fel
                  do n=1,3
                     stress_amber(n,m)=stress_amber(n,m)+
     .                    stress_fact*dr(n)*dr(m)*De
                  enddo
               enddo
            enddo
         enddo

C     Substract self-energy
         do i=1,na_mm
            Eelec_amber=Eelec_amber-const1*pc(i)*pc(i)*
     .           unidades/epsilon
            Erecip_amber=Erecip_amber-const1*pc(i)*pc(i)*
     .           unidades/epsilon
         enddo

C     end of conditional for ewald sum of reciprocal space
      endif

      end

c******************************************************************
C subroutine that calculates the coulomb interaction using the Ewald method
        subroutine coulomb_ewald_real(pc_i,pc_j,distancia,E2,
     .        fel,ewald_alpha,sqrt_ewald_alpha,const1)

       use precision, only: dp
        real(dp) pc_i, pc_j
        real(dp) distancia, E2, fel
        real(dp) ewald_alpha,sqrt_ewald_alpha,const1
        real(dp) distancia2

        distancia2=distancia*distancia

        E2=((pc_i*pc_j)/distancia)*erfc(sqrt_ewald_alpha*distancia)
        fel=-(pc_i*pc_j)/distancia2*(erfc(sqrt_ewald_alpha
     .   *distancia)/distancia+2.0*const1*exp(-ewald_alpha*distancia2))

        return
        end

c******************************************************************
C subroutine that calculates the coulomb interaction using a cut-off distance
        subroutine coulomb_cutoff(pc_i,pc_j,distancia,distancia2,E2,fel)

        use precision, only: dp
        real(dp) pc_i, pc_j
        real(dp) distancia, distancia2, E2, fel

        E2=((pc_i*pc_j)/distancia)
        fel=-E2/distancia2

        return
        end


c******************************************************************
c subrutina q calcula el water restarin potential

	subroutine waters(na_qm,na_mm,natot,rclas,masst,
     $     noaa,noat,ewat,fwat,cell,kcell)

       use precision, only: dp
        integer na_qm,na_mm,natot,watlist(2000),watlistnum
        real(dp) rclas(3,natot),ewat,rwat,masscenter(3),
     .  rt(3),rij,E,fwat(3,na_mm),dx,dy,dz,masst(natot),
     .  kte,ramber(3,natot),dist,dist2,mdist
        character noat(na_mm)*4,noaa(na_mm)*4
        integer i,j,k,l,m,n
        real(dp) pi
        real(dp)  cell(3,3), kcell(3,3)
        character lattice_type

        pi=DACOS(-1.d0)
        kte=200.0
        dist2=0.0
        mdist=0.0

c     calcula el masscenter del sistema y el rwat
        ramber(1:3,1:natot)=rclas(1:3,1:natot)*0.529177
        masscenter=0.0 
        do i=1,natot
           masscenter(1:3)=masscenter(1:3)+masst(i)*ramber(1:3,i)
        enddo
        masscenter=masscenter/natot
        do i=1,natot
           dx=ramber(1,i)-masscenter(1)
           dy=ramber(2,i)-masscenter(2)
           dz=ramber(3,i)-masscenter(3)
           dist2=dx**2+dy**2+dz**2
           if(dist2.gt.mdist) mdist=dist2
        enddo

	rwat=sqrt(mdist) - 2.5        
	write(6,'(a,2x,f8.4)') 'Water Cutoff Radius:', rwat

c     calculo la matrix con las aguas xa los at MM
        ramber=0.0
        ramber(1:3,1:na_mm)=rclas(1:3,na_qm+1:natot)*0.529177
        k=1
        do i=1,na_mm
           if(noaa(i).eq.'HOH'.and.noat(i).eq.'O') then
              rij=dist(ramber(1,i),ramber(2,i),ramber(3,i),
     .             masscenter(1),masscenter(2),masscenter(3),cell)
              if(rij.gt.rwat) then
                 watlist(k)=i
                 k=k+1
c     fin de si esta en la zona buffer
              endif
c     fin de si es agua
           endif
        enddo
        watlistnum=k-1
        ewat=0.0
        fwat=0.0
c     calculo de la Ene y la fza para el potencial de las aguas Ewat = 0.0
        do j=1,watlistnum
           i = watlist (j)
           rij = dist(ramber (1,i),ramber(2,i),ramber(3,i),
     .          masscenter(1),masscenter(2), masscenter(3),cell)
           if(rij.gt.rwat) then
              E = kte*((rij-rwat)**2)
              ewat = ewat + E
              dx=ramber(1,i)-masscenter(1)
              dx = (1.0/rij)*dx
              dx = 2.*kte*(rij-rwat)*dx
              dy=ramber(2,i)-masscenter(2)
              dy = (1.0/rij)*dy
              dy = 2.*kte*(rij-rwat)*dy
              dz=ramber(3,i)-masscenter(3)
              dz = (1.0/rij)*dz
              dz =  2.*kte*(rij-rwat)*dz
              fwat(1,i)=fwat(1,i)-dx
              fwat(2,i)=fwat(2,i)-dy
              fwat(3,i)=fwat(3,i)-dz
           endif
        enddo
 
        end
c********************************************************************************

