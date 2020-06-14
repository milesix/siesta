c     funciones q utiliza el programa

      function dist(x1,y1,z1,x2,y2,z2,cell)
      use precision, only: dp
      implicit none
      real(dp) dist,x1,y1,z1,x2,y2,z2
      real(dp) cell(3,3)
      real(dp) dx, dy, dz
      dx = x1-x2
C     Count for PBC along X
      dx = dx - ANINT(dx/cell(1,1))*cell(1,1)
      dy = y1-y2
C     Count for PBC along Y
      dy = dy - ANINT(dy/cell(2,2))*cell(2,2)
      dz = z1-z2
C     Count for PBC along Z
      dz = dz - ANINT(dz/cell(3,3))*cell(3,3)
      dist=dx**2+dy**2+dz**2
      dist = sqrt(dist)
      return
      end function dist

      function dist_v2(dx,dy,dz)
      use precision, only: dp
      implicit none
      real(dp) dx, dy, dz
      real(dp) dist_v2
      dist_v2=dx**2+dy**2+dz**2
      dist_v2 = sqrt(dist_v2)
      return
      end function dist_v2

      function dist2(x1,y1,z1,x2,y2,z2,cell)
      use precision, only: dp
      implicit none
      real(dp) dist2,x1,y1,z1,x2,y2,z2
      real(dp) cell(3,3)
      real(dp) dx, dy, dz
      dx = x1-x2
C     Count for PBC along X
      dx = dx - ANINT(dx/cell(1,1))*cell(1,1)
      dy = y1-y2
C     Count for PBC along Y
      dy = dy - ANINT(dy/cell(2,2))*cell(2,2)
      dz = z1-z2
C     Count for PBC along Z
      dz = dz - ANINT(dz/cell(3,3))*cell(3,3)
      dist2=dx**2+dy**2+dz**2
      return
      end function dist2

      function dist2_v2(dx,dy,dz)
      use precision, only: dp
      implicit none
      real(dp) dx, dy, dz
      real(dp) dist2_v2
      dist2_v2=dx**2+dy**2+dz**2
      return
      end function dist2_v2

      function angle(x1,y1,z1,x2,y2,z2,x3,y3,z3,cell)
      use precision, only: dp
      implicit none
      real(dp) scalar,dist,angle,pi
      real(dp) x1,y1,z1,x2,y2,z2,x3,y3,z3
      real(dp) cell(3,3)
      real(dp) dx12, dy12, dz12, dx32, dy32, dz32

      pi=DACOS(-1.d0)
c     calculo del producto escalar
      dx12 = x1-x2
C     Count for PBC along X
      dx12 = dx12 - ANINT(dx12/cell(1,1))*cell(1,1)
      dy12 = y1-y2
C     Count for PBC along Y
      dy12 = dy12 - ANINT(dy12/cell(2,2))*cell(2,2)
      dz12 = z1-z2
C     Count for PBC along Z
      dz12 = dz12 - ANINT(dz12/cell(3,3))*cell(3,3)
c     calculo del producto escalar
      dx32 = x3-x2
C     Count for PBC along X
      dx32 = dx32 - ANINT(dx32/cell(1,1))*cell(1,1)
      dy32 = y3-y2
C     Count for PBC along Y
      dy32 = dy32 - ANINT(dy32/cell(2,2))*cell(2,2)
      dz32 = z3-z2
C     Count for PBC along Z
      dz32 = dz32 - ANINT(dz32/cell(3,3))*cell(3,3)
      scalar = dx12*dx32+dy12*dy32+dz12*dz32
      angle = dist(x1,y1,z1,x2,y2,z2,cell)*
     $     dist(x3,y3,z3,x2,y2,z2,cell)
      angle = scalar/angle
      if (angle.ge.1.0) then
         angle = 1.0
      elseif(angle.le.-1.0) then
         angle=-1.0
      endif
      angle = dACOS(angle)*180/pi 

      end function angle

      function angle_v2(dx12,dy12,dz12,dx32,dy32,dz32)
      use precision, only: dp
      implicit none
      real(dp) scalar,dist_v2,pi
      real(dp) dx12,dy12,dz12,dx32,dy32,dz32
      real(dp) angle_v2

      pi=DACOS(-1.d0)
      scalar = dx12*dx32+dy12*dy32+dz12*dz32
      angle_v2 = dist_v2(dx12,dy12,dz12)*dist_v2(dx32,dy32,dz32)
      angle_v2 = scalar/angle_v2
      if (angle_v2.ge.1.0) then
         angle_v2 = 1.0
      elseif(angle_v2.le.-1.0) then
         angle_v2=-1.0
      endif
      angle_v2 = dACOS(angle_v2)*180/pi 

      end function angle_v2

      function scalar(x1,y1,z1,x2,y2,z2,x3,y3,z3,cell)
      use precision, only: dp
      implicit none
      real(dp) scalar,pi
      real(dp) x1,y1,z1,x2,y2,z2,x3,y3,z3
      real(dp) cell(3,3)
      real(dp) dx12, dy12, dz12, dx32, dy32, dz32
c     calculo del producto escalar
      dx12 = x1-x2
C     Count for PBC along X
      dx12 = dx12 - ANINT(dx12/cell(1,1))*cell(1,1)
      dy12 = y1-y2
C     Count for PBC along Y
      dy12 = dy12 - ANINT(dy12/cell(2,2))*cell(2,2)
      dz12 = z1-z2
C     Count for PBC along Z
      dz12 = dz12 - ANINT(dz12/cell(3,3))*cell(3,3)
c     calculo del producto escalar
      dx32 = x3-x2
C     Count for PBC along X
      dx32 = dx32 - ANINT(dx32/cell(1,1))*cell(1,1)
      dy32 = y3-y2
C     Count for PBC along Y
      dy32 = dy32 - ANINT(dy32/cell(2,2))*cell(2,2)
      dz32 = z3-z2
C     Count for PBC along Z
      dz32 = dz32 - ANINT(dz32/cell(3,3))*cell(3,3)
      pi=DACOS(-1.d0)
c     calculo del producto escalar
      scalar = dx12*dx32+dy12*dy32+dz12*dz32
      end function scalar

      function scalar_v2(x1,y1,z1,x2,y2,z2)
      use precision, only: dp
      implicit none
      real(dp) x1, y1, z1,x2,y2,z2
      real(dp) scalar_v2
c     calculo del producto escalar
      scalar_v2=x1*x2+y1*y2+z1*z2
      return
      end function scalar_v2

      function dihedro(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,cell)
      use precision, only: dp
      implicit none
      real(dp) dihedro,dist,angle,l1,l2,pi,arg
      real(dp) x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4
      real(dp) l,A,B,C,D,m,n,nx,ny,nz,mx,my,mz,scalar
      real(dp) cell(3,3)
      real(dp) dx12, dy12, dz12, dx32, dy32, dz32
      real(dp) dx23, dy23, dz23, dx43, dy43, dz43
      pi=DACOS(-1.d0)

      dx12 = x1-x2
C     Count for PBC along X
      dx12 = dx12 - ANINT(dx12/cell(1,1))*cell(1,1)
      dy12 = y1-y2
C     Count for PBC along Y
      dy12 = dy12 - ANINT(dy12/cell(2,2))*cell(2,2)
      dz12 = z1-z2
C     Count for PBC along Z
      dz12 = dz12 - ANINT(dz12/cell(3,3))*cell(3,3)
      dx32 = x3-x2
C     Count for PBC along X
      dx32 = dx32 - ANINT(dx32/cell(1,1))*cell(1,1)
      dy32 = y3-y2
C     Count for PBC along Y
      dy32 = dy32 - ANINT(dy32/cell(2,2))*cell(2,2)
      dz32 = z3-z2
C     Count for PBC along Z
      dz32 = dz32 - ANINT(dz32/cell(3,3))*cell(3,3)
      mx = dy12*dz32 - dz12*dy32
      my = -(dx12*dz32 - dz12*dx32)
      mz = dx12*dy32-dy12*dx32

      dx23 = x2-x3
C     Count for PBC along X
      dx23 = dx23 - ANINT(dx23/cell(1,1))*cell(1,1)
      dy23 = y2-y3
C     Count for PBC along Y
      dy23 = dy23 - ANINT(dy23/cell(2,2))*cell(2,2)
      dz23 = z2-z3
C     Count for PBC along Z
      dz23 = dz23 - ANINT(dz23/cell(3,3))*cell(3,3)
      dx43 = x4-x3
C     Count for PBC along X
      dx43 = dx43 - ANINT(dx43/cell(1,1))*cell(1,1)
      dy43 = y4-y3
C     Count for PBC along Y
      dy43 = dy43 - ANINT(dy43/cell(2,2))*cell(2,2)
      dz43 = z4-z3
C     Count for PBC along Z
      dz43 = dz43 - ANINT(dz43/cell(3,3))*cell(3,3)
      nx = dy23*dz43 - dz23*dy43
      ny = -(dx23*dz43 - dz23*dx43)
      nz = dx23*dy43-dy23*dx43
c     calculo del prod scalar n*m
      scalar = mx*nx + my*ny + mz*nz
      m = mx**2 + my**2 + mz**2
      m = m**(0.5)
      n = nx**2 + ny**2 + nz**2
      n = n**(0.5)
c     si el argumento da mal (m*n)=0.0 avisa y sale el dihe vale 500.0
      if(m*n.eq.0) then
         dihedro=500.0
         go to 10
      endif
      arg=scalar/(m*n)
      if(arg.ge.1.0) then
         dihedro= 0.0
      elseif (arg.le.-1.0) then
         dihedro=180.0
      else
         dihedro = dACOS(arg)
      endif
      dihedro = dihedro*180/pi
 10   continue
      end function dihedro

      function dihedro_v2(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,cell,
     .kcell,lattice_type)
      use precision, only: dp
      implicit none
      real(dp) dihedro_v2,dist,angle,l1,l2,pi,arg
      real(dp) l,A,B,C,D,m,n,nx,ny,nz,mx,my,mz,scalar
      real(dp) x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4
      real(dp) dx12, dy12, dz12, dx32, dy32, dz32
      real(dp) dx43, dy43, dz43
      real(dp) cell(3,3),kcell(3,3)
      character lattice_type
      pi=DACOS(-1.d0)

      dx12=x1-x2
      dy12=y1-y2
      dz12=z1-z2
      call pbc_displ_vector(lattice_type,cell,kcell,dx12,
     .     dy12,dz12)

      dx32=x3-x2
      dy32=y3-y2
      dz32=z3-z2
      call pbc_displ_vector(lattice_type,cell,kcell,dx32,
     .     dy32,dz32)
      mx = dy12*dz32 - dz12*dy32
      my = -(dx12*dz32 - dz12*dx32)
      mz = dx12*dy32-dy12*dx32

      dx43=x4-x3
      dy43=y4-y3
      dz43=z4-z3
      call pbc_displ_vector(lattice_type,cell,kcell,dx43,
     .     dy43,dz43)
      nx = -dy32*dz43 + dz32*dy43
      ny = dx32*dz43 - dz32*dx43
      nz = -dx32*dy43+dy32*dx43

c     calculo del prod scalar n*m
      scalar = mx*nx + my*ny + mz*nz
      m = mx**2 + my**2 + mz**2
      m = m**(0.5)
      n = nx**2 + ny**2 + nz**2
      n = n**(0.5)
c     si el argumento da mal (m*n)=0.0 avisa y sale el dihe vale 500.0
      if(m*n.eq.0) then
         dihedro_v2=500.0
         go to 10
      endif
      arg=scalar/(m*n)
      if(arg.ge.1.0) then
         dihedro_v2= 0.0
      elseif (arg.le.-1.0) then
         dihedro_v2=180.0
      else
         dihedro_v2 = dACOS(arg)
      endif
      dihedro_v2 = dihedro_v2*180/pi
 10   continue
      end function dihedro_v2

      function dihedro2(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,cell)
      use precision, only: dp
      use sys, only : die
      implicit none
      real(dp) dihedro2,dist,angle,l1,l2,pi
      real(dp) x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4
      real(dp) l,A,B,C,D,m,n,nx,ny,nz,mx,my,mz,scalar
      real(dp) cell(3,3)
      real(dp) dx12, dy12, dz12, dx32, dy32, dz32
      real(dp) dx23, dy23, dz23, dx43, dy43, dz43
      real(dp) dx21,dy21, dz21, dx31, dy31, dz31
      real(dp) cos_alpha
      
c     vectors n and m
c     m = vectorial product 21 and 23
c     n = vectorial product 23 and 34
      
      pi=DACOS(-1.d0)
c     which cuadrant? 
c     plane generation with points 1 2 3
      dx12 = x1-x2
C     Count for PBC along X
      dx12 = dx12 - ANINT(dx12/cell(1,1))*cell(1,1)
      dy12 = y1-y2
C     Count for PBC along Y
      dy12 = dy12 - ANINT(dy12/cell(2,2))*cell(2,2)
      dz12 = z1-z2
C     Count for PBC along Z
      dz12 = dz12 - ANINT(dz12/cell(3,3))*cell(3,3)
      dx32 = x3-x2
C     Count for PBC along X
      dx32 = dx32 - ANINT(dx32/cell(1,1))*cell(1,1)
      dy32 = y3-y2
C     Count for PBC along Y
      dy32 = dy32 - ANINT(dy32/cell(2,2))*cell(2,2)
      dz32 = z3-z2
C     Count for PBC along Z
      dz32 = dz32 - ANINT(dz32/cell(3,3))*cell(3,3)
      mx = dy12*dz32 - dz12*dy32
      my = -(dx12*dz32 - dz12*dx32)
      mz = dx12*dy32-dy12*dx32
      

c     plane generation with points 1 2 3
      dx23 = x2-x3
C     Count for PBC along X
      dx23 = dx23 - ANINT(dx23/cell(1,1))*cell(1,1)
      dy23 = y2-y3
C     Count for PBC along Y
      dy23 = dy23 - ANINT(dy23/cell(2,2))*cell(2,2)
      dz23 = z2-z3
C     Count for PBC along Z
      dz23 = dz23 - ANINT(dz23/cell(3,3))*cell(3,3)
      dx43 = x4-x3
C     Count for PBC along X
      dx43 = dx43 - ANINT(dx43/cell(1,1))*cell(1,1)
      dy43 = y4-y3
C     Count for PBC along Y
      dy43 = dy43 - ANINT(dy43/cell(2,2))*cell(2,2)
      dz43 = z4-z3
C     Count for PBC along Z
      dz43 = dz43 - ANINT(dz43/cell(3,3))*cell(3,3)
      nx = dy23*dz43 - dz23*dy43
      ny = -(dx23*dz43 - dz23*dx43)
      nz = dx23*dy43-dy23*dx43
      
c     scalar product n*m
      scalar = mx*nx + my*ny + mz*nz
      m = mx**2 + my**2 + mz**2
      m = m**(0.5)
      n = nx**2 + ny**2 + nz**2
      n = n**(0.5)
      cos_alpha=scalar/(m*n)
      if (abs(cos_alpha)>1.00001) then 
         call die('dihedro2: abs(cos_alpha)>1.00001')
      else if (abs(cos_alpha)>1.0d0) then
         if (cos_alpha<0.0d0) then
            cos_alpha=-1.0d0 
         else
            cos_alpha=1.0d0
         endif
      endif
      dihedro2 = ACOS(cos_alpha)
      dihedro2 = dihedro2*180/pi
      
c     which cuadrant? 
c     plane generation with points 1 2 3
      dx21 = x2-x1
C     Count for PBC along X
      dx21 = dx21 - ANINT(dx21/cell(1,1))*cell(1,1)
      dy21 = y2-y1
C     Count for PBC along Y
      dy21 = dy21 - ANINT(dy21/cell(2,2))*cell(2,2)
      dz21 = z2-z1
C     Count for PBC along Z
      dz21 = dz21 - ANINT(dz21/cell(3,3))*cell(3,3)
      dx31 = x3-x1
C     Count for PBC along X
      dx31 = dx31 - ANINT(dx31/cell(1,1))*cell(1,1)
      dy31 = y3-y1
C     Count for PBC along Y
      dy31 = dy31 - ANINT(dy31/cell(2,2))*cell(2,2)
      dz31 = z3-z1
C     Count for PBC along Z
      dz31 = dz31 - ANINT(dz31/cell(3,3))*cell(3,3)
c     calculo del producto escalar
      A = dy21*dz31-dz21*dy31
      B = dz21*dx31-dx21*dz31
      C = dx21*dy31-dy21*dx31
      D = -x1*A-y1*B-z1*C
      
c     distance(l) at4 to ABCD=0 plane
      l1 = A*x4+B*y4+C*z4+D
      l2 = (A**2+B**2+C**2)
      l2 = l2**0.5
      l= l1/l2
      
c     if l>0 -> dihe<0 , else >0
      if(l.lt.0) then
         dihedro2 = 360-dihedro2
      endif
      
      end function dihedro2

      function dihedro2_v2(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,cell,
     . kcell,lattice_type)
      use precision, only: dp
      use sys, only : die
      implicit none
      real(dp) dihedro2_v2,dist,angle,l1,l2,pi
      real(dp) l,A,B,C,D,m,n,nx,ny,nz,mx,my,mz,scalar
      real(dp) x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4
      real(dp) dx12, dy12, dz12, dx32, dy32, dz32
      real(dp) dx31, dy31, dz31, dx43, dy43, dz43 
      real(dp) cos_alpha
      real(dp) cell(3,3),kcell(3,3)
      character lattice_type
      
c     vectors n and m
c     m = vectorial product 21 and 23
c     n = vectorial product 23 and 34
      
      pi=DACOS(-1.d0)
c     which cuadrant? 
      dx12=x1-x2
      dy12=y1-y2
      dz12=z1-z2
      call pbc_displ_vector(lattice_type,cell,kcell,dx12,
     .     dy12,dz12)

      dx32=x3-x2
      dy32=y3-y2
      dz32=z3-z2
      call pbc_displ_vector(lattice_type,cell,kcell,dx32,
     .     dy32,dz32)
c     plane generation with points 1 2 3
      mx = dy12*dz32 - dz12*dy32
      my = -(dx12*dz32 - dz12*dx32)
      mz = dx12*dy32-dy12*dx32      

      dx43=x4-x3
      dy43=y4-y3
      dz43=z4-z3
      call pbc_displ_vector(lattice_type,cell,kcell,dx43,
     .     dy43,dz43)
c     plane generation with points 1 2 3
      nx = -dy32*dz43 + dz32*dy43
      ny = dx32*dz43 - dz32*dx43
      nz = -dx32*dy43+dy32*dx43
      
c     scalar product n*m
      scalar = mx*nx + my*ny + mz*nz
      m = mx**2 + my**2 + mz**2
      m = m**(0.5)
      n = nx**2 + ny**2 + nz**2
      n = n**(0.5)
      cos_alpha=scalar/(m*n)
      if (abs(cos_alpha)>1.00001) then 
         call die('dihedro2_v2: abs(cos_alpha)>1.00001')
      else if (abs(cos_alpha)>1.0d0) then
         if (cos_alpha<0.0d0) then
            cos_alpha=-1.0d0 
         else
            cos_alpha=1.0d0
         endif
      endif
      dihedro2_v2 = ACOS(cos_alpha)
      dihedro2_v2 = dihedro2_v2*180/pi
      
c     which cuadrant? 
c     plane generation with points 1 2 3
c     calculo del producto escalar
      A = -dy12*dz31+dz12*dy31
      B = -dz12*dx31+dx12*dz31
      C = -dx12*dy31+dy12*dx31
      D = -x1*A-y1*B-z1*C
      
c     distance(l) at4 to ABCD=0 plane
      l1 = A*x4+B*y4+C*z4+D
      l2 = (A**2+B**2+C**2)
      l2 = l2**0.5
      l= l1/l2
      
c     if l>0 -> dihe<0 , else >0
      if(l.lt.0) then
         dihedro2_v2 = 360-dihedro2_v2
      endif
      
      end function dihedro2_v2

c************************************************************************************
c     subroutine that calculates the force of a dihedral angle in amber

      subroutine diheforce(nac,ramber,i1,i2,i3,i4,
     .     atom,kd,eq,per,mult,fce,cell)
      use precision, only: dp
c     variables asoc a funcion como subrutina
      integer i1,i2,i3,i4
      real(dp) ramber(3,nac),fce(12),dtot
c     i1-i4 numero de atomos 1 a 4  atom es que derivada tiene que calcular 
c     atom=1 indica el primer atomo del dihe atomo=2 el segundo etc...
c     fce(12) es la fuerza (derivada) para at1x,at1y,at1z,at2x....at4z
      integer  mult
      real(dp)  kd,eq,per
c     variables asoc al calculo de la derivada
      real(dp) dih,dihedro,dist,rm,rn,mx,my,mz,nx,ny,nz,scal
      real(dp) dscalar,dm,dn,dmx,dmy,dmz,dnx,dny,dnz,dmn,
     .     fx,fy,fz,E(3),fpx,fpxmasd,fpxmenosd,fr,step,prue
c     variables generales
      integer i,j,k,l,m,n,nac,vez,coord,atom
      character exp
      real(dp) pi
      real(dp) cell(3,3)
      pi=DACOS(-1.d0)
      call dihevars(   ramber(1,i1),ramber(2,i1),ramber(3,i1),
     .     ramber(1,i2),ramber(2,i2),ramber(3,i2),
     .     ramber(1,i3),ramber(2,i3),ramber(3,i3),
     .     ramber(1,i4),ramber(2,i4),ramber(3,i4),
     .     mx,my,mz,rm,nx,ny,nz,rn,dih,cell)
      fce=0.0       
      scal=mx*nx + my*ny + mz*nz    
      dtot = (kd/mult)*(-dSIN((pi/180)*
     .     (per*dih-eq)))*(per)
      prue=scal/(rn*rm)
      prue=(1.0-(prue)**2) 
      if (prue.lt.1.0E-15.and.prue.gt.-1.0E-15) go to 10  
      prue=dsqrt(prue) 
      dtot = -dtot/prue
      do j=1,3
         i=(atom-1)*3+j
         dmx=0.0
         dmy=0.0
         dmz=0.0
         dnx=0.0
         dny=0.0
         dnz=0.0
         if(i.eq.1) then
            dmy=ramber(3,i2)-ramber(3,i3)
            dmz=ramber(2,i3)-ramber(2,i2)
         elseif(i.eq.4) then
            dmy=ramber(3,i3)-ramber(3,i1)
            dmz=ramber(2,i1)-ramber(2,i3)
            dny=ramber(3,i3)-ramber(3,i4)
            dnz=ramber(2,i4)-ramber(2,i3)
         elseif(i.eq.7) then
            dmy=ramber(3,i1)-ramber(3,i2)
            dmz=ramber(2,i2)-ramber(2,i1)
            dny=ramber(3,i4)-ramber(3,i2)
            dnz=ramber(2,i2)-ramber(2,i4)
         elseif(i.eq.10) then
            dny=ramber(3,i2)-ramber(3,i3)
            dnz=ramber(2,i3)-ramber(2,i2)
         elseif(i.eq.2) then
            dmx=ramber(3,i3)-ramber(3,i2)
            dmz=ramber(1,i2)-ramber(1,i3)
         elseif(i.eq.5) then
            dmx=ramber(3,i1)-ramber(3,i3)
            dmz=ramber(1,i3)-ramber(1,i1)
            dnx=ramber(3,i4)-ramber(3,i3)
            dnz=ramber(1,i3)-ramber(1,i4)
         elseif(i.eq.8) then
            dmx=ramber(3,i2)-ramber(3,i1)
            dmz=ramber(1,i1)-ramber(1,i2)
            dnx=ramber(3,i2)-ramber(3,i4)
            dnz=ramber(1,i4)-ramber(1,i2)
         elseif(i.eq.11) then
            dnx=ramber(3,i3)-ramber(3,i2)
            dnz=ramber(1,i2)-ramber(1,i3)
         elseif(i.eq.3) then
            dmx=ramber(2,i2)-ramber(2,i3)
            dmy=ramber(1,i3)-ramber(1,i2)
         elseif(i.eq.6) then
            dmx=ramber(2,i3)-ramber(2,i1)
            dmy=ramber(1,i1)-ramber(1,i3)
            dnx=ramber(2,i3)-ramber(2,i4)
            dny=ramber(1,i4)-ramber(1,i3)
         elseif(i.eq.9) then
            dmx=ramber(2,i1)-ramber(2,i2)
            dmy=ramber(1,i2)-ramber(1,i1)
            dnx=ramber(2,i4)-ramber(2,i2)
            dny=ramber(1,i2)-ramber(1,i4)
         elseif(i.eq.12) then
            dnx=ramber(2,i2)-ramber(2,i3)
            dny=ramber(1,i3)-ramber(1,i2)
         endif
         dmx=dmx-ANINT(dmx/cell(1,1))*dmx
         dmy=dmy-ANINT(dmy/cell(2,2))*dmy
         dmz=dmz-ANINT(dmz/cell(3,3))*dmz
         dnx=dnx-ANINT(dnx/cell(1,1))*dnx
         dny=dny-ANINT(dny/cell(2,2))*dny
         dnz=dnz-ANINT(dnz/cell(3,3))*dnz
         dm=(mx*dmx+my*dmy+mz*dmz)/rm
         dn=(nx*dnx+ny*dny+nz*dnz)/rn
         dmn=rm*dn+rn*dm
         dscalar=nx*dmx+mx*dnx+ny*dmy+my*dny+nz*dmz+mz*dnz
         fce(i)= dtot*(dscalar*rm*rn-dmn*scal)/(rn*rm)**2
      enddo
 10   end

c     subroutine that calculates the force of a dihedral angle in amber

      subroutine diheforce_v2(nac,ramber,i1,i2,i3,i4,
     .     atom,kd,eq,per,mult,fce,cell,kcell,lattice_type)
      use precision, only: dp
c     variables asoc a funcion como subrutina
      integer i1,i2,i3,i4
      real(dp) ramber(3,nac),fce(12),dtot
c     i1-i4 numero de atomos 1 a 4  atom es que derivada tiene que calcular 
c     atom=1 indica el primer atomo del dihe atomo=2 el segundo etc...
c     fce(12) es la fuerza (derivada) para at1x,at1y,at1z,at2x....at4z
      integer  mult
      real(dp)  kd,eq,per
c     variables asoc al calculo de la derivada
      real(dp) dih,dihedro,dist,rm,rn,mx,my,mz,nx,ny,nz,scal
      real(dp) dscalar,dm,dn,dmx,dmy,dmz,dnx,dny,dnz,dmn,
     .     fx,fy,fz,E(3),fpx,fpxmasd,fpxmenosd,fr,step,prue
c     variables generales
      integer i,j,k,l,m,n,nac,vez,coord,atom
      character exp
      real(dp) pi
      real(dp) cell(3,3),kcell(3,3)
      character lattice_type
      pi=DACOS(-1.d0)

      call dihevars_v2(   ramber(1,i1),ramber(2,i1),ramber(3,i1),
     .     ramber(1,i2),ramber(2,i2),ramber(3,i2),
     .     ramber(1,i3),ramber(2,i3),ramber(3,i3),
     .     ramber(1,i4),ramber(2,i4),ramber(3,i4),
     .     mx,my,mz,rm,nx,ny,nz,rn,dih,cell,kcell,lattice_type)

      fce=0.0       
      scal=mx*nx + my*ny + mz*nz    
      dtot = (kd/mult)*(-dSIN((pi/180)*
     .     (per*dih-eq)))*(per)
      prue=scal/(rn*rm)
      prue=(1.0-(prue)**2) 
      if (prue.lt.1.0E-15.and.prue.gt.-1.0E-15) go to 10  
      prue=dsqrt(prue) 
      dtot = -dtot/prue
      do j=1,3
         i=(atom-1)*3+j
         dmx=0.0
         dmy=0.0
         dmz=0.0
         dnx=0.0
         dny=0.0
         dnz=0.0
         if(i.eq.1) then
            dmy=ramber(3,i2)-ramber(3,i3)
            dmz=ramber(2,i3)-ramber(2,i2)
         elseif(i.eq.4) then
            dmy=ramber(3,i3)-ramber(3,i1)
            dmz=ramber(2,i1)-ramber(2,i3)
            dny=ramber(3,i3)-ramber(3,i4)
            dnz=ramber(2,i4)-ramber(2,i3)
         elseif(i.eq.7) then
            dmy=ramber(3,i1)-ramber(3,i2)
            dmz=ramber(2,i2)-ramber(2,i1)
            dny=ramber(3,i4)-ramber(3,i2)
            dnz=ramber(2,i2)-ramber(2,i4)
         elseif(i.eq.10) then
            dny=ramber(3,i2)-ramber(3,i3)
            dnz=ramber(2,i3)-ramber(2,i2)
         elseif(i.eq.2) then
            dmx=ramber(3,i3)-ramber(3,i2)
            dmz=ramber(1,i2)-ramber(1,i3)
         elseif(i.eq.5) then
            dmx=ramber(3,i1)-ramber(3,i3)
            dmz=ramber(1,i3)-ramber(1,i1)
            dnx=ramber(3,i4)-ramber(3,i3)
            dnz=ramber(1,i3)-ramber(1,i4)
         elseif(i.eq.8) then
            dmx=ramber(3,i2)-ramber(3,i1)
            dmz=ramber(1,i1)-ramber(1,i2)
            dnx=ramber(3,i2)-ramber(3,i4)
            dnz=ramber(1,i4)-ramber(1,i2)
         elseif(i.eq.11) then
            dnx=ramber(3,i3)-ramber(3,i2)
            dnz=ramber(1,i2)-ramber(1,i3)
         elseif(i.eq.3) then
            dmx=ramber(2,i2)-ramber(2,i3)
            dmy=ramber(1,i3)-ramber(1,i2)
         elseif(i.eq.6) then
            dmx=ramber(2,i3)-ramber(2,i1)
            dmy=ramber(1,i1)-ramber(1,i3)
            dnx=ramber(2,i3)-ramber(2,i4)
            dny=ramber(1,i4)-ramber(1,i3)
         elseif(i.eq.9) then
            dmx=ramber(2,i1)-ramber(2,i2)
            dmy=ramber(1,i2)-ramber(1,i1)
            dnx=ramber(2,i4)-ramber(2,i2)
            dny=ramber(1,i2)-ramber(1,i4)
         elseif(i.eq.12) then
            dnx=ramber(2,i2)-ramber(2,i3)
            dny=ramber(1,i3)-ramber(1,i2)
         endif
         call pbc_displ_vector(lattice_type,cell,kcell,dmx,
     .     dmy,dmz)
         call pbc_displ_vector(lattice_type,cell,kcell,dnx,
     .     dny,dnz)
         dm=(mx*dmx+my*dmy+mz*dmz)/rm
         dn=(nx*dnx+ny*dny+nz*dnz)/rn
         dmn=rm*dn+rn*dm
         dscalar=nx*dmx+mx*dnx+ny*dmy+my*dny+nz*dmz+mz*dnz
         fce(i)= dtot*(dscalar*rm*rn-dmn*scal)/(rn*rm)**2
      enddo
 10   end

c*************************************************************************
c     subroutine that calculates the force of a dihedral angle in general

      subroutine diheforce2(nac,ramber,i1,i2,i3,i4,
     .     atom,kd,fce,cell)

      use precision, only: dp      
c     variables asoc a funcion como subrutina
      integer i1,i2,i3,i4
      real(dp) ramber(3,nac),fce(12),dtot
c     i1-i4 numero de atomos 1 a 4  atom es que derivada tiene que calcular
c     atom=1 indica el primer atomo del dihe atomo=2 el segundo etc...
c     fce(12) es la fuerza (derivada) para at1x,at1y,at1z,at2x....at4z
      integer  mult
      real(dp)  kd,eq,per
c     variables asoc al calculo de la derivada
      real(dp) dih,dihedro,dist,rm,rn,mx,my,mz,nx,ny,nz,scal
      real(dp) dscalar,dm,dn,dmx,dmy,dmz,dnx,dny,dnz,dmn,
     .     fx,fy,fz,E(3),fpx,fpxmasd,fpxmenosd,fr,step,prue
c     variables generales
      integer i,j,k,l,m,n,nac,vez,coord,atom
      character exp
      real(dp) pi
      real(dp) cell(3,3)
      pi=DACOS(-1.d0)
      call dihevars(   ramber(1,i1),ramber(2,i1),ramber(3,i1),
     .     ramber(1,i2),ramber(2,i2),ramber(3,i2),
     .     ramber(1,i3),ramber(2,i3),ramber(3,i3),
     .     ramber(1,i4),ramber(2,i4),ramber(3,i4),
     .     mx,my,mz,rm,nx,ny,nz,rn,dih,cell)
      fce=0.0
      scal=mx*nx + my*ny + mz*nz
      dtot=kd
      prue=scal/(rn*rm)
      prue=(1.0-(prue)**2)
      if (prue.lt.1.0E-15.and.prue.gt.-1.0E-15) go to 10
      prue=dsqrt(prue)
      dtot = -dtot/prue
      do j=1,3
         i=(atom-1)*3+j
         dmx=0.0
         dmy=0.0
         dmz=0.0
         dnx=0.0
         dny=0.0
         dnz=0.0
         if(i.eq.1) then
            dmy=ramber(3,i2)-ramber(3,i3)
            dmz=ramber(2,i3)-ramber(2,i2)
         elseif(i.eq.4) then
            dmy=ramber(3,i3)-ramber(3,i1)
            dmz=ramber(2,i1)-ramber(2,i3)
            dny=ramber(3,i3)-ramber(3,i4)
            dnz=ramber(2,i4)-ramber(2,i3)
         elseif(i.eq.7) then
            dmy=ramber(3,i1)-ramber(3,i2)
            dmz=ramber(2,i2)-ramber(2,i1)
            dny=ramber(3,i4)-ramber(3,i2)
            dnz=ramber(2,i2)-ramber(2,i4)
         elseif(i.eq.10) then
            dny=ramber(3,i2)-ramber(3,i3)
            dnz=ramber(2,i3)-ramber(2,i2)
         elseif(i.eq.2) then
            dmx=ramber(3,i3)-ramber(3,i2)
            dmz=ramber(1,i2)-ramber(1,i3)
         elseif(i.eq.5) then
            dmx=ramber(3,i1)-ramber(3,i3)
            dmz=ramber(1,i3)-ramber(1,i1)
            dnx=ramber(3,i4)-ramber(3,i3)
            dnz=ramber(1,i3)-ramber(1,i4)
         elseif(i.eq.8) then
            dmx=ramber(3,i2)-ramber(3,i1)
            dmz=ramber(1,i1)-ramber(1,i2)
            dnx=ramber(3,i2)-ramber(3,i4)
            dnz=ramber(1,i4)-ramber(1,i2)
         elseif(i.eq.11) then
            dnx=ramber(3,i3)-ramber(3,i2)
            dnz=ramber(1,i2)-ramber(1,i3)
         elseif(i.eq.3) then
            dmx=ramber(2,i2)-ramber(2,i3)
            dmy=ramber(1,i3)-ramber(1,i2)
         elseif(i.eq.6) then
            dmx=ramber(2,i3)-ramber(2,i1)
            dmy=ramber(1,i1)-ramber(1,i3)
            dnx=ramber(2,i3)-ramber(2,i4)
            dny=ramber(1,i4)-ramber(1,i3)
         elseif(i.eq.9) then
            dmx=ramber(2,i1)-ramber(2,i2)
            dmy=ramber(1,i2)-ramber(1,i1)
            dnx=ramber(2,i4)-ramber(2,i2)
            dny=ramber(1,i2)-ramber(1,i4)
         elseif(i.eq.12) then
            dnx=ramber(2,i2)-ramber(2,i3)
            dny=ramber(1,i3)-ramber(1,i2)
         endif
         dmx=dmx-ANINT(dmx/cell(1,1))*dmx
         dmy=dmy-ANINT(dmy/cell(2,2))*dmy
         dmz=dmz-ANINT(dmz/cell(3,3))*dmz
         dnx=dnx-ANINT(dnx/cell(1,1))*dnx
         dny=dny-ANINT(dny/cell(2,2))*dny
         dnz=dnz-ANINT(dnz/cell(3,3))*dnz
         dm=(mx*dmx+my*dmy+mz*dmz)/rm
         dn=(nx*dnx+ny*dny+nz*dnz)/rn
         dmn=rm*dn+rn*dm
         dscalar=nx*dmx+mx*dnx+ny*dmy+my*dny+nz*dmz+mz*dnz
         fce(i)= dtot*(dscalar*rm*rn-dmn*scal)/(rn*rm)**2
      enddo
 10   end

c*******************************************************************
      subroutine diheforce2_v2(nac,ramber,i1,i2,i3,i4,
     .     atom,kd,fce,cell,kcell,lattice_type)

      use precision, only: dp      
c     variables asoc a funcion como subrutina
      integer i1,i2,i3,i4
      real(dp) ramber(3,nac),fce(12),dtot
c     i1-i4 numero de atomos 1 a 4  atom es que derivada tiene que calcular
c     atom=1 indica el primer atomo del dihe atomo=2 el segundo etc...
c     fce(12) es la fuerza (derivada) para at1x,at1y,at1z,at2x....at4z
      integer  mult
      real(dp)  kd,eq,per
c     variables asoc al calculo de la derivada
      real(dp) dih,dihedro,dist,rm,rn,mx,my,mz,nx,ny,nz,scal
      real(dp) dscalar,dm,dn,dmx,dmy,dmz,dnx,dny,dnz,dmn,
     .     fx,fy,fz,E(3),fpx,fpxmasd,fpxmenosd,fr,step,prue
c     variables generales
      integer i,j,k,l,m,n,nac,vez,coord,atom
      character exp
      real(dp) pi
      real(dp) cell(3,3),kcell(3,3)
      character lattice_type
      pi=DACOS(-1.d0)
      call dihevars_v2(   ramber(1,i1),ramber(2,i1),ramber(3,i1),
     .     ramber(1,i2),ramber(2,i2),ramber(3,i2),
     .     ramber(1,i3),ramber(2,i3),ramber(3,i3),
     .     ramber(1,i4),ramber(2,i4),ramber(3,i4),
     .     mx,my,mz,rm,nx,ny,nz,rn,dih,cell,kcell,lattice_type)
      fce=0.0
      scal=mx*nx + my*ny + mz*nz
      dtot=kd
      prue=scal/(rn*rm)
      prue=(1.0-(prue)**2)
      if (prue.lt.1.0E-15.and.prue.gt.-1.0E-15) go to 10
      prue=dsqrt(prue)
      dtot = -dtot/prue
      do j=1,3
         i=(atom-1)*3+j
         dmx=0.0
         dmy=0.0
         dmz=0.0
         dnx=0.0
         dny=0.0
         dnz=0.0
         if(i.eq.1) then
            dmy=ramber(3,i2)-ramber(3,i3)
            dmz=ramber(2,i3)-ramber(2,i2)
         elseif(i.eq.4) then
            dmy=ramber(3,i3)-ramber(3,i1)
            dmz=ramber(2,i1)-ramber(2,i3)
            dny=ramber(3,i3)-ramber(3,i4)
            dnz=ramber(2,i4)-ramber(2,i3)
         elseif(i.eq.7) then
            dmy=ramber(3,i1)-ramber(3,i2)
            dmz=ramber(2,i2)-ramber(2,i1)
            dny=ramber(3,i4)-ramber(3,i2)
            dnz=ramber(2,i2)-ramber(2,i4)
         elseif(i.eq.10) then
            dny=ramber(3,i2)-ramber(3,i3)
            dnz=ramber(2,i3)-ramber(2,i2)
         elseif(i.eq.2) then
            dmx=ramber(3,i3)-ramber(3,i2)
            dmz=ramber(1,i2)-ramber(1,i3)
         elseif(i.eq.5) then
            dmx=ramber(3,i1)-ramber(3,i3)
            dmz=ramber(1,i3)-ramber(1,i1)
            dnx=ramber(3,i4)-ramber(3,i3)
            dnz=ramber(1,i3)-ramber(1,i4)
         elseif(i.eq.8) then
            dmx=ramber(3,i2)-ramber(3,i1)
            dmz=ramber(1,i1)-ramber(1,i2)
            dnx=ramber(3,i2)-ramber(3,i4)
            dnz=ramber(1,i4)-ramber(1,i2)
         elseif(i.eq.11) then
            dnx=ramber(3,i3)-ramber(3,i2)
            dnz=ramber(1,i2)-ramber(1,i3)
         elseif(i.eq.3) then
            dmx=ramber(2,i2)-ramber(2,i3)
            dmy=ramber(1,i3)-ramber(1,i2)
         elseif(i.eq.6) then
            dmx=ramber(2,i3)-ramber(2,i1)
            dmy=ramber(1,i1)-ramber(1,i3)
            dnx=ramber(2,i3)-ramber(2,i4)
            dny=ramber(1,i4)-ramber(1,i3)
         elseif(i.eq.9) then
            dmx=ramber(2,i1)-ramber(2,i2)
            dmy=ramber(1,i2)-ramber(1,i1)
            dnx=ramber(2,i4)-ramber(2,i2)
            dny=ramber(1,i2)-ramber(1,i4)
         elseif(i.eq.12) then
            dnx=ramber(2,i2)-ramber(2,i3)
            dny=ramber(1,i3)-ramber(1,i2)
         endif
         call pbc_displ_vector(lattice_type,cell,kcell,dmx,
     .     dmy,dmz)
         call pbc_displ_vector(lattice_type,cell,kcell,dnx,
     .     dny,dnz)
         dm=(mx*dmx+my*dmy+mz*dmz)/rm
         dn=(nx*dnx+ny*dny+nz*dnz)/rn
         dmn=rm*dn+rn*dm
         dscalar=nx*dmx+mx*dnx+ny*dmy+my*dny+nz*dmz+mz*dnz
         fce(i)= dtot*(dscalar*rm*rn-dmn*scal)/(rn*rm)**2
      enddo
 10   end

c*******************************************************************
c     subroutine that calculates the variables of a dihedral angle

      subroutine dihevars(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,
     .     mx,my,mz,m,nx,ny,nz,n,dihedro,cell)
      use precision, only: dp
      implicit none
      real(dp) dihedro,dist,angle,l1,l2,pi,arg
      real(dp) x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4
      real(dp) l,A,B,C,D,m,n,nx,ny,nz,mx,my,mz,scalar
      real(dp) cell(3,3)
      real(dp) dx12, dy12, dz12, dx32, dy32, dz32
      real(dp) dx23, dy23, dz23, dx43, dy43, dz43

      pi=DACOS(-1.d0)
      dx12 = x1-x2
C     Count for PBC along X
      dx12 = dx12 - ANINT(dx12/cell(1,1))*cell(1,1)
      dy12 = y1-y2
C     Count for PBC along Y
      dy12 = dy12 - ANINT(dy12/cell(2,2))*cell(2,2)
      dz12 = z1-z2
C     Count for PBC along Z
      dz12 = dz12 - ANINT(dz12/cell(3,3))*cell(3,3)
      dx32 = x3-x2
C     Count for PBC along X
      dx32 = dx32 - ANINT(dx32/cell(1,1))*cell(1,1)
      dy32 = y3-y2
C     Count for PBC along Y
      dy32 = dy32 - ANINT(dy32/cell(2,2))*cell(2,2)
      dz32 = z3-z2
C     Count for PBC along Z
      dz32 = dz32 - ANINT(dz32/cell(3,3))*cell(3,3)
      mx = dy12*dz32 - dz12*dy32
      my = -(dx12*dz32 - dz12*dx32)
      mz = dx12*dy32-dy12*dx32

      dx23 = x2-x3
C     Count for PBC along X
      dx23 = dx23 - ANINT(dx23/cell(1,1))*cell(1,1)
      dy23 = y2-y3
C     Count for PBC along Y
      dy23 = dy23 - ANINT(dy23/cell(2,2))*cell(2,2)
      dz23 = z2-z3
C     Count for PBC along Z
      dz23 = dz23 - ANINT(dz23/cell(3,3))*cell(3,3)
      dx43 = x4-x3
C     Count for PBC along X
      dx43 = dx43 - ANINT(dx43/cell(1,1))*cell(1,1)
      dy43 = y4-y3
C     Count for PBC along Y
      dy43 = dy43 - ANINT(dy43/cell(2,2))*cell(2,2)
      dz43 = z4-z3
C     Count for PBC along Z
      dz43 = dz43 - ANINT(dz43/cell(3,3))*cell(3,3)
      nx = dy23*dz43 - dz23*dy43
      ny = -(dx23*dz43 - dz23*dx43)
      nz = dx23*dy43-dy23*dx43

c     calculo del prod scalar n*m
      scalar = mx*nx + my*ny + mz*nz
      m = mx**2 + my**2 + mz**2
      m = m**(0.5)
      n = nx**2 + ny**2 + nz**2
      n = n**(0.5)
      arg=scalar/(m*n)
      if(arg.ge.1.0) then 
         dihedro=0.0 
      elseif(arg.le.-1.0) then
         dihedro=180.0
      else
         dihedro = dACOS(arg)
      endif
      dihedro = dihedro*180/pi 
 10   continue 
      end

      subroutine dihevars_v2(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,
     .     mx,my,mz,m,nx,ny,nz,n,dihedro,cell,kcell,lattice_type)
      use precision, only: dp
      implicit none
      real(dp) dihedro,dist,angle,l1,l2,pi,arg
      real(dp) x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4
      real(dp) l,A,B,C,D,m,n,nx,ny,nz,mx,my,mz,scalar
      real(dp) dx12, dy12, dz12, dx32, dy32, dz32
      real(dp) dx43, dy43, dz43
      real(dp) cell(3,3), kcell(3,3)
      character lattice_type

      pi=DACOS(-1.d0)
      dx12=x1-x2
      dy12=y1-y2
      dz12=z1-z2
      call pbc_displ_vector(lattice_type,cell,kcell,dx12,
     .     dy12,dz12)
      dx32=x3-x2
      dy32=y3-y2
      dz32=z3-z2
      call pbc_displ_vector(lattice_type,cell,kcell,dx32,
     .     dy32,dz32)
      mx = dy12*dz32 - dz12*dy32
      my = -(dx12*dz32 - dz12*dx32)
      mz = dx12*dy32-dy12*dx32

      dx43=x4-x3
      dy43=y4-y3
      dz43=z4-z3
      call pbc_displ_vector(lattice_type,cell,kcell,dx43,
     .     dy43,dz43)
      nx = -dy32*dz43 + dz32*dy43
      ny = dx32*dz43 - dz32*dx43
      nz = -dx32*dy43+dy32*dx43

c     calculo del prod scalar n*m
      scalar = mx*nx + my*ny + mz*nz
      m = mx**2 + my**2 + mz**2
      m = m**(0.5)
      n = nx**2 + ny**2 + nz**2
      n = n**(0.5)
      arg=scalar/(m*n)
      if(arg.ge.1.0) then 
         dihedro=0.0 
      elseif(arg.le.-1.0) then
         dihedro=180.0
      else
         dihedro = dACOS(arg)
      endif
      dihedro = dihedro*180/pi 
 10   continue 
      end

c     *************************************************************************
c     subroutine that calculates position of 4 from 1-3 

      subroutine pos4(x1,y1,z1,x2,y2,z2,x3,y3,z3,r4,a4,d4,x4,y4,z4,
     .     cell)
      use precision, only: dp
      implicit none
      real(dp) x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4
      real(dp) r4,a4,d4,xejx,yejx,zejx,xejz,yejz,zejz
      real(dp) pi,rejx,r23,rejz,x,y,z
      real(dp) l1,m1,n1,l2,m2,n2,l3,m3,n3
      real(dp) cell(3,3)
      real(dp) dx32, dy32, dz32, dx12, dy12, dz12
      
      pi=DACOS(-1.d0)
      
c     directrices cosines l1,m1,n1,etc
c     l1: x' axis respect to x 
c     m1 x'axis respect to y
c     n1 x' axis respect to z
c     2 and 3 respect to y and z axis
c     x' axis ejx
      
      dx32 = x3-x2
C     Count for PBC along X
      dx32 = dx32 - ANINT(dx32/cell(1,1))*cell(1,1)
      dy32 = y3-y2
C     Count for PBC along Y
      dy32 = dy32 - ANINT(dy32/cell(2,2))*cell(2,2)
      dz32 = z3-z2
C     Count for PBC along Z
      dz32 = dz32 - ANINT(dz32/cell(3,3))*cell(3,3)
      dx12 = x1-x2
C     Count for PBC along X
      dx12 = dx12 - ANINT(dx12/cell(1,1))*cell(1,1)
      dy12 = y1-y2
C     Count for PBC along Y
      dy12 = dy12 - ANINT(dy12/cell(2,2))*cell(2,2)
      dz12 = z1-z2
C     Count for PBC along Z
      dz12 = dz12 - ANINT(dz12/cell(3,3))*cell(3,3)

      xejx = dy32*dz12 - dz32*dy12
      yejx = -(dx32*dz12 - dz32*dx12)
      zejx = dx32*dy12-dy32*dx12
      
      rejx = dsqrt(xejx**2+yejx**2+zejx**2)
      
      l1 = xejx/rejx
      
      r23=dsqrt(dx32**2+dy32**2+dz32**2)
      m1 = yejx/rejx
      n1= zejx/rejx
      
      l2 = dx32/r23
      m2 = dy32/r23
      n2 = dz32/r23
      
      
      xejz = yejx*dz32 - zejx*dy32
      yejz = -(xejx*dz32 - zejx*dx32)
      zejz = xejx*dy32- yejx*dx32
      
      rejz = dsqrt(xejz**2+yejz**2+zejz**2)
      
      l3 = xejz/rejz
      m3 = yejz/rejz
      n3 = zejz/rejz

c     at4 coords in x' y' z' system
c     dihedral equivalent to angle in z' x' plane
c     180-angle equivalent to angle with y' axis
c     r4=distance at3 - at4 
      d4 = d4*pi/180.0
      a4 =180.0-a4
      a4 = a4*pi/180.0
      
      z = r4*dSIN(a4)*dCOS(d4)
      x = r4*dSIN(a4)*dSIN(d4)
      y = r4*dCOS(a4)
      
      y = y + r23
      
c     translating 
      
      x4 = (l1*x + l2*y + l3*z + x2)
      y4 = (m1*x + m2*y + m3*z + y2)
      z4 = (n1*x + n2*y + n3*z + z2)
      
      return
      end subroutine pos4

c******************************************************************************
