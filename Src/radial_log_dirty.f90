!     This file is part of the SIESTA package.
!     
!     Copyright (c) Fundacion General Universidad Autonoma de Madrid:
!     E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal
!     and J.M.Soler, 1996-2006.
!     
!     Use of this software constitutes agreement with the full conditions
!     given in the SIESTA license, as signed by all legitimate users.
!     

! radial log functions with an interface which has not been updated.

module radial_log_dirty
  use precision, only: dp
  use radial_log
  use radial_logGrid, only: logGrid_t, log_grid_get_number_of_points, &
       log_grid_get_r
  use radial_log_low_level, only: log_rad_func_t
  use radial_util, only: fit_parabola,findp
  implicit none
  public

  contains

     !---------------------------------------------------------------------------

  subroutine log_rad_smooth(func,zval,r_match, smooth, integral)
    !    Corresponds to the old vloc1 subroutine
    !    Local-potential size parameter 'rgauss'
    !    We choose as a smooth pseudopotential the one generated 
    !    by a 'Vanderbilt-function' charge distribution. We have to select 
    !    the size of this distribution somehow.
    !    'Vanderbilt-functions' are of the form :
    !    p(r)=N*exp(-(sinh(van*r)/sinh(van)**2)
    !    when van---> 0 we will obtain a 'gaussian'
    !    when van---> Inf. we will obtain a step function
    !    Some test has revealed that the best election to achieve 
    !    a good convergence in real and reciprocal space is b in the 
    !    range 0.5-1.0 .
    !    *
    
    !    So, the 'gaussian' charge distribution 
    !    must go to zero at a distance 'rgauss'.
    
    type(log_rad_func_t), intent(in) :: func
    real(dp), intent(in) :: zval,r_match

    type(log_rad_func_t), intent(out)::smooth
    type(log_rad_func_t), intent(out):: integral

    !double precision, intent(in) :: Zval, a
    !integer, intent(in) :: nrval
    !real(dp), intent(in) :: rofi(:), drdi(:), s(:)
    !real(dp), intent(out) :: vlocal(:)
    !real(dp), intent(out) :: chlocal(:)
    !integer, intent(out) :: nchloc
    !real(dp), intent(inout) :: rgauss      !!???

    !     *Internal variables* 

    real(dp) van, factor, alp, cutoff1, cutoff2, &
         gexp, qtot, pi, chc, r, Rchloc, rhor1, rhor, rpb, ea
    real(dp), dimension(:), pointer :: chlocal,vlocal,s
    integer ir,nrval,nchloc
    type(logGrid_t), pointer :: grid
    real(dp) :: eps
    parameter(eps=1.0d-4)
    !     AG       Maximum exponent for exp(-x)
    real(dp), parameter :: exp_range=60.d0
    !     AG

    !**   Usual local potential 
    !     (generated with an optimum Vanderbilt function)**


    !***  The very first local potential used by SIESTA was 
    !     the electrostatic potential generated by a gaussian 
    !     distribution ===> loctype='old' 
    !     loctype='old'
    !***  
    pi=acos(-1.0d0)

    !    Local-potential size parameter 'rgauss'
    !    We choose as a smooth pseudopotential the one generated 
    !    by a 'Vanderbilt-function' charge distribution. We have to select 
    !    the size of this distribution somehow.
    !    'Vanderbilt-functions' are of the form :
    !    p(r)=N*exp(-(sinh(van*r)/sinh(van)**2)
    !    when van---> 0 we will obtain a 'gaussian'
    !    when van---> Inf. we will obtain a step function
    !    Some test has revealed that the best election to achieve 
    !    a good convergence in real and reciprocal space is b in the 
    !    range 0.5-1.0 .
    !    *

    !    So, the 'gaussian' charge distribution 
    !    must go to zero at a distance 'rgauss'.



    !    We take a 'Vanderbilt-function' as local potential
    !    van=1.0d0 all the parameter have optimized for this value 

    grid => func%grid
    van=1.0d0
    cutoff1=3.63d0
    cutoff2=5.48d0
    !**   99% of charge inside Rgauss**
    !     factor=1.627d0

    !**   99.9% of charge inside Rgauss
    factor=1.815d0

    !    * Scaling factor for local-pseudopot. charge**
    alp=factor/r_match

    write(6,'(/,a,f10.3,a)') 'VLOCAL1: 99.0% of the norm of Vloc inside ',&
         (alp*cutoff1)**2,' Ry'
    write(6,'(a,f10.3,a)') 'VLOCAL1: 99.9% of the norm of Vloc inside ',&
         (alp*cutoff2)**2,' Ry'     
    !the above message shouldn't be inside radial_log"

    !--------------------

    nrval = log_grid_get_number_of_points(grid)
    allocate(chlocal(1:nrval),vlocal(1:nrval),s(1:nrval))
   
    qtot=0.0d0 
    do ir=1,nrval
       r=grid%r(ir) 
       gexp=sinh(van*alp*r)/sinh(van)
       gexp=gexp*gexp
       !    AG  ---- prevent underflow
       rhor = 0.d0
       if (gexp .lt. exp_range) rhor=exp(-gexp)  
       !    AG
       !     if(ir.eq.1) rhor1=-rhor
       chlocal(ir)=(-4.0d0)*pi*rhor*r*r
       qtot=qtot+rhor*grid%drdi(ir)*r*r
    enddo

    gexp=sinh(van*alp*grid%r(1))/sinh(van)
    gexp=gexp*gexp
    rhor1 = 0.0_dp
    if (gexp .lt. exp_range) rhor1=exp(-gexp)

    qtot=4.0_dp*pi*qtot 
    nchloc=0  
    do ir=nrval,1,-1
       chc=zval*chlocal(ir)/qtot
       chlocal(ir)=chc
       !print *, ir,chlocal(ir)
       if((abs(chc).gt.eps).and.(nchloc.eq.0)) then    
          nchloc=ir+1
       endif
    enddo
   
    Rchloc=grid%r(nchloc)


    rpb=grid%b
    ea=exp(grid%a)
    do ir=1,nrval
       s(ir)=sqrt(grid%a*rpb)
       rpb=rpb*ea       
    enddo

    call vhrtre(chlocal,vlocal,grid%r,grid%drdi,s,nrval,grid%a)
    !print *, "v=", vlocal(1),vlocal(20),vlocal(100),vlocal(300)

    do ir=2,nrval
       r=grid%r(ir)  
       chlocal(ir)=chlocal(ir)/(4.0d0*pi*r*r)
      
       if (r.gt.1.1d0*Rchloc) then
          vlocal(ir)=(-2.0_dp)*zval/grid%r(ir)
       endif
    enddo

    chlocal(1)=-1.0_dp*rhor1*zval/qtot

    call log_rad_alloc(smooth,chlocal(1:nchloc),grid)
    call log_rad_alloc(integral,vlocal,grid)

    deallocate(chlocal,vlocal,s)
  end subroutine log_rad_smooth

  !----------------------------------------------------------------

    subroutine log_rad_parab_params(func,l, splnorm,rmatch, cons1,cons2,dpbug) 

!!$C
!!$C For a value of the SplitNorm parameter equal
!!$C to splnorm, this routine returns 
!!$C the parameters Cons1, Cons2 and Nm, such that 
!!$C the doble-Z function would be defined as 
!!$C
!!$C     phi(r)=r^l (cons1*r^2+cons2) if r < Rm
!!$C
!!$C     phi(r)= rphi(ir)/r           if r > Rm
!!$C    
!!$C with  Rm= b*[exp(a(nm-1)) -1 ] 
!!$C Continuity in the wavefunction and its derivative
!!$C is imposed.   
!!$C The arrays rphi and rnrm belong to the input
!!$C rphi(nrmax): The original PAO function multiply
!!$C   by the radius.
!!$C rnrm(nrmax): PAOs norm as a function of the 
!!$C   radius. 
!!$C
!!$C  Written by D. Sanchez-Portal, July 1997.
!!$C  
!!$C Algorithm based on routine for Golden Section Search
!!$C from Numerical Recipes.
!!$C

     type(log_rad_func_t), intent(in) :: func
    !real(dp), intent(in)    ::  a, b
    integer, intent(in)     ::  l
    real(dp), intent(in)    ::  splnorm
    real(dp), intent(out)   :: rmatch
    real(dp), intent(out)   ::  cons1
    real(dp), intent(out)   :: cons2
    !integer, intent(out)    ::  nm
    logical, intent(in) :: dpbug

    real(dp),  parameter  :: Ratio=0.61803399D0       

    real(dp):: slopold, slop, rmin, gmin, cmin, rnrmin
    real(dp):: gmax, cmax, rmax, rnrmax, valmin, valmax, gmed
    real(dp):: cmed, rmed, rnrmed, valmed, g1, c1, r, rn1, val1
    real(dp):: g2, c2, rn2, val2, dnrm
    integer :: n0, n1, n2, n3, nm
    integer :: ir, nr_max, nmin, nmax, nmed, iter, nlast, nrc

    !real(dp), pointer :: rphi(:)
    real(dp), pointer :: norm(:)
    real(dp) :: a,b
    
    a=func%grid%a
    b=func%grid%b
    !rphi=>func%f

    nrc = size(func%f)
    !print *, "nrc=",nrc
    allocate(norm(1:nrc))

    norm = 0.0_dp
    dnrm = 0.0_dp
    do ir=1,nrc
       dnrm = dnrm + func%Grid%drdi(ir)*func%f(ir)*func%f(ir)
       norm(ir)= dnrm
    enddo


    !         Find the last maximum of the wavefunction

    nlast=nrc-2
    slopold=0.0d0
    do ir=nlast,2,-1
       slop=func%f(ir)-func%f(ir-1)
       if(slop*slopold.lt.0.0d0) exit
       slopold=slop
    enddo

    nr_max=ir-1
    rmin=b*(exp(a*(nr_max-1))-1.0d0)
    rmin=1.01d0*rmin
    nmin=nint(dlog(rmin/b+1.0d0)/a)+1
    nmin=max(nmin,2)
    nmax=nrc-1

    !
    !         Initial brackets
    !
    call findp(nrc,nmin,func%f,a,b,l,cmin,gmin,dpbug)

    rmin=b*(exp(a*(nmin-1))-1.0d0)
    call nrmpal(cmin,gmin,rmin,l,rnrmin)
    rnrmin=1.0d0+rnrmin-norm(nmin)
    call findp(nrc,nmax,func%f,a,b,l,cmax,gmax,dpbug)
    
    rmax=b*(exp(a*(nmax-1))-1.0d0)
    call nrmpal(cmax,gmax,rmax,l,rnrmax)
    rnrmax=1.0d0+rnrmax-norm(nmax)
    !
    !         Start the algorithm to find the matching point at
    !         which the *total* norm of the parabola+tail = splitnorm
    !         (compare with the JPC paper: there it appears that only
    !         the tail should have a norm=splitnorm.
    !

    ! Under certain circunstances the algorithm is not going to work
    if(rnrmin.gt.splnorm.and.rnrmax.gt.splnorm) then
       write(6,'(/,A,/,A)')'parabola: The program failed in finding a SPLIT orbital ',&
            'parabola: with the desired splitnorm'
       call die()
    endif


    valmin=(splnorm-rnrmin)**2
    valmax=(splnorm-rnrmax)**2
    nmed=(nmin+nmax)/2

    do iter=1,nrc
      
       call findp(nrc,nmed,func%f,a,b,l,cmed,gmed,dpbug)

       rmed=b*(exp(a*(nmed-1))-1.0d0)
       call nrmpal(cmed,gmed,rmed,l,rnrmed)
       rnrmed=1.0d0+rnrmed-norm(nmed)

       valmed=(splnorm-rnrmed)**2

       if((valmed.lt.valmin).and.(valmed.lt.valmax)) goto 20
       nmed=nmed+1
       if(nmed.eq.nmax) goto 15
    enddo
15  continue
    nmed=(nmin+nmax)/2
    do iter=1,nrc
       call findp(nrc,nmed,func%f,a,b,l,cmed,gmed,dpbug)
       
       rmed=b*(exp(a*(nmed-1))-1.0d0)
       call nrmpal(cmed,gmed,rmed,l,rnrmed)
       rnrmed=1.0d0+rnrmed-norm(nmed)

       valmed=(splnorm-rnrmed)**2


       if((valmed.lt.valmin).and.(valmed.lt.valmax)) goto 20
       nmed=nmed-1
       if(nmed.eq.nmin) goto  20
    enddo
20  continue

    if(nmed.eq.nmin) then
       if(valmin.lt.valmax) then
          nm=nmin
          cons1=cmin
          cons2=gmin
       elseif(valmax.le.valmin) then
          nm=nmax
          cons1=cmax
          cons2=gmax
       endif
       return
    endif

    ! The interval has been bracketed.   
    n0=nmin
    n3=nmax
    if(abs(nmed-nmax).gt.abs(nmed-nmin)) then
       n1=nmed
       n2=nmed+nint((1.0d0-ratio)*(nmax-nmed))
    else
       n2=nmed
       n1=nmed-nint((1.0d0-ratio)*(nmed-nmin))
    endif

    call findp(nrc,n1,func%f,a,b,l,c1,g1,dpbug)


    r=b*(exp(a*(n1-1))-1.0d0)
    call nrmpal(c1,g1,r,l,rn1)
    rn1=1.0d0+rn1-norm(n1)
    val1=(splnorm-rn1)**2

    call findp(nrc,n2,func%f,a,b,l,c2,g2,dpbug)
    
    r=b*(exp(a*(n2-1))-1.0d0)
    call nrmpal(c2,g2,r,l,rn2)
    rn2=1.0d0+rn2-norm(n2)
    val2=(splnorm-rn2)**2

1   if(abs(n3-n0).gt.1) then
       if(val2.lt.val1) then
          n0=n1
          n1=n2
          n2=nint(ratio*n1+(1-ratio)*n3)
          !              val0=val1
          val1=val2
          call findp(nrc,n2,func%f,a,b,l,c2,g2,dpbug)
          r=b*(exp(a*(n2-1))-1.0d0)
          call nrmpal(c2,g2,r,l,rn2)
          rn2=1.0d0+rn2-norm(n2)
          val2=(splnorm-rn2)**2
       else
          n3=n2
          n2=n1
          n1=nint(ratio*n2+(1-ratio)*n0)
          !              val3=val2
          val2=val1
          call findp(nrc,n1,func%f,a,b,l,c1,g1,dpbug)
          r=b*(exp(a*(n1-1))-1.0d0)
          call nrmpal(c1,g1,r,l,rn1)
          rn1=1.0d0+rn1-norm(n1)
          val1=(splnorm-rn1)**2
       endif
       goto 1
    endif
    if(val1.lt.val2) then
       cons2=g1
       cons1=c1
       nm=n1
    else
       cons2=g2
       cons1=c2
       nm=n2
    endif

    if(restricted_grid) nm=nm+1-mod(nm,2)
    rmatch = log_rad_get_r(func,nm)
    deallocate(norm)
  end subroutine log_rad_parab_params

  !--------------------------------------------------------

   subroutine log_rad_smooth_large(f,zval,r_match,smooth,integral)
    type(log_rad_func_t), intent(in) :: f
    real(dp), intent(in) :: zval,r_match

    type(log_rad_func_t), intent(out)::smooth
    type(log_rad_func_t), intent(out):: integral

  !subroutine vlocal2(Zval, nrval, a, rofi, drdi, s, f,&
  !     nrgauss,vlocal,nchloc,chlocal) 

    !    This routine generates the local pseudopotential appropiate 
    !    for species with  a large core.
    !    Written by D. Sanchez-Portal, Aug. 1998

    !real(dp), intent(in) :: Zval, a
    !integer, intent(in) :: nrval
    !integer, intent(inout) :: nrgauss
    !real(dp), intent(in) :: rofi(:), drdi(:), s(:),f(:)
    !real(dp), intent(out) :: vlocal(:), chlocal(:)
    !integer, intent(out) :: nchloc


    !    Internal variables****

    real(dp)  vlc, r, dev, dev2, dev3, var1, var2, var3, v1, v2, v3, v4,&
         dm11, dm12, dm13, dm21, dm22, dm23, dm31, dm32, dm33, &
         g0, g1, g2, g3, g4, d2g, d2u, cons, a2b4, qtot, pi, a,b, rpb, ea   
    integer::  ndevfit, ir, nrgauss, nrval, nchloc

    real(dp), parameter  :: eps=1.0d-5

    real(dp), dimension(:), pointer :: vlocal, chlocal,s

    type(logGrid_t), pointer :: grid

    pi=acos(-1.0d0)        

    !    Continuity up to second derivative***
    ndevfit=2
    !    Continuity up to third derivative****
    !    ndevfit=3

    grid => f%grid
    a = grid%a
    b = grid%b

    nrval = log_rad_get_length(f)
    nrgauss = log_rad_get_ir(f,r_match) 

    nrgauss=nrgauss+3        !! For good measure...

    allocate(vlocal(1:nrval),chlocal(1:nrval),s(1:nrval))
    do ir=1,nrval
       vlocal(ir)=f%f(ir)*grid%r(ir)
    enddo
    
    rpb=b
    ea=exp(a)
    do ir=1,nrval
       s(ir)=sqrt(a*rpb)
       rpb=rpb*ea       
    enddo

    ir=nrgauss
    dev=(vlocal(ir+1)-vlocal(ir-1))*0.5d0
    dev2=(vlocal(ir+1)+vlocal(ir-1)-2.0d0*vlocal(ir))
    dev3=(vlocal(ir+2)-2.0d0*vlocal(ir+1) &
         +2.0d0*vlocal(ir-1)-vlocal(ir-2))*0.5d0
    dev3=(dev3-3.0d0*a*dev2+2.0d0*(a**2)*dev)&
         /(grid%drdi(ir)**3)
    dev2=(dev2-a*dev)/(grid%drdi(ir)**2)
    dev=dev/grid%drdi(ir)

    !    Local potential is Vloc(r)=v3*exp(v1*r^2+v2*r^3) 
    !    inside Rgauss and equals the 
    !    all-electron atomic potential outside Rgauss
    !    We impose the continuity up to second        

    if(ndevfit.eq.2) then               
       vlc=vlocal(nrgauss)
       r=grid%r(nrgauss)

       var1=dev/vlc-1.0d0/r
       var2=dev2/vlc-2.0d0*var1/r -(var1**2)

       dm11=2.0d0*r
       dm12=3.0d0*r*r
       dm21=2.0d0
       dm22=6.0d0*r

       v1=(dm22*var1-dm12*var2)/(6.0d0*r*r)
       v2=(dm11*var2-dm21*var1)/(6.0d0*r*r)
       v3=vlc/(r*exp((v1+v2*r)*r*r))


       !     elseif(ndevfit.eq.3) then 
    else

       !    We can also construct a local potential 
       !    Vloc(r)=v4*exp(v1*r^2+v2*r^3+v3*r^4),
       !    this new coefficient allows us to impose the continuity 
       !    of the potential up  to the third derivative.

       vlc=vlocal(nrgauss)
       r=grid%r(nrgauss)

       var1=dev/vlc-1.0d0/r
       var2=dev2/vlc-2.0d0*var1/r-(var1**2)
       var3=dev3/vlc-3.0d0*var1*var2-(var1**3) &
            -3.0d0*(var1**2+var2)/r

       dm11=2.0d0*r
       dm12=3.0d0*r*r
       dm13=4.0d0*r*r*r
       dm21=2.0d0
       dm22=6.0d0*r
       dm23=12.0d0*r*r
       dm31=0.0d0
       dm32=6.0d0
       dm33=24.0d0*r

       v1=((var1*dm22*dm33+var2*dm13*dm32+var3*dm12*dm23) &
            -(var3*dm22*dm13+var1*dm32*dm23+var2*dm12*dm33))/(48.0d0*r*r*r) 
       v2=((var2*dm11*dm33+var3*dm21*dm13+var1*dm23*dm31) &
            -(var2*dm31*dm13+var3*dm23*dm11+var1*dm21*dm33))/(48.0d0*r*r*r)
       v3=((var3*dm11*dm22+var2*dm12*dm31+var1*dm32*dm21) &
            -(var1*dm22*dm31+var3*dm21*dm12+var2*dm11*dm32))/(48.0d0*r*r*r)
       v4=vlc/(r*exp((v1+v2*r+v3*r*r)*r*r))

    endif

    do ir=1,nrval
       r=grid%r(ir)
       if(ir.le.nrgauss) then 

          !**   If second derivative fit***
          if(ndevfit.eq.2) then 
             vlocal(ir)=v3*exp((v1+v2*r)*r*r)

             !**   If third derivative fit****
          elseif(ndevfit.eq.3) then 

             vlocal(ir)=v4*exp((v1+v2*r+v3*r*r)*r*r)
             !**** 
          endif

       else
          vlocal(ir)=f%f(ir)
       endif

    enddo

    !    Once we have the local potential we define the 'local-pseudopotential 
    !    charge' which help us to calculate the electrostatic interation 
    !    between the ions

    a2b4=0.25d0*a*a 
    qtot=0.d0 

    do ir=1,nrval-1

       g2=vlocal(ir)*grid%r(ir)
       if(abs(g2+2.0d0*zval).lt.eps) goto 10

       if(ir.gt.nrgauss) then  

          if((ir.gt.2).and.(ir.lt.(nrval-1))) then 
             g0=vlocal(ir-2)*grid%r(ir-2)/s(ir-2)
             g1=vlocal(ir-1)*grid%r(ir-1)/s(ir-1)
             g2=vlocal(ir)*grid%r(ir)/s(ir)
             g3=vlocal(ir+1)*grid%r(ir+1)/s(ir+1)
             g4=vlocal(ir+2)*grid%r(ir+2)/s(ir+2)

             d2g=(16.0d0*(g1+g3)-(g0+g4)-30.0d0*g2)/12.0d0

          else
             g1=vlocal(ir-1)*grid%r(ir-1)/s(ir-1)
             g2=vlocal(ir)*grid%r(ir)/s(ir) 
             g3=vlocal(ir+1)*grid%r(ir+1)/s(ir+1)  

             d2g=g1+g3-2.0d0*g2

          endif

          d2u=d2g-a2b4*g2

          r=grid%r(ir)
          cons=8.0d0*pi*r*grid%drdi(ir)*s(ir)
          chlocal(ir)=(-d2u)/cons
          qtot= qtot + 0.5d0*d2u*r/s(ir)

       else

          !    If second derivative fit*** 
          if(ndevfit.eq.2)  then
             r=grid%r(ir)

             g0=v3*exp((v1+v2*r)*r**2)
             g1=(2.0d0*v1+3.0d0*v2*r)
             g2=2.0d0*v1+6.0d0*v2*r
             g3=(g2+g1*g1*r*r+2.0d0*g1)*g0

             cons=8.0d0*pi
             chlocal(ir)= (-g3)/cons
             qtot= qtot  + 0.5d0*g3*r*r*grid%drdi(ir)

             !**** If third derivative fit

          elseif(ndevfit.eq.3)  then

             r=grid%r(ir)

             g0=v4*exp((v1+v2*r+v3*r*r)*r*r)     
             g1=(2.0d0*v1+3.0d0*v2*r+4.0d0*v3*r*r)
             g2=(2.0d0*v1+6.0d0*v2*r+12.0d0*v3*r*r)    
             g3=(g2+g1*g1*r*r+2.0d0*g1)*g0   

             cons=8.0d0*pi
             chlocal(ir)= -g3/cons
             qtot= qtot  + 0.5d0*g3*r*r*grid%drdi(ir)
          endif

       endif
    enddo

10  continue
    nchloc=ir          
    
    do ir=1,nchloc-1
       chlocal(ir)=zval*chlocal(ir)/qtot
    enddo
    do ir=nchloc,nrval
       chlocal(ir)=0.0d0
    enddo

    call log_rad_alloc(smooth,chlocal(1:nchloc),grid)
    call log_rad_alloc(integral,vlocal,grid)

    deallocate(vlocal,chlocal,s)

  endsubroutine log_rad_smooth_large

  !----------------------------------------------------------------

    subroutine log_rad_ghost(func,vps_f,vlocal_f,ve_f,l,eigenl,zval,ighost)
!      subroutine ghost(Zval,rofi,vps,vlocal,ve,s,drdi,nrval,l,a,b, &
!           nrc, eigenl,rphi,ighost)

!    This routine checks the possible existence of ghost states. 
!    Input:
!    vps(*)    : pseudopotential for angular momentum l
!    vlocal(*) : local pseudopotential
!    rphi (*)  : first radial pseudowavefunctions for Vps.
!    eigenl    : eigenvalue 
!    
!    Output:
!    ighost:   if ighost=0 no ghost states
!    if ighost=1, ghost states exist
!    
!    Written by D. Sanchez-Portal, Aug. 1998
!    

    type (log_rad_func_t), intent(in) :: func ! The function to be tested.
    type (log_rad_func_t), intent(in) :: vps_f, vlocal_f, ve_f
    integer, intent(in) :: l
    real(dp) :: zval, eigenl
    integer, intent(out) :: ighost

    !    * Internal variables***
    !integer, intent(inout)   :: nrc
    
    real(dp), pointer :: a, b
    real(dp), pointer :: rofi(:), vps(:), vlocal(:), ve(:)
    real(dp), pointer :: s(:), drdi(:), g(:)
    real(dp) ::   dnrm, avgv, phi, e, elocal(2), vl, vphi, dkbcos
    
    integer ::  ir, nprin, nnode, nrval, nrc
    logical ::  called
    !**   SAVE FOR NEXT CALL
    !    
    save called
    data  called /.false./
    !    
    if(.not.called) then 
       ighost=0 
       called=.true.
    endif

    a=> func%grid%a
    b=> func%grid%b
    rofi => func%grid%r
    vps => vps_f%f
    ve => ve_f%f
    vlocal => vlocal_f%f
    drdi => func%grid%drdi
    
    nrc = log_rad_get_ir(func,log_rad_cutoff(func))
    nrval = size(func%f)
    if(nrval .gt. maximum_length) nrval=maximum_length

    allocate(s(1:nrval), g(1:nrval))
    s(2:nrval) = drdi(2:nrval)*drdi(2:nrval)
    s(1) = s(2)
    !***  CALCULATE EIGENVALUES OF LOCAL POTENTIAL FOR GHOST ANALYSIS*
    ! ATTENTION , 'Ve' is the screenig potential generated from valence*
    ! pseudo-charge given by the pseudopotential generation program ****

    do nnode=1,2
       nprin=l+1

       call schro_eq(Zval,rofi,vlocal,ve,s,drdi,&
            nrval,l,a,b,nnode,nprin,e,g)

       elocal(nnode)=e

    enddo

    !    write(6,*) 'GHOST: Ground state vlocal for L=',l,elocal(1)
    !    write(6,*) 'GHOST: First excited state for L=',l,elocal(2)
    
    !    *CALCULATE KB-COSINE****

    nrc=min(nrc,nrval)
    dnrm=0.0d0
    Avgv=0.0d0

    do ir=2,nrc
       vl=(vps(ir)-vlocal(ir))
       phi=log_rad_get_value_from_ir(func,ir)
       vphi=vl*phi
       dnrm=dnrm+vphi*vphi*drdi(ir)
       avgv=avgv+vphi*phi*drdi(ir)
    enddo


    dkbcos=dnrm/(avgv+1.0d-20)
    !    dknrm=1.0d0/(sqrt(dnrm)+1.0d-20)

    !    GHOST ANALYSIS*


    if(dkbcos.gt.0.0d0) then

       if(eigenl.gt.elocal(2)) then
          write(6,"(a,i3)") 'GHOST: WARNING: Ghost state for L =', l
          ighost=1
       else
          write(6,'(a,i3)') '     GHOST: No ghost state for L =',l
          ighost=0
       endif

    elseif(dkbcos.lt.0d0) then
       
       if(eigenl.gt.elocal(1)) then
          write(6,"(a,i3)") 'GHOST: WARNING: Ghost state for L =', l
          ighost=1
       else
          write(6,'(a,i3)') '     GHOST: No ghost state for L =',l
          ighost=0
       endif
             
    elseif(dkbcos.eq.0.0d0) then
             
       write(6,"('     GHOST: vps = vlocal, no ghost for L =',i3)") l
       ighost=0 
       
    endif
    deallocate(s,g)
    nullify(a,b,rofi,vps,ve,vlocal,drdi)
    
  end subroutine log_rad_ghost

!-------------------------------------------------------------------------

  function log_rad_kb_proj(func,pseudo,pseudolocal,l,dkbcos,ekb) result(kbproj)

    !    This routine calculates the Kleinman-Bylander projector
    !    with angular momentum l.
    !    Written by D. Sanchez-Portal, Aug. 1998
    !    Modified by DSP to allow more than one projector per l, July 1999.

    !real(dp), intent(in) :: rofi(:), drdi(:), vps(:)
    !real(dp), intent(in) :: vlocal(:), rphi(:)
    !integer, intent(in) :: nrwf, l
    !integer, intent(out) :: nrc

    type(log_rad_func_t), intent(in) :: func
    type(log_rad_func_t), intent(in) :: pseudo
    type(log_rad_func_t), intent(in) :: pseudolocal
    integer, intent(in) :: l
    real(dp), intent(out) :: dkbcos
    real(dp), intent(out) :: ekb

    type(log_rad_func_t) :: kbproj
    !    * Internal variables***

    
    integer, parameter  :: nkbmx  =    2
    integer, parameter  :: nrmax  = 6000

    real(dp)  dnrm, vl, vphi, avgv, r, phi, dknrm, &
         dincv, rc, rphi2(nrmax,nkbmx), vii(nkbmx), &
         sum, vij(nkbmx)

    real(dp), pointer :: rofi(:),drdi(:), vps(:), vlocal(:),rphi(:),proj(:)
    integer ir, l_last, nkb_last, jkb, ikb,nrwf,nrc

    logical  called
    
    !parameter (eps=1.0d-6)

    save called, l_last, nkb_last, rphi2, vii
    data called /.false./

    rphi=>func%f
    rofi=>func%grid%r
    drdi=>func%grid%drdi
    vps=>pseudo%f
    vlocal=>pseudolocal%f

    !nrwf=log_rad_get_ir_from_r(func,log_rad_cutoff(func))
    !if(restricted_grid) nrwf=nrwf+1-mod(nrwf,2)
    nrwf = size(func%f)

    if(called) then 
       if(l.ne.l_last) then 
          l_last=l 
          nkb_last=0
       endif
    else
       called=.true.
       l_last=l
       nkb_last=0
    endif


    !    We need to orthogonalize to the previous projectors. 
    !    We follow the scheme proposed by Blochl, PRB 41, 5414 (1990)
    ikb=nkb_last+1
    if(nkb_last.eq.0) then 
       do ir=1,nrwf
          rphi2(ir,1)=rphi(ir)
       enddo
    else
       do jkb=1,nkb_last
          sum=0.0d0
          do ir=1,nrwf
             vl=vps(ir)-vlocal(ir)
             sum=sum+rphi(ir)*vl*rphi2(ir,jkb)*drdi(ir)
          enddo
          vij(jkb)=sum
       enddo
       do ir=1,nrwf
          sum=0.0d0
          do jkb=1,nkb_last
             sum=sum+vij(jkb)*rphi2(ir,jkb)/(vii(jkb)+1.0d-20)
          enddo
          rphi2(ir,ikb)=rphi(ir)-sum
       enddo
    endif
    nkb_last=ikb
    !    Normalize the new function            
    dnrm=0.0d0
    do ir=1,nrwf
       dnrm=dnrm+ drdi(ir)*(rphi2(ir,ikb)**2)
    enddo
    dnrm=sqrt(dnrm)
    if(dnrm.lt.eps) then 
       do ir=1,nrwf
          rphi2(ir,ikb)=0.0d0
       enddo
    else
       do ir=1,nrwf
          rphi2(ir,ikb)=rphi2(ir,ikb)/dnrm
       enddo
    endif

    dnrm=0.0d0
    avgv=0.0d0
    do ir=2,nrwf
       r=rofi(ir)
       vl=(vps(ir)-vlocal(ir))
       phi=rphi2(ir,ikb)
       vphi=vl*phi
       dnrm=dnrm+vphi*vphi*drdi(ir)
       avgv=avgv+vphi*phi*drdi(ir)
      
    enddo
    
    vii(ikb)=avgv

    ekb=dnrm/(avgv+1.0d-20)
    dknrm=1.0d0/(sqrt(dnrm)+1.0d-20)
    dkbcos=avgv*dknrm

    !print *, "ekb,dkbcos=", ekb, dkbcos

    !    DEFINE THE CUT-OFF RADII FOR THE KB PROJECTORS**
    !    Warning these radii should be quite short, if it is not the case
    !    something is probably wrong in this part of the program.
    !    It will display a warning if Rc>4.5 a.u.or Rc < 0.5a.u.!!!!!!!!!!!!

    do 20 ir=nrwf,2,-1 
       phi=(rphi2(ir,ikb)/rofi(ir))*dknrm
       dincv=abs((vps(ir)-vlocal(ir))*phi)
       if (dincv.gt.eps) then 
          if (ir.ge.nrwf-1) then
             write(6,"(2a,/,2a)") 'KBproj: WARNING: ',&
                  'KB projector does not decay to zero',&
                  'KBproj: WARNING: ', 'parameter Rmax in routine KBgen should be increased'
          endif
          exit   
       endif
20  enddo
 
    nrc=ir+1
    rc=rofi(nrc)
    
    if(rc.lt.0.5d0) then
       write(6,"('KBproj: WARNING: Rc(',i2,')=',f12.4)")l,rc
       write(6,"(2a)") 'KBproj: WARNING: ',&
            'Cut of radius for the KB projector too small' 
    elseif(rc.gt.4.5d0) then
       write(6,"('KBproj: WARNING: Rc(',i2,')=',f12.4)")l,rc
       write(6,"(2a)") 'KBproj: WARNING: ',&
            'Cut of radius for the KB projector too big'
       write(6,"(2a)") 'KBproj: WARNING: ',&
            'Increasing the tolerance parameter eps'
       write(6,"(a)") 'KBproj: WARNING: might be a good idea'
    endif

    allocate(proj(1:nrc))

    do 30 ir=2,nrc
       r=rofi(ir)
       vl=(vps(ir)-vlocal(ir))
       phi=rphi2(ir,ikb)/r
       vphi=vl*phi*dknrm
       proj(ir)=vphi/r**l
       !write(38,'(5f10.6)'),r,vl,phi,vphi,proj(ir)
30  end do

    proj(1)= ( proj(2)*rofi(3)**2 - proj(3)*rofi(2)**2 ) / &
         ( rofi(3)**2 - rofi(2)**2 )

    call log_rad_alloc(kbproj,proj,func%grid)
    deallocate(proj)
    nullify(rphi,rofi,drdi,vps,vlocal)
  endfunction log_rad_kb_proj

  !----------------------------------------------------------

  function log_rad_parabolic_split(func,rmatch,c1,c2,l,nsm,shells,lambdas) result(split)
    type(log_rad_func_t),intent(in) :: func
    real(dp), intent(in) :: c1,c2,rmatch
    integer, intent(in) :: l, nsm
    type(log_rad_func_t), dimension(1:nsm), intent(in) :: shells
    real(dp), dimension(:), intent(in) :: lambdas

    type(log_rad_func_t) :: split
    real(dp), pointer :: g(:), over(:), norm(:)
    real(dp)::rc1,dlt,d,dn,dnrm,r,const1,const2, phi
    integer :: nrc,ir,ism, nrc1, nrc2,nrc3,nrc4
    type(logGrid_t), pointer :: grid

    parameter(dlt=0.6_dp)

    grid => func%grid
    nrc = log_rad_get_ir(func,rmatch)
    allocate(g(1:nrc))
    if (nsm > 1) allocate(over(1:nsm),norm(1:nsm))

    do ir=2,nrc-1
       r=Grid%r(ir)         
       g(ir)=-(c1*r**2+c2)*r**(l+1)+func%f(ir) 
    enddo
    
    g(1)  =g(2)
    g(nrc)=0.0_dp
    

    if(nsm.gt.1) then
       do ism=1,nsm-1
          rc1=log_rad_cutoff(shells(ism))/lambdas(ism)

          nrc4=nint(dlog(rc1/Grid%b+1.0_dp)/Grid%a)+1
          rc1=rc1-dlt
          nrc3=nint(dlog(rc1/Grid%b+1.0_dp)/Grid%a)+1
          rc1=Grid%r(nrc3)
          const1=log_rad_get_value_from_ir(shells(ism),nrc3)
          const1=const1*log_grid_get_r(Grid,nrc3)
          d=dsqrt(const1**2+dlt**2)
          rc1=rc1+d+dlt
          nrc2=nint(dlog(rc1/Grid%b+1.0d0)/Grid%a)+1
          d=Grid%r(nrc2)-Grid%r(nrc3)
          rc1=Grid%r(nrc2)
          const2=0.5d0*(const1**2+(d)**2)/const1
          dnrm=0.0d0
          dn=0.0d0

          do ir=1, nrc2
             r=Grid%r(ir)
             if(ir.gt.nrc3.and.nrc.gt.nrc4) then
                phi=const2-const1*dsqrt(const2**2-(r-rc1)**2)/ &
                     dabs(const1)
             elseif(ir .ge. size(shells(ism)%f))then
                phi=0.0_dp
             else
                phi=log_rad_get_value_from_ir(shells(ism),ir)
                phi=phi*log_grid_get_r(grid,ir) 
             endif
             if (ir <= nrc) then
                dnrm=dnrm+Grid%drdi(ir)*phi*g(ir)
                dn=dn+Grid%drdi(ir)*phi*phi
             endif
          enddo
          over(ism)=dnrm
          norm(ism)=dn          
       enddo

       nrc1=nrc
       do ism=1,nsm-1
          rc1=log_rad_cutoff(shells(ism))/lambdas(ism)
          !
          nrc4=nint(dlog(rc1/Grid%b+1.0d0)/Grid%a)+1
          rc1=rc1-dlt
          rc1=log_rad_cutoff(func)-dlt !/shell%lambda(1)
          nrc3=nint(dlog(rc1/Grid%b+1.0d0)/Grid%a)+1
          rc1=Grid%r(nrc3)
          !print *, "nrc3=",nrc3
          !print *, "f size=",size(shells(ism)%f)
          if(nrc3 > size(shells(ism)%f)) nrc3=size(shells(ism)%f)
          const1=shells(ism)%f(nrc3)*Grid%r(nrc3)
          d=dsqrt(const1**2+dlt**2)
          rc1=rc1+d+dlt
          nrc2=nint(dlog(rc1/Grid%b+1.0d0)/Grid%a)+1
          d=Grid%r(nrc2)-Grid%r(nrc3)
          rc1=Grid%r(nrc2)
          const2=0.5d0*(const1**2+(d)**2)/const1
          if(nrc4.lt.nrc) nrc1=max(nrc1,nrc2)

          do ir=1, nrc2
             r=Grid%r(ir)
             if(ir.gt.nrc3.and.nrc.gt.nrc4) then
                phi=const2-const1*dsqrt(const2**2-(r-rc1)**2)/dabs(const1)
             elseif(ir .ge. size(shells(ism)%f))then
                phi=0.0_dp
             else
                phi=shells(ism)%f(ir)*Grid%r(ir)
             endif
             if(ir < nrc) g(ir)=g(ir)-over(ism)*phi/(norm(ism)+1.0d-20)
          enddo
       enddo

       if(nrc.ne.nrc1) then
          if (nrc1 < nrc) then
             nrc=nrc1
          else
             print *, "radial: multiple z nrc1 > nrc"
          endif
       endif
    end if
   
    !Store the orb divided by r**l
    do ir=2,nrc
       g(ir)=g(ir)/(Grid%r(ir)**(l+1))
    enddo

    g(1)=g(2)
    g(nrc)=0.0_dp

    call log_rad_alloc(split,g,grid)
    deallocate(g)
    if (nsm > 1) deallocate(over,norm)
  end function log_rad_parabolic_split

  !-------------------------------------------------------------------------

  function log_rad_scan_tail_parabola(rphi,l,fix_split_table,label) result(split)
!(nrc,r, drdi, l, rphi, rnrm, &
!       split_table,fix_split_table,label) 
    !
    ! Scans the whole orbital, fitting parabolas and
    ! reporting effective norms to file
    
    type(log_rad_func_t), intent(in) :: rphi
    integer, intent(in)   :: l
    logical, intent(in)   :: fix_split_table
    character(len=*), intent(in), optional :: label

    !integer, intent(in)   :: nrc  ! Index of last point of orbital
    !real(dp), intent(in),dimension(1:nrc)  :: r, drdi
    !real(dp), intent(in)  :: rphi(*), rnrm(*) ! orb, norm(r)
    

    integer :: ir, iu, jr, nrc
    character(len=50) :: fname
    real(dp) :: cons1, cons2, parab_norm, spln, rmin, rc, factor, dnrm
   
    real(dp), pointer :: rnrm(:), split_table(:), split_table_raw(:)
    real(dp), pointer :: drdi(:),r(:)

    type(log_rad_func_t) :: split
   
    drdi => rphi%grid%drdi
    r=>rphi%grid%r

    nrc = size(rphi%f)

    allocate(rnrm(1:nrc),split_table(1:nrc),split_table_raw(1:nrc))
    rnrm = 0.0_dp
    dnrm = 0.0_dp

    do ir = 1, nrc
       dnrm = dnrm + drdi(ir)*rphi%f(ir)**2
       rnrm(ir)=dnrm
    enddo

    do ir = 3, nrc-1          ! Have to avoid /0 
       call fit_parabola(ir, r, drdi,rphi%f, l, cons1,cons2,parab_norm)
       spln = 1.0_dp - rnrm(ir)
       split_table_raw(ir) =  (spln+parab_norm)
       !print ('(i,5f6.3)'), ir, rphi%f(ir), rnrm(ir), parab_norm, cons1,cons2
    enddo

    split_table_raw(2) = split_table_raw(3)
    split_table_raw(1) = split_table_raw(2)
    split_table_raw(nrc) = split_table_raw(nrc-1)

    if (fix_split_table) then
       jr = nrc - 20         ! heuristic
       rmin = r(jr)
       rc = r(nrc)
       split_table(1:jr) = split_table_raw(1:jr)
       do ir = jr+1, nrc
          factor = dampen(rmin,rc,r(ir))
          split_table(ir) = factor*split_table_raw(ir)
       enddo
    else
       split_table(1:nrc) = split_table_raw(1:nrc)
    endif

    if (present(label))then
       write(fname,"(3a,i1)") "SPLIT_SCAN.", trim(label), ".", l
       call io_assign(iu)
       open(iu,file=trim(fname),form="formatted", &
            status="unknown",action="write", position="rewind")
       do ir = 1, nrc
          write(iu,"(i4,5f14.8)") ir, r(ir), rphi%f(ir), &
               (1.0_dp - rnrm(ir)), split_table_raw(ir), split_table(ir)
       enddo
       call io_close(iu)
    endif

    call log_rad_alloc(split,split_table,rphi%grid)
    deallocate(rnrm,split_table,split_table_raw)

  end function log_rad_scan_tail_parabola

  !----------------------------------------------------------

  function dampen(a,b,r) result (y)
    real(dp), intent(in) :: a, b, r
    real(dp)             :: y

    real(dp)             :: x
    real(dp), parameter  :: tiny = 1.0e-12_dp

    x = (r-a)/(b-a)
    y = tanh(1.0_dp/(x+tiny) - 1.0_dp)
  end function dampen

  !---------------------------------------------------------

  function log_rad_self_energy(func) result(slf)
    type(log_rad_func_t), intent(in) :: func

    !C Calculates the self-energy associated to the local-pseudopotential
    !C charge density.
    !C Written by D. Sanchez-Portal, Aug. 1998.

    !real(dp), intent(out)     :: slfe
    !real(dp), intent(in)      :: vlocal(:), rofi(:), drdi(:)
    !real(dp), intent(in)      :: a
    !integer, intent(in)     :: nVna 

    !C     *Internal variables**

    real(dp) :: slf, a2b4, g0, g1, g2, g3, g4, d2g, d2u,a
    integer  :: ir,nVna
    real(dp), pointer :: vlocal(:), rofi(:),drdi(:),s(:)

    a = func%grid%a
    nVna = size(func%f)
    vlocal => func%f
    rofi => func%grid%r
    drdi => func%grid%drdi

    allocate(s(1:nVna))
    do ir=1,nVna
       s(ir)=sqrt(drdi(ir))
    enddo
    

    a2b4=0.25d0*a*a
    slf=0.0d0
    do ir=2,nVna-1
       if((ir.gt.2).and.(ir.lt.(nVna-1))) then
          g0=vlocal(ir-2)*rofi(ir-2)/s(ir-2)
          g1=vlocal(ir-1)*rofi(ir-1)/s(ir-1)
          g2=vlocal(ir)*rofi(ir)/s(ir)
          g3=vlocal(ir+1)*rofi(ir+1)/s(ir+1)
          g4=vlocal(ir+2)*rofi(ir+2)/s(ir+2)

          d2g=(16.0d0*(g1+g3)-(g0+g4)-30.0d0*g2)/12.0d0

       else
          g1=vlocal(ir-1)*rofi(ir-1)/s(ir-1)
          g2=vlocal(ir)*rofi(ir)/s(ir)
          g3=vlocal(ir+1)*rofi(ir+1)/s(ir+1)
          d2g=g1+g3-2.0d0*g2

       endif

       d2u=d2g-a2b4*g2

       slf=slf - g2*d2u*0.25d0

    enddo

    deallocate(s)
  end function log_rad_self_energy

end module radial_log_dirty
