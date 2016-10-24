      module m_localgen
      ! Generates vlocal (and chlocal)

      implicit none

      public :: compute_vlocal_chlocal, chlocal_from_vlocal

      private

      integer, parameter :: dp = selected_real_kind(10,100)

      CONTAINS

      subroutine compute_vlocal_chlocal(rofi,nrval,drdi,s,Zval,
     $                                  lmxkb, vps,
     $                                  a, b, nicore,
     $                                  nchloc,chlocal,
     $                                  vlocal,
     $                                  force_chlocal_method,
     $                                  rmax_ps_check,
     $                                  fit_3derivs,
     $                                  use_charge_cutoff,
     $                                  method_used)

      real(dp), intent(in) :: rofi(:), drdi(:), s(:)
      integer, intent(in)  :: nrval
      integer, intent(in)  :: lmxkb
      real(dp), intent(in) :: vps(:,0:)
      real(dp), intent(in) :: Zval, a, b
      character(len=4), intent(in)     :: nicore

      integer, intent(out) :: nchloc
      real(dp), intent(out):: vlocal(:)
      real(dp), intent(out):: chlocal(:)

      ! Explicitly passed options
      logical, intent(in)  :: force_chlocal_method
      real(dp), intent(in) :: rmax_ps_check
      logical, intent(in)  :: fit_3derivs
      logical, intent(in)  :: use_charge_cutoff

      character(len=*), intent(out)     :: method_used

      ! local variables
      real(dp) :: rgauss, rgauss2
      integer  :: nrgauss

! Rgauss is approximately the maximum cut-off radius used in the 
! pseudopotential generation. 
! Rgauss is determined by comparison of the pseudopot.   
! Corresponding to different l ( this is not possible if we have     
! just one pseudopotential)                                          
! Rgauss2 is the radius where the pseudopotentials reach  the        
! asymptotic behaviour 2*Zval/r.                                     
! For just one pseudopotential Rgauss is taken equal to Rgauss2       
! 
      call radii_ps(vps,rofi,Zval,nrval,lmxkb,
     $     nrgauss, rgauss, rgauss2, rmax_ps_check)
! 
!  
! Calculate local pseudopotential
! 
      if ((rgauss2.gt.1.3d0*rgauss) .and.
     $     .not. force_chlocal_method) then

              write(6,'(a)') "Using large-core scheme (fit) for Vlocal"

! In this case the atom core is so big that we do not have an asymptotic
! of 2*Zval/r until Rgauss2 (> Rc) . To retain the same asymptotic
! behaviour as in the pseudopotentials we modify the definition
! of the local potential, making it join the Vps's smoothly at rgauss.
! 
              write(6,'(/,a,f10.5)') 'atom: Estimated core radius ',
     .           rgauss2
              if (nicore.eq.'nc ')
     .          write(6,'(/,2a)') 
     .            'atom: Including non-local core corrections',
     .            ' could be a good idea'
!
! As all potentials are equal beyond rgauss, we can just use the
! s-potential here.
!
              call vlocal_as_fit(Zval,nrval,a,rofi,drdi,s,vps(:,0),
     .                     nrgauss,vlocal,nchloc,chlocal,
     $                     fit_3derivs, use_charge_cutoff) 
              method_used = "siesta-fit"
! 
            else
! 
! In this case the pseudopotential reach to an asymptotic behaviour 2*Zval/r  
! for a radius approximately equal to Rc. We build a generalized-gaussian
! "local charge density" and set Vlocal as the potential generated by
! it. Note that chlocal is negative.
! 
              write(6,'(a)') "Using small-core scheme" //
     $                       " (local charge) for Vlocal"

              call vlocal_from_chlocal(Zval, nrval, a, rofi, drdi, s,
     $                                 rgauss,vlocal,nchloc,chlocal)

              method_used = "siesta-charge"
            endif 

       end subroutine compute_vlocal_chlocal


        subroutine radii_ps(vps,rofi,Zval,nrval,lmxkb,
     .                      nrgauss, rgauss, rgauss2,
     $                      rmax_ps_check)
C
C     This routine returns the maximum radius for the
C     Kleinman-Bylander projectors with a standard choice
C     of the local potential.
C     Check also at which radius the asymptotic 2*Zval/r
C     behaviour is achieved.
C     D. Sanchez-Portal, Aug. 1998
C
        real(dp), intent(in)    ::  vps(:,0:), rofi(:)
        real(dp), intent(in)    ::  Zval
        integer,  intent(in)    ::  nrval, lmxkb
        integer,  intent(out)   ::  nrgauss
        real(dp), intent(out)   ::  rgauss, rgauss2
        real(dp), intent(in)    ::  rmax_ps_check

        real(dp) dincv, r
        integer ir, l, nrgauss2

        real(dp), parameter     ::  eps=1.0d-4
C
C Iterate over the possible local potentials
C
        rgauss=0.0_dp
        rgauss2=0.0_dp
        nrgauss=0
        nrgauss2=0


!       Optionally choose a large but not huge "infinity" for
!       the check on potentials, say 12 bohr, as passed in the
!       rmax_ps_check argument
!
        do l=0,lmxkb-1
          do ir=nrval,2,-1
            if (rofi(ir) > Rmax_ps_check) cycle
            dincv=abs(vps(ir,l)-vps(ir,lmxkb))
            if (dincv.gt.eps) exit
          enddo
          rgauss=max(rofi(ir),rgauss)
          nrgauss=max(ir,nrgauss)
        enddo
!
!       New: Use all potentials, not just l=0, since
!       potentials with larger l can converge later...
!
        do l=0,lmxkb
          do ir=nrval,2,-1
            if (rofi(ir) > Rmax_ps_check) cycle
            r=rofi(ir)
            dincv=abs(vps(ir,l)*r+2.0_dp*zval)
            if (dincv.gt.eps) exit
          enddo
          write(6,'(a,i1,a,f8.4)')
     .     'V l=', l,' = -2*Zval/r beyond r=', rofi(ir)
          rgauss2=max(rofi(ir),rgauss2)
          nrgauss2=max(ir,nrgauss2)
        enddo

        if (lmxkb.eq.0) then
          rgauss=rgauss2
          nrgauss=nrgauss2
        endif

        write(6,'(a,f8.4)') 
     .    'All V_l potentials equal beyond r=', rgauss
        write(6,'(a)')
     .    "This should be close to max(r_c) in ps generation"
        write(6,'(a,f8.4)')
     .    'All pots = -2*Zval/r beyond r=', rgauss2

        end subroutine radii_ps
            
!-----------------------------------------------------

        subroutine vlocal_from_chlocal(Zval, nrval, a, rofi, drdi, s,
     $                                 rgauss,vlocal, nchloc, chlocal)
C
C     This routine generates a smooth local pseudopotential.
C     Written by D. Sanchez-Portal, Aug. 1998
C
        real(dp), intent(in)    :: Zval, a
        integer,  intent(in)    :: nrval
        real(dp), intent(in)    :: rofi(:), drdi(:), s(:)
        real(dp), intent(out)   :: vlocal(:)
        real(dp), intent(out)   :: chlocal(:)
        integer,  intent(out)   :: nchloc
        real(dp), intent(inout) :: rgauss  ! r beyond which V_l's are equal
C
C Internal variables 
C
        real(dp) van, factor, alp, cutoff1, cutoff2,
     .    qtot, eps, pi, chc, r, Rchloc, rhor1, rhor
        integer ir
        character loctype*3
        real(dp), allocatable :: alt_chlocal(:)
        real(dp), allocatable :: chlocal_spline(:)

        parameter(eps=1.0d-4)
C
C     Usual local potential 
C     (generated with an optimum Vandebilt function)
C
        allocate(alt_chlocal(size(chlocal)))
        allocate(chlocal_spline(size(chlocal)))

        loctype='new' 
C
C     The very first local potential used by SIESTA was 
C     the electrostatic potential generated by a gaussian 
C     distribution ===> loctype='old' 
C     loctype='old'
C
        pi=acos(-1.0_dp)
C
C     Local-potential size parameter 'rgauss'
C     We choose as a smooth pseudopotential the one generated 
C     by a 'Vanderbilt-function' charge distribution. We have to select 
C     the size of this distribution somehow.
C     'Vanderbilt-functions' are of the form :
C     p(r)=N*exp(-(sinh(van*r)/sinh(van))**2)
C     when van---> 0 we will obtain a 'gaussian'
C     when van---> Inf. we will obtain a step function
C     Some test has revealed that the best election to achieve 
C     a good convergence in real and reciprocal space is b in the 
C     range 0.5-1.0 .
C
C     So, the 'gaussian' charge distribution 
C     must go to zero at a distance 'rgauss'.
C
        if (loctype.eq.'new') then          
C
C     We take a 'Vanderbilt-function' as local potential
C     van=1.0_dp all the parameter have optimized for this value 
C
          van=1.0_dp
          cutoff1=3.63_dp
          cutoff2=5.48_dp
C     99% of charge inside Rgauss
c     factor=1.627_dp

C     99.9% of charge inside Rgauss
          factor=1.815_dp
         
C       Scaling factor for local-pseudopot. charge
          alp=factor/rgauss

          write(6,'(/,a,f10.3,a)')
     .      'VLOCAL_FROM_CHLOCAL: 99.0% of the norm of Vloc inside ',
     .      (alp*cutoff1)**2,' Ry'
          write(6,'(a,f10.3,a)')
     .      'VLOCAL_FROM_CHLOCAL: 99.9% of the norm of Vloc inside ',
     .      (alp*cutoff2)**2,' Ry'

        else 

C     This is just a gaussian !!!!!!!!!!!!!!!!!                 
          van=0.00001_dp 
          rgauss=0.80_dp
          factor=2.0_dp
C     Scaling factor for local-pseudopot. charge
          alp=factor/rgauss  
        endif 
!--------------------

        ! Note that qtot is for the whole range
        qtot=0.0_dp 
        rhor1 = vander(van,alp*rofi(1))     ! This is 1...
        do ir=1,nrval
          r=rofi(ir) 
          rhor = vander(van,alp*r)
          chlocal(ir)=(-4.0_dp)*pi*rhor*r*r
          qtot=qtot+rhor*drdi(ir)*r*r
        enddo

        qtot=4.0_dp*pi*qtot 
        nchloc=0  
        do ir=nrval,1,-1
          chc=zval*chlocal(ir)/qtot
          chlocal(ir)=chc
          if ((abs(chc).gt.eps).and.(nchloc.eq.0)) then    
            nchloc=ir+1
          endif
        enddo 
        Rchloc=rofi(nchloc)

!     Note that the above cutoff is for 4*pi*r*r*rho_local(r)...
!
        ! If we now integrate chlocal up to rchloc we will not
        ! get the full charge

        call vhrtre(chlocal,vlocal,rofi,drdi,s,nrval,a)

        open(44,file="r_ch_v.charge",form="formatted")
        write(44,"(a,f12.6)") "# r  chlocal r*vlocal rchloc: ", rchloc
        do ir=2,nrval 
          r=rofi(ir)  
          chlocal(ir)=chlocal(ir)/(4.0_dp*pi*r*r)
          write(44,"(3f16.10)") r, chlocal(ir), r*vlocal(ir)
!
!     Poor man's cutoff!! Largely irrelevant?
!     This might introduce discontinuities. It is better
!     to use one of the standard tapering functions. **AG
!      
          if (r.gt.1.1_dp*Rchloc) then
            vlocal(ir)=(-2.0_dp)*zval/rofi(ir)
          endif

        enddo 
        close(44)
        chlocal(1)= -rhor1* zval/qtot

        !
        ! Test the reversibility
        !
        call ch_of_vlocal(rofi, vlocal, drdi, s, a,
     $                    nrval, alt_chlocal,
     $                    quadratic=.true.)

        open(44,file="ch_ch.charge",form="formatted")
        write(44,"(a)") "# r  chlocal_Siesta, computed"
        do ir=2,nrval 
           r=rofi(ir)  
           write(44,"(3f16.10)") r, chlocal(ir), alt_chlocal(ir)
        enddo
        close(44)

        end subroutine vlocal_from_chlocal


        subroutine vlocal_as_fit(Zval, nrval, a, rofi, drdi, s, vps,
     .                     nrgauss,vlocal,nchloc,chlocal,
     $                     fit_3derivs, use_charge_cutoff) 
C
C     This routine generates the local pseudopotential appropiate 
C     for species with  a large core.
C     Written by D. Sanchez-Portal, Aug. 1998
C     
        real(dp), intent(in)    :: Zval, a
        integer,  intent(in)    :: nrval
        integer,  intent(inout) :: nrgauss
        real(dp), intent(in)    :: rofi(:), drdi(:), s(:), vps(:)
        real(dp), intent(out)   :: vlocal(:), chlocal(:)
        integer,  intent(out)   :: nchloc

        logical, intent(in)     :: fit_3derivs
        logical, intent(in)     :: use_charge_cutoff
C
C Internal variables
C
        integer :: other_nchloc, nchloc_vlocal, nchloc_charge
        real(dp), allocatable :: other_chlocal(:) ! for debugging
        real(dp), allocatable :: alt_chlocal(:)
        real(dp), allocatable :: chlocal_spline(:)

        real(dp) 
     .    vlc, r, dev, dev2, dev3, var1, var2, var3, v1, v2, v3, v4,
     .    dm11, dm12, dm13, dm21, dm22, dm23, dm31, dm32, dm33, 
     .    g0, g1, g2, g3, g4, d2g, d2u, cons, a2b4, qtot, pi, ch
        integer 
     .    ndevfit, ir  

        real(dp), parameter  :: eps_vlocal=1.0d-5  ! this is eps
        real(dp), parameter  :: eps_charge=1.0d-4  ! for charge criterion
        real(dp) :: q1, q2
        allocate(other_chlocal(size(chlocal)))
        allocate(alt_chlocal(size(chlocal)))
        allocate(chlocal_spline(size(chlocal)))

        pi=acos(-1.0_dp)        

C Continuity up to second derivative
        if (fit_3derivs) then
           ! Continuity up to third derivative
           write(6,"(a)") "Fit of Vlocal with continuous 3rd derivative"
           ndevfit = 3
        else
           ! Continuity up to second derivative
           write(6,"(a)") "Fit of Vlocal with continuous 2nd derivative"
           ndevfit=2
        endif

        nrgauss=nrgauss+3        !! For good measure...

        do ir=1,nrval
          vlocal(ir)=vps(ir)*rofi(ir)
        enddo 

        write(6,"(a,f12.4)") "Fitting vlocal at ", rofi(nrgauss)

        ir=nrgauss
        dev=(vlocal(ir+1)-vlocal(ir-1))*0.5_dp
        dev2=(vlocal(ir+1)+vlocal(ir-1)-2.0_dp*vlocal(ir))
        dev3=(vlocal(ir+2)-2.0_dp*vlocal(ir+1)
     .     +2.0_dp*vlocal(ir-1)-vlocal(ir-2))*0.5_dp
        dev3=(dev3-3.0_dp*a*dev2+2.0_dp*(a**2)*dev)
     .     /(drdi(ir)**3)
        dev2=(dev2-a*dev)/(drdi(ir)**2)
        dev=dev/drdi(ir)

C     Local potential is Vloc(r)=v3*exp(v1*r^2+v2*r^3) 
C     inside Rgauss and equals the 
C     all-electron atomic potential outside Rgauss
C     We impose the continuity up to second derivative
      
        if (ndevfit.eq.2) then               
          vlc=vlocal(nrgauss)
          r=rofi(nrgauss)

          var1=dev/vlc-1.0_dp/r
          var2=dev2/vlc-2.0_dp*var1/r -(var1**2)

          dm11=2.0_dp*r
          dm12=3.0_dp*r*r
          dm21=2.0_dp
          dm22=6.0_dp*r
 
          v1=(dm22*var1-dm12*var2)/(6.0_dp*r*r)
          v2=(dm11*var2-dm21*var1)/(6.0_dp*r*r)
          v3=vlc/(r*exp((v1+v2*r)*r*r))

c     elseif(ndevfit.eq.3) then 
        else

C     We can also construct a local potential 
C     Vloc(r)=v4*exp(v1*r^2+v2*r^3+v3*r^4),
C     this new coefficient allows us to impose the continuity 
C     of the potential up  to the third derivative.

          vlc=vlocal(nrgauss)
          r=rofi(nrgauss)
         
          var1=dev/vlc-1.0_dp/r
          var2=dev2/vlc-2.0_dp*var1/r-(var1**2)
          var3=dev3/vlc-3.0_dp*var1*var2-(var1**3)
     .        -3.0_dp*(var1**2+var2)/r

          dm11=2.0_dp*r
          dm12=3.0_dp*r*r
          dm13=4.0_dp*r*r*r
          dm21=2.0_dp
          dm22=6.0_dp*r
          dm23=12.0_dp*r*r
          dm31=0.0_dp
          dm32=6.0_dp
          dm33=24.0_dp*r

          v1=((var1*dm22*dm33+var2*dm13*dm32+var3*dm12*dm23)
     . -(var3*dm22*dm13+var1*dm32*dm23+var2*dm12*dm33))/(48.0_dp*r*r*r)
          v2=((var2*dm11*dm33+var3*dm21*dm13+var1*dm23*dm31)
     . -(var2*dm31*dm13+var3*dm23*dm11+var1*dm21*dm33))/(48.0_dp*r*r*r)
          v3=((var3*dm11*dm22+var2*dm12*dm31+var1*dm32*dm21)
     . -(var1*dm22*dm31+var3*dm21*dm12+var2*dm11*dm32))/(48.0_dp*r*r*r)
          v4=vlc/(r*exp((v1+v2*r+v3*r*r)*r*r))
         
        endif 
      
        do ir=1,nrval
          r=rofi(ir)
          if (ir.le.nrgauss) then 
C   If second derivative fit
            if (ndevfit.eq.2) then 
              vlocal(ir)=v3*exp((v1+v2*r)*r*r)
C   If third derivative fit
            elseif(ndevfit.eq.3) then 
              vlocal(ir)=v4*exp((v1+v2*r+v3*r*r)*r*r)
            endif 
          else
            vlocal(ir)=vps(ir)
          endif 

        enddo 


C     Once we have the local potential we define the 'local-pseudopotential 
C     charge' which help us to calculate the electrostatic interation 
C     between the ions
!
!     Poisson's eq.:
!
!           1/r* d2(rV)/dr2 = -8*pi*rho
!
        a2b4=0.25_dp*a*a 
        qtot=0._dp 
        do ir=1,nrval-1
!
!        To determine the chlocal cutoff, use the reduced_vlocal cutoff
!        Comment this out for clarity (see below)
!          g2=vlocal(ir)*rofi(ir)
!          if (abs(g2+2.0_dp*zval).lt.eps) exit   !exit loop  !eps_vlocal

          if (ir.gt.nrgauss) then  

            if ((ir.gt.2).and.(ir.lt.(nrval-1))) then 
              g0=vlocal(ir-2)*rofi(ir-2)/s(ir-2)
              g1=vlocal(ir-1)*rofi(ir-1)/s(ir-1)
              g2=vlocal(ir)*rofi(ir)/s(ir)
              g3=vlocal(ir+1)*rofi(ir+1)/s(ir+1)
              g4=vlocal(ir+2)*rofi(ir+2)/s(ir+2)

              d2g=(16.0_dp*(g1+g3)-(g0+g4)-30.0_dp*g2)/12.0_dp
               
            else
              g1=vlocal(ir-1)*rofi(ir-1)/s(ir-1)
              g2=vlocal(ir)*rofi(ir)/s(ir) 
              g3=vlocal(ir+1)*rofi(ir+1)/s(ir+1)  
               
              d2g=g1+g3-2.0_dp*g2

            endif  

            d2u=d2g-a2b4*g2

            r=rofi(ir)
            cons=8.0_dp*pi*r*drdi(ir)*s(ir)
            chlocal(ir)=(-d2u)/cons
            qtot= qtot + 0.5_dp*d2u*r/s(ir)

          else

C     If second derivative fit
            if (ndevfit.eq.2)  then
              r=rofi(ir)

              g0=v3*exp((v1+v2*r)*r**2)
              g1=(2.0_dp*v1+3.0_dp*v2*r)
              g2=2.0_dp*v1+6.0_dp*v2*r
              g3=(g2+g1*g1*r*r+2.0_dp*g1)*g0
             
              cons=8.0_dp*pi
              chlocal(ir)= (-g3)/cons
              qtot= qtot  + 0.5_dp*g3*r*r*drdi(ir)
               
C If third derivative fit
               
            elseif(ndevfit.eq.3)  then

              r=rofi(ir)
               
              g0=v4*exp((v1+v2*r+v3*r*r)*r*r)     
              g1=(2.0_dp*v1+3.0_dp*v2*r+4.0_dp*v3*r*r)
              g2=(2.0_dp*v1+6.0_dp*v2*r+12.0_dp*v3*r*r)    
              g3=(g2+g1*g1*r*r+2.0_dp*g1)*g0   

              cons=8.0_dp*pi
              chlocal(ir)= -g3/cons
              qtot= qtot  + 0.5_dp*g3*r*r*drdi(ir)
            endif 

          endif
        enddo              

        ! get simply chlocal with an alternative method
        call ch_of_vlocal(rofi, vlocal, drdi, s, a,
     $                    nrval, alt_chlocal,
     $                    quadratic=.false.)

        open(44,file="ch_ch.fit",form="formatted")
        write(44,"(a)") "# r  chlocal_Siesta, computed"
        do ir=2,nrval 
           r=rofi(ir)  
           write(44,"(3f16.10)") r, chlocal(ir), alt_chlocal(ir)
        enddo
        close(44)
        

!     Decouple the different operations performed in the
!     above loop:
!      1. computation of the -Z/r-based cutoff
!      2. calculation of chlocal
!      3. computation of qtot

!     2. has now been performed above, for the whole range
!
!     1. 
!     "vlocal" option for cutoff
!     This sets the cutoff point for chlocal in a rather
!     arbitrary way, as that in which Vlocal is close to
!     -Z/r (see test inside the loop).
!
!     For consistency with the other method, the cutoff
!     should be based on the charge itself.

        do ir=1,nrval-1
           g2=vlocal(ir)*rofi(ir)
           if (abs(g2+2.0_dp*zval).lt.eps_vlocal) exit !exit loop
        enddo
        nchloc_vlocal = ir          

!
!     "charge" option for cutoff
!
        do ir = nrval, 1, -1
           r = rofi(ir)
           ch = 4*pi*r*r*chlocal(ir)
           if (abs(ch) .gt. eps_charge) exit
        enddo
        nchloc_charge = ir

        ! Choose
        if (use_charge_cutoff) then
           nchloc = nchloc_charge
           write(6,"(a,i3,f10.6)") "Choosing charge chloc cutoff:",
     $          nchloc, rofi(nchloc)
        else
           ! classic behavior
           nchloc = nchloc_vlocal
           write(6,"(a,i3,f10.6)") "Choosing vlocal chloc cutoff:",
     $          nchloc, rofi(nchloc)
        endif

        ! Now compute qtot only up to nchloc
        qtot = 0.0_dp
        do ir=1,nchloc-1
           r = rofi(ir)
           ch = 4*pi*r*r*chlocal(ir)*drdi(ir)
           qtot = qtot - ch
        enddo
        write(6,"(a,f14.8)") "qtot up to nchloc:", qtot

        do ir=1,nchloc-1
          chlocal(ir)=zval*chlocal(ir)/qtot
        enddo  
        do ir=nchloc,nrval
          chlocal(ir)=0.0_dp
        enddo 

        end subroutine vlocal_as_fit

      subroutine chlocal_from_vlocal(rofi, vlocal, drdi, s, a, nrval,
     $                               Zval, nchloc, chlocal,
     $                               use_charge_cutoff,quadratic)

      ! This is a general version, which does not use
      ! the details of the fit for r < rgauss
      ! It should work for an arbitrary vlocal

        real(dp), intent(in)    :: rofi(:), drdi(:), s(:), vlocal(:)
        real(dp), intent(in)    :: a, Zval
        integer,  intent(in)    :: nrval
        real(dp), intent(out)   :: chlocal(:)
        integer,  intent(out)   :: nchloc
        ! Option for compatibility
        logical, intent(in)  :: use_charge_cutoff
        logical, intent(in)  :: quadratic

        ! local variables

        real(dp), parameter  :: eps=1.0d-5
        real(dp), parameter  :: eps_charge=1.0d-4 ! To match
        real(dp) :: qtot, charge
        integer  :: ir
        real(dp), parameter :: pi = 3.14159265358979323846_dp

        integer :: status

C     Once we have the local potential we define the 'local-pseudopotential 
C     charge' which help us to calculate the electrostatic interation 
C     between the ions


        call ch_of_vlocal(rofi, vlocal, drdi, s, a,
     $                  nrval, chlocal, quadratic)

!
!       This is the classic way of computing an effective cutoff
!       for chlocal, but it is only a heuristic based on when
!       rVlocal equals 2*Zval. There are other criteria, including
!       a specific test for smallness of 4pi*r*r*chlocal, as in
!       the "vlocal from a localized charge" routine
!
        do ir = 1, nrval
           if (abs(rofi(ir)*vlocal(ir)+2.0_dp*zval).lt.eps) exit   !exit loop
        enddo
        nchloc = ir
        write(6,"(a,i4,f10.6)") "nchloc from rV+2Z: ",
     $       nchloc, rofi(nchloc)

        if (use_charge_cutoff) then
         !     Alternative way to compute nchloc,
         !     used by Siesta when method="charge"
           do ir = nrval, 1, -1
              charge =  4*pi*rofi(ir)**2 * chlocal(ir)
              if (abs(charge).gt.eps_charge) exit
           enddo
           nchloc=ir+1
           write(6,"(a,i4,f10.6)") "nchloc from 4*pi*r*r*chloc: ",
     $          nchloc, rofi(nchloc)
        endif

        qtot = 0.0_dp
        do ir = 1, nchloc-1
           qtot = qtot - 4*pi*rofi(ir)**2 * chlocal(ir) * drdi(ir)
        enddo
        print *, "qtot in chlocal_from_vlocal up to nchloc: ", qtot

        do ir=1,nchloc-1
          chlocal(ir)=zval*chlocal(ir)/qtot
        enddo  
        do ir=nchloc,nrval
          chlocal(ir)=0.0_dp
        enddo 

        end subroutine chlocal_from_vlocal

      subroutine ch_of_vlocal(rofi, vlocal, drdi, s, a,
     $                        nrval, chlocal, quadratic)

      use flib_spline
      use m_interpol, only: interpolate

      ! This is a general version, which
      ! should work for an arbitrary vlocal

      ! It uses some help from the user, if she knows that the
      ! behavior of chlocal is quadratic near the origin

        real(dp), intent(in)    :: rofi(:), drdi(:), s(:), vlocal(:)
        real(dp), intent(in)    :: a
        integer,  intent(in)    :: nrval
        real(dp), intent(out)   :: chlocal(:)
        logical, intent(in)     :: quadratic

        ! local variables

        real(dp) :: pi, r
        real(dp) :: rmax, delta, a1, b1, x2, splval
        real(dp) :: delta_1, factor
        integer  :: ir, npts, ninitial, nstep

        real(dp), allocatable :: f(:), f2(:)
        real(dp), allocatable :: ralt(:), drd1(:), s1(:),
     $                           chlocal_ralt(:)
        real(dp), allocatable :: chlocal_ralt_spl(:)
        real(dp), allocatable :: g(:), g2der(:), g2der_spl(:), ch2(:)
        real(dp), allocatable :: xaux(:)

        integer :: status

C     Once we have the local potential we define the 'local-pseudopotential 
C     charge' which help us to calculate the electrostatic interation 
C     between the ions
!
!     Poisson's eq.:
!
!           1/r* d2(rV)/dr2 = -8*pi*rho
!
!     Important note:
!   
!     Near the origin, rho goes like -2*V'/r - V'', so V' must be zero at
!     the origin in order for V to be rho-representable!
!
!     If V(r) = V_0 + a*r^2 + b*r^3 + c*r^4  near the origin,
!
!     8*pi*rho(r) = -6a - 18b*r -20c*r^2 ...
!
!     We could check:
!
!     1. Whether V' is indeed zero near the origin
!     2. Try to estimate a (and maybe b) to help with the extrapolation
!        of rho near the origin.
!

        pi=acos(-1.0_dp)        

        allocate (f(nrval), f2(nrval))
        do ir = 1, nrval
           f(ir) = rofi(ir)*vlocal(ir)
           write(44,*) rofi(ir), f(ir)
        enddo

        ! re-sample the range to avoid very small differences
!        Example for linear grid:
!        rmax = 10.0_dp
!        delta= 0.0002_dp
!        npts = rmax/delta + 1
!        npts = 1000
!        ----
!        Empirically it is better to use a logarithmic grid, not so
!        fine as the standard atomic one, and with reduced range.
!        
!        I want Rmax=10 au,
         rmax = 10.0_dp
!        Initial delta_r
         delta_1 = 2e-4_dp
!        expansion factor: 1 to 500:  1.e-1 at the end of the range
         factor = 500.0_dp

!        (note exchanged a, b with respect to Siesta usage...)
!        We have, approximately:  exp(bN) = factor
!                                 ab = delta_1
!                                 a*exp(bN) = Rmax
!        so
!
         a1 = rmax / factor
         b1 = delta_1 / a1
         npts = log(factor) / b1

        allocate(ralt(npts), drd1(npts), s1(npts),
     $           chlocal_ralt(npts))
        allocate(g(npts), g2der(npts), ch2(npts))
        allocate(g2der_spl(npts), chlocal_ralt_spl(npts))


        ! Two ways of resampling, basically equivalent
        ! Splines, and polynomial interpolation

        call generate_spline(rofi,f,nrval,f2,stat=status)
        if (status/=0) call die("error in generate_spline")

        open(55,file="r.rV.rVspl",form="formatted")
        write(55,"(a)") "# ralt, r*Vlocal (interp), r*Vlocal (spline)"
        do ir = 1, npts
           ralt(ir) = a1*(exp(b1*(ir-1)) -1)
           drd1(ir) = a1*b1*(exp(b1*(ir-1)))
           s1(ir) = sqrt(drd1(ir))
!           ralt(ir) = (ir-1)*delta   ! for linear grid
           call interpolate(rofi(1:nrval),f,ralt(ir),g(ir),npoint=2)
           call evaluate_spline(rofi,f,f2,nrval,ralt(ir),splval)
           write(55,*) ralt(ir), g(ir), splval
        enddo
        close(55)

        ! Two ways of computing the second derivative
        call compute_2nd_der(ralt,b1,drd1,s1,g,npts,g2der,g2der_spl)

        ! Should we trust the second derivative for the points
        ! near r=0? No.

        open(55,file="r.ch.ch",form="formatted")
        write(55,"(a)") "# ralt, chlocal_ralt_fd, chlocal_ralt_spline"
        do ir = 2, npts
           chlocal_ralt(ir) = - g2der(ir) / (ralt(ir)*8*pi)
           chlocal_ralt_spl(ir) = - g2der_spl(ir) / (ralt(ir)*8*pi)
           write(55,*) ralt(ir), chlocal_ralt(ir), chlocal_ralt_spl(ir)
        enddo
        close(55)

        ! Now we have chlocal on the alternate grid
        ! re-sample back to the original logarithmic grid

        call generate_spline(ralt(2:npts),chlocal_ralt_spl(2:npts),
     $                       npts-1,ch2,stat=status)
        if (status/=0) call die("error in generate_spline")

        if (quadratic) then
           allocate(xaux(npts))
           xaux(:) = ralt(:)**2
        endif

        ! Experimentally, it is better to use the spline data
        !
        do ir = 1, nrval
           if (rofi(ir) > rmax) then
              chlocal(ir) = 0.0_dp
           else if (rofi(ir) > ralt(5)) then
              call evaluate_spline(ralt(2:npts),
     $              chlocal_ralt_spl(2:npts),ch2,npts-1,
     $              rofi(ir),chlocal(ir))
           else
              ! The region near r=0 is tricky
              ! We have to use some heuristics
              if (quadratic) then
                 ninitial = 4
                 nstep = 2
                 x2 = rofi(ir)**2
                 call interpolate(xaux(ninitial:npts:nstep),
     $                         chlocal_ralt_spl(ninitial:npts:nstep),
     $                         x2,chlocal(ir),npoint=2)
              else
                 ninitial = 5
                 nstep = 3
                 call interpolate(ralt(ninitial:npts:nstep),
     $                            chlocal_ralt_spl(ninitial:npts:nstep),
     $                            rofi(ir),chlocal(ir),npoint=3)
              endif
           endif
        enddo

        CONTAINS

        subroutine compute_2nd_der(ralt,b1,drd1,s1,g,npts,
     $                             g2der,g2der_spl)
        real(dp), intent(in) :: ralt(:), drd1(:), s1(:), g(:)
        real(dp), intent(in) :: b1
        integer, intent(in)  :: npts
        real(dp), intent(out):: g2der(:)
        real(dp), intent(out):: g2der_spl(:)

        integer :: status, ir
        real(dp) :: d2g, d2u, g1, g2, g3, a2b4

        call generate_spline(ralt,g,npts,g2der_spl,stat=status)
        if (status/=0) call die("error in generate_spline")

        ! This section for finite-difference derivative
        ! (it seems to have problems near the origin)

        a2b4 = 0.25_dp * b1*b1  ! note exchange a-b with respect to Siesta

        do ir = 2, npts-1
           g1=g(ir-1)/s1(ir-1)
           g2=g(ir)/s1(ir) 
           g3=g(ir+1)/s1(ir+1)  
               
           d2g=g1+g3-2.0_dp*g2
           d2u=d2g-a2b4*g2
           g2der(ir) = d2u / (drd1(ir)*s1(ir))
        enddo
        g2der(npts) = 0.0_dp

        end subroutine compute_2nd_der

        end subroutine ch_of_vlocal

!--------------------------------------------------------------
!
!      The famous "Vanderbilt generalized cutoff"
!
       function vander(a,x) result(f)
       real(dp), intent(in) :: a    ! Generalized gaussian shape
       real(dp), intent(in) :: x    
       real(dp)             :: f

       real(dp), parameter :: log10_e = 0.4343
       real(dp), parameter :: exp_range = (range(1.0_dp)-1)/log10_e

!!     real(dp), parameter :: exp_range = 40.0_dp
       real(dp)   :: gexp

       gexp = sinh(a*x)/sinh(a)
       gexp = gexp*gexp

       if (gexp .lt. exp_range) then
          f=exp(-gexp)  
       else
          f = 0.0_dp
       endif

       end function vander


      end module m_localgen