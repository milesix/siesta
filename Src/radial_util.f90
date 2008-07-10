module radial_util
use precision
implicit none

public

contains

  subroutine findp(nrc,nm,rphi,a,b,l,cons1,cons2,bug)
    integer, intent(in)   ::  nrc, nm, l
    real(dp), intent(in)  ::  a, b
    real(dp), intent(out) ::  cons1, cons2
    real(dp), intent(in)  ::  rphi(nrc)
    logical, intent(in) :: bug

    !  This routine provides the constants Cons1 and 
    !  Cons2 and described in subroutine 'parabola' 
    !  for fitting at point of index nm

  

    real(dp) rm, rm1, rm2, drdi_local, frsp
    real(dp) dfrsp

     if (bug) then
       ! old, wrong code
       rm=b*(exp(a*(nm-1)) + 1.0d0) 
       rm1=b*(exp(a*(nm-2)) + 1.0d0)
       rm2=b*(exp(a*nm) + 1.0d0)
    else
       ! correct code
       rm=b*(exp(a*(nm-1)) - 1.0d0)     
       rm1=b*(exp(a*(nm-2)) -  1.0d0)   
       rm2=b*(exp(a*nm) - 1.0d0)        
    endif

    drdi_local =a*b*exp(a*(nm-1))

    frsp=rphi(nm)/rm
    dfrsp=0.5d0*(rphi(nm+1)/rm2-rphi(nm-1)/rm1)
    dfrsp=dfrsp/drdi_local

    cons1= 0.5d0*(dfrsp*rm-l*frsp)/(rm**(l+2))
    cons2= frsp/(rm**l)-cons1*(rm**2)

  end subroutine findp

  !--------------------------------------------------------

    subroutine nrmpal(c1,c2,r,l,dnrm)
    real(dp), intent(in)   ::  c1, c2, r
    integer, intent(in)  ::  l
    real(dp), intent(out)  ::  dnrm

    !    returns the norm of a parabolic function
    !    f(r')= r'^l (c1*r'^2 + c2)  r'< r
    !    0 otherwise
    !    

    dnrm=(c1**2)*r**(2*l+7)/(2*l+7) & 
         + (2.0d0*c1*c2)*r**(2*l+5)/(2*l+5) &
         + (c2**2)*r**(2*l+3)/(2*l+3)

  end subroutine nrmpal

  !--------------------------------------------------------------------------

    subroutine fit_parabola(nsp,r, drdi,rphi, l, cons1,cons2,parab_norm)

    ! Fits C1 and C2 to match a parabola to rphi at
    ! point of index nsp

    integer, intent(in)   :: nsp
    real(dp), intent(in)  :: r(*), drdi(*)
    real(dp), intent(in)  :: rphi(*)
    integer, intent(in)   :: l
    real(dp), intent(out) :: cons1, cons2
    real(dp), intent(out) :: parab_norm


    real(dp) :: rsp, frsp, dfrsp

    rsp = r(nsp)
    frsp=rphi(nsp)/rsp
    dfrsp=0.5d0*(rphi(nsp+1)/r(nsp+1) &
         -rphi(nsp-1)/r(nsp-1))
    dfrsp=dfrsp/drdi(nsp)

    cons1= 0.5d0*(dfrsp*rsp-l*frsp)/(rsp**(l+2))
    cons2= frsp/(rsp**l)-cons1*rsp**2
    call nrmpal(cons1,cons2,rsp,l,parab_norm)

  end subroutine fit_parabola


end module radial_util
