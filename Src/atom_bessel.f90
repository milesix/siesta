!     This file is part of the SIESTA package.
!     
!     Copyright (c) Fundacion General Universidad Autonoma de Madrid:
!     E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal
!     and J.M.Soler, 1996-2006.
!     
!     Use of this software constitutes agreement with the full conditions
!     given in the SIESTA license, as signed by all legitimate users.
!     

!Generation of multiple-z using bessel functions

module atom_bessel
use precision
use atom_generation_types , only : basis_def_t, basis_parameters,lshell_t, shell_t
use atomparams
use radial
use atom_types 
use hilbert_vector_collection
use atom_pao_util, only : number_of_orbs
use schro, only : schro_eq
implicit none

public bessel
private

contains

  subroutine bessel(isp) 
    !****
    ! Caculates Bessel functions as a floating basis
    ! Written by D. Sanchez-Portal, Aug. 1998.
    ! Modify by DSP, July 1999

    implicit none

    integer, intent(in) :: isp

    type(basis_def_t), pointer :: basp     !Structure with the parameters corresponding to this
    !basis set.
    
    type(shell_t),     pointer :: shell    !Pointer to a shell.
    type(lshell_t),    pointer :: lshell   !Pointer to a l-shell.


    type(rad_func_t)       :: rad_orb 
    type(hilbert_vector_t) :: orb_vector
    type(rad_grid_t)        :: grid


    !***Internal variables**

    integer ::l,nprin, nnodes, nodd, nrc, ir,iorb,izeta,lmxo, nsm, norbs 

    real(dp) :: rc,dnrm, phi,g(nrmax), eorb, eps,a,b

    real(dp) :: v(nrmax)=0.0_dp,s(nrmax)=0.0_dp,drdi(nrmax)=0.0_dp,rofi(nrmax)=0.0_dp

    basp => basis_parameters(isp)

    !

    a=0.247875217667E-02 
    b=0.125000000000E-01

    ! ARRAY S FOR THE SCHRODINGER EQ. INTEGRATION**
    ! (Needed for the calculation of the basis set )
    call set_mesh(a,b,rofi,drdi,s)

    grid = rad_grid_alloc(size(rofi),a=a,b=b,r=rofi)

    s(2:nrmax) = drdi(2:nrmax)*drdi(2:nrmax)
    s(1)       = s(2)

    norbs = number_of_orbs(isp)

    call init_orbs(isp,norbs)
        
    iorb=1
    lmxo=get_lmax_orbs(isp) 
    do l=0,lmxo 

       lshell => basp%lshell(l)
       do nsm=1,lshell%nn          
          shell => lshell%shell(nsm)
          if(shell%nzeta.gt.0) then

             write(6,'(/2A,I2)') 'Bessel: floating Bessel functions ', &
                  'with angular momentum L=',l
             
             do izeta=1, shell%nzeta 

                !**Cut-off radius for Bessel functions must be an explicit input***

                if (shell%rc(izeta).lt.1.0d-5) then 
                   write(6,'(a)')  'Bessel: ERROR Zero cut-off radius with Z=0 option'
                   write(6,'(a)')  'Bessel: ERROR Cut-off radius must be explicitely specified'
                   write(6,'(a)')  'Bessel: ERROR using Z=0 (Floating Bessel functions)'
                   call die

                endif

                if(dabs(shell%lambda(izeta)) < 1.0d-3) shell%lambda(izeta)=1.0d0
                if( (dabs(shell%lambda(izeta))-1.0d0) >1.0d-3) then
                   write(6,'(/,a)')  'Bessel: WARNING Scale factor is not active with Z=0 option' 
                endif
                shell%lambda(izeta)=1.0d0
                !C*

                rc=shell%rc(izeta)
                nrc=nint(log(rc/b+1.0d0)/a)+1
                nodd=mod(nrc,2)
                if(nodd.eq.0) then
                   nrc=nrc+1
                endif
                rc=b*(exp(a*(nrc-1))-1.0d0)
                shell%rc(izeta)=rc

                nnodes=izeta
                nprin=l+1
                call schro_eq(1.0d0,rofi,v,v,s,drdi, &
                     nrc,l,a,b,nnodes,nprin,eorb,g) 
                dnrm=0.0d0
                do ir=2,nrc
                   phi=g(ir)
                   dnrm=dnrm+drdi(ir)*phi*phi
                   g(ir)=phi/rofi(ir)
                enddo
                g(1)=g(2)        

                !****Checking normalization of the wavefunctions**
                eps=1.0d-4
                if(abs(dnrm-1.0d0).gt.eps) then
                   do ir=1,nrc
                      g(ir)=g(ir)/sqrt(dnrm)
                   enddo
                endif
                !****  

                call rad_alloc(rad_orb,g(1:nrc),grid)
                

                write(6,'(/,(3x,a,i2),2(/,a25,f12.6))') 'izeta =',izeta, &
                     'rc=',shell%rc(izeta),'energy =',eorb  

                !We are done. Fill species with the orbs.

                call init_vector(orb_vector,rad_orb,shell%n,l,izeta,0.0_dp,0.0_dp,.false.)               
                call set_orb(isp,orb_vector,iorb)
                
                iorb=iorb+1

             enddo

          endif

       enddo

    enddo

    call set_orbs_deg(isp)

    write(6,'(/a,i3)') 'ATOM: Total number of floating Bessel orbitals:', get_number_of_orbs(isp)

    

  end subroutine bessel

!------------------------------------------------------------------------------------

  !---------------------------------------------------------------
  !
  subroutine set_mesh(a,b,rofi,drdi,s)

    !**
    !   Setting up mesh points an its derivatives from standard
    !   values
    !   D. Sanchez-Portal, Aug. 98

    implicit none

    real(dp) rofi(nrmax), drdi(nrmax), s(nrmax), a, b


    !**** Internal variables**
    real(dp) aa, bb, zt, rpb, ea, ea2
    integer ir

    parameter(zt=1.0d0)
    parameter(aa=80.0d0)
    parameter(bb=6.0d0) 


    !*STANDART VALUES FOR MESH PARAMETERS

    b=exp(-bb)/zt
    a=1.0d0/aa

    !*SET UP THE MESH POINTS AND ITS DERIVATIVE***

    rpb=b
    ea=exp(a)
    ea2=1.0d0
    do ir=1,nrmax
       drdi(ir)=a*rpb
       rofi(ir)=b*(ea2-1.0d0)
       s(ir)=(a*rpb)**2
       rpb=rpb*ea
       ea2=ea2*ea
    enddo
  end subroutine set_mesh

  !---------------------------------------------------------------------

end module atom_bessel
