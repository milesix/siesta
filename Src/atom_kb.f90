!     This file is part of the SIESTA package.
!     
!     Copyright (c) Fundacion General Universidad Autonoma de Madrid:
!     E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal
!     and J.M.Soler, 1996-2006.
!     
!     Use of this software constitutes agreement with the full conditions
!     given in the SIESTA license, as signed by all legitimate users.
!     

module atom_kb
  !
  !    Generates the KB projectors
  use precision, only : dp
  use atom_types, only : get_number_of_kb_projs, &
       init_kb_proj,set_kb_projs_deg,set_lmax_kb_proj, set_kb_proj
  use hilbert_vector_m, only : hilbert_vector_t,  init_vector, destroy_vector
  use radial, only : rad_func_t, rad_schro, rad_kbproj,rad_integral_vs_value, &
       rad_energy_deriv, rad_cutoff, rad_ghost, rad_dump_file, rad_log_to_linear, &
       rad_dealloc
  use pseudopotential_new, only : pseudopotential_new_t, get_pseudo_down, get_ve_val,&
       get_vlocal, get_pseudo_valence_charge
  use atom_generation_types, only : basis_def_t, kbshell_t, basis_parameters
  use sys, only : die
  implicit none

  private
  public kbgen

CONTAINS

  subroutine KBgen(isp)
    

    !    Call routines for 1) the generation of the Kleinman-Bylander projectors,
    !    2) Cheking for the presence of ghost states and 3), the storage of 
    !    all the information in the corresponding common blocks.
    !    
    !    Written D. Sanchez-Portal, Aug. 1998.
    !    Modified by DSP to allow more than one projector per l, July 1999.

    integer,      intent(in) :: isp

    type (basis_def_t),pointer :: basp
    type (kbshell_t),pointer   :: k
    type (pseudopotential_new_t), pointer :: vps
    type (hilbert_vector_t)      :: kb_vector
    type (rad_func_t)            :: kb_proj
    type (rad_func_t), pointer   :: rphi(:)
    type (rad_func_t)            :: vdown, ve_val, vlocal

    integer :: l,nkb,nnodes,nprin, ighost, ikb,i, nkb_deg, idx

    real(dp) :: rmax, pop, r_c, zval
    real(dp), pointer :: ekb(:), dkbcos(:), rc(:)

    !   The atomic wavefunctions and/or its energy derivatives are* 
    !  calculated only inside a sphere of radius Rmax. To define the***  
    !   KB projectors they will not be need very far from the nucleus,** 
    !   and this limitation simplifies the handling of not bound states*
    !    
    parameter (Rmax=6.0d0)
    !    
    ighost = 0

    basp => basis_parameters(isp)
    vps => basp%pseudopotential

    nkb=0
    nkb_deg=0
    do l=0,basp%lmxkb
       k => basp%kbshell(l)
       do ikb=1,k%nkbl
          nkb_deg=nkb+(2*l+1) 
          nkb=nkb+1
       enddo
    enddo

    allocate(ekb(1:nkb),dkbcos(1:nkb),rc(1:nkb))

    call init_kb_proj(isp,nkb)
    call set_lmax_kb_proj(isp,basp%lmxkb)

    print '(/a)', "--------------------------------------------"
    print *, "KB: Generation of KB projectors"

    idx=1
    do l=0,basp%lmxkb
       print *, "KB:  L=",l
       k => basp%kbshell(l)
       if (k%nkbl > 0) allocate(rphi(1:k%nkbl))

       print '(a,i2)', " KB:    Number of Kleinman-Bylander projectors:",k%nkbl 
       do ikb= 1, k%nkbl
          print '(A,i2)', " KB:       Generating projector:",ikb
          !    Atomic wavefunctions and eigenvalues
          !    for the construction of the KB projectors

            vdown = get_pseudo_down(vps,l)
            ve_val = get_ve_val(vps) 
            zval = get_pseudo_valence_charge(vps)
            vlocal = get_vlocal(vps)

          if(k%erefkb(ikb).ge.1.0d3) then           
             !    If the reference energies have not been specifed, the eigenstates
             !    with the condition of being zero at r(nrval) will be used.
             !   
             print *, "KB:          Projector kind: standard"
             nnodes=ikb
             nprin=l+1
                        
             r_c = rad_cutoff(vdown)
           
             call rad_dump_file(vdown,"vdown.dat")
             call rad_dump_file(ve_val,"ve_val.dat")
             print *, "rc,l,nnodes,nprint,zval,e=",r_c,l,nnodes,nprin,zval,k%erefkb(ikb)

             rphi(ikb) = rad_schro(vdown,ve_val,r_c,l,nnodes,nprin,zval,k%erefkb(ikb))    
             

          elseif((k%erefkb(ikb).le.-1.0d3).and.(ikb.gt.1) ) then 
             print *, " KB:          Projector generated using the energy derivative"
             print *, " KB:          of the previous (L-1) projector"

             !    If the energy is specified to be 1000 Ry, the energy derivative
             !    of the previous wavefunction will be used

             rphi(ikb) = rad_energy_deriv(rphi(ikb-1),vdown,ve_val,l,k%erefkb(ikb-1))
             k%erefkb(ikb)=0.0d0

          else 

             print *, " KB:          Projector generated using the specified reference energy"
             !    If the reference energies have been specified, we just use them 

             rphi(ikb) = rad_integral_vs_value(vdown,ve_val,l,k%erefkb(ikb),rmax)

          endif

          !    
          !***  GHOST ANALYSIS****
          !    
          if(k%nkbl.eq.1) then
             call rad_ghost(rphi(ikb),vdown,vlocal,ve_val,l,k%erefkb(ikb),zval,ighost)

                if (ighost.eq.1) then
                    write(6,"(2a)")'KBgen: WARNING: ','Ghost states have been detected'
                    write(6,"(2a)")'KBgen: WARNING: ','Some parameter should be changed in the '
                    write(6,"(2a)")'KBgen: WARNING: ','pseudopotential generation procedure.'
                    call die()
                endif

          else 
             if (ikb.eq.1) &
                  write(6,'(a,i3,/a)') '      KBgen: More than one KB projector for l=',l,&
                  '      KBgen: ghost states analysis will be not performed'
          endif

          kb_proj=rad_kbproj(rphi(ikb),vdown,vlocal,l,dkbcos(idx),ekb(idx))
          
          rc(idx) = rad_cutoff(kb_proj)

          pop=0.0_dp
          call init_vector(kb_vector,kb_proj,0,l,ikb,pop,ekb(idx),.false.) 
          call set_kb_proj(isp,kb_vector,idx)
          call rad_dealloc(kb_proj)
          call destroy_vector(kb_vector)
          idx=idx+1

          call rad_dealloc(vdown)
          call rad_dealloc(vlocal)
          call rad_dealloc(ve_val)
       enddo !ikb

       if (k%nkbl > 0) then
          do ikb=1,size(rphi)
             call rad_dealloc(rphi(ikb))
          end do
          deallocate(rphi)
       endif
    enddo
    
    if (ighost.eq.1) then
       write(6,"(2a)")'KBgen: WARNING: ','Ghost states have been detected'
       write(6,"(2a)")'KBgen: WARNING: ','Some parameter should be changed in the '
       write(6,"(2a)")'KBgen: WARNING: ','pseudopotential generation procedure.'
       call die
    endif


    write(6,'(/,a)')'KBgen: Kleinman-Bylander projectors: '
    ikb=1
    do l=0,basp%lmxkb
       k => basp%kbshell(l)
       do i=1, k%nkbl
          write(6,'(3x,a,i2,4(3x,a,f10.6))')&
               'l=',l, 'rc=',rc(ikb), 'el=',k%erefkb(i),& 
               'Ekb=',ekb(ikb),'kbcos=',dkbcos(ikb)
          ikb=ikb+1
       enddo
    enddo
    call set_kb_projs_deg(isp)

    write(6,'(/,a, i4)') 'KBgen: Total number of  Kleinman-Bylander projectors: ', &
         get_number_of_kb_projs(isp)

     print '(/a)', "--------------------------------------------"

    deallocate(ekb,dkbcos,rc)
  end  subroutine KBgen

  !--------------------------------------------------

endmodule atom_kb
