module ldauproj
  use precision,        only: dp
  use fdf
  use hilbert_vector_collection
  use hilbert_vector_m, only:hilbert_vector_t
  use radial,           only: rad_func_t, rad_grid_t, rad_copy, &
        rad_get_value_from_ir, rad_get_r_from_ir, rad_get_ir_from_r, &
        rad_get_length, rad_get_norm_r_l_1_from_ir, rad_get_grid, &
        rad_multiply, rad_alloc
  use atom_types,        only:  get_lmax_ldau_proj, &
        set_lmax_ldau_proj, set_lmax_orbs, init_ldau_proj, &
        get_atomic_number, set_ldau_proj, set_ldau_projs_deg, & 
        get_number_of_ldau_proj, get_number_of_ldau_proj_non_deg, &
        set_switch_ldau, get_switch_ldau, init_UJ_ldau_proj, &
        set_U_ldau_proj, set_J_ldau_proj, get_U_ldau_proj, get_J_ldau_proj
  
  use pseudopotential_new,     only:pseudopotential_new_t, get_ve_val, get_ve_val_scaled, &
       get_pseudo_down
  use atom_generation_types,   only:basis_def_t,basis_parameters, ldaushell_t, l_ldaushell_t
  use sys,              only: die
  use atom_pao_util,         only: sym , generate_vsoft, soft_confinement
  use atom_options, only: write_ion_plot_files
  use units, only: eV
  implicit none

  private
  public ldauproj_gen

  real(dp), parameter             :: deltmax=0.05d0

  !     Default energy-shift to define the cut off radius of orbitals
  ! In Rydbergs
  real(dp), parameter             :: eshift_default=0.05d0

contains



  subroutine ldauproj_gen(isp)
    !Generation of LDAU projectors
    ! LDAU projectors are, basically, PAOs with artificially small radii 
    !    Written by D. Sanchez-Portal, Aug. 2008 after module
    !    basis_gen



    integer, intent(in)        :: isp
    

    !***  Internal variables**

    integer :: l,iorb,norbs
    integer :: izeta, nsm, nvalence, maxn, nmin
    integer :: nldaupj
    logical :: switch
    integer :: method
 
    type(basis_def_t), pointer :: basp     !Parameters corresponding to this basis set.
    type(pseudopotential_new_t),pointer :: vps !Psuedopotential info.
    type(ldaushell_t),     pointer :: shell    !Pointer to a shell.
    type(l_ldaushell_t),    pointer :: lshell   !Pointer to a l-shell.
    
    type(hilbert_vector_t) :: orb_vector
    type(rad_func_t)   :: rad_tmp,  vdown
    

    basp => basis_parameters(isp)
    vps  => basp%pseudopotential

    nldaupj = number_of_ldau_proj(isp)
    call init_ldau_proj(isp,nldaupj)
    call set_lmax_ldau_proj(isp,basp%lmxldaupj)

    switch=.false.
    if(nldaupj.gt.0) switch=.true. 
    call set_switch_ldau(isp,switch)
    if(.not.switch) return
    maxn=1
    do l=0,basp%lmxldaupj
        lshell => basp%l_ldaushell(l)
        do nsm=1,lshell%nn
            shell => lshell%ldaushell(nsm)
   !Determine max principal quantum number to allocate
   ! enough space to store all possible different U and J parameters
            maxn=max(maxn,shell%n)
        enddo
    enddo

    call init_UJ_ldau_proj(isp,maxn)
    
    print ('(/,a)'), "LDAUprojgen begin"
    
    !LOOP over angular momenta    
    iorb = 1

    do l=0,get_lmax_ldau_proj(isp) ! species(isp)%lmax_basis
       lshell => basp%l_ldaushell(l)
       nmin=100
       do nsm=1,lshell%nn
          shell => lshell%ldaushell(nsm)
          nmin=min(nmin,shell%n)
       enddo 
       do nsm=1,lshell%nn
          
          shell => lshell%ldaushell(nsm)
          shell%i_sm = shell%n-nmin+1

          call set_U_ldau_proj(isp,shell%U,l,shell%n)
          call set_J_ldau_proj(isp,shell%J,l,shell%n)

            write(6,'(/A,I2)')'LDAUprojs with angular momentum L=',l

          !IF THE COMPRESSION FACTOR IS NEGATIVE OR ZERO THE ORBITALS ARE
          !LEFT UNTOUCHED
          if(shell%lambda.le.0.0d0) shell%lambda=1.0d0
          vdown = get_pseudo_down(vps,l)
          call rad_copy(vdown, shell%pseudo)

          shell%z_valence  = basp%ground_state%z_valence
          !Set pseudo for schroed. eq.
          call pseudopotential_setup(shell,vps)

          method=shell%method
          write(6,'(a,i4)')'LDAUproj generation method', method
          if(method.eq.1) then 
            !Automatic rc?
            if (shell%rc .eq. 0.0_dp) call auto_rc(shell,vps,l,shell%i_sm)
            !Generate the orbital    
            call generate_first_zeta(shell,basp%label)
          elseif(method.eq.2) then 
            call generate_long_first_zeta(shell,basp%label)    
            call fermicutoff(shell)
            write(6,'(a)')'LDAUproj is an extended PAO orbital cut off with a'
            write(6,'(a)')'Fermi function 1/[1+exp(r-rc)/w] with' 
            write(6,'(a,f12.6)') ' rc=',shell%rc
            write(6,'(a,f12.6)') ' w =',shell%width
            write(6,'(a,f12.6)')     &
            'LDAUproj cutoff radious ',rad_cutoff(shell%ldaupj) 
          else
            stop "LDAU: Projector Generation Method Incorrect" 
          endif
           
          call init_vector(orb_vector,shell%ldaupj,shell%n,l,1,&
                  0.0_dp,0.0_dp,.false.)


          call set_ldau_proj(isp,orb_vector,iorb)
          call destroy_vector(orb_vector)

          iorb = iorb + 1

      enddo
    enddo  
    call set_ldau_projs_deg(isp)

    print ('(/,a)'), "LDAUprojgen end"

  end subroutine ldauproj_gen
  !---------------------------------------------------------------

  function number_of_ldau_proj(isp) result(nldaupj)
    integer, intent(in) :: isp
    integer             :: nldaupj

    integer             :: nsm,l

    type(l_ldaushell_t),    pointer :: lshell
    type(ldaushell_t),     pointer :: shell
    type(basis_def_t), pointer :: basp     !Parameters corresponding to this basis set.

    basp => basis_parameters(isp)

    !Find total number of orbs, including degeneracy
    nldaupj=0
    do l=0, basp%lmxldaupj
       lshell => basp%l_ldaushell(l)
       do nsm=1,lshell%nn
          shell => lshell%ldaushell(nsm)
          nldaupj = nldaupj+1
       end do
    enddo
  end function number_of_ldau_proj

   !---------------------------------------------------


  !-------------------------------------------------------------------------

  subroutine pseudopotential_setup(shell,vps)
    type(ldaushell_t), intent(inout)         :: shell
    type(pseudopotential_new_t),intent(in)   :: vps    !Pseudo information

    type(rad_func_t) :: ve_val, ve_val_scaled

    ve_val = get_ve_val(vps)
    ve_val_scaled = get_ve_val_scaled(vps)

    call rad_copy(ve_val,shell%ve)
    call rad_copy(ve_val_scaled,shell%ve_PAO)

  end subroutine pseudopotential_setup

  !-------------------------------------------------------------------------
 
  subroutine auto_rc(shell,vps,l,nsm)
    type(ldaushell_t) :: shell
    Type(pseudopotential_new_t) :: vps
    integer, intent(in) :: nsm !which semicore shell we're generating
    integer, intent(in) :: l

    integer  :: nnodes, nprin
    real(dp) :: eigen, eshift, el,rc

    type(rad_func_t) :: rad_orb,ve_val
    !Automatic determination of the cut off radius for the PAOs
    !***  READING THE ENERGY-SHIFT TO DEFINE THE CUT-OFF RADIUS OF ORBITALS***
    eshift=fdf_physical('LDAU.EnergyShift',eshift_default,'Ry')

    nnodes=nsm
    nprin=l+nsm

    ve_val = get_ve_val(vps)
    rc = rad_cutoff(ve_val)
    rad_orb = rad_schro(shell%pseudo,ve_val,rc,l,nnodes,nprin,shell%z_valence,eigen)

    !Rc given by eshift   
    if(eigen.gt.0.0d0) then 
       write(6,'(/2A,I2,A)') 'ERROR ', &
            'LDAUproj with angular momentum L=',l,' not bound in the atom'
       write(6,'(2A)') 'ERROR a cut off radius must be explicitely given' 
       call die
    endif

    if(abs(eshift).gt.1.0d-5) then
       el=eigen+eshift
       shell%rc = rad_rc_vs_e(rad_orb,shell%pseudo,ve_val,l,el,nnodes)
    else
       shell%rc = rad_default_length(rad_orb) - 0.2_dp
       !rad_get_r_from_ir(shell%pseudo, rad_default_length(shell%pseudo)-2)
    endif
    write(6,'(/,A,/,A,f10.6,A)') &
         'LDAUproj cut-off radius determinated from an', &
         'energy shift=',eshift,' Ry'
    call rad_dealloc(rad_orb)
    call rad_dealloc(ve_val)
  end subroutine auto_rc

  !--------------------------------------------------------------------

  subroutine generate_first_zeta(shell,label) 
    type(ldaushell_t),intent(inout)    :: shell
    character(len=*), intent(in)   :: label

    !Internal vars
    integer                        :: nprin,l
    real(dp)                       :: eorb !, dnrm
    real(dp)                       :: rc
    real(dp)                       :: rcsan,vcte,rinn
    type(rad_func_t)               :: vtot


    l=shell%l
    write(6,'(/A,I1,A)') 'LDAUproj corresponding to state ',shell%n, sym(l)

    rc=shell%rc/shell%lambda    
    if (rc < 0.0_dp)  call die("rc < 0 for LDAU proj")

    nprin=l+shell%i_sm

    !SOFT-CONFINEMENT            
    !Calculate the soft confinement potential for this shell  
    rcsan = shell%rc+1.0E-4
    vcte = shell%vcte
   
    if (shell%rinn < 0.0_dp) then
       rinn = -shell%rinn*rcsan
    else
       rinn = shell%rinn
    endif
      write(6,*) 'rinn, rcsan ', rinn, rcsan
    vtot=rad_sum_function(shell%ve_pao,soft_confinement,rinn,rcsan)


    shell%rldaupj = rad_schro(shell%pseudo,vtot,rc,l,shell%i_sm,nprin,shell%z_valence,eorb)
    
    !The orb is multiplied by r and the spherical harmonics routines
    !expect the orb divided by r**l, so we divide it by r**(l+1)
    shell%ldaupj = rad_divide_by_r_l_1(shell%rldaupj,l,shell%lambda)

    ! Normalize
    l=shell%l
    call rad_normalize_r_l_1(shell%ldaupj,l)
    shell%rc = rad_cutoff(shell%ldaupj)

    call rad_dealloc(vtot)

  end subroutine generate_first_zeta


  !--------------------------------------------------------------------

  subroutine generate_long_first_zeta(shell,label) 
    type(ldaushell_t),intent(inout)    :: shell
    character(len=*), intent(in)   :: label

    !Internal vars
    integer                        :: nprin,l
    real(dp)                       :: eorb !, dnrm
    real(dp)                       :: rc
    real(dp)                       :: rcsan,vcte,rinn
! Arbitrary long localization radious for these orbitals
    real(dp), parameter            :: Rmax=60.0_dp
    type(rad_func_t)               :: vtot


    l=shell%l
    write(6,'(/A,I1,A)') 'LDAUproj corresponding to state ',shell%n, sym(l)
    
    rc=Rmax/shell%lambda   
    nprin=l+shell%i_sm
! Since Rc is very long here we do not use soft confinement 
! potential here
    rinn=0.0d0
    rcsan=Rmax
    vtot=rad_sum_function(shell%ve_pao,soft_confinement,rinn,rcsan)
     
    shell%rldaupj = rad_schro(shell%pseudo,shell%ve_pao,rc,l,shell%i_sm,nprin,shell%z_valence,eorb)
    
    !The orb is multiplied by r and the spherical harmonics routines
    !expect the orb divided by r**l, so we divide it by r**(l+1)
    shell%ldaupj = rad_divide_by_r_l_1(shell%rldaupj,l,shell%lambda)

    ! Normalize
    l=shell%l
    call rad_normalize_r_l_1(shell%ldaupj,l)
    call rad_dealloc(vtot)

  end subroutine generate_long_first_zeta

  subroutine fermicutoff(shell)

    type(ldaushell_t),intent(inout)    :: shell

    !Internal vars
    integer :: nrc, ir, l
    real(dp) ::r, dnrm, rc, a, dnrm_rc, w
    real(dp), parameter :: gexp=60.0_dp
    real(dp), parameter :: eps=1.0e-4_dp
    real(dp), parameter :: dnrm_rc_default=0.90_dp
    real(dp), parameter :: width_default=0.05_dp
    type(rad_func_t)          ::  fermi_func, multiply
    type(rad_grid_t)          ::  grid
    real(dp), allocatable :: aux(:)
     
           l=shell%l
           nrc=rad_get_length(shell%ldaupj)

           if(dabs(shell%rc).lt.eps) then 
            dnrm_rc=fdf_double("LDAU.CutoffNorm",dnrm_rc_default)
            write(6,'(a)')'LDAUproj cut-off radious determined from a'
            write(6,'(a,f12.6)')'cutoff norm parameter = ',dnrm_rc
            do ir=1,nrc
               r=rad_get_r_from_ir(shell%ldaupj,ir)
               dnrm=rad_get_norm_r_l_1_from_ir(shell%ldaupj,l,ir)
               if(dnrm.gt.dnrm_rc) exit
             enddo
             shell%rc=r
           endif
           rc=shell%rc

           if(dabs(shell%width).lt.eps) shell%width=width_default
           allocate(aux(1:nrc))
           w=shell%width

           do ir=1,nrc
               r=rad_get_r_from_ir(shell%ldaupj,ir)
               a=(r-rc)/w
               if(a.lt.-gexp) then 
                 aux(ir)=1.0_dp
               else if(a.gt.gexp) then 
                 aux(ir)=0.0_dp
               else
                  aux(ir)=1.0_dp/(1.0_dp+dexp(a))
               endif
           enddo

           grid=rad_get_grid(shell%ldaupj)
           call rad_alloc(fermi_func,aux,grid)
           multiply=rad_multiply(shell%ldaupj,fermi_func)
           call rad_copy(multiply,shell%ldaupj)
           call rad_normalize_r_l_1(shell%ldaupj,l)
           do ir=nrc,1,-1
              a=rad_get_value_from_ir(shell%ldaupj,ir)
              if(dabs(a).gt.eps) exit
           enddo
           nrc=ir+1
           deallocate(aux)
           allocate(aux(1:nrc))
           do ir=1,nrc
               aux(ir)=rad_get_value_from_ir(shell%ldaupj,ir)
           enddo 
           call rad_dealloc(shell%ldaupj)
           call rad_alloc(shell%ldaupj,aux,grid)
           call rad_normalize_r_l_1(shell%ldaupj,l)

           call rad_dealloc(multiply)
           call rad_dealloc(fermi_func)
           deallocate(aux)
     
   end  subroutine fermicutoff

   
end module ldauproj



