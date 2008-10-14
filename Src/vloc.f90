!     This file is part of the SIESTA package.
!     
!     Copyright (c) Fundacion General Universidad Autonoma de Madrid:
!     E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal
!     and J.M.Soler, 1996-2006.
!     
!     Use of this software constitutes agreement with the full conditions
!     given in the SIESTA license, as signed by all legitimate users.
!     

module vloc
  use precision
  use atom_types, only:set_pseudo_local_charge, set_reduced_vlocal
  use radial
  use atom_generation_types, only: basis_parameters
  use pseudopotential, only : pseudopotential_t

  implicit none

  private
  public gen_vlocal 

contains

  subroutine gen_vlocal(isp)
    !Generate the local part of the pseudopotential for a given species.

    integer, intent(in)     ::  isp !The species we are interested in

    !Internal vars
    type(pseudopotential_t),pointer  ::  vps

    ! Rgauss is the maximum cut-off radius used in the pseudopotential
    ! generation, Rgauss is determined by comparison of the pseudopot.
    ! Corresponding to different l ( this is not possible is we have
    ! just one pseudopotential)
    ! Rgauss2 is the radius where the pseudopotentials reach  the
    ! asymptotic behaviour 2*Zval/r.
    ! For just one pseudopotential Rgauss is taken equal to Rgauss2
    !

    type(Rad_Func_t) :: chlocal, reduced_vlocal,rad_tmp

    integer              :: lmxkb,i
    real(dp)             :: zval,r_match, r_2zval, r_2zval_match,rchloc, match_value, r_i

    vps => basis_parameters(isp)%pseudopotential

    zval =  vps%zval

    lmxkb = basis_parameters(isp)%lmxkb

    r_match = 0.0_dp
    do i=0,lmxkb-1
       r_i     = rad_matching_radius(vps%vdown(i),vps%vdown(lmxkb))
       r_match = max(r_i,r_match)
    enddo

    match_value = 2.0_dp*zval

    r_2zval_match = 0.0_dp

    do i=0, lmxkb
       r_2zval = rad_get_radius_from_value(vps%vdown(i),match_value)
       print '(a,i1,a,f8.4)','V l=',i,' =-2*Zval/r beyond r=', r_2zval
       r_2zval_match = max(r_2zval,r_2zval_match)
    enddo
 
    write(6,'(a,f8.4)') 'All V_l potentials equal beyond r=', r_match
    write(6,'(a)') "This should be close to max(r_c) in ps generation"
    write(6,'(a,f8.4)') 'All pots = -2*Zval/r beyond r=', r_2zval_match



    ! Calculate local pseudopotential*
    !      
    if (r_2zval.gt.1.30_dp*r_match) then

       write(6,'(a)') "Using large-core scheme for Vlocal"

       ! In this case the atom core is so big that we do not have an asymptotic
       ! of 2*Zval/r until Rgauss2 > Rc . To retain the same asymptotic
       ! behaviour as in the pseudopotentials we modified the definition
       ! of the local potential, making it join the Vps's smoothly at rgauss
       !
       write(6,'(/,a,f10.5)') 'ATOM: Estimated core radius ',r_2zval
       if (basis_parameters(isp)%pseudopotential%nicore.eq.'nc ')&
            write(6,'(/,2a)') 'ATOM: To include non-local',&
            ' core corrections could be a good idea'

       !As all potentials are equal beyond rgauss, we can just use the
       ! s-potential here.

       call rad_smooth_large(vps%vdown(0),zval,r_match,chlocal,vps%vlocal)
       
       !!!!call vlocal2(zval,vps%nrval,vps%logGrid%a,vps%logGrid%r,vps%logGrid%drdi,vps%s%f,vps%vdown(0)%f(:),&
       !     nrgauss,vps%vlocal%f,nchloc,chlocal%f)
       !
    else

       !
       ! In this case the pseudopotential reach to an asymptotic 
       ! behaviour 2*Zval/r for a radius approximately equal to Rc.
       ! We build a generalized-gaussian
       ! "local charge density" and set Vlocal as the potential generated by
       ! it. Note that chlocal is negative.

        call rad_smooth(vps%vdown(0),zval,r_match,chlocal,vps%vlocal)

    endif


    ! Local-pseudopotential charge
    rchloc = rad_cutoff(chlocal)
   
    write(6,'(2a,f10.5)') 'ATOM: Maximum radius for' ,&
         '  4*pi*r*r*local-pseudopot. charge ',Rchloc

    !     Fill species(isp)%reduced_vlocal and species(isp)%chlocal
    call set_pseudo_local_charge(isp,chlocal)
    call rad_dealloc(chlocal)

    rad_tmp = rad_multiply_by_rl(vps%vlocal,1)
    reduced_vlocal = rad_sum_value(rad_tmp,2.0_dp*zval)
    call rad_set_origin(reduced_vlocal,2.0_dp*zval)
    !call rad_dump_file(reduced_vlocal,"reduced_vlocal.dat")

    write(6,'(2a,f10.5)') 'atom: Maximum radius for' , &
         ' r*vlocal+2*Zval: ', rad_cutoff(reduced_vlocal) 
  
    call set_reduced_vlocal(isp,reduced_vlocal)
    call rad_dealloc(reduced_vlocal)
    call rad_dealloc(rad_tmp)
  end subroutine gen_vlocal

end module vloc
