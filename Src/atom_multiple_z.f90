!     This file is part of the SIESTA package.
!     
!     Copyright (c) Fundacion General Universidad Autonoma de Madrid:
!     E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal
!     and J.M.Soler, 1996-2006.
!     
!     Use of this software constitutes agreement with the full conditions
!     given in the SIESTA license, as signed by all legitimate users.
!     

!Generation of multiple z orbitals, not the first one!

module atom_multiple_z
use precision
use fdf
use atom_generation_types,     only:shell_t,lshell_t
use atom_pao_util, only : normalize_orbital, calculate_energies
use radial, only: rad_func_t, rad_fit_parabola,rad_find_parabola_parameters, &
     rad_split_gauss, rad_get_r_from_value,rad_cutoff, rad_min, rad_parabolic_split, &
     rad_multiply_by_rl,rad_get, rad_schro,rad_dump_file
use sys, only: die
implicit none

public  generate_multiple_zeta
private


contains

 
  subroutine generate_multiple_zeta(lshell,izeta,nsm) 

    type(lshell_t),intent(inout)  :: lshell  
    integer, intent(in)           :: izeta
    integer, intent(in)           :: nsm

    !internal vars   
    real(dp)                       :: spln, rmatch
    integer                        :: i,l, ism
    real(dp)                       :: parab_norm

    real(dp)                       :: const1,const2, rc
    real(dp)                       :: spln_min, df

    type(shell_t),pointer          :: shell
    type(rad_func_t), pointer      :: same_n(:)
    real(dp), pointer              :: lambdas(:)
    real(dp)                       :: eorb

    logical :: new_split_code, fix_split_table, split_tail_norm,keep_findp_bug
    integer :: nprin, nnodes

    shell=>lshell%shell(nsm)

    !Define kind of splitting.
    split_tail_norm=fdf_boolean('PAO.SplitTailNorm',.false.)
    fix_split_table=fdf_boolean('PAO.FixSplitTable',.false.)

    if (split_tail_norm .or. fix_split_table) then
       new_split_code = .true.
    else
       new_split_code=fdf_boolean('PAO.NewSplitCode',.false.)
    endif

    if (shell%rc(izeta) < 0.0_dp) then       
       shell%rc(izeta) = shell%rc(1) * (-shell%rc(izeta))       
    endif
    rc = shell%rc(izeta)

    l=shell%l
    if(shell%rc(izeta).gt.shell%rc(1)) then
       write(6,'(/,A)') 'SPLIT: ERROR: SPLIT OPTION FOR BASIS SET '
       write(6,'(A)')   'SPLIT: ERROR: Rc FOR DOUBLE-Z, TRIPLE-Z,... SHOULD BE SMALLER '
       write(6,'(A)') 'SPLIT: ERROR:  THAN THAT OF THE INITIAL PAO!'
       call die
    endif

    !Cut-off radius for double-Z, triple-Z,..., if it is set to
    !zero in the input then it is calculated from the splitnorm
    !parameter 

    if(shell%rc(izeta).gt.1.0d-5) then
       rmatch = rc
       !Find parabola's params
       call rad_fit_parabola(shell%rphi(1),rc, l, const1, const2, parab_norm)

       call rad_get(shell%split_table,rc,spln,df)

       do i=1,izeta-1
          if(Abs(shell%rc(izeta)-shell%rc(i)).lt.1.0d-5)then 
             write(6,'(/,A,I2,A,I2,A,I2)') 'SPLIT: WARNING: Split-orbital', &
                  'with zeta=',izeta,' and zeta=',i,' are identical for l=',l
          endif
       enddo

    else !Automatic rc             

       spln=shell%split_norm

       if(izeta.gt.2) then
          spln=spln/(2.0_dp*(izeta-2) )
       endif

       
       if (new_split_code) then
          call find_split_location(shell%rphi(1),shell%split_table,l,spln,rmatch,const1,const2)
       else
          spln_min = rad_min(shell%split_table)
          if (spln < spln_min) then
             write(6,"(a,f8.5,a,f8.5)") &
                  "WARNING: Minimum split_norm parameter: ", &
                  spln_min, ". Will not be able to generate " &
                  // "orbital with split_norm = ", spln 
             !call die("See manual for new split options")
          endif          
          keep_findp_bug = fdf_boolean("PAO.Keep.Findp.Bug", .false.)
          call rad_find_parabola_parameters(shell%rphi(1),l,spln,rmatch,const1,const2,keep_findp_bug)
       endif

       !Cut-off radius for the split orbital with a desired norm
       !nrc=nsp
       !rc = logGrid%b*(exp(logGrid%a*(nsp-1))-1.0d0)*shell%lambda(izeta)
       rc = rmatch
       shell%rc(izeta) = rmatch

       do i=1,izeta-1
          if(abs(shell%rc(izeta)-shell%rc(i)).lt.1.0d-5) then             
             write(6,'(/,A,I2,A,I2,A,I2)') &
                  'SPLIT: WARNING: Split-orbital with', &
                  'zeta=',izeta,' and zeta=',i,' are identicals for l=',l
          endif
       enddo

    endif

    shell%spln(izeta) = spln

    allocate(same_n(1:nsm),lambdas(1:nsm))

    do ism=1,nsm
       same_n(ism)=lshell%shell(ism)%orb(1)
       lambdas(ism)=lshell%shell(ism)%lambda(1)
    end do
    !Parabolic split
    if (shell%multiple_z_kind == 'split') then
       shell%orb(izeta)=rad_parabolic_split(shell%rphi(1),rmatch,const1,const2,l,nsm,same_n,lambdas)
    elseif(shell%multiple_z_kind == 'splitgauss') then
       shell%orb(izeta)=rad_split_gauss(shell%rphi(1),rmatch,shell%lambda(izeta),l)
    elseif(shell%multiple_z_kind == 'nodes') then
       nprin = l+nsm+izeta
       nnodes = izeta+nsm
       shell%orb(izeta)=rad_schro(shell%pseudo,shell%ve_pao,rc,l,nnodes,nprin,shell%z_valence,eorb)
    elseif(shell%multiple_z_kind == 'nonodes') then
       nprin = l+nsm
       nnodes = nsm
       shell%orb(izeta)=rad_schro(shell%pseudo,shell%ve_pao,rc,l,nnodes,nprin,shell%z_valence,eorb)
    else
       call die("multiple_z: unknown generation scheme!")
    endif

    call normalize_orbital(shell,izeta)
    shell%rc(izeta)=rad_cutoff(shell%orb(izeta))
    shell%rphi(izeta)=rad_multiply_by_rl(shell%orb(izeta),l+1)
    call calculate_energies(shell,izeta)    
    
    deallocate(same_n,lambdas)
  end subroutine generate_multiple_zeta

  !--------------------------------------------------------------------------

  subroutine find_split_location(rphi, split_table, l,spln,rmatch,cons1,cons2)

    !
    ! Searches for the right point to fit, and
    ! fits a parabola.

    type(rad_func_t), intent(in) :: rphi, split_table

    integer, intent(in)   :: l
    real(dp), intent(in)  :: spln
    !real(dp), intent(in)  :: split_table(*)  ! Precomputed table

    real(dp), intent(out) :: cons1, cons2, rmatch
    real(dp) :: parab_norm

    rmatch=rad_get_r_from_value(split_table,spln)

    if (rmatch == rad_cutoff(split_table)) then
       write(6,"(a,f10.6)") "Split-norm parameter is too small, " &
            // "(degenerate 2nd zeta): ", spln
       call die()
    endif
    if (rmatch <= 0.1_dp) then
       call die("Cannot find split_valence match point")
    endif

    
    call rad_fit_parabola(rphi,rmatch, l, cons1,cons2,parab_norm)

  end subroutine find_split_location

  !--------------------------------------------------------------

end module atom_multiple_z
