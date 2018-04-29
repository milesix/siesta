!> An example of the use of the PSML processing library
!!
!! Siesta uses internally the old "Froyen" ps structure. This module provides
!! the functionality to fill that structure using calls to the PSML 
!! processing library.
!!
!> @author Alberto Garcia
!
module m_ncps_translators

  public :: ncps_psml2froyen

CONTAINS

  subroutine ncps_psml2froyen( ps, p, new_grid, a, b, rmax )

! Translate the more complete xml-adapted data structure 
! into the 'Froyen' ps type used in Atom and Siesta.
! Use accessors, and assume nothing about the grid in 
! the XML file. For this,
!
! p%nr, p%a, and p%b should be set on entry if we want
! a new grid. By default, a logarithmic grid with the
! vanilla ATOM parameters is used.
!

        use m_psml
        use m_ncps_psml_plugins
        use m_ncps_froyen_ps_t,  only: froyen_ps_t
        use m_libxc_compat, only: xc_id_t, get_xc_id_from_libxc
        use m_libxc_compat, only: xc_id_to_string

        implicit none 

        integer, parameter  :: dp = selected_real_kind(14)


        type(ps_t), intent(in)               :: ps
        type(froyen_ps_t), intent(inout)     :: p
        logical, intent(in), optional        :: new_grid 
        real(dp), intent(in), optional       :: a
        real(dp), intent(in), optional       :: b
        real(dp), intent(in), optional       :: rmax


        logical          :: want_new_grid
        integer          :: position, i, il, ir, l, n, nval_shells, ii
        character(len=1) :: ispp, lshell
        logical          :: polarized
        real(dp)         :: zeld, zelu, occupation, jval
        real(dp)         :: r2, rc, rmax_grid, znuc
        character(len=64):: xc_string
        character(len=40):: method_string
        character(len=1), dimension(0:4) :: &
                         sym = (/ "s", "p", "d", "f", "g" /)
        integer          :: libxc_ids(2), n_xcfuncs
        integer          :: status

        ! These are the "current" ATOM parameters
        ! There is another set turned on by the UCB_COMPAT flag in ATOM
        real(dp), parameter :: aa_def = 6.0_dp
        real(dp), parameter :: bb_def = 80.0_dp    ! UCB_COMPAT: 40.0
        real(dp), parameter :: rmax_def = 120.0_dp ! UCB_COMPAT: 80.0

        type(xc_id_t)                        :: xc_id
        type(ps_annotation_t)                :: grid_annotation, annot
        character(len=40)                    :: strvalue
        logical                              :: log_grid_in_file

        integer, allocatable, dimension(:)   :: idxd, idxu, idxlj
        integer, allocatable, dimension(:)   :: nn, ll
        real(dp), allocatable, dimension(:)  :: rrc

        character(len=10) :: psml_version
        character(len=10) :: relativity
        logical           :: spin_polarized, core_corrections
        integer           :: l_shell
        real(dp)          :: zup_shell, zdown_shell, occ_shell
        
        integer :: iu, id, li, npotd, npotu, nscalar, nlj, lmax
        logical :: has_lj, has_sr, has_sr_so, has_nonrel
        logical :: has_up_down, has_spin_ave
        real(dp) :: v

        call ps_RootAttributes_Get(ps,version=psml_version)
        write(6,"(a)") "PSML file version: " // &
                        trim(psml_version)
        call ps_PseudoAtomSpec_Get(ps,atomic_symbol=p%name,z_pseudo=p%zval,&
             atomic_number=znuc, pseudo_flavor=method_string,&
             relativity=relativity,spin_dft=spin_polarized,&
             core_corrections=core_corrections)
        
        call ps_ValenceConfiguration_Get(ps,nshells=nval_shells, &
                                         charge=p%gen_zval)

        !
        call ps_ExchangeCorrelation_Get(ps,annotation=annot,&
                                           n_libxc_functionals=n_xcfuncs)
        
!       Partial support for libxc functionals
!       (no single-functional cases, no 'cocktails')
!
        if (n_xcfuncs == 2) then
           do i = 1, n_xcfuncs
              call ps_LibxcFunctional_Get(ps,i,code=libxc_ids(i))
!              libxc_ids(i) = ps_LibxcId(ps,i)
           enddo
           call get_xc_id_from_libxc(libxc_ids,xc_id,status)
        else
           write(6,"(a,2i4)") "Cannot handle libxc cases with nfunc/=2..."
           status = -1
        endif

        if (status == 0) then
              write(6,"(a,2i4)") "Using libxc ids: ", libxc_ids(:)
              write(6,"(a)") trim(xc_id_to_string(xc_id))
              p%icorr = xc_id%atom_id
        else
           ! Fall back to querying a possible XC annotation
           call get_annotation_value(annot,  &
                                      "atom-xc-code",p%icorr,status)
           if (status == 0) then
              write(6,"(a)") "Atom-xc-code from annotation: " // p%icorr
           else
              call die("Cannot get atom-xc code")
           endif
        endif
!
!       Note that most (all?) formats include only "scalar-relativistic" plus
!       maybe spin-orbit components, but never "polarized" pseudos.
!       This is a shortcoming of the "Froyen" format: there is no support
!       for "scalar_relativistic" calculations in the "irel" label.

        if (trim(relativity) == "dirac" ) then
           p%irel    = 'rel'
           ispp      = 'r'
           polarized = .false.
        else
           if (spin_polarized) then
              p%irel    = 'isp'
              ispp      = 's'
              polarized = .true.
           else
              p%irel    = 'nrl'
              ispp      = ' '
              polarized = .false.
           end if
        endif

        if (core_corrections) then
            p%nicore = 'pcec'
         else
            p%nicore = 'nc'
         endif
!
!       Grid handling. We do not assume that the file
!       uses a logarithmic grid.

        want_new_grid = .false.
        if (present(new_grid)) then
           want_new_grid = new_grid
        endif

        ! We want to check for grid annotations, in case
        ! the grid is already of the "atom" type
        ! This is an example of a "PSML plugin"
        call check_atom_grid(ps,log_grid_in_file,p%nrval,p%a,p%b)

        if (log_grid_in_file) then
           
           if (p%nrval > 0) then  ! Normal ATOM grid
              print *, "Using ATOM log grid in PSML file ..."
              p%nr = p%nrval - 1  ! For backwards compatibility
           else    ! Sampled ATOM grid; nrval not available               
              print *, "Using ATOM log grid {a,b} parameters in PSML file ..."
              rmax_grid = rmax_def
              p%nr = nint(log(rmax_grid/(p%b)+1.0d0)/(p%a))
              p%nrval = p%nr + 1  ! Count also r=0
           endif

        else if (want_new_grid) then
           if (.not. present(a)) call die("new grid: a not present")
           if (.not. present(b)) call die("new grid: b not present")
           if (.not. present(rmax)) call die("new grid: rmax not present")
           rmax_grid = rmax
           p%a = a
           p%b = b
           p%nr = nint(log(rmax_grid/b+1.0d0)/a)
           p%nrval = p%nr + 1  ! Count also r=0

        else 
              print *, "Using ATOM defaults for log grid..."
              ! use the ATOM defaults 
              ! Note that a and b are interchanged in Siesta!
              p%a = 1.0_dp / bb_def
              p%b = exp(-aa_def)/znuc
              rmax_grid = rmax_def
              p%nr = nint(log(rmax_grid/(p%b)+1.0d0)/(p%a))
              p%nrval = p%nr + 1  ! Count also r=0
        endif

        ! Assume first action is ps generation...
        call ps_Provenance_Get(ps,level=1,creator=p%method(1),&
                               date=p%method(2))
        read(method_string,'(4a10)') (p%method(i),i=3,6) 

        has_nonrel = .false.
        has_sr = .false.
        has_sr_so = .false.
        has_up_down = .false.
        has_spin_ave = .false.
        has_lj = .false.

        select case (trim(relativity))
        case ("dirac")

          call ps_Potential_Filter(ps,set=SET_SREL,number=nscalar,indexes=idxd)

           if (nscalar == 0) then

              ! Will get the scalar-relativistic SL potentials
              ! from the lj set

              call ps_Potential_Filter(ps,set=SET_LJ,number=nlj,indexes=idxlj)
              if (nlj == 0) call die("Cannot find srel SL potentials for dirac case")
              has_lj = .true.
              npotd = 0
              npotu = 0
              do i = 1, nlj
                 call ps_Potential_Get(ps,idxlj(i),l=l,j=jval)
                 if ( (l==0) .or. (jval>l)) then
                    npotd = npotd + 1
                 else
                    npotu = npotu + 1
                 endif
              enddo

           else  ! nscalar > 0

              ! We have a scalar-relativistic set
              npotd = nscalar
              ! Get the information on the companion spin-orbit set
              call ps_Potential_Filter(ps,set=SET_SO,number=npotu,indexes=idxu)
              has_sr_so = .true.

           endif

        case ("scalar")

           call ps_Potential_Filter(ps,set=SET_SREL,number=npotd,indexes=idxd)
           if (npotd == 0) call die("Cannot find srel SL potentials for srel case")
           npotu = 0
           has_sr = .true.

           ! We assume that srel calculations are not polarized...
        case ("no")

           if (spin_polarized) then
              call ps_Potential_Filter(ps,set=SET_UP,number=npotu,indexes=idxu)
              call ps_Potential_Filter(ps,set=SET_DOWN,number=npotd,indexes=idxd)
              if (  (npotu > 0) .and. (npotd > 0) ) then

                 ! We have spin_up and spin_down potentials
                 ! Will get the average later

                 has_up_down = .true.

              else
                 ! We must have (at least) the spin_average
                 call ps_Potential_Filter(ps,set=SET_SPINAVE,number=npotd,indexes=idxd)
                 call ps_Potential_Filter(ps,set=SET_SPINDIFF,number=npotu,indexes=idxu)
                 if (npotd == 0) call die("Cannot get spin-averaged SL potentials")
                 has_spin_ave = .true.
              endif

           else ! not polarized

              call ps_Potential_Filter(ps,set=SET_NONREL,number=npotd,indexes=idxd)
              if (npotd == 0) call die("Cannot get non-relativistic SL potentials")
              npotu = 0
              has_nonrel = .true.

           endif

        case default

           call die("Wrong relativity scheme in PSML file")

        end select

! Allocate the radial variables and semilocal potentials

        p%npotd = npotd
        p%npotu = npotu

        allocate(p%r(1:p%nrval))
        allocate(p%chcore(1:p%nrval))
        allocate(p%chval(1:p%nrval))

        if (p%npotd.gt.0) then
           allocate(p%vdown(1:p%npotd,1:p%nrval))
           allocate(p%ldown(1:p%npotd))
           allocate(nn(1:npotd),ll(1:npotd),rrc(1:npotd))
        endif

        if (p%npotu.gt.0) then
           allocate(p%vup(1:p%npotu,1:p%nrval))
           allocate(p%lup(1:p%npotu))
        endif
! ---

! Calculate the points of the logarithmic radial grid 
        ! The first point here (NOT in the classic files
        ! produced by ATOM (vps,psf) is r=0
        do ir = 1, p%nrval
           p%r(ir) = p%b * (exp(p%a*(ir-1))-1)
        enddo
! ---

! Translate the valence charge density and the pseudo-core charge density,
! and define the value at the first point of the logarithmic grid
        if (core_corrections) then
           do ir = 2, p%nrval
              p%chcore(ir) = ps_CoreCharge_Value(ps,p%r(ir))
              p%chcore(ir) = p%chcore(ir) * (p%r(ir))**2
           enddo
        else
           p%chcore(2:p%nrval) = 0.0_dp
        endif

        do ir = 2, p%nrval
           p%chval(ir) = ps_ValenceCharge_Value(ps,p%r(ir))
           p%chval(ir) = p%chval(ir) * (p%r(ir))**2
        enddo
        !
        ! Note this re-scaling to comply with the psf (Froyen)
        ! convention of writing a neutral-atom total valence charge
        !
        if (abs(p%zval-p%gen_zval) > 1.0e-3_dp) then
           print "(a)", "Rescaling valence charge in psml file to neutral-atom"
           print "(a,2f10.4)", "Zval, GenerationZval:", p%zval, p%gen_zval
           do ir = 2, p%nrval
              p%chval(ir) = (p%zval/p%gen_zval) * p%chval(ir)
           enddo
        endif
!
        ! This is no longer necessary
        ! We can directly evaluate the PSML radfuncs at r=0
        ! without an explicit extrapolation

        r2=p%r(2)/(p%r(3)-p%r(2))
        p%chcore(1) = p%chcore(2) - r2*(p%chcore(3)-p%chcore(2))
        p%chval(1) = p%chval(2) - r2*(p%chval(3)-p%chval(2))
! ---

        if ( (has_sr_so) .or. (has_spin_ave) .or. (has_nonrel) .or. &
             (has_sr)   ) then
         ! No need for any extra computations
         do il = 1, p%npotd
           call ps_Potential_Get(ps,idxd(il),&
                     l=p%ldown(il),n=nn(il),rc=rrc(il))
           ll(il) = p%ldown(il)
                     
           do ir = 2, p%nrval
             p%vdown(il,ir) = p%r(ir) * &
                           ps_Potential_Value(ps,idxd(il),p%r(ir))
             p%vdown(il,ir) = p%vdown(il,ir) * 2.0_dp   ! rydberg

           enddo
           p%vdown(il,1) = p%vdown(il,2) - r2*(p%vdown(il,3)-p%vdown(il,2))
         enddo

         do il = 1, p%npotu
           call ps_Potential_Get(ps,idxu(il),l=p%lup(il))
           do ir = 2, p%nrval
              p%vup(il,ir) = p%r(ir) * &
                           ps_Potential_Value(ps,idxu(il),p%r(ir))
              p%vup(il,ir) = p%vup(il,ir) * 2.0_dp   ! rydberg
           enddo
           p%vup(il,1) = p%vup(il,2) - r2*(p%vup(il,3)-p%vup(il,2))
         enddo

        ! if not, we need to get the right averages
        else if (has_lj) then

           id = 0
           iu = 0
           do i = 1, nlj
              call ps_Potential_Get(ps,idxlj(i),l=l,j=jval)
              if ( (l==0) .or. (jval>l)) then
                 id = id + 1
                 ! If the lj slpots are not ordered by l in the psml
                 ! file this array will not be monotonic
                 p%ldown(id) = l
                 ! get some extra info needed later
                 ll(id) = l
                 call ps_Potential_Get(ps,idxlj(i),n=nn(id),rc=rrc(id))
              else
                 iu = iu + 1   
                 p%lup(iu) = l
              endif
           enddo
           call assert((id == npotd),"Wrong check on number of down pots")
           call assert((iu == npotu),"Wrong check on number of up pots")

           print *, "l down: ", p%ldown(1:npotd)
           print *, "l up: ", p%lup(1:npotu)

           lmax = maxval(p%ldown(1:npotd))

           p%vdown(:,:) = 0.0_dp
           p%vup(:,:) = 0.0_dp

           do l = 0, lmax   ! Loop on l

              ! determine corresponding indexes in the arrays

              id = 0
              do ii = 1, npotd
                 if (p%ldown(ii) == l) then
                    id = ii
                    exit
                 endif
              enddo
              call assert((id>0),"l mismatch in down array")

              iu = 0   ! This will be zero for l=0, and skipped below
              do ii = 1, npotu
                 if (p%lup(ii) == l) then
                    iu = ii
                    exit
                 endif
              enddo

              do i = 1, nlj
                 call ps_Potential_Get(ps,idxlj(i),l=li,j=jval)
                 if (li /= l) cycle

                 ! Process the two (except for l=0) j channels for this l

                 if ((l==0) .or. jval > l) then  ! j=l+1/2 or l=0,j=0
                    !print *, "l,j+, i, id, iu ", l, jval, i, id, iu
                    do ir = 2, p%nrval
                       v = ps_Potential_Value(ps,idxlj(i),p%r(ir))
                       p%vdown(id,ir) =  p%vdown(id,ir) + (l+1)*v / dble(2*l+1)
                       if (iu>0) p%vup(iu,ir) = p%vup(iu,ir) + 2*v / dble(2*l+1)
                    enddo
                 else   ! j=l-1/2
                    !print *, "l,j=-, i, id, iu ", l, jval, i, id, iu
                    do ir = 2, p%nrval
                       v = ps_Potential_Value(ps,idxlj(i),p%r(ir))
                       p%vdown(id,ir) =  p%vdown(id,ir) + l*v/dble(2*l+1)
                       if (iu>0) p%vup(iu,ir) = p%vup(iu,ir) - 2*v/dble(2*l+1)
                    enddo
                 endif  ! j+ or j-
              enddo   ! over lj set

              ! rydberg and rV
              do ir = 2, p%nrval
                 p%vdown(id,ir) = 2.0_dp * p%r(ir)*p%vdown(id,ir) 
                 if (iu>0) p%vup(iu,ir) = 2.0_dp * p%r(ir)* p%vup(iu,ir)
              enddo
              ! extrapolate to r=0
              p%vdown(id,1) = p%vdown(id,2) - r2*(p%vdown(id,3)-p%vdown(id,2))
              if (iu>0) p%vup(iu,1) = p%vup(iu,2) - r2*(p%vup(iu,3)-p%vup(iu,2))

           enddo        ! over l values

           ! Do we need to re-sort the arrays if l is not monotonic?

        else if (has_up_down) then
           ! to be implemented
           call die("up/down case not implemented yet")
           !
        endif

        ! Encode generation configuration and cutoffs

        p%text = ' '
        position = 1
        ! We deal with "down" potentials only, as they are enough
        ! for this particular bookeeping
        do il = 1, p%npotd
           n = nn(il)  
           l = ll(il)  
           rc = rrc(il)

           if ( .not. polarized) then
              occupation = 0.0_dp
              do i = 1, nval_shells
                 call ps_ValenceShell_Get(ps,i,l=l_shell,occupation=occ_shell)
                 if (l_shell == l) then
                    occupation = occ_shell
                    exit
                 endif
              enddo
              write(p%text(position:),9070)     &
                   n, sym(l), occupation, ispp, rc
 9070         format(i1,a1,f5.2,a1,' r=',f5.2,'/')
              position = position + 17
           else
              zeld = 0.0_dp
              zelu = 0.0_dp
              do i = 1, nval_shells
                 call ps_ValenceShell_Get(ps,i,l=l_shell,occ_up=zup_shell,occ_down=zdown_shell)
                 if (l_shell == l) then
                    zeld = zdown_shell
                    zelu = zup_shell
                    exit
                 endif
              enddo
              write(p%text(position:),9090)    &
                   n, sym(l), zeld,zelu,ispp, rc
 9090         format(i1,a1,f4.2,',',f4.2,a1,f4.2,'/')
              position = position + 17
           end if
        enddo

      end subroutine ncps_psml2froyen

      subroutine assert(cond,message)
        logical, intent(in) :: cond
        character(len=*) message

        if (.not. cond) call die(message)
      end subroutine assert

    end module m_ncps_translators