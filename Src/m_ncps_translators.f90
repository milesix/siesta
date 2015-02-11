!> An example of the use of the PSML processing library
!!
!! Siesta uses internally the old "Froyen" ps structure. This module provides
!! the functionality to fill that structure using calls to the PSML 
!! processing library.
!!
!> @author Alberto Garcia
!
module m_ncps_translators

  public :: ncps_xml2froyen_new
  public :: ncps_xml2froyen
  interface ncps_xml2froyen
     module procedure ncps_xml2froyen_new
  end interface

CONTAINS

  subroutine ncps_xml2froyen_new( ps, p, new_grid, a, b, rmax )

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
        use m_ncps_froyen_ps_t,  only: froyen_ps_t
        use SiestaXC, only: xc_id_t, get_xc_id_from_libxc
        use SiestaXC, only: xc_id_to_string

        implicit none 

        integer, parameter  :: dp = selected_real_kind(14)


        type(ps_t), intent(in)               :: ps
        type(froyen_ps_t), intent(inout)     :: p
        logical, intent(in), optional        :: new_grid 
        real(dp), intent(in), optional       :: a
        real(dp), intent(in), optional       :: b
        real(dp), intent(in), optional       :: rmax


        logical          :: want_new_grid
        integer          :: position, i, il, ir, l, n, nval_shells
        character(len=1) :: ispp, lshell
        logical          :: polarized
        real(dp)         :: zeld, zelu, occupation
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
        type(ps_annotation_t)                :: grid_annotation
        character(len=40)                    :: strvalue
        logical                              :: log_grid_in_file
        

        p%name = ps_AtomicSymbol(ps)
        p%zval         = ps_ZPseudo(ps)
        znuc           = ps_AtomicNumber(ps)
        ! This needs to be generalized
        p%gen_zval     = ps_GenerationZval(ps)
        
!
!       Partial support for libxc functionals
!       (no single-functional cases, no 'cocktails')
!
        n_xcfuncs = ps_NLibXCFunctionals(ps)
        if (n_xcfuncs == 2) then
           do i = 1, n_xcfuncs
              libxc_ids(i) = ps_LibxcId(ps,i)
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
           call get_annotation_value(ps_XCAnnotation(ps),  &
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
!
        if (ps_IsRelativistic(ps)) then
           p%irel    = 'rel'
           ispp      = 'r'
           polarized = .false.
        else
           if (ps_IsSpinPolarized(ps)) then
              p%irel    = 'isp'
              ispp      = 's'
              polarized = .true.
           else
              p%irel    = 'nrl'
              ispp      = ' '
              polarized = .false.
           end if
        endif

        if (ps_HasCoreCorrections(ps)) then
            p%nicore = 'pcec'
         else
            p%nicore = 'nc'
         endif
!
!       Grid handling. We do not assume that the file
!       uses a logarithmic grid.

           ! We want to check for grid annotations, in case
           ! the grid is already of the "atom" type

         grid_annotation = ps_GridAnnotation(ps)
         call get_annotation_value(grid_annotation,  &
              "type",strvalue,status)

         log_grid_in_file = .false.
         if (status == 0) then
            log_grid_in_file = (trim(strvalue) == "log-atom")
         endif

        want_new_grid = .false.
        if (present(new_grid)) then
           want_new_grid = new_grid
        endif

        if (want_new_grid) then
           if (.not. present(a)) call die("new grid: a not present")
           if (.not. present(b)) call die("new grid: b not present")
           if (.not. present(rmax)) call die("new grid: rmax not present")
           rmax_grid = rmax
           p%a = a
           p%b = b
           p%nr = nint(log(rmax_grid/b+1.0d0)/a)
           p%nrval = p%nr + 1  ! Count also r=0

        else if (log_grid_in_file) then
           
           print *, "Using ATOM log grid already in PSML file ..."

           call get_annotation_value(grid_annotation,  &
                                     "nrval",strvalue,status)
           if (status /= 0) call die("Cannot read nrval")
           read(strvalue,*) p%nrval

              ! Note that a and b are interchanged in Siesta!
           call get_annotation_value(grid_annotation,  &
                                      "scale",strvalue,status)
           if (status /= 0) call die("Cannot read log grid scale")
           read(strvalue,*) p%b
           call get_annotation_value(grid_annotation,  &
                                      "step",strvalue,status)
           if (status /= 0) call die("Cannot read log grid step")
           read(strvalue,*) p%a

           p%nr = p%nrval - 1  ! For backwards compatibility

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

        p%method(1)    = ps_Creator(ps)
        p%method(2)    = ps_Date(ps)
        method_string = ps_PseudoFlavor(ps)
        read(method_string,'(4a10)') (p%method(i),i=3,6) 

        p%npotd        = ps_NPotentials(ps,set="major")
        p%npotu        = ps_NPotentials(ps,set="minor")

! Allocate the radial variables and semilocal potentials


        allocate(p%r(1:p%nrval))
        allocate(p%chcore(1:p%nrval))
        allocate(p%chval(1:p%nrval))

        if (p%npotd.gt.0) then
           allocate(p%vdown(1:p%npotd,1:p%nrval))
           allocate(p%ldown(1:p%npotd))
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

        if (ps_HasCoreCorrections(ps)) then
           do ir = 2, p%nrval
              p%chcore(ir) = ps_EvaluateCoreCharge(ps,p%r(ir))
              p%chcore(ir) = p%chcore(ir) * (p%r(ir))**2
           enddo
        else
           p%chcore(2:p%nrval) = 0.0_dp
        endif

        do ir = 2, p%nrval
           p%chval(ir) = ps_EvaluateValenceCharge(ps,p%r(ir))
           p%chval(ir) = p%chval(ir) * (p%r(ir))**2
        enddo

        ! This is no longer necessary
        ! We can directly evaluate the PSML radfuncs at r=0
        ! without an explicit extrapolation

        r2=p%r(2)/(p%r(3)-p%r(2))
        p%chcore(1) = p%chcore(2) - r2*(p%chcore(3)-p%chcore(2))
        p%chval(1) = p%chval(2) - r2*(p%chval(3)-p%chval(2))
! ---


        do il = 1, p%npotd
          p%ldown(il) = ps_PotentialL(ps,il)
          do ir = 2, p%nrval
             p%vdown(il,ir) = p%r(ir) * &
                           ps_EvaluatePotential(ps,il,p%r(ir))
             p%vdown(il,ir) = p%vdown(il,ir) * 2.0_dp   ! rydberg

          enddo
          p%vdown(il,1) = p%vdown(il,2) - r2*(p%vdown(il,3)-p%vdown(il,2))
        enddo

        do il = 1, p%npotu
           p%lup(il) = ps_PotentialL(ps,il,set="minor")
           do ir = 2, p%nrval
              p%vup(il,ir) = p%r(ir) * &
                          ps_EvaluatePotential(ps,il,p%r(ir),set="minor")
              p%vup(il,ir) = p%vup(il,ir) * 2.0_dp   ! rydberg
           enddo
          p%vup(il,1) = p%vup(il,2) - r2*(p%vup(il,3)-p%vup(il,2))
        enddo

        ! Encode generation configuration and cutoffs

        nval_shells = ps_NValenceShells(ps)
        p%text = ' '
        position = 1
        do il = 1, p%npotd
           n = ps_PotentialN(ps,il)
           l = ps_PotentialL(ps,il)
           rc = ps_PotentialRc(ps,il)

           if ( .not. polarized) then
              occupation = 0.0_dp
              do i = 1, nval_shells
                 if (ps_ValenceShellL(ps,i) == l) then
                    occupation = ps_ValenceShellOccupation(ps,i)
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
                 if (ps_ValenceShellL(ps,i) == l) then
                    zeld = ps_ValenceShellOccupation(ps,i,channel="d")
                    zelu = ps_ValenceShellOccupation(ps,i,channel="u")
                    exit
                 endif
              enddo
              write(p%text(position:),9090)    &
                   n, sym(l), zeld,zelu,ispp, rc
 9090         format(i1,a1,f4.2,',',f4.2,a1,f4.2,'/')
              position = position + 17
           end if
        enddo

      end subroutine ncps_xml2froyen_new
    end module m_ncps_translators
