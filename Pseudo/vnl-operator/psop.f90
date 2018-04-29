      program psop

! Stand alone program to:
! 1. Read the pseudopotential files
! 2. Generate the local part and the Kleynman-Bylander projector
!
! Written by Javier Junquera using code from atom.F in Siesta.
! Re-written by A. Garcia to use the PSML library
! Re-written again by A. Garcia to use the standalone psop lib.
!
! Input required:
!
!     The name of the file where the pseudopotential is stored 
!     The angular momentum cutoff for KB non-local pseudopotential
!     (here lmxkb) is fixed by the maximum l in the pseudo file.

!     The number of KB projectors per angular momentum shell is
!     based on whether there are semicore states or not. There is
!     currently no support for extra "precision" projectors.

!     The reference energy for the calculation of the KB projectors
!     is fixed to the Siesta defaults (i.e.: eigenvalue of a bound
!     state)
!
!     Implementation note:
!       
!     Even though this program is able to deal with .psml files, the
!     information in them is first converted to the old "Froyen" 
!     form, which is able to offer the lowest common functionality.
!     Eventually, support for .psf and .vps files might be discontinued.

      use m_ncps, only: pseudopotential_t => froyen_ps_t, pseudo_read
      use m_ncps, only: ncps_psml2froyen
      use m_psml, psml_t => ps_t
      use m_uuid

      use m_psop, only: kbgen, compute_vlocal_chlocal
      use m_psop, only: nrmax, nkbmx

      use m_semicore_info_froyen, only: get_n_semicore_shells
      use m_libxc_sxc_translation, only: xc_id_t, get_xc_id_from_atom_id
      use GridXC, only: atomxc => gridxc_atomxc
      use GridXC, only: setXC => gridxc_setXC
      

      use psop_options
      use m_getopts

      use xmlf90_wxml, only: str
      use m_kb, only: kb_t, nprojs, kbprojs

      implicit none

      character(len=*), parameter :: PSML_VERSION = "1.1"
      character(len=*), parameter :: PSML_CREATOR = "psop-1.1"

      
      integer, parameter :: dp = selected_real_kind(10,100)

      integer, parameter :: POLY_ORDER_EXTRAPOL = 7
      real(dp), parameter :: ryd_to_hartree = 0.5_dp
      
      type(kb_t), pointer :: kb
      character(len=1), dimension(0:4) ::  lsymb = (/'s','p','d','f','g'/)

!
!     INPUT VARIABLES REQUIRED TO GENERATE THE LOCAL PART AND KB PROJECTORS
!     THEY MUST BE PROVIDED BY THE USER
!  
      character(len=200)      :: work_string
      character(len=200)      :: filename
      character(len=20)       :: file_version
      character(len=100)      :: creator

      character*10            :: functl      ! Exchange and correlation function
      character*10            :: author      ! Exchange and correlation parametr
      integer                 :: lmxkb       ! Angular momentum cutoff for 
                                             !   Kleinman-Bylander nonlocal 
                                             !   pseudopotential
      integer,  allocatable    :: nkbl(:)    ! Number of KB projectors for
                                             !   each angular momentum
      real(dp), allocatable    :: erefkb(:,:)! Reference energies (in Ry) for
                                             !   the calculation of the KB
                                             !   projectors
      logical, allocatable     :: shifted_erefkb(:,:)
      
      integer                 :: is          ! Species index               

!
!     DERIVED VARIABLE WHERE THE PSEUDO IS STORED
!
      type(pseudopotential_t) :: psr
      type(psml_t), target    :: psml_handle
      logical                 :: has_psml, proj_file_exists
      type(xc_id_t)           :: xc_id
      integer                 :: status
      type(ps_annotation_t)   :: ann, vlocal_ann

      integer  :: max_l_ps
      integer  :: nsemic(0:3)
!
!     VARIABLES RELATED WITH THE RADIAL LOGARITHMIC GRID 
!
      integer                  :: nrval      ! Number of points required to 
                                             !   store the pseudopotential and
                                             !   the wave functions in the
                                             !   logarithmic grid
                                             !   (directly read from the 
                                             !   pseudopotential file)
      real(dp)                 :: a          ! Step parameter of log. grid
                                             !   (directly read from the 
                                             !   pseudopotential file)
      real(dp)                 :: b          ! Scale parameter of log. grid
                                             !   (directly read from the 
                                             !   pseudopotential file)
      real(dp), allocatable    :: rofi(:)    ! Radial points of the 
                                             !   logarithmic grid 
                                             !   rofi(r)=b*[exp(a*(i-1)) - 1]
                                             !   (directly read from the 
                                             !   pseudopotential file)
      real(dp), allocatable    :: drdi(:)    ! Derivative of the radial 
                                             !   distance respect the mesh index
                                             !   Computed after the radial mesh
                                             !    is read
      real(dp), allocatable    :: s(:)       ! Metric array
                                             !   Computed after the radial mesh
                                             !    is read
      real(dp)                 :: rpb, ea    ! Local variables used in the 
                                             !   calculation of the log. grid

!
!     VARIABLE USED TO STORE THE SEMILOCAL COMPONENTS OF THE PSEUDOPOTENTIAL
!
      real(dp), allocatable    :: vps(:,:)   ! Semilocal components of the
                                             !   pseudopotentials 
                                             !   (directly read from the 
                                             !   pseudopotential file)
!
!     Local variables
!  
      real(dp)                 :: Zval       ! Valence charge of the atom   
                                             !   (directly read from the 
                                             !   pseudopotential file)
!
!     VARIABLES READ FROM THE PSEUDO FILE
!
      character*4              ::  nicore    ! Flag that determines whether
                                             !   non-linear core corrections
                                             !   are included
      character*3              ::  irel      ! Flag that determines whether
                                             !   the atomic calculation is 
                                             !   relativistic or not
!
!     LOCAL VARIABLES TO DEFINE THE LOCAL PART OF THE PSEUDOPOTENTIAL
!
      integer                  :: nchloc     ! Number of radial points required
                                             !   to describe chlocal
      real(dp)                 :: rchloc     ! Radius where all the semilocal
                                             !   pseudopotentials get the 
                                             !   asymptotic behaviour
      real(dp), allocatable    :: vlocal(:)  ! Local component of the pseudopot.
                                             !   Output of vlocal1 or vlocal2.
      real(dp), allocatable    :: chlocal(:) ! Charge distribution that 
                                             !   generates vlocal
                                             !   Output of vlocal1 or vlocal2.
!
!     LOCAL VARIABLES TO DEFINE THE KLEINMAN-BYLANDER PROJECTORS 
!
      real(dp), allocatable    :: rho(:)     ! Valence charge density 
                                             !   As read from the pseudo file,
                                             !   it is angularly integrated
                                             !   (i.e. multiplied by 4*pi*r^2).
      real(dp), allocatable    :: ve(:)      ! Electrostatic potential
                                             !   generated by the valence charge
                                             !   density, readed from the 
                                             !   pseudo file
      real(dp), allocatable    :: chcore(:)  ! Core charge density 
                                             !   As read from the pseudo file,
                                             !   it is angularly integrated
                                             !   (i.e. multiplied by 4*pi*r^2).
      real(dp), allocatable    :: auxrho(:)  !  Sum of the valence charge and 
                                             !   core charge (if NLCC included)
                                             !   densities to compute the 
                                             !   atomic exchange and correl.
                                             !   potential. 
                                             !   auxrho is NOT angularly integr.
                                             !   (not multiplied by 4*pi*r^2)
      integer                  :: irelt      ! Flag that determines whether the
                                             !   atomic calculation to
                                             !   generate the pseudopotential
                                             !   was relativistic (irelt = 1)
                                             !   or no relativistic (irelt = 0)
      real(dp), allocatable    :: vxc(:) ! Exchange and correlation potential
!     ps%gen_zval                            ! Valence charge of the pseudoion 
                                             !   for which the pseudo was 
                                             !   generated in the ATM code
                                             !   (it might not correspond with 
                                             !   the nominal valence charge
                                             !   of the atom if the pseudo 
                                             !   has been generated for an ionic
                                             !   configuration, for instance 
                                             !   when semicore has been
                                             !   explicitly included in the 
                                             !   valence).
                                             !   For instance, for Ba with 
                                             !   the semicore in valence,
                                             !   (atomic reference configuration
                                             !   5s2 5p6 5d0 4f0),
                                             !   chgvps = 8  (two in the 5s 
                                             !                and six in the 5p)
                                             !   zval   = 10 (two in the 5s, 
                                             !                six in the 5p, 
                                             !                and two in the 6s.
                                             !   These two last not included in
                                             !   reference atomic configuration)
      real(dp)                 :: ex         ! Total exchange energy 
      real(dp)                 :: ec         ! Total correlation energy 
      real(dp)                 :: dx         ! IntegralOf( rho * (eps_x - v_x) )
      real(dp)                 :: dc         ! IntegralOf( rho * (eps_c - v_c) )
!
!     LOCAL VARIABLES
!
      integer                  :: ndown, ir  ! Counters for the do loops
      integer                  :: l          ! Angular momentum of the channel
                                             !   that is read 
      real(dp)                 :: pi         ! pi
      real(dp)                 :: r2         ! Local variables
      integer                  :: nkb        ! Number of KB projectors

      integer                  :: set
      character(len=40)        :: method_used  !  "siesta-fit" or "-charge"

      character(len=512) :: cmd_line
      character(len=36)  :: uuid
      character(len=200) :: opt_arg, mflnm, ref_line
      character(len=10)  :: opt_name 
      integer :: nargs, iostat, n_opts, nlabels, iorb, ikb, i

      real(dp) :: rlmax
      character(len=10)     :: datestr
      integer               :: dtime(8)

      integer  :: nrl
      real(dp) :: rmax, delta
      integer, allocatable  :: isample(:)
      real(dp), allocatable :: r0(:), f0(:)
      real(dp), allocatable :: fvlocal0(:), fchlocal0(:)
      real(dp), pointer :: r(:) => null()

      real(dp) :: delta_e
      integer  :: nkb_l, nlines
      character(len=200) proj_filename, output_filename
      character(len=40) keyname, line
      

      external :: store_proj_psml, check_grid
!
!     Process options
!
      n_opts = 0

      restricted_grid = .true.
      new_kb_reference_orbitals = .false.
      debug_kb_generation = .false.
      ignore_ghosts = .false.
      kb_rmax       = 0.0_dp

      force_chlocal_method = .false.
      fit_3derivs = .false.
      use_charge_cutoff = .false.

      write_ion_plot_files = .false.
      rmax_ps_check = 0.0_dp

      proj_file_exists = .false.
      output_filename = "PSOP_PSML"

      do
         call getopts('hdgpKF:R:C:3fcvo:',opt_name,opt_arg,n_opts,iostat)
         if (iostat /= 0) exit
         select case(opt_name)
           case ('d')
              debug_kb_generation = .true.
           case ('g')
              restricted_grid = .false.
           case ('p')
              write_ion_plot_files = .true.
           case ('K')
              new_kb_reference_orbitals = .true.
           case ('F')
              proj_file_exists = .true.
              read(opt_arg,*) proj_filename
           case ('R')
              read(opt_arg,*) kb_rmax
           case ('C')
              read(opt_arg,*) rmax_ps_check
           case ('3')
              fit_3derivs = .true.
           case ('f')
              force_chlocal_method = .true.
           case ('c')
              use_charge_cutoff = .true.
           case ('o')
              read(opt_arg,*) output_filename
           case ('v')
              write(6,"(a)") "Version: " // PSML_CREATOR
              STOP
           case ('h')
            call manual()
           case ('?',':')
             write(0,*) "Invalid option: ", opt_arg(1:1)
             write(0,*) "Usage: psop [ options ] FILE"
             write(0,*) "Use -h option for manual"
             STOP
          end select
       enddo

       nargs = command_argument_count()
       nlabels = nargs - n_opts + 1
       if (nlabels /= 1)  then
          write(0,*) "Usage: psop [ options ] FILE"
          write(0,*) "Use -h option for manual"
          STOP
       endif

       call get_command_argument(n_opts,value=filename,status=iostat)
       if (iostat /= 0) then
          STOP "Cannot get filename"
       endif

      pi = dacos(-1.0_dp)
!
!     DEFINE THE INPUT VARIABLES:
!     This must be done by the user before compiling and running the code
!
      is          = 1

!
!     READ THE PSML FILE
!
      call psml_reader(filename,psml_handle)
      call ncps_psml2froyen(psml_handle,psr)

      call ps_RootAttributes_Get(psml_handle,version=file_version,uuid=uuid)
      call ps_Provenance_Get(psml_handle,level=1,creator=creator)
      write(6,"(a)") "Processing a PSML-" // trim(file_version) // " file"
      write(6,"(a)") "Last processed by: " // trim(creator)

!
!     STORE IN LOCAL VARIABLES SOME OF THE PARAMETERS READ IN THE 
!     PSEUDOPOTENTIAL, AND DEFINITION OF THE RADIAL LOGARITHMIC GRID
!
      nrval         = psr%nrval
      a             = psr%a
      b             = psr%b

      allocate( rofi(nrmax) )
      allocate( drdi(nrmax) )
      allocate( s(nrmax)    )

      rofi(1:nrval) = psr%r(1:nrval)

!     Calculate drdi and s
!     drdi is the derivative of the radial distance respect to the mesh index
!     i.e. rofi(ir)= b*[ exp( a*(i-1) ) - 1 ] and therefore
!     drdi=dr/di =a*b*exp(a*(i-1))= a*[rofi(ir)+b]

      rpb = b
      ea  = dexp(a)
      do ir = 1, nrval
        drdi(ir) = a * rpb
        s(ir)    = dsqrt( a * rpb )
        rpb      = rpb * ea
      enddo

!     Define the angular momentum cutoff for Kleinman-Bylander nonlocal pseudopo
!     In this example, we will expand up to the f-shell (l=3).
!     Therefore, we will include (lmxkb + 1) shells
!     l = 0 (s)
!     l = 1 (p)
!     l = 2 (d)
!     l = 3 (f)

!AG:  This should be configurable via the command line, and checked
!     with the maximum l in the pseudo file...

      max_l_ps = psr%ldown(psr%npotd)
      lmxkb    = max_l_ps

!     Define the number of KB projectors for each angular momentum
      allocate( nkbl(0:lmxkb) )
!     Define the reference energies (in Ry) for the calculation of the KB proj.
      allocate( erefkb(nkbmx,0:lmxkb) )
      allocate( shifted_erefkb(nkbmx,0:lmxkb) )
      call get_n_semicore_shells(psr,nsemic)
      nkbl(:) = 1
      erefkb(:,:) = huge(1.0_dp)   ! defaults 'a la Siesta'
      shifted_erefkb(:,:) = .false.

      nlines = 0
      if (proj_file_exists) then
         ! Specify number of projectors for each l, and energy shifts
         ! Maximum 2 projectors...

         nkbl(:) = 0
         open(unit=1,file=trim(proj_filename),position="rewind")
         do
            read(1,fmt=*,iostat=status) l, nkb_l, delta_e
            if (status < 0) exit
            nlines = nlines + 1
            if (l > lmxkb) call die("l> lmxkb in projs spec file")
            nkbl(l) = nkb_l
            if (nkbl(l) > 2) call die("More than 2 projectors requested")
            ! For semicore states, make sure that we get at least two
            if (nkbl(l) < (nsemic(l)+1)) then
               nkbl(l) = (nsemic(l)+1)
               print *, "Number of projectors for l=",l," increased to", nkbl(l), " (semicore states)"
            endif
            ! For non-semicore states, honor the energy shift specified
            if (nsemic(l)==0) then
               erefkb(2,l) = delta_e
               shifted_erefkb(2,l) = .true.
            endif
         enddo
         close(1)
      else
         ! Just generate multiple projectors for channels with semicore
         nkbl(0:lmxkb) =  1 + nsemic(0:lmxkb)
      endif
         
! 
!     STORE THE IONIC PSEUDOPOTENTIALS IN A LOCAL VARIABLE 
!     Only the 'down'/major component is used in this version
!     (This will change when we implement "lj" operation)
!
      allocate( vps(nrmax,0:lmxkb) )

      do ndown = 1, lmxkb + 1
         l = psr%ldown(ndown)
         if( l .ne. ndown-1) then
           write(6,'(a)') &
             'atom: Unexpected angular momentum  for pseudopotential'
           write(6,'(a)') &
             'atom: Pseudopotential should be ordered by increasing l'
         endif
         vps(1:nrval,l) = psr%vdown(ndown,1:nrval)
!       vps contains r*V...
!       Here we compute the pseudopotential part dividing by r.
         do ir = 2, nrval
           vps(ir,l) = vps(ir,l) / rofi(ir)
         enddo
         vps(1,l) = vps(2,l)     ! AG
      enddo

!     Storing locally other variables
      Zval   = psr%zval
      nicore = psr%nicore
      irel   = psr%irel

      ! This is not sufficiently general. If a libxc functional
      ! is used, psr%icorr will not record the actual id's
      ! We have to rewrite psop forsaking the Froyen data structure
      
      call get_xc_id_from_atom_id(psr%icorr,xc_id,status)
      if (status == 0) then
         functl = xc_id%siestaxc_id%family
         author = xc_id%siestaxc_id%authors
      else
         call die("**** Cannot process XC info ***")
         functl = "LDA"
         author = "PZ"
      endif
      call setxc(1,(/functl/),(/author/),(/1.0_dp/),(/1.0_dp/))

      if (rmax_ps_check == 0.0_dp) rmax_ps_check = rofi(nrval)

      allocate (vlocal(nrmax), chlocal(nrmax))

      call compute_vlocal_chlocal(rofi,nrval,drdi,s,Zval,  &
                                  lmxkb, vps,        &
                                  a, b, nicore,      &
                                  nchloc,chlocal,    &
                                  vlocal,                &
                                  force_chlocal_method,  &
                                  rmax_ps_check,         &
                                  fit_3derivs,           &
                                  use_charge_cutoff,     &
                                  method_used)

      rchloc = rofi(nchloc)

!
!     COMPUTE THE NON-LOCAL KLEINMAN-BYLANDER PROJECTORS
!
!     Allocate internal variables
!
      allocate( ve(nrmax)     )
      allocate( vxc(nrmax)    )
      allocate( rho(nrmax)    )
      allocate( auxrho(nrmax) )
      allocate( chcore(nrmax) )

!     Read the valence charge density from the pseudo file
!     and scale it if the ionic charge of the reference configuration
!     is not the same as the nominal valence charge of the atom

!     This will set rho to the actual "descreening charge" of the
!     pseudopotential generation run (total charge will be less than       
!     zval if an ionic configuration was used)
      
      do ir = 1, nrval
        rho(ir) = psr%gen_zval * psr%chval(ir)/zval
      enddo

!     Find the Hartree potential created by a radial electron density
!     using the Numerov's method to integrate the radial Poisson equation.
!     The input charge density at this point has to be angularly integrated.
      call vhrtre( rho, ve, rofi, drdi, s, nrval, a )

!     Read the core charge density from the pseudo file
      chcore(1:nrval) = psr%chcore(1:nrval)

!     Compute the exchange and correlation potential in the atom
!     Note that atomxc expects true rho(r), not 4 * pi * r^2 * rho(r)
!     We use auxrho for that
!
      do ir = 2, nrval
        r2 = rofi(ir)**2
        r2 = 4.0_dp * pi * r2
        dc = rho(ir) / r2
        if( nicore .ne. 'nc  ')  dc = dc + chcore(ir) / r2
        auxrho(ir) = dc
      enddo

      r2        = rofi(2) / (rofi(3)-rofi(2))
      auxrho(1) = auxrho(2) - ( auxrho(3) - auxrho(2) ) * r2

!     Determine whether the atomic calculation to generate the pseudopotential
!     is relativistic or not
      if (irel.eq.'rel') then
         irelt=1
      else
         irelt=0
      endif

      call atomxc( irelt, nrval, nrmax, rofi, &
                   1, auxrho, ex, ec, dx, dc, vxc )

!     Add the exchange and correlation potential to the Hartree potential
      ve(1:nrval) = ve(1:nrval) + vxc(1:nrval)
      
!
!     Redefine the array s for the Schrodinger equation integration 
!
      s(2:nrval) = drdi(2:nrval) * drdi(2:nrval)
      s(1) = s(2)

!
!     Calculation of the Kleinman-Bylander projector functions
!

!        call get_rmax(mmax,rr,nv,uua,rmax)
        ! Choose a maximum range which is large enough for Vlocal to
        ! go to the Coulomb form

        delta = 0.005d0
        rmax = 12.0d0  ! For now
        call get_sampled_grid(nrval,rofi,rmax,delta,nrl,isample,r0)
        allocate(f0(nrl))

      ! This is also not sufficiently general. The Froyen format cannot
      ! tell apart a scalar-relativistic calculation from a fully-relativistic one.
        
      if (psr%irel == "rel") then
         ! For now, psop will generate only the scalar-relativistic set 
         set = SET_SREL
      else if (psr%irel == "isp") then
         ! Use the spin-averaged pseudos
         set = SET_SPINAVE
      else
         set = SET_NONREL
      endif

      call KBgen( is, a, b, rofi, drdi, s, &
                 vps, vlocal, ve, nrval, Zval, lmxkb, &
                 nkbl, erefkb, nkb,     &
                 new_kb_reference_orbitals, &
                 restricted_grid,   &
                 debug_kb_generation, &
                 ignore_ghosts,       &
                 kb_rmax,             &
                 process_proj=store_proj_psml,&
                 shifted_erefkb=shifted_erefkb)

      !
         call reset_annotation(ann)
         call insert_annotation_pair(ann,"source-uuid",uuid,status)
         if (status /= 0) call die("Cannot insert source-uuid")
         call get_command(cmd_line)
         call insert_annotation_pair(ann,"command-line",trim(cmd_line),status)
         if (status /= 0) call die("Cannot insert options")
         if (proj_file_exists) then
            i = 0
            open(unit=1,file=trim(proj_filename),position="rewind")
            do
               read(1,fmt="(a)",iostat=status) line
               if (status < 0) exit
               i = i + 1
               write(keyname,"(a,i1)") "proj-spec-file-",i
               call insert_annotation_pair(ann,trim(keyname),trim(line),status)
               if (status /= 0) call die("Cannot insert proj-spec-file line")
            enddo
            close(1)
         endif
         !
         if (ps_HasLocalPotential(psml_handle)) then
            call ps_Delete_LocalPotential(psml_handle)
            call insert_annotation_pair(ann,"action","replaced-local-potential",status)
            if (status /= 0) call die("Cannot insert lpot record")
         else
            call insert_annotation_pair(ann,"action","inserted-local-potential",status)
            if (status /= 0) call die("Cannot insert lpot record")
         endif
         
         if (ps_HasProjectors(psml_handle)) then
            call ps_Delete_NonLocalProjectors(psml_handle)
            call insert_annotation_pair(ann,"action-cont","replaced-nonlocal-projectors",status)
            if (status /= 0) call die("Cannot insert nl record")
         else
            call insert_annotation_pair(ann,"action-cont","inserted-nonlocal-projectors",status)
            if (status /= 0) call die("Cannot insert nl record")
         endif

         call get_uuid(uuid)
         call ps_SetUUID(psml_handle,uuid)
         call ps_SetPSMLVersion(psml_handle,PSML_VERSION)
         call date_and_time(VALUES=dtime)
         write(datestr,"(i4,'-',i2.2,'-',i2.2)") dtime(1:3)
         call ps_AddProvenanceRecord(psml_handle,creator=PSML_CREATOR, &
              date=trim(datestr), annotation=ann)

         allocate(fvlocal0(nrl), fchlocal0(nrl))
         call resample(rofi,vlocal,nrval,r0,isample,fvlocal0,nrl)
         fvlocal0 = ryd_to_hartree*fvlocal0
         call resample(rofi,chlocal,nrval,r0,isample,fchlocal0,nrl)
         where (abs(fchlocal0) < 1.0e-98_dp) fchlocal0 = 0.0_dp

         ! Grid annotation
         call reset_annotation(ann)
           !   Note interchanged a, b
           !   r(i) = b*(exp(a*(i-1))-1)
         call insert_annotation_pair(ann,"type","sampled-log-atom",status)
         if (status /= 0) call die("Cannot insert grid type")
         call insert_annotation_pair(ann,"scale",str(b),status)
         if (status /= 0) call die("Cannot insert grid scale")
         call insert_annotation_pair(ann,"step",str(a),status)
         if (status /= 0) call die("Cannot insert grid step")
         call insert_annotation_pair(ann,"delta",str(delta),status)
         if (status /= 0) call die("Cannot insert delta of sampled grid")
         call insert_annotation_pair(ann,"rmax",str(rmax),status)
         if (status /= 0) call die("Cannot insert rmax of sampled grid")

         call reset_annotation(vlocal_ann)
         call insert_annotation_pair(vlocal_ann,"chlocal-cutoff",str(rchloc),status)
         if (status /= 0) call die("Cannot insert chocal-cutoff in vlocal annotation")

         call ps_AddLocalPotential(psml_handle,grid=r0(1:nrl), &
              grid_annotation=ann, &
              vlocal=fvlocal0,vlocal_type=trim(method_used), &
              chlocal = fchlocal0, annotation=vlocal_ann)

         call ps_AddNonLocalProjectors(psml_handle,grid=r0(1:nrl), &
              grid_annotation=ann, annotation=EMPTY_ANNOTATION, set=set, &
              nprojs=nprojs,kbprojs=kbprojs)

         call ps_DumpToPSMLFile(psml_handle,trim(output_filename))

      deallocate( rofi    )
      deallocate( drdi    )
      deallocate( s       )
      deallocate( vps     )
      deallocate( vlocal  )
      deallocate( chlocal )
      deallocate( ve      )
      deallocate( vxc     )
      deallocate( rho     )
      deallocate( auxrho  )
      deallocate( erefkb  )
      deallocate( nkbl    )
      
      deallocate(r0,f0,isample)

CONTAINS

  subroutine ps_AddLocalPotential(ps,grid,grid_annotation,vlocal,vlocal_type,chlocal,annotation)
    type(psml_t), intent(inout) :: ps
    real(dp), intent(in) :: grid(:)
    type(ps_annotation_t), intent(in)   :: grid_annotation
    real(dp), intent(in) :: vlocal(:)
    real(dp), intent(in) :: chlocal(:)
    character(len=*), intent(in) :: vlocal_type
    type(ps_annotation_t), intent(in)   :: annotation

    integer :: npts
    real(dp), pointer :: gdata(:) 
    type(ps_annotation_t), pointer   :: gannot

    npts = size(grid)
    call newGrid(ps%local%grid,npts)
    gdata => valGrid(ps%local%grid)
    gdata(1:npts) = grid(1:npts)
    gannot => annotationGrid(ps%local%grid)
    gannot = grid_annotation

    ps%local%vlocal_type = trim(vlocal_type)
    ps%local%annotation = annotation
    allocate(ps%local%Vlocal%data(npts))
    ps%local%Vlocal%grid = ps%local%grid
    allocate(ps%local%Chlocal%data(npts))
    ps%local%Chlocal%grid = ps%local%grid
    ps%local%Vlocal%data(1:npts) = vlocal(1:npts)
    ps%local%Chlocal%data(1:npts) = chlocal(1:npts)
  end subroutine ps_AddLocalPotential
  
  subroutine ps_AddNonLocalProjectors(ps,grid,grid_annotation,annotation,set,nprojs,kbprojs)
    type(psml_t), intent(inout) :: ps
    real(dp), intent(in) :: grid(:)
    type(ps_annotation_t), intent(in)   :: grid_annotation
    type(ps_annotation_t), intent(in)   :: annotation
    integer, intent(in)   :: set
    integer, intent(in)   :: nprojs
    type(kb_t), intent(in), target :: kbprojs(:)

    type(nonlocal_t), pointer :: nlp, qnlp
    type(nlpj_t), pointer :: nlpp, qnlpp
    type(kb_t), pointer   :: kb
    
    integer :: npts, status, i
    real(dp), pointer :: gdata(:) 
    type(ps_annotation_t), pointer   :: gannot

    ! Allocate new node and add to the end of the linked list
    allocate(nlp)

    if (associated(ps%nonlocal)) then
       qnlp => ps%nonlocal
       do while (associated(qnlp%next))
          qnlp => qnlp%next
       enddo
       qnlp%next => nlp
    else
       ps%nonlocal => nlp
    endif

    nlp%set = set
    
    npts = size(grid)
    call newGrid(nlp%grid,npts)
    gdata => valGrid(nlp%grid)
    gdata(1:npts) = grid(1:npts)
    gannot => annotationGrid(nlp%grid)
    gannot = grid_annotation

    nlp%annotation = annotation

    do i = 1, nprojs
       kb => kbprojs(i)
       allocate(nlpp)

       ! Append to end of list  !! call append(nlp%proj,nlpp)
       if (associated(nlp%proj)) then
          qnlpp => nlp%proj
          do while (associated(qnlpp%next))
             qnlpp => qnlpp%next
          enddo
          qnlpp%next => nlpp
       else
          !First link
          nlp%proj => nlpp
       endif

       nlpp%parent_group => nlp   ! current nonlocal-projectors element

       nlpp%l = lsymb(kb%l)
       nlpp%j = kb%j
       nlpp%seq = kb%seq
       nlpp%ekb = ryd_to_hartree*kb%ekb
       nlpp%eref = ryd_to_hartree*kb%erefkb
       nlpp%type  = "KB"

       ! radfunc 'proj'
       
       allocate(nlpp%proj%data(npts))
       nlpp%proj%grid = nlp%grid
       call resample(kb%r,kb%proj,kb%nrval,r0,isample,f0,nrl)
       if (kb%l /= 0)  f0(1) = 0.0_dp
       nlpp%proj%data = f0

    enddo

  end subroutine ps_AddNonLocalProjectors
  
subroutine get_label(str,label,stat)
 character(len=*), intent(in)   :: str
 character(len=*), intent(out)  :: label
 integer, intent(out)           :: stat

 integer n, i, lo, hi

 n = len_trim(str)
 stat = -1
 lo = 1
 hi = -1
 do i = n, 1, -1
!    print *, "i, c:", i, "|",str(i:i),"|"
    if (str(i:i) == ".") then
       hi = i-1
!       print *, "hi set to: ", hi
       exit
    endif
 enddo

 if (hi>=lo) then
    stat = 0
    label=str(lo:hi)
 endif

end subroutine get_label

  subroutine get_sampled_grid(mmax,rr,rmax,delta,nrl,isample,r0)
   integer, intent(in)   :: mmax
   real(dp), intent(in)  :: rr(:)
   real(dp), intent(in)  :: rmax
   real(dp), intent(in)  :: delta
   integer, intent(out)  :: nrl
   integer, allocatable, intent(out) :: isample(:)
   real(dp), allocatable, intent(out) :: r0(:)

   integer  :: is, j
   real(dp) :: rs

   ! First scan to get size of sampled grid
   is = 1
   rs = 0.0_dp
   do j = 1, mmax
      if (rr(j) > rmax) exit
      if ((rr(j)-rs) < delta) cycle
      is = is + 1
      rs = rr(j)
   enddo
   
   nrl = is
   allocate(isample(nrl),r0(nrl))

   is = 1
   r0(is) = 0.0_dp
   isample(is) = 0
   do j = 1, mmax
      if (rr(j) > rmax) exit
      if ((rr(j)-r0(is)) < delta) cycle
      is = is + 1
      r0(is) = rr(j)
      isample(is) = j
   enddo
 end subroutine get_sampled_grid

  subroutine resample(rr,ff,mmax,r0,isample,f0,nrl)
   integer, intent(in)   :: mmax
   real(dp), intent(in)  :: rr(:)
   real(dp), intent(in)  :: ff(:)
   integer, intent(in)   :: nrl
   real(dp), intent(in)  :: r0(:)
   integer, intent(in)   :: isample(:)
   real(dp), intent(out) :: f0(:)

   integer  :: is
   real(dp) :: val
   
   do is = 2, nrl
      f0(is) = ff(isample(is))
   enddo
   ! Choice of treatments of point at r=0
   ! Polynomial extrapolation with sampled points
   call dpnint1(POLY_ORDER_EXTRAPOL,r0(2:),f0(2:),nrl-1,0.0_dp,val,.false.)
   f0(1) = val
   ! Simply set f0(r=0) = ff(r=r1)
   !...
   ! Others
   ! ...
   
 end subroutine resample

!
! Copyright (c) 1989-2014 by D. R. Hamann, Mat-Sim Research LLC and Rutgers
! University
! 
! Modified by Alberto Garcia, March 2015
! This routine is included in this module with permission from D.R. Hamann.
!
 subroutine dpnint1(npoly, xx, yy, nn, r, val, debug)

! Modified by Alberto Garcia, March 2015 from routine
! dpnint by D.R. Hamann. 
! Changes:
!   -- A single value is returned
!   -- It can extrapolate, instead of stopping,
!      when called with an abscissa outside the
!      data range.
!   -- If the number of data points is less than
!      npoly+1, npoly is implicitly reduced, without
!      error, and without warning.
!   -- Debug interface 
!
! local polynomial interpolation of data yy on nn points xx
! giving value val on point r
! npoly sets order of polynomial
! xx must be ordered in ascending order
! output interpolated value val on point r

 implicit none

 integer, parameter :: dp=kind(1.0d0)

!Input variables
 real(dp), intent(in) :: xx(*),yy(*)
 real(dp), intent(in) :: r
 real(dp), intent(out) :: val
 integer, intent(in)   ::  nn,npoly
 logical, intent(in)   ::  debug

!Local variables
 real(dp) :: sum,term,zz
 integer ii,imin,imax,iprod,iy,istart,kk,iend

! interval halving search for xx(ii) points bracketing r

   imin = 1
   imax = nn
   do kk = 1, nn
     ii = (imin + imax) / 2
     if(r>xx(ii)) then
       imin = ii
     else
       imax = ii
     end if
     if(imax - imin .eq. 1) then
       exit
     end if
   end do


   zz=r

!   if (debug) print *, "imin, imax: ", imin, imax

   if(mod(npoly,2)==1) then
    istart=imin-npoly/2
   else if(zz-xx(imin) < xx(imax)-zz) then
     istart=imin-npoly/2
   else
     istart=imax-npoly/2
   end if

   istart = min(istart, nn - npoly)
   istart = max(istart, 1)
   iend = min(istart+npoly,nn)

 !  if (debug) print *, "istart, iend: ", istart, iend
   sum=0.0d0
   do iy=istart,iend
    if(yy(iy)==0.0d0) cycle
    term=yy(iy)
    do iprod=istart, iend
     if(iprod==iy) cycle
     term=term*(zz-xx(iprod))/(xx(iy)-xx(iprod))
    end do
    sum=sum+term
   end do
   val=sum

 end subroutine dpnint1

subroutine manual()
  write(0,*) "Usage: psop [ options ] FILE"
  write(0,*) " FILE: A PSML file"
  write(0,*) " Options: "
  write(0,*) " -d                                debug"
  write(0,*) " -p                 write_ion_plot_files"
  write(0,*) " -F  PROJ_SPEC"
  write(0,*) "     Read number of projectors and      "
  write(0,*) "     reference energies from PROJ_SPEC  "
  write(0,*) " -K  use NEW-style KB reference orbitals"
  write(0,*) " -g  use unrestricted log grid"
  write(0,*) " -R  Rmax_kb (bohr)    for KB generation"
  write(0,*) " -C rmax_ps_Check (bohr) for tail checks"
  write(0,*) " -3  fit Vlocal with continuous 3rd derivative"
  write(0,*) " -f  force 'gaussian charge' method for Vlocal"
  write(0,*) " -c  force the use of 'charge cutoff' for chlocal"
  write(0,*) " -o OUTPUT_FILENAME (default PSOP_PSML)"

end subroutine manual

end program psop