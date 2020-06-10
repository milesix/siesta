!     
!     This file is part of the SIESTA package.
!     
!     Copyright (c) Fundacion General Universidad Autonoma de Madrid:
!     E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal
!     and J.M.Soler, 1996- .
!     
!     Use of this software constitutes agreement with the full conditions
!     given in the SIESTA license, as signed by all legitimate users.
!     
      module m_genq_broyden_optim
      public :: genq_broyden_optimizer
      private

      CONTAINS
!     
      subroutine genq_broyden_optimizer( na, xa, fa, cell,
     .     stress, dxmax, 
     .     tp, ftol, strtol, varcel, relaxd)

c     ******************************** INPUT *******************************
c     integer na            : Number of atoms in the simulation cell
c     real*8 fa(3,na)       : Atomic forces
c     real*8 stress(3,3)    : Stress tensor components
c     real*8 dxmax          : Maximum atomic (or lattice vector) displacement
c     real*8 tp             : Target pressure
c     real*8 ftol           : Maximum force tolerance
c     real*8 strtol         : Maximum stress tolerance
c     logical varcel        : true if variable cell optimization
c     *************************** INPUT AND OUTPUT **************************
c     real*8 xa(3,na)       : Atomic coordinates
c     input: current step; output: next step
c     real*8 cell(3,3)      : Matrix of the vectors defining the unit cell 
c     input: current step; output: next step
c     cell(ix,ja) is the ix-th component of ja-th vector
c     real*8 stress(3,3)    : Stress tensor components
c     real*8 strtol         : Maximum stress tolerance
c     ******************************** OUTPUT *******************************
c     logical relaxd        : True when converged
c     ***********************************************************************

C     
C     Modules
C     
      use precision,   only : dp
      use sys,         only : die
      use parallel,    only : Node, IONode
      use fdf,         only : fdf_get

      use units, only: kBar, Ang, eV
      use m_broyddj_nocomm, only: broyden_t
      use m_broyddj_nocomm, only: broyden_reset, broyden_step,
     $     broyden_destroy, broyden_init,
     $     broyden_is_setup

      use m_genq,      only : genq
      use m_genq,      only : geom_from_genq
      use m_genq,      only : cartF_to_genqF

      implicit none

C     Subroutine arguments:
      integer,  intent(in)   :: na
      logical,  intent(in)   :: varcel
      logical,  intent(out)  :: relaxd
      real(dp), intent(in)   :: fa(3,na), dxmax, tp, ftol
      real(dp), intent(inout):: xa(3,na), stress(3,3), strtol, cell(3,3)

C     Internal variables and arrays
      integer             :: iu, ia, i, j, n, indi,indi1,vi,k, lun
      real(dp), dimension(:), allocatable :: gxa, gfa, epsgfa
      real(dp)            :: force

C     Saved internal variables:
      integer,       save :: ndeg

      real(dp), dimension(:), allocatable    :: ftoln, dxmaxn
      real(dp), dimension(:), allocatable    :: rold, rnew, rdiff

      type(broyden_t), save  :: br
      logical, save           :: initialization_done = .false.

      real(dp), save :: jinv0
      integer, save  :: maxit
      logical, save  :: cycle_on_maxit, variable_weight
      logical, save  :: broyden_debug
      real(dp)       :: weight

      integer        :: ndegi, ndi
      real(dp)       :: xafrac(3), max_force, jinv0_def
      real(dp), save :: max_ucell, ftol_genq, dxmax_genq
      real(dp), save :: genq_tolerance
      logical, save  :: frstme = .true.
      logical        :: relaxd_genq, want_coord_change_criterion
      logical        :: use_cart_ftol
c ---------------------------------------------------------

      if ( varcel ) then
         call die("Cannot do variable_cell optim with genq yet")
      endif     

      if (.not. initialization_done) then

         maxit = fdf_get("MD.Broyden.History.Steps",5)
         cycle_on_maxit = fdf_get("MD.Broyden.Cycle.On.Maxit",.true.)
         variable_weight = fdf_get("MD.Broyden.Variable.Weight",.false.)
         broyden_debug = fdf_get("MD.Broyden.Debug",.false.)

         if (ionode) then
            write(6,*) 
            write(6,'(a,i3)')
     $           "Broyden_optim: max_history for broyden: ", maxit
            write(6,'(a,l1)')
     $           "Broyden_optim: cycle on maxit: ", cycle_on_maxit
            if (variable_weight) write(6,'(a)')
     $           "Broyden_optim: Variable weight not implemented yet"
            write(6,*) 
         endif

         call broyden_init(br,debug=broyden_debug)
         initialization_done = .true.
      endif

      if ( frstme ) then
         ! Find number of variables
         ndeg = genq%nvars
         if (Ionode) then
            write(6,'(a,i6)') "Broyden Genq optim: No of elements: ",
     $           ndeg
         endif
      endif                     ! First time

      allocate(gfa(ndeg))
      allocate(epsgfa(ndeg))
      allocate(gxa(ndeg))

      gxa(1:ndeg) = genq%q(1:ndeg)   

      call cartF_to_genqF(ndeg,genq%q,cell,na,fa,gfa,ftol,epsgfa)

      if (ionode) then
         print "(a)",
     $        "Generalized forces, estimated tolerances: (eV)"
         do i = 1, ndeg
            print "(a,2f15.4)", trim(genq%qname(i)),
     $           gfa(i)/eV, epsgfa(i)/eV
         enddo
            call io_assign(lun)
            open(lun,file="GENQ.latest",form="formatted",
     $           status="unknown",
     $           action="write", position="rewind")
            write(lun,"(i4)") ndeg
            do i = 1, ndeg
               write(lun,"(f10.6,2x,a)"), genq%q(i), trim(genq%qname(i))
            enddo
            close(lun)
      endif

      if (frstme) then

         ! NOTE: the use of cell vectors is limited to "crystallographic"
         ! applications in which the generalized coordinates are typically
         ! fractional coordinates. A better heuristic should be used.

         max_force = maxval(abs(gfa))
         max_ucell = maxval(abs(cell))

         use_cart_ftol =  fdf_get("MD.Genq.UseCartForceTol",.false.)

         if (ionode) then
            if (use_cart_ftol) then
               print "(/,a)",
     $           "Using cartesian force tolerance for all atoms (ftol)"

            else
               print "(/,a)",
     $           "Using gen force tolerance (see above)"
            endif
         endif

         want_coord_change_criterion =
     $         fdf_get("MD.Genq.UseGenqChangeTol",.false.)

         if (want_coord_change_criterion) then

            genq_tolerance = fdf_get("MD.Genq.Tolerance",0.0001_dp)

            if (ionode) then
               print "(/,a,f10.6)",
     $              "Effective tolerance for " //
     $              "Gen. Coords change: ", genq_tolerance
               print "(a)", "Use MD.Genq.Tolerance fdf var to change"
            endif
         endif

!  Heuristics to avoid moving the atoms too much
         dxmax_genq = dxmax / max_ucell

         if (ionode) then
            print "(/,a,f12.5)",
     $      "Effective max displacement for Gen. Coords: ", dxmax_genq
            print "(a,f12.5)",
     $         "(standard dxmax) / max cell size (Ang):", max_ucell/Ang
            print "(a)",
     $           "Use MD.MaxCGDispl (physical fdf var) to change dxmax"
         endif

! Note: A very small jinv0 is needed because of
! the "cell factors" involved in the calculation
! of the generalized forces. Need to explore further
! (See comment above regarding crystallographic applications)

         jinv0 = dxmax_genq / max_force 

         if (ionode) then
            print "(/,a,f12.6)",
     $     "Suggested effective Initial Inverse J (eV^-1): ", jinv0*eV
            print "(a,f12.6)",
     $     "Suggested effective Initial Inverse J (Ry^-1): ", jinv0
            print "(a)",
     $    "Use MD.Broyden.Genq.Initial.Inverse.Jacobian" //
     $           " to change (*implicitly in Ry^-1 *)"
         endif

         jinv0_def = jinv0
         jinv0 = fdf_get("MD.Broyden.Genq.Initial.Inverse.Jacobian",
     $                   jinv0_def)
         if (jinv0 /= jinv0_def) then
            print "(a,f12.6)",
     $           "Requested Initial Inverse J (Ry^-1): ", jinv0
         endif

         frstme = .false.
      endif

      if (use_cart_ftol) then
         ! This assumes that the forces on all atoms must be zero
         relaxd = all(abs(fa)<ftol)
         if (relaxd) then
            if (ionode) then
               print *,
     $              "Relaxed according to cart force tol criterion"
            endif
         endif
      else
         ! Look only at the generalized forces

         relaxd = .true.
         do i = 1, ndeg
            if (abs(gfa(i)) > abs(epsgfa(i))) relaxd = .false.
         enddo

         if (relaxd) then
            if (ionode) then
               print *,
     $           "Relaxed according to generalized force tol criterion"
            endif
         endif
      endif

         if (relaxd) then
            deallocate(gxa,gfa,epsgfa)
            RETURN                
         endif


         if (.not. broyden_is_setup(br)) then
            allocate(dxmaxn(ndeg))
            dxmaxn = dxmax_genq
            call broyden_reset(br,ndeg,maxit,cycle_on_maxit,
     $           jinv0,0.01_dp,dxmaxn)
            deallocate(dxmaxn)
         endif

         allocate (rnew(ndeg))
         weight = 1.0_dp
         call broyden_step(br,gxa,gfa,w=weight,newx=rnew)

         if (want_coord_change_criterion) then

            relaxd_genq = all(abs(gxa-rnew) < genq_tolerance) 
         

            if (relaxd_genq) then
               if (ionode) then
                  print *,
     $            "Relaxed according to coordinate change criterion"
               endif
               deallocate(gxa,gfa)
               relaxd = .true.
               RETURN
            endif
         endif

         genq%q = rnew
         deallocate(rnew)

! Transform back

         call geom_from_genq(ndeg,genq%q,xa)

         if (ionode) then
            print "(/,a)", "New values for gen_coords:"
            do i = 1, ndeg
               print "(a,f10.6)", trim(genq%qname(i)), genq%q(i)
            enddo
            call io_assign(lun)
            open(lun,file="GENQ.next",form="formatted",
     $           status="unknown",
     $           action="write", position="rewind")
            write(lun,"(i4)") ndeg
            do i = 1, ndeg
               write(lun,"(f10.6,2x,a)"), genq%q(i), trim(genq%qname(i))
            enddo
            close(lun)
            
            print "(/,a)", " New Fractional coordinates: "
            do ia = 1, na
               print "(i45,3f10.6)", ia, xa(1:3,ia)
            enddo
         endif

         do ia = 1, na
            xafrac(:) = xa(:,ia)
            xa(:,ia) = matmul(cell,xafrac(:))
         enddo

      deallocate(gxa,gfa)

      end subroutine genq_broyden_optimizer

      end module m_genq_broyden_optim
