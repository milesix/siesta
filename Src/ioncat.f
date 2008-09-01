      program ioncat
!
!     Looks inside .ion files (which hold the PAO, KB, and Vna tables)
!     Generates data files for plotting, optionally with zoom feature
!     near rc.

      use precision
      use atm_types, only: species, nspecies, species,species_info_t,
     $     set_label, set_read_from_file
      use atmfuncs
      use atmfuncs_types
      use m_getopts

      use basis_io, only: read_ion_ascii

      implicit none

C Internal variables ...................................................
      integer is, i, j
      integer, parameter :: n=1000
      real(dp) :: delta, rc, rmin, rmax, range

      real(dp) :: x, val, grad

      type(species_info_t), pointer   :: spp

      character(len=200) :: opt_arg
      character(len=10)  :: opt_name 
      integer :: nargs, iostat, n_opts, nlabels, iorb, ikb
      integer :: no, nkb

      logical :: show_all    = .false.,
     $           show_shells = .false.,
     $           show_kbshells = .false.,
     $           zoom        = .false.,
     $           process_orb = .false.,
     $           process_kbp = .false.,
     $           process_vna = .false.
      character(len=1000) nameat
!
!     Process options
!
      n_opts = 0
      do
         call getopts('sijo:k:vZh',opt_name,opt_arg,n_opts,iostat)
         if (iostat /= 0) exit
         select case(opt_name)
         case ('s', '+s')
            show_all = .true.
            write(0,*) "Will show everything ", trim(opt_arg)
         case ('i', '+i')
            show_shells = .true.
         case ('j', '+j')
            show_kbshells = .true.
         case ('o', '+o')
            process_orb = .true.
            write(0,*) "Will process orb ", trim(opt_arg)
            read(opt_arg,*) iorb
         case ('k', '+k')
            process_kbp = .true.
            write(0,*) "Will process KB proj ", trim(opt_arg)
            read(opt_arg,*) ikb
         case ('v', '+v')
            process_vna = .true.
            write(0,*) "Will process Vna "
         case ('Z', '+Z')
            zoom = .true.
            write(0,*) "Will use zoom around rc"
         case ('h', '+h')
            call print_help()
         case ('?',':')
            write(0,*) "Invalid option: ", opt_arg(1:1)
            call print_help()
            STOP
         end select
      enddo

      nargs = command_argument_count()
      nlabels = nargs - n_opts + 1
      if (nlabels /= 1)  then
         call print_help()
         STOP
      endif

      call get_command_argument(n_opts,value=nameat,status=iostat)
      if (iostat /= 0) then
          STOP "Cannot get Species_Label"
      endif

      nspecies = 1
      allocate(species(1))
      spp => species(1)
      call set_label(spp,trim(nameat))
      call set_read_from_file(spp,.true.)
      call read_ion_ascii(spp)

      no = nofis(1)
      nkb = nkbfis(1)

      if (show_all) then
         write(0,*) "Atomic number, norbs, nkbs: ", izofis(1), no, nkb
         write(0,*) "Orbitals (#, l, z, m, rc):"
         do i= 1, no
            write(6,*) i,
     $           lofio(1,orb_f,i), zetafio(1,i), 
     $           mofio(1,orb_f,i), rcut(1,orb_f,i)
         enddo
         write(6,*) "KB projs (#, l, m, rc):"
         do i= 1, nkb
            write(6,*) i, lofio(1,kbpj_f,i), mofio(1,kbpj_f,i),
     $           rcut(1,kbpj_f,i)
         enddo
         write(6,*) "Vna rcut: ", rcut(1,vna_f,0)
      endif
!
      if (show_shells) then
         i = 0
         do 
            if (i == no) exit
            i = i + 1
            write(6,fmt="(i3)",advance="no") i
            i = i +  2* lofio(1,orb_f,i)           ! Skip other m copies
         enddo
         write(6,*)
      endif
!
      if (show_kbshells) then
         i = 0
         do 
            if (i == nkb) exit
            i = i + 1
            write(6,fmt="(i3)",advance="no") i
            i = i +  2* lofio(1,kbpj_f,i)           ! Skip other m copies
         enddo
         write(6,*)
      endif
!
      if (process_orb) then
         if (iorb > no) STOP "no such orbital"
         write(6,*) "# Orbital (#, l, z, m, rc):",
     $        lofio(1,orb_f,iorb), zetafio(1,iorb), 
     $        mofio(1,orb_f,iorb), rcut(1,orb_f,iorb)
         rc = rcut(1,orb_f,iorb)
         rmin = 0.0_dp
         rmax = 1.05_dp * rc
         if (zoom) then
            rmin = 0.95_dp * rc
            rmax = 1.05_dp * rc
         endif
         range = rmax - rmin
         delta = range/n
         do i=0,n
            x = rmin + delta*i
            call rphiatm(1,orb_f,iorb,x,val,grad)
            write(*,"(3f14.8)") x, val, grad
         enddo
      endif

      if (process_kbp) then
         if (ikb > nkb) STOP "no such KB projector"
         write(6,*) "# KB proj (#, l, m, rc):",
     $        i, lofio(1,kbpj_f,ikb), mofio(1,kbpj_f,ikb), 
     $        rcut(1,kbpj_f,ikb)
         rc = rcut(1,kbpj_f,ikb)
         rmin = 0.0_dp
         rmax = 1.05_dp * rc
         if (zoom) then
            rmin = 0.95_dp * rc
            rmax = 1.05_dp * rc
         endif
         range = rmax - rmin
         delta = range/n
         do i=0,n
            x = rmin + delta*i
            call rphiatm(1,kbpj_f,ikb,x,val,grad)
            write(*,"(3f14.8)") x, val, grad
         enddo
      endif

      if (process_vna) then
         write(6,*) "# Vna, rcut: ", rcut(1,vna_f,0)
         rc = rcut(1,vna_f,0)
         rmin = 0.0_dp
         rmax = 1.05_dp * rc
         if (zoom) then
            rmin = 0.95_dp * rc
            rmax = 1.05_dp * rc
         endif
         range = rmax - rmin
         delta = range/n
         do i=0,n
            x = rmin + delta*i
            call rphiatm(1,vna_f,i,x,val,grad)
            write(*,"(3f14.8)") x, val, grad
         enddo
      endif

      end program ioncat

      subroutine print_help()
         print *, "Usage: ioncat [options] Species_Label"
         print *, "Options:"
         print *, " "
         print *, " -s          : Show header information"
         print *, " -i          : Print indexes of unique orbitals"
         print *, " -j          : Print indexes of unique KB projectors"
         print *, " -o ORBINDEX : Generate table for orbital ORBINDEX"
         print *, " -k KBINDEX  : Generate table for KB proj KBINDEX"
         print *, " -v          : Generate table for Vna potential"
         print *, " -Z          : Zoom in near rc for table generation"
         print *, " -h          : Print this help message"
         print *, " "
      end subroutine print_help



