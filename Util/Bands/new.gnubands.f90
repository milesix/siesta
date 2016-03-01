! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!

program gnubands
!
! Processes a "SystemLabel.bands" file to produce a data file suitable for
! plotting by Gnuplot
!
! This is an update of the venerable "gnubands.f" program, adding options
! to specify a subset of the bands contained in the .bands file, select
! spin, and more.  
!
! A. Garcia, May 2012, 2015
!
! See also the "fatbands" functionality provided by the program "fat" in
! directory Util/COOP and the program "eigfat2plot" in this directory.
!
! Updated to be able to extract an energy range and automatic Ef shift
! Bugfix for max_bands. It was required that one supplied the -B
! option, else it was not set.
!
! Nick Papior Andersen, April 2013, 2015

  use m_getopts
  use f2kcli

  implicit none

  integer, parameter :: dp = selected_real_kind(10,100)

  integer :: nk, nspin, nband, ik, is, ib
  
  real(dp), allocatable ::  e(:,:,:), k(:)
  real(dp) :: ef, kmin, kmax, emin, emax
  real(dp) :: dummy, delta

  integer :: spin_idx
  integer :: min_band = 1
  integer :: max_band = huge(1)
  logical :: min_band_set = .false.
  logical :: max_band_set = .false.

  integer :: bands_u = 100 ! Some compilers use 1 for special things
  integer :: n_opts, iostat, nargs, nlabels
  character(len=132) :: opt_name, opt_arg
  character(len=132) :: bandfile
  integer :: nbands
  logical :: add_new_line, gnu_ticks
  logical :: Fermi_shift, emin_set, emax_set

  ! Number of lines
  integer :: nlines
  real(dp), allocatable :: listk(:)
  character(len=8), allocatable :: labels(:)

  Fermi_shift = .false.
  emin_set = .false.
  emax_set = .false.
  spin_idx = 0
  gnu_ticks = .false.

  ! Process options
  n_opts = 0
  do
     call getopts('hb:GB:Fe:E:s:',opt_name,opt_arg,n_opts,iostat)
     if (iostat /= 0) exit
     select case(opt_name)
     case ('G')
        gnu_ticks = .true.
     case ('s')
        read(opt_arg,*) spin_idx
     case ('F') 
        Fermi_shift = .true.
     case ('e')
        emin_set = .true.
        read(opt_arg,*) emin
     case ('E')
        emax_set = .true.
        read(opt_arg,*) emax
     case ('b')
        read(opt_arg,*) min_band
        min_band_set = .true.
     case ('B')
        read(opt_arg,*) max_band
        max_band_set = .true.
     case ('h')
        call manual()
        STOP
     case ('?',':')
        write(0,*) "Invalid option: ", opt_arg(1:1)
        call manual()
        STOP
     end select
  enddo

  nargs = command_argument_count()
  nlabels = nargs - n_opts + 1
  if      (nlabels == 0) then
     ! Reads in from pipe
     bands_u = 5
  else if (nlabels == 1) then
     call get_command_argument(n_opts,value=bandfile,status=iostat)
     if (iostat /= 0) then
        STOP "Cannot get bands file"
     end if
  else if (nlabels /= 1)  then
     call manual()
     STOP
  end if

  ! If the file-unit is now std-in, then we should open the file
  if ( bands_u /= 5 ) then
     open(bands_u,file=trim(bandfile), &
          status='old',action="read",position="rewind")
  end if

  read(bands_u,*) ef
  read(bands_u,*) kmin, kmax
  if ( emin_set .and. emax_set ) then
     read(bands_u,*) dummy, dummy
  else if ( emin_set ) then
     read(bands_u,*) dummy, emax
  else if ( emax_set ) then
     read(bands_u,*) emin, dummy
  else
     read(bands_u,*) emin, emax
  end if
  read(bands_u,*) nband, nspin, nk
  if ( spin_idx /= 0 ) then
     if ( spin_idx < 1 .and. nspin < spin_idx ) then
        write(0,"(a)") " ** Selected spin does not exist..."
        stop
     end if
  end if

  if (min_band_set .and. (min_band < 1)) then
     write(0,"(a)") " ** Min_band implicitly reset to 1..."
     min_band = 1
  endif
  if (min_band_set .and. (min_band > nband)) then
     write(0,"(a,2i5)") " ** Min_band is too large  (min_band, nband):", min_band, nband
     STOP
  endif
  if ( max_band > nband ) then
     if ( max_band_set ) then
        write(0,"(a,2i5)") " ** Max_band is too large (max_band, nband):", max_band, nband
        write(0,"(a)") " ** Max_band will be effectively reset to its maximum allowed value"
     end if
     max_band = nband
  endif
  if (max_band_set .and. (max_band < max(1,min_band))) then
     write(0,"(a,2i5)") " ** Max_band is less than min_band: (max_band, eff min_band):", &
                        max_band, max(1,min_band)
     STOP
  endif

  allocate(k(nk))
  allocate(e(nband,nspin,nk))

  read(bands_u,*) (k(ik),((e(ib,is,ik),ib=1,nband), is=1,nspin), ik=1,nk)

  ! Read in optional strings for the band lines
  if ( gnu_ticks ) then
     read(bands_u,*,iostat=is) nlines
     if ( nlines > 0 ) then
        allocate(listk(nlines))
        allocate(labels(nlines))
     else
        ! Errors in reading the ticks
        write(0,'(a)') '*** Could not read number of labels used for GNUplot...'
        gnu_ticks = .false.
     end if
  end if
  if ( gnu_ticks ) then
     do ik = 1 , nlines
        read(bands_u,*,iostat=is) listk(ik), labels(ik)
     end do
  end if
  
  ! We can not close std-in
  if ( bands_u /= 5 ) then
     close(bands_u)
  end if

  if ( Fermi_shift ) then
     e = e - ef
  end if

  nbands = max_band - min_band + 1

  write(6,"(2a)") '# GNUBANDS: Utility for SIESTA to transform ',  &
       'bands output into Gnuplot format'
  write(6,"(a)") '#'
  write(6,"(2a)") '#                                           ',  &
       '       Emilio Artacho, Feb. 1999'
  write(6,"(2a)") '#                                           ',  &
       '        Alberto Garcia, May 2012'
  write(6,"(2a)") '#                                           ',  &
       'Nick Papior Andersen, April 2013'
  write(6,"(2a)") '# ------------------------------------------', &
       '--------------------------------'
  if ( Fermi_shift ) then
     write(6,"(a)")  '# Energies are subtracted Ef (zero is Ef)'
  end if
  if ( spin_idx > 0 ) then
     write(6,"(a,i0)")  '# Only bands for spin ',spin_idx
  else if ( nspin > 1 ) then
     write(6,"(a)")  '# Bands for all spins'
  end if
  write(6,"(a,f10.4)")  '# E_F               = ', Ef
  write(6,"(a,2f10.4)") '# k_min, k_max      = ', kmin, kmax
  write(6,"(a,2f10.4)") '# E_min, E_max      = ', emin, emax
  write(6,"(a,3i6)")    '# Nbands, Nspin, Nk = ', nband, nspin, nk
  write(6,"(a,2i6)")    '# Using min_band, max_band = ', min_band, max_band
  write(6,"(a,i6)")     '# Total number of bands = ', nbands
  write(6,"(a)") '#'
  write(6,"(a)") '#        k            E'
  write(6,"(2a,/)") '# ------------------------------------------',   &
       '--------------------------------'

  delta = 1.0e-5_dp
  if ( spin_idx > 0 ) then
     is = spin_idx
     do ib = min_band, max_band
        add_new_line = .false.
        do ik = 1 , nk
           ! We will only write out in an energy range
           if ( emin-delta <= e(ib,is,ik) .and. e(ib,is,ik) <= emax+delta ) then
              add_new_line = .true.
              ! write spin variable to differentiate
              write(6,"(2f14.6,i3)") k(ik), e(ib,is,ik), is
           end if
        end do
        ! If the energy range has no contribution of this
        ! band, then do not add new-lines
        if ( add_new_line ) write(6,'(/)')
     end do
  else
     do is = 1, nspin
        do ib = min_band, max_band
           add_new_line = .false.
           do ik = 1 , nk
              ! We will only write out in an energy range
              if ( emin-delta <= e(ib,is,ik) .and. e(ib,is,ik) <= emax+delta ) then
                 add_new_line = .true.
                 ! write spin variable to differentiate
                 write(6,"(2f14.6,i3)") k(ik), e(ib,is,ik), is
              end if
           end do
           ! If the energy range has no contribution of this
           ! band, then do not add new-lines
           if ( add_new_line ) write(6,'(/)')
        end do
     end do
  end if

  if ( gnu_ticks ) then
     ! Print out the tick-marks on stderr
     write(0,'(a)',advance='no') 'set xtics ('
     do ik = 1 , nlines
        write(0,'(a,tr1,f9.6)',advance='no') &
             '"'//trim(labels(ik))//'"', listk(ik)
        if ( ik < nlines ) write(0,'(a)',advance='no') ', '
     end do
     write(0,'(a)') ')'
     write(0,'(a)') 'plot "bands.dat" using 1:2:3 with lines lc variable'
     write(0,'(a)') '# -- Use line below for single-color'
     write(0,'(a)') '#plot "bands.dat" with lines'
  end if

contains

  subroutine manual()
    write(0,'(a)') ' Usage: gnubands [options] [bandsfile|PIPE]'
    write(0,*) ! new line
    write(0,'(a)') '   bandsfile   : SystemLabel.bands'
    write(0,'(a)') '        PIPE   : < SystemLabel.bands'
    write(0,*) ! new line
    write(0,'(a)') ' Options:'
    write(0,*) ! new line
    write(0,'(a)') '     -h        : print help'
    write(0,'(a)') '     -G        : print GNUplot commands for correct labels to stderr'
    write(0,'(a)') '                 Suggested usage: prog options 2> bands.gplot 1> bands.dat'
    write(0,'(a)') '                  and then: gnuplot -persist bands.gplot'
    write(0,'(a)') '     -s arg    : only plot selected spin bands [1,nspin]'
    write(0,'(a)') '     -F        : shift energy to Fermi-level'
    write(0,'(a)') '     -b arg    : first band to write'
    write(0,'(a)') '     -B arg    : last band to write'
    write(0,'(a)') '     -e arg    : minimum energy to write'
    write(0,'(a)') '               :   If -F set, will be with respect'
    write(0,'(a)') '               :   to Fermi level'
    write(0,'(a)') '     -E arg    : maximum energy to write'
    write(0,'(a)') '               :   Note, see -e'
  end subroutine manual

end program gnubands
