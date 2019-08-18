! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!

program plstm


! Program PLSTM reads a charge density or local density of states
! generated by SIESTA, and simulates STM images at constant-current
! or constant-height mode using the Tersoff-Hamann approximation.

! Written by P. Ordejon, June 2001
!     (basic structure from Plrho package of J.M.Soler)
! Spin capability and restructuring by Alberto García (March 2019)      

! CAVEATS:
!
! This version works assuming that the scanning plane is the XY plane.
! This plane must be perpendicular to the third lattice vector of
! the supercell (the Z direction).

! USAGE:

! This program reads grid files generated from SIESTA or other programs
! (such as ol-stm/wfs2ldos), with information on the local density of
! states (filename.LDOS) and computes a simulated STM image, in the
! Tersoff-Hamann approximation. Two modes are available:
! constant-current (simulated by computing a constant ldos surface
! z=z(x,y)) and constant-height (obtaining the ldos at the tip position
! at a given height).

! The program uses command-line options. Type 'plstm -h' for usage.
  
! It is possible to analyze just the 'charge', or arbitrary components
! of the spin, by using the '-s' or '-v' options.  
! '-v', followed by a line with three real numbers,
! serves to simulate a spin-polarized tip. The scalar product of the
! "tip spin vector" and the "local spin vector" is recorded at each point.

! The program generates some informative output on standard output,
! including valid ranges, and writes one file. In the case of
! 'constant-current' mode, this file's name is SystemLabel.*.CC.STM,
! and contains the X,Y,Z values of the isosurface (a grid of X,Y and
! the value of Z(X,Y) of the isosurface). In the case of
! 'constant-height' mode, the name is SystemLabel.*.CH.STM, and
! contains the values X,Y,RHO for each X,Y of the grid, where RHO is
! the charge (or spin component) computed at the point X,Y,Z (Z being
! the height specified in the input).

! In the filenames, '*' stands for the 'spin code'
! (as entered with the '-s' flag).

  use m_gridfunc, only: monoclinic_z
  use m_gridfunc, only: gridfunc_t, read_gridfunc
  use m_getopts
      
  implicit none

  type(gridfunc_t) :: gf

  integer, parameter :: dp = selected_real_kind(10,100)

  character(len=200) :: opt_arg
  character(len=10)  :: opt_name 
  integer :: nargs, iostat, n_opts, nlabels

  logical ::  debug    = .false.

  real, allocatable :: rho(:), f2d(:,:)

  character  name*75, fform*12, fname*80, oname*80, task*15, &
             mode*25, spin_code*1
  integer :: i, ip, is, j, mesh(3), np, nspin, nt, Ind, iz, iy
  integer :: k1, k2, ic, ix
  real ::    fvalue, sx, sy, sz, x, y
  real ::    fmin, fmax
  real(dp) ::    zmin, zmax, stepz
  integer :: n3, q1, q2
  real(dp) :: cell(3,3), origin(3)
  real(dp) :: dxdm(2,2)  ! 2D plane grid vectors
  real(dp) :: tip_spin(3)
  logical  :: get_current_range, get_height_range
  integer  :: nx, ny
!
!     Process options
!
  ! defaults  
  get_current_range = .false.
  get_height_range = .false.
  spin_code = 'q'   ! default value
  nx = 1
  ny = 1
  n_opts = 0
  do
     call getopts('hdz:i:s:v:X:Y:IH',opt_name,opt_arg,n_opts,iostat)
     if (iostat /= 0) exit
     select case(opt_name)
     case ('d')
        debug = .true.
     case ('z')
        mode = 'constant-height'
        read(opt_arg,*) fvalue
     case ('i')
        mode = 'constant-current'
        read(opt_arg,*) fvalue
     case ('s')
        read(opt_arg,*) spin_code
     case ('v')
        spin_code = 'v'
        read(opt_arg,*) tip_spin(1:3)
     case ('X')
        read(opt_arg,*) nx
     case ('Y')
        read(opt_arg,*) ny
     case ('I')
        get_current_range = .true.
     case ('H')
        get_height_range = .true.
     case ('h')
        call manual()
        STOP
     case ('?',':')
        write(0,*) "Invalid option: ", opt_arg(1:1)
        write(0,*) "Use -h option for manual"
        write(0,*) ""
        call manual()
        STOP
     end select
  enddo

  nargs = command_argument_count()
  nlabels = nargs - n_opts + 1
  if (nlabels /= 1)  then
     write(0,*) "Use -h option for manual"
     write(0,*) ""
     call manual()
     STOP
  endif

  call get_command_argument(n_opts,value=fname,status=iostat)
  if ( iostat /= 0 ) then
     stop "Cannot get LDOS file"
  end if
  
  if (spin_code == 'q') then
     write(0,*) "Using 'total charge' ('q') mode"
  endif

  write(6,*) 'Reading grid data from file ',trim(fname)

  !! SHOULD MAKE SURE OF THE UNITS... although for plots
  !  it should not make much difference

  ! We might want to try to detect file extension and use
  ! the proper reader.
  call read_gridfunc(fname,gf)
         
  np = product(gf%n(1:3))
  nspin = gf%nspin
  cell = gf%cell
  mesh = gf%n
  origin = gf%origin
       
       write(6,*)
       write(6,*) 'Cell vectors (bohr)'
       write(6,*)
       write(6,*) cell(:,1)
       write(6,*) cell(:,2)
       write(6,*) cell(:,3)
       write(6,*)
       write(6,*) 'Grid mesh: ',mesh(1),'x',mesh(2),'x',mesh(3)
       write(6,*)
       write(6,*) 'nspin = ',nspin
       write(6,*)
       write(6,"(a,3f10.5)") 'Box origin (bohr): ', origin(1:3)

       if (.not. monoclinic_z(cell)) then
          write(0,*) 'The cell is not monoclinic with ' // &
                     'c lattice vector along z...'
          stop " *** Unsuitable cell for STM program"
       endif

          ! Note that the last slice with information in the file is
          ! the (n3-1)th plane
          n3 = mesh(3)
          zmax = (n3-1) * cell(3,3) / n3  + origin(3)
          zmin = origin(3)
          if (n3==1) then
             write(6,"(/,a)") "The file contains a single plane of data"
             write(6,"(/,a)") "Only constant-height mode is allowed"
             stepz = 0.0
          else
             stepz = cell(3,3) / n3
          endif
          write(6,"(a,3f12.5)") "Zmin, Zmax (bohr): ", zmin, zmax


          if (mode .eq. 'constant-current') then

             if (n3==1) then
                STOP "Constant-current mode not available for n3=1"
             endif

          ! It does not make physical sense to get topography
          ! for spin components
          if (spin_code /= "q") then
             STOP "Constant-current mode only available " // &
                 "for 'q' (charge)"
          endif

       else if (mode .eq. 'constant-height') then
         
          if (n3==1) then
             write(6,"(/,a)") "The file contains a single plane of data"
             write(6,"(/,a)") "Height set to plane's height"
             fvalue = zmin
          else

             ! Check that value is within range
             if ((fvalue > zmax) .or. (fvalue < zmin)) then
                STOP 'Z requested is beyond box limits'
             endif
          endif

       else 
         write(6,*) 'ERROR: mode must be either constant-current'
         write(6,*) '       or constant-height (in lower case)'
         stop
       endif
       write(6,*)

       allocate(rho(product(mesh(1:3))))
       
       call get_function(nspin, spin_code, mesh, gf%val, &
            rho, tip_spin, fmin, fmax)
       
       select case (spin_code)
       case ( 'q' )
          write(6,*) "Using 'total charge' ('q') mode"
       case ( 'x' )
          write(6,*) "Using 'x-component of spin' ('x') mode"
       case ( 'y' )
          write(6,*) "Using 'y-component of spin' ('y') mode"
       case ( 'z' )
          write(6,*) "Using 'z-component of spin' ('z') mode"
       case ( 's' )
          write(6,*) "Using 'total spin' ('s') mode"
       case ( 'v' )
          write(6,"(a,3f10.5)") "Using a 'polarized tip' ('-v') " // &
               "with spin: ", tip_spin(1:3)
       end select

       write(6,"(a,3g20.8,/)") "Range of values of processed function: ", fmin, fmax

       allocate(f2d(0:mesh(1)-1,0:mesh(2)-1))

       if (mode .eq. 'constant-current') then
          write(6,"(/a)") 'Calculating STM image in Constant Current mode'
          write(6,"(a)") 'The STM image is obtained as the isosurface of'
          write(6,"(a)") 'constant function value =', fvalue,' e/Bohr**3'
          ! Check that value is within range
          if ((fvalue > fmax) .or. (fvalue < fmin)) then
             STOP 'Function value outside range of values in box'
          endif
          oname = trim(fname)// "." // spin_code // '.CC.STM'
         
         ! Generate x,y,z surface (dump to output file)
          call isocharge( stepz, origin, mesh, rho, fvalue, f2d)
          
       else

          ! Dump slice to output file
          write(6,"(/a)") 'Calculating STM image in Constant Height mode'
          write(6,"(a)") 'The STM image is obtained as the value of the'
          write(6,"(a,f9.4,a)") 'charge at a given tip height Z = ', fvalue, 'Bohr'
          oname = trim(fname)// "." // spin_code // '.CH.STM'

          call isoz( stepz, origin, mesh, rho, fvalue, f2d )
       endif

       ! Write 2D file info
       
       OPEN( unit=2, file=oname )
       write(6,*) 'Writing STM image in file ', trim(oname)

          DO IC = 1,2
             DO IX = 1,2
                DXDM(IX,IC) = CELL(IX,IC) / MESH(IC)
             ENDDO
          ENDDO

          ! Possibly use nx, ny multipliers here
                do k1 = 0, nx*mesh(1) - 1
                   do k2 = 0, ny*mesh(2) - 1
                      x = dxdm(1,1) * k1 + dxdm(1,2) * k2
                      y = dxdm(2,1) * k1 + dxdm(2,2) * k2
                      q1 = mod(k1,mesh(1))      
                      q2 = mod(k2,mesh(2))      
                      write(2,*) x, y, f2d(q1,q2)
                   enddo
                   write(2,*)
                enddo
                close(2)
       deallocate(rho,f2d)

  CONTAINS
  subroutine manual()
    write(0,"(a)") " -------------------"
    write(0,"(a)") " Usage: plstm [options] LDOSfile"
    write(0,"(a)") "  "
    write(0,"(a)") " "
    write(0,"(a)") " OPTIONS: "
    write(0,"(a)") " "
    write(0,"(a)") " -h             Print this help"
    write(0,"(a)") " -d             Print debugging info"
    write(0,"(a)") " "
    write(0,"(a)") " -i current     Constant-current calculation"
    write(0,"(a)") "                with 'current' in e/bohr**3 "
    write(0,"(a)") " -z height      Constant-height calculation"
    write(0,"(a)") "                with 'height' in bohr      "
    write(0,"(a)") " "
    write(0,"(a)") " -s {q,x,y,z,s} Spin code (default 'q' for total 'charge')"
    write(0,"(a)") "                (x|y|z) select cartesian components of spin "
    write(0,"(a)") "                s selects total spin magnitude            "
    write(0,"(a)") " -v 'ux uy uz'  Tip spin direction for selection of spin component"
    write(0,"(a)") "                The vector (ux,uy,uz) should be normalized"
    write(0,"(a)") " "
    write(0,"(a)") " -X NX          Request multiple copies of plot domain along X"
    write(0,"(a)") " -Y NY          Request multiple copies of plot domain along Y"
    write(0,"(a)") " "
    write(0,"(a)") " -H             Return range of height (not implemented yet)"
    write(0,"(a)") " -I             Return range of current (not implemented yet)"
    write(0,"(a)") " -------------------"

  end subroutine manual
end program plstm

      subroutine get_function(nspin, spin_code, mesh, f, rho,  &
                             tip_spin, fmin, fmax)

       integer, parameter :: dp = selected_real_kind(10,100)
       integer, parameter :: sp = kind(1.0)

       integer, intent(in)          :: nspin, mesh(3)
       real(sp), intent(in)         :: f(mesh(1)*mesh(2)*mesh(3),nspin)
       character(len=1), intent(in) :: spin_code
       real(dp), intent(in)         :: tip_spin(3)

       real(sp), intent(out) :: rho(mesh(1)*mesh(2)*mesh(3))
       real(sp), intent(out) :: fmin, fmax

       real ::    sx, sy, sz
       integer :: i, np

       np = product(mesh(1:3))
       
       if (nspin == 1) then
          select case (spin_code)
          case ( 'x', 'y', 'z', 's')
             stop "Cannot choose spin for spin-less file"
          case ( 'v')
             stop "Cannot use tip spin for spin-less file"
          case ( 'q')
             rho(:) = f(:,1)
          end select
       else if (nspin == 2) then
          select case (spin_code)
          case ( 'x', 'y')
             stop "Cannot choose x, y comps for collinear file"
          case ( 'v')
             stop "Cannot use tip spin for collinear file"
          case ( 'q')
             rho(:) = f(:,1) + f(:,2)
          case ( 'z')
             rho(:) = f(:,1) - f(:,2)
          case ( 's')
             rho(:) = abs(f(:,1) - f(:,2))
          end select
       else if (nspin == 4) then
          select case (spin_code)
          case ( 'x' )
             rho(:) = 2.0 * f(:,3)
          case ( 'y' )
             rho(:) = 2.0 * f(:,4)
          case ( 'z')
             rho(:) = f(:,1) - f(:,2)
          case ( 'q')
             rho(:) = f(:,1) + f(:,2)
          case ( 's')
             do i = 1, np
                sx = 2.0 * f(i,3)
                sy = 2.0 * f(i,4)
                sz = f(i,1) - f(i,2)
                rho(i) = sqrt(sx**2 + sy**2 + sz**2)
             enddo
          case ( 'v')
             do i = 1, np
                sx = 2.0 * f(i,3)
                sy = 2.0 * f(i,4)
                sz = f(i,1) - f(i,2)
                rho(i) = tip_spin(1) * sx +  &
                         tip_spin(2) * sy +  &
                         tip_spin(3) * sz
             enddo
          end select
       endif
       fmin = minval(rho)
       fmax = maxval(rho)
     end subroutine get_function
      
      ! Some extra optimizations and clarifications are possible in these routines
      ! for a monoclinic cell

      SUBROUTINE ISOCHARGE( stepz, origin, NMESH, F, FVALUE, f2d)
 
! *******************************************************************
! Calculates the surface z=z(x,y) with constant function value.
! The surface is determined by the condition function=value, and
! it is printed in a file as x,y,z. The function must
! be given in a regular 3-D grid of points.
! Notice single precision in this version

! Written by P. Ordejon. June 2001.
!     from plsurf.f (written by J. M. Soler)
! Modified by A. Garcia, March 2019      
! ************************* INPUT ***********************************
! REAL(DP)  STEPZ
! INTEGER NMESH(3)     : Number of mesh divisions of each vector
! REAL    F(:,:,:)     : Function such that F=FVALUE determines
!                        the shape of the solid surface.
! REAL    FVALUE       : Value such that F=FVALUE
!                        determines the shape of the solid surface.
! ************************* OUTPUT **********************************
! f2d function
! *******************************************************************

       IMPLICIT NONE
       integer, parameter :: dp = selected_real_kind(10,100)

       INTEGER, intent(in) ::  NMESH(3)
       REAL(DP), intent(in)  ::  STEPZ
       REAL(DP), intent(in)  ::  origin(3)
       REAL, intent(in)  ::  F(0:nmesh(1)-1,0:nmesh(2)-1,0:nmesh(3)-1)
       REAL, intent(in)  ::  FVALUE
       REAL, intent(out) ::  F2D(0:nmesh(1)-1,0:nmesh(2)-1)

! Local variables and arrays
       LOGICAL  HIGH, ZERO
       INTEGER  IC, IX, K1, K2, K3
       REAL(dp) ZK3
       REAL     f_cur, f_up

       ZERO = .FALSE.

       DO K1 = 0,NMESH(1)-1
          DO K2 = 0,NMESH(2)-1
             ! z-direction is scanned from top to bottom
             ! We assume that the function decreases as z increases
             ! This is reasonable for the charge density, but not
             ! necessarily true for spin components.
             ! We might need to work with absolute values
             DO K3 = NMESH(3)-1,0,-1

                !Find if this point is above FVALUE
                f_cur = f(k1,k2,k3)     
                HIGH = ( f_cur .GT. FVALUE)
 
                if (HIGH) then   ! our point is between k3 and k3+1
                   if (K3 .eq. NMESH(3)-1) then
                      stop 'Surface above box boundary!!!'
                   endif
                   ! Linear interpolation to find z-coordinate of surface
                   f_up = f(k1,k2,k3+1)
                   ! This is a real number in [k3,k3+1]
                   ZK3 = (K3+1) - (FVALUE - f_up) / (f_cur - f_up)
                   f2d(k1,k2) = stepz * zk3 + origin(3)
                   goto 10
                endif

             ENDDO
             ! Note that K3 is zero here, at the end of the loop
             f2d(k1,k2) = origin(3)
             ZERO = .TRUE.

 10       ENDDO
       ENDDO

       IF (ZERO) THEN
          write(6,*) 'WARNING: I could not find the isosurface'
          write(6,*) '   for some X,Y points. For these, Z = ZMIN (0)'
       ENDIF

     END SUBROUTINE ISOCHARGE

     SUBROUTINE ISOZ( stepz, origin, NMESH, F, ZVALUE, f2d)
 
! Calculates the value of the function F at the plane Z=ZVALUE
! The function must be given in a regular 3-D grid of points.
! Notice single precision in this version

! Written by P. Ordejon. June 2001.
! from plsurf.f (written by J. M. Soler)
! ************************* INPUT ***********************************
! REAL(DP)  CELL(3,3)  : Unit cell vectors CELL(ixyz,ivector)
! INTEGER NMESH(3)     : Number of mesh divisions of each vector
! REAL    F(*)         : Function such that F=FVALUE determines
!                        the shape of the solid surface.
! REAL    ZVALUE       : Z level where the function is written
! CHARACTER*80 ONAME   : Output file name
! ************************* OUTPUT **********************************
! real f2d(:,:)
! *******************************************************************

       IMPLICIT NONE

       integer, parameter :: dp = selected_real_kind(10,100)

       INTEGER, intent(in)  :: NMESH(3)
       REAL(dp), intent(in) :: STEPZ
       REAL(dp), intent(in) :: origin(3)
       REAL, intent(in)     :: F(0:nmesh(1)-1,0:nmesh(2)-1,0:nmesh(3)-1)
       REAL, intent(in)     :: ZVALUE
       REAL, intent(out)    :: f2d(0:nmesh(1)-1,0:nmesh(2)-1)

       
       INTEGER  IC, IP, IPM, IX, K1, K2, K3
       REAL(dp) zm, z
       REAL    f_cur, f_up, fv

       if (nmesh(3) == 1) then
          f2d(:,:) = f(:,:,0)
       else
!     Loop on mesh points
        DO K1 = 0,NMESH(1)-1
           DO K2 = 0,NMESH(2)-1
             ! z-direction is scanned from top to bottom
             DO K3 = NMESH(3)-1,0,-1

               ! Calculate Z coordinate of this point:
               Z = stepz * K3 +  origin(3)

               IF (Z .LT. ZVALUE) THEN

                ! Linear interpolation to find the value of F at ZVALUE

                  f_cur = f(k1,k2,k3)
                  f_up = f(k1,k2,k3+1)
                  ZM = stepz * (K3+1) + origin(3)

                  FV = f_up + (f_cur-f_up) * (ZVALUE - ZM) / (Z - ZM)
                  f2d(k1,k2) = FV
             
                  GOTO 10

               ENDIF

            ENDDO

            WRITE(6,*)  'Z = ',ZVALUE,  &
                 ' not found. It is probably outside your cell'
            STOP
 10      ENDDO
      ENDDO

      endif

    END SUBROUTINE ISOZ
      