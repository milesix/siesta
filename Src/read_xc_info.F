! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      subroutine read_xc_info()
C
C     Gathers species' exchange-correlation functional information,
C     optionally reads extra settings from the fdf file,
C     and calls setXC to store it
C
      use precision, only : dp
      use gridXC,  only : setXC=>gridxc_setXC
      use fdf
      use parallel,  only : Node
      use sys,       only : die
      use atm_types, only: species, species_info
      use fdf
      
      implicit none

C     LOCAL variables

      integer, parameter      :: MaxFunc = 10
      integer                 :: nXCfunc
      character(len=50)       :: XCauth(MaxFunc)
      character(len=20)       :: XCfunc(MaxFunc)
      real(dp)                :: XCweightX(MaxFunc)
      real(dp)                :: XCweightC(MaxFunc)

        integer            :: n
        integer            :: ni
        integer            :: nn
        integer            :: nr

        type(block_fdf)            :: bfdf
        type(parsed_line), pointer :: pline

        integer isp, nsp
        character(len=2) :: atom_id
        integer :: libxc_packed_code
        type(species_info), pointer :: spp
        
        logical :: same_xc, force_xc_consistency
        logical :: allow_xc_control

        external :: set_xc_species
        
        nsp = size(species)
        spp => species(1)
        atom_id = spp%xc_info%atom_id
        libxc_packed_code = spp%xc_info%libxc_packed_code
        same_xc = .true.
        do isp = 2, nsp
           spp => species(isp)
           if (spp%xc_info%atom_id /= atom_id) same_xc = .false.
           if (spp%xc_info%libxc_packed_code /= libxc_packed_code)
     $          same_xc = .false.
        enddo

        force_xc_consistency = fdf_get("XC-Force-Consistency",.true.)
        allow_xc_control = fdf_get("XC-Allow-Control",.false.)
        
        if (Node==0) write(6,'(/,a,/)')
     $                  "Exchange-correlation information"
        if (same_xc) then
           ! all species have the same XC. Use it globally
           call set_xc_species(spp)
           if (.not. allow_xc_control) then

              if ( (fdf_defined('xc.functional')) .or.
     $             (fdf_defined('xc.authors')) .or.
     $             (fdf_defined('xc.cocktail')) )  then
                 call message("WARNING",
     $                "XC settings in fdf file disregarded")
              endif
              
              RETURN            ! Do not read more xc input
           endif
        else
           if (force_xc_consistency) then
              call die("Inconsistent XC functionals. See manual")
           else
              call message("WARNING",
     $             "Inconsistent species XC functionals")
              call message("WARNING","Will set XC from fdf file")
           endif
        endif
           
!     Read XC functional information from the fdf file, to override
!     the species' setting.        
        
        if (fdf_block('xc.hybrid',bfdf)) then
           call die('XC.hybrid block is deprecated. Use XC.cocktail')
        endif
        
        if (fdf_block('xc.cocktail',bfdf)) then
          if (.not. fdf_bline(bfdf,pline)) then
            call die('setXC: ERROR no data in xc.cocktail block')
          endif
          ni = fdf_bnintegers(pline)

          if (ni .eq. 0) then
            call die('setXC: Number of functionals missing in ' //
     .               'xc.cocktail')
          endif
          nXCfunc = abs(fdf_bintegers(pline,1))
          if (nXCfunc .gt. MaxFunc) then
            call die('setXC: Too many functionals in xc.cocktail')
          endif
          do n= 1, nXCfunc
            if (.not. fdf_bline(bfdf,pline)) then
              call die('setXC: Number of XC functionals does not match')
            endif
            nn = fdf_bnnames(pline)
            nr = fdf_bnreals(pline)

            if (nn .gt. 0) then
              XCfunc(n) = fdf_bnames(pline,1)
            else
              XCfunc(n) = 'LDA'
            endif
            if (nn .gt. 1) then
              XCauth(n) = fdf_bnames(pline,2)
            else
              XCauth(n) = 'PZ'
            endif
            if (nr .gt. 1) then
              XCweightX(n) = fdf_breals(pline,1)
              XCweightC(n) = fdf_breals(pline,2)
            elseif (nr .eq. 1) then
              XCweightX(n) = fdf_breals(pline,1)
              XCweightC(n) = fdf_breals(pline,1)
            else
              XCweightX(n) = 1.0_dp
              XCweightC(n) = 1.0_dp
            endif
          enddo
        else
          nXCfunc = 1 
          XCfunc(1) = fdf_string('xc.functional','LDA')
          XCauth(1) = fdf_string('xc.authors','PZ')
          XCweightX(1) = 1.0_dp
          XCweightC(1) = 1.0_dp
        endif

        ! Set only if explicitly defined in the fdf file
        if ( (fdf_defined('xc.functional')) .or.
     $       (fdf_defined('xc.authors')) .or.
     $       (fdf_defined('xc.cocktail')) )  then

           if (same_xc) then  ! allow_xc_control must be true
              call message("WARNING",
     $             "Species-based XC setting overridden")
           endif

         if (Node .eq. 0) then
          if (nXCfunc > 1) then
             write(6,'(a,/)') "Cocktail exchange-correlation functional"
          else
             write(6,'(a,/)') "Functional information"
          endif
          write(6,'(a,a)') " Number     Functional     Authors  ",
     .      "            Weight(Ex)   Weight(Ec)"
          do n = 1,nXCfunc
          write(6,'(i4,3x,a20,2x,a20,3x,f5.3,8x,f5.3)')
     .        n,XCfunc(n),XCauth(n),XCweightX(n),XCweightC(n)
          enddo
          write(6,*)
        endif

        ! Store information in module

        call setXC (n=nXCfunc, func=XCfunc, auth=XCauth,
     $           wx=XCweightX, wc=XCweightC)

      else
         call die("Need to set XC functional info explicitly")
      endif

      end subroutine read_xc_info
!
!>     Set xc record for gridxc based on the information
!>     for the species
!      
      subroutine set_xc_species(spp)
      use gridXC,  only : setXC => gridxc_setXC
#ifndef SIESTA__NO_LIBXC
      use gridXC,  only : setXC_libxc => gridxc_setXC_libxc
#endif      
      use atm_types, only: species_info
      use precision, only: dp
      use parallel, only: node

      implicit none
      
      type(species_info), intent(in) :: spp
      
      integer :: libxc_id(2), xc_code, x_code, c_code, nfuncs, i
      character(len=10) :: xc_family, xc_authors
      character(len=100) :: xc_str
      character(len=25) :: xc_type
      
      xc_family = spp%xc_info%xc_family
      xc_authors = spp%xc_info%xc_authors
      xc_code = spp%xc_info%libxc_packed_code
      
      if ((xc_family /= "--") .and. (xc_authors /= "--")) then
         call setXC( 1 , (/ xc_family /), (/ xc_authors /),
     $        wx= (/ 1.0_dp /), wc= (/ 1.0_dp /) )
         if (node == 0) then
            write(6,"(a,a,'/',a)") "Using Siesta-lineage XC: ",
     $                      trim(xc_family),
     $                      trim(xc_authors)
         endif
      else
         ! Fall back to libxc
         if (xc_code <= 0) then
            call die("Cannot fall back to libxc codes")
         else
            if (xc_code < 10000) then
               ! Special syntax for single functional
               ! (For example, Teter exch-corr functional: xc_code=0020
               nfuncs = 1
               libxc_id(1) = xc_code
            else
               x_code = xc_code/10000
               c_code = xc_code - 10000*x_code
               nfuncs = 2
               libxc_id = (/ x_code, c_code /)
            endif
#ifndef SIESTA__NO_LIBXC
            if (node == 0) then
               write(6,"(a)") "Using libXC functionals: "
               do i = 1, nfuncs
                  call libxc_info(libxc_id(i),xc_str, xc_type)
                  write(6,"(i4,':',a,2x,'(',a,')')")
     $                 libxc_id(i), trim(xc_str), trim(xc_type)
               enddo
            endif
            call setxc_libxc(nfuncs, libxc_id)
#else
            if (node == 0) then
               write(6,"(a,2(1x,i3))") "Will need libXC functionals: ",
     $                                 (libxc_id(i), i=1,nfuncs)
            endif
            call die("GridXC does not have LibXC support")
#endif            
         endif
      endif
      if (node == 0) write(6,*)
      
      end subroutine set_xc_species

#ifndef SIESTA__NO_LIBXC
      subroutine libxc_info(id,libxc_str,libxc_type)
      use xc_f90_types_m
      use xc_f90_lib_m

      implicit none

      integer, intent(in) :: id
      character(len=*), intent(out) :: libxc_str, libxc_type

      TYPE(xc_f90_pointer_t) :: xc_func
      TYPE(xc_f90_pointer_t) :: xc_info
!     integer :: i
      character(len=120) :: s1, s2
!     type(xc_f90_pointer_t) :: str

      call xc_f90_func_init(xc_func, xc_info, id, XC_UNPOLARIZED)

      select case(xc_f90_info_kind(xc_info))
      case(XC_EXCHANGE)
         libxc_type = 'exchange'
      case(XC_CORRELATION)
         libxc_type = 'correlation'
      case(XC_EXCHANGE_CORRELATION)
         libxc_type =  'exchange-correlation'
      case(XC_KINETIC)
         libxc_type =  'kinetic'
      end select

      call xc_f90_info_name(xc_info, s1)
      select case(xc_f90_info_family(xc_info))
      case (XC_FAMILY_LDA);       write(s2,'(a)') "LDA"
      case (XC_FAMILY_GGA);       write(s2,'(a)') "GGA"
      case (XC_FAMILY_HYB_GGA);   write(s2,'(a)') "Hybrid GGA"
      case (XC_FAMILY_MGGA);      write(s2,'(a)') "MGGA"
      case (XC_FAMILY_HYB_MGGA);  write(s2,'(a)') "Hybrid MGGA"
      case (XC_FAMILY_LCA);       write(s2,'(a)') "LCA"
      end select
      write(libxc_str, '(4a)') trim(s1), ' (', trim(s2), ')'

!     i = 0
!     call xc_f90_info_refs(xc_info, i, str, s1)
!     do while(i >= 0)
!     write(*, '(a,i1,2a)') '[', i, '] ', trim(s1)
!     call xc_f90_info_refs(xc_info, i, str, s1)
!     end do

      call xc_f90_func_end(xc_func)

      end subroutine libxc_info
#endif
