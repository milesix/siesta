
module m_qmmm_fdf

  private

  public :: fdf_block_qmmm, chrlen_qmmm

  !
  ! Copyright Alberto Garcia, Jose Soler, 1996, 1997, 1998
  !---
  !     I/O variables for the fdf package. 
  !
  !     In Fortran 90, all this should go in a module...
  !
  !     ndepth: number of open files (maximum maxdepth)
  !     fdf_stack holds their unit numbers.

  integer, parameter :: maxdepth=5

  integer ndepth, fdf_stack(maxdepth)
  !
  !     Unit numbers for input, output, error notification, and
  !     debugging output (the latter active if fdf_debug is true)
  !
  integer, save :: fdf_in, fdf_out, fdf_err, fdf_log
  ! Unit numbers for input, output, error notification, and
  ! debugging output (the latter active if fdf_debug is true)
  logical, save              :: fdf_debug   = .TRUE.,             &
       fdf_debug2  = .TRUE.,             &
       fdf_started = .FALSE.,             &
       fdf_donothing = .FALSE.
  !
  !     Line just read and parsing info
  !
  character(len=132) :: line
  integer, parameter :: maxntokens=50

  integer ntokens
  integer first(maxntokens), last(maxntokens)
  !---
contains
  !-----------------------------------------------------------------------
  subroutine qmmm_fdf_init()
    !
    !     New initialization for fdf. Simplified user interface using
    !     the io package.
    !

    implicit none

    integer debug_level

    integer lun_tmp
    character(len=20) string

!  Internal variables .................................................
    character(len=20) :: filein, fileout
    character(len=150) :: line
    integer :: count, length

    logical debug_input, file_exists

    !
    if (fdf_donothing) return
    !
    !     Prevent the user from opening two head files
    !
    if (fdf_started) then 
       write(fdf_err,'(a)') 'FDF: Head file already set...'
       stop 'HEAD'
    endif

    call io_geterr(fdf_err)

    ndepth = 0

    call io_assign(fdf_out)
    open(unit=fdf_out,file='qmmm_fdf.log',form='formatted', &
         status='unknown')
    rewind(fdf_out)

    filein='qmmm_block_input'
    call fdf_open(filein)
    write(fdf_out,'(/,a,a,a,i3,/)') &
         '#FDF: Opened ',filein, ' for input. Unit:',fdf_in

    fdf_started = .true.

    debug_level = fdf_integer('fdf-debug',0)
    call fdf_setdebug(debug_level)

  end subroutine qmmm_fdf_init
  !
  !---------------------------------------------------
  subroutine fdf_shutdown
    !
    !     Closes the 'head' file
    !
    implicit none

    if (.not. fdf_started) return

    call fdf_refresh
    call io_close(fdf_in)
    call io_close(fdf_out)
    fdf_started = .false.

  end subroutine fdf_shutdown
  !
  logical function fdf_block_qmmm(label,unit)

    !
    !     Returns "true" and the unit number of the file from which to read
    !     the contents of a block if "label" is associated with a block, and
    !     false if not (unit is set to -1 in this case).
    !
    implicit none
    character*(*) label
    integer unit

    character*256 token1, filename
    integer iless

    if (.not.fdf_started) call qmmm_fdf_init()

    fdf_block_qmmm = .false.
    unit = -1

    if (.not. fdf_locate(label)) return

    token1 = line(first(1):last(1))
    if (.not. leqi(token1,'%block')) then
       write(fdf_err,*) 'FDF_BLOCK_QMMM: Not a block:',label
       !
       !        Return instead of stopping
       !
       return
    endif

    iless = fdf_search('<')
    if ((iless .ne. 0) .and. (ntokens .gt. iless)) then
       !
       !           Read block from file
       !
       filename = line(first(iless+1):last(iless+1))
       if (fdf_debug) write(fdf_log,'(2a)') &
            '*Reading block from file', filename
       call fdf_open(filename)
       if (fdf_search('%dump') .ne. 0) call fdf_dumpfile(label)
       fdf_block_qmmm = .true.
       unit = fdf_in
       return
    endif
    !
    !     Standard block in fdf file. Dump contents
    !
    call fdf_dumpblock(label)
    fdf_block_qmmm = .true.
    unit = fdf_in

  end function fdf_block_qmmm
  !
  !-------------------------------------------------------------------
  !
  subroutine fdf_dumpblock(label)
    !     
    !     Dumps block contents starting at the current line
    !     
    implicit none

    character*(*) label

    integer i, lblock
    character*128 token1

    write(fdf_out,'(/,a79)') line
    lblock = 0
120 continue
    if (fdf_getline()) then
       lblock = lblock + 1
       write(fdf_out,'(a79)') line
       token1 = line(first(1):last(1))
       if (.not. leqi(token1,'%endblock')) goto 120
    else
       write(fdf_err,'(a,a,a)') &
            'FDF_LOCATE: Block ', label, ' does not end!'
       stop 'FDF'
    endif
    write(fdf_out,*)
    !     
    !     Sanity check (optional construct %endblock [ Label [Label] ])
    !     
    if ((ntokens .gt. 1) .and.(fdf_search(label) .eq. 0)) then
       write(fdf_err,'(a,a,a)') &
            'FDF_LOCATE: Block ', label, ' does not end!'
       stop 'FDF'
    endif
    !
    !     Backspace the lines read
    !
    do i=1,lblock
       backspace(fdf_in)
    enddo

    return
  end subroutine fdf_dumpblock
  !  
  !-------------------------------------------------------------------
  !     
  subroutine fdf_dumpfile(label)
    !     
    !     Dumps the contents of a file to fdf_out.
    !     The lines are embedded in a %block ... %endblock pair.
    !     
    implicit none

    character*(*) label
    character form*30
    integer length
    !     
    !     Build the right format
    !     
    call chrlen_qmmm(label,0,length)
    write(form,'(a,i2.2,a)') '(a,a',length,',10x,a)'

    write(fdf_out,*)
    write(fdf_out,form) '%block ', label, &
         '# Originally in include file' 
    !     
    rewind(fdf_in)
10  continue
    if (fdf_getline()) then
       write(fdf_out,'(a79)') line
       goto 10
    endif

    write(fdf_out,form) '%endblock ', label, &
         '# Originally in include file' 
    write(fdf_out,*)
    rewind(fdf_in)
    !     
  end subroutine fdf_dumpfile
  !
  subroutine fdf_parse
    !
    !     Processes the input line looking for meaningful tokens.
    !
    implicit none
    !
    logical intoken, instring

    integer c
    integer stringdel
    !
    !     Character statement functions
    !
    integer i
    logical isdigit, isupper, islower, isalpha, &
         isalnum, isextra, istokch
    logical iscomment, isdelstr, isspecial
    !
    isdigit(i) = (i .ge. 48) .and. (i .le. 57)
    isupper(i) = (i .ge. 65) .and. (i .le. 90)
    islower(i) = (i .ge. 97) .and. (i .le. 122)
    isalpha(i) = isupper(i) .or. islower(i)
    isalnum(i) = isdigit(i) .or. isalpha(i)

    !     Extra characters allowed in tokens:  $ % * + & - . / @ ^ _ ~
    isextra(i) = ((i .ge. 36) .and. (i .le. 38)) &
         .or. (i .eq. 42) .or. (i .eq. 43) &
         .or. (i .eq. 45) .or. (i .eq. 46) &
         .or. (i .eq. 47) .or. (i .eq. 64) .or. (i .eq. 94) &
         .or. (i .eq. 95) .or. (i .eq. 126)

    istokch(i) = isalnum(i) .or. isextra(i)
    !
    !     Comments are signaled by:  !  #  ; 
    iscomment(i) = (i.eq.33) .or. (i.eq.35) .or. (i.eq.59)
    !
    !     String delimiters: "  '  `
    isdelstr(i) = (i.eq.34) .or. (i.eq.39) .or. (i.eq.96)
    !
    !     Special characters which are tokens by themselves: <
    isspecial(i) = (i.eq.60)
    !
    !========================================================
    !
    intoken = .false.
    instring = .false.
    ntokens = 0
    stringdel = 0

    do i = 1, len(line)
       c = ichar(line(i:i))

       if (iscomment(c)) then
          ! possible comment...
          if (instring) then
             last(ntokens) = i
          else
             goto 1000
          endif

       else if (istokch(c)) then
          ! character allowed in a token...
          if (.not. intoken) then
             intoken = .true.
             ntokens = ntokens+1
             first(ntokens) = i
          endif
          last(ntokens) = i

       else if (isspecial(c)) then
          ! character that forms a token by itself...
          if (.not. instring) then
             ntokens=ntokens+1
             first(ntokens) = i
             intoken = .false.
          endif
          last(ntokens) = i

       else if (isdelstr(c)) then
          ! string delimiter... make sure it is the right one before
          ! closing the string.
          ! If we are currently in a token, the delimiter is appended to it.

          if (instring) then
             if (c.eq.stringdel) then
                instring = .false.
                intoken = .false.
                stringdel = 0
             else
                last(ntokens) = i
             endif
          else
             if (intoken) then
                last(ntokens) = i
             else
                instring = .true.
                stringdel = c
                intoken = .true.
                ntokens = ntokens+1
                first(ntokens) = i+1
                last(ntokens) = i+1
             endif
          endif

       else
          ! token delimiter...

          if (instring) then
             last(ntokens) = i
          else
             if (intoken) intoken=.false.
          endif
       endif

    enddo

1000 continue

    if (fdf_debug2) then
       write(fdf_log,*) '            ',  ntokens, ' tokens:'
       do i=1,ntokens
          write(fdf_log,*) '                 ', &
               '|',line(first(i):last(i)),'|'
       enddo
    endif

  end subroutine fdf_parse

  !------------------------------------------------------------------
  !
  integer function fdf_search(label)
    !
    !     Performs a case-and-punctuation-insensitive search for 'label'
    !     among the tokens in a line.
    !
    implicit none

    character*(*) label

    integer i

    fdf_search = 0
    do i = 1, ntokens
       if (labeleq(label,line(first(i):last(i)))) then
          fdf_search = i
          return
       endif
    enddo

  end function fdf_search
  !
  !----------------------------------------------------------------------
  !
  integer function fdf_integer(label,default)
    !
    !     Returns an integer associated with label, or default if label
    !     is not found in the fdf file.
    !
    implicit none

    character*(*) label
    integer default
    !
    character*10 fmtstr
    !
    fdf_integer = default

    if (.not. fdf_locate(label)) then
       write(fdf_out,'(a,5x,i10,5x,a)') &
            label, default, '# Default value'
       return
    endif

    if (ntokens.eq.1) then
       write(fdf_err,*) 'FDF_INTEGER: No value for ', label
       stop
    endif

    write(fmtstr,9000) last(2)-first(2)+1
9000 format('(i',i2.2,')')
    read(line(first(2):last(2)),fmt=fmtstr) fdf_integer
    write(fdf_out,'(a,5x,i20)') label, fdf_integer

  end function fdf_integer
  !
  !----------------------------------------------------------------
  !
  logical function labeleq(s1,s2)
    !
    !     Compares s1 and s2 without regard for case, or appearance
    !     of '_', '.', '-'.
    !
    implicit none

    character*(*) s1, s2
    character*256 n1, n2

    call fdf_pack(s1,n1)
    call fdf_pack(s2,n2)
    labeleq=leqi(n1,n2)
    if (fdf_debug) then
       if (labeleq .and. .not. leqi(s1,s2)) &
            write(fdf_log,'(a,/,a,/,a)') &
            '--------- Considered equivalent:', s1, s2
    endif

  end function labeleq

  !-----------------------------

  subroutine fdf_pack(s,n)
    implicit none
    character*(*) s, n
    !
    !     Removes occurrences of '_ .-'  from s1
    !
    character*1 c
    integer i, j
    logical issep
    issep(i) = (i.eq.95) .or. (i.eq.46) .or. (i.eq.45)

    n = ' '
    j = 0
    do i = 1, len(s)
       c = s(i:i)
       if (.not.issep(ichar(c))) then
          j = j+1
          n(j:j) = c
       endif
    enddo

  end subroutine fdf_pack
  !
  logical function fdf_getline()
    implicit none

    read(fdf_in,end=100,err=100,fmt='(a)') line
    fdf_getline = .true.
    if (fdf_debug2) &
         write(fdf_log,'(a,a76)') '> ', line
    call fdf_parse
    return

100 continue
    fdf_getline = .false.

  end function fdf_getline
  !-----------------------------------------------------------------------
  !
  logical function fdf_locate(label)
    !
    !     Searches for label in the fdf hierarchy. If it appears and it
    !     is not part of a comment, the function returns .true. and leaves
    !     the file positioned at the next line. Otherwise, it returns .false.
    !
    !     It supports two kinds of "include" files:
    !
    !     %include filename  
    !     Indicates an unconditional opening of filename for 
    !     further fdf processing.
    !
    !     Label1 Label2 ... < filename  
    !     Indicates that filename should be opened only when 
    !     searching for any of the labels indicated.
    !     'filename' should be an fdf file.
    !
    implicit none

    character*(*) label

    character*256 token1, filename
    integer ilabel, iless
    !
    fdf_locate = .false.
    if (fdf_donothing) return
    !
    call fdf_refresh
    if (fdf_debug) write(fdf_log,'(/,a,1x,a)') &
         'Looking for ', label

    rewind(fdf_in)

10  continue

    if (.not. fdf_getline()) then
       if (ndepth .gt. 1) then
          call fdf_close
          goto 10
       endif
       if (fdf_debug) write(fdf_log,'(a,1x,a)') &
            '*Did not find ', label
       return
    endif
    !
    if (ntokens .eq. 0) goto 10
    !
    token1 = line(first(1):last(1))
    !
    if (leqi(token1,'%include')) then
       !
       !        Include file
       !
       if (ntokens .eq. 1) then
          write(fdf_err,*) 'FDF: No valid filename after %include'
          stop
       endif
       filename = line(first(2):last(2))
       call fdf_open(filename)
       goto 10
    endif

    ilabel = fdf_search(label)

    if (ilabel .ne. 0) then
       !
       !        Label found...
       !
       if (leqi(token1,'%block')) then
          fdf_locate = .true.
          if (fdf_debug) write(fdf_log,'(a,1x,a)')  &
               '*Found ', label
          return
       endif

       iless = fdf_search('<')
       if ((iless .ne. 0) .and. (ntokens .gt. iless)) then
          !
          !           Continue search in other file
          !
          filename = line(first(iless+1):last(iless+1))
          call fdf_open(filename)
          goto 10
       endif
       !
       !        If we reach this point we must be dealing with a line
       !        of the form 'Label Value'. But we are not interested if
       !        the string appears in the "Value" section
       !
       if (ilabel .eq. 1) then
          fdf_locate = .true.
          if (fdf_debug) write(fdf_log,'(a,1x,a)') '*Found ', label
          return
       else
          goto 10
       endif

    else

       goto 10

    endif

  end function fdf_locate
!
!---------------------------------------------------
      subroutine fdf_setdebug(level)
!
!     Debugging levels: 
!     level <=0: nothing
!     level  =1: standard
!     level >=2: exhaustive
!
      implicit none

      integer level

      if (level .le. 0) then

         if (fdf_debug) then
            call io_close(fdf_log)
            fdf_debug = .false.
         endif

      else

         if (.not. fdf_debug) then
            call io_assign(fdf_log)
            open(fdf_log,file='FDF.debug',form='formatted', &
                 status='unknown')
            rewind(fdf_log)
            fdf_debug = .true.
         endif
      endif
      
      fdf_debug2 = (level .ge. 2)

      end subroutine fdf_setdebug
  !
  !---------------------------------------
  !
  subroutine fdf_open(filename)
    implicit none
    !
    !     Opens a file for fdf processing.
    !
    character*(*) filename

    integer lun
    logical file_exists
    !
    ndepth = ndepth + 1
    if (ndepth .gt. maxdepth) then
       write(fdf_err,'(a)') 'FDF: Too many nested fdf files...'
       stop 'DEPTH'
    endif

    if (leqi(filename,'stdin')) then
       lun = 5
       if (fdf_debug) write(fdf_log,'(a,i1,a)')  &
            '--->Reading from Standard Input [depth:', ndepth,'] '

    else

       call io_assign(lun)

       inquire(file=filename,exist=file_exists)
       if (file_exists) then
          open(unit=lun,file=filename,status='old',form='formatted')
          rewind(lun)
          if (fdf_debug) write(fdf_log,'(a,i1,a,a50)') &
               '--->Opened [depth:', ndepth,'] ', filename
       else
          write(fdf_err,'(a,a60)')  &
               'FDF: Cannot open ',filename
       endif
    endif

    fdf_stack(ndepth) = lun
    fdf_in = lun

  end subroutine fdf_open
  !-----------------------------------
  subroutine fdf_close
    implicit none
    !
    !     Closes currently opened fdf file, except if it is the original one.
    !

    if (ndepth .gt. 1) then
       call io_close(fdf_in)
       if (fdf_debug) write(fdf_log,'(a,i1,a)')  &
            '--->Closed [depth:', ndepth,']'

       ndepth = ndepth -1
       fdf_in = fdf_stack(ndepth)
    endif

  end subroutine fdf_close
  !-------------------------------------
  subroutine fdf_refresh
    !
    !     Closes all the open files in the stack (except the first).
    !     Failure to do so would imply that the next Label is searched 
    !     first in the 'deeper' files. fdf_locate calls fdf_refresh 
    !     before doing anything. 
    !
    implicit none

    integer i

    do i = ndepth, 1 , -1
       call fdf_close
    enddo

  end subroutine fdf_refresh
  !
  !-------------
  !
  SUBROUTINE CHRLEN_QMMM(STRING,NCHAR,LCHAR)
    !
    !  CHRLEN_QMMM accepts a STRING of NCHAR characters and returns LCHAR,
    !  the length of the string up to the last nonblank, nonnull.
    !     
    implicit none

    CHARACTER CHAR*1
    CHARACTER STRING*(*)
    integer nchar,lchar
    !
    integer ncopy, i
    !
    NCOPY=NCHAR
    IF(NCOPY.LE.0)NCOPY=LEN(STRING)
    !
    DO 10 I=1,NCOPY
       LCHAR=NCOPY+1-I
       IF(STRING(LCHAR:LCHAR).NE.' '.AND. &
            STRING(LCHAR:LCHAR).NE.CHAR(0))RETURN
10  ENDDO
    LCHAR=0

  END SUBROUTINE CHRLEN_QMMM
  !-------------
  !
  SUBROUTINE CHRCAP(STRING,NCHAR)
    !
    !  CHRCAP accepts a STRING of NCHAR characters and replaces
    !  any lowercase letters by uppercase ones.
    !
    implicit none

    CHARACTER CHAR*1
    integer nchar, ncopy, i, itemp
    LOGICAL   LGE
    LOGICAL   LLE
    CHARACTER STRING*(*)
    !
    NCOPY=NCHAR
    IF(NCOPY.LE.0)NCOPY=LEN(STRING)
    DO 10 I=1,NCOPY
       !
       IF(LGE(STRING(I:I),'a').AND.LLE(STRING(I:I),'z'))THEN
          ITEMP=ICHAR(STRING(I:I))+ICHAR('A')-ICHAR('a')
          STRING(I:I)=CHAR(ITEMP)
       ENDIF
10  ENDDO

  END SUBROUTINE CHRCAP
  !
  !
  LOGICAL FUNCTION LEQI(STRNG1,STRNG2)
    !
    !  Case-insensitive lexical equal-to comparison
    !
    implicit none
    !
    CHARACTER*1   S1,S2
    CHARACTER*(*) STRNG1
    CHARACTER*(*) STRNG2
    !
    integer len1, len2, lenc, i
    LEN1=LEN(STRNG1)
    LEN2=LEN(STRNG2)
    LENC=MIN(LEN1,LEN2)
    !
    LEQI=.FALSE.
    DO 10 I=1,LENC
       S1=STRNG1(I:I)
       S2=STRNG2(I:I)
       CALL CHRCAP(S1,1)
       CALL CHRCAP(S2,1)
       IF(S1.NE.S2)RETURN
10  ENDDO
    ! 
    IF(LEN1.GT.LENC.AND.STRNG1(LENC+1:LEN1).NE.' ')RETURN
    IF(LEN2.GT.LENC.AND.STRNG2(LENC+1:LEN2).NE.' ')RETURN
    LEQI=.TRUE.
    RETURN
  END FUNCTION LEQI
  !
end module m_qmmm_fdf
