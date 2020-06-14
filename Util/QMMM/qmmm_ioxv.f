      subroutine qmmm_ioxv(slabel, task, istep, n, u, vu, r, v, fxv ,fv)

c *******************************************************************
c Saves positions and velocities.
c ********** INPUT **************************************************
C  Modules
C
      use precision
      use fdf

      implicit          none

      character         task*(*), slabel*30, paste*33
      logical           fxv,fv
      integer           n, istep
      real(dp)          r(3,n), v(3,n), u(3,3), vu(3,3) 
      external          io_assign, io_close, paste

c Internal variables and arrays
      character  fname*40
      integer    ia,iu,iv,ix,na
      logical    frstme
      save       frstme, fname
      data       frstme /.true./

c Find name of file
        if (frstme) then
c Read slabel.fdf
          fname = paste( slabel, '.qmmm.XV' )
          frstme = .false.
        endif

c Choose between read or write
        if (task.eq.'read' .or. task.eq.'READ') then

c       Check if input file exists
        fxv = .false.
          inquire( file=fname, exist=fxv )
          if (fxv) then

c         Open file
            call io_assign( iu )
            open( iu, file=fname, status='old' )      

c         Read data
             write(6,'(/,a)') 
     .      'qmmm_ioxv: Reading coordinates and velocities from file'
            do iv = 1,3
              read(iu,*,err=1000,end=1000) (u(ix,iv),ix=1,3),
     .                           (vu(ix,iv),ix=1,3)
            enddo
            read(iu,*,err=1000,end=1000) na,istep
            if (na.ne.n) stop 'qmmm_ioxv: Wrong number of atoms!!'
            do ia = 1,n
              read(iu,*,err=1000,end=1000)
     .          (r(ix,ia),ix=1,3),(v(ix,ia),ix=1,3)
            enddo

c         Close file
            call io_close( iu )

c         Check for velocities
            do ix=1,3
              do ia=1,n
                if(v(ix,ia).ne.0.0) fv=.true.
              enddo
            enddo

          else
c         If input file not found, go to exit point
          write(6,'(/,a)') 'qmmm_ioxv: WARNING: XV file not found'
            goto 999
          endif

        elseif (task.eq.'write' .or. task.eq.'WRITE') then

c       Open file
          call io_assign( iu )
          open( iu, file=fname, form='formatted', status='unknown' )

c       Write data on file
          write(iu,'(2(3x,3f18.9))') ((u(ix,iv),ix=1,3),
     .                        (vu(ix,iv),ix=1,3),iv=1,3) 
          write(iu,*) n,istep
          do ia = 1,n
            write(iu,'(3f18.9,3x,3f18.9)')
     .        (r(ix,ia),ix=1,3),(v(ix,ia),ix=1,3)
          enddo

c       Close file
          call io_close( iu )

        endif

  999 continue

      return

 1000 stop 'qmmm_ioxv: problem reading from file'

      end
