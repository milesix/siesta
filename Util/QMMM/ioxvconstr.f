c====================================================================
      subroutine ioxvconstr( n, u, r, v, s )

c *******************************************************************
c SaveS positions and velocities at each constr step.
c ********** INPUT **************************************************
C  Modules
C
      use precision
      use fdf
      implicit          none

      character         paste*33
      logical           f
      integer           n,s
      real(dp)          r(3,n), v(3,n), u(3,3) 
      external          io_assign, io_close, paste  

c Internal variables and arrays
      character  sname*30, fname*33,count(1:100)*2
      integer    ia,iu,iv,ix
      data  count/ 
     .     '1','2','3','4','5','6','7','8','9','10',
     .     '11','12','13','14','15','16','17','18','19','20',
     .     '21','22','23','24','25','26','27','28','29','30',
     .     '31','32','33','34','35','36','37','38','39','40',
     .     '41','42','43','44','45','46','47','48','49','50',
     .     '51','52','53','54','55','56','57','58','59','60',
     .     '61','62','63','64','65','66','67','68','69','70',
     .     '71','72','73','74','75','76','77','78','79','80',
     .     '81','82','83','84','85','86','87','88','89','90',
     .     '91','92','93','94','95','96','97','98','99','0'/

      if(s.gt.100) stop 'ioxvcosntr: Increase count variable'
c Find name of file
          sname = fdf_string( 'SystemLabel', 'siesta' )
          fname = paste( sname, '.XV.' )
          fname = paste( fname, count(s) )

c       Open file
          call io_assign( iu )
          open( iu, file=fname, form='formatted', status='unknown' )

c       Write data on file
          write(iu,'(3x,3f18.9)') ((u(ix,iv),ix=1,3),iv=1,3) 
          write(iu,*) n
          do ia = 1,n
            write(iu,'(3f18.9,3x,3f18.9)')
     .        (r(ix,ia),ix=1,3),(v(ix,ia),ix=1,3)
          enddo

c       Close file
          call io_close( iu )

      return
      end

