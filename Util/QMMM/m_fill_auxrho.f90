
module m_fill_auxrho

  use precision
  use sys, only : die

contains

  subroutine fill_auxrho(na_siesta,qmcell,rclas,ntm,ntpl,atm_rcut,auxrho)

    implicit none

    logical :: auxrho_debug = .true.

    integer, intent(in) :: na_siesta
    real(dp), dimension(3,3), intent(in) :: qmcell
    real(dp), dimension(3,na_siesta) :: rclas
    integer, dimension(3), intent(in) :: ntm
    integer, intent(in) :: ntpl
    real(dp), dimension(:), intent(in) :: atm_rcut
    real(grid_p), intent(inout) :: auxrho(ntpl)

    integer :: i, j, k, ia, imesh
    integer :: ix, iy, iz, ixc, iyc, izc
    integer :: pix, piy, piz, nx, ny, nz
    real(dp) :: rix, riy, riz
    real(dp), dimension(3,3) :: A
    real(dp), dimension(3) :: b
    integer, dimension(3) :: ipiv

    real(dp) :: xm, ym, zm

    integer :: info

    do ia = 1,na_siesta
       do i=1,3
          b(i)=rclas(i,ia)
          do j=1,3
             A(i,j)=qmcell(i,j)/ntm(j)
          enddo
       enddo
       call dgesv(3,1,A,3,ipiv,b,3,info);
       if (info/=0) &
             call die('Failed to solve system of equations')
       ixc=int(b(1))
       iyc=int(b(2))
       izc=int(b(3))
       do i=1,3
          b(i)=atm_rcut(ia)
          do j=1,3
             A(i,j)=qmcell(i,j)/ntm(j)
          enddo
       enddo
       call dgesv(3,1,A,3,ipiv,b,3,info);
       if (info/=0) &
             call die('Failed to solve system of equations')
       nx=int(b(1))+1
       ny=int(b(2))+1
       nz=int(b(3))+1
       do iz=-nz,nz
         do iy=-ny,ny
            do ix=-nx,nx
               rix=real(ix,kind=dp)/real(nx,kind=dp)
               riy=real(iy,kind=dp)/real(ny,kind=dp)
               riz=real(iz,kind=dp)/real(nz,kind=dp)
               if (rix*rix+riy*riy+riz*riz>1.0_dp) cycle
               pix=ixc+ix
               pix=mod(pix,ntm(1))
               piy=iyc+iy
               piy=mod(piy,ntm(2))
               piz=izc+iz
               piz=mod(piz,ntm(3))
               imesh = 1 + pix + ntm(1)*piy + ntm(1)*ntm(2)*piz
               auxrho(imesh)=1.0
            enddo
         enddo
      enddo
    enddo

!!$    if (auxrho_debug) then
!!$       open(1333,file='auxrho.xyz')
!!$       k=0
!!$       do iz=0,ntm(3)-1
!!$          do iy=0,ntm(2)-1
!!$             do ix=0,ntm(1)-1
!!$                imesh = 1 + ix + ntm(1)*iy + ntm(1)*ntm(2)*iz
!!$                if (auxrho(imesh).gt.0.0) k=k+1
!!$             enddo
!!$          enddo
!!$       enddo
!!$       write(1333,*)k+3
!!$       write(1333,*)
!!$       write(1333,FMT='(A,3(1X,F8.3))')'O ',rclas(1:3,1)
!!$       write(1333,FMT='(A,3(1X,F8.3))')'H ',rclas(1:3,2)
!!$       write(1333,FMT='(A,3(1X,F8.3))')'H ',rclas(1:3,3)
!!$       do iz=0,ntm(3)-1
!!$          do iy=0,ntm(2)-1
!!$             do ix=0,ntm(1)-1
!!$                imesh = 1 + ix + ntm(1)*iy + ntm(1)*ntm(2)*iz
!!$                if (auxrho(imesh).gt.0.0) then
!!$                   xm = qmcell(1,1) * ix/ntm(1) + qmcell(1,2) * &
!!$                        iy/ntm(2) + qmcell(1,3) * iz/ntm(3)
!!$                   ym = qmcell(2,1) * ix/ntm(1) + qmcell(2,2) * &
!!$                        iy/ntm(2) + qmcell(2,3) * iz/ntm(3)
!!$                   zm = qmcell(3,1) * ix/ntm(1) + qmcell(3,2) * &
!!$                        iy/ntm(2) + qmcell(3,3) * iz/ntm(3)
!!$                   write(1333,FMT='(A,3(1X,F8.3),1X,F8.3)')'X  ',xm,ym,zm,auxrho(imesh)
!!$                endif
!!$             enddo
!!$          enddo
!!$       enddo
!!$       call flush(1333)
!!$
!!$       close(1333)
!!$
!!$       stop
!!$
!!$    endif

  end subroutine fill_auxrho

end module m_fill_auxrho
