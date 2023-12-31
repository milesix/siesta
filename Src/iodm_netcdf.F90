! ---
! Copyright (C) 1996-2014	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
! 
!
module iodm_netcdf

#ifdef CDF
use netcdf
use sys, only: die
use precision, only: dp
use parallel, only: Node, Nodes
use alloc, only: re_alloc, de_alloc

#ifdef MPI
use mpi_siesta
use parallelsubs, only : LocalToGlobalOrb
#endif

implicit none

! These module variables should be put into a derived type, and maybe
! not all of them are really necessary
!
integer norbs_id, nspin_id, nnzs_id, scf_step_id
integer no_s_id, indxuo_id
integer numd_id, row_pointer_id, column_id, dm_id
!

integer  :: ncid = -1        ! Initialize handle to file descriptor to flag value

public :: setup_dm_netcdf_file, write_dm_netcdf
private

CONTAINS

subroutine setup_dm_netcdf_file( maxnd, nbasis, nspin,    &
                                 no_s, indxuo,            &
                                 numd,  listdptr, listd, geom_step)

      integer, intent(in)   ::    maxnd  ! First dimension of listd and dm
      integer, intent(in)   ::    nbasis ! Number of atomic orbitals
      integer, intent(in)   ::    nspin  ! Number of spins 
      integer, intent(in)   ::    no_s   ! Number of orbitals in supercell

      integer, intent(in)   :: indxuo(no_s) ! Mapping of sc to unit cell orbs

      integer, intent(in)   :: numd(nbasis)
      integer, intent(in)   :: listdptr(nbasis)
      integer, intent(in)   :: listd(maxnd)

      integer, intent(in), optional   :: geom_step

      integer               :: norbs, nnzs
      character(len=12)     :: fname

#ifdef MPI
      integer, dimension(:), pointer  :: norbs_node => null()
      integer, dimension(:), pointer  :: nnzs_node  => null()

      integer, dimension(:), pointer  :: numd_buf  => null()
      integer, dimension(:), pointer  :: listd_buf  => null()

      integer, dimension(:), pointer  :: numd_global  => null()
      integer, dimension(:), pointer  :: global_row_pointer  => null()

      integer  :: MPIerror, stat(MPI_STATUS_SIZE), count, BNode
      integer  :: max_norbs, max_nnzs, ipt, io, iog

      nullify ( norbs_node, nnzs_node)
      call re_alloc( norbs_node, 0, Nodes-1, name='norbs_node', routine='iodm_netcdf' )
      call re_alloc( nnzs_node, 0, Nodes-1,  name='nnzs_node', routine='iodm_netcdf' )
      call mpi_gather(nbasis,1,MPI_Integer, norbs_node(:),1,MPI_integer,0,MPI_Comm_World, mpierror)
      if (Node == 0) then
         norbs = sum(norbs_node(0:Nodes-1))
         max_norbs = maxval(norbs_node(0:Nodes-1))
      endif
!
      call mpi_gather(maxnd,1,MPI_Integer, nnzs_node(:) ,1,MPI_integer,0,MPI_Comm_World, mpierror)
      if (Node == 0) then
         call pxfflush(6)
         nnzs = sum(nnzs_node(0:Nodes-1))
         max_nnzs = maxval(nnzs_node(0:Nodes-1))
      endif
#else
      norbs = nbasis
      nnzs  = maxnd
#endif

      if (Node == 0) then
         if (present(geom_step)) then
            write(fname,"(a3,i4.4,a3)") "DM-", geom_step, ".nc"
         else
            write(fname,"(a)") "DM.nc"
         endif
 
         !  The netCDF file descriptor is kept open during the SCF cycle, but should be
         !  closed at the end of it to avoid resource exhaustion in, say, MD runs. 
         !  It is now actually closed at the beginning of the next SCF cycle. 
         !  (Thanks to Nick Papior Andersen)

         if (ncid >= 0) then
            call check( nf90_close(ncid))
         endif
         call check( nf90_create(fname,NF90_CLOBBER,ncid))
!
!     Dimensions
!
      call check( nf90_def_dim(ncid,'norbs',norbs,norbs_id))  !"Number of basis orbitals"
      call check( nf90_def_dim(ncid,'no_s',no_s,no_s_id))      !"Number of orbitals in supercell"
      call check( nf90_def_dim(ncid,'nspin',nspin,nspin_id))   !"Number of spin components"
      call check( nf90_def_dim(ncid,'nnzs',nnzs,nnzs_id))     !"Number of non-zero interactions"
      call check( nf90_def_dim(ncid,'scf_step',NF90_UNLIMITED,scf_step_id)) !"Index of SCF step"
!
!     Variables
!
      call check( nf90_def_var(ncid,'numd',nf90_int,(/norbs_id/),numd_id))
      call check( nf90_put_att(ncid,numd_id,'Description',"Number of interactions of a given orbital"))
      call check( nf90_def_var(ncid,'row_pointer',nf90_int,(/norbs_id/),row_pointer_id))
      call check( nf90_put_att(ncid,row_pointer_id,'Description',"Index (minus 1) of the start of a given row"))
      call check( nf90_def_var(ncid,'column',nf90_int,(/nnzs_id/),column_id))
      call check( nf90_put_att(ncid,column_id,'Description',"Column index of a given element"))
      call check( nf90_def_var(ncid,'dm',nf90_float,(/nnzs_id,nspin_id,scf_step_id/),dm_id))
      call check( nf90_put_att(ncid,dm_id,'Description',"Density matrix"))

      if (norbs /= no_s) then
         call check( nf90_def_var(ncid,'indxuo',nf90_int,(/no_s_id/),indxuo_id))
         call check( nf90_put_att(ncid,indxuo_id,'Description',"Index of equivalent orb in unit cell"))
      endif
      call check( nf90_enddef(ncid))

      endif    ! Node == 0
!
!     Fill-in unchanging values
!
#ifdef MPI
      if (Node == 0) then
         nullify( numd_buf )
         call re_alloc( numd_buf, 1, max_norbs, name='numd_buf', routine='iodm_netcdf' )
         nullify( listd_buf )
         call re_alloc( listd_buf, 1, max_nnzs, name='listd_buf', routine='iodm_netcdf' )
!
         nullify( numd_global, global_row_pointer )
         call re_alloc( numd_global, 1, norbs, name='numd_global', routine='iodm_netcdf' )
         call re_alloc( global_row_pointer, 1, norbs, name='row_poiner_global', routine='iodm_netcdf' )
      endif

      do BNode = 0, Nodes - 1
         if (Node == 0) then
            if (BNode == 0) then
               numd_buf(1:norbs_node(0)) = numd(1:nbasis)     ! Could do with pointer
            else
               call MPI_Recv(numd_buf,norbs_node(BNode),MPI_integer,BNode,BNode, &
                                                        MPI_Comm_World,stat,mpierror)
              !! call MPI_GET_COUNT(stat, MPI_CHARACTER, count, mpierror)
              !! print *, 'Task ', Node ,': Received', count, 'char(s) from task',  &
              !!          stat(MPI_SOURCE), 'with tag',stat(MPI_TAG)
            endif
            do io = 1, norbs_node(BNode)
               call LocalToGlobalOrb(io,BNode,Nodes,iog)
               numd_global(iog) = numd_buf(io)
            enddo
         else      ! Non-master nodes
            if (Node == BNode) then
               call MPI_Send(numd(1:nbasis),nbasis,MPI_integer,0,BNode,MPI_Comm_World,mpierror)
            endif
         endif
      enddo

      if (Node == 0) then
         call check( nf90_put_var(ncid,numd_id,numd_global(1:norbs),count=(/norbs/)))

         ! Compute global row pointer
         global_row_pointer(1) = 0
         do iog=2,norbs
            global_row_pointer(iog) = global_row_pointer(iog-1) + numd_global(iog-1)
         enddo

         call check( nf90_put_var(ncid,row_pointer_id,global_row_pointer,count=(/norbs/)))
      endif

#else
      call check( nf90_put_var(ncid,numd_id,numd,count=(/norbs/)))
      call check( nf90_put_var(ncid,row_pointer_id,listdptr,count=(/norbs/)))
#endif
!
!   Column information
!
#ifdef MPI

      do BNode = 0, Nodes - 1
         if (Node == 0) then
            if (BNode == 0) then
               listd_buf(1:nnzs_node(0)) = listd(1:maxnd)     ! Could do with pointer
            else
               call MPI_Recv(listd_buf,nnzs_node(BNode),MPI_integer,BNode,BNode,    &
                                                        MPI_Comm_World,stat,mpierror)
               !!call MPI_GET_COUNT(stat, MPI_CHARACTER, count, mpierror) 
               !!print *, 'Task ', Node ,': Received', count, 'char(s) from task',    &
               !!         stat(MPI_SOURCE), 'with tag',stat(MPI_TAG)
            endif
            !
            ! Fill in the column information using the proper offsets
            ipt = 0
            do io = 1, norbs_node(BNode)
               call LocalToGlobalOrb(io,BNode,Nodes,iog)
               call check( nf90_put_var(ncid,column_id,listd_buf(ipt+1:ipt+numd_global(iog)),         &
                                                    start=(/global_row_pointer(iog)+1/),        &
                                                    count=(/numd_global(iog)/) ) )
               ipt = ipt + numd_global(iog)
            enddo
         else      ! Non-master nodes
            if (Node == BNode) then
               call MPI_Send(listd(1:maxnd),maxnd,MPI_integer,0,BNode,MPI_Comm_World,mpierror)
            endif
         endif
      enddo
      call de_alloc(norbs_node,name="norbs_node", routine="iodm_netcdf")
      call de_alloc(nnzs_node,name="nnzs_node", routine="iodm_netcdf")
      call de_alloc(numd_global,name="numd_global", routine="iodm_netcdf")
      call de_alloc(global_row_pointer,name="global_row_pointer", routine="iodm_netcdf")
      call de_alloc(numd_buf,name="numd_buf", routine="iodm_netcdf")
      call de_alloc(listd_buf,name="listd_buf", routine="iodm_netcdf")
#else
      call check( nf90_put_var(ncid,column_id,listd,count=(/nnzs/)))
#endif

      if (Node == 0) then
         if (norbs /= no_s) then
            call check( nf90_put_var(ncid,indxuo_id,indxuo,count=(/no_s/)))
         endif

         call check( nf90_sync(ncid))
      endif

!
!     Should we close the file at this point?
!     (... It might be clearer, but see descriptor fix above)
!
end subroutine setup_dm_netcdf_file

subroutine write_dm_netcdf( nbasis, maxnd, nspin, dm, overwrite )

use precision, only : dp

integer, intent(in)   ::    nbasis ! Number of basis orbitals (in this node)
integer, intent(in)   ::    maxnd  ! First dimension of listd and dm
integer, intent(in)   ::    nspin  ! Number of spins (1 or 2)
logical, intent(in), optional  :: overwrite    ! Overwrite info along scf_step dimension

real(dp), intent(in)  :: dm(maxnd, nspin)

integer               :: norbs, nnzs
integer               :: step_no, step_location

#ifdef MPI
      integer, dimension(:), pointer  :: norbs_node => null()
      integer, dimension(:), pointer  :: nnzs_node  => null()

      real(dp), dimension(:), pointer  ::dm_buf  => null()

      integer, dimension(:), pointer  :: numd_global  => null()
      integer, dimension(:), pointer  :: global_row_pointer  => null()

      integer  :: MPIerror, stat(MPI_STATUS_SIZE), count, BNode
      integer  :: max_norbs, max_nnzs, ipt, io, iog, ispin

      nullify ( norbs_node, nnzs_node)
      call re_alloc( norbs_node, 0, Nodes-1, name='norbs_node', routine='iodm_netcdf' )
      call re_alloc( nnzs_node, 0, Nodes-1,  name='nnzs_node', routine='iodm_netcdf' )

      call mpi_gather(nbasis,1,MPI_Integer, norbs_node(:),1,MPI_integer,0,MPI_Comm_World, mpierror)
      call mpi_gather(maxnd,1,MPI_Integer, nnzs_node(:) ,1,MPI_integer,0,MPI_Comm_World, mpierror)
      if (Node == 0) then
         norbs = sum(norbs_node(0:Nodes-1))
         max_norbs = maxval(norbs_node(0:Nodes-1))
         nnzs = sum(nnzs_node(0:Nodes-1))
         max_nnzs = maxval(nnzs_node(0:Nodes-1))
      endif
#else
      norbs = nbasis
      nnzs  = maxnd
#endif

if (Node == 0) then
   call check(nf90_inquire_dimension(ncid, dimid=scf_step_id, len=step_no))
   step_location = step_no + 1
   if (present(overwrite)) then
      if (overwrite) then
         step_location = 1
      endif
   endif
endif

#ifdef MPI
if (Node == 0) then

   nullify( numd_global, global_row_pointer )
   call re_alloc( numd_global, 1, norbs, name='numd_global', routine='iodm_netcdf' )
   call re_alloc( global_row_pointer, 1, norbs, name='row_poiner_global', routine='iodm_netcdf' )

   nullify( dm_buf )
   call re_alloc( dm_buf, 1, max_nnzs, name='listd_buf', routine='iodm_netcdf' )

   call check( nf90_get_var(ncid,numd_id,numd_global(1:norbs),count=(/norbs/)))
   call check( nf90_get_var(ncid,row_pointer_id,global_row_pointer,count=(/norbs/)))
endif

   do ispin = 1, nspin               ! Outer loop to simplify the logic
                                     ! Cannot send non-contiguous arrays
      do BNode = 0, Nodes - 1
         if (Node == 0) then
            if (BNode == 0) then
               dm_buf(1:nnzs_node(0)) = dm(1:maxnd,ispin)     ! Could do with pointer
            else
               call MPI_Recv(dm_buf,nnzs_node(BNode),MPI_double_precision,BNode,ispin,  &
                                                        MPI_Comm_World,stat,mpierror)
               !! call MPI_GET_COUNT(stat, MPI_CHARACTER, count, mpierror)
               !! print *, 'Task ', Node ,': Received', count, 'char(s) from task',        &
               !!         stat(MPI_SOURCE), 'with tag',stat(MPI_TAG)
            endif
            !
            ! Fill in the column information using the proper offsets
            ipt = 0
            do io = 1, norbs_node(BNode)
               call LocalToGlobalOrb(io,BNode,Nodes,iog)
               call check( nf90_put_var(ncid,dm_id,dm_buf(ipt+1:ipt+numd_global(iog)),         &
                                                    start=(/global_row_pointer(iog)+1, ispin, step_location/),        &
                                                    count=(/numd_global(iog), 1, 1 /) ) )
               ipt = ipt + numd_global(iog)
            enddo

         else      ! Non-master nodes

            if (Node == BNode) then
               call MPI_Send(dm(1:maxnd,ispin),maxnd,MPI_double_precision,0,ispin,MPI_Comm_World,mpierror)
            endif

         endif

      enddo           ! Bnode
   enddo        ! ispin

      call de_alloc(norbs_node,name="norbs_node", routine="iodm_netcdf")
      call de_alloc(nnzs_node,name="nnzs_node", routine="iodm_netcdf")
      call de_alloc(numd_global,name="numd_global", routine="iodm_netcdf")
      call de_alloc(global_row_pointer,name="global_row_pointer", routine="iodm_netcdf")
      call de_alloc(dm_buf,name="dm_buf", routine="iodm_netcdf")

#else
      call check( nf90_put_var(ncid, dm_id, dm(1:maxnd,1:nspin),  & 
                              start = (/1, 1, step_location /), &
                              count = (/maxnd, nspin, 1 /) ))
#endif

if (Node == 0) then
   call check( nf90_sync(ncid))
endif

end subroutine write_dm_netcdf

subroutine check(code)
use netcdf, only: nf90_noerr, nf90_strerror
integer, intent(in) :: code
if (code /= nf90_noerr) call die("netCDF error: " // NF90_STRERROR(code))
end subroutine check



#endif
end module iodm_netcdf
