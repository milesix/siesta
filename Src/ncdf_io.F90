! Module for easy writing of all information
! in a single netCDF file.

! We rely on the NetCDF4 library and can
! utilise parallel IO or standard IO

module ncdf_io_m

  use precision, only : sp, dp, grid_p, i8b
  use parallel, only : Node, Nodes, Blocksize

  use class_OrbitalDistribution
  use class_Sparsity
  use io_sparse_m, only: max_consecutive_sum, Node_Sp_gncol
  use io_sparse_m, only: count_consecutive, count_blocks

#ifdef NCDF_4
  use netcdf_ncdf, ncdf_parallel => parallel
#endif
#ifdef MPI
  use mpi_siesta
#endif
  use alloc

  implicit none

  private

#ifdef NCDF_4
  public :: cdf_init_mesh
  public :: cdf_r_grid
  public :: cdf_w_grid

  public :: cdf_w_sp
  public :: cdf_w_d1D, cdf_w_d2D

  public :: cdf_r_sp
  public :: cdf_r_d1D, cdf_r_d2D

#endif

#ifndef NEW_NCDF
# define NEW_NCDF 1
# endif

# if NEW_NCDF==1
  type :: t_ncdfSparse
    logical          :: IOnode             ! True if I am a writer
    logical          :: memset             ! True if memory is set
    integer, pointer :: wrtrs(:)           ! List of IOnodes
    integer          :: nwriters           ! Number of writers
    integer, pointer :: g_ncol(:)          ! Number of columns of every global orbital
    integer, pointer :: sndnz(:)           ! Number of elements to send to every writer
    integer, pointer :: rcvnz(:)           ! Number of elements to receive from every process
    integer, pointer :: rcvds(:)           ! Receive displacement of elements
    integer          :: l_nnz              ! Local Number of Non zeros
    integer          :: wroff              ! Writer offset
    integer          :: wr_fp              ! Writer first process
    integer          :: wr_fo              ! Writer first orbital
    integer          :: wr_no              ! Writer number of orbitals
    contains
    procedure :: init_Sp
    procedure :: write_Sp
    procedure :: write_d1D
    procedure :: write_d2D
    procedure :: delete_Sp
  end type t_ncdfSparse
  public :: t_ncdfSparse
# elif NEW_NCDF<5
  type :: t_ncdfSparse
    integer          :: nwriters           ! Number of writers
    logical          :: memset             ! True if memory is set
    logical          :: IOnode             ! True if I am a writer
    integer          :: dsize              ! Max data size
    integer          :: frst_p
    integer, pointer :: ncol(:)            ! New distribution of columns
    integer, pointer :: xnco(:)            ! Offset of the received columns
    integer          :: no                 ! New distribution number of orbitals
    integer          :: ocol               ! offset of the new distribution
    integer, pointer :: nsend(:)           ! Number of nonzeros that we have to send
    integer, pointer :: xsend(:)           ! Offset of nonzeros that we have to send
    integer, pointer :: nrecv(:)           ! Number of nonzeros that we have to receive
    integer, pointer :: xrecv(:)           ! Offset of nonzeros that we have to receive
    integer          :: nnz                ! nnz of the new distribution
    integer          :: onnz               ! Offset for nnz of the new distribution
    integer          :: bsize
# if NEW_NCDF==4
    integer          :: g_no
    integer, pointer :: nno(:)
    integer, pointer :: xno(:)
    integer, pointer :: gnz(:)
#endif
    contains
    procedure :: init_Sp
    procedure :: write_Sp
    procedure :: write_d1D
    procedure :: write_d2D
    procedure :: delete_Sp
  end type t_ncdfSparse
  public :: t_ncdfSparse
# elif NEW_NCDF==5
# warning NEW_NCDF==5
  type :: t_ncdfSparse
    integer :: nwriters
    logical :: IOnode               ! True if I am a writer
    integer :: dsize                ! Max data size
    integer :: p_first              ! Owner (in original distribution of the first block in the new distribution
    integer, pointer :: p_nrows(:)  ! Number of rows of every process in the new distribution
    integer, pointer :: p_xrows(:)  ! Accumulated p_xrows
    !integer, pointer :: ibuff(:)    ! Buffer to receive integer data
    integer :: nb
    integer :: nnz                  ! NZs in new distribution
    integer :: maxnnz               ! max NZs
    integer, pointer :: p_nbloc(:)  ! Number of blocks of every process in the new distribution
    integer, pointer :: p_xbloc(:)  ! Accumulated p_nbloc
    integer, pointer :: blonz(:)    ! Number of NZs per block in original distribution
    integer, pointer :: p_blonz(:)  ! Number of NZs per block in new distribution
    integer, pointer :: p_nnz(:)    ! Number of NZs of every process in the new distribution
    !real*8,  pointer :: dbuff(:)    ! Buffer to receive real data

    integer, pointer :: sreq(:)     ! Buffer to contain MPI send request
    integer, pointer :: rreq(:)     ! Buffer to contain MPI receive request
    contains
    procedure :: init_Sp
    procedure :: write_Sp
    procedure :: write_d1D
    procedure :: write_d2D
    procedure :: delete_Sp
  end type t_ncdfSparse
  public :: t_ncdfSparse
#endif


  ! Create the type to hold the mesh distribution data
  type :: tMeshDist
    private
    ! number of mesh divisions in each axis
    integer :: m(3)
    ! mesh box bounds of each node in each direction
    ! box(1,iAxis,iNode)=lower bounds
    ! box(2,iAxis,iNode)=upper bounds
    integer, allocatable :: box(:,:,:)
  end type tMeshDist
  type(tMeshDist), save :: distr

contains

#ifdef NCDF_4

  subroutine cdf_err( ierr )
    implicit none
    integer :: ierr
    if (ierr /= NF90_NOERR) then
      print *, node, trim(nf90_strerror(ierr))
      call MPI_Abort( MPI_COMM_WORLD, 0, ierr )
    end if
  end subroutine cdf_err



  subroutine cdf_init_mesh(mesh,nsm)

    use parallel, ProcY => ProcessorY

    ! mesh divisions (fine points), number of fine points per big points
    integer, intent(in) :: mesh(3), nsm

    ! Local quantities
    integer :: nm(3), nmesh, ntot
    integer :: ProcZ, blocY, blocZ, nremY, nremZ
    integer :: dimX, dimY, dimZ
    integer :: PP, iniY, iniZ, PY, PZ

    if ( .not. allocated(distr%box) ) then
      allocate(distr%box(2,3,0:Nodes-1))
    end if

    nm(1:3) = mesh(1:3) / nsm
    nmesh = product(mesh)
    distr%m(1:3) = mesh(1:3)

    ProcZ = Nodes/ProcY

    blocY = nm(2)/ProcY
    nremY = nm(2) - blocY*ProcY
    blocZ = nm(3)/ProcZ
    nremZ = nm(3) - blocZ*ProcZ

    dimX = nm(1) * nsm

    ntot = 0

    PP   = 0
    iniY = 1
    do PY = 1, ProcY

      dimY = blocY
      if ( PY <= nremY ) dimY = dimY + 1  ! Add extra points starting from the first nodes
      dimY = dimY * nsm                 ! For fine points

      iniZ = 1
      do PZ = 1, ProcZ
        dimZ = blocZ
        if ( PZ <= nremZ ) dimZ = dimZ + 1
        dimZ = dimZ*nsm                 ! For fine points

        distr%box(1,1,PP) = 1
        distr%box(2,1,PP) = dimX
        distr%box(1,2,PP) = iniY
        distr%box(2,2,PP) = iniY + dimY - 1
        distr%box(1,3,PP) = iniZ
        distr%box(2,3,PP) = iniZ + dimZ - 1

        ntot = ntot + dimX * dimY * dimZ

        iniZ = iniZ + dimZ
        PP   = PP + 1

      end do

      iniY = iniY + dimY

    end do

    if (ntot /= nmesh) then
      if (Node == 0) then
        write(6,*) "Nominal npt: ", nmesh, " /= assigned npt:", ntot
      end if
      call die()
    end if

  end subroutine cdf_init_mesh

  ! Reads in a sparsity pattern at the
  ! current position in the file (iu)
  ! The sparsity pattern "sp" will be returned
  ! as populated.
  ! If dist is supplied it will distribute
  ! the sparsity pattern as supplied (this implies Bcast = .true.)
  ! Else if Bcast is true it will b-cast the sparsity 
  ! pattern fully.
  subroutine cdf_r_Sp(ncdf, no, sp, tag, dit, Bcast, gncol)

    ! File handle
    type(hNCDF), intent(inout) :: ncdf
    ! Number of orbitals readed
    integer, intent(in) :: no
    ! Sparsity pattern
    type(Sparsity), intent(inout) :: sp
    ! The tag of the sparsity pattern
    character(len=*), intent(in) :: tag
    ! distribution if needed to be b-cast in a non-global
    ! fashion
    type(OrbitalDistribution), intent(inout), optional :: dit
    ! Bcast the values?
    logical, intent(in), optional :: Bcast
    integer, intent(inout), target, optional :: gncol(no)

    ! Local variables for reading the values
    integer, pointer :: ncol(:) => null()
    integer, pointer :: l_ptr(:) => null()
    integer, pointer :: l_col(:) => null()

    integer :: io, n_nzs, gind, ind, nl, n, i, nb
    logical :: ldit, lBcast, lIO
    integer, pointer :: lgncol(:) => null()
#ifdef MPI
    integer, allocatable :: ibuf(:)
    integer :: gio, max_n
    integer :: MPIerror, BNode
#endif

#ifdef MPI
    ldit = present(dit)
#else
    ldit = .false.
#endif
    lBcast = .false.
    if ( present(Bcast) ) lBcast = Bcast
    lIO = .not. (lBcast .or. ldit)
    if ( .not. lIO ) lIO = (Node == 0)

    if ( lIO ) then

      if ( present(gncol) ) then
        lgncol => gncol
      else
        allocate(lgncol(no))
      end if

      ! First read in number of non-zero 
      ! entries per orbital
      call ncdf_get_var(ncdf,'n_col',lgncol)

    end if

    nl = no

    ! If a distribution is present, then do something
#ifdef MPI
    if ( ldit ) then

      ! Number of local elements
      nl = num_local_elements(dit,no,Node)
      allocate(ncol(nl))

      ! allocate all requests
      nb = count_blocks(dit,no)
      allocate(ibuf(nb))

      ! Distribute it
      gio = 1
      nb = 0
      do while ( gio <= no ) 

        BNode = node_handling_element(dit,gio)

        ! Get number of consecutive orbitals
        n = count_consecutive(dit,no,gio)

        if ( BNode == Node ) then

          io = index_global_to_local(dit,gio,Node)

          if ( Node == 0 ) then
            ncol(io:io-1+n) = lgncol(gio:gio-1+n)
          else
            nb = nb + 1
            call MPI_IRecv(ncol(io), n, MPI_Integer, &
                0, gio, MPI_Comm_World, ibuf(nb), MPIerror)
          end if

        else if ( Node == 0 ) then
          nb = nb + 1
          call MPI_ISend(lgncol(gio), n, MPI_Integer, &
              BNode, gio, MPI_Comm_World, ibuf(nb), MPIerror)

        end if

        gio = gio + n

      end do

      if ( nb > 0 ) then
        call MPI_WaitAll(nb,ibuf,MPI_STATUSES_IGNORE,MPIerror)
      end if
      deallocate(ibuf)

    else if ( lBcast ) then

      ! Everything should be b-casted
      if ( Node == 0 ) then
        ncol => lgncol
      else
        allocate(ncol(nl))
      end if

      ! Bcast everything
      call MPI_Bcast(ncol(1), nl, MPI_Integer, &
          0, MPI_Comm_World, MPIError)

    else if ( lIO ) then

      ncol => lgncol

    end if
#else
    ! Point to the buffer
    ncol => lgncol
#endif

    ! Allocate pointer
    allocate(l_ptr(nl))

    l_ptr(1) = 0
    do io = 2 , nl
      l_ptr(io) = l_ptr(io-1) + ncol(io-1)
    end do

    ! Number of local non-zero elements
    ! (also works for any bcast methods)
    n_nzs = l_ptr(nl) + ncol(nl)

    ! Allocate space
    allocate(l_col(n_nzs))

#ifdef MPI
    if ( ldit ) then

      ! We have a distributed read
      if ( Node == 0 ) then
        max_n = max_consecutive_sum(dit,no,lgncol)
        allocate(ibuf(max_n))
      else
        allocate(ibuf(nb))
      end if

      ! Read in columns
      gio = 1
      gind = 1
      ind = 0
      nb = 0
      do while ( gio <= no ) 

        BNode = node_handling_element(dit,gio)

        ! Get number of consecutive orbitals
        n = count_consecutive(dit,no,gio)

        if ( BNode == Node ) then

          ! Get the local orbital
          io = index_global_to_local(dit,gio,Node)

          if ( Node == 0 ) then

            i = sum(ncol(io:io-1+n))
            call ncdf_get_var(ncdf, 'list_col', l_col(ind+1:ind+i), &
                start=(/gind/), count=(/i/))
            ind = ind + i
            gind = gind + i

          else

            ! count the number of received entities
            i = sum(ncol(io:io-1+n))
            nb = nb + 1
            call MPI_IRecv(l_col(ind+1), i, MPI_Integer, &
                0, gio, MPI_Comm_World, ibuf(nb), MPIerror)
            ind = ind + i

          end if

        else if ( Node == 0 ) then

          i = sum(lgncol(gio:gio-1+n))
          call ncdf_get_var(ncdf, 'list_col', ibuf(1:i), &
              start=(/gind/), count=(/i/))

          call MPI_Send(ibuf(1), i, MPI_Integer, &
              BNode, gio, MPI_Comm_World, MPIerror)

          gind = gind + i

        end if

        gio = gio + n

      end do

      if ( Node == 0 ) then
        if ( .not. present(gncol) ) deallocate(lgncol)
      else
        if ( nb > 0 ) then
          call MPI_Waitall(nb, ibuf, MPI_STATUSES_IGNORE, MPIerror)
        end if
      end if
      deallocate(ibuf)

    else if ( lBcast ) then

      if ( Node == 0 ) then

        ind = sum(ncol(1:no))
        call ncdf_get_var(ncdf,'list_col',l_col(1:ind), &
            count=(/ind/))
        if ( ind /= n_nzs ) then
          call die('Error in reading sparsity pattern, &
              &size not equivalent.')
        end if
      end if

      ! Bcast
      call MPI_Bcast(l_col(1), n_nzs, MPI_Integer, &
          0, MPI_Comm_World, MPIError)

    else
#endif

      ind = sum(ncol(1:no))
      call ncdf_get_var(ncdf, 'list_col', l_col(1:ind), &
          count=(/ind/) )

#ifdef MPI       
    end if
#endif

    ! Create the sparsity pattern
    call newSparsity(sp,nl,no, n_nzs, ncol, l_ptr, l_col, trim(tag))

    ! de-allocate
    deallocate(l_ptr,l_col)
    if ( ldit ) deallocate(ncol)
    if ( lBcast .and. Node /= 0 ) deallocate(ncol)
    if ( lIO .and. .not. present(gncol) ) deallocate(lgncol)

  end subroutine cdf_r_Sp

# if NEW_NCDF==1
  subroutine delete_Sp( this )
  implicit none
  class(t_ncdfSparse) :: this
  if (this%memset) then
    deallocate(this%wrtrs)
    call de_alloc( this%g_ncol, "g_ncol", "t_ncdfSparse" )
    call de_alloc( this%sndnz, "sndnz", "t_ncdfSparse" )
    call de_alloc( this%rcvnz, "rcvnz", "t_ncdfSparse" )
    call de_alloc( this%rcvds, "rcvds", "t_ncdfSparse" )
  endif
  end subroutine delete_Sp

  subroutine init_Sp( this, sp, max_dsize )
  use m_host, only : getAvailMem, getCHost
  implicit none
  class(t_ncdfSparse) :: this
  type(Sparsity), intent(inout) :: sp
  integer :: max_dsize
  integer :: no, g_no, nnzs, nhost, ierr, i, j, n, p
  integer :: psize, off, nblocks, di, mo, re, ex, lblocks, rlim
  integer, pointer :: ncol(:), count(:), displ(:), recv(:), idx(:), &
                      wrinb(:), wridb(:), sndnb(:)
  integer(i8b) :: l_nnzs, g_nnzs, req_mem, ava_mem
  character*32 :: nw
  this%memset =.false.
  this%nwriters = 1
  this%IOnode = Node==0
  if (max_dsize==0) return
  if (Nodes==1) return
#ifdef MPI
  nullify(ncol,count,displ,recv,idx,wrinb,wridb,sndnb)
  call attach( sp, nrows=no, nrows_g=g_no, n_col=ncol, nnzs=nnzs )
  this%l_nnz = nnzs
  ! Number of computing nodes => All processes in the same computing node
  ! share the available memory
  call getCHost( nhost, this%wrtrs )
  ! Required memory
  l_nnzs = nnzs
  call MPI_AllReduce( l_nnzs, g_nnzs, 1, MPI_INTEGER8, MPI_SUM, &
      MPI_COMM_WORLD, ierr )
  req_mem = g_nnzs*max_dsize*2  ! recv(?:nnzs) reor(?:nnz)
 	req_mem = req_mem + g_no*4    ! g_ncol(g_no)
 	req_mem = req_mem + nodes*4*2 ! rcvnz(nodes)  rcvds(nodes)
  req_mem = req_mem + nodes*4   ! sndnz(nwriter)*cpusXhost
  ava_mem = getAvailMem( )
  this%nwriters = req_mem/ava_mem + 1
  call getenv( "SIESTA_NWRITERS", nw )
  if (trim(nw)/='') read(nw,*) this%nwriters
  this%nwriters = min(nhost,max(1,this%nwriters))
  do i=1, this%nwriters
    if (node==this%wrtrs(i)) then
      this%IOnode = .TRUE.
      exit
    endif
  enddo
  ! Allocate "t_ncdfSparse" arrays
  nullify(this%g_ncol,this%sndnz,this%rcvnz,this%rcvds)
  psize = merge(g_no,1,this%IOnode)
  call re_alloc( this%g_ncol, 1, psize, "g_ncol", "t_ncdfSparse" )
  call re_alloc( this%sndnz, 1, this%nwriters, "sndnz", "t_ncdfSparse" )
  psize = merge(Nodes,1,this%IOnode)
  call re_alloc( this%rcvnz, 1, psize, "rcvnz", "t_ncdfSparse" )
  call re_alloc( this%rcvds, 1, psize, "rcvds", "t_ncdfSparse" )
  ! Receive local no of every process
  psize = merge(Nodes,1,this%IOnode)
  call re_alloc( count, 1, psize, "count", "t_ncdfSparse" )
  call re_alloc( displ, 1, psize, "displ", "t_ncdfSparse" )
  nblocks = (g_no+BlockSize-1)/BlockSize
  di = g_no/(BlockSize*Nodes)
  re = mod(g_no,BlockSize*Nodes)
  if (this%IOnode) then
    off = 0
    do i=1, Nodes
      ex = min(re,BlockSize)
      displ(i) = off
      count(i) = di*BlockSize + ex
      re = re - ex
      off = off + count(i)
    enddo
  endif

  ! IOnodes receive global 'ncol'
  psize = merge(g_no,1,this%IOnode)
  call re_alloc( recv, 1, psize, "recv", "t_ncdfSparse" )
  do i=1, this%nwriters
    p = this%wrtrs(i)
    call MPI_Gatherv( ncol, no, MPI_Integer, recv, count, &
          displ, MPI_Integer, p, MPI_COMM_WORLD, ierr )
  enddo
  ! Reorder recv => g_ncol
  if (this%IOnode) then
    call re_alloc( idx, 1, Nodes, "idx", "t_ncdfSparse" )
    p=1
    idx = displ + 1
    do i=1, g_no, BlockSize
      n = min(BlockSize,g_no-i+1)
      j = idx(p)
      this%g_ncol(i:i+n-1) = recv(j:j+n-1)
      idx(p) = j + n
      p = p+1
      if (p>Nodes) p = 1
      !ib = ib + 1
    enddo
    call de_alloc( idx, "idx", "t_ncdfSparse" )
  endif
  ! Split blocks among writers
  call re_alloc( wrinb, 1, this%nwriters, "wrinb", "t_ncdfSparse" )
  call re_alloc( wridb, 1, this%nwriters, "wridb", "t_ncdfSparse" )
  call re_alloc( sndnb, 1, this%nwriters, "sndnb", "t_ncdfSparse" )
  off = 0
  di = nblocks/this%nwriters
  mo = mod(nblocks,this%nwriters)
  do i=1, this%nwriters
    p = this%wrtrs(i)
    wrinb(i) = di + merge(1,0,i<=mo)
    wridb(i) = off
    if (Node==p) then
      this%wroff = sum(this%g_ncol(1:off*BlockSize))+1
      this%wr_fp = 1+mod(off,nodes)
      this%wr_fo = off*BlockSize+1
      this%wr_no = min(wrinb(i)*BlockSize,g_no-off*BlockSize)
    endif
    off = off + wrinb(i)
  enddo
  ! Number of blocks to send to every writer
  sndnb = 0
  off=Node+1
  p=1
  rlim = wridb(p)+wrinb(p)
  lblocks = (no+BlockSize-1)/BlockSize
  do i=1, lblocks
    do while(off>rlim)
      p = p+1
      rlim = wridb(p)+wrinb(p)
    enddo
    sndnb(p) = sndnb(p)+1
    off = off + Nodes
  enddo
  ! Number of elements to send/receive from every process
  this%sndnz = 0
  off = 0
  do i=1, this%nwriters
    n = min(off+sndnb(i)*BlockSize,no)
    this%sndnz(i) = sum(ncol(off+1:n))
    p = this%wrtrs(i)
    call MPI_Gather( this%sndnz(i), 1, MPI_Integer, this%rcvnz, 1, &
        MPI_integer, p, MPI_COMM_WORLD, ierr )
    off = n
  enddo

  if (Node==0) then
    call re_alloc( idx, 1, Nodes, "idx", "t_ncdfSparse" )
    idx = 0
    p = 1
    do i= 1, this%wr_no, Blocksize
      idx(p) = idx(p) + sum(this%g_ncol(i:min(i+Blocksize-1,this%wr_no)))
      p = merge(1,p+1,p==Nodes)
    enddo
    do p=0, Nodes-1
      n = 0
      do i= 1+p*Blocksize, this%wr_no, Nodes*Blocksize
        n = n + sum(this%g_ncol(i:min(i+Blocksize-1,this%wr_no)))
      enddo
    enddo
    call de_alloc( idx, "idx", "t_ncdfSparse" )
  endif
  if (this%IOnode) then
    off = 0
    do i=1, Nodes
      this%rcvds(i) = off
      off = off + this%rcvnz(i)
    enddo
  endif

  call de_alloc( sndnb, "sndnb", "t_ncdfSparse" )
  call de_alloc( wridb, "wridb", "t_ncdfSparse" )
  call de_alloc( wrinb, "wrinb", "t_ncdfSparse" )
  call de_alloc( recv, "recv", "t_ncdfSparse" )
  call de_alloc( count, "count", "t_ncdfSparse" )
  call de_alloc( displ, "displ", "t_ncdfSparse" )
  this%memset =.true.
#endif
  end subroutine init_Sp
  
  subroutine write_Sp( this, grp, sp )
  implicit none
  class(t_ncdfSparse) :: this
  integer :: grp
  type(Sparsity), intent(inout) :: sp
  integer :: var, nnz, isz, i, i_max, ir, p, psize, off, n, ierr, of2
  integer, pointer :: ncol(:), l_col(:), recv(:), reor(:), idx(:)
  ! Write the sparsity to the file...
  if (nodes==1) then
    nullify(ncol,l_col)
    call attach( sp, n_col=ncol, list_col=l_col )
    call cdf_err( nf90_inq_varid( grp, 'n_col', var ) )
    call cdf_err( nf90_put_var( grp, var, ncol ) )
    call cdf_err( nf90_inq_varid( grp, 'list_col', var ) )
    call cdf_err( nf90_put_var( grp, var, l_col ) )
#	ifdef MPI
  else
    if (Node==0) then
      call cdf_err( nf90_inq_varid( grp, 'n_col', var ) )
      call cdf_err( nf90_put_var( grp, var, this%g_ncol ) )
    endif
    nullify(recv,reor,idx,l_col)
    nnz = merge(this%rcvds(nodes)+this%rcvnz(nodes),1,this%IOnode)
    isz = merge(Nodes,1,this%IOnode)
    call re_alloc( recv, 1, nnz, "recv", "t_ncdfSparse" )
    call re_alloc( reor, 1, nnz, "reor", "t_ncdfSparse" )
    call re_alloc( idx, 1, isz, "idx", "t_ncdfSparse" )

    call attach( sp, list_col=l_col )
    off = 1
    do i=1, this%nwriters
      p=this%wrtrs(i)
      n = this%sndnz(i)
      call MPI_Gatherv( l_col(off), n, MPI_Integer, recv, this%rcvnz, &
          this%rcvds, MPI_Integer, p, MPI_COMM_WORLD, ierr )
      off = off + n
    enddo
    if (this%IOnode) then
      idx = this%rcvds
      p = this%wr_fp
      i_max = this%wr_fo+this%wr_no-1
      off = 0
      do i = this%wr_fo, i_max, Blocksize
        ir = min(i_max,i+Blocksize-1)
        n = sum(this%g_ncol(i:ir))
        of2 = idx(p)
        reor(off+1:off+n) = recv(of2+1:of2+n)
        idx(p) = of2+n
        off = off+n
        p = merge(1,p+1,p==Nodes)
      enddo
      call cdf_err( nf90_inq_varid( grp, 'list_col', var ) )
      !if (this%nwriters>1) then
      !  call cdf_err(nf90_var_par_access( grp, var, NF90_COLLECTIVE ))
      !endif
      call cdf_err( nf90_put_var( grp, var, reor, &
          start=(/ this%wroff /), count=(/ nnz /) ) )
    endif
    call de_alloc( recv, "recv", "t_ncdfSparse" )
    call de_alloc( reor, "reor", "t_ncdfSparse" )
    call de_alloc( idx, "idx", "t_ncdfSparse" )
# endif
  endif
  end subroutine write_Sp

  subroutine write_d1D( this, grp, vname, dSp1D )
  use class_dSpData1D
  implicit none
  class(t_ncdfSparse) :: this
  integer :: grp
  character(len=*), intent(in) :: vname
  type(dSpData1D), intent(inout) :: dSp1D
  integer :: var, nnz, isz, off, i, i_max, p, n, ierr, ir, of2
  real(dp), pointer :: a(:), recv(:), reor(:)
  integer,  pointer :: idx(:)
  a => val(dSp1D)
  ! Write the matrix to the file...
  if (nodes==1) then
    call cdf_err( nf90_inq_varid( grp, vname, var ) )
    call cdf_err( nf90_put_var( grp, var, a ) )
#	ifdef MPI
  else
    nullify(recv,reor,idx)
    nnz = merge(this%rcvds(nodes)+this%rcvnz(nodes),1,this%IOnode)
    isz = merge(Nodes,1,this%IOnode)
    call re_alloc( recv, 1, nnz, "recv", "t_ncdfSparse" )
    call re_alloc( reor, 1, nnz, "reor", "t_ncdfSparse" )
    call re_alloc( idx, 1, isz, "idx", "t_ncdfSparse" )

    ! All writers gathers data
    off = 1
    do i=1, this%nwriters
      p=this%wrtrs(i)
      n = this%sndnz(i)
      call MPI_Gatherv( a(off), n, MPI_Double_Precision, recv, &
          this%rcvnz, this%rcvds, MPI_Double_Precision, p, MPI_COMM_WORLD, ierr )
      off = off + n
    enddo
    ! Reords data
    if (this%IOnode) then
      idx = this%rcvds
      p = this%wr_fp
      i_max = this%wr_fo+this%wr_no-1
      off = 0
      do i = this%wr_fo, i_max, Blocksize
        ir = min(i_max,i+Blocksize-1)
        n = sum(this%g_ncol(i:ir))
        of2 = idx(p)
        reor(off+1:off+n) = recv(of2+1:of2+n)
        idx(p) = of2+n
        off = off+n
        p = merge(1,p+1,p==Nodes)
      enddo
      call cdf_err( nf90_inq_varid( grp, vname, var ) )
      !if (this%nwriters>1) then
      !  call cdf_err(nf90_var_par_access( grp, var, NF90_COLLECTIVE ))
      !endif
      call cdf_err( nf90_put_var( grp, var, reor, &
          start=(/ this%wroff /), count=(/ nnz /) ) )
    endif
    call de_alloc( recv, "recv", "t_ncdfSparse" )
    call de_alloc( reor, "reor", "t_ncdfSparse" )
    call de_alloc( idx, "idx", "t_ncdfSparse" )
# endif
  endif
  end subroutine write_d1D

  subroutine write_d2D( this, grp, vname, dSp2D )
  use class_dSpData2D
  implicit none
  class(t_ncdfSparse) :: this
  integer :: grp
  character(len=*), intent(in) :: vname
  type(dSpData2D), intent(inout) :: dSp2D
  integer :: var, nnz, isz, w, it, off, i, p, n, ierr, i_max, ir, of2
  real(dp), pointer :: a(:,:), recv(:), reor(:)
  integer,  pointer :: idx(:)
  a => val(dSp2D)
  ! Write the sparsity to the file...
  if (nodes==1) then
    call cdf_err( nf90_inq_varid( grp, vname, var ) )
    call cdf_err( nf90_put_var( grp, var, a ) )
#	ifdef MPI
  else
    if (this%IOnode) call cdf_err( nf90_inq_varid( grp, vname, var ) )
    nullify(recv,reor,idx)
    nnz = merge(this%rcvds(nodes)+this%rcvnz(nodes),1,this%IOnode)
    isz = merge(Nodes,1,this%IOnode)
    w = merge(1,size(A,1),size(A,1)==this%l_nnz)
    call re_alloc( recv, 1, nnz*w, "recv", "t_ncdfSparse" )
    call re_alloc( reor, 1, nnz*w, "reor", "t_ncdfSparse" )
    call re_alloc( idx, 1, isz, "idx", "t_ncdfSparse" )
    ! Reduce array to MASTER
    if (size(A,1)==this%l_nnz) then
      do it= 1, size(A,2)
        ! All writers gathers data
        off = 1
        do i=1, this%nwriters
          p=this%wrtrs(i)
          n = this%sndnz(i)
          call MPI_Gatherv( a(off,it), n, MPI_Double_Precision, recv, &
              this%rcvnz, this%rcvds, MPI_Double_Precision, p, MPI_COMM_WORLD, ierr )
          off = off + n
        enddo
        ! Reords data
        if (this%IOnode) then
          idx = this%rcvds
          p = this%wr_fp
          i_max = this%wr_fo+this%wr_no-1
          off = 0
          do i = this%wr_fo, i_max, Blocksize
            ir = min(i_max,i+Blocksize-1)
            n = sum(this%g_ncol(i:ir))
            of2 = idx(p)
            reor(off+1:off+n) = recv(of2+1:of2+n)
            idx(p) = of2+n
            off = off+n
            p = merge(1,p+1,p==Nodes)
          enddo
          !if (this%nwriters>1) then
          !  call cdf_err(nf90_var_par_access( grp, var, NF90_COLLECTIVE ))
          !endif
          call cdf_err( nf90_put_var( grp, var, reor, &
              start=(/ this%wroff, it /), count=(/ nnz, 1 /) ) )
        endif
      enddo
    else
      ! All writers gathers data
      off = 1
      do i=1, this%nwriters
        p = this%wrtrs(i)
        n = this%sndnz(i)
        call MPI_Gatherv( a(1,off), n*w, MPI_Double_Precision, recv, &
            this%rcvnz*w, this%rcvds*w, MPI_Double_Precision, p, &
            MPI_COMM_WORLD, ierr )
        off = off + n
      enddo
      ! Reords data
      if (this%IOnode) then
        idx = this%rcvds*w
        p = this%wr_fp
        i_max = this%wr_fo+this%wr_no-1
        off = 0
        do i = this%wr_fo, i_max, Blocksize
          ir = min(i_max,i+Blocksize-1)
          n = sum(this%g_ncol(i:ir))*w
          of2 = idx(p)
          reor(off+1:off+n) = recv(of2+1:of2+n)
          idx(p) = of2+n
          off = off+n
          p = merge(1,p+1,p==Nodes)
        enddo
        !if (this%nwriters>1) then
        !  call cdf_err(nf90_var_par_access( grp, var, NF90_COLLECTIVE ))
        !endif
        call cdf_err( nf90_put_var( grp, var, reor, &
            start=(/ 1, this%wroff /), count=(/ w, nnz /) ) )
      endif
    endif

    call de_alloc( recv, "recv", "t_ncdfSparse" )
    call de_alloc( reor, "reor", "t_ncdfSparse" )
    call de_alloc( idx, "idx", "t_ncdfSparse" )
# endif
  endif
  end subroutine write_d2D

# elif NEW_NCDF<5
  subroutine delete_Sp( this )
  implicit none
  class(t_ncdfSparse) :: this
  if (this%memset) then
    call de_alloc( this%ncol,  "this%ncol",  "t_ncdfSparse" )
    call de_alloc( this%xnco,  "this%xnco",  "t_ncdfSparse" )
    call de_alloc( this%nsend, "this%nsend", "t_ncdfSparse" )
    call de_alloc( this%xsend, "this%xsend", "t_ncdfSparse" )
    call de_alloc( this%nrecv, "this%nrecv", "t_ncdfSparse" )
    call de_alloc( this%xrecv, "this%xrecv", "t_ncdfSparse" )
  endif
  end subroutine delete_Sp

  subroutine init_Sp( this, sp, max_dsize )
  use m_host, only : getAvailMem, getCHost
  implicit none
  class(t_ncdfSparse) :: this
  type(Sparsity), intent(inout) :: sp
  integer :: max_dsize
  integer :: no, g_no, nbl, g_nbl, di, mo, off, p, nn, i
  integer :: ini, fin, nnz
  integer, pointer :: ncol(:), ndist(:), xdist(:), nsend(:), xsend(:)
  integer, pointer :: nrecv(:), xrecv(:)
  integer :: ierr
  if (max_dsize==0 .or. Nodes==1) then
    this%nwriters = 1
    this%memset = .false.
    this%IOnode = Node==0
    return
  endif
# ifdef MPI
  this%memset = .true.
# if NEW_NCDF==4
  this%nwriters = 1
  this%IOnode = Node==0
# else
  this%nwriters = Nodes
  this%IOnode = .true.
# endif
  nullify(ncol)
  call attach( sp, nrows=no, nrows_g=g_no, n_col=ncol )

  nbl = (no+BlockSize-1)/BlockSize
  g_nbl = (g_no+BlockSize-1)/BlockSize

  ! Compute the new data distribution. Let's use multiples of BlockSize
  ! Only last node would have a non multiple of BlockSize
  nullify(ndist,xdist)
  call re_alloc( ndist, 1, Nodes, "ndist", "t_ncdfSparse" )
  call re_alloc( xdist, 1, Nodes, "xdist", "t_ncdfSparse" )
  di = g_nbl/Nodes
  mo = g_nbl-di*Nodes
  off = 0
  do p=1, Nodes
    nn = merge(di+1,di,p<=mo)*BlockSize
    xdist(p) = off
    ndist(p) = min(nn,g_no-off)
    off = off + nn
  enddo

  ! Compute the columns that we have to send
  nullify(nsend,xsend)
  call re_alloc( nsend, 1, Nodes, "nsend", "t_ncdfSparse" )
  call re_alloc( xsend, 1, Nodes, "xsend", "t_ncdfSparse" )
  nsend(1:Nodes) = 0
  off = Node*BlockSize
  p   = 1
  xsend(p) = 0
  do i=1, no, BlockSize
    do while(off+1>xdist(p)+ndist(p))
      p = p+1
      xsend(p) = i-1
    enddo
    nn = min(BlockSize,no-i+1)
    nsend(p) = nsend(p) + nn
    off = off+Nodes*BlockSize
  enddo
  if (p<Nodes) xsend(p+1:Nodes) = no

  ! Compute the columns that we have to receive
  nullify(nrecv,this%xnco)
  call re_alloc( nrecv, 1, Nodes, "nrecv", "t_ncdfSparse" )
  call re_alloc( this%xnco, 1, Nodes, "this%xnco", "t_ncdfSparse" )
  xrecv => this%xnco
  nrecv(1:Nodes) = 0
  ini = xdist(node+1)
  fin = ini + ndist(node+1)
  this%frst_p = mod(ini/Blocksize,Nodes)+1
  p = this%frst_p
  do i=ini+1, fin, Blocksize
    nn = min(Blocksize,fin-i+1)
    nrecv(p) = nrecv(p)+nn
    p = merge(p+1,1,p<Nodes)
  enddo
  off = 0
  do p=1, Nodes
    xrecv(p) = off
    off = off + nrecv(p)
  enddo

  ! Send/Receive NCOL
  nullify(this%ncol)
  this%no = ndist(node+1)
  call re_alloc( this%ncol, 1, this%no, "this%ncol", "t_ncdfSparse" )
  call MPI_Alltoallv( ncol, nsend, xsend, MPI_Integer, &
      this%ncol, nrecv, xrecv, MPI_Integer, MPI_COMM_WORLD, ierr )

  ! Global Offset of the new distribution
  this%ocol = xdist(node+1)+1

  ! Compute the number of nonzeros that we have to send
  nullify(this%nsend,this%xsend)
  call re_alloc( this%nsend, 1, Nodes, "this%nsend", "t_ncdfSparse" )
  call re_alloc( this%xsend, 1, Nodes, "this%xsend", "t_ncdfSparse" )
  off = 0
  nnz = 0
  do p=1, Nodes
    nn = sum(ncol(off+1:off+nsend(p)))
    this%xsend(p) = nnz
    this%nsend(p) = nn
    off = off + nsend(p)
    nnz = nnz + nn
  enddo

  ! Compute the number of nonzeros that we have to receive
  nullify(this%nrecv,this%xrecv)
  call re_alloc( this%nrecv, 1, Nodes, "this%nrecv", "t_ncdfSparse" )
  call re_alloc( this%xrecv, 1, Nodes, "this%xrecv", "t_ncdfSparse" )
  off = 0
  nnz = 0
  do p=1, Nodes
    nn = sum(this%ncol(off+1:off+nrecv(p)))
    this%xrecv(p) = nnz
    this%nrecv(p) = nn
    off = off + nrecv(p)
    nnz = nnz + nn
  enddo

  ! NNZ and offset of the new distribution
  this%nnz = nnz
  call MPI_Scan( nnz, this%onnz, 1, MPI_Integer, MPI_SUM, &
      MPI_COMM_WORLD, ierr )
  this%onnz = this%onnz-nnz+1
  this%bsize = max(this%no,this%nnz)

# if NEW_NCDF==4
    nullify(this%nno,this%xno,this%gnz)
    if (node==0) then
      call re_alloc( this%nno, 1, Nodes, "this%nno", "t_ncdfSparse" )
      call re_alloc( this%xno, 1, Nodes, "this%xno", "t_ncdfSparse" )
      call re_alloc( this%gnz, 1, Nodes, "this%gnz", "t_ncdfSparse" )
    endif
    call MPI_Gather( this%no, 1, MPI_Integer, this%nno, 1, &
        MPI_integer, 0, MPI_COMM_WORLD, ierr )
    call MPI_Gather( this%nnz, 1, MPI_Integer, this%gnz, 1, &
        MPI_integer, 0, MPI_COMM_WORLD, ierr )
    if (node==0) then
      this%g_no = g_no
      off = 0
      do p=1, Nodes
        this%xno(p) = off
        off = off + this%nno(p)
      enddo
      this%bsize = max(maxval(this%gnz),g_no)
    endif
# endif

  !call de_alloc( xrecv, "xrecv", "t_ncdfSparse" )
  call de_alloc( nrecv, "nrecv", "t_ncdfSparse" )
  call de_alloc( xsend, "xsend", "t_ncdfSparse" )
  call de_alloc( nsend, "nsend", "t_ncdfSparse" )
  call de_alloc( xdist, "xdist", "t_ncdfSparse" )
  call de_alloc( ndist, "ndist", "t_ncdfSparse" )
#endif
  end subroutine init_Sp

  subroutine reord_iarr( frst, indx, indz, no, ncol, rbuf, wbuf )
    integer,    intent(in) :: frst, no, ncol(no), rbuf(*)
    integer, intent(inout) :: indx(nodes), indz(nodes)
    integer,   intent(out) :: wbuf(*)
    integer :: p, o, i, nn, j, nz, k
    p = frst
    o = 1
    do i=1, no, Blocksize
      nn = min(Blocksize,no-i+1)
      j = indx(p)
      nz = sum(ncol(j:j+nn-1))
      k = indz(p)
      wbuf(o:o+nz-1) = rbuf(k:k+nz-1)
      indx(p) = indx(p)+nn
      indz(p) = indz(p)+nz
      o = o+nz
      p = merge(p+1,1,p<Nodes)
    enddo
  end subroutine reord_iarr

  subroutine reord_darr( frst, indx, indz, no, ncol, rbuf, wbuf, sc )
    integer,    intent(in) :: frst, no, ncol(no)
    real(dp),   intent(in) :: rbuf(*)
    integer, intent(inout) :: indx(nodes), indz(nodes)
    real(dp),  intent(out) :: wbuf(*)
    integer, optional :: sc
    integer :: p, o, i, nn, j, nz, k, d
    d = 1
    if (present(sc)) d = sc
    p = frst
    o = 1
    do i=1, no, Blocksize
      nn = min(Blocksize,no-i+1)
      j = indx(p)
      nz = sum(ncol(j:j+nn-1))*d
      k = indz(p)
      wbuf(o:o+nz-1) = rbuf(k:k+nz-1)
      indx(p) = indx(p)+nn
      indz(p) = indz(p)+nz
      o = o+nz
      p = merge(p+1,1,p<Nodes)
    enddo
  end subroutine reord_darr

#	ifdef MPI
# if NEW_NCDF==4
  subroutine write_iarr( grp, var, nz, gnz, wbuf )
  implicit none
  integer, intent(in) :: grp, var, nz, gnz(*)
  integer, intent(inout) ::wbuf(*)
  integer :: o, p, nn, ierr
  integer :: MPIstatus(MPI_STATUS_SIZE)
  if (Node==0) then
    o = 1
    do p=1, Nodes
      nn = gnz(p)
      if (p/=1) then
        call MPI_Recv( wbuf, max(nn,32768), MPI_Integer, &
            p-1, 0, MPI_Comm_World, MPIstatus, ierr )
      endif
      call cdf_err( nf90_put_var( grp, var, wbuf(:nn), &
        start=(/ o /), count=(/ nn /) ) )
      o = o + nn
    enddo
  else
    call MPI_Send( wbuf, nz, MPI_Integer, &
            0, 0, MPI_Comm_World, ierr )
  endif
  end subroutine write_iarr

  subroutine write_darr( grp, var, nz, gnz, wbuf )
  implicit none
  integer, intent(in) :: grp, var, nz, gnz(Nodes)
  real(dp), intent(inout) ::wbuf(*)
  integer :: o, p, nn, ierr
  integer :: MPIstatus(MPI_STATUS_SIZE)
  if (Node==0) then
    o = 1
    do p=1, Nodes
      nn = gnz(p)
      if (p/=1) then
        call MPI_Recv( wbuf, max(nn,16384), MPI_Double_Precision, &
            p-1, 0, MPI_Comm_World, MPIstatus, ierr )
      endif
      call cdf_err( nf90_put_var( grp, var, wbuf(:nn), &
        start=(/ o /), count=(/ nn /) ) )
      o = o + nn
    enddo
  else
    call MPI_Send( wbuf, nz, MPI_Double_Precision, &
          0, 0, MPI_Comm_World, ierr )
  endif
  end subroutine write_darr

  subroutine write_2darr( grp, var, nz, gnz, id, wbuf )
  implicit none
  integer, intent(in) :: grp, var, nz, id, gnz(*)
  real(dp), intent(inout) ::wbuf(*)
  integer :: o, p, nn, ierr, w, st(2), co(2)
  integer :: MPIstatus(MPI_STATUS_SIZE)
  w = merge(-id,1,id<0)
  if (Node==0) then
    o = 1
    do p=1, Nodes
      nn = gnz(p)
      if (p/=1) then
        call MPI_Recv( wbuf, max(w*nn,16384), MPI_Double_Precision, &
            p-1, 0, MPI_Comm_World, MPIstatus, ierr )
      endif
      if (id<0) then
        st = (/ 1, o /)
        co = (/ w, nn /)
      else
        st = (/ o, id /)
        co = (/ nn, 1 /)
      endif
      call cdf_err( nf90_put_var( grp, var, wbuf(:w*nn), &
        start=st, count=co ) )
      o = o + nn
    enddo
  else
    call MPI_Send( wbuf, w*nz, MPI_Double_Precision, &
            0, 0, MPI_Comm_World, ierr )
  endif
  end subroutine write_2darr
# elif NEW_NCDF<4
  subroutine write_iarr( grp, var, o, nz, wbuf )
  implicit none
  integer, intent(in) :: grp, var, o, nz, wbuf(nz)
# if NEW_NCDF==3
  call cdf_err( nf90_var_par_access( grp, var, nf90_collective) )
#	endif
  call cdf_err( nf90_put_var( grp, var, wbuf, start=(/ o /), count=(/ nz /) ) )
  end subroutine write_iarr

  subroutine write_darr( grp, var, o, nz, wbuf )
  implicit none
  integer, intent(in) :: grp, var, o, nz
  real(dp), intent(in) ::wbuf(nz)
# if NEW_NCDF==3
  call cdf_err( nf90_var_par_access( grp, var, nf90_collective) )
#	endif
  call cdf_err( nf90_put_var( grp, var, wbuf, start=(/ o /), count=(/ nz /) ) )
  end subroutine write_darr

  subroutine write_2darr( grp, var, o, nz, id, wbuf )
  implicit none
  integer, intent(in) :: grp, var, o, nz, id
  real(dp), intent(in) ::wbuf(nz)
  integer :: st(2), co(2)
  if (id<0) then
    st = (/ 1, o /)
    co = (/ -id, nz /)
  else
    st = (/ o, id /)
    co = (/ nz, 1 /)
  endif
# if NEW_NCDF==3
  call cdf_err( nf90_var_par_access( grp, var, nf90_collective) )
#	endif
  call cdf_err( nf90_put_var( grp, var, wbuf, start=st, count=co ) )
  end subroutine write_2darr


#	endif /* NEW_NCDF<4 */
# endif /* ifdef MPI */

  subroutine write_Sp( this, grp, sp )
  implicit none
  class(t_ncdfSparse) :: this
  integer :: grp
  type(Sparsity), intent(inout) :: sp
  integer :: ierr, var, p, i, j, k, o, nn, nz, bsize
  integer, pointer :: ncol(:), l_col(:), r_col(:)
  integer, pointer :: wbuf(:), rbuf(:), indx(:), indz(:)
  ! Write the sparsity to the file...
  if (Nodes==1) then
    nullify(ncol,l_col)
    call attach( sp, n_col=ncol, list_col=l_col )
    call cdf_err( nf90_inq_varid( grp, 'n_col', var ) )
    call cdf_err( nf90_put_var( grp, var, ncol ) )
    call cdf_err( nf90_inq_varid( grp, 'list_col', var ) )
    call cdf_err( nf90_put_var( grp, var, l_col ) )
#	ifdef MPI
  else
    bsize = max(this%bsize,32768)
    ! Reorder n_col
    nullify(rbuf,wbuf,indx,indz)
    call re_alloc( wbuf, 1, bsize, "wbuf", "t_ncdfSparse" )
    call re_alloc( rbuf, 1, bsize, "rbuf", "t_ncdfSparse" )
    call re_alloc( indx, 1, Nodes, "indx", "t_ncdfSparse" )
    call re_alloc( indz, 1, Nodes, "indz", "t_ncdfSparse" )
    p = this%frst_p
    indx = this%xnco+1
    do i=1, this%no, Blocksize
      nn = min(Blocksize,this%no-i+1)
      j = indx(p)
      wbuf(i:i+nn-1) = this%ncol(j:j+nn-1)
      indx(p) = indx(p)+nn
      p = merge(p+1,1,p<Nodes)
    enddo
# if NEW_NCDF==4
    call MPI_Gatherv( wbuf, this%no, MPI_Integer, rbuf, this%nno, &
          this%xno, MPI_Integer, 0, MPI_COMM_WORLD, ierr )
    if (node==0) then
      call cdf_err( nf90_inq_varid( grp, 'n_col', var ) )
      call cdf_err( nf90_put_var( grp, var, rbuf(:this%g_no) ) )
    endif
# else  /* NEW_NCDF==4 */
    call cdf_err( nf90_inq_varid( grp, 'n_col', var ) )
#   if NEW_NCDF==3
    call cdf_err( nf90_var_par_access( grp, var, nf90_collective) )
#   endif
    call cdf_err( nf90_put_var( grp, var, wbuf(:this%no), &
          start=(/ this%ocol /), count=(/ this%no /) ) )
# endif
    ! Get the list of columns of current distribution
    nullify(l_col)
    call attach( sp, list_col=l_col )
    ! Send/Receive the list of columns with the new distribution
    call MPI_Alltoallv( l_col, this%nsend, this%xsend, MPI_Integer, &
      rbuf, this%nrecv, this%xrecv, MPI_Integer, MPI_COMM_WORLD, ierr )
    ! Reorder Array
    indx = this%xnco+1
    indz = this%xrecv+1
    call reord_iarr( this%frst_p, indx, indz, &
        this%no, this%ncol, rbuf, wbuf )
    ! Write list_col to NCDF file
    if (this%IOnode) call cdf_err( nf90_inq_varid( grp, 'list_col', var ) )
# if NEW_NCDF==4
    call write_iarr( grp, var, bsize, this%gnz, wbuf )
# else
    call write_iarr( grp, var, this%onnz, this%nnz, wbuf )
# endif
    ! Free local memory
    call de_alloc( indz, "indz", "t_ncdfSparse" )
    call de_alloc( indx, "indx", "t_ncdfSparse" )
    call de_alloc( rbuf, "rbuf", "t_ncdfSparse" )
    call de_alloc( wbuf, "wbuf", "t_ncdfSparse" )
  endif
#	endif
  end subroutine write_Sp

  subroutine write_d1D( this, grp, vname, dSp1D )
  use class_dSpData1D
  implicit none
  class(t_ncdfSparse) :: this
  integer :: grp
  character(len=*), intent(in) :: vname
  type(dSpData1D), intent(inout) :: dSp1D
  integer :: ierr, var, p, o, nn, bsize
  real(dp), pointer :: a(:), wbuf(:), rbuf(:)
  integer, pointer :: indx(:), indz(:)
  a => val(dSp1D)
  ! Write the sparsity to the file...
  if (Nodes==1) then
    call cdf_err( nf90_inq_varid( grp, vname, var ) )
    call cdf_err( nf90_put_var( grp, var, a ) )
#	ifdef MPI
  else
    bsize = max(this%bsize,16384)
    ! Send/Receive the new data distribution
    nullify(rbuf,wbuf,indx,indz)
    call re_alloc( wbuf, 1, bsize, "wbuf", "t_ncdfSparse" )
    call re_alloc( rbuf, 1, bsize, "rbuf", "t_ncdfSparse" )
    call re_alloc( indx, 1, Nodes, "indx", "t_ncdfSparse" )
    call re_alloc( indz, 1, Nodes, "indz", "t_ncdfSparse" )
    call MPI_Alltoallv( a, this%nsend, this%xsend, MPI_Double_Precision, &
      rbuf, this%nrecv, this%xrecv, MPI_Double_Precision, MPI_COMM_WORLD, ierr )
    ! Reorder Array
    indx = this%xnco+1
    indz = this%xrecv+1
    call reord_darr( this%frst_p, indx, indz, &
        this%no, this%ncol, rbuf, wbuf )
    ! Write data to NCDF file
    if (this%IOnode) call cdf_err( nf90_inq_varid( grp, vname, var ) )
#   if NEW_NCDF==4
    call write_darr( grp, var, bsize, this%gnz, wbuf )
#   else
    call write_darr( grp, var, this%onnz, this%nnz, wbuf )
#   endif
    ! Free local memory
    call de_alloc( indz, "indz", "t_ncdfSparse" )
    call de_alloc( indx, "indx", "t_ncdfSparse" )
    call de_alloc( rbuf, "rbuf", "t_ncdfSparse" )
    call de_alloc( wbuf, "wbuf", "t_ncdfSparse" )
  endif
#	endif
  end subroutine write_d1D

  subroutine write_d2D( this, grp, vname, dSp2D )
  use class_dSpData2D
  implicit none
  class(t_ncdfSparse) :: this
  integer :: grp
  character(len=*), intent(in) :: vname
  type(dSpData2D), intent(inout) :: dSp2D
  integer :: ierr, i, var, spdim, d1, nspin, p, o, nn, bsize
  real(dp), pointer :: a(:,:), wbuf(:), rbuf(:)
  integer, pointer :: indx(:), indz(:)
#	ifdef MPI
  integer :: MPIstatus(MPI_STATUS_SIZE)
# endif
  a => val(dSp2D)
  ! Write the sparsity to the file...
  if (Nodes==1) then
    call cdf_err( nf90_inq_varid( grp, vname, var ) )
    call cdf_err( nf90_put_var( grp, var, a ) )
#	ifdef MPI
  else
    spdim=spar_dim(dSp2D)
    d1=merge(size(a,1),1,spdim==2)
    bsize=max(d1*this%bsize,16384)
    nullify(rbuf,wbuf,indx,indz)
    call re_alloc( wbuf, 1, bsize, "wbuf", "t_ncdfSparse" )
    call re_alloc( rbuf, 1, bsize, "rbuf", "t_ncdfSparse" )
    call re_alloc( indx, 1, Nodes, "indx", "t_ncdfSparse" )
    call re_alloc( indz, 1, Nodes, "indz", "t_ncdfSparse" )
    if (this%IOnode) call cdf_err( nf90_inq_varid( grp, vname, var ) )
    if (spdim==2) then
      ! Send/Receive the new data distribution
      call MPI_Alltoallv( a(1,1), d1*this%nsend, d1*this%xsend, MPI_Double_Precision, &
        rbuf, d1*this%nrecv, d1*this%xrecv, MPI_Double_Precision, MPI_COMM_WORLD, ierr )
      ! Reorder Array
      indx = this%xnco+1
      indz = d1*this%xrecv+1
      call reord_darr( this%frst_p, indx, indz, &
          this%no, this%ncol, rbuf, wbuf, d1 )
      ! Write data to NCDF file
#     if NEW_NCDF==4
      call write_2darr( grp, var, bsize/d1, this%gnz, -d1, wbuf )
#     else
      call write_2darr( grp, var, this%onnz, this%nnz, -d1, wbuf )
#     endif
    else
      nspin = size(A,2)
      do i=1, nspin
        ! Send/Receive the new data distribution
        call MPI_Alltoallv( a(:,i), this%nsend, this%xsend, MPI_Double_Precision, &
          rbuf, this%nrecv, this%xrecv, MPI_Double_Precision, MPI_COMM_WORLD, ierr )
        ! Reorder Array
        indx = this%xnco+1
        indz = this%xrecv+1
        call reord_darr( this%frst_p, indx, indz, &
            this%no, this%ncol, rbuf, wbuf )
        ! Write data to NCDF file
#       if NEW_NCDF==4
        call write_2darr( grp, var, bsize, this%gnz, i, wbuf )
#       else
        call write_2darr( grp, var, this%onnz, this%nnz, i, wbuf )
#       endif
      enddo
    endif
    ! Free local memory
    call de_alloc( indz, "indz", "t_ncdfSparse" )
    call de_alloc( indx, "indx", "t_ncdfSparse" )
    call de_alloc( rbuf, "rbuf", "t_ncdfSparse" )
    call de_alloc( wbuf, "wbuf", "t_ncdfSparse" )
  endif
#	endif
  end subroutine write_d2D
# elif NEW_NCDF==5

  subroutine init_Sp( this, sp, dsize )
  use m_host, only : getAvailMem, getCHost
  implicit none
  class(t_ncdfSparse) :: this
  type(Sparsity), intent(inout) :: sp
  integer :: dsize
  integer :: p_nb, g_nb, no, g_no, nnzs
  integer :: ierr, gsize, nw, di, mo, b, gb, o, w, p, n
  integer, pointer :: ncol(:)

  this%nwriters = 1
  this%IOnode   = Node==0
  this%dsize    = dsize
  if (dsize==0 .or. Nodes==1) return

#	ifdef MPI
  call attach( sp, nrows=no, nrows_g=g_no, n_col=ncol, nnzs=nnzs )
  ! Compute the number of nodes that will do the reduction of data
  ! This depends on the infrastructure: Problem size, Number of processors per host, speed of network, etc
  !  We optimize it for Marenostrum IV
  this%nb = (no+Blocksize-1)/Blocksize     ! Number of blocks of current process
  g_nb = (g_no+Blocksize-1)/Blocksize      ! Total number of blocks
  if (Nodes<= 48) then
    gsize = Nodes
  elseif (Nodes<=384 .and. g_nb<1000) then
    gsize = Nodes
  elseif (this%nb>20) then
    gsize = 4
  else
    gsize = 24
  endif
  nw = (Nodes+gsize-1)/gsize
  ! allocate memory
  nullify(this%p_nrows,this%p_xrows,this%p_nbloc,this%p_xbloc)
  call re_alloc( this%p_nrows, 1, Nodes, "p_nrows", "t_ncdfSparse" )
  call re_alloc( this%p_xrows, 1, Nodes, "p_xrows", "t_ncdfSparse" )
  call re_alloc( this%p_nbloc, 1, Nodes, "p_nbloc", "t_ncdfSparse" )
  call re_alloc( this%p_xbloc, 1, Nodes, "p_xbloc", "t_ncdfSparse" )
  ! Compute new distribution. Number of rows for every process
  di = g_nb/nw
  mo = g_nb-nw*di
  o = 0
  w =1
  do p=1, Nodes
    if (mod(p-1,gsize)==0) then
      n = (di + merge(1,0,w<=mo))*Blocksize
      if (o+n>g_no) n = g_no-o
      this%p_nrows(p) = n
      this%p_xrows(p) = o
      o = o+n
      w = w+1
    else
      this%p_nrows(p) = 0
      this%p_xrows(p) = o
    endif
  enddo
  ! Origin process of the new distribution (Fortran numbering)
  this%p_first = mod(this%p_xrows(node+1)/Blocksize,Nodes)+1
  ! Number of blocks of every process in the new distribution
  this%p_nbloc(:) = (this%p_nrows(:)+Blocksize-1)/Blocksize
  o = 0
  do p=1, Nodes
    this%p_xbloc(p) = o
    o = o+this%p_nbloc(p)
  enddo
  ! Compute NZ of every block of local distribution
  p_nb = this%p_nbloc(node+1)
  nullify(this%blonz,this%p_blonz)
  call re_alloc( this%blonz, 1, this%nb, "blonz", "t_ncdfSparse" )
  call re_alloc( this%p_blonz, 1, p_nb, "p_blonz", "t_ncdfSparse" )
  call re_alloc( this%sreq, 1, this%nb, "sreq", "t_ncdfSparse" )
  call re_alloc( this%rreq, 1, p_nb, "rreq", "t_ncdfSparse" )
  o = 1
  do b=1, this%nb
    n = min(Blocksize,no-o+1)
    this%blonz(b) = sum(ncol(o:o+n-1))
    o = o+n
  enddo
  ! Redistribute blonz in the new distribution
  p = 1
  gb = Node+1  ! First global block ID
  do b=1, this%nb
    do while(gb>this%p_nbloc(p)+this%p_xbloc(p))
      p=p+1
    enddo
    call MPI_Isend( this%blonz(b), 1, MPI_Integer, p-1, 0, &
        MPI_COMM_WORLD, this%sreq(b), ierr )
    o = o+n
    gb  = gb+Nodes
  enddo
  if (p_nb>0) then
    p = this%p_first
    do b=1, p_nb
      call MPI_Irecv( this%p_blonz(b), 1, MPI_Integer, &
          p-1, 0, MPI_COMM_WORLD, this%rreq(b), ierr )
      p = merge(1,p+1,p==Nodes)
    enddo
    call MPI_WAitAll( p_nb, this%rreq, MPI_STATUSES_IGNORE, ierr )
  endif
  call MPI_WAitAll( this%nb, this%sreq, MPI_STATUSES_IGNORE, ierr )
  ! Number of Nzs of current distribution
  this%nnz = sum(this%p_blonz)
  nullify(this%p_nnz)
  n = merge(Nodes,0,Node==0)
  call re_alloc( this%p_nnz, 1, n, "p_nnz", "t_ncdfSparse" )
  call MPI_Gather( this%nnz, 1, MPI_Integer, this%p_nnz, 1, &
      MPI_Integer, 0, MPI_COMM_WORLD, ierr )
  this%maxnnz = merge(maxval(this%p_nnz),this%nnz,Node==0)
#	endif /* MPI */
  end subroutine init_Sp

  subroutine delete_Sp( this )
  implicit none
  class(t_ncdfSparse) :: this
  call de_alloc( this%p_nrows, "p_nrows", "t_ncdfSparse" )
  call de_alloc( this%p_xrows, "p_xrows", "t_ncdfSparse" )
  call de_alloc( this%p_nbloc, "p_nbloc", "t_ncdfSparse" )
  call de_alloc( this%p_xbloc, "p_xbloc", "t_ncdfSparse" )
  call de_alloc( this%blonz, "blonz", "t_ncdfSparse" )
  call de_alloc( this%p_blonz, "p_blonz", "t_ncdfSparse" )
  call de_alloc( this%sreq, "sreq", "t_ncdfSparse" )
  call de_alloc( this%rreq, "rreq", "t_ncdfSparse" )
  call de_alloc( this%p_nnz, "p_nnz", "t_ncdfSparse" )
  end subroutine delete_Sp

  subroutine write_Sp( this, grp, sp )
  implicit none
  class(t_ncdfSparse) :: this
  integer :: grp
  type(Sparsity), intent(inout) :: sp
  integer :: var, no, p_no, p_nb
  integer :: o, go, p, b, gb, n, nz
  integer, pointer :: ncol(:), l_col(:), ibuff(:)
# ifdef MPI
  integer :: MPIstatus(MPI_STATUS_SIZE), ierr
# endif
  if (nodes==1) then
    nullify(ncol,l_col)
    call attach( sp, n_col=ncol, list_col=l_col )
    call cdf_err( nf90_inq_varid( grp, 'n_col', var ) )
    call cdf_err( nf90_put_var( grp, var, ncol ) )
    call cdf_err( nf90_inq_varid( grp, 'list_col', var ) )
    call cdf_err( nf90_put_var( grp, var, l_col ) )
#	ifdef MPI
  else
    nullify(ibuff)
    call re_alloc( ibuff, 1, this%maxnnz, "ibuff", "t_ncdfSparse" )
    call attach( sp, nrows=no, n_col=ncol, list_col=l_col )

    o = 1                   ! Local  Offset
    go = Node*Blocksize+1   ! Global offset
    p = 1
    do b=1, this%nb
      do while(go>this%p_nrows(p)+this%p_xrows(p))
        p=p+1
      enddo
      n = merge(Blocksize,no-o+1,b<this%nb)
      call MPI_Isend( ncol(o), n, MPI_Integer, p-1, 0, &
          MPI_COMM_WORLD, this%sreq(b), ierr )
      o = o+n
      go = go+Nodes*Blocksize
    enddo

    p_no = this%p_nrows(Node+1)
    p_nb = (p_no+Blocksize-1)/Blocksize
    if (p_no>0) then
      ncol => ibuff
      p = this%p_first
      o = 1
      do b=1, p_nb
        n = merge(Blocksize,p_no-o+1,b<p_nb)
        call MPI_Irecv( ncol(o), n, MPI_Integer, &
            p-1, 0, MPI_COMM_WORLD, this%rreq(b), ierr )
        o = o+n
        p = merge(1,p+1,p==Nodes)
      enddo
      call MPI_WAitAll( p_nb, this%rreq, MPI_STATUSES_IGNORE, ierr )
    endif
    call MPI_WAitAll( this%nb, this%sreq, MPI_STATUSES_IGNORE, ierr )
    if (node==0) then
      call cdf_err( nf90_inq_varid( grp, 'n_col', var ) )
      do p=1, Nodes
        n = this%p_nrows(p)
        if (n==0) cycle
        o = this%p_xrows(p)+1
        if (p>1) then
          call MPI_Recv( ncol, n, MPI_Integer, &
                 p-1, 0, MPI_Comm_World, MPIstatus, ierr )
        endif
        call cdf_err( nf90_put_var( grp, var, ncol, &
            start=(/o/), count=(/n/) ) )
      enddo
    else if (p_no>0) then
        call MPI_Send( ncol, p_no, MPI_Integer, &
            0, 0, MPI_Comm_World, ierr )
    endif

    o = 1
    gb = Node+1 ! Global block
    p = 1
    do b=1, this%nb
      do while(gb>this%p_nbloc(p)+this%p_xbloc(p))
        p=p+1
      enddo
      nz = this%blonz(b)
      call MPI_Isend( l_col(o), nz, MPI_Integer, p-1, 0, &
          MPI_COMM_WORLD, this%sreq(b), ierr )
      o = o+nz
      gb = gb+Nodes
    enddo
    p_nb = this%p_nbloc(Node+1)
    if (p_nb>0) then
      p = this%p_first
      o = 1
      do b=1, p_nb
        nz = this%p_blonz(b)
        call MPI_Irecv( ibuff(o), nz, MPI_Integer, &
            p-1, 0, MPI_COMM_WORLD, this%rreq(b), ierr )
        o = o+nz
        p = merge(1,p+1,p==Nodes)
      enddo
      call MPI_WAitAll( p_nb, this%rreq, MPI_STATUSES_IGNORE, ierr )
    endif
    call MPI_WAitAll( this%nb, this%sreq, MPI_STATUSES_IGNORE, ierr )
    if (node==0) then
      call cdf_err( nf90_inq_varid( grp, 'list_col', var ) )
      o = 1
      do p=1, Nodes
        nz = this%p_nnz(p)
        if (nz==0) cycle
        if (p>1) then
          call MPI_Recv( ibuff, nz, MPI_Integer, &
                 p-1, 0, MPI_Comm_World, MPIstatus, ierr )
        endif
        call cdf_err( nf90_put_var( grp, var, ibuff, &
            start=(/o/), count=(/nz/) ) )
        o = o+nz
      enddo
    else if (p_nb>0) then
      nz = this%nnz
      call MPI_Send( ibuff, nz, MPI_Integer, &
          0, 0, MPI_Comm_World, ierr )
    endif
    call de_alloc( ibuff, "ibuff", "t_ncdfSparse" )
  endif
#	endif /* MPI */
  end subroutine write_Sp


  subroutine write_d1D( this, grp, vname, dSp1D )
  use class_dSpData1D
  implicit none
  class(t_ncdfSparse) :: this
  integer :: grp
  character(len=*), intent(in) :: vname
  type(dSpData1D), intent(inout) :: dSp1D
  integer :: o, gb, p, b, p_nb, var, nz
  real(dp), pointer :: a(:), dbuff(:)
# ifdef MPI
  integer :: MPIstatus(MPI_STATUS_SIZE), ierr
# endif
  a => val(dSp1D)
  if (nodes==1) then
    call cdf_err( nf90_inq_varid( grp, vname, var ) )
    call cdf_err( nf90_put_var( grp, var, a ) )
#	ifdef MPI
  else
    nullify(dbuff)
    call re_alloc( dbuff, 1, this%maxnnz, "dbuff", "t_ncdfSparse" )

    o = 1
    gb = Node+1 ! Global block
    p = 1
    do b=1, this%nb
      do while(gb>this%p_nbloc(p)+this%p_xbloc(p))
        p=p+1
      enddo
      nz = this%blonz(b)
      call MPI_Isend( a(o), nz, MPI_Double_precision, p-1, 0, &
          MPI_COMM_WORLD, this%sreq(b), ierr )
      o = o+nz
      gb = gb+Nodes
    enddo
    p_nb = this%p_nbloc(Node+1)
    if (p_nb>0) then
      p = this%p_first
      o = 1
      do b=1, p_nb
        nz = this%p_blonz(b)
        call MPI_Irecv( dbuff(o), nz, MPI_Double_precision, &
            p-1, 0, MPI_COMM_WORLD, this%rreq(b), ierr )
        o = o+nz
        p = merge(1,p+1,p==Nodes)
      enddo
      call MPI_WAitAll( p_nb, this%rreq, MPI_STATUSES_IGNORE, ierr )
    endif
    call MPI_WAitAll( this%nb, this%sreq, MPI_STATUSES_IGNORE, ierr )
    if (node==0) then
      call cdf_err( nf90_inq_varid( grp, vname, var ) )
      o = 1
      do p=1, Nodes
        nz = this%p_nnz(p)
        if (nz==0) cycle
        if (p>1) then
          call MPI_Recv( dbuff, nz, MPI_Double_precision, &
                 p-1, 0, MPI_Comm_World, MPIstatus, ierr )
        endif
        call cdf_err( nf90_put_var( grp, var, dbuff, &
            start=(/o/), count=(/nz/) ) )
        o = o+nz
      enddo
    else if (p_nb>0) then
      nz = this%nnz
      call MPI_Send( dbuff, nz, MPI_Double_precision, &
          0, 0, MPI_Comm_World, ierr )
    endif
    call de_alloc( dbuff, "dbuff", "t_ncdfSparse" )
  endif
#	endif /* MPI */
  end subroutine write_d1D

  subroutine write_d2D( this, grp, vname, dSp2D )
  use class_dSpData2D
  implicit none
  class(t_ncdfSparse) :: this
  integer :: grp
  character(len=*), intent(in) :: vname
  type(dSpData2D), intent(inout) :: dSp2D
  integer :: o, gb, p, b, p_nb, var, nz, width, niters, it, spdim
  real(dp), pointer :: a(:,:), dbuff(:,:)
# ifdef MPI
  integer :: MPIstatus(MPI_STATUS_SIZE), ierr
# endif
  a => val(dSp2D)
  if (nodes==1) then
    call cdf_err( nf90_inq_varid( grp, vname, var ) )
    call cdf_err( nf90_put_var( grp, var, a ) )
#	ifdef MPI
  else
    spdim=spar_dim(dSp2D)
    nullify(dbuff)
    if (spdim==1) then
      niters = size(a,2)
      width  = 1
      call re_alloc( dbuff, 1, 1, 1, this%maxnnz,"dbuff", "t_ncdfSparse" )
    else
      niters = 1
      width  = size(a,1)
      call re_alloc( dbuff, 1, width, 1, this%maxnnz, "dbuff", "t_ncdfSparse" )
    endif
    do it=1, niters
      o = 1
      gb = Node+1 ! Global block
      p = 1
      do b=1, this%nb
        do while(gb>this%p_nbloc(p)+this%p_xbloc(p))
          p=p+1
        enddo
        nz = this%blonz(b)
        if (spdim==1) then
          call MPI_Isend( a(o,it), nz*width, MPI_Double_precision, &
              p-1, 0, MPI_COMM_WORLD, this%sreq(b), ierr )
        else
          call MPI_Isend( a(1,o), nz*width, MPI_Double_precision, &
              p-1, 0, MPI_COMM_WORLD, this%sreq(b), ierr )
        endif
        o = o+nz
        gb = gb+Nodes
      enddo
      p_nb = this%p_nbloc(Node+1)
      if (p_nb>0) then
        p = this%p_first
        o = 1
        do b=1, p_nb
          nz = this%p_blonz(b)
          call MPI_Irecv( dbuff(1,o), nz*width, MPI_Double_precision, &
              p-1, 0, MPI_COMM_WORLD, this%rreq(b), ierr )
          o = o+nz
          p = merge(1,p+1,p==Nodes)
        enddo
        call MPI_WAitAll( p_nb, this%rreq, MPI_STATUSES_IGNORE, ierr )
      endif
      call MPI_WAitAll( this%nb, this%sreq, MPI_STATUSES_IGNORE, ierr )
      if (node==0) then
        call cdf_err( nf90_inq_varid( grp, vname, var ) )
        o = 1
        do p=1, Nodes
          nz = this%p_nnz(p)
          if (nz==0) cycle
          if (p>1) then
            call MPI_Recv( dbuff(1,1), nz*width, MPI_Double_precision, &
                   p-1, 0, MPI_Comm_World, MPIstatus, ierr )
          endif
          if (spdim==1) then
            call cdf_err( nf90_put_var( grp, var, dbuff, &
                start=(/o,it/), count=(/nz,1/) ) )
          else
            call cdf_err( nf90_put_var( grp, var, dbuff, &
                start=(/1,o/), count=(/width,nz/) ) )
          endif
          o = o+nz
        enddo
      else if (p_nb>0) then
        nz = this%nnz
        call MPI_Send( dbuff(1,1), nz*width, MPI_Double_precision, &
            0, 0, MPI_Comm_World, ierr )
      endif
    enddo
    call de_alloc( dbuff, "dbuff", "t_ncdfSparse" )
  endif
#	endif /* MPI */
  end subroutine write_d2D

#endif /* NEW_NCDF>2 */

  subroutine cdf_w_Sp(ncdf, sp, dit, gncol)

    type(hNCDF), intent(inout) :: ncdf
    type(Sparsity), intent(inout) :: sp
    type(OrbitalDistribution), intent(in), optional :: dit
    integer, intent(inout), target, optional :: gncol(:)

    ! Local variables for reading the values
    integer, pointer :: lgncol(:) => null()
    integer, pointer :: ncol(:), l_col(:) => null()

    logical :: ldit
    integer :: lno, no, n_nzs
#ifdef MPI
    integer :: MPIerror
#endif
    integer :: n_nnzs, g_nzs

    ! Write the sparsity to the file...
    call attach(sp,nrows=lno, nrows_g=no, n_col=ncol,list_col=l_col,nnzs=n_nzs)

    ldit = present(dit)
    if ( ldit ) ldit = lno /= no

    if ( ldit ) then
#ifdef MPI
      call MPI_Reduce(n_nzs,g_nzs,1,MPI_Integer, &
          MPI_Sum, 0, MPI_Comm_World, MPIerror)

      if ( present(gncol) ) then
        lgncol => gncol
      else
        allocate(lgncol(no))
        lgncol(1) = -1
      end if
#else
      call die('Error in code, non-full contained sp')
#endif
    else
      ! serial (no distribution)
      lgncol => ncol
      g_nzs = n_nzs
    end if

    ! Check dimension nnzs
    call ncdf_inq_dim(ncdf,'nnzs',len=n_nnzs)
    if ( Node == 0 .and. n_nnzs /= g_nzs ) then
      call die('Number of non-zero elements is not equivalent.')
    end if

#ifdef MPI
    if ( ldit ) then
      if ( parallel_io(ncdf) ) then
        ! Write out in parallel
        call write_parallel()
      else
        call write_parallel_single()
      end if

      if ( .not. present(gncol) ) deallocate(lgncol)
    else
      call write_sequential()
    end if
#else
    call write_sequential()
#endif

  contains

#ifdef MPI
    subroutine write_parallel()

      integer :: BNode
      integer :: gio, io, ind, gind, n, i, nb

      if ( lgncol(1) < 0 ) then
        ! this means broadcast lgncol on all nodes
        call Node_Sp_gncol(-1, sp, dit, no, lgncol)
      end if

      ! we write it using MPI
      call ncdf_par_access(ncdf,name='n_col',access=NF90_COLLECTIVE)

      ! Calculate indices for this node
      io = no / Nodes
      ind = io * Node + 1
      if ( Node == Nodes - 1 ) then
        ! correct the last node to write the rest, probably not
        ! balanced, but shouldn't be too problematic for
        ! very large matrices
        io = io + mod(no, Nodes)
      end if

      ! Write in parallel (faster than in below loop)
      call ncdf_put_var(ncdf,'n_col',lgncol(ind:ind+io-1), &
          start=(/ind/), count=(/io/))

      call ncdf_par_access(ncdf,name='list_col',access=NF90_COLLECTIVE)

      ! Determine the number of blocks
      gio = 1
      gind = 1
      ind = 0
      nb = 1
      do while ( gio <= no )

        BNode = node_handling_element(dit,gio)

        ! Get number of consecutive orbitals
        n = count_consecutive(dit,no,gio)

        ! Number of consecutive elements
        i = sum(lgncol(gio:gio-1+n))

        if ( BNode == Node ) then
          ! Count number of actual writes
          nb = nb + 1
          call ncdf_put_var(ncdf,'list_col',l_col(ind+1:ind+i), &
              count=(/i/),start=(/gind/))
          ind = ind + i
        end if

        gind = gind + i
        gio = gio + n

      end do

      ind = count_blocks(dit, no)
      do i = nb, ind ! only runned if some iterations are missed
        ! Fake access
        call ncdf_put_var(ncdf,'list_col',l_col(1:1), &
            start=(/1/), count=(/0/) )
      end do

    end subroutine write_parallel
    
    subroutine write_parallel_single()

      integer :: BNode, MPIerror, MPIstatus(MPI_STATUS_SIZE)
      integer :: gio, io, max_n, ind, gind, n, i, nb
      integer, allocatable :: ibuf(:)

      if ( lgncol(1) < 0 ) then
        ! this means gather lgncol on IO-node
        call Node_Sp_gncol(0, sp, dit, no, lgncol)
      end if

      ! Write gncol
      call ncdf_put_var(ncdf,'n_col',lgncol)

      nb = count_blocks(dit,no)

      if ( Node == 0 ) then
        max_n = max_consecutive_sum(dit,no,lgncol)
        allocate(ibuf(max_n))
      else
        allocate(ibuf(nb))
      end if

      ! Loop size
      gio = 1
      gind = 1
      ind = 0
      nb = 0
      do while ( gio <= no )

        BNode = node_handling_element(dit,gio)

        ! Get number of consecutive orbitals
        n = count_consecutive(dit,no,gio)

        if ( Node == BNode ) then
          
          io = index_global_to_local(dit,gio,Node)
          i = sum(ncol(io:io-1+n))

          if ( Node == 0 ) then
            call ncdf_put_var(ncdf,'list_col',l_col(ind+1:ind+i), &
                count=(/i/),start=(/gind/))
            gind = gind + i
          else
            nb = nb + 1
            call MPI_ISend(l_col(ind+1), i, MPI_Integer, &
                0, gio, MPI_Comm_World, ibuf(nb), MPIerror)
          end if
          ind = ind + i
          
        else if ( Node == 0 ) then
          call MPI_Recv(ibuf(1), max_n, MPI_Integer, &
              BNode, gio, MPI_Comm_World, MPIstatus, MPIerror)
          if ( MPIerror /= MPI_Success ) &
              call die('Error in code: cdf_w_Sp')
          
          call MPI_Get_Count(MPIstatus, MPI_Integer, i, MPIerror)
          call ncdf_put_var(ncdf,'list_col',ibuf(1:i), &
              count=(/i/),start=(/gind/))
          gind = gind + i
        end if

        gio = gio + n
      end do

      if ( Node /= 0 .and. nb > 0 ) then
        call MPI_WaitAll(nb, ibuf, MPI_STATUSES_IGNORE, MPIerror)
      end if

      deallocate(ibuf)

    end subroutine write_parallel_single
#endif

    subroutine write_sequential()
      call ncdf_put_var(ncdf,'n_col',ncol)
      call ncdf_put_var(ncdf,'list_col',l_col)
    end subroutine write_sequential
    
  end subroutine cdf_w_Sp


  subroutine cdf_r_d1D(ncdf,vname,sp, dSp1D, tag, dit, Bcast, gncol)

    use class_dSpData1D

    type(hNCDF), intent(inout) :: ncdf
    character(len=*), intent(in) :: vname
    type(Sparsity), intent(inout) :: sp
    type(dSpData1D), intent(inout) :: dSp1D
    ! The tag of the data format
    character(len=*), intent(in) :: tag

    ! distribution if needed to be b-cast in a non-global
    ! fashion
    type(OrbitalDistribution), intent(in), optional :: dit
    ! Bcast the values?
    logical, intent(in), optional :: Bcast
    integer, intent(inout), target, optional :: gncol(:)

    ! Local variables for reading the values
    type(OrbitalDistribution) :: fdit
    real(dp), pointer :: a(:) => null()
    integer, pointer :: lgncol(:) => null(), ncol(:)

    integer :: lno, no, n_nzs
    logical :: ldit, lBcast, lIO
#ifdef MPI
    real(dp), allocatable :: buf(:)
    integer, allocatable :: ibuf(:)
    integer :: MPIerror, BNode
    integer :: max_n
#endif

    integer :: gind, gio, ind, n, nb, i, io

    ldit = present(dit)
    lBcast = .false.
    if ( present(Bcast) ) lBcast = Bcast
    lIO = .not. (lBcast .or. ldit)
    if ( .not. lIO ) lIO = (Node == 0)

    call attach(sp,nrows=lno,nrows_g=no, n_col=ncol,nnzs=n_nzs)
    if ( ldit ) ldit = lno /= no

    if ( ldit ) then
      call newdSpData1D(sp,dit,dSp1D,name=trim(tag))

#ifdef MPI
      if ( present(gncol) ) then
        lgncol => gncol
      else
        allocate(lgncol(no))
        lgncol(1) = -1
      end if
      if ( lgncol(1) < 0 ) then
        call Node_Sp_gncol(0, sp, dit, no, lgncol)
      end if
#else
      call die('Error in distribution, cdf_r_d1D')
#endif

    else
      ! Create the Fake distribution
#ifdef MPI
      call newDistribution(no,MPI_Comm_Self,fdit,name='Fake dist')
#else
      call newDistribution(no,-1           ,fdit,name='Fake dist')
#endif
      call newdSpData1D(sp,fdit,dSp1D,name=trim(tag))
      ! Clean up the distribution again
      call delete(fdit)
    end if

    a => val(dSp1D)

    if ( ldit ) then
#ifdef MPI

      nb = count_blocks(dit,no)

      ! Allocate the maximum number of entries
      if ( Node == 0 ) then
        max_n = max_consecutive_sum(dit,no,lgncol)
        allocate(buf(max_n))
      else
        allocate(ibuf(nb))
      end if

      ! Loop size
      gio = 1 
      gind = 1
      ind = 0
      nb = 0
      do while ( gio <= no )

        BNode = node_handling_element(dit,gio)

        ! Get number of consecutive orbitals
        n = count_consecutive(dit,no,gio)

        if ( BNode == Node ) then

          ! Get the local orbital
          io = index_global_to_local(dit,gio,Node)

          if ( Node == 0 ) then

            ! Count number of elements
            i = sum(ncol(io:io-1+n))
            call ncdf_get_var(ncdf, vname, a(ind+1:ind+i), &
                start=(/gind/), count=(/i/))
            ind = ind + i
            gind = gind + i

          else

            ! count the number of received entities
            i = sum(ncol(io:io-1+n))
            nb = nb + 1
            call MPI_IRecv(a(ind+1), i, MPI_Double_Precision, &
                0, gio, MPI_Comm_World, ibuf(nb), MPIerror)
            ind = ind + i

          end if

        else if ( Node == 0 ) then

          i = sum(lgncol(gio:gio-1+n))
          call ncdf_get_var(ncdf, vname, buf(1:i), &
              start=(/gind/), count=(/i/))

          call MPI_Send(buf(1), i, MPI_Double_Precision, &
              BNode, gio, MPI_Comm_World, MPIerror)

          gind = gind + i

        end if

        gio = gio + n

      end do

      if ( .not. present(gncol) ) deallocate(lgncol)
      if ( Node == 0 ) then
        deallocate(buf)
      else
        if ( nb > 0 ) then
          call MPI_WaitAll(nb, ibuf, MPI_STATUSES_IGNORE, MPIerror)
        end if
        deallocate(ibuf)
      end if
#else
      call die('Error in distribution for, cdf_r_d1D')
#endif
    else if ( lIO ) then

      i = sum(ncol(1:no))
      call ncdf_get_var(ncdf, vname, a(:), count=(/i/))

    end if

    ! If a distribution is present, then do something
#ifdef MPI
    if ( lBcast ) then

      call MPI_Bcast(a(1),n_nzs,MPI_Double_Precision, &
          0, MPI_Comm_World, MPIError)

    end if
#endif

  end subroutine cdf_r_d1D
  
  subroutine cdf_w_d1D(ncdf, vname, dSp1D, gncol)

    use class_dSpData1D

    type(hNCDF), intent(inout) :: ncdf
    character(len=*), intent(in) :: vname
    type(dSpData1D), intent(inout) :: dSp1D
    integer, intent(inout), target, optional :: gncol(:)

    ! Local variables
    type(OrbitalDistribution), pointer :: dit
    type(Sparsity), pointer :: sp
    integer, pointer :: ncol(:), l_col(:), lgncol(:)
    integer :: lno, no, n_nzs
    real(dp), pointer :: a(:)

    logical :: ldit

    dit => dist(dSp1D)
    sp => spar(dSp1D)
    ! Write the sparsity to the file...
    call attach(sp,nrows=lno, nrows_g=no, n_col=ncol,list_col=l_col,nnzs=n_nzs)

    ! If they are different we should 
    ! use the distribution setting
    ldit = lno /= no

    if ( ldit ) then

      if ( present(gncol) ) then
        lgncol => gncol
      else
        allocate(lgncol(no))
        lgncol(1) = -1
      end if

    else
      ! serial (no distribution)
      lgncol => ncol
    end if
    
    a => val(dSp1D)

#ifdef MPI
    if ( ldit ) then
      if ( parallel_io(ncdf) ) then
        ! Write out in parallel
        call write_parallel()
      else
        call write_parallel_single()
      end if
      
      if ( .not. present(gncol) ) deallocate(lgncol)
    else
      call write_sequential()
    end if
#else
    call write_sequential()
#endif

  contains

#ifdef MPI
    subroutine write_parallel()

      integer :: BNode
      integer :: gio, ind, gind, n, i, nb

      if ( lgncol(1) < 0 ) then
        ! this means broadcast lgncol on all nodes
        call Node_Sp_gncol(-1, sp, dit, no, lgncol)
      end if

      call ncdf_par_access(ncdf, name=vname, access=NF90_COLLECTIVE)

      ! Loop blocks
      gio = 1
      gind = 1
      ind = 0
      nb = 1
      do while ( gio <= no )

        BNode = node_handling_element(dit,gio)

        ! Get number of consecutive orbitals
        n = count_consecutive(dit,no,gio)

        ! Number of consecutive elements
        i = sum(lgncol(gio:gio-1+n))

        if ( BNode == Node ) then
          ! Count number of actual writes
          nb = nb + 1
          call ncdf_put_var(ncdf, vname, a(ind+1:ind+i), &
              count=(/i/), start=(/gind/))
          ind = ind + i
        end if

        gind = gind + i
        gio = gio + n

      end do

      ind = count_blocks(dit, no)
      do i = nb, ind ! only runned if some iterations are missed
        ! Fake access
        call ncdf_put_var(ncdf, vname, a(1:1), start=(/1/), count=(/0/))
      end do

    end subroutine write_parallel
    
    subroutine write_parallel_single()

      integer :: BNode, MPIerror, MPIstatus(MPI_STATUS_SIZE)
      integer :: gio, io, max_n, ind, gind, n, i, nb
      real(dp), allocatable :: buf(:)
      integer, allocatable :: ibuf(:)

      if ( lgncol(1) < 0 ) then
        ! this means gather lgncol on IO-node
        call Node_Sp_gncol(0, sp, dit, no, lgncol)
      end if

      nb = count_blocks(dit,no)

      if ( Node == 0 ) then
        max_n = max_consecutive_sum(dit,no,lgncol)
        allocate(buf(max_n))
      else
        allocate(ibuf(nb))
      end if

      ! Loop size
      gio = 1
      gind = 1
      ind = 0
      nb = 0
      do while ( gio <= no )

        BNode = node_handling_element(dit,gio)

        ! Get number of consecutive orbitals
        n = count_consecutive(dit,no,gio)

        if ( Node == BNode ) then
          
          io = index_global_to_local(dit,gio,Node)
          i = sum(ncol(io:io-1+n))

          if ( Node == 0 ) then
            call ncdf_put_var(ncdf, vname, a(ind+1:ind+i), &
                count=(/i/), start=(/gind/))
            gind = gind + i
          else
            nb = nb + 1
            call MPI_ISend(a(ind+1), i, MPI_Double_Precision, &
                0, gio, MPI_Comm_World, ibuf(nb), MPIerror)
          end if
          ind = ind + i
          
        else if ( Node == 0 ) then
          call MPI_Recv(buf(1), max_n, MPI_Double_Precision, &
              BNode, gio, MPI_Comm_World, MPIstatus, MPIerror)
          if ( MPIerror /= MPI_Success ) &
              call die('Error in code: cdf_w_d1D')
          
          call MPI_Get_Count(MPIstatus, MPI_Double_Precision, i, MPIerror)
          call ncdf_put_var(ncdf, vname, buf(1:i), &
              count=(/i/), start=(/gind/))
          gind = gind + i
        end if

        gio = gio + n
      end do

      if ( Node == 0 ) then
        deallocate(buf)
      else
        if ( nb > 0 ) then
          call MPI_WaitAll(nb, ibuf, MPI_STATUSES_IGNORE, MPIerror)
        end if
        deallocate(ibuf)
      end if

    end subroutine write_parallel_single
#endif

    subroutine write_sequential()
      call ncdf_put_var(ncdf, vname, a)
    end subroutine write_sequential
    
  end subroutine cdf_w_d1D


  subroutine cdf_r_d2D(ncdf, vname, sp, dSp2D, dim2, tag, &
      sparsity_dim, dit, Bcast, gncol)

    use class_dSpData2D

    ! File handle
    type(hNCDF), intent(inout) :: ncdf
    ! Variable name
    character(len=*), intent(in) :: vname
    ! Sparsity pattern
    type(Sparsity), intent(inout) :: sp
    ! Data array with sparsity pattern
    type(dSpData2D), intent(inout) :: dSp2D
    ! The non-sparse dimension
    integer, intent(in) :: dim2
    ! The tag of the sparsity pattern
    character(len=*), intent(in) :: tag
    ! This denotes the sparsity dimension (either 1 or 2) 1=default
    integer, intent(in), optional :: sparsity_dim
    ! distribution if needed to be b-cast in a non-global
    ! fashion
    type(OrbitalDistribution), intent(in), optional :: dit
    ! Bcast the values?
    logical, intent(in), optional :: Bcast
    integer, intent(inout), target, optional :: gncol(:)

    ! Local variables for reading the values
    type(OrbitalDistribution) :: fdit
    real(dp), pointer :: a(:,:) => null()
    integer, pointer :: lgncol(:) => null(), ncol(:)

    integer :: lno, no, n_nzs
    integer :: sp_dim
    logical :: ldit, lBcast, lIO
#ifdef MPI
    integer :: MPIerror
#endif

    ldit = present(dit)
    lBcast = .false.
    if ( present(Bcast) ) lBcast = Bcast
    lIO = .not. (lBcast .or. ldit)
    if ( .not. lIO ) lIO = (Node == 0)

    sp_dim = 1
    if ( present(sparsity_dim) ) sp_dim = sparsity_dim

    call attach(sp,nrows=lno, nrows_g=no, n_col=ncol,nnzs=n_nzs)
    if ( ldit ) ldit = lno /= no

    if ( ldit ) then
      call newdSpData2D(sp,dim2,dit,dSp2D,name=trim(tag), &
          sparsity_dim=sp_dim)

#ifdef MPI
      if ( present(gncol) ) then
        lgncol => gncol
      else
        allocate(lgncol(no))
        lgncol(1) = -1
      end if
      if ( lgncol(1) < 0 ) then
        call Node_Sp_gncol(0, sp, dit, no, lgncol)
      end if
#else
      call die('Error in distribution, cdf_r_d2D')
#endif

    else
      ! Create the Fake distribution
#ifdef MPI
      call newDistribution(no,MPI_Comm_Self,fdit,name='Fake dist')
#else
      call newDistribution(no,-1           ,fdit,name='Fake dist')
#endif
      call newdSpData2D(sp,dim2,fdit,dSp2D,name=trim(tag), &
          sparsity_dim=sp_dim)
    end if

    a => val(dSp2D)

#ifdef MPI
    if ( ldit ) then
      if ( sp_dim == 1 ) then
        call read_sp_dim1()
      else
        call read_sp_dim2()
      end if

      if ( .not. present(gncol) ) deallocate(lgncol)
    else
      ! this is both lio and serial read
      call ncdf_get_var(ncdf, vname, a)
    end if
#else
    call ncdf_get_var(ncdf, vname, a)
#endif

    ! If a distribution is present, then do something
#ifdef MPI
    if ( lBcast ) then
      
      call MPI_Bcast(a(1,1), dim2*n_nzs, MPI_Double_Precision, &
          0, MPI_Comm_World, MPIError)

    end if
#endif

#ifdef MPI
  contains

    subroutine read_sp_dim1()

      integer :: io, ind, n, i, nb, s, gio, gind, max_n
      real(dp), allocatable :: buf(:)
      integer, allocatable :: ibuf(:)
      integer :: MPIerror, BNode

      nb = count_blocks(dit,no)

      ! Allocate maximum number of entries
      if ( Node == 0 ) then
        max_n = max_consecutive_sum(dit,no,lgncol)
        allocate(buf(max_n))
      else
        allocate(ibuf(nb))
      end if

      do s = 1 , dim2

        gio = 1 
        gind = 1
        ind = 0
        nb = 0
        do while ( gio <= no )

          BNode = node_handling_element(dit,gio)

          ! Get number of consecutive orbitals
          n = count_consecutive(dit,no,gio)

          if ( BNode == Node ) then

            ! Get the local orbital
            io = index_global_to_local(dit,gio,Node)

            i = sum(ncol(io:io-1+n))
            if ( Node == 0 ) then

              call ncdf_get_var(ncdf, vname, a(ind+1:ind+i,s), &
                  start=(/gind,s/), count=(/i/))
              gind = gind + i

            else

              ! count the number of received entities
              nb = nb + 1
              call MPI_IRecv(a(ind+1,s), i, MPI_Double_Precision, &
                  0, gio, MPI_Comm_World, ibuf(nb), MPIerror)

            end if
            ind = ind + i

          else if ( Node == 0 ) then

            i = sum(lgncol(gio:gio-1+n))
            call ncdf_get_var(ncdf, vname, buf(1:i), &
                start=(/gind,s/), count=(/i/))
            call MPI_Send(buf(1), i, MPI_Double_Precision, &
                BNode, gio, MPI_Comm_World, MPIerror)
            gind = gind + i
          end if
          
          gio = gio + n
        end do

        if ( Node /= 0 .and. nb > 0 ) then
          call MPI_WaitAll(nb, ibuf, MPI_STATUSES_IGNORE, MPIerror)
        end if

      end do

      if ( Node == 0 ) then
        deallocate(buf)
      else
        deallocate(ibuf)
      end if

    end subroutine read_sp_dim1

    subroutine read_sp_dim2()
      
      integer :: io, ind, n, i, nb, gio, gind, max_n
      real(dp), allocatable :: buf(:)
      integer, allocatable :: ibuf(:)
      integer :: MPIerror, BNode

      nb = count_blocks(dit,no)

      ! Allocate maximum number of entries
      if ( Node == 0 ) then
        max_n = max_consecutive_sum(dit,no,lgncol) * dim2
        allocate(buf(max_n))
      else
        allocate(ibuf(nb))
      end if

      ! Read sparse blocks and distribute
      gind = 1
      ind = 0
      gio = 1 
      nb = 0
      do while ( gio <= no )

        BNode = node_handling_element(dit,gio)

        ! Get number of consecutive orbitals
        n = count_consecutive(dit,no,gio)
        if ( BNode == Node ) then

          ! Get the local orbital
          io = index_global_to_local(dit,gio,Node)

          i = sum(ncol(io:io-1+n))
          if ( Node == 0 ) then
            call ncdf_get_var(ncdf, vname, a(1:dim2,ind+1:ind+i), &
                start=(/1,gind/), count=(/dim2,i/) )
            gind = gind + i

          else

            ! count the number of received entities
            nb = nb + 1
            call MPI_IRecv(a(1,ind+1), dim2*i, MPI_Double_Precision, &
                0, gio, MPI_Comm_World, ibuf(nb), MPIerror)
          end if
          ind = ind + i

        else if ( Node == 0 ) then

          i = sum(lgncol(gio:gio-1+n))
          call ncdf_get_var(ncdf, vname, buf(1:i), &
              start=(/1,gind/), count=(/dim2,i/))
          call MPI_Send(buf(1), i, MPI_Double_Precision, &
                BNode, gio, MPI_Comm_World, MPIerror)
          gind = gind + i

        end if

        gio = gio + n
      end do

      if ( Node == 0 ) then
        deallocate(buf)
      else
        if ( nb > 0 ) then
          call MPI_WaitAll(nb, ibuf, MPI_STATUSES_IGNORE, MPIerror)
        end if
        deallocate(ibuf)
      end if

    end subroutine read_sp_dim2

#endif
    
  end subroutine cdf_r_d2D

  subroutine cdf_w_d2D(ncdf, vname, dSp2D, gncol)

    use class_dSpData2D

    type(hNCDF), intent(inout) :: ncdf
    character(len=*), intent(in) :: vname
    type(dSpData2D), intent(inout) :: dSp2D
    integer, intent(inout), target, optional :: gncol(:)

    ! Local variables
    type(OrbitalDistribution), pointer :: dit
    type(Sparsity), pointer :: sp
    integer, pointer :: ncol(:), l_col(:), lgncol(:)
    integer :: lno, no, n_nzs
    integer :: dim2, sp_dim
    real(dp), pointer :: a(:,:)

    logical :: ldit

    dit => dist(dSp2D)
    sp => spar(dSp2D)
    ! Write the sparsity to the file...
    call attach(sp,nrows=lno, nrows_g=no, n_col=ncol,list_col=l_col,nnzs=n_nzs)

    ! If they are different we should 
    ! use the distribution setting
    ldit = lno /= no

    if ( ldit ) then

      if ( present(gncol) ) then
        lgncol => gncol
      else
        allocate(lgncol(no))
        lgncol(1) = -1
      end if

    else
      ! serial (no distribution)
      lgncol => ncol
    end if

    a => val(dSp2D)
    sp_dim = spar_dim(dSp2D)
    if ( sp_dim == 1 ) then
      dim2 = size(a, dim=2)
    else
      dim2 = size(a, dim=1)
    end if

#ifdef MPI
    if ( ldit ) then
      if ( parallel_io(ncdf) ) then
        ! Write out in parallel
        call write_parallel()
      else
        call write_parallel_single()
      end if

      if ( .not. present(gncol) ) deallocate(lgncol)
    else
      call write_sequential()
    end if
#else
    call write_sequential()
#endif

  contains

#ifdef MPI
    subroutine write_parallel()
      integer :: BNode
      integer :: gio, ind, gind, n, i, nb

      if ( lgncol(1) < 0 ) then
        ! this means broadcast lgncol on all nodes
        call Node_Sp_gncol(-1, sp, dit, no, lgncol)
      end if

      call ncdf_par_access(ncdf, name=vname, access=NF90_COLLECTIVE)

      gio = 1
      gind = 1
      ind = 0
      nb = 1
      do while ( gio <= no )

        BNode = node_handling_element(dit,gio)

        ! Get number of consecutive orbitals
        n = count_consecutive(dit,no,gio)

        ! Number of consecutive elements
        i = sum(lgncol(gio:gio-1+n))

        if ( BNode == Node ) then
          ! Count number of actual writes
          nb = nb + 1
          if ( sp_dim == 1 ) then
            ! this will make a temporary copy... :(
            call ncdf_put_var(ncdf, vname, a(ind+1:ind+i,:), &
                count=(/i,dim2/), start=(/gind,1/))
          else
            call ncdf_put_var(ncdf, vname, a(:,ind+1:ind+i), &
                count=(/dim2,i/), start=(/1,gind/))
          end if
          ind = ind + i
        end if

        gind = gind + i
        gio = gio + n

      end do

      ind = count_blocks(dit, no)
      do i = nb, ind ! only runned if some iterations are missed
        ! Fake access
        call ncdf_put_var(ncdf, vname, a(1:1,1:1), start=(/1,1/), count=(/0,0/))
      end do

    end subroutine write_parallel

    subroutine write_parallel_single()
      if ( lgncol(1) < 0 ) then
        ! this means gather lgncol on IO-node
        call Node_Sp_gncol(0, sp, dit, no, lgncol)
      end if

      if ( sp_dim == 1 ) then
        call write_parallel_single_sp_dim1()
      else
        call write_parallel_single_sp_dim2()
      end if

    end subroutine write_parallel_single

    subroutine write_parallel_single_sp_dim1()
      integer :: BNode, MPIstatus(MPI_STATUS_SIZE), MPIerror
      integer :: gio, io, max_n, ind, gind, n, i, nb, is
      real(dp), allocatable :: buf(:)
      integer, allocatable :: ibuf(:)

      nb = count_blocks(dit,no)

      if ( Node == 0 ) then
        max_n = max_consecutive_sum(dit,no,lgncol)
        allocate(buf(max_n))
      else
        allocate(ibuf(nb))
      end if

      ! Loop size
      do is = 1, dim2
        
        gio = 1
        gind = 1
        ind = 0
        nb = 0
        do while ( gio <= no )

          BNode = node_handling_element(dit,gio)

          ! Get number of consecutive orbitals
          n = count_consecutive(dit,no,gio)

          if ( Node == BNode ) then

            io = index_global_to_local(dit,gio,Node)
            i = sum(ncol(io:io-1+n))

            if ( Node == 0 ) then
              call ncdf_put_var(ncdf, vname, a(ind+1:ind+i,is), &
                  count=(/i,1/), start=(/gind,is/))
              gind = gind + i
            else
              nb = nb + 1
              call MPI_ISend(a(ind+1,is), i, MPI_Double_Precision, &
                  0, gio, MPI_Comm_World, ibuf(nb), MPIerror)
            end if
            ind = ind + i
            
          else if ( Node == 0 ) then
            call MPI_Recv(buf(1), max_n, MPI_Double_Precision, &
                BNode, gio, MPI_Comm_World, MPIstatus, MPIerror)
            if ( MPIerror /= MPI_Success ) &
                call die('Error in code: cdf_w_d2D')
          
            call MPI_Get_Count(MPIstatus, MPI_Double_Precision, i, MPIerror)
            call ncdf_put_var(ncdf, vname, buf(1:i), &
                count=(/i,1/), start=(/gind,is/))
            gind = gind + i
          end if
          
          gio = gio + n
        end do

        if ( Node /= 0 .and. nb > 0 ) then
          call MPI_WaitAll(nb, ibuf, MPI_STATUSES_IGNORE, MPIerror)
        end if

      end do

      if ( Node == 0 ) then
        deallocate(buf)
      else
        deallocate(ibuf)
      end if

    end subroutine write_parallel_single_sp_dim1

    subroutine write_parallel_single_sp_dim2()
      integer :: BNode, MPIstatus(MPI_STATUS_SIZE), MPIerror
      integer :: gio, io, max_n, ind, gind, n, i, nb
      real(dp), allocatable :: buf(:)
      integer, allocatable :: ibuf(:)

      nb = count_blocks(dit,no)

      if ( Node == 0 ) then
        max_n = max_consecutive_sum(dit,no,lgncol) * dim2
        allocate(buf(max_n))
      else
        allocate(ibuf(nb))
      end if

      gio = 1
      gind = 1
      ind = 0
      nb = 0
      do while ( gio <= no )

        BNode = node_handling_element(dit,gio)

        ! Get number of consecutive orbitals
        n = count_consecutive(dit,no,gio)
        
        if ( Node == BNode ) then

          io = index_global_to_local(dit,gio,Node)
          i = sum(ncol(io:io-1+n))

          if ( Node == 0 ) then
            call ncdf_put_var(ncdf, vname, a(:,ind+1:ind+i), &
                count=(/dim2,i/), start=(/1,gind/))
            gind = gind + i
          else
            nb = nb + 1
            call MPI_ISend(a(1,ind+1), i*dim2, MPI_Double_Precision, &
                0, gio, MPI_Comm_World, ibuf(nb), MPIerror)
          end if
          ind = ind + i

        else if ( Node == 0 ) then
          call MPI_Recv(buf(1), max_n, MPI_Double_Precision, &
              BNode, gio, MPI_Comm_World, MPIstatus, MPIerror)
          if ( MPIerror /= MPI_Success ) &
              call die('Error in code: cdf_w_d2D')

          call MPI_Get_Count(MPIstatus, MPI_Double_Precision, i, MPIerror)
          i=i/dim2
          call ncdf_put_var(ncdf, vname, buf(1:i), &
              count=(/dim2,i/), start=(/1,gind/))
          gind = gind + i
        end if
        
        gio = gio + n
      end do

      if ( Node == 0 ) then
        deallocate(buf)
      else
        if ( nb > 0 ) then
          call MPI_WaitAll(nb, ibuf, MPI_STATUSES_IGNORE, MPIerror)
        end if
        deallocate(ibuf)
      end if

    end subroutine write_parallel_single_sp_dim2
#endif
    
    subroutine write_sequential()
      call ncdf_put_var(ncdf, vname, a)
    end subroutine write_sequential

  end subroutine cdf_w_d2D
  
  subroutine cdf_w_grid(ncdf,name,nmeshl,grid,idx)

    type(hNCDF), intent(inout) :: ncdf
    character(len=*), intent(in) :: name
    integer, intent(in) :: nmeshl(3)
    real(grid_p), intent(in) :: grid(:)
    integer, intent(in), optional :: idx

#ifdef MPI
    integer :: MPIstat(MPI_STATUS_SIZE)
    integer :: MPIerror, mnpt
    integer :: lb(3), nel(3), iN, inpt
    real(grid_p), allocatable :: gb(:)
#endif

#ifdef MPI

    if ( parallel_io(ncdf) ) then

      ! Ensure collective writing
      call ncdf_par_access(ncdf,name=name,access=NF90_COLLECTIVE)

      lb(:)  = distr%box(1,:,Node)
      nel(:) = distr%box(2,:,Node) - lb(:) + 1
      if ( .not. all(nel == nmeshl) ) then
        call die('cdf_r_grid: cannot assert the grid-size from the &
            &stored grid and the denoted size')
      end if

      if ( present(idx) ) then
        call ncdf_put_var(ncdf,name,grid, &
            start=(/lb(1),lb(2),lb(3),idx/), &
            count=(/nel(1),nel(2),nel(3),1/) )          
      else
        call ncdf_put_var(ncdf,name,grid, start=lb, count=nel )
      end if

    else

      mnpt = 0
      do iN = 0 , Nodes - 1
        nel(:) = distr%box(2,:,iN) - distr%box(1,:,iN) + 1
        mnpt = max(mnpt,product(nel))
      end do

      ! The main node can safely write the data...
      if ( Node == 0 ) then

        allocate(gb(mnpt))

        ! First save it's own data
        lb(:) = distr%box(1,:,0) 
        nel(:) = distr%box(2,:,0) - lb(:) + 1
        if ( present(idx) ) then
          call ncdf_put_var(ncdf,name,grid, &
              start=(/lb(1),lb(2),lb(3),idx/), &
              count=(/nel(1),nel(2),nel(3),1/) )
        else
          call ncdf_put_var(ncdf,name,grid, &
              start=lb, count=nel )
        end if

        ! Loop on the remaining nodes
        do iN = 1 , Nodes - 1

          ! we retrieve data from the iN'th node.
          call MPI_Recv(gb,mnpt,MPI_grid_real,iN,iN, &
              MPI_Comm_World,MPIstat, MPIerror)

          ! Just make sure we only pass the correct size
          call MPI_Get_Count(MPIstat, MPI_Grid_Real, inpt, MPIerror)

          lb(:) = distr%box(1,:,iN) 
          nel(:) = distr%box(2,:,iN) - lb(:) + 1
          if ( inpt /= product(nel) ) then
            call die('Error when receiving the distributed &
                &grid data for writing to the NetCDF file.')
          end if

          if ( present(idx) ) then
            call ncdf_put_var(ncdf,name,gb(1:inpt), &
                start=(/lb(1),lb(2),lb(3),idx/), &
                count=(/nel(1),nel(2),nel(3),1/) )
          else
            call ncdf_put_var(ncdf,name,gb(1:inpt), &
                start=lb, count=nel )
          end if

        end do

        deallocate(gb)
      else
        mnpt = product(nmeshl)
        call MPI_Send(grid,mnpt,MPI_grid_real,0, &
            Node,MPI_Comm_World,MPIerror)
      end if

    end if

#else
    if ( present(idx) ) then
      call ncdf_put_var(ncdf,name,grid,start=(/1,1,1,idx/), &
          count=nmeshl)
    else
      call ncdf_put_var(ncdf,name,grid, &
          count=nmeshl)
    end if
#endif

  end subroutine cdf_w_grid

  subroutine cdf_r_grid(ncdf,name,nmeshl,grid,idx)

    type(hNCDF), intent(inout) :: ncdf
    character(len=*), intent(in) :: name
    integer, intent(in) :: nmeshl(3)
    real(grid_p), intent(inout), target :: grid(:)
    integer, intent(in), optional :: idx

#ifdef MPI
    integer :: MPIstat(MPI_STATUS_SIZE)
    integer :: MPIerror, mnpt
    integer :: lb(3), nel(3), iN, inpt
    real(grid_p), pointer :: gb(:) => null()
#endif

#ifdef MPI

    if ( parallel_io(ncdf) ) then

      ! Ensure collective writing
      call ncdf_par_access(ncdf,name=name,access=NF90_COLLECTIVE)

      lb(:)  = distr%box(1,:,Node)
      nel(:) = distr%box(2,:,Node) - lb(:) + 1
      if ( .not. all(nel == nmeshl) ) then
        call die('cdf_r_grid: cannot assert the grid-size from the &
            &stored grid and the denoted size')
      end if

      if ( present(idx) ) then
        call ncdf_get_var(ncdf,name,grid, &
            start=(/lb(1),lb(2),lb(3),idx/), &
            count=(/nel(1),nel(2),nel(3),1/) )          
      else
        call ncdf_get_var(ncdf,name,grid, start=lb, count=nel )
      end if

    else

      mnpt = 0
      do iN = 0 , Nodes - 1
        nel(:) = distr%box(2,:,iN) - distr%box(1,:,iN) + 1
        mnpt = max(mnpt,product(nel))
      end do

      ! The main node can safely read the data...
      if ( Node == 0 ) then

        if ( mnpt > product(nmeshl) ) then
          ! allocate, the IO node have a smaller
          ! space than the others
          allocate(gb(mnpt))
        else
          gb => grid
        end if

        ! Loop on the remaining nodes (we may possibly
        ! reuse the current grid array, and hence
        ! we first read the other nodes...)
        do iN = 1 , Nodes - 1

          lb(:) = distr%box(1,:,iN) 
          nel(:) = distr%box(2,:,iN) - lb(:) + 1

          ! get total number of points
          inpt = product(nel)

          if ( present(idx) ) then
            call ncdf_get_var(ncdf,name,gb(1:inpt), &
                start=(/lb(1),lb(2),lb(3),idx/), &
                count=(/nel(1),nel(2),nel(3),1/) )
          else
            call ncdf_get_var(ncdf,name,gb(1:inpt), &
                start=lb, count=nel )
          end if

          ! we send data to the iN'th node.
          call MPI_Send(gb,inpt,MPI_grid_real,iN,iN, &
              MPI_Comm_World,MPIerror)

        end do

        if ( mnpt > product(nmeshl) ) then
          deallocate(gb)
        end if

        ! Retrieve it's own data
        lb(:) = distr%box(1,:,0) 
        nel(:) = distr%box(2,:,0) - lb(:) + 1

        if ( present(idx) ) then
          call ncdf_get_var(ncdf,name,grid, &
              start=(/lb(1),lb(2),lb(3),idx/), &
              count=(/nel(1),nel(2),nel(3),1/) )
        else
          call ncdf_get_var(ncdf,name,grid, &
              start=lb, count=nel )
        end if

      else
        mnpt = product(nmeshl)
        call MPI_Recv(grid,mnpt,MPI_grid_real,0, &
            Node,MPI_Comm_World,MPIstat,MPIerror)
      end if

    end if

#else
    if ( present(idx) ) then
      call ncdf_get_var(ncdf,name,grid,start=(/1,1,1,idx/), &
          count=nmeshl)
    else
      call ncdf_get_var(ncdf,name,grid,start=(/1,1,1/), &
          count=nmeshl)
    end if
#endif

  end subroutine cdf_r_grid

#else
  subroutine dummy()
    !pass
  end subroutine dummy
#endif

end module ncdf_io_m
