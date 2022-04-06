! Module for easy writing of all information
! in a single netCDF file.

! We rely on the NetCDF4 library and can
! utilise parallel IO or standard IO

module ncdf_io_m

  use precision, only : sp, dp, grid_p, i8b
  use parallel, only : Node, Nodes

  use class_OrbitalDistribution
  use class_Sparsity
  use io_sparse_m, only: max_consecutive_sum, Node_Sp_gncol
  use io_sparse_m, only: count_consecutive, count_blocks

# ifdef NCDF_4
  use netcdf_ncdf, ncdf_parallel => parallel
# endif
# ifdef MPI
  use mpi_siesta
  use mpi, only : MPI_COMM_TYPE_SHARED, MPI_INFO_NULL
# endif
  use alloc

  implicit none

  private

# ifdef NCDF_4
  public :: cdf_init_mesh
  public :: cdf_r_grid
  public :: cdf_w_grid

  public :: cdf_w_sp
  public :: cdf_w_d1D, cdf_w_d2D

  public :: cdf_r_sp
  public :: cdf_r_d1D, cdf_r_d2D
# endif

  type :: t_ncdfSparse
    logical          :: IOnode             ! True if I am a writer
    logical          :: memset             ! True if memory is set
    integer          :: nwriters           ! Number of writers
    integer          :: nb                 ! Local number of blocks of rows
    integer          :: nnz                ! Number of Nzs of current distribution
    integer          :: bsize              ! Buffer size
    integer          :: g_no               ! Global number of orbitals
# if NEW_NCDF<=3
    integer          :: wroff              ! Writer offset
#endif
    integer, pointer :: p_nrows(:)         ! Number of rows of the new distribution
    integer, pointer :: p_xrows(:)         ! Offset of p_nrows
    integer, pointer :: p_nbloc(:)         ! Number of blocks (group of rows) of the new distribution
    integer, pointer :: p_xbloc(:)         ! Offset of  p_nbloc
    integer, pointer :: blonz(:)           ! Number of nonzeros of every block (original distribution)
    integer, pointer :: p_blonz(:)         ! Number of nonzeros of every block (new distribution)
    integer, pointer :: g_ncol(:)          ! Number of columns of every global orbital
    integer, pointer :: nsend(:)           ! Number of NZs to send to every process
    integer, pointer :: xsend(:)           ! Offset of nsend
    integer, pointer :: nrecv(:)           ! Number of NZs to receive from every process
    integer, pointer :: xrecv(:)           ! Offset of nrecv
    integer, pointer :: p_nnz(:)           ! Number of NZs of every process in the new distribution
# if NEW_NCDF==5
    integer, pointer :: sreq(:)            ! Buffer to store the MPI send request
    integer, pointer :: rreq(:)            ! Buffer to store the MPI receive request
# endif
    type(OrbitalDistribution), pointer :: dit  ! A pointer to the orbital distribution
    contains
    procedure :: init_Sp
    procedure :: write_Sp
    procedure :: write_d1D
    procedure :: write_d2D
    procedure :: delete_Sp
  end type t_ncdfSparse
  public :: t_ncdfSparse

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






  subroutine delete_Sp( this )
  implicit none
  class(t_ncdfSparse) :: this
  if (this%memset) then
    call de_alloc( this%p_nrows, "p_nrows", "t_ncdfSparse" )
    call de_alloc( this%p_xrows, "p_xrows", "t_ncdfSparse" )
    call de_alloc( this%p_nbloc, "p_nbloc", "t_ncdfSparse" )
    call de_alloc( this%p_xbloc, "p_xbloc", "t_ncdfSparse" )

    call de_alloc( this%blonz, "blonz", "t_ncdfSparse" )
    call de_alloc( this%p_blonz, "p_blonz", "t_ncdfSparse" )
# 	if NEW_NCDF==5
    call de_alloc( this%sreq, "sreq", "t_ncdfSparse" )
    call de_alloc( this%rreq, "rreq", "t_ncdfSparse" )
#   else
    call de_alloc( this%g_ncol, "g_ncol", "t_ncdfSparse" )
    call de_alloc( this%nsend, "nsend", "t_ncdfSparse" )
    call de_alloc( this%nrecv, "nrecv", "t_ncdfSparse" )
    call de_alloc( this%xsend, "xsend", "t_ncdfSparse" )
    call de_alloc( this%xrecv, "xrecv", "t_ncdfSparse" )
#   endif
#	  if NEW_NCDF==4 || NEW_NCDF==5
    call de_alloc( this%p_nnz, "p_nnz", "t_ncdfSparse" )
#   endif
  endif
  end subroutine delete_Sp

  subroutine init_Sp( this, sp, dit, max_dsize )
  implicit none
  class(t_ncdfSparse) :: this
  type(Sparsity), intent(inout) :: sp
  type(OrbitalDistribution), target :: dit
  integer :: max_dsize
  integer :: p, nw, gio, gnb, BNode, n, no, g_no, nnzs, gsize, ngroups
  integer :: lnb, lno, i, io, b, gb, COMM_SM, ierr
  integer, pointer :: ncol(:), ibuf(:)
  integer, pointer :: nsend(:), nrecv(:), xsend(:), xrecv(:)

  this%memset =.false.
  this%nwriters = 1
  this%IOnode = Node==0
  if (max_dsize==0) return
  if (Nodes==1) return
#ifdef MPI
# if NEW_NCDF==1
  ! Compute the number of nodes that will write to the CDF file
  ! One reductor for every host
  CALL MPI_Comm_SPLIT_TYPE( MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, COMM_SM, ierr )
  CALL MPI_Comm_size( COMM_SM, gsize, ierr )
  call MPI_Bcast( gsize, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )
  ! Number of writers
  this%nwriters = (Nodes+gsize-1)/gsize
  this%IOnode = mod(node,gsize) == 0
  ngroups = this%nwriters
# elif NEW_NCDF==2 || NEW_NCDF==3
  this%nwriters = Nodes
  this%IOnode = .true.
  gsize = 1
  ngroups = Nodes
#	elif NEW_NCDF==4
  this%nwriters = 1
  this%IOnode = Node==0
  gsize = 1
  ngroups = Nodes
#	elif NEW_NCDF==5
  ! Compute the number of nodes that will reduce the data
  ! One reductor for every host
  CALL MPI_Comm_SPLIT_TYPE( MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, COMM_SM, ierr )
  CALL MPI_Comm_size( COMM_SM, gsize, ierr )
  call MPI_Bcast( gsize, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )
  ! Number of writers
  ngroups = (Nodes+gsize-1)/gsize
  this%nwriters = 1
  this%IOnode = Node == 0
# endif /* NEW_NCDF==1,2,3,4,5 */
  ! Get information from Sparse matrix
  nullify(ncol)
  call attach( sp, nrows=no, nrows_g=g_no, n_col=ncol, nnzs=nnzs )
  this%dit => dit
  ! allocate memory
  nullify(this%p_nrows,this%p_xrows,this%p_nbloc,this%p_xbloc)
  call re_alloc( this%p_nrows, 1, Nodes, "p_nrows", "t_ncdfSparse" )
  call re_alloc( this%p_xrows, 1, Nodes, "p_xrows", "t_ncdfSparse" )
  call re_alloc( this%p_nbloc, 1, Nodes, "p_nbloc", "t_ncdfSparse" )
  call re_alloc( this%p_xbloc, 1, Nodes, "p_xbloc", "t_ncdfSparse" )
  ! Compute the number of rows and the number of blocks of rows
  ! that will handle every process
  p = 1
  nw = 1
  gio = 1
  gnb = 0
  this%p_nrows(p) = 0
  this%p_nbloc(p) = 0
  this%p_xrows(p) = 0
  this%p_xbloc(p) = 0
  this%nb = 0
  do while( gio<=g_no )
    BNode = node_handling_element(dit,gio)
    if (node==BNode) this%nb = this%nb + 1
    ! Get number of consecutive orbitals
    n = count_consecutive( dit, g_no, gio )
    this%p_nrows(p) = this%p_nrows(p) + n
    this%p_nbloc(p) = this%p_nbloc(p) + 1
    gio = gio + n
    gnb = gnb + 1
    if (gio>(nw*g_no)/ngroups) then
      do while (p<min(nw*gsize+1,Nodes))
        p = p+1
        this%p_nrows(p) = 0
        this%p_nbloc(p) = 0
        this%p_xrows(p) = gio-1
        this%p_xbloc(p) = gnb
      enddo
      nw = nw + 1
    endif
  enddo
# if NEW_NCDF==5
  ! Compute NZ of every block of local distribution
  lnb = this%p_nbloc(node+1)
  nullify(this%blonz,this%p_blonz,this%sreq,this%rreq)
  call re_alloc( this%blonz, 1, this%nb, "blonz", "t_ncdfSparse" )
  call re_alloc( this%p_blonz, 1, lnb, "p_blonz", "t_ncdfSparse" )
  call re_alloc( this%sreq, 1, this%nb, "sreq", "t_ncdfSparse" )
  call re_alloc( this%rreq, 1, lnb, "rreq", "t_ncdfSparse" )
  io = 1
  do b=1, this%nb
    gio = index_local_to_global( dit, io )
    n = count_consecutive( dit, g_no, gio )
    this%blonz(b) = sum(ncol(io:io+n-1))
    io = io+n
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
    gb  = gb+Nodes
  enddo
  if (lnb>0) then
    !p = this%p_first
    gio = this%p_xrows(Node+1)+1
    p = node_handling_element(dit,gio)+1
    do b=1, lnb
      call MPI_Irecv( this%p_blonz(b), 1, MPI_Integer, &
          p-1, 0, MPI_COMM_WORLD, this%rreq(b), ierr )
      p = merge(1,p+1,p==Nodes)
    enddo
    call MPI_WAitAll( lnb, this%rreq, MPI_STATUSES_IGNORE, ierr )
  endif
  call MPI_WAitAll( this%nb, this%sreq, MPI_STATUSES_IGNORE, ierr )
  ! Number of Nzs of current distribution
  this%nnz = sum(this%p_blonz)
# else /* NEW_NCDF==5 */
  lnb = this%p_nbloc(node+1)
  lno = this%p_nrows(node+1)
  nullify(this%g_ncol,this%blonz,this%p_blonz)
  call re_alloc( this%g_ncol, 1, lno, "g_ncol", "t_ncdfSparse" )
  call re_alloc( this%blonz, 1, this%nb, "blonz", "t_ncdfSparse" )
  call re_alloc( this%p_blonz, 1, lnb, "p_blonz", "t_ncdfSparse" )
  ! Compute NZ of every block of local distribution
  io = 1
  do b=1, this%nb
    gio = index_local_to_global( dit, io )
    n = count_consecutive( dit, g_no, gio )
    this%blonz(b) = sum(ncol(io:io+n-1))
    io = io+n
  enddo
  nullify(nsend,nrecv,xsend,xrecv,ibuf)
  call re_alloc( nsend, 1, Nodes, "nsend", "t_ncdfSparse" )
  call re_alloc( nrecv, 1, Nodes, "nrecv", "t_ncdfSparse" )
  call re_alloc( xsend, 1, Nodes, "xsend", "t_ncdfSparse" )
  call re_alloc( xrecv, 1, Nodes, "xrecv", "t_ncdfSparse" )
  call re_alloc( ibuf, 1, lno, "ibuf", "t_ncdfSparse" )
  nsend(1:Nodes) = 0
  io = 1
  p = 1
  do b=1, this%nb
    gio = index_local_to_global( dit, io )
    do while( gio>this%p_xrows(p)+this%p_nrows(p))
      p = p+1
    enddo
    n = count_consecutive( dit, g_no, gio )
    nsend(p) = nsend(p) + n
    io = io+n
  enddo
  n = 0
  do p=1, Nodes
    xsend(p) = n
    n = n+nsend(p)
  enddo
  nrecv(1:Nodes) = 0
  gio = this%p_xrows(Node+1)+1
  do b=1, lnb
    p = node_handling_element( dit, gio )+1
    n = count_consecutive( dit, g_no, gio )
    nrecv(p) = nrecv(p)+n
    gio = gio+n
  enddo
  n = 0
  do p=1, Nodes
    xrecv(p) = n
    n = n+nrecv(p)
  enddo
  call MPI_Alltoallv( ncol, nsend, xsend, MPI_INTEGER, &
                      ibuf, nrecv, xrecv, MPI_INTEGER, &
                      MPI_COMM_WORLD, ierr )
  io = 1
  gio = this%p_xrows(Node+1)+1
  do b=1, lnb
    p = node_handling_element( dit, gio )+1
    n = count_consecutive( dit, g_no, gio )
    i = xrecv(p)+1
    this%g_ncol(io:io+n-1) = ibuf(i:i+n-1)
    io = io+n
    gio = gio+n
    xrecv(p) = xrecv(p)+n
  enddo
  call de_alloc( ibuf, "ibuf", "t_ncdfSparse" )
  ! Compute the number of non zeros to send/receive
  this%nsend => nsend
  this%xsend => xsend
  this%nrecv => nrecv
  this%xrecv => xrecv
  this%nsend(1:Nodes) = 0
  io = 1
  p = 1
  do b=1, this%nb
    gio = index_local_to_global( dit, io )
    do while( gio>this%p_xrows(p)+this%p_nrows(p))
      p = p+1
    enddo
    n = count_consecutive( dit, g_no, gio )
    this%nsend(p) = this%nsend(p) + sum(ncol(io:io+n-1))
    io = io+n
  enddo
  n = 0
  do p=1, Nodes
    this%xsend(p) = n
    n = n+this%nsend(p)
  enddo
  this%nrecv(1:Nodes) = 0
  gio = this%p_xrows(Node+1)+1
  io = 1
  do b=1, lnb
    p = node_handling_element( dit, gio )+1
    n = count_consecutive( dit, g_no, gio )
    this%nrecv(p) = this%nrecv(p)+sum(this%g_ncol(io:io+n-1))
    gio = gio+n
    io = io+n
  enddo
  n = 0
  do p=1, Nodes
    this%xrecv(p) = n
    n = n+this%nrecv(p)
  enddo
  ! Number of Nzs of current distribution
  this%nnz = this%nrecv(Nodes)+this%xrecv(Nodes)
# endif /* NEW_NCDF==5 */
  this%bsize = this%nnz
#	if NEW_NCDF==4 || NEW_NCDF==5
  nullify(this%p_nnz)
  call re_alloc( this%p_nnz, 1, merge(Nodes,0,this%IOnode), &
                 "p_nnz", "t_ncdfSparse" )
  call MPI_Gather( this%nnz, 1, MPI_INTEGER, this%p_nnz, &
                   1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )
  this%bsize = max(maxval(this%p_nnz),this%bsize)
# else
  call MPI_Scan( this%nnz, n, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr )
  this%wroff = n-this%nnz+1
# endif
  this%g_no = g_no
  this%memset =.true.
# endif   /* MPI */
  end subroutine init_Sp

  subroutine write_Sp( this, grp, sp )
  implicit none
  class(t_ncdfSparse) :: this
  integer :: grp
  type(Sparsity), intent(inout) :: sp
  integer :: var, i, j, o, go, p, b, gb, n, nz, nb, p_nb
  integer :: no, g_no, p_no, ierr
  integer, pointer :: ncol(:), l_col(:), rbuf(:), wbuf(:), ioff(:)
# ifdef MPI
  integer :: MPIstatus(MPI_STATUS_SIZE)
# endif
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
    ! allocate buffers
#   if NEW_NCDF<5
    nullify(rbuf,ioff)
    call re_alloc( rbuf, 1, this%bsize, "rbuf", "t_ncdfSparse" )
    call re_alloc( ioff, 1, Nodes, "ioff", "t_ncdfSparse" )
#   endif
    nullify(wbuf)
    call re_alloc( wbuf, 1, this%bsize, "wbuf", "t_ncdfSparse" )
    ! Reduce and reorder ncol in every reductor process (if necessary)
#	  if NEW_NCDF==5
    call attach( sp, nrows=no, nrows_g=g_no, n_col=ncol, list_col=l_col )
    o = 1                   ! Local  Offset
    p = 1
    do b=1, this%nb
      go = index_local_to_global( this%dit, o )
      n = count_consecutive( this%dit, g_no, go )
      do while(go>this%p_nrows(p)+this%p_xrows(p))
        p=p+1
      enddo
      call MPI_Isend( ncol(o), n, MPI_Integer, p-1, 0, &
          MPI_COMM_WORLD, this%sreq(b), ierr )
      o = o+n
    enddo
    p_no = this%p_nrows(Node+1)
    p_nb = this%p_nbloc(Node+1)
    this%g_ncol => wbuf
    if (p_no>0) then
      o = 1
      go = this%p_xrows(Node+1)+1
      p = node_handling_element(this%dit,go)+1
      do b=1, p_nb
        n = count_consecutive( this%dit, g_no, go )
        call MPI_Irecv( this%g_ncol(o), n, MPI_Integer, &
            p-1, 0, MPI_COMM_WORLD, this%rreq(b), ierr )
        o = o+n
        go = go+n
        p = merge(1,p+1,p==Nodes)
      enddo
      call MPI_WAitAll( p_nb, this%rreq, MPI_STATUSES_IGNORE, ierr )
    endif
    call MPI_WAitAll( this%nb, this%sreq, MPI_STATUSES_IGNORE, ierr )
#	  endif

    ! Reduce data to IOnodes and write to disk
#	  if NEW_NCDF<4
    if (this%IOnode) then
      o = this%p_xrows(node+1)+1
      n = this%p_nrows(node+1)
      call cdf_err( nf90_inq_varid( grp, 'n_col', var ) )
#     if NEW_NCDF==3
  	  call cdf_err( nf90_var_par_access( grp, var, nf90_collective) )
#	    endif
      call cdf_err( nf90_put_var( grp, var, this%g_ncol, &
          start=(/ o /), count=(/ n /) ) )
    endif
#   else
    if (this%IOnode) then
      call cdf_err( nf90_inq_varid( grp, 'n_col', var ) )
      o = 1
      n = this%p_nrows(1)
      call cdf_err( nf90_put_var( grp, var, this%g_ncol, &
            start=(/ o /), count=(/ n /) ) )
      do p=2, Nodes
        o = o+n
        n = this%p_nrows(p)
        if (n/=0) then
          call MPI_Recv( wbuf, n, MPI_Integer, &
                   p-1, 0, MPI_Comm_World, MPIstatus, ierr )
          call cdf_err( nf90_put_var( grp, var, wbuf, &
              start=(/ o /), count=(/ n /) ) )
        endif
      enddo
    else
      n = this%p_nrows(Node+1)
      if (n/=0) call MPI_Send( this%g_ncol, n, MPI_Integer, 0, 0, &
                                MPI_Comm_World, ierr )
    endif
#   endif   /* if NEW_NCDF<4 */

    ! Reduce and reorder l_col in every reductor process (if necessary)
#	  if NEW_NCDF<5
    call attach( sp, nrows_g=g_no, list_col=l_col )
    call MPI_Alltoallv( l_col, this%nsend, this%xsend, MPI_INTEGER, &
                        rbuf, this%nrecv, this%xrecv, MPI_INTEGER, &
                        MPI_COMM_WORLD, ierr )

    ioff(1:Nodes) = this%xrecv(1:Nodes)
    nb = this%p_nbloc(node+1)
    o = 1
    i = 1
    go = this%p_xrows(Node+1)+1
    do b=1, nb
      p = node_handling_element( this%dit, go )+1
      n = count_consecutive( this%dit, g_no, go )
      nz = sum(this%g_ncol(o:o+n-1))
      j = ioff(p)+1
      wbuf(i:i+nz-1) = rbuf(j:j+nz-1)
      o = o+n
      go = go+n
      i = i+nz
      ioff(p) = ioff(p)+nz
    enddo
#   else /* NEW_NCDF==5 */
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
      go = this%p_xrows(Node+1)+1
      p = node_handling_element(this%dit,go)+1
      o = 1
      do b=1, p_nb
        nz = this%p_blonz(b)
        call MPI_Irecv( wbuf(o), nz, MPI_Integer, &
            p-1, 0, MPI_COMM_WORLD, this%rreq(b), ierr )
        o = o+nz
        p = merge(1,p+1,p==Nodes)
      enddo
      call MPI_WAitAll( p_nb, this%rreq, MPI_STATUSES_IGNORE, ierr )
    endif
    call MPI_WAitAll( this%nb, this%sreq, MPI_STATUSES_IGNORE, ierr )
#   endif /* NEW_NCDF==5 */

    ! Reduce data to IOnodes and write to disk
#	  if NEW_NCDF<4
    if (this%nnz>0) then
      call cdf_err( nf90_inq_varid( grp, 'list_col', var ) )
#     if NEW_NCDF==3
  	  call cdf_err( nf90_var_par_access( grp, var, nf90_collective) )
#	    endif
      call cdf_err( nf90_put_var( grp, var, wbuf, &
          start=(/ this%wroff /), count=(/ this%nnz /) ) )
    endif
#   else  /* NEW_NCDF<4 */
    if (this%IOnode) then
      call cdf_err( nf90_inq_varid( grp, 'list_col', var ) )
      o = 1 
      do p=1, Nodes
        n = this%p_nnz(p)
        if (n>0) then
          if (p/=1) call MPI_Recv( wbuf, n, MPI_Integer, &
                        p-1, 0, MPI_Comm_World, MPIstatus, ierr )
          call cdf_err( nf90_put_var( grp, var, wbuf, &
              start=(/ o /), count=(/ n /) ) )
        endif
        o = o+n
      enddo
    else
      n = this%nnz
      if (n>0) call MPI_Send( wbuf, n, MPI_Integer, 0, 0, &
                              MPI_Comm_World, ierr )
    endif
#   endif /* NEW_NCDF<4 */
    call de_alloc( wbuf, "wbuf", "t_ncdfSparse" )
#	  if NEW_NCDF<5
    call de_alloc( ioff, "ioff", "t_ncdfSparse" )
    call de_alloc( rbuf, "rbuf", "t_ncdfSparse" )
#   endif /* NEW_NCDF<5 */
# endif /* MPI */
  endif
  end subroutine write_Sp

  subroutine write_d1D( this, grp, vname, dSp1D )
  use class_dSpData1D
  implicit none
  class(t_ncdfSparse) :: this
  integer :: grp
  character(len=*), intent(in) :: vname
  type(dSpData1D), intent(inout) :: dSp1D
  integer :: var, ierr, nb, o, i, go, b, gb, p_nb, p, n, nz, j
  integer,  pointer :: ioff(:)
  real(dp), pointer :: a(:), rbuf(:), wbuf(:)
# ifdef MPI
  integer :: MPIstatus(MPI_STATUS_SIZE)
# endif
  a => val(dSp1D)
  ! Write the matrix to the file...
  if (nodes==1) then
    call cdf_err( nf90_inq_varid( grp, vname, var ) )
    call cdf_err( nf90_put_var( grp, var, a ) )
#	ifdef MPI
  else
    ! allocate buffers
#   if NEW_NCDF<5
    nullify(rbuf,ioff)
    call re_alloc( rbuf, 1, this%bsize, "rbuf", "t_ncdfSparse" )
    call re_alloc( ioff, 1, Nodes, "ioff", "t_ncdfSparse" )
#   endif
    nullify(wbuf)
    call re_alloc( wbuf, 1, this%bsize, "wbuf", "t_ncdfSparse" )
#	  if NEW_NCDF<5
    call MPI_Alltoallv( a, this%nsend, this%xsend, MPI_Double_Precision, &
                        rbuf, this%nrecv, this%xrecv, MPI_Double_Precision, &
                        MPI_COMM_WORLD, ierr )

    ioff(1:Nodes) = this%xrecv(1:Nodes)
    nb = this%p_nbloc(node+1)
    o = 1
    i = 1
    go = this%p_xrows(Node+1)+1
    do b=1, nb
      p = node_handling_element( this%dit, go )+1
      n = count_consecutive( this%dit, this%g_no, go )
      nz = sum(this%g_ncol(o:o+n-1))
      j = ioff(p)+1
      wbuf(i:i+nz-1) = rbuf(j:j+nz-1)
      o = o+n
      go = go+n
      i = i+nz
      ioff(p) = ioff(p)+nz
    enddo
#   else /* NEW_NCDF==5 */
    o = 1
    gb = Node+1 ! Global block
    p = 1
    do b=1, this%nb
      do while(gb>this%p_nbloc(p)+this%p_xbloc(p))
        p=p+1
      enddo
      nz = this%blonz(b)
      call MPI_Isend( a(o), nz, MPI_Double_Precision, p-1, 0, &
          MPI_COMM_WORLD, this%sreq(b), ierr )
      o = o+nz
      gb = gb+Nodes
    enddo
    p_nb = this%p_nbloc(Node+1)
    if (p_nb>0) then
      go = this%p_xrows(Node+1)+1
      p = node_handling_element(this%dit,go)+1
      o = 1
      do b=1, p_nb
        nz = this%p_blonz(b)
        call MPI_Irecv( wbuf(o), nz, MPI_Double_Precision, &
            p-1, 0, MPI_COMM_WORLD, this%rreq(b), ierr )
        o = o+nz
        p = merge(1,p+1,p==Nodes)
      enddo
      call MPI_WAitAll( p_nb, this%rreq, MPI_STATUSES_IGNORE, ierr )
    endif
    call MPI_WAitAll( this%nb, this%sreq, MPI_STATUSES_IGNORE, ierr )
#   endif /* NEW_NCDF==5 */

#   if NEW_NCDF<4
    if (this%nnz>0) then
      call cdf_err( nf90_inq_varid( grp, vname, var ) )
#     if NEW_NCDF==3
  	  call cdf_err( nf90_var_par_access( grp, var, nf90_collective) )
#	    endif
      call cdf_err( nf90_put_var( grp, var, wbuf, &
          start=(/ this%wroff /), count=(/ this%nnz /) ) )
    endif
#	  else
    if (this%IOnode) then
      call cdf_err( nf90_inq_varid( grp, vname, var ) )
      o = 1
      do p=1, Nodes
        n = this%p_nnz(p)
        if (n>0) then
          if (p/=1) call MPI_Recv( wbuf, n, MPI_Double_Precision, &
                        p-1, 0, MPI_Comm_World, MPIstatus, ierr )
          call cdf_err( nf90_put_var( grp, var, wbuf, &
              start=(/ o /), count=(/ n /) ) )
        endif
        o = o+n
      enddo
    else
      n = this%nnz
      if (n>0) call MPI_Send( wbuf, n, MPI_Double_Precision, &
                    0, 0, MPI_Comm_World, ierr )
    endif
#   endif /* NEW_NCDF>=4 */
#	  if NEW_NCDF<5
    call de_alloc( ioff, "ioff", "t_ncdfSparse" )
    call de_alloc( rbuf, "rbuf", "t_ncdfSparse" )
#	  endif /* NEW_NCDF==5 */
    call de_alloc( wbuf, "wbuf", "t_ncdfSparse" )
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
  integer :: var, ierr, spdim, niters, width, it, nb, o, i
  integer :: go, b, gb, p_nb, p, n, nz, j
  real(dp), pointer :: a(:,:), rbuf(:,:), wbuf(:,:)
  integer,  pointer :: ioff(:)
# ifdef MPI
  integer :: MPIstatus(MPI_STATUS_SIZE)
# endif
  a => val(dSp2D)
  ! Write the sparsity to the file...
  if (nodes==1) then
    call cdf_err( nf90_inq_varid( grp, vname, var ) )
    call cdf_err( nf90_put_var( grp, var, a ) )
#	ifdef MPI
  else
    spdim=spar_dim(dSp2D)
    nullify(rbuf,wbuf)
    if (spdim==1) then
      niters = size(a,2)
      width  = 1
    else
      niters = 1
      width  = size(a,1)
    endif
#   if NEW_NCDF<5
    nullify(rbuf,ioff)
    call re_alloc( rbuf, 1, width, 1, this%bsize, "rbuf", "t_ncdfSparse" )
    call re_alloc( ioff, 1, Nodes, "ioff", "t_ncdfSparse" )
#   endif
    nullify(wbuf)
    call re_alloc( wbuf, 1, width, 1, this%bsize, "wbuf", "t_ncdfSparse" )
    if (this%IOnode) call cdf_err( nf90_inq_varid( grp, vname, var ) )
    do it=1, niters
#	    if NEW_NCDF<5
      if (spdim==1) then
        call MPI_Alltoallv( a(1,it), this%nsend, this%xsend, &
          MPI_Double_Precision, rbuf(1,1), this%nrecv, this%xrecv, &
          MPI_Double_Precision, MPI_COMM_WORLD, ierr )
      else
        call MPI_Alltoallv( a(1,1), this%nsend*width, this%xsend*width, &
          MPI_Double_Precision, rbuf(1,1), this%nrecv*width, this%xrecv*width, &
          MPI_Double_Precision, MPI_COMM_WORLD, ierr )
      endif
      ioff(1:Nodes) = this%xrecv(1:Nodes)
      nb = this%p_nbloc(node+1)
      o = 1
      i = 1
      go = this%p_xrows(Node+1)+1
      do b=1, nb
        p = node_handling_element( this%dit, go )+1
        n = count_consecutive( this%dit, this%g_no, go )
        nz = sum(this%g_ncol(o:o+n-1))
        j = ioff(p)+1
        wbuf(:,i:i+nz-1) = rbuf(:,j:j+nz-1)
        o = o+n
        go = go+n
        i = i+nz
        ioff(p) = ioff(p)+nz
      enddo
#     else
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
        go = this%p_xrows(Node+1)+1
        p = node_handling_element(this%dit,go)+1
        o = 1
        do b=1, p_nb
          nz = this%p_blonz(b)
          call MPI_Irecv( wbuf(1,o), nz*width, MPI_Double_precision, &
              p-1, 0, MPI_COMM_WORLD, this%rreq(b), ierr )
          o = o+nz
          p = merge(1,p+1,p==Nodes)
        enddo
        call MPI_WAitAll( p_nb, this%rreq, MPI_STATUSES_IGNORE, ierr )
      endif
      call MPI_WAitAll( this%nb, this%sreq, MPI_STATUSES_IGNORE, ierr )
#     endif
#	    if NEW_NCDF<4
      if (this%nnz>0) then
        if (spdim==1) then
          call cdf_err( nf90_put_var( grp, var, wbuf, &
              start=(/ this%wroff, it /), count=(/ this%nnz,1 /) ) )
        else
          call cdf_err( nf90_put_var( grp, var, wbuf, &
              start=(/1,this%wroff/), count=(/width,this%nnz/) ) )
        endif
      endif
#     else /* NEW_NCDF>=4 */
      if (this%IOnode) then
        o = 1
        do p=1, Nodes
          n = this%p_nnz(p)
          if (p/=1) then
            call MPI_Recv( wbuf(1,1), n*width, MPI_Double_Precision, &
                             p-1, 0, MPI_Comm_World, MPIstatus, ierr )
          endif
          if (spdim==1) then
            call cdf_err( nf90_put_var( grp, var, wbuf, &
                start=(/ o, it /), count=(/ n,1 /) ) )
          else
            call cdf_err( nf90_put_var( grp, var, wbuf, &
                start=(/1,o/), count=(/width,n/) ) )
          endif
          o = o+n
        enddo
      else
        n = this%nnz
        call MPI_Send( wbuf(1,1), n*width, MPI_Double_Precision, 0, 0, &
                       MPI_Comm_World, ierr )
      endif
#     endif /* NEW_NCDF>=4 */
    enddo
#	  if NEW_NCDF<5
    call de_alloc( rbuf, "rbuf", "t_ncdfSparse" )
    call de_alloc( wbuf, "wbuf", "t_ncdfSparse" )
#   endif /* NEW_NCDF<5 */
    call de_alloc( ioff, "ioff", "t_ncdfSparse" )
# endif /* MPI */
  endif
  end subroutine write_d2D

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
