! Module for easy writing of all information
! in a single netCDF file.

! We rely on the NetCDF4 library and can
! utilise parallel IO or standard IO

module m_ncdf_siesta

  use precision, only : sp, dp, grid_p
  use parallel, only : Node, Nodes
#ifdef NCDF_4
  use siesta_options, only: cdf_comp_lvl, cdf_w_parallel
  use variable
  use dictionary
  use netcdf_ncdf, ncdf_parallel => parallel
  use ncdf_io_m
#endif
#ifdef MPI
  use mpi_siesta, only : MPI_INTEGER, MPI_INTEGER8, MPI_SUM, &
                         MPI_GRID_REAL, MPI_STATUS_SIZE
  use mpi, only : MPI_COMM_WORLD, MPI_INFO_NULL, &
                  MPI_STATUSES_IGNORE
#endif
  use alloc
  use m_host, only : getCHost, getAvailMem
  implicit none

  private

#ifdef NCDF_4

  public :: cdf_init_file
  public :: cdf_init_grid

  public :: cdf_save_settings
  public :: cdf_save_basis
  public :: cdf_save_fc
  public :: cdf_save_state
  public :: cdf_save_grid

  public :: cdf_init_file_new
  public :: cdf_init_grid_new
  public :: cdf_save_settings_new
  public :: cdf_save_basis_new
  public :: cdf_save_state_new
  public :: cdf_save_grid_new

#endif

# if NEW_NCDF==5
  type t_ncdfGrid
    integer :: ntm(3)
    integer :: grsize
    integer :: nred
    logical :: Iamred
    integer :: nsend
    integer :: nrecv
    integer :: bsize
    integer :: redid
    integer, pointer :: ndist(:)
    integer, pointer :: scount(:)
    integer, pointer :: rcount(:)
    contains
    procedure :: init_ng, getDims, getArrays
  end type t_ncdfGrid
# else
  type t_ncdfGrid
    logical :: arr_init
    integer :: ntm(3), nwriters, IamW, nzoff
    integer, pointer :: wrtrs(:), nzs(:), sndnz(:), count(:), displ(:)
    contains
    procedure :: init_ng, getDims, getArrays
  end type t_ncdfGrid
# endif

  type(t_ncdfGrid) :: cdfGrid

contains

# if NEW_NCDF==5
  subroutine init_ng( this, ntm )
  use moreMeshSubs, only : getLocalBox
  use fdf, only : fdf_get
  implicit none
  class(t_ncdfGrid) :: this
  integer :: ntm(3)
  integer :: nproc_width
  integer :: grsize, nred, nz, di, mo, i, j, o, l, r, bsize
  integer :: lbox(2,3)
  !integer :: myred
  logical :: Iamred
  integer, pointer :: ndist(:), scount(:), rcount(:)
  this%ntm = ntm

# ifdef MPI
  if (nodes>1) then
     nproc_width = fdf_get("CDF.NPROC.WIDTH",12)
     print *, "nproc_width: ", nproc_width
    grsize = nproc_width
    Iamred = mod(node,grsize)==0
    ! Split Z dimension among the reductor processors
    nz = ntm(3)
    nred = (Nodes+grsize-1)/grsize
    ! Save the offset of every cut in Z simension
    nullify(ndist)
    allocate(ndist(nred+1))
    di = nz/nred
    mo = mod(nz,nred)
    o = 1
    do i=1, nred
      ndist(i) = o
      o = o+di+merge(1,0,i<=mo)
    enddo
    ndist(nred+1) = nz+1
    ! Number of XY planes to send to every reductor
    nullify(scount)
    allocate(scount(nred))
    call getLocalBox( 1, Node+1, lbox ) ! Distribution == 1
    do i=1, nred
      l = max(ndist(i),lbox(1,3))
      r = min(ndist(i+1)-1,lbox(2,3))
      scount(i) = max(r-l+1,0)
    enddo
    ! Number of XY planes to receive from every process
    nullify(rcount)
    if (Iamred) then
      allocate(rcount(Nodes))
      j = Node/grsize+1
      do i=1, Nodes
        call getLocalBox( 1, i, lbox )
        l = max(ndist(j),lbox(1,3))
        r = min(ndist(j+1)-1,lbox(2,3))
        rcount(i) = max(r-l+1,0)
      enddo
      bsize = ntm(1)*ntm(2)*(ndist(j+1)-ndist(j))
    else
      allocate(rcount(0))
      bsize = 0
    endif
    !myred = (Node/grsize)*grsize
    ! Fill data type
    this%grsize = grsize
    this%nred   = nred
    this%Iamred = Iamred
    this%redid  = j
    this%nsend  = sum(scount)
    this%nrecv  = sum(rcount)
    this%bsize  = bsize
    this%ndist  => ndist
    this%scount => scount
    this%rcount => rcount
  endif
# endif /* MPI */
  end subroutine init_ng

# else /* NEW_NCDF/=5 */
  subroutine init_ng( this, ntm )
  implicit none
  class(t_ncdfGrid) :: this
  integer :: ntm(3)

  integer :: nhost, n, i
  integer*8 :: req_mem, ava_mem
  character*128 :: nw
  this%ntm = ntm
  ! Number of computing nodes => All processes in the same computing node
  ! share the available memory
  nullify(this%wrtrs,this%nzs,this%sndnz,this%count,this%displ)
  call getCHost( nhost, this%wrtrs )
  ! Compute required memory to gather and write the data
  n = ntm(1)*ntm(2)*ntm(3)
  req_mem = n*8*2  ! rbuf(1:n) wbuf(1:n)
  req_mem = req_mem + nodes*4*2 ! count(nodes) displ(nodes)
  req_mem = req_mem + nodes*4*2 ! (sndnz(nwriter)+nzs(nwriter))*cpusXhost
  req_mem = int(1.05*req_mem)   ! Let's use a 5% of safety margin
  ava_mem = getAvailMem( )
  ! Set the number of writers
  this%nwriters = req_mem/ava_mem + 1
  call getenv( "SIESTA_NWRITERS", nw )
  if (trim(nw)/='') read(nw,*) this%nwriters
  this%nwriters = min(nhost,max(1,this%nwriters))

  this%IamW = 0
  do i=1, this%nwriters
    if (node==this%wrtrs(i)) then
      this%IamW = i
      exit
    endif
  enddo

  n = merge(nodes,1,this%IamW>0)
  call re_alloc( this%nzs, 1, this%nwriters, 'nzs', 't_ncdfGrid' )
  call re_alloc( this%sndnz, 1, this%nwriters, 'sndnz', 't_ncdfGrid' )
  call re_alloc( this%count, 1, n, 'count', 't_ncdfGrid' )
  call re_alloc( this%displ, 1, n, 'displ', 't_ncdfGrid' )
  this%arr_init = .false.
  end subroutine init_ng
# endif

  subroutine getDims( this, nx, ny, nz )
  implicit none
  class(t_ncdfGrid) :: this
  integer :: nx, ny, nz
  nx = this%ntm(1)
  ny = this%ntm(2)
  nz = this%ntm(3)
  end subroutine getDims

  subroutine getArrays( this, nwriters, IamW, nzoff, wrtrs, nzs, &
                        sndnz, count, displ )
  use moreMeshSubs, only : getLocalBox
  implicit none
  class(t_ncdfGrid) :: this
  integer :: nwriters, IamW, nzoff
  integer, pointer :: wrtrs(:), nzs(:), sndnz(:), count(:), displ(:)
# if NEW_NCDF==5
# else

  integer :: lbox(2,3), di, mo, j, nx, ny, nz, lx, ly, lz, i, p
  integer :: ierr, n
  nwriters = this%nwriters
  IamW     = this%IamW
  wrtrs   => this%wrtrs
  nzs     => this%nzs
  sndnz   => this%sndnz
  count   => this%count
  displ   => this%displ
  if (.not. this%arr_init) then
    this%arr_init = .true.
    call getLocalBox( 1, Node+1, lbox ) ! Distribution == 1
    lx = lbox(2,1)-lbox(1,1)+1
    ly = lbox(2,2)-lbox(1,2)+1

    nx = this%ntm(1)
    ny = this%ntm(2)
    nz = this%ntm(3)
    di = nz/nwriters
    mo = mod(nz,nwriters)
    j = 0
    do i=1, nwriters
      p = wrtrs(i)
      if (node==p) this%nzoff = j
      nzs(i) = di + merge(1,0,i<=mo)
      sndnz(i) = min(j+nzs(i),lbox(2,3))-max(j,lbox(1,3)-1)
      sndnz(i) = max(sndnz(i),0)*lx*ly
      call MPI_Gather( sndnz(i), 1, MPI_Integer, count, 1, &
          MPI_integer, p, MPI_COMM_WORLD, ierr )
      j = j + nzs(i)
    enddo
    if (IamW>0) then
      j=0
      do i=1, nodes
        displ(i) = j
        j = j + count(i)
      enddo
    endif
  endif
  nzoff = this%nzoff
# endif
  end subroutine getArrays

#ifdef NCDF_4

  subroutine cdf_err( ierr )
    implicit none
    integer :: ierr
    if (ierr /= NF90_NOERR) then
      print *, node, trim(nf90_strerror(ierr))
      call MPI_Abort( MPI_COMM_WORLD, 0, ierr )
    end if
  end subroutine cdf_err

  subroutine cdf_init_file(fname, is_fc)
    use class_Sparsity
    use files, only : slabel
    use atomlist, only: no_u, no_s, lasto, Qtot
    use siesta_geom, only: na_u, nsc
    use sparse_matrices, only: sparse_pattern
    use m_spin, only : spin
    use siesta_options, only: sname, isolve
    use siesta_options, only: SOLVE_DIAGON, SOLVE_ORDERN
    use siesta_options, only: SOLVE_MINIM, SOLVE_TRANSI
    use siesta_options, only: SOLVE_PEXSI
    use siesta_options, only: savehs
    use siesta_options, only: fixspin, total_spin
    use siesta_options, only: ia1, ia2, dx ! FC information
    use m_timestamp, only: datestring
    use ts_electrode_m, only: electrode_t
    use m_ts_options, only : Volt, N_Elec, Elecs
    use m_ts_options, only : TS_HS_Save

    character(len=*), intent(in) :: fname
    logical, intent(in), optional :: is_fc

    ! Local variables
    type(hNCDF) :: ncdf, grp, grp2
    type(dictionary_t) :: dic, d
    integer :: n_nzs, tmp, i, chks(3), iEl
    integer, allocatable :: ibuf(:)
    logical :: lis_fc
#ifdef MPI
    integer :: MPIerror
#endif

    lis_fc = .false.
    if ( present(is_fc) ) lis_fc = is_fc

    ! We always re-write the file...
    call ncdf_create(ncdf,fname, &
        mode=ior(NF90_WRITE,NF90_NETCDF4), overwrite=.true.)

#ifdef MPI
    tmp = nnzs(sparse_pattern)
    call MPI_Reduce(tmp,n_nzs,1,MPI_Integer, &
        MPI_Sum, 0, MPI_Comm_World, MPIerror)
#else
    n_nzs = nnzs(sparse_pattern)
#endif

    ! First we create all the dimensions
    ! necessary
    d = ('DIMna_u'.kv.na_u)//('DIMno_u'.kv.no_u)
    d = d//('DIMno_s'.kv.no_s)//('DIMspin'.kv.spin%H)
    d = d//('DIMxyz'.kv. 3) // ('DIMn_s'.kv. product(nsc))
    d = d// ('DIMone'.kv. 1)
    call ncdf_crt(ncdf,d)
    call delete(d)

    dic = dic//('info'.kv.'Last orbital of equivalent atom')
    call ncdf_def_var(ncdf,'lasto',NF90_INT,(/'na_u'/), &
        compress_lvl=0, atts=dic)
    call ncdf_put_var(ncdf,'lasto',lasto(1:na_u))

    dic = dic//('info'.kv.'Total charge')
    call ncdf_def_var(ncdf,'Qtot',NF90_DOUBLE,(/'one'/), &
        compress_lvl=0, atts=dic)
    call ncdf_put_var(ncdf,'Qtot',Qtot)

    if ( fixspin ) then
      dic = dic//('info'.kv.'Total spin')
      call ncdf_def_var(ncdf,'Qspin',NF90_DOUBLE,(/'one'/), &
          compress_lvl=0, atts=dic)
      call ncdf_put_var(ncdf,'Qspin',total_spin)
    end if

    ! Create all necessary containers...
    dic = dic//('info'.kv.'Number of supercells in each unit-cell direction') 
    call ncdf_def_var(ncdf,'nsc',NF90_INT,(/'xyz'/), &
        compress_lvl=0, atts=dic)

    dic = dic//('info'.kv.'Fermi level')//('unit'.kv.'Ry')
    if ( fixspin ) then
      call ncdf_def_var(ncdf,'Ef',NF90_DOUBLE,(/'spin'/), &
          compress_lvl=0, atts=dic)
    else
      call ncdf_def_var(ncdf,'Ef',NF90_DOUBLE,(/'one'/), &
          compress_lvl=0, atts=dic)
    end if

    dic = dic//('info'.kv.'Atomic coordinates')//('unit'.kv.'Bohr')
    call ncdf_def_var(ncdf,'xa',NF90_DOUBLE,(/'xyz ','na_u'/), &
        compress_lvl=0, atts=dic)

    dic = dic//('info'.kv.'Unit cell')
    call ncdf_def_var(ncdf,'cell',NF90_DOUBLE,(/'xyz','xyz'/), &
        compress_lvl=0, atts=dic)

    dic = dic//('info'.kv.'Atomic forces')//('unit'.kv.'Ry/Bohr')
    call ncdf_def_var(ncdf,'fa',NF90_DOUBLE,(/'xyz ','na_u'/), &
        compress_lvl=0, atts=dic)

    dic = dic//('info'.kv.'Cell stress')//('unit'.kv.'Ry/Bohr**3')
    call ncdf_def_var(ncdf,'stress',NF90_DOUBLE,(/'xyz','xyz'/), &
        compress_lvl=0, atts=dic)

    call delete(dic)

    ! Create matrix group
    call ncdf_def_grp(ncdf,'SPARSE',grp)

    call ncdf_def_dim(grp,'nnzs',n_nzs)

    ! EDM is required to have its own spin (because of spin-orbit coupling
    ! where the spin == 4, and not 8)
    call ncdf_def_dim(grp,'spin_EDM',spin%EDM)


    dic = dic//('info'.kv.'Index of supercell coordinates')
    call ncdf_def_var(grp,'isc_off',NF90_INT,(/'xyz','n_s'/), &
        compress_lvl=0, atts=dic)

    dic = dic//('info'.kv.'Number of non-zero elements per row')
    call ncdf_def_var(grp,'n_col',NF90_INT,(/'no_u'/), &
        compress_lvl=0, atts=dic)

    chks = (/n_nzs,1,1/)

    dic = dic//('info'.kv.'Supercell column indices in the sparse format')
    call ncdf_def_var(grp,'list_col',NF90_INT,(/'nnzs'/), &
        compress_lvl=cdf_comp_lvl,atts=dic,chunks=chks)

    ! Create the overlap matrix (we know it will not change)
    dic = dic//('info'.kv.'Overlap matrix')
    call ncdf_def_var(grp,'S',NF90_DOUBLE,(/'nnzs'/), &
        compress_lvl=cdf_comp_lvl,atts=dic,chunks=chks)

    dic = dic//('info'.kv.'Overlap matrix gradient')//('unit'.kv.'1/Bohr')
    call ncdf_def_var(grp,'S_gradient',NF90_DOUBLE,(/'xyz ', 'nnzs'/), &
        compress_lvl=cdf_comp_lvl,atts=dic,chunks=(/1,n_nzs/))
    call delete(dic)

    dic = dic//('info'.kv.'Density matrix')
    call ncdf_def_var(grp,'DM',NF90_DOUBLE,(/'nnzs','spin'/), &
        compress_lvl=cdf_comp_lvl,atts=dic,chunks=chks)

    dic = dic//('info'.kv.'Energy density matrix')
    dic = dic//('unit'.kv.'Ry')
    call ncdf_def_var(grp,'EDM',NF90_DOUBLE,(/'nnzs    ','spin_EDM'/), &
        compress_lvl=cdf_comp_lvl,atts=dic,chunks=chks)

    if ( savehs .or. TS_HS_save ) then
      ! Unit is already present in dictionary
      dic = dic//('info'.kv.'Hamiltonian')
      call ncdf_def_var(grp,'H',NF90_DOUBLE,(/'nnzs','spin'/), &
          compress_lvl=cdf_comp_lvl,atts=dic,chunks=chks)
    end if

    ! Even though I think we could do without, I add the
    ! xij array to the file
    ! Note that xij can be re-created using
    !    nsc, isc_off and xa
    !    if ( .not. Gamma ) then
    !       dic = dic//('info'.kv. &
    !            'Distance between orbital i and j')
    !       dic = dic//('unit'.kv.'Bohr')
    !       call ncdf_def_var(grp,'xij',NF90_DOUBLE,(/'xyz ','nnzs'/), &
    !            compress_lvl=cdf_comp_lvl,atts=dic)
    !    end if

    ! Delete the dictionary
    call delete(dic)

    ! Create grid group
    call ncdf_def_grp(ncdf,'GRID',grp)

    ! Note that there are now 2 spin dimensions
    ! 1. in the top level NC file and one in the GRID group
    ! This is because for spin-orbit coupling the grid and
    ! matrix dimensions are not the same (4 vs 8, respectively)
    d = ('DIMspin' .kv. spin%Grid)
    call ncdf_crt(grp, d)

    call delete(d)
    call delete(dic)

    call ncdf_def_grp(ncdf,'SETTINGS',grp)

    dic = ('info'.kv.'Tolerance for converging the density matrix')
    call ncdf_def_var(grp,'DMTolerance',NF90_DOUBLE,(/'one'/), &
        compress_lvl=0, atts=dic)

    dic = dic//('info'.kv.'Tolerance for converging the Hamiltonian')
    call ncdf_def_var(grp,'HTolerance',NF90_DOUBLE,(/'one'/), &
        compress_lvl=0, atts=dic)

    dic = dic//('info'.kv.'Net charge of the system')
    call ncdf_def_var(grp,'NetCharge',NF90_DOUBLE,(/'one'/), &
        compress_lvl=0, atts=dic)

    dic = dic//('info'.kv.'Mixing weight')
    call ncdf_def_var(grp,'MixingWeight',NF90_DOUBLE,(/'one'/), &
        compress_lvl=0, atts=dic)

    dic = dic//('info'.kv.'Grid used for the Brillouin zone integration')
    call ncdf_def_var(grp,'BZ',NF90_INT,(/'xyz','xyz'/), &
        compress_lvl=0, atts=dic)

    dic = dic//('info'.kv.'Grid displacement used in Brillouin zone') &
        //('unit'.kv.'b**-1')
    call ncdf_def_var(grp,'BZ_displ',NF90_DOUBLE,(/'xyz'/), &
        compress_lvl=0, atts=dic)

    dic = dic//('info'.kv.'Temperature for electrons')// &
        ('unit'.kv.'Ry')
    call ncdf_def_var(grp,'ElectronicTemperature',NF90_DOUBLE,(/'one'/), &
        compress_lvl=0, atts=dic)

    dic = dic//('info'.kv.'Mesh cutoff for real space grid')
    call ncdf_def_var(grp,'MeshCutoff',NF90_DOUBLE,(/'one'/), &
        compress_lvl=0, atts=dic)

    call delete(dic)

    ! Create FC group
    if ( lis_FC ) then

      call ncdf_def_grp(ncdf,'FC',grp)

      ! Create dimension arrays
      ! Create arrays for containers
      allocate(ibuf(ia2-ia1+1))
      do i = ia1, ia2
        ibuf(i-ia1+1) = i
      end do
      call ncdf_def_dim(grp,'na_fc',ia2-ia1+1)
      call ncdf_def_dim(grp,'m_p',2)

      dic = dic//('info'.kv.'Displacement length')//('unit'.kv.'Bohr')
      call ncdf_def_var(grp,'disp',NF90_DOUBLE,(/'one'/), &
          compress_lvl=0, atts=dic)
      call ncdf_put_var(grp,'disp',dx)

      call delete(dic)
      dic = dic//('info'.kv.'Displaced atoms')
      call ncdf_def_var(grp,'ia_fc',NF90_INT,(/'na_fc'/), &
          compress_lvl=cdf_comp_lvl, atts=dic)

      call ncdf_put_var(grp,'ia_fc',ibuf)
      deallocate(ibuf)

      dic = dic//('info'.kv.'Undisplaced atomic forces')//('unit'.kv.'Ry/Bohr')
      call ncdf_def_var(grp,'fa0',NF90_DOUBLE,(/'xyz ','na_u'/), &
          compress_lvl=cdf_comp_lvl, atts=dic)

      dic = dic//('info'.kv.'Displaced atomic forces')//('unit'.kv.'Ry/Bohr')
      call ncdf_def_var(grp,'fa',NF90_DOUBLE,(/'xyz  ','na_u ','m_p  ','xyz  ', 'na_fc'/), &
          compress_lvl=cdf_comp_lvl, atts=dic)

    end if

    call delete(dic)

    if ( isolve == SOLVE_TRANSI ) then

      ! Save all information about the transiesta method
      call ncdf_def_grp(ncdf,'TRANSIESTA',grp)

      dic = ('info'.kv.'Grid used for the Brillouin zone integration')
      call ncdf_def_var(grp,'BZ',NF90_INT,(/'xyz','xyz'/), &
          compress_lvl=0, atts=dic)

      dic = dic//('info'.kv.'Grid displacement used in Brillouin zone') &
          //('unit'.kv.'b**-1')
      call ncdf_def_var(grp,'BZ_displ',NF90_DOUBLE,(/'xyz'/), &
          compress_lvl=0, atts=dic)

      dic = dic//('info'.kv.'Applied voltage')//('unit'.kv.'Ry')
      call ncdf_def_var(grp,'Volt',NF90_DOUBLE,(/'one'/), &
          compress_lvl=0,atts=dic)

      dic = dic//('info'.kv.'Hartree correction for FFT Poisson solver')
      call ncdf_def_var(grp,'dHartree',NF90_DOUBLE,(/'one'/), &
          compress_lvl=0,atts=dic)
      call delete(dic)

       ! Add all the electrodes
      do iEl = 1 , N_Elec

        call ncdf_def_grp(grp,trim(Elecs(iEl)%name),grp2)

        tmp = Elecs(iEl)%device_atoms()
        call ncdf_def_dim(grp2,'na',tmp)

        dic = ('info'.kv.'Atoms belonging to electrode')
        call ncdf_def_var(grp2,'a_idx',NF90_INT,(/'na'/), &
            compress_lvl=0,atts=dic)

        dic = dic//('info'.kv.'Chemical potential')//('unit'.kv.'Ry')
        call ncdf_def_var(grp2,'mu',NF90_DOUBLE,(/'one'/), &
            compress_lvl=0,atts=dic)

        dic = dic//('info'.kv.'Electronic temperature')//('unit'.kv.'Ry')
        call ncdf_def_var(grp2,'kT',NF90_DOUBLE,(/'one'/), &
            compress_lvl=0,atts=dic)

        call delete(dic)

      end do

    end if

    ! Save all things necessary here
    dic = ('time'.kv.datestring())
    dic = dic//('name'.kv.trim(sname))
    dic = dic//('label'.kv.trim(slabel))

    if ( isolve == SOLVE_DIAGON ) then
      dic = dic//('method'.kv.'diagon')
    else if ( isolve == SOLVE_ORDERN ) then
      dic = dic//('method'.kv.'order-n')
    else if ( isolve == SOLVE_MINIM ) then
      dic = dic//('method'.kv.'omm')
    else if ( isolve == SOLVE_TRANSI ) then
      dic = dic//('method'.kv.'transiesta')
    else if ( isolve == SOLVE_PEXSI ) then
      dic = dic//('method'.kv.'pexsi')
    end if

    ! Attributes are collective
    call ncdf_put_gatt(ncdf,atts=dic)

    call delete(dic)

    call ncdf_put_var(ncdf,'lasto',lasto(1:na_u))

    if ( isolve == SOLVE_TRANSI ) then

      ! Save all information about the transiesta method
      call ncdf_open_grp(ncdf,'TRANSIESTA',grp)

      call ncdf_put_var(grp,'Volt',Volt)

      ! Add all the electrodes
      do iEl = 1 , N_Elec

        call ncdf_open_grp(grp,trim(Elecs(iEl)%name),grp2)

        tmp = Elecs(iEl)%device_atoms()

        allocate(ibuf(tmp))
        do i = 1 , tmp
          ibuf(i) = Elecs(iEl)%idx_a + i - 1
        end do
        call ncdf_put_var(grp2,'a_idx',ibuf)
        deallocate(ibuf)

        call ncdf_put_var(grp2,'mu',Elecs(iEl)%mu%mu)
        call ncdf_put_var(grp2,'kT',Elecs(iEl)%mu%kT)

      end do

    end if

    ! Close the file
    call ncdf_close(ncdf)

  end subroutine cdf_init_file

  subroutine cdf_init_file_new(fname, is_fc)
    use class_Sparsity
    use files, only : slabel
    use atomlist, only: no_u, no_s, lasto, Qtot
    use siesta_geom, only: na_u, nsc
    use sparse_matrices, only: sparse_pattern
    use m_spin, only : spin
    use siesta_options, only: sname, isolve
    use siesta_options, only: SOLVE_DIAGON, SOLVE_ORDERN
    use siesta_options, only: SOLVE_MINIM, SOLVE_TRANSI
    use siesta_options, only: SOLVE_PEXSI
    use siesta_options, only: savehs
    use siesta_options, only: fixspin, total_spin
    use siesta_options, only: ia1, ia2, dx ! FC information
    use m_timestamp, only: datestring
    use m_ts_options, only : Volt, N_Elec, Elecs
    use m_ts_options, only : TS_HS_Save

    character(len=*), intent(in) :: fname
    logical, intent(in), optional :: is_fc

    ! Local variables
    integer :: ncid, n_nzs, tmp, cmode, i, iEl
    integer :: n_s_, na_u_, spin_, one, no_u_, no_s_, xyz, nnzs_, &
               spin_EDM, spin_Grid, na_fc, m_p, na
    integer :: var, vlasto, vQtot, vQspin, vdisp, vVolt, via_fc
    integer :: grp, grp2, grp_tr
    integer, pointer :: ibuf(:)
    logical :: lis_fc
#ifdef MPI
    integer :: MPIerror
#endif

    lis_fc = .false.
    if ( present(is_fc) ) lis_fc = is_fc

#ifdef MPI
    tmp = nnzs(sparse_pattern)
    call MPI_Reduce(tmp,n_nzs,1,MPI_Integer, &
        MPI_Sum, 0, MPI_Comm_World, MPIerror)
#else
    n_nzs = nnzs(sparse_pattern)
#endif

    if (node/=0) return

    ! We always re-write the file...
    cmode = IOR(NF90_WRITE,NF90_NETCDF4)
    cmode = IOR(cmode,NF90_CLOBBER)
    call cdf_err( nf90_create( fname, cmode, ncid ) )

    ! First we create all the dimensions
    ! necessary
    call cdf_err( nf90_def_dim( ncid, 'n_s', product(nsc), n_s_ ) )
    call cdf_err( nf90_def_dim( ncid, 'na_u', na_u, na_u_ ) )
    call cdf_err( nf90_def_dim( ncid, 'spin', spin%H, spin_ ) )
    call cdf_err( nf90_def_dim( ncid, 'one', 1, one ) )
    call cdf_err( nf90_def_dim( ncid, 'no_u', no_u, no_u_ ) )
    call cdf_err( nf90_def_dim( ncid, 'no_s', no_s, no_s_ ) )
    call cdf_err( nf90_def_dim( ncid, 'xyz', 3, xyz ) )

    call cdf_err( nf90_def_var( ncid, 'lasto', &
          NF90_INT, (/ na_u_ /), vlasto ) )
  	call cdf_err( nf90_put_att( ncid, vlasto, &
          'info', 'Last orbital of equivalent atom' ) )

    call cdf_err( nf90_def_var( ncid, 'Qtot', &
          NF90_DOUBLE, (/ one /), vQtot ) )
  	call cdf_err( nf90_put_att( ncid, vQtot, &
          'info', 'Total charge' ) )

    if ( fixspin ) then
      call cdf_err( nf90_def_var( ncid, 'Qspin', &
            NF90_DOUBLE, (/ one /), vQspin ) )
      call cdf_err( nf90_put_att( ncid, vQspin, &
            'info', 'Total spin' ) )
    end if

    ! Create all necessary containers...
    call cdf_err( nf90_def_var( ncid, 'nsc', &
          NF90_INT, (/ xyz /), var ) )
    call cdf_err( nf90_put_att( ncid, var, 'info', &
          'Number of supercells in each unit-cell direction' ) )

    call cdf_err( nf90_def_var( ncid, 'Ef', &
          NF90_DOUBLE, (/ merge(spin_,one,fixspin) /), var ) )
    call cdf_err( nf90_put_att( ncid, var, 'unit', 'Ry' ) )
    call cdf_err( nf90_put_att( ncid, var, 'info', 'Fermi level' ) )

    call cdf_err( nf90_def_var( ncid, 'xa', &
          NF90_DOUBLE, (/ xyz, na_u_ /), var ) )
    call cdf_err( nf90_put_att( ncid, var, 'unit', 'Bohr' ) )
    call cdf_err( nf90_put_att( ncid, var, &
          'info', 'Atomic coordinates' ) )

    call cdf_err( nf90_def_var( ncid, 'cell', &
          NF90_DOUBLE, (/ xyz, xyz /), var ) )
    call cdf_err( nf90_put_att( ncid, var, 'unit', 'Bohr' ) )
    call cdf_err( nf90_put_att( ncid, var, &
          'info', 'Unit cell' ) )

    call cdf_err( nf90_def_var( ncid, 'fa', &
          NF90_DOUBLE, (/ xyz, na_u_ /), var ) )
    call cdf_err( nf90_put_att( ncid, var, 'unit', 'Ry/Bohr' ) )
    call cdf_err( nf90_put_att( ncid, var, &
          'info', 'Atomic forces' ) )

    call cdf_err( nf90_def_var( ncid, 'stress', &
          NF90_DOUBLE, (/ xyz, xyz /), var ) )
    call cdf_err( nf90_put_att( ncid, var, 'unit', 'Ry/Bohr**3' ) )
    call cdf_err( nf90_put_att( ncid, var, &
          'info', 'Cell stress' ) )

    ! Create matrix group
  	call cdf_err( nf90_def_grp( ncid, 'SPARSE', grp ) )

    call cdf_err( nf90_def_dim( grp, 'nnzs', n_nzs, nnzs_ ) )
    ! EDM is required to have its own spin (because of spin-orbit coupling
    ! where the spin == 4, and not 8)
    call cdf_err( nf90_def_dim( grp, 'spin_EDM', spin%EDM, spin_EDM ) )

    call cdf_err( nf90_def_var( grp, 'isc_off', &
          NF90_INT, (/ xyz, n_s_ /), var ) )
    call cdf_err( nf90_put_att( grp, var, &
          'info', 'Index of supercell coordinates' ) )

    call cdf_err( nf90_def_var( grp, 'n_col', &
          NF90_INT, (/ no_u_ /), var ) )
    call cdf_err( nf90_put_att( grp, var, &
          'info', 'Number of non-zero elements per row' ) )

    call cdf_err( nf90_def_var( grp, 'list_col', &
          NF90_INT, (/ nnzs_ /), var ) )
    call cdf_err( nf90_put_att( grp, var, &
          'info', 'Supercell column indices in the sparse format' ) )

    ! Create the overlap matrix (we know it will not change)
    call cdf_err( nf90_def_var( grp, 'S', &
          NF90_DOUBLE, (/ nnzs_ /), var ) )
    call cdf_err( nf90_put_att( grp, var, &
          'info', 'Overlap matrix' ) )

    call cdf_err( nf90_def_var( grp, 'S_gradient', &
          NF90_DOUBLE, (/ xyz, nnzs_ /), var ) )
    call cdf_err( nf90_put_att( grp, var, &
          'unit', '1/Bohr' ) )
    call cdf_err( nf90_put_att( grp, var, &
          'info', 'Overlap matrix gradient' ) )

    call cdf_err( nf90_def_var( grp, 'DM', &
          NF90_DOUBLE, (/ nnzs_, spin_ /), var ) )
    call cdf_err( nf90_put_att( grp, var, &
          'info', 'Density matrix' ) )

    call cdf_err( nf90_def_var( grp, 'EDM', &
          NF90_DOUBLE, (/ nnzs_, spin_EDM /), var ) )
    call cdf_err( nf90_put_att( grp, var, &
          'unit', 'Ry' ) )
    call cdf_err( nf90_put_att( grp, var, &
          'info', 'Energy density matrix' ) )

    if ( savehs .or. TS_HS_save ) then
      ! Unit is already present in dictionary
      call cdf_err( nf90_def_var( grp, 'H', &
            NF90_DOUBLE, (/ nnzs_, spin_ /), var ) )
      call cdf_err( nf90_put_att( grp, var, &
          'unit', 'Ry' ) )
      call cdf_err( nf90_put_att( grp, var, &
            'info', 'Hamiltonian' ) )
    end if

    ! Even though I think we could do without, I add the
    ! xij array to the file
    ! Note that xij can be re-created using
    !    nsc, isc_off and xa
    !    if ( .not. Gamma ) then
    !       dic = dic//('info'.kv. &
    !            'Distance between orbital i and j')
    !       dic = dic//('unit'.kv.'Bohr')
    !       call ncdf_def_var(grp,'xij',NF90_DOUBLE,(/'xyz ','nnzs'/), &
    !            compress_lvl=cdf_comp_lvl,atts=dic)
    !    end if


    ! Create grid group
  	call cdf_err( nf90_def_grp( ncid, 'GRID', grp ) )

    ! Note that there are now 2 spin dimensions
    ! 1. in the top level NC file and one in the GRID group
    ! This is because for spin-orbit coupling the grid and
    ! matrix dimensions are not the same (4 vs 8, respectively)
    call cdf_err( nf90_def_dim( grp, 'spin', spin%Grid, spin_Grid ) )

  	call cdf_err( nf90_def_grp( ncid, 'SETTINGS', grp ) )

    call cdf_err( nf90_def_var( grp, 'DMTolerance', &
          NF90_DOUBLE, (/ one /), var ) )
    call cdf_err( nf90_put_att( grp, var, &
          'info', 'Tolerance for converging the density matrix' ) )

    call cdf_err( nf90_def_var( grp, 'HTolerance', &
          NF90_DOUBLE, (/ one /), var ) )
    call cdf_err( nf90_put_att( grp, var, &
          'info', 'Tolerance for converging the Hamiltonian' ) )

    call cdf_err( nf90_def_var( grp, 'NetCharge', &
          NF90_DOUBLE, (/ one /), var ) )
    call cdf_err( nf90_put_att( grp, var, &
          'info', 'Net charge of the system' ) )

    call cdf_err( nf90_def_var( grp, 'MixingWeight', &
          NF90_DOUBLE, (/ one /), var ) )
    call cdf_err( nf90_put_att( grp, var, &
          'info', 'Mixing weight' ) )

    call cdf_err( nf90_def_var( grp, 'BZ', &
          NF90_INT, (/ xyz, xyz /), var ) )
    call cdf_err( nf90_put_att( grp, var, &
          'info', 'Grid used for the Brillouin zone integration' ) )

    call cdf_err( nf90_def_var( grp, 'BZ_displ', &
          NF90_DOUBLE, (/ xyz /), var ) )
    call cdf_err( nf90_put_att( grp, var, 'unit', 'b**-1' ) )
    call cdf_err( nf90_put_att( grp, var, &
          'info', 'Grid displacement used in Brillouin zone' ) )

    call cdf_err( nf90_def_var( grp, 'ElectronicTemperature', &
          NF90_DOUBLE, (/ one /), var ) )
    call cdf_err( nf90_put_att( grp, var, 'unit', 'Ry' ) )
    call cdf_err( nf90_put_att( grp, var, &
          'info', 'Temperature for electrons' ) )

    call cdf_err( nf90_def_var( grp, 'MeshCutoff', &
          NF90_DOUBLE, (/ one /), var ) )
    call cdf_err( nf90_put_att( grp, var, 'unit', 'Ry' ) )
    call cdf_err( nf90_put_att( grp, var, &
          'info', 'Mesh cutoff for real space grid' ) )

    ! Create FC group
    if ( lis_FC ) then
      call cdf_err( nf90_def_grp( ncid, 'FC', grp ) )

      ! Create dimension arrays
      call cdf_err( nf90_def_dim( grp, 'na_fc', ia2-ia1+1, na_fc ) )
      call cdf_err( nf90_def_dim( grp, 'm_p', 2, m_p ) )

      call cdf_err( nf90_def_var( grp, 'disp', &
            NF90_DOUBLE, (/ one /), vdisp ) )
      call cdf_err( nf90_put_att( grp, vdisp, 'unit', 'Bohr' ) )
      call cdf_err( nf90_put_att( grp, vdisp, &
            'info', 'Displacement length' ) )

      call cdf_err( nf90_def_var( grp, 'ia_fc', &
            NF90_INT, (/ na_fc /), via_fc ) )
      call cdf_err( nf90_put_att( grp, via_fc, &
            'info', 'Displaced atoms' ) )

      call cdf_err( nf90_def_var( grp, 'fa0', &
            NF90_DOUBLE, (/ xyz, na_fc /), var ) )
      call cdf_err( nf90_put_att( grp, var, 'unit', 'Ry/Bohr' ) )
      call cdf_err( nf90_put_att( grp, var, &
            'info', 'Undisplaced atomic forces' ) )

      call cdf_err( nf90_def_var( grp, 'fa', &
            NF90_DOUBLE, (/ xyz, na_fc, m_p, xyz, na_fc /), var ) )
      call cdf_err( nf90_put_att( grp, var, 'unit', 'Ry/Bohr' ) )
      call cdf_err( nf90_put_att( grp, var, &
            'info', 'Displaced atomic forces' ) )
    endif

    if ( isolve == SOLVE_TRANSI ) then
      ! Save all information about the transiesta method
      call cdf_err( nf90_def_grp( ncid, 'TRANSIESTA', grp_tr ) )

      call cdf_err( nf90_def_var( grp_tr, 'BZ', &
            NF90_INT, (/ xyz, xyz /), var ) )
      call cdf_err( nf90_put_att( grp_tr, var, &
            'info', 'Grid used for the Brillouin zone integration' ) )

      call cdf_err( nf90_def_var( grp_tr, 'BZ_displ', &
            NF90_DOUBLE, (/ xyz /), var ) )
      call cdf_err( nf90_put_att( grp_tr, var, 'unit', 'b**-1' ) )
      call cdf_err( nf90_put_att( grp_tr, var, &
            'info', 'Grid displacement used in Brillouin zone' ) )

      call cdf_err( nf90_def_var( grp_tr, 'Volt', &
            NF90_DOUBLE, (/ one /), vVolt ) )
      call cdf_err( nf90_put_att( grp_tr, vVolt, 'unit', 'Ry' ) )
      call cdf_err( nf90_put_att( grp_tr, vVolt, &
            'info', 'Applied voltage' ) )

      ! Add all the electrodes
      do iEl = 1, N_Elec
        call cdf_err( nf90_def_grp( grp_tr, trim(Elecs(iEl)%name), grp2 ) )
        tmp = Elecs(iEl)%device_atoms()
        call cdf_err( nf90_def_dim( grp2, 'na', tmp, na ) )

        call cdf_err( nf90_def_var( grp2, 'a_idx', &
              NF90_INT, (/ na /), var ) )
        call cdf_err( nf90_put_att( grp2, var, &
              'info', 'Atoms belonging to electrode' ) )

        call cdf_err( nf90_def_var( grp2, 'mu', &
              NF90_DOUBLE, (/ one /), var ) )
        call cdf_err( nf90_put_att( grp2, var, 'unit', 'Ry' ) )
        call cdf_err( nf90_put_att( grp2, var, &
              'info', 'Chemical potential' ) )

        call cdf_err( nf90_def_var( grp2, 'kT', &
              NF90_DOUBLE, (/ one /), var ) )
        call cdf_err( nf90_put_att( grp2, var, 'unit', 'Ry' ) )
        call cdf_err( nf90_put_att( grp2, var, &
              'info', 'Electronic temperature' ) )
      end do
    end if

    ! Save all things necessary here
    if ( isolve == SOLVE_DIAGON ) then
  	  call cdf_err( nf90_put_att( ncid, NF90_GLOBAL, &
            'method', 'diagon' ) )
    else if ( isolve == SOLVE_ORDERN ) then
  	  call cdf_err( nf90_put_att( ncid, NF90_GLOBAL, &
            'method', 'order-n' ) )
    else if ( isolve == SOLVE_MINIM ) then
  	  call cdf_err( nf90_put_att( ncid, NF90_GLOBAL, &
            'method', 'omm' ) )
    else if ( isolve == SOLVE_TRANSI ) then
  	  call cdf_err( nf90_put_att( ncid, NF90_GLOBAL, &
            'method', 'transiesta' ) )
    else if ( isolve == SOLVE_PEXSI ) then
  	  call cdf_err( nf90_put_att( ncid, NF90_GLOBAL, &
            'method', 'pexsi' ) )
    end if
  	call cdf_err( nf90_put_att( ncid, NF90_GLOBAL, &
          'name', trim(sname) ) )
  	call cdf_err( nf90_put_att( ncid, NF90_GLOBAL, &
          'time', datestring() ) )
  	call cdf_err( nf90_put_att( ncid, NF90_GLOBAL, &
          'label', trim(slabel) ) )

    call cdf_err( nf90_enddef( ncid ) )

    call cdf_err( nf90_put_var( ncid, vlasto, lasto(1:na_u) ) )
    call cdf_err( nf90_put_var( ncid, vQtot, Qtot ) )
    if ( fixspin ) then
      call cdf_err( nf90_put_var( ncid, vQspin, total_spin ) )
    endif
    if ( lis_FC ) then
      call cdf_err( nf90_put_var( ncid, vdisp, dx ) )
      nullify(ibuf)
      call re_alloc( ibuf, 1, ia2-ia1+1, 'ibuf', 'cdf_init' )
      do i = ia1, ia2
        ibuf(i-ia1+1) = i
      end do
      call cdf_err( nf90_put_var( ncid, via_fc, ibuf ) )
      call de_alloc( ibuf, 'ibuf', 'cdf_init' )
    endif

    if ( isolve == SOLVE_TRANSI ) then
      ! Save all information about the transiesta method
      call cdf_err( nf90_put_var( grp_tr, vVolt, Volt ) )

      ! Add all the electrodes
      do iEl = 1 , N_Elec
        ! Inquire the group ID
        call cdf_err( nf90_inq_ncid( grp_tr,trim(Elecs(iEl)%name), grp2 ) )
        tmp = Elecs(iEl)%device_atoms()
        nullify(ibuf)
        call re_alloc( ibuf, 1, tmp, 'ibuf', 'cdf_init' )
        do i = 1 , tmp
          ibuf(i) = Elecs(iEl)%idx_a + i - 1
        end do
        call cdf_err( nf90_inq_varid( grp2, 'a_idx', var ) )
        call cdf_err( nf90_put_var( grp2, var, ibuf ) )
        call de_alloc( ibuf, 'ibuf', 'cdf_init' )

        call cdf_err( nf90_inq_varid( grp2, 'mu', var ) )
        call cdf_err( nf90_put_var( grp2, var, Elecs(iEl)%mu%mu ) )

        call cdf_err( nf90_inq_varid( grp2, 'kT', var ) )
        call cdf_err( nf90_put_var( grp2, var, Elecs(iEl)%mu%kT ) )
      end do
    end if
    ! Close the file
    call cdf_err( nf90_close( ncid ) )
  end subroutine cdf_init_file_new

  subroutine cdf_init_grid(fname, ntm)
    use fdf, only : fdf_get, leqi, fdf_defined
    use siesta_options, only: saverho, savedrho, savevh, savevna
    use siesta_options, only: savevt, savepsch, savetoch, saverhoxc
    use siesta_options, only: savebader
    use siesta_options, only: save_initial_charge_density

    character(len=*), intent(in) :: fname
    integer, intent(in) :: ntm(3)

    ! Local variables
    type(hNCDF) :: ncdf, grp
    integer :: prec, chks(3)
    character(len=64) :: key
    type(dictionary_t) :: dic

    ! We always re-write the file...
    call ncdf_open(ncdf,fname, &
        mode=ior(NF90_WRITE,NF90_NETCDF4))

    ! Create grid group
    call ncdf_open_grp(ncdf,'GRID', grp)

    call ncdf_def_dim(grp,'nx',ntm(1))
    call ncdf_def_dim(grp,'ny',ntm(2))
    call ncdf_def_dim(grp,'nz',ntm(3))

    ! Create the grid functions...

    ! all grids are using the grid_p precision
    if ( grid_p == dp ) then
      prec = NF90_DOUBLE
      ! In case the user thinks the double precision
      ! is needed, but only single precision
      ! is needed saving, we allow that.
      key = fdf_get('CDF.Grid.Precision','double')
      if ( leqi(key,'single') ) prec = NF90_FLOAT
      if ( leqi(key,'float') )  prec = NF90_FLOAT
    else
      ! The grid is in single precision, so
      ! we save it in that precision.
      prec = NF90_FLOAT
    end if
    
    chks = (/ntm(1),ntm(2),1/)
    
    ! Define the units for all charge variables
    dic = ('unit'.kv.'e/Bohr**3')

    if ( save_initial_charge_density ) then
      dic = dic//('info'.kv.'Initial charge density')
      call ncdf_def_var(grp,'RhoInit',prec,(/'nx  ','ny  ','nz  ','spin'/), &
          compress_lvl=cdf_comp_lvl,atts=dic,chunks=chks)
    end if

    if ( saverho ) then
      dic = dic//('info'.kv.'Charge density')
      call ncdf_def_var(grp,'Rho',prec,(/'nx  ','ny  ','nz  ','spin'/), &
          compress_lvl=cdf_comp_lvl,atts=dic,chunks=chks)
    end if

    if ( savepsch ) then
      dic = dic//('info'.kv.'Diffuse ionic charge')
      call ncdf_def_var(grp,'Chlocal',prec,(/'nx','ny','nz'/), &
          compress_lvl=cdf_comp_lvl,atts=dic,chunks=chks)
    end if

    if ( savetoch ) then
      dic = dic//('info'.kv.'Total charge')
      call ncdf_def_var(grp,'RhoTot',prec,(/'nx','ny','nz'/), &
          compress_lvl=cdf_comp_lvl,atts=dic,chunks=chks)
    end if

    if ( fdf_defined("LocalDensityOfStates") ) then
      dic = dic//('info'.kv.'Local Density of States')
      call ncdf_def_var(grp,'LDOS',prec,(/'nx  ','ny  ','nz  ','spin'/), &
          compress_lvl=cdf_comp_lvl,atts=dic,chunks=chks)
    end if

    if ( savedrho ) then
      dic = dic//('info'.kv.'Density difference from atomic densities')
      call ncdf_def_var(grp,'RhoDelta',prec,(/'nx  ','ny  ','nz  ','spin'/), &
          compress_lvl=cdf_comp_lvl,atts=dic,chunks=chks)
    end if

    if ( saverhoxc ) then
      dic = dic//('info'.kv.'Density used to calculate XC functional')
      call ncdf_def_var(grp,'RhoXC',prec,(/'nx  ','ny  ','nz  ','spin'/), &
          compress_lvl=cdf_comp_lvl,atts=dic,chunks=chks)
    end if

    if ( savebader ) then
      dic = dic//('info'.kv.'Bader charge')
      call ncdf_def_var(grp,'RhoBader',prec,(/'nx','ny','nz'/), &
          compress_lvl=cdf_comp_lvl,atts=dic,chunks=chks)
    end if

    ! Define the units for the potential grids
    dic = dic//('unit'.kv.'Ry')

    if ( savevna ) then
      dic = dic//('info'.kv.'Neutral atom potential')
      call ncdf_def_var(grp,'Vna',prec,(/'nx','ny','nz'/), &
          compress_lvl=cdf_comp_lvl,atts=dic,chunks=chks)
    end if

    if ( savevh ) then
      dic = dic//('info'.kv.'Electrostatic potential')
      call ncdf_def_var(grp,'Vh',prec,(/'nx','ny','nz'/), &
          compress_lvl=cdf_comp_lvl,atts=dic,chunks=chks)
    end if

    if ( savevt ) then
      dic = dic//('info'.kv.'Total potential')
      call ncdf_def_var(grp,'Vt',prec,(/'nx  ','ny  ','nz  ','spin'/), &
          compress_lvl=cdf_comp_lvl,atts=dic,chunks=chks)
    end if

    call delete(dic)

    ! Close the file
    call ncdf_close(ncdf)

  end subroutine cdf_init_grid

  subroutine cdf_init_grid_new(fname, ntm)
    use fdf, only : fdf_get, leqi
    use siesta_options, only: saverho, savedrho, savevh, savevna
    use siesta_options, only: savevt, savepsch, savetoch, saverhoxc
    use siesta_options, only: savebader
    use siesta_options, only: save_initial_charge_density
    implicit none
    character(len=*), intent(in) :: fname
    integer, intent(in) :: ntm(3)
    integer :: cmode, ncid, grp, nx_, ny_, nz_, spin_, prec, var, i
    character(len=64) :: key
    call cdfGrid%init_ng( ntm )

    if (node/=0) return
    ! We always re-write the file...
    cmode = IOR(NF90_WRITE,NF90_NETCDF4)
    call cdf_err( nf90_open( fname, cmode, ncid ) )

    ! Open GRID group
    call cdf_err( nf90_inq_ncid( ncid, 'GRID', grp ) )
    call cdf_err( nf90_inq_dimid( grp, 'spin', spin_ ) )

    call cdf_err( nf90_def_dim( grp, 'nx', ntm(1), nx_ ) )
    call cdf_err( nf90_def_dim( grp, 'ny', ntm(2), ny_ ) )
    call cdf_err( nf90_def_dim( grp, 'nz', ntm(3), nz_ ) )

    ! Create the grid functions...
    ! all grids are using the grid_p precision
    if ( grid_p == dp ) then
      prec = NF90_DOUBLE
      ! In case the user thinks the double precision
      ! is needed, but only single precision
      ! is needed saving, we allow that.
      key = fdf_get('CDF.Grid.Precision','double')
      if ( leqi(key,'single') ) prec = NF90_FLOAT
      if ( leqi(key,'float') )  prec = NF90_FLOAT
    else
      ! The grid is in single precision, so
      ! we save it in that precision.
      prec = NF90_FLOAT
    end if

    if ( save_initial_charge_density ) then
      call cdf_err( nf90_def_var( grp, 'RhoInit', &
          prec, (/ nx_, ny_, nz_, spin_ /), var ) )
      call cdf_err( nf90_put_att( grp, var, 'unit', 'e/Bohr**3' ) )
      call cdf_err( nf90_put_att( grp, var, &
          'info', 'Initial charge density' ) )
    end if

    if ( saverho ) then
      call cdf_err( nf90_def_var( grp, 'Rho', &
          prec, (/ nx_, ny_, nz_, spin_ /), var ) )
      call cdf_err( nf90_put_att( grp, var, 'unit', 'e/Bohr**3' ) )
      call cdf_err( nf90_put_att( grp, var, &
          'info', 'Charge density' ) )
    end if

    if ( savepsch ) then
      call cdf_err( nf90_def_var( grp, 'Chlocal', &
          prec, (/ nx_, ny_, nz_ /), var ) )
      call cdf_err( nf90_put_att( grp, var, 'unit', 'e/Bohr**3' ) )
      call cdf_err( nf90_put_att( grp, var, &
          'info', 'Diffuse ionic charge' ) )
    end if

    if ( savetoch ) then
      call cdf_err( nf90_def_var( grp, 'RhoTot', &
          prec, (/ nx_, ny_, nz_ /), var ) )
      call cdf_err( nf90_put_att( grp, var, 'unit', 'e/Bohr**3' ) )
      call cdf_err( nf90_put_att( grp, var, &
          'info', 'Total charge' ) )
    end if

    if ( savedrho ) then
      call cdf_err( nf90_def_var( grp, 'RhoDelta', &
          prec, (/ nx_, ny_, nz_, spin_ /), var ) )
      call cdf_err( nf90_put_att( grp, var, 'unit', 'e/Bohr**3' ) )
      call cdf_err( nf90_put_att( grp, var, &
          'info', 'Density difference from atomic densities' ) )
    end if

    if ( saverhoxc ) then
      call cdf_err( nf90_def_var( grp, 'RhoXC', &
          prec, (/ nx_, ny_, nz_, spin_ /), var ) )
      call cdf_err( nf90_put_att( grp, var, 'unit', 'e/Bohr**3' ) )
      call cdf_err( nf90_put_att( grp, var, &
          'info', 'Density used to calculate XC functional' ) )
    end if

    if ( savebader ) then
      call cdf_err( nf90_def_var( grp, 'RhoBader', &
          prec, (/ nx_, ny_, nz_ /), var ) )
      call cdf_err( nf90_put_att( grp, var, 'unit', 'e/Bohr**3' ) )
      call cdf_err( nf90_put_att( grp, var, &
          'info', 'Bader charge' ) )
    end if

    if ( savevna ) then
      call cdf_err( nf90_def_var( grp, 'Vna', &
          prec, (/ nx_, ny_, nz_ /), var ) )
      call cdf_err( nf90_put_att( grp, var, 'unit', 'Ry' ) )
      call cdf_err( nf90_put_att( grp, var, &
          'info', 'Neutral atom potential' ) )
    end if

    if ( savevh ) then
      call cdf_err( nf90_def_var( grp, 'Vh', &
          prec, (/ nx_, ny_, nz_ /), var ) )
      call cdf_err( nf90_put_att( grp, var, 'unit', 'Ry' ) )
      call cdf_err( nf90_put_att( grp, var, &
          'info', 'Electrostatic potential' ) )
    end if

    if ( savevt ) then
      call cdf_err( nf90_def_var( grp, 'Vt', &
          prec, (/ nx_, ny_, nz_, spin_ /), var ) )
      call cdf_err( nf90_put_att( grp, var, 'unit', 'Ry' ) )
      call cdf_err( nf90_put_att( grp, var, &
          'info', 'Total potential' ) )
    end if

    ! Close the file
    call cdf_err( nf90_close( ncid ) )

  end subroutine cdf_init_grid_new

  subroutine cdf_save_settings(fname)

    use kpoint_scf_m,   only: kpoint_scf
    use siesta_options, only: dDtol, dHtol, charnet, wmix, temp, g2cut
    use siesta_options, only: isolve
    use siesta_options, only: SOLVE_DIAGON, SOLVE_ORDERN, SOLVE_TRANSI
    use ts_kpoint_scf_m, only: ts_kpoint_scf

    character(len=*), intent(in) :: fname
    
    type(hNCDF) :: ncdf, grp

    ! We just open it (prepending)
    call ncdf_open(ncdf,fname, &
        mode=ior(NF90_WRITE,NF90_NETCDF4))

    call ncdf_open_grp(ncdf,'SETTINGS',grp)

    ! Save settings
    call ncdf_put_var(grp,'BZ',kpoint_scf%k_cell)
    call ncdf_put_var(grp,'BZ_displ',kpoint_scf%k_displ)
    call ncdf_put_var(grp,'DMTolerance',dDtol)
    call ncdf_put_var(grp,'HTolerance',dHtol)
    call ncdf_put_var(grp,'NetCharge',charnet)
    call ncdf_put_var(grp,'MixingWeight',wmix)
    call ncdf_put_var(grp,'ElectronicTemperature',Temp)
    call ncdf_put_var(grp,'MeshCutoff',g2cut)

    ! I suggest that all MD settings are 
    ! put in the MD-group

    if ( isolve == SOLVE_TRANSI ) then
       call ncdf_open_grp(ncdf,'TRANSIESTA',grp)

       call ncdf_put_var(grp,'BZ',ts_kpoint_scf%k_cell)
       call ncdf_put_var(grp,'BZ_displ',ts_kpoint_scf%k_displ)

    end if

    call ncdf_close(ncdf)

  end subroutine cdf_save_settings

  subroutine cdf_save_settings_new(fname)
    use kpoint_scf_m,   only: kpoint_scf
    use siesta_options, only: dDtol, dHtol, charnet, wmix, temp, g2cut
    use siesta_options, only: isolve
    use siesta_options, only: SOLVE_DIAGON, SOLVE_ORDERN, SOLVE_TRANSI
    use ts_kpoint_scf_m, only: ts_kpoint_scf
    implicit none
    character(len=*), intent(in) :: fname
    integer :: ncid, grp, cmode, var
    integer :: MPIerror
    if ( Node /= 0 ) return
    ! We just open it (prepending)
    cmode = IOR(NF90_WRITE,NF90_NETCDF4)
    call cdf_err( nf90_open( fname, cmode, ncid ) )

    ! Open SETTINGS group
    call cdf_err( nf90_inq_ncid( ncid, 'SETTINGS', grp ) )

    ! Save settings
    call cdf_err( nf90_inq_varid( grp, 'BZ', var ) )
    call cdf_err( nf90_put_var( grp, var, kpoint_scf%k_cell ) )

    call cdf_err( nf90_inq_varid( grp, 'BZ_displ', var ) )
    call cdf_err( nf90_put_var( grp, var, kpoint_scf%k_displ ) )

    call cdf_err( nf90_inq_varid( grp, 'DMTolerance', var ) )
    call cdf_err( nf90_put_var( grp, var, dDtol ) )

    call cdf_err( nf90_inq_varid( grp, 'HTolerance', var ) )
    call cdf_err( nf90_put_var( grp, var, dHtol ) )

    call cdf_err( nf90_inq_varid( grp, 'NetCharge', var ) )
    call cdf_err( nf90_put_var( grp, var, charnet ) )

    call cdf_err( nf90_inq_varid( grp, 'MixingWeight', var ) )
    call cdf_err( nf90_put_var( grp, var, wmix ) )

    call cdf_err( nf90_inq_varid( grp, 'ElectronicTemperature', var ) )
    call cdf_err( nf90_put_var( grp, var, Temp ) )

    call cdf_err( nf90_inq_varid( grp, 'MeshCutoff', var ) )
    call cdf_err( nf90_put_var( grp, var, g2cut ) )

    ! I suggest that all MD settings are 
    ! put in the MD-group
    if ( isolve == SOLVE_TRANSI ) then
      ! Open TRANSIESTA group
      call cdf_err( nf90_inq_ncid( ncid, 'TRANSIESTA', grp ) )
      call cdf_err( nf90_inq_varid( grp, 'BZ', var ) )
      call cdf_err( nf90_put_var( grp, var, ts_kpoint_scf%k_cell ) )
      call cdf_err( nf90_inq_varid( grp, 'BZ_displ', var ) )
      call cdf_err( nf90_put_var( grp, var, ts_kpoint_scf%k_displ ) )
    end if
    call cdf_err( nf90_close( ncid ) )
  end subroutine cdf_save_settings_new

  subroutine cdf_save_state(fname,dic_save)
    !    use m_gamma, only : Gamma
    use m_energies, only: Ef, Efs, NEGF_Vha
    use atomlist, only : Qtot
    use siesta_options, only: fixspin, total_spin, isolve, SOLVE_TRANSI
    use siesta_geom, only: na_u, ucell, xa, va
    use siesta_geom, only: nsc, isc_off
    use class_Sparsity, only: nrows_g
    use sparse_matrices, only: sparse_pattern, block_dist
    use sparse_matrices, only: S_1D, xij_2D
    use sparse_matrices, only: DM_2D, EDM_2D, H_2D
    use sparse_matrices, only: gradS_2D
    use m_stress, only : stress
    use m_forces, only: fa

    character(len=*), intent(in) :: fname
    ! Dictionary containing keys that we will save
    type(dictionary_t), intent(in) :: dic_save
    type(hNCDF) :: ncdf, grp
    integer :: no
    integer, allocatable :: gncol(:)

    call timer('CDF',1)

    ! We just open it (prepending)
#ifdef MPI
    if ( cdf_w_parallel ) then
      call ncdf_open(ncdf,fname, &
          mode=ior(NF90_WRITE,NF90_MPIIO), comm=MPI_Comm_World)
    else
#endif
      call ncdf_open(ncdf,fname,mode=ior(NF90_WRITE,NF90_NETCDF4))
#ifdef MPI
    end if
#endif

    ! Add attribute of current Fermi-level
    if ( fixspin ) then
      if ( ('Ef' .in. dic_save) .and. Node == 0 ) &
          call ncdf_put_var(ncdf,'Ef',Efs)
      if ( ('Qspin' .in. dic_save) .and. Node == 0 ) &
          call ncdf_put_var(ncdf,'Qspin',total_spin)
    else
      if ( ('Ef' .in. dic_save) .and. Node == 0 ) &
          call ncdf_put_var(ncdf,'Ef',Ef)
    end if
    if ( ('Qtot' .in. dic_save) .and. Node == 0 ) &
        call ncdf_put_var(ncdf,'Qtot',Qtot)
    ! Save nsc, xa, fa, lasto, ucell
    if ( ('nsc' .in. dic_save) .and. Node == 0 ) &
        call ncdf_put_var(ncdf,'nsc',nsc)
    if ( ('cell' .in. dic_save) .and. Node == 0 ) &
        call ncdf_put_var(ncdf,'cell',ucell)
    if ( ('xa' .in. dic_save) .and. Node == 0 ) &
        call ncdf_put_var(ncdf,'xa',xa(:,1:na_u))
    if ( ('va' .in. dic_save) .and. Node == 0 ) &
        call ncdf_put_var(ncdf,'va',va(:,1:na_u))
    if ( ('fa' .in. dic_save) .and. Node == 0 ) &
        call ncdf_put_var(ncdf,'fa',fa(:,1:na_u))
    if ( ('stress' .in. dic_save) .and. Node == 0 ) &
        call ncdf_put_var(ncdf,'stress',stress)

    ! Sparsity format
    call ncdf_open_grp(ncdf,'SPARSE',grp)

    no = nrows_g(sparse_pattern)
    allocate(gncol(no))
    ! Signal that it needs to be filled (first element is negative)
    gncol(1) = -1

    if ( 'isc_off'.in. dic_save) &
        call ncdf_put_var(grp,'isc_off',isc_off)
    if ( 'sp' .in. dic_save ) &
        call cdf_w_Sp(grp, sparse_pattern, block_dist, gncol=gncol)
    if ( 'S' .in. dic_save ) &
        call cdf_w_d1D(grp,'S',S_1D, gncol=gncol)
    if ( 'gradS' .in. dic_save ) &
        call cdf_w_d2D(grp,'S_gradient', gradS_2D, gncol=gncol)
    ! we don't have a flag for saving xij
    !if ( 'xij' .in. dic_save ) then
    !  call cdf_w_d2D(grp,'xij', xij_2D, gncol=gncol)
    !end if
    if ( 'H' .in. dic_save ) &
        call cdf_w_d2D(grp,'H',H_2D, gncol=gncol)
    if ( 'DM' .in. dic_save ) &
        call cdf_w_d2D(grp,'DM',DM_2D, gncol=gncol)
    if ( 'EDM' .in. dic_save ) &
        call cdf_w_d2D(grp,'EDM',EDM_2D, gncol=gncol)
    
    if ( isolve == SOLVE_TRANSI ) then
      call ncdf_open_grp(ncdf,'TRANSIESTA',grp)

      if ( 'NEGF_Vha' .in. dic_save ) &
          call ncdf_put_var(grp,'dHartree',NEGF_Vha)
    end if

    deallocate(gncol)

    call ncdf_close(ncdf)

    call timer('CDF',2)
  end subroutine cdf_save_state

  subroutine cdf_save_state_new(fname,dic_save)
    !    use m_gamma, only : Gamma
    use m_energies, only: Ef, Efs
    use atomlist, only : Qtot
    use siesta_options, only: fixspin, total_spin
    use siesta_geom, only: na_u, ucell, xa, va
    use siesta_geom, only: nsc, isc_off
    use class_Sparsity, only: nrows_g
    use sparse_matrices, only: sparse_pattern, block_dist
    use sparse_matrices, only: S_1D, xij_2D
    use sparse_matrices, only: DM_2D, EDM_2D, H_2D
    use sparse_matrices, only: gradS_2D
    use m_stress, only : stress
    use m_forces, only: fa
    implicit none
    character(len=*), intent(in) :: fname
    type(dictionary_t), intent(in) :: dic_save
    integer :: cmode, ncid, var, grp, no, i
    type(t_ncdfSparse) :: ncdfSp
    integer :: max_dsize
    call timer('CDF',1)
    max_dsize = 0
    if ( 'sp' .in. dic_save ) max_dsize = max(max_dsize,4)
    if ( ('S' .in. dic_save) .or. ('H' .in. dic_save) .or. &
         ('DM' .in. dic_save) .or. ('EDM' .in. dic_save) ) &
         max_dsize = max(max_dsize,dp)
    if ( 'gradS' .in. dic_save ) max_dsize = max(max_dsize,dp*3)
    call ncdfSp%init_Sp( sparse_pattern, block_dist, max_dsize )

    cmode = IOR(NF90_WRITE,NF90_NETCDF4)
    if (ncdfSp%nwriters==1) then
      if (Node==0) call cdf_err( nf90_open( fname, cmode, ncid ) )
#ifdef MPI
    else
      cmode = IOR(cmode,NF90_MPIIO)
      call cdf_err( nf90_open_par( fname, cmode, MPI_COMM_WORLD, &
     &               MPI_INFO_NULL, ncid ) )
#endif
    endif
    ! We just open it (prepending)
    if (node==0) then
      ! Add attribute of current Fermi-level
      if ( fixspin ) then
        if ('Ef' .in. dic_save ) then
          call cdf_err( nf90_inq_varid( ncid, 'Ef', var ) )
          call cdf_err( nf90_put_var( ncid, var, Efs ) )
        endif
        if ('Qspin' .in. dic_save ) then
          call cdf_err( nf90_inq_varid( ncid, 'Qspin', var ) )
          call cdf_err( nf90_put_var( ncid, var, total_spin ) )
        endif
      else
        if ('Ef' .in. dic_save ) then
          call cdf_err( nf90_inq_varid( ncid, 'Ef', var ) )
          call cdf_err( nf90_put_var( ncid, var, Ef ) )
        endif
      end if
      if ('Qtot' .in. dic_save ) then
        call cdf_err( nf90_inq_varid( ncid, 'Qtot', var ) )
        call cdf_err( nf90_put_var( ncid, var, Qtot ) )
      endif
      if ('nsc' .in. dic_save ) then
        call cdf_err( nf90_inq_varid( ncid, 'nsc', var ) )
        call cdf_err( nf90_put_var( ncid, var, nsc ) )
      endif
      if ('cell' .in. dic_save ) then
        call cdf_err( nf90_inq_varid( ncid, 'cell', var ) )
        call cdf_err( nf90_put_var( ncid, var, ucell ) )
      endif
      if ('xa' .in. dic_save ) then
        call cdf_err( nf90_inq_varid( ncid, 'xa', var ) )
        call cdf_err( nf90_put_var( ncid, var, xa(:,1:na_u) ) )
      endif
      if ('va' .in. dic_save ) then
        call cdf_err( nf90_inq_varid( ncid, 'va', var ) )
        call cdf_err( nf90_put_var( ncid, var, va(:,1:na_u) ) )
      endif
      if ('fa' .in. dic_save ) then
        call cdf_err( nf90_inq_varid( ncid, 'fa', var ) )
        call cdf_err( nf90_put_var( ncid, var, fa(:,1:na_u) ) )
      endif
      if ('stress' .in. dic_save ) then
        call cdf_err( nf90_inq_varid( ncid, 'stress', var ) )
        call cdf_err( nf90_put_var( ncid, var, stress ) )
      endif
    endif
    grp = -1
    if (ncdfSp%IOnode) &
        call cdf_err( nf90_inq_ncid( ncid, 'SPARSE', grp ) )
    if (('isc_off' .in. dic_save).and. Node==0) then
      call cdf_err( nf90_inq_varid( grp, 'isc_off', var ) )
      call cdf_err( nf90_put_var( grp, var, isc_off ) )
    endif
    if ( 'sp' .in. dic_save ) &
      call ncdfSp%write_Sp( grp, sparse_pattern )
    if ( 'S' .in. dic_save ) &
      call ncdfSp%write_d1D( grp, 'S', S_1D )
    if ( 'gradS' .in. dic_save ) &
      call ncdfSp%write_d2D( grp, 'S_gradient', gradS_2D )
    if ( 'H' .in. dic_save ) &
      call ncdfSp%write_d2D( grp, 'H', H_2D )
    if ( 'DM' .in. dic_save ) &
      call ncdfSp%write_d2D( grp, 'DM', DM_2D )
    if ( 'EDM' .in. dic_save ) &
      call ncdfSp%write_d2D( grp, 'EDM', EDM_2D )

    if (ncdfSp%nwriters>1 .or. Node==0) call cdf_err( nf90_close( ncid ) )
    call ncdfSp%delete_Sp( )
    call timer('CDF',2)
  end subroutine cdf_save_state_new

  subroutine cdf_save_grid(fname,vname,nspin,nmeshl,grid)
    character(len=*), intent(in) :: fname, vname
    integer, intent(in) :: nspin, nmeshl(3)
    real(grid_p), intent(in) :: grid(product(nmeshl),nspin)

    type(hNCDF) :: ncdf
    integer :: is

    call timer('CDF-grid',1)

    ! We just open it (prepending)
#ifdef MPI
    if ( cdf_w_parallel ) then
      call ncdf_open(ncdf,fname, group='GRID', &
          mode=ior(NF90_WRITE,NF90_MPIIO), &
          comm=MPI_Comm_World)
    else
#endif
      call ncdf_open(ncdf,fname, group='GRID', &
          mode=ior(NF90_WRITE,NF90_NETCDF4))
#ifdef MPI
    end if
#endif

    ! Save the grid
    if ( nspin > 1 ) then
      do is = 1 , nspin 
        call cdf_w_grid(ncdf,vname,nmeshl,grid(:,is),idx=is)
      end do
    else
      call cdf_w_grid(ncdf,vname,nmeshl,grid(:,1))
    end if

    call ncdf_close(ncdf)

    call timer('CDF-grid',2)

  end subroutine cdf_save_grid

  subroutine cdf_save_grid_new( fname, vname, nspin, nmeshl, grid )
    !use moreMeshSubs, only : prepareMeshBuffers, reduceMeshToMaster
    use moreMeshSubs, only : getLocalBox
    implicit none
    character(len=*), intent(in) :: fname, vname
    integer, intent(in) :: nspin, nmeshl(3)
    real(grid_p), intent(in) :: grid(product(nmeshl),nspin)
#   if NEW_NCDF==5
    integer :: cmode, ierr, ncid, grp, var, lbox(2,3), fz
    integer :: is, i, o, nx, ny, xy, p, r, n, rq, gnx, gny, gxy, id, nz
    integer :: MPIstatus(MPI_STATUS_SIZE)
    integer, pointer :: srequest(:), rrequest(:)
    real(grid_p), pointer :: wbuff(:)

    call timer('CDF-grid',1)
    if (node==0) then
      cmode = IOR(NF90_WRITE,NF90_NETCDF4)
      call cdf_err( nf90_open( fname, cmode, ncid ) )
      call cdf_err( nf90_inq_ncid( ncid, 'GRID', grp ) )
      call cdf_err( nf90_inq_varid( grp, vname, var ) )
    endif
    if (nodes==1) then
      call cdf_err( nf90_put_var( grp, var, grid, &
          count=(/ cdfGrid%ntm, nspin /) ) )
    else
      nullify(srequest,rrequest,wbuff)
      allocate(srequest(cdfGrid%nsend),rrequest(cdfGrid%nrecv),wbuff(cdfGrid%bsize))
      call getLocalBox( 1, Node+1, lbox )
      do is=1, nspin
        nx = lbox(2,1)-lbox(1,1)+1
        ny = lbox(2,2)-lbox(1,2)+1
        xy = nx*ny
        o = 1
        rq = 1
        do r=1, cdfGrid%nred
          p = (r-1)*cdfGrid%grsize
          do i= 1, cdfGrid%scount(r)
            call MPI_Isend( grid(o,is), xy, MPI_grid_real, p, &
              0, MPI_Comm_World, srequest(rq), ierr )
            o = o+xy
            rq = rq+1
          enddo
        enddo
        if (cdfGrid%Iamred) then
          id = cdfGrid%redid
          fz = cdfGrid%ndist(id)
          nz = cdfGrid%ndist(id+1)-fz
          gnx = cdfGrid%ntm(1)
          gny = cdfGrid%ntm(2)
          gxy = gnx*gny
          rq = 1
          do p=1, Nodes
            if (cdfGrid%rcount(p)==0) cycle
            call getLocalBox( 1, p, lbox )
            nx = lbox(2,1)-lbox(1,1)+1
            ny = lbox(2,2)-lbox(1,2)+1
            o = 1+(lbox(1,2)-1)*gnx
            if (lbox(1,3)>fz) &
                o = o+gxy*(lbox(1,3)-fz)
            xy = nx*ny
            do i= 1, cdfGrid%rcount(p)
              call MPI_Irecv( wbuff(o), xy, MPI_grid_real, p-1, &
                  0, MPI_Comm_World, rrequest(rq), ierr )
              o = o+gxy
              rq = rq+1
            enddo
          enddo
        	call MPI_WAITALL( cdfGrid%nrecv, rrequest, MPI_STATUSES_IGNORE, ierr )
          if (node==0) then
            do r=1, cdfGrid%nred
              p = (r-1)*cdfGrid%grsize
              o = cdfGrid%ndist(r)
              nz = cdfGrid%ndist(r+1)-o
              if (nz==0) cycle
              if (p/=0) then
                call MPI_Recv( wbuff, nz*gxy, MPI_grid_real, &
                    p, 0, MPI_Comm_World, MPIstatus, ierr )
              endif
              call cdf_err( nf90_put_var( grp, var, wbuff, &
                  start=(/1,1,o,is/), count=(/gnx,gny,nz,1/) ) )
            enddo
          else
            call MPI_Send( wbuff, nz*gxy, MPI_grid_real, &
                    0, 0, MPI_Comm_World, ierr )
          endif
        endif
        call MPI_WAITALL( cdfGrid%nsend, srequest, MPI_STATUSES_IGNORE, ierr )
      enddo
      deallocate(srequest,rrequest,wbuff)
    endif

    if (node==0) call cdf_err( nf90_close( ncid ) )
# 	else  /* NEW_NCDF/=5 */
    integer :: cmode, ncid, grp, var, n, i, j, p, di, mo
    integer :: ierr, ndims, nx, ny, nz, lz, ix, iy, iz
    integer :: off, ns, lbox(2,3), nzoff
    integer :: nwriters, IamW, status(MPI_STATUS_SIZE)
    integer*8 :: mem
    integer, pointer :: wrtrs(:), nzs(:), sndnz(:), count(:), displ(:)
    real(grid_p), pointer :: rbuf(:), wbuf(:,:,:), wbuf2(:,:,:,:)
    logical :: init
    call timer('CDF-grid',1)
    call cdfGrid%getDims( nx, ny, nz )
    if (Nodes==1) then
      cmode = IOR(NF90_WRITE,NF90_NETCDF4)
      call cdf_err( nf90_open( fname, cmode, ncid ) )
      ! Open GRID group
      call cdf_err( nf90_inq_ncid( ncid, 'GRID', grp ) )
      call cdf_err( nf90_inq_varid( grp, vname, var ) )
      call cdf_err( nf90_inquire_variable( grp, var, ndims=ndims ) )
      call cdf_err( nf90_put_var( grp, var, grid, &
          count=(/ nx, ny, nz, nspin /) ) )
#ifdef TMP
      if (ndims==3) then
        call Reshape3D( grid, wbuf, nx, ny, nz )
        call cdf_err( nf90_put_var( grp, var, wbuf ) )
      else
        call Reshape4D( grid, wbuf2, nx, ny, nz, nspin )
        call cdf_err( nf90_put_var( grp, var, wbuf2 ) )
      endif
#endif
    else
      call cdfGrid%getArrays( nwriters, IamW, nzoff, wrtrs, nzs, &
                              sndnz, count, displ )
      cmode = IOR(NF90_WRITE,NF90_NETCDF4)
      if (nwriters==1) then
        if (Node==0) call cdf_err( nf90_open( fname, cmode, ncid ) )
      else
        cmode = IOR(cmode,NF90_MPIIO)
        call cdf_err( nf90_open_par( fname, cmode, MPI_COMM_WORLD, &
                      MPI_INFO_NULL, ncid ) )
      endif
      nullify(rbuf,wbuf)
      if (IamW>0) then
        lz = nzs(IamW)
        n = nx*ny*lz
        call re_alloc( rbuf, 1, n, 'rbuf', 'cdf_grid' )
        call re_alloc( wbuf, 1, nx, 1, ny, 1, lz, 'wbuf', 'cdf_grid' )
        ! Open GRID group
        call cdf_err( nf90_inq_ncid( ncid, 'GRID', grp ) )
        call cdf_err( nf90_inq_varid( grp, vname, var ) )
        call cdf_err( nf90_inquire_variable( grp, var, ndims=ndims ) )
      else
        call re_alloc( rbuf, 1, 1, 'rbuf', 'cdf_grid' )
        call re_alloc( wbuf, 1, 1, 1, 1, 1, 1, 'wbuf', 'cdf_grid' )
      endif

      do ns=1, nspin
        off = 1
        do i=1, nwriters
          p = wrtrs(i)
          call MPI_Gatherv( grid(off,ns), sndnz(i), MPI_grid_real, rbuf, &
              count, displ, MPI_grid_real, p, MPI_COMM_WORLD, ierr )
          off = off + sndnz(i)
        enddo
        if (IamW>0) then
          off = 1
          do i=1, nodes
            if (count(i)==0) cycle
            call getLocalBox( 1, i, lbox ) ! Distribution == 1
            lbox(:,3) = (/ max(lbox(1,3)-nzoff,1), &
                           min(lbox(2,3)-nzoff,nzs(IamW)) /)
            do iz= lbox(1,3), lbox(2,3)
              do iy= lbox(1,2), lbox(2,2)
                do ix= lbox(1,1), lbox(2,1)
                  wbuf(ix,iy,iz) = rbuf(off)
                  off = off + 1
                enddo
              enddo
            enddo
          enddo
          if (ndims==3) then
            call cdf_err( nf90_put_var( grp, var, wbuf, &
                start=(/1,1,nzoff+1/), count=(/nx,ny,nzs(IamW)/)))
          else
            call cdf_err( nf90_put_var( grp, var, wbuf, &
                start=(/1,1,nzoff+1,ns/), count=(/nx,ny,nzs(IamW),1/)))
          endif
        endif
      enddo
      call de_alloc( wbuf, 'wbuf', 'cdf_grid' )
      call de_alloc( rbuf, 'rbuf', 'cdf_grid' )
      ! Close the file
      if (nwriters>1 .or. Node==0) call cdf_err( nf90_close( ncid ) )
    endif
#	  endif /* NEW_NCDF/=5 */
    call timer('CDF-grid',2)

  contains
    subroutine Reshape3D( x, p, n1, n2, n3 )
    integer::         n1, n2, n3
    real(grid_p), target::    x(n1,n2,n3)
    real(grid_p), pointer::   p(:,:,:)
    p => x
    end subroutine Reshape3D

    subroutine Reshape4D( x, p, n1, n2, n3, nn )
    integer::         n1, n2, n3, nn
    real(grid_p), target::    x(n1,n2,n3,nn)
    real(grid_p), pointer::   p(:,:,:,:)
    p => x
    end subroutine Reshape4D

  end subroutine cdf_save_grid_new


  subroutine cdf_save_basis(fname)

    use siesta_geom, only: na_u, isa

    use atmparams, only : nt => NTBMAX
    use atm_types, only : species_info, species, nspecies
    use radial, only : rad_func

    character(len=*), intent(in) :: fname

    ! Local variables
    type(species_info), pointer :: spp
    type(rad_func), pointer :: p

    type(hNCDF) :: nf, ncdf, grp
    type(dictionary_t) :: dic, d
    character(len=DICTIONARY_KEY_LENGTH) :: key
    type(variable_t) :: v
    integer :: is, i

    ! Used for saving variables
    integer :: no, nk
    integer, allocatable :: iaux(:)

    ! Unluckily the new basis saves only
    ! saved on the IO node. 
    ! No other node must therefore access this routine
    if ( Node /= 0 ) return

    call timer('CDF-basis',1)

    call ncdf_open(nf,fname,mode=ior(NF90_WRITE,NF90_NETCDF4))

    ! create the BASIS group
    call ncdf_def_grp(nf,'BASIS',ncdf)

    ! Create a list of the species associated with each atom
    dic = ('info'.kv.'Basis of each atom by ID')
    call ncdf_def_var(ncdf,'basis',NF90_INT,(/'na_u'/), atts= dic)
    call delete(dic)
    call ncdf_put_var(ncdf,'basis',isa(1:na_u))

    do is = 1 , nspecies

      ! Get current specie
      spp => species(is)

      ! Get array sizes
      no = spp%n_orbnl
      nk = spp%n_pjnl
      if ( no == 0 ) cycle

      ! Create the group
      call ncdf_def_grp(ncdf,trim(spp%label),grp)

      call ncdf_def_dim(grp,'norbs',no)
      if ( nk > 0 ) then
        call ncdf_def_dim(grp,'nkbs',nk)
      end if
      call ncdf_def_dim(grp,'ntb',nt)

      ! Save the orbital global attributes
      dic = ('Element'.kv.trim(spp%symbol))
      dic = dic//('Label'.kv.trim(spp%label))
      dic = dic//('Atomic_number'.kv.spp%z)
      dic = dic//('Valence_charge'.kv.spp%zval)
      dic = dic//('Mass'.kv.spp%mass)
      dic = dic//('Self_energy'.kv.spp%self_energy)
      dic = dic//('Number_of_orbitals'.kv.spp%norbs)
      dic = dic//('L_max_basis'.kv.spp%lmax_basis)
      dic = dic//('Number_of_projectors'.kv.spp%nprojs)
      dic = dic//('L_max_projs'.kv.spp%lmax_projs)
      dic = dic//('ID'.kv.is)
      call ncdf_put_gatt(grp,atts=dic)
      call delete(dic)

      ! Create all orbital variables...
      call ncdf_def_var(grp,'orbnl_l',NF90_INT,(/'norbs'/), &
          compress_lvl=0,chunks=(/no/))
      call ncdf_def_var(grp,'orbnl_n',NF90_INT,(/'norbs'/), &
          compress_lvl=0,chunks=(/no/))
      call ncdf_def_var(grp,'orbnl_z',NF90_INT,(/'norbs'/), &
          compress_lvl=0,chunks=(/no/))
      call ncdf_def_var(grp,'orbnl_ispol',NF90_INT,(/'norbs'/), &
          compress_lvl=0,chunks=(/no/))
      call ncdf_def_var(grp,'orbnl_pop',NF90_DOUBLE,(/'norbs'/), &
          compress_lvl=0,chunks=(/no/))
      call ncdf_def_var(grp,'cutoff',NF90_DOUBLE,(/'norbs'/), &
          compress_lvl=0,chunks=(/no/))
      call ncdf_def_var(grp,'delta',NF90_DOUBLE,(/'norbs'/), &
          compress_lvl=0,chunks=(/no/))

      ! Create all projector variables...
      if ( nk > 0 ) then
        call ncdf_def_var(grp,'proj',NF90_DOUBLE, (/'ntb ','nkbs'/), &
            compress_lvl=cdf_comp_lvl,chunks=(/nt,1/))
        call ncdf_def_var(grp,'pjnl_l',NF90_INT,(/'nkbs'/), &
            compress_lvl=0,chunks=(/nk/))
        call ncdf_def_var(grp,'pjnl_n',NF90_INT,(/'nkbs'/), &
            compress_lvl=0,chunks=(/nk/))
        call ncdf_def_var(grp,'pjnl_ekb',NF90_DOUBLE,(/'nkbs'/), &
            compress_lvl=0,chunks=(/nk/))
        call ncdf_def_var(grp,'kbcutoff',NF90_DOUBLE,(/'nkbs'/), &
            compress_lvl=0,chunks=(/nk/))
        call ncdf_def_var(grp,'kbdelta',NF90_DOUBLE,(/'nkbs'/), &
            compress_lvl=0,chunks=(/nk/))
        if ( spp%lj_projs ) then
          call ncdf_def_var(grp,'pjnl_j',NF90_INT,(/'nkbs'/), &
              compress_lvl=0,chunks=(/nk/))
        end if
      end if
      call delete(v)

      ! Create orbital projector
      call ncdf_def_var(grp,'orb',NF90_DOUBLE,(/'ntb  ','norbs'/), &
          compress_lvl=cdf_comp_lvl,chunks=(/nt,1/))

      if ( spp%z > 0 ) then ! negative are floating orbitals

        ! Local potential
        call save_rad_func('vna', spp%vna)

        ! Local potential charge density
        call save_rad_func('chlocal', spp%chlocal)

        ! Reduced local potential (rV+2*Zval)
        call save_rad_func('reduced_vlocal', spp%reduced_vlocal)

        if ( spp%there_is_core ) then
          ! Core charge, if it exists, the variable will be created
          ! Hence, the old way of designating whether core is present
          ! or not is removed.
          ! I.e. no Core_flag [1|0] will be saved, Core_flag == .true. is 
          ! apt if 'core' variable exists.
          call save_rad_func('core', spp%core)
        end if

      end if

      call delete(dic)

      ! Save all variables to the group

      ! Save orbital
      call ncdf_put_var(grp,'orbnl_l',spp%orbnl_l(1:no))
      call ncdf_put_var(grp,'orbnl_n',spp%orbnl_n(1:no))
      call ncdf_put_var(grp,'orbnl_z',spp%orbnl_z(1:no))

      allocate(iaux(no))
      do i = 1, no
        if ( spp%orbnl_ispol(i) ) then
          iaux(i) = 1
        else
          iaux(i) = 0 
        end if
      end do
      call ncdf_put_var(grp,'orbnl_ispol',iaux(1:no))
      deallocate(iaux)
      call ncdf_put_var(grp,'orbnl_pop',spp%orbnl_pop(1:no))

      ! Save projector
      if ( nk > 0 ) then
        call ncdf_put_var(grp,'pjnl_l',spp%pjnl_l(1:nk))
        if ( spp%lj_projs ) then
          call ncdf_put_var(grp,'pjnl_j',spp%pjnl_j(1:nk))
        end if
        call ncdf_put_var(grp,'pjnl_n',spp%pjnl_n(1:nk))
        call ncdf_put_var(grp,'pjnl_ekb',spp%pjnl_ekb(1:nk))
      end if

      do i = 1, nk
        p => spp%pjnl(i)
        call ncdf_put_var(grp,'proj',p%f(1:nt),start=(/1,i/))
        call ncdf_put_var(grp,'kbcutoff',p%cutoff,start=(/i/))
        call ncdf_put_var(grp,'kbdelta',p%delta,start=(/i/))
      end do

      do i = 1, no
        p => spp%orbnl(i)
        call ncdf_put_var(grp,'orb',p%f(1:nt),start=(/1,i/))
        call ncdf_put_var(grp,'cutoff',p%cutoff,start=(/i/))
        call ncdf_put_var(grp,'delta',p%delta,start=(/i/))
      end do

    end do

    call ncdf_close(nf)

    call timer('CDF-basis',2)

  contains

    subroutine save_rad_func(name, rfunc)
      use radial, only: rad_func

      character(len=*), intent(in) :: name
      type(rad_func), intent(in) :: rfunc
      type(dictionary_t) :: dic

      ! Only create it if it exists in the pseudo
      if ( rfunc%n <= 0 ) return

      dic = ('cutoff'.kv.rfunc%cutoff) // ('delta'.kv.rfunc%delta)
      call ncdf_def_var(grp, name, NF90_DOUBLE, (/'ntb'/), atts=dic, &
          compress_lvl=cdf_comp_lvl,chunks=(/nt/))
      ! Save immediately
      call ncdf_put_var(grp, name, rfunc%f(1:nt))
      call delete(dic)

    end subroutine save_rad_func

  end subroutine cdf_save_basis

  subroutine cdf_save_basis_new(fname)

    use siesta_geom, only: na_u, isa

    use atmparams, only : nt => NTBMAX
    use atm_types, only : species_info, species, nspecies
    use radial, only : rad_func
    implicit none

    character(len=*), intent(in) :: fname
    ! Local variables
    type(species_info), pointer :: spp
    type(rad_func), pointer :: p
    integer :: cmode, ncid, grp, grp2, is, no, nk, i
    integer :: na_u_, no_, nk_, nt_
    integer :: var
    integer, pointer :: iaux(:)
    real*8, pointer :: raux1(:), raux2(:), raux2D(:,:)

    ! Unluckily the new basis saves only
    ! saved on the IO node. 
    ! No other node must therefore access this routine
    if ( Node /= 0 ) return

    call timer('CDF-basis',1)

    cmode = IOR(NF90_WRITE,NF90_NETCDF4)
    call cdf_err( nf90_open( fname, cmode, ncid ) )

    ! create the BASIS group
  	call cdf_err( nf90_def_grp( ncid, 'BASIS', grp ) )

    ! Create a list of the species associated with each atom
    call cdf_err( nf90_inq_dimid( ncid, 'na_u', na_u_ ) )
    call cdf_err( nf90_def_var( grp, 'basis', &
          NF90_INT, (/ na_u_ /), var ) )
    call cdf_err( nf90_put_att( grp, var, &
          'info', 'Basis of each atom by ID' ) )

    do is = 1 , nspecies
      ! Get current specie
      spp => species(is)

      ! Get array sizes
      no = spp%n_orbnl
      nk = spp%n_pjnl
      if ( no == 0 ) cycle

      ! Create the group
      call cdf_err( nf90_def_grp( grp, trim(spp%label), grp2 ) )

      call cdf_err( nf90_def_dim( grp2, 'norbs', no, no_ ) )
      if ( nk > 0 ) then
        call cdf_err( nf90_def_dim( grp2, 'nkbs', nk, nk_ ) )
      endif
      call cdf_err( nf90_def_dim( grp2, 'ntb', nt, nt_ ) )

      ! Save the orbital global attributes
  	  call cdf_err( nf90_put_att( grp2, NF90_GLOBAL, 'ID', is ) )
  	  call cdf_err( nf90_put_att( grp2, NF90_GLOBAL, &
          'Atomic_number', spp%z ) )
  	  call cdf_err( nf90_put_att( grp2, NF90_GLOBAL, &
          'Mass', spp%mass ) )
  	  call cdf_err( nf90_put_att( grp2, NF90_GLOBAL, &
          'Element', trim(spp%symbol) ) )
  	  call cdf_err( nf90_put_att( grp2, NF90_GLOBAL, &
          'Self_energy', spp%self_energy ) )
  	  call cdf_err( nf90_put_att( grp2, NF90_GLOBAL, &
          'Number_of_projectors', spp%nprojs ) )
  	  call cdf_err( nf90_put_att( grp2, NF90_GLOBAL, &
          'Valence_charge', spp%zval ) )
  	  call cdf_err( nf90_put_att( grp2, NF90_GLOBAL, &
          'Number_of_orbitals', spp%norbs ) )
  	  call cdf_err( nf90_put_att( grp2, NF90_GLOBAL, &
          'L_max_basis', spp%lmax_basis ) )
  	  call cdf_err( nf90_put_att( grp2, NF90_GLOBAL, &
          'L_max_projs', spp%lmax_projs ) )
  	  call cdf_err( nf90_put_att( grp2, NF90_GLOBAL, &
          'Label', trim(spp%label) ) )

      ! Create all orbital variables...
      call cdf_err( nf90_def_var( grp2, 'orbnl_l', &
            NF90_INT, (/ no_ /), var ) )
      call cdf_err( nf90_def_var( grp2, 'orbnl_n', &
            NF90_INT, (/ no_ /), var ) )
      call cdf_err( nf90_def_var( grp2, 'orbnl_z', &
            NF90_INT, (/ no_ /), var ) )
      call cdf_err( nf90_def_var( grp2, 'orbnl_ispol', &
            NF90_INT, (/ no_ /), var ) )
      call cdf_err( nf90_def_var( grp2, 'orbnl_pop', &
            NF90_DOUBLE, (/ no_ /), var ) )
      call cdf_err( nf90_def_var( grp2, 'cutoff', &
            NF90_DOUBLE, (/ no_ /), var ) )
      call cdf_err( nf90_def_var( grp2, 'delta', &
            NF90_DOUBLE, (/ no_ /), var ) )

      ! Create all projector variables...
      if ( nk > 0 ) then
        call cdf_err( nf90_def_var( grp2, 'proj', &
              NF90_DOUBLE, (/ nt_, nk_ /), var ) )
        call cdf_err( nf90_def_var( grp2, 'pjnl_l', &
              NF90_INT, (/ nk_ /), var ) )
        call cdf_err( nf90_def_var( grp2, 'pjnl_n', &
              NF90_INT, (/ nk_ /), var ) )
        call cdf_err( nf90_def_var( grp2, 'pjnl_ekb', &
              NF90_DOUBLE, (/ nk_ /), var ) )
        call cdf_err( nf90_def_var( grp2, 'kbcutoff', &
              NF90_DOUBLE, (/ nk_ /), var ) )
        call cdf_err( nf90_def_var( grp2, 'kbdelta', &
              NF90_DOUBLE, (/ nk_ /), var ) )
        if ( spp%lj_projs ) then
          call cdf_err( nf90_def_var( grp2, 'pjnl_j', &
                NF90_INT, (/ nk_ /), var ) )
        end if
      end if

      ! Create orbital projector
      call cdf_err( nf90_def_var( grp2, 'orb', &
            NF90_DOUBLE, (/ nt_, no_ /), var ) )

      if ( spp%z > 0 ) then ! negative are floating orbitals
        ! Local potential
        if (spp%vna%n > 0) then
          call cdf_err( nf90_def_var( grp2, 'vna', &
                NF90_DOUBLE, (/ nt_ /), var ) )
          call cdf_err( nf90_put_att( grp2, var, &
                'cutoff', spp%vna%cutoff ) )
          call cdf_err( nf90_put_att( grp2, var, &
                'delta', spp%vna%delta ) )
        endif

        ! Local potential charge density
        if (spp%chlocal%n > 0) then
          call cdf_err( nf90_def_var( grp2, 'chlocal', &
                NF90_DOUBLE, (/ nt_ /), var ) )
          call cdf_err( nf90_put_att( grp2, var, &
                'cutoff', spp%chlocal%cutoff ) )
          call cdf_err( nf90_put_att( grp2, var, &
                'delta', spp%chlocal%delta ) )
        endif

        ! Reduced local potential (rV+2*Zval)
        if (spp%reduced_vlocal%n > 0) then
          call cdf_err( nf90_def_var( grp2, 'reduced_vlocal', &
                NF90_DOUBLE, (/ nt_ /), var ) )
          call cdf_err( nf90_put_att( grp2, var, &
                'cutoff', spp%reduced_vlocal%cutoff ) )
          call cdf_err( nf90_put_att( grp2, var, &
                'delta', spp%reduced_vlocal%delta ) )
        endif

        if ( spp%there_is_core ) then
          ! Core charge, if it exists, the variable will be created
          ! Hence, the old way of designating whether core is present
          ! or not is removed.
          ! I.e. no Core_flag [1|0] will be saved, Core_flag == .true. is 
          ! apt if 'core' variable exists.
          if (spp%core%n > 0) then
            call cdf_err( nf90_def_var( grp2, 'core', &
                  NF90_DOUBLE, (/ nt_ /), var ) )
            call cdf_err( nf90_put_att( grp2, var, &
                  'cutoff', spp%core%cutoff ) )
            call cdf_err( nf90_put_att( grp2, var, &
                  'delta', spp%core%delta ) )
          endif
        end if
      end if
    enddo

    call cdf_err( nf90_enddef( ncid ) )

    ! Save all variables to the group
    call cdf_err( nf90_inq_varid( grp, 'basis', var ) )
    call cdf_err( nf90_put_var( grp, var, isa(1:na_u) ) )

    do is = 1 , nspecies
      ! Get current specie
      spp => species(is)

      ! Get array sizes
      no = spp%n_orbnl
      nk = spp%n_pjnl

      if ( no == 0 ) cycle

      ! Inquire the group ID
      call cdf_err( nf90_inq_ncid( grp, trim(spp%label), grp2 ) )

      ! Save orbital
      call cdf_err( nf90_inq_varid( grp2, 'orbnl_l', var ) )
      call cdf_err( nf90_put_var( grp2, var, spp%orbnl_l(1:no) ) )
      call cdf_err( nf90_inq_varid( grp2, 'orbnl_n', var ) )
      call cdf_err( nf90_put_var( grp2, var, spp%orbnl_n(1:no) ) )
      call cdf_err( nf90_inq_varid( grp2, 'orbnl_z', var ) )
      call cdf_err( nf90_put_var( grp2, var, spp%orbnl_z(1:no)) )


      nullify(iaux)
      call re_alloc( iaux, 1, no, 'iaux', 'cdf_basis' )
      do i = 1, no
        if ( spp%orbnl_ispol(i) ) then
          iaux(i) = 1
        else
          iaux(i) = 0 
        end if
      end do
      call cdf_err( nf90_inq_varid( grp2, 'orbnl_ispol', var ) )
      call cdf_err( nf90_put_var( grp2, var, iaux ) )
      call de_alloc( iaux, 'iaux', 'cdf_basis' )

      call cdf_err( nf90_inq_varid( grp2, 'orbnl_pop', var ) )
      call cdf_err( nf90_put_var( grp2, var, spp%orbnl_pop(1:no) ) )

      ! Save projector
      if ( nk > 0 ) then
        call cdf_err( nf90_inq_varid( grp2, 'pjnl_l', var ) )
        call cdf_err( nf90_put_var( grp2, var, spp%pjnl_l(1:nk) ) )
        if ( spp%lj_projs ) then
          call cdf_err( nf90_inq_varid( grp2, 'pjnl_j', var ) )
          call cdf_err( nf90_put_var( grp2, var, spp%pjnl_j(1:nk) ) )
        endif
        call cdf_err( nf90_inq_varid( grp2, 'pjnl_n', var ) )
        call cdf_err( nf90_put_var( grp2, var, spp%pjnl_n(1:nk) ) )
        call cdf_err( nf90_inq_varid( grp2, 'pjnl_ekb', var ) )
        call cdf_err( nf90_put_var( grp2, var, spp%pjnl_ekb(1:nk) ) )

        nullify(raux1,raux2,raux2D)
      	call re_alloc( raux1, 1, nk, 'raux1', 'cdf_basis' )
      	call re_alloc( raux2, 1, nk, 'raux2', 'cdf_basis' )
      	call re_alloc( raux2D, 1, nt, 1, nk, 'raux2D', 'cdf_basis' )
        do i = 1, nk
          p => spp%pjnl(i)
          raux1(i) = p%cutoff
          raux2(i) = p%delta
          raux2D(:,i) = p%f
        enddo
        call cdf_err( nf90_inq_varid( grp2, 'proj', var ) )
        call cdf_err( nf90_put_var( grp2, var, raux2D ) )
        call cdf_err( nf90_inq_varid( grp2, 'kbcutoff', var ) )
        call cdf_err( nf90_put_var( grp2, var, raux1 ) )
        call cdf_err( nf90_inq_varid( grp2, 'kbdelta', var ) )
        call cdf_err( nf90_put_var( grp2, var, raux2 ) )
      	call de_alloc( raux1, 'raux1', 'cdf_basis' )
      	call de_alloc( raux2, 'raux2', 'cdf_basis' )
      	call de_alloc( raux2D, 'raux2D', 'cdf_basis' )
      endif
      nullify(raux1,raux2,raux2D)
      call re_alloc( raux1, 1, no, 'raux1', 'cdf_basis' )
      call re_alloc( raux2, 1, no, 'raux2', 'cdf_basis' )
      call re_alloc( raux2D, 1, nt, 1, no, 'raux2D', 'cdf_basis' )
      do i = 1, no
        p => spp%orbnl(i)
        raux1(i) = p%cutoff
        raux2(i) = p%delta
        raux2D(:,i) = p%f
      end do
      call cdf_err( nf90_inq_varid( grp2, 'orb', var ) )
      call cdf_err( nf90_put_var( grp2, var, raux2D ) )
      call cdf_err( nf90_inq_varid( grp2, 'cutoff', var ) )
      call cdf_err( nf90_put_var( grp2, var, raux1 ) )
      call cdf_err( nf90_inq_varid( grp2, 'delta', var ) )
      call cdf_err( nf90_put_var( grp2, var, raux2 ) )
      call de_alloc( raux1, 'raux1', 'cdf_basis' )
      call de_alloc( raux2, 'raux2', 'cdf_basis' )
      call de_alloc( raux2D, 'raux2D', 'cdf_basis' )
      if ( spp%z > 0 ) then ! negative are floating orbitals
        ! Local potential
        if (spp%vna%n > 0) then
          call cdf_err( nf90_inq_varid( grp2, 'vna', var ) )
          call cdf_err( nf90_put_var( grp2, var, spp%vna%f ) )
        endif
        ! Local potential charge density
        if (spp%chlocal%n > 0) then
          call cdf_err( nf90_inq_varid( grp2, 'chlocal', var ) )
          call cdf_err( nf90_put_var( grp2, var, spp%chlocal%f ) )
        endif

        ! Reduced local potential (rV+2*Zval)
        if (spp%reduced_vlocal%n > 0) then
          call cdf_err( nf90_inq_varid( grp2, 'reduced_vlocal', var ) )
          call cdf_err( nf90_put_var( grp2, var, spp%reduced_vlocal%f ) )
        endif

        if ( spp%there_is_core ) then
          ! Core charge, if it exists, the variable will be created
          ! Hence, the old way of designating whether core is present
          ! or not is removed.
          ! I.e. no Core_flag [1|0] will be saved, Core_flag == .true. is 
          ! apt if 'core' variable exists.
          if (spp%core%n > 0) then
            call cdf_err( nf90_inq_varid( grp2, 'core', var ) )
            call cdf_err( nf90_put_var( grp2, var, spp%core%f ) )
          endif
        end if
      endif
    enddo
    ! Close the file
    call cdf_err( nf90_close( ncid ) )

    call timer('CDF-basis',2)
  end subroutine cdf_save_basis_new

  subroutine cdf_save_fc(fname,istep)
    use siesta_geom, only: na_u
    use m_forces, only: fa

    character(len=*), intent(in) :: fname
    integer, intent(in) :: istep

    type(hNCDF) :: ncdf
    integer :: ia, dir, pm

    if ( Node /= 0 ) return

    ! open the file...
    call ncdf_open(ncdf,fname,mode=ior(NF90_WRITE,NF90_NETCDF4),group='FC')

    if ( istep == 0 ) then
      call ncdf_put_var(ncdf,'fa0',fa)

    else
      ! Atomic index (with respect to ia1)
      ia = mod(istep - 1, 6) + 1
      dir = (ia - 1) / 2 + 1
      pm = mod(ia-1, 2) + 1
      ia = (istep - 1) / 6 + 1

      call ncdf_put_var(ncdf,'fa',fa,count=(/3,na_u/),start=(/1,1,pm,dir,ia/))
    end if

    call ncdf_close(ncdf)

  end subroutine cdf_save_fc
  
#else
  subroutine dummy()
    !pass
  end subroutine dummy
#endif

end module m_ncdf_siesta
