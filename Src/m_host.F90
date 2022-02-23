module m_host
implicit none

contains

  function getNProcessors( ) result(nproc)
  use parallel, only : Node
# ifdef MPI
  use mpi
# endif
  implicit none
  integer :: nproc
  integer   :: iu, io, ierr
  character*128 :: line
  if (Node==0) then
    nproc = 0
    call io_assign( iu )
  	open(unit=iu, file='/proc/cpuinfo')
    do while( .TRUE. )
      read(iu,'(A)',IOSTAT=io) line
      if (io/=0)  exit
      if (index(line,"processor")==1) nproc = nproc+1
    enddo
  endif
# ifdef MPI
  call MPI_Bcast( nproc, 1, MPI_integer, 0, MPI_Comm_World, ierr )
# endif
  end function getNProcessors

  function getAvailMem( ) result(mem)
  use parallel, only : Node
# ifdef MPI
  use mpi
# endif
  implicit none
  integer*8 :: mem
  integer   :: iu, io, ierr
  character*128 :: line
  character*32 :: tag, aux
  if (node==0) then
    call io_assign( iu )
  	open(unit=iu, file='/proc/meminfo')
    do while( .TRUE. )
      read(iu,'(A)',IOSTAT=io) line
      if (io/=0)  exit
      read(line,*) tag, mem, aux
      if (tag=="MemFree:") exit
    enddo
    mem = mem*1024
  endif
# ifdef MPI
  call MPI_Bcast( mem, 1, MPI_integer8, 0, MPI_Comm_World, ierr )
# endif
  end function getAvailMem
  
  subroutine getCHost( nhost, lhost )
  use parallel, only : Node, Nodes
  use alloc
# ifdef MPI
  use mpi
# endif
  implicit none
  integer :: nhost
  integer, pointer :: lhost(:)
  character*(MPI_MAX_PROCESSOR_NAME) :: name, jid
  character*64 :: nw
  character*(MPI_MAX_PROCESSOR_NAME), pointer :: names(:)
  integer :: len, ierr, i, n, nextc, j, nwriters
  logical, pointer :: check(:)
  integer, pointer :: ln(:)

# ifdef MPI
  nullify(lhost,names,check,ln)
  call re_alloc( names, 1, Nodes, "names", "getNumCores" )
  call re_alloc( check, 1, Nodes, "check", "getNumCores" )
  call re_alloc( ln, 1, Nodes, "ln", "getNumCores" )
  check = .TRUE.

  call MPI_Get_processor_name( name, len, ierr )
# define _DEBUG_
# ifdef _DEBUG_
  call getenv( "SLURM_JOB_ID", jid )
  if (jid=="") call getenv( "LSB_JOBID", jid )
  if (jid=="") then
    write(name,"(A,'.',I0.3)") trim(name), Node/8
  else
    call getenv( "SIESTA_NWRITERS", nw )
    if (trim(nw)/='') then
      read(nw,*) nwriters
      nwriters = min(nodes,max(1,nwriters))
      if (nwriters>nodes/getNProcessors( )) then
        write(name,"(A,'.',I0.3)") trim(name), Node/(Nodes/nwriters)
      endif
    endif
  endif
# endif
  call MPI_AllGather( name, MPI_MAX_PROCESSOR_NAME, MPI_CHARACTER, &
      names, MPI_MAX_PROCESSOR_NAME, MPI_CHARACTER, MPI_Comm_World, &
      ierr )

  nhost = 0
  nextc  = 1
  do while(nextc<=Nodes)
    i  = nextc
    check(i) = .FALSE.
    nhost = nhost+1
    ln(nhost) = i-1
    nextc = nextc+1
    do while (nextc<=Nodes)
      if (check(nextc)) then
        if (trim(names(nextc))/=trim(names(i))) exit
        check(nextc) = .FALSE.
      endif
      nextc = nextc+1
    enddo
    j = nextc+1
    do while (j<=Nodes)
      if (check(j) .and. trim(names(j))==trim(names(i))) then
        check(j) = .FALSE.
      endif
      j = j+1
    enddo
  enddo
  call de_alloc( check, "check", "getNumCores" )
  call de_alloc( names, "names", "getNumCores" )
  allocate(lhost(nhost))
  lhost = ln(1:nhost)
  call de_alloc( names, "ln", "getNumCores" )
# else
  nhost = 1
  allocate(lhost(nhost))
  lhost(1) = Node
# endif
  end subroutine getCHost


end module m_host
