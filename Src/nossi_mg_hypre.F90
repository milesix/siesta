!------------------------------------------------------------------------------------
!  Copyright  (c) 2010,2011,2012 INRIA    
!
! author O. Coulaud - Olivier.Coulaud@inria.fr 
!
! This software is a computer program whose purpose is to [descr1ibe
! functionalities and technical features of your software].
!
! This software is governed by the CeCILL C license under French law and
! abiding by the rules of distribution of free software.  You can  use, 
! modify and/ or redistribute the software under the terms of the CeCILL
! license as circulated by CEA, CNRS and INRIA at the following URL
! "http://www.cecill.info". 
!
! As a counterpart to the access to the source code and  rights to copy,
! modify and redistribute granted by the license, users are provided only
! with a limited warranty  and the software's author,  the holder of the
! economic rights,  and the successive licensors  have only  limited
! liability. 
! 
! In this respect, the user's attention is drawn to the risks associated
! with loading,  using,  modifying and/or developing or reproducing the
! software by the user in light of its specific status of free software,
! that may mean  that it is complicated to manipulate,  and  that  also
! therefore means  that it is reserved for developers  and  experienced
! professionals having in-depth computer knowledge. Users are therefore
! encouraged to load and test the software's suitability as regards their
! requirements in conditions enabling the security of their systems and/or 
! data to be ensured and,  more generally, to use and operate it in the 
! same conditions as regards security. 
!------------------------------------------------------------------------------------
!
! The fact that you are presently reading this means that you have had
! knowledge of the CeCILL license and that you accept its terms.                               
!

module nossi_mg_hypre
  !
  ! Solve Poisson equation with Dirichlet boundary conditions in [0,L1]x[0,L2]x[0,L3]
  !   - Delta V = 4 pi rho  in  Omega 
  !      u = g on the boundary
  ! g is approximate by a multipole expansion 
  ! The system is discretize by finite difference method with a second order scheme. 
  ! The matrix system is solved with the HYPRE library
  ! For more details, we invite you to read the HYPRE's manuel 
  ! https://computation.llnl.gov/casc/linear_solvers/sls_hypre.html
  ! The multipole expansions to compute 
  !
  !
  ! We consider the classical second order discretisation of the Laplacian
  ! stencil  :   Discretisation stencil in 3d
  !                                 X   X
  !                                 | /
  !                             X - X - X
  !                               / |
  !                             X   X
  !
  !  Solve  Ax = b 
  !    As A is symmetric we consider the symetric storage of the stencil. 
  !      Then we only have four terms
  !            X                   
  !            |                   
  !            X - X    or      compact scheme with first neigbors 
  !           /                 
  !          X                  
  !          x                  
  !=======================================================================================
  !                                      FDF parameters
  !
  !  MG.UseMG      false use FFT solver for Poisson, true use non full periodic 
  !                  solver based on hypre library.
  !  MG.Solver     0 PFMG multigrid solver (default)
  !	  	   1 PCG  Precontioned conjugate gradient.  The preconditioner is PFMG
  !
  !  MG.BCY        0 Dirichlet boundary conditions in Y direction ; 1 periodic condition in Y
  !  MG.BCY        1 Dirichlet boundary conditions in Z direction ; 1 periodic condition in Z
  !
  !  MG.Order      discretisation order of Posisson operator by compact schemes 
  !                  order 2 or order 4
  !  MG.NumberMultipole l order of the spherical multipole expansion (default 4)
  !
  !  MG.Verbose     1 print the L2 residual norm after the Poisson solver, 0 none (default 0)
  !  MG.Nlevel      l (integer) set MG.N{X,Y,Z} = a *2^l ! best values are l = 4,3 
  !                         default value is 2
  !  MG.Tolerance   Tolerance to decide if the multigrid has converged (default 0.001)
  !  MG.MaxIter     Maximun iterations number for the solver to reach the convergence
  !  MG.PrecondIter Number of iteration for the preconditionner.
  !
  !
  !  MG.Distrib    0 (default) classical UNIFORM distribution of the grid used in parallel Siesta
  !                1  3D distribution UNIFORMSPLITX splt the grid in the three directions according ti bisection algorithm
  !                    NODES SHOULD BE a power of 2  8, 64, ...  -EXPERIMENTAL 
  !  MG.SetMGPoints  if true we assign NX,NY, NZ with the keywords MG.N{X,Y,Z}
  !
  !=======================================================================================
  !   subroutine hypre_solver(rho,energ,V,stress)
  !  input
  !      rho    density function (1D array)
  !  output
  !      energ  local contribution of the electrostatic energy
  !      V      Electrostatic potential
  !      stress stress from electrostatic energy (NOT COMPUTED YET)
  !
  use ISO_C_BINDING, only: C_INT
  use precision
  use nossi_harmonic, only : init_harmonic, evalylm
  implicit none
#ifdef MPI
  include 'mpif.h'
#endif
  ! Public variables
  ! useMG    : if true we consider multigrid method to solve the Poisson problem.
  !
  !----------------------------------------------------------
  ! Private variables
  ! solver_id     : 0 PFMG solver , 1 PCG 
  ! precond_id    : preconditionning method only if solver_id =1
  !                   0 none, 1 PFMG, 2 diagonal
  !
  ! dim                : Space dimension
  ! NGPoints(3)        : Global number of points in each direction
  ! NLPoints(3)        : Local number of points in each direction of the unknown 
  ! NLSiestaPoints(3)  : Local number of points in each direction of the siesta Grid
  ! hypreIndex(3,2)    : The Hypre Grid index space owned by the processor
  ! siestaIndex(3,2)   : The Siesta Grid index space owned by the processor
  !
  ! BoundaryConditions array to specify the boundary condition type 
  !                      in each direction  
  !        BoundaryConditions(i)  0  dirichlet conditions in the i-direction
  !                               1  periodic conditions  in the i-direction
  ! ------------------------------------------------------------------------------------------
  ! pos_x               Coordinate of the points in the x direction
  ! pos_y               Coordinate of the points in the y direction
  ! pos_z               Coordinate of the points in the z direction
  ! h(3)                Space discretisation step in each direction
  ! centre(3)           Centre of the simulation box
  ! ------------------------------------------------------------------------------------------
  ! P_multipole         Degre of the expansion 2 in cartesian or 2,3,4 in spherical
  ! num_coeffMultipole  number of multipole ciefficients (or harmonics)
  ! ylm (:)             value of the harmonics at the current point 
  ! multipoles          multipole coefficients - array of size num_coeffMultipole
  ! ------------------------------------------------------------------------------------------
  !  Hypre Parameter
  !
  ! Matrix
  !    nentries :  Number of values in the stencil (finites differences) here 4 rather than 9
  !    grid     :  Grid object for Hypre
  !    stencil  :  Discretisation stencil in 3d
  !    A        : matrix
  !    b        : right en side
  !    x        : the solution  of the Poisson's system
  ! solver
  !   num_pre num_post : number of pre and post relaxations
  !  relax         : relaxation type
  !                        0 - Jacobi
  !                        1 - Weighted Jacobi (default)
  !                        2 - R/B Gauss-Seidel
  !                        3 - R/B Gauss-Seidel (nonsymmetric)
  ! tolerance      : Parameter to control MG convergence
  ! ------------------------------------------------------------------------------------------
  !
  public ::  nmg_hypre, nossi_MGreadData, hypre_init,  hypre_createMatrix
  public ::  hypre_checkBox, hypre_cleanGridAndMatrix, hypre_solver
  !
  
  logical,public   :: useMG, setMGPoints
  private
  !
  integer, parameter, private  :: dim = 3
  integer, private    :: poissonDist = 0
  integer, private    :: nentries, num_ghost
  logical, private    :: firstTime, firstMultip,firstSolver
  logical, private    :: sphericalMultipole
  !
  integer,private      :: NGPoints(3),NLPoints(3),NLSiestaPoints(3),MGPoints(3)
  integer,private      :: myIndex(3,2),siestaIndex(3,2),hypreIndex(3,2),shiftIndex(3)
  logical,private      :: bord(2,3)
  integer,private      :: BoundaryConditions(3),nbBord(3)

  real(8),private      :: h(3),centre(3)
  real(8), private, allocatable :: pos_x(:),pos_y(:),pos_z(:)
  integer, private     :: P_multipole, num_coeffMultipole
!  real(dp), private,allocatable :: multipoles(:),ylm(:),Anm(:),values(:)
  real(dp), private,allocatable :: multipoles(:),ylm(:),values(:,:)
  integer(8), private  :: grid
  integer(8), private  :: stencil
  integer(8), private  :: A, rhs, sol,B ,Brhs
  real(dp), pointer     :: rho(:)

  !
  integer, private     :: MGOrder
  integer , private    :: stencil_indices(14)  
  integer, private     :: offsets(3,14),nbSiestaValues,nbHypreValues,A_num_ghost(6)

  integer(8), private  :: solver,precond
  integer   , private  :: num_pre,num_post,rap,relax,skip,maxIterations,mg_level,mgVerbose,precondIter
  real(8),    private  :: tolerance
  integer              :: solver_id, precond_id
  real(8), private     :: lx,ly,lz
  !
  integer, private     :: myDistr
  !
  !  Density and potential so the Multigrid Grid and not on the SIESTA one
  !
  real(8), target, allocatable :: rhoMG(:),VMG(:)

  integer, private     :: nf ! To remove after degug
#ifdef DEBUG_MG
!  integer, private     :: nf
  logical, private              :: not_assign=.true.
#endif
  real(dp), private,allocatable :: hypreBuffer(:)
  real(dp), private,allocatable :: BoundaryX0(:,:),BoundaryY0(:,:),BoundaryZ0(:,:)
  integer(C_INT), private       :: hypre_mpi_comm
#ifdef TEST_MG
  real :: start, finish   
#endif
! 
contains
  !
  !------------------------------------------------------------
  ! Read the parameter of the multigrid module
  !
  subroutine nossi_MGreadData()
    !
    use fdf,            only : fdf_get
    use parallel,       only : IOnode,node
    use sys,            only : die
    use m_io,           only : io_assign
    use moreMeshSubs,   only : UNIFORM,UNIFORMSPLITX
    implicit none
    !
#ifdef DEBUG_MG
    character(3)  :: cnode    ! For debug purpose
    character(20) :: logName  ! For debug purpose
#endif

    !
    useMG = fdf_get('MG.UseMG', .false.)
    if( .NOT. useMG ) return 
    ! Read the parameters of the multigrid method if useMG is true
    !
#ifdef DEBUG_MG
    if( not_assign) then
       call io_assign( nf )
       logName = 'debugMG'
#ifdef MPI
       write(cnode,'(i3.3)') node
       logName = trim(logName)//".node"//cnode
#endif
       open( unit=nf, file=logName, form='formatted', status='unknown' )
       not_assign =.false.
    end if
#endif
    
    setMGPoints  = fdf_get('MG.SetMGPoints', .false.)
    MGPoints(1)  = fdf_get('MG.NX', 64 ) 
    MGPoints(2)  = fdf_get('MG.NY', 64 ) 
    MGPoints(3)  = fdf_get('MG.NZ', 64 ) 
    !
    !
    BoundaryConditions(1) = 0 
    BoundaryConditions(2) = fdf_get('MG.BCY', 0 ) 
    BoundaryConditions(3) = fdf_get('MG.BCZ', 0 ) 
    !
    ! Mesh distribution for hypre
    !
    poissonDist  =  fdf_get('MG.Distrib', 0 )
    !
#ifndef MPI
    poissonDist = UNIFORM
#else
    if ( poissonDist /= 0 .AND. poissonDist /= 1 ) then 
       call die( 'redata: Invalid value for MG.Distrib. The value must be O (UNIFORM) or 1 (UNIFORMSPLITX).' )
    end if
    if( poissonDist == 0 ) then
       poissonDist = UNIFORM
    else
       poissonDist = UNIFORMSPLITX
    endif
#endif
    !
    ! Finite difference scheme order
    !       
    MGOrder =  fdf_get('MG.Order', 4 )
    if ( MGOrder /= 2 .AND. MGOrder /= 4 ) then 
       call die( 'redata: Invalid value for MG.Order. The value must be 2 or 4.  ' )
    end if
    firstTime   = .true.
    firstMultip = .true.
    firstSolver = .true.
    ! Minimum number of level in each direction N_i = a_i 2**mg_level.
    ! rebember that a_i is a multible of nsm
    solver_id           = fdf_get('MG.Solver', 0 )
    if( (solver_id < 0) .OR. (solver_id > 1) ) then
       call die( 'redata: Invalid value for MG.Solver. The value must be 0 (PFMG), 1 (PCG).  ' )
    end if
    precond_id          = fdf_get('MG.Precond', 1 ) 
    !
    mg_level            = fdf_get('MG.Nlevel', 2 ) 
    mgVerbose           = fdf_get('MG.Verbose', 0 ) 
    ! Tolerance to decide if the multigrid has converged
    tolerance          = fdf_get( 'MG.Tolerance', 0.001_dp )
    ! Maximum Number of iterations for the solver
    maxIterations       = fdf_get('MG.MaxIter', 40 )
    ! Maximum Number of iterations for the preconditioner
    precondIter         = fdf_get('MG.PrecondIter', 2 )
    !
    !========================================================================================
    ! Advanced parameters
    !
    ! Number of prerelaxation steps:
    num_pre            = fdf_get( 'MG.PreRelaxSteps', 2 )
    ! Number of postrelaxation steps:
    num_post           = fdf_get( 'MG.PostRelaxSteps', 2 )
    !========================================================================================
    !
#ifdef NOSSI_MG
    sphericalMultipole = fdf_get( 'MG.SphericalMultipole', .true. )
#else
    sphericalMultipole = .true.
#endif
    ! relaxation type:
    if (sphericalMultipole ) then
       P_multipole =  fdf_get(  'MG.NumberMultipole', 4 )
    else
       P_multipole = 2
    end if
    if (ionode) then
       write(6,'(a,L1)') "Multigrid points given: ",setMGPoints
       
       write(6,'(a)') &
            "readata: Use MG for Poisson's equation with dirichlet conditions."
       if( poissonDist == UNIFORM ) then 
          write(6,'(a,i5)')  &
               'redata: MG - MG Distrib                                = UNIFORM'
       else
          write(6,'(a,i5)') &
               'redata: MG - MG Distrib                                = UNIFORMSPLIX'
       end if
       write(6,'(a,i5)') &
            'redata: MG - MG Boundary in Y                          = ',BoundaryConditions(2)
       write(6,'(a,i5)') &
            'redata: MG - MG Boundary in Z                          = ',BoundaryConditions(3)
       write(6,'(a,i5)') &
            'redata: MG - MG Order                                  = ',MGOrder
       write(6,'(a,i5)') &
            'redata: MG - Minimum number of level in each direction = ',mg_level
       write(6,'(a,i5)') &
            'redata: MG - Maximum Number of iterations              = ',maxIterations
       write(6,'(a,f10.4)') &
            'redata: MG - Tolerance for the convergence             = ',tolerance
       write(6,'(a,I5)') &
            'redata: MG - Verbose parameter                         = ',mgVerbose
       write(6,'(a,i5)') &
            'redata: MG - Maximum Number of iterations precondioner = ',precondIter
       write(6,'(a,i5)') &
            'redata: MG - Number of prerelaxation steps             = ',num_pre
       write(6,'(a,i5)') &
            'redata: MG - Number of postrelaxation steps            = ',num_post
       write(6,'(a,4x,L1)') &
            'redata: MG Spherical Multipole                         = ',sphericalMultipole
       if(sphericalMultipole) then 
          write(6,'(a,3x,I2)') &
               'redata: MG Spherical Multipole order                   = ',P_multipole
       else
          write(6,'(a,3x,I2)') &
               'redata: MG Cartesian Multipole order                   = ',P_multipole
       endif
    endif
    !
 end subroutine nossi_MGreadData
 !
 ! ----------------------------------------------------------
 ! This method compute the number of points n in dim direction 
 !   according to:
 !     n = a * 2^mg_level           if mod(nsm,2) == 0
 !       = a * nsm * 2^mg_level      otherwise
 !
 subroutine nmg_hypre( dim, n,nsm)
   ! dim   x, y or z directions (1,2,3)
   ! n     number of points given by siesta in dim direction
   ! nsm   number of subdivision
   !
    use parallel, only : ProcessorY,Nodes,Node,IOnode
    use sys,      only : die
    !
    implicit none
    integer, intent(in)    :: dim,nsm
    integer, intent(inout) :: n
    !
    integer :: p , r, q,n1,am,l
    
#ifdef  DEBUG_MG
    integer oldn
  !
    oldn = n
#endif
    am = 13 ;
    !
    if(.NOT. setMGPoints) then 
       ! mg_level mustn't be 0
       if(mod(nsm,2) == 0 ) then 
          p = 2**mg_level
       else
          p = nsm*(2**mg_level) ; am = am/nsm
       endif
       q = n/p ; r = n -q*p
       n1 = n 
       if(r/=0) n1 = n1 + p - r
       n1 = n1/p
       l  = 0
       do while( n1 > am )
          if(mod(n1,2) == 0 ) then 
             l = l + 1 ; n1 = n1/2
          else
             n1 = n1 + 1
          endif
       end do
       n = p*n1*2**l
    else
       if (n <  MGPoints(dim) ) then 
          n = MGPoints(dim)
       else
          if(IOnode) then 
              write(*,*) "$$$$$"
             write(*,*)  "$$$$$"
             write(*,'(a,I2,I5,a,I5)') '$$ MG In direction ',dim, "the number of points ",MGPoints(dim),  &
                  " is < the number of points given by MeshCutoff: ", n
             write(*,*) "$$$$$"
             write(*,*) "$$$$$"
          end if
          call die('Invalid value for MG.N?')
       end if
    end if
    !   
#ifdef  DEBUG_MG
    if(IOnode) then    
       write(nf,'(A,I2,A,I5,A,I5)') 'nmg_hypre --> dim:',dim, &
            '  Old number of points N =', oldn,' New points number      :',n
    end if
#endif
  end subroutine nmg_hypre
  !
  ! ----------------------------------------------------------------------------
  !
  subroutine hypre_init(box, maxOrbitalRadius)
!  subroutine hypre_init(box,myBoxIn, maxOrbitalRadius)
    !
    ! box              Size of the simulation domain 
    ! myBoxIn          Indexes of the Siesta grid used the current processor
    ! maxOrbitalRadius Maximal cutoff radius of orbital functions
    ! 
    use precision,    only: grid_p
    use parallel,     only: node, ProcessorY,Nodes
    use sys,          only: die
    use siesta_geom,  only: ucell,xa
    use mesh,         only: nsm, cmesh, nmeshg,nsp
    use moreMeshSubs, only: meshDistr,setMeshDistr,UNIFORMSPLITX
    use siestaXC,     only: myMeshBox
    !
    implicit none 
    real(8), intent(in) :: box(3,3)
    real(dp),intent(in) :: maxOrbitalRadius
    ! Internal variables
    integer  :: i, ierr, ind,n,m
    integer  :: ilower(3), iupper(3),myBox(2,3)
    integer  :: nml(3), nmpl, ntml(3), ntpl   
    integer, pointer :: myBoxIn(:,:)

#ifndef MPI
    integer  :: MPI_COMM_WORLD=0
#endif
    !
#ifdef DEBUG
      call write_debug( '    PRE nossi_hypre_init' )
#endif
#ifdef DEBUG_MG
  write(nf,*) "Begin nossi_hypre_init"
#endif
    CALL timer ('HYPRE_INIT',1)
    !
    call setMeshDistr( poissonDist, nsm, nsp,  nml, nmpl, ntml, ntpl )
    myBoxIn => meshDistr(poissonDist)%box(1:2,1:3,node+1)
#ifdef MPI 
    if( (Nodes/ProcessorY)*ProcessorY /= Nodes) then
       call die( 'redata: Invalid value for Nodes or ProcessorY. Nodes HAVE TO BE a multiple of ProcessorY.  ' )
    end if
#endif       
    if(sphericalMultipole ) then
       num_coeffMultipole = ( P_multipole+1)**2
    else
       num_coeffMultipole = 10
    endif
    if(.NOT. allocated(multipoles) ) then 
       allocate(multipoles(num_coeffMultipole))
       if(sphericalMultipole ) then 
          allocate (ylm(num_coeffMultipole) ) !,values(num_coeffMultipole), Anm( ( P_multipole+1)* ( P_multipole+2)/2))
       end if
    endif
    call  init_harmonic(num_coeffMultipole,P_multipole)
    !
    centre(1)   = box(1,1)*0.5_dp
    centre(2)   = box(2,2)*0.5_dp
    centre(3)   = box(3,3)*0.5_dp
    !
    lx = box(1,1) ; ly = box(2,2); lz =  box(3,3)
    !
    myBox(:,:)  = myBoxIn(:,:)
    NGpoints(:) = nmeshg(:)
    IF ((box(1,2) /= 0.0_dp).OR.(box(1,3) /= 0.0_dp).OR. &
         (box(2,1) /= 0.0_dp).OR.(box(2,3) /= 0.0_dp).OR. &
         (box(3,1) /= 0.0_dp).OR.(box(3,2) /= 0.0_dp)) THEN
       WRITE(*,"('MG solver:  ******ERROR******')")
       WRITE(*,"('MG solver: The multigrid solver only deals with orthorhombic cells.   ')")
       STOP 'stop hypre_init'
    END IF

#ifdef MG_FULLGRID
    ! We construct the full grid (not only the grid where we have only unknowns in periodic) 
    NGpoints(:) = NGpoints(:) + 1
#endif
!
#ifdef DEBUG_MG
    write(nf,*) 'Npoints: ',NGpoints
    WRITE(NF,*) 'nmeshg   ',nmeshg  
    write(nf,*) 'cmesh:  '
    write(nf,*)' ',cmesh(1:3,1)
    write(nf,*)' ',cmesh(1:3,2)
    write(nf,*)' ',cmesh(1:3,3)
    write(nf,*) 'box:  '
    write(nf,*)' ',box(1:3,1)
    write(nf,*)' ',box(1:3,2)
    write(nf,*)' ',box(1:3,3)
    write(nf,*)'centre: '
    write(nf,*)centre(1:3)
    write(nf,*) 'myBox 1:  ',myBox(1:2,1)
    write(nf,*) 'myBox 2:  ',myBox(1:2,2)
    write(nf,*) 'myBox 3:  ',myBox(1:2,3)
#endif
    !
    siestaIndex(1:3,1)  = 1 + (myBox(1,1:3)-1)*nsm
    siestaIndex(1:3,2)  = 1 + (myBox(2,1:3)-1)*nsm
    siestaIndex(1:3,2)  = siestaIndex(1:3,2) + nsm - 1
#ifdef MG_FULLGRID
    myIndex(:,:)        = siestaIndex(:,:) 
    ! We construct the full grid (not only the grid where we have only unknowns in preriodic)
    if (myIndex(1,2) ==  nmeshg(1) ) myIndex(1,2) = myIndex(1,2) + 1
    if (myIndex(2,2) ==  nmeshg(2) ) myIndex(2,2) = myIndex(2,2) + 1
    if (myIndex(3,2) ==  nmeshg(3) ) myIndex(3,2) = myIndex(3,2) + 1
    !
    hypreIndex(1:3,:)  = siestaIndex(1:3,:) - 1 ! ???? WRONG !!!!
#else
    ! We consider only interior points where there is no periodic boundary contionons
    hypreIndex(1:3,:)  = siestaIndex(1:3,:) - 1  
#endif
    !
    if ( BoundaryConditions(2) /= 0 )  hypreIndex(2,:)     = hypreIndex(2,:) + 1  
    if ( BoundaryConditions(3) /= 0 )  hypreIndex(3,:)     = hypreIndex(3,:) + 1  
    ! 
    bord(:,:)      = .false.
    shiftIndex(:)  = 0 
    nbbord(:)      = 0
    do i = 1, 3 
       if( (siestaIndex(i,1) == 1) .and.  (BoundaryConditions(i) == 0) )   then 
          bord(1,i)        = .true. 
          nbbord(i)        =  nbbord(i) + 1     
          shiftIndex(i)    = 1
          hypreIndex(i,1)  = 1  ! otherwise we have -1
       end if
       if( (siestaIndex(i,2) == NGpoints(i)) .and.  (BoundaryConditions(i) == 0)) then
           nbbord(i) =  nbbord(i) + 1 
           bord(2,i) = .true. 
        end if
        if ( BoundaryConditions(i) > 0 ) then 
           BoundaryConditions(i)  = NGPoints(i) 
        end if
    end do
    NLPoints(:)        = hypreIndex(:,2)  - hypreIndex(:,1)  + 1
    NLSiestaPoints(:)  = siestaIndex(:,2) - siestaIndex(:,1) + 1
    nbSiestaValues     = NLSiestaPoints(1)*NLSiestaPoints(2)*NLSiestaPoints(3)
    nbHypreValues      = NLPoints(1)*NLPoints(2)*NLPoints(3)
    !
    H(1) = cmesh(1,1)/nsm
    H(2) = cmesh(2,2)/nsm
    H(3) = cmesh(3,3)/nsm
    !
#ifdef DEBUG_MG
    write(nf,'(A)')  '   hypreIndex:       Start     End'
    write(nf,'(A,I5,4X,I5)')  'hypreIndex X:      ',hypreIndex(1,1),hypreIndex(1,2) 
    write(nf,'(A,I5,4X,I5)')  'hypreIndex Y:      ',hypreIndex(2,1),hypreIndex(2,2) 
    write(nf,'(A,I5,4X,I5)')  'hypreIndex Z:      ',hypreIndex(3,1),hypreIndex(3,2) 
    write(nf,'(A)')  '   siestaIndex:       Start     End'
    write(nf,'(A,I5,4X,I5)')  'siestaIndex X:      ',siestaIndex(1,1),siestaIndex(1,2) 
    write(nf,'(A,I5,4X,I5)')  'siestaIndex Y:      ',siestaIndex(2,1),siestaIndex(2,2) 
    write(nf,'(A,I5,4X,I5)')  'siestaIndex Z:      ',siestaIndex(3,1),siestaIndex(3,2)

    write(nf,'(A,3I5)')  'Local points:   ',NLPoints
    write(nf,'(A,I8,A,I8)')  ' Siesta points:   ',nbSiestaValues,' MG points:   ',nbHypreValues
    !
    write(nf,'(A)')             '       R  L  shiftIndex nbbord BoundaryConditions '
    write(nf,'(A,2(L,1X),2x,2(I4,5X),2x,I5)')   ' x:   ',bord(:,1),shiftIndex(1),nbbord(1),BoundaryConditions(1)
    write(nf,'(A,2(L,1X),2x,2(I4,5X),2x,I5)')   ' y:   ',bord(:,2),shiftIndex(2),nbbord(2),BoundaryConditions(2)
    write(nf,'(A,2(L,1X),2x,2(I4,5X),2x,I5)')   ' z:   ',bord(:,3),shiftIndex(3),nbbord(3),BoundaryConditions(3)
    write(nf,'(A,3E14.5)') "Mesh Step size (Bohr): ",H(:)
#endif
    !
    call hypre_checkBox(maxOrbitalRadius)
!
    allocate( pos_x(NGPoints(1)+1),pos_y(shiftIndex(2):NGPoints(2)+1),pos_z(shiftIndex(3):NGPoints(3)+1)) 

    Do i = 0, NGPoints(1)
       pos_x(i+1) = H(1) * i 
    end do
    Do i = 0, NGPoints(2) 
       pos_y(i+1) = H(2) * i 
    end do
    Do i = 0, NGPoints(3)
       pos_z(i+1) = H(3) * i 
    end do
    ! 
    if (MGOrder ==  2) then     ! second order Laplacian stencil 
       nentries     = 4
       offsets(:,1) = (/0,0,0/) 
       offsets(:,2) = (/1,0,0/)
       offsets(:,3) = (/0,1,0/)
       offsets(:,4) = (/0,0,1/)
    else     ! fourth order Laplacian stencil 
       if (BoundaryConditions(3) > 0 ) then
          pos_z(0)             = pos_z(NGPoints(3))
          pos_z(NGPoints(3)+1) = pos_z(1)
       end if
       if (BoundaryConditions(2) > 0 ) then
          pos_y(0)             = pos_y(NGPoints(2))
          pos_y(NGPoints(3)+1) = pos_y(1)
       end if

       nentries     = 14
       offsets(:,1) = (/0,0,0/) 
       offsets(:,2) = (/1,0,0/)
       offsets(:,3) = (/0,1,0/)
       offsets(:,4) = (/0,0,1/)
       offsets(:,5) = (/1,0,1/)
       offsets(:,6) = (/1,0,-1/)
       offsets(:,7) = (/0,-1,1/)
       offsets(:,8) = (/1,-1,0/)
       offsets(:,9) = (/0,1,1/)
       offsets(:,10) = (/1,1,0/)
       offsets(:,11) = (/1,-1,1/)
       offsets(:,12) = (/-1,1,1/)
       offsets(:,13) = (/1,1,-1/)
       offsets(:,14) = (/1,1,1/)
    endif
    num_ghost    = 1
    !
    grid    = 0
    stencil = 0
    do i = 1, nentries
       stencil_indices(i) = i-1
    enddo
    !
#ifdef MPI
    call MPI_COMM_DUP(MPI_COMM_WORLD, hypre_mpi_comm,ierr)
#else
    hypre_mpi_comm = MPI_COMM_WORLD
#endif
    !
    ! Set up a grid. Each processor describes the piece
    !     of the grid that it owns. 
    ! Create an empty 3D grid object 
    call HYPRE_StructGridCreate(hypre_mpi_comm, dim, grid,ierr)    
    !
    call  HYPRE_StructGridSetExtents(grid, hypreIndex(1,1), hypreIndex(1,2),ierr) 
    ! 
    call HYPRE_StructGridSetPeriodic (grid,BoundaryConditions(1), ierr)
    !
    ! This is a collective call finalizing the grid assembly.
    !   The grid is now ``ready to be used'' 
    call HYPRE_StructGridAssemble(grid,ierr) 
    !
    ! Create an empty vector object 
    call HYPRE_StructVectorCreate(hypre_mpi_comm, grid, rhs,ierr)
    call HYPRE_StructVectorCreate(hypre_mpi_comm, grid, sol,ierr)
    !
    ! Indicate that the vector coefficients are ready to be set 
    call HYPRE_StructVectorInitialize(rhs,ierr)
    call HYPRE_StructVectorInitialize(sol,ierr)
    !
    allocate(hypreBuffer( nbHypreValues ) )
    if( bord(1,3) .OR.  bord(2,3) ) then  ! z = 0 
       allocate( &
        BoundaryZ0(max(1,siestaIndex(1,1)-1):(siestaIndex(1,2)+1), &
                   max(shiftIndex(2),siestaIndex(2,1)-1):(siestaIndex(2,2)+1)))
    end if
    if( bord(1,2) .OR.  bord(2,2) ) then  ! y = 0 
       allocate( &
        BoundaryY0(max(1,siestaIndex(1,1)-1):(siestaIndex(1,2)+1), &
                  max(shiftIndex(3),siestaIndex(3,1)-1):(siestaIndex(3,2)+1)))
    end if
    if( bord(1,1) .OR.  bord(2,1)) then  ! x = 0 
       allocate( &
      BoundaryX0(max(shiftIndex(2),siestaIndex(2,1)-1):(siestaIndex(2,2)+1), &
                 max(shiftIndex(3),siestaIndex(3,1)-1):(siestaIndex(3,2)+1)))
    end if
    !
    !-----------------------------------------------------------------------
    ! Initialize 4th order Right Hand Side and buffer for Boundary Conditions
    !-----------------------------------------------------------------------
    if (MGOrder ==  4) then 
      call HYPRE_StructVectorCreate(hypre_mpi_comm, grid, Brhs,ierr)
      call HYPRE_StructVectorInitialize(Brhs,ierr)
      ! /!\ 2+.... because NLPoints comes from hypreIndex and Periodic conditions 
    endif
    !-----------------------------------------------------------------
    if( poissonDist == UNIFORMSPLITX) then
       ! We construct the full grid (not only the grid where we have only unknowns in preriodic
       allocate(rhoMG(nbSiestaValues),VMG(nbSiestaValues) )
    endif
    !
    CALL timer ('HYPRE_INIT',2)
#ifdef DEBUG_MG
    write(nf,*) "End nossi_hypre_init"
#endif
#ifdef DEBUG
      call write_debug( '    POS nossi_hypre_init' )
#endif
    !
    end subroutine hypre_init
  !
  ! -------------------------------------------------------------------------
  ! nossi_hypre_createMatrix
  !    Build the matrix of the Poisson'sequation according to
  !           [A_ii 0; 0 I] 
  ! -------------------------------------------------------------------------
  !
    subroutine hypre_createMatrix()
    !
    use alloc,       only : re_alloc, de_alloc
    use parallel,    only : IOnode 
    implicit none
    integer i,j,k,nvalues,ierr,n_max
    integer  bc_ilower(3), bc_iupper(3),ilower(3), iupper(3)
    real(8)              ::invhx2,invhy2,invhz2,coeffMatA(14)
    real(dp)             :: coeffMatB(14)
    !
#ifdef DEBUG
      call write_debug( '    PRE nossi_hypre_createMatrix' )
#endif
    CALL timer ('HYPRE_MATRIX',1)
#ifdef DEBUG_MG
    write(nf,*) "Begin nossi_hypre_createMatrix"
#endif
    !
    ! Define the discretization stencil 
    !------------------------------------------
    !
    ! Create an empty stencil object (We use symetric storage) 
    !
     call HYPRE_StructStencilCreate(dim, nentries , stencil, ierr)
    !
    ! Define the geometry of the stencil. Each represents a
    !        relative offset (in the index space). 
    do i = 1, nentries
       call HYPRE_StructStencilSetElement(stencil, i-1, offsets(1,i),ierr)
    end do

    !
    ! Set up a Struct Matrix 
    !------------------------------------------
    !
    ! Create an empty matrix object 
     call HYPRE_StructMatrixCreate(hypre_mpi_comm, grid, stencil, A,ierr)
!    write(nf,*) 'NOSSI_StructMatrixCreate ierr',ierr
     call HYPRE_StructMatrixSetConstantEn(A, nentries, stencil_indices,ierr)
    !
    ! Use symmetric storage 
     call HYPRE_StructMatrixSetSymmetric(A, 1,ierr)
     do i = 1, 6
        A_num_ghost(i) = num_ghost
     end do
 !   write(nf,*) 'HYPRE_StructMatrixSetNumGhost ierr',ierr, A_num_ghost
     call HYPRE_StructMatrixSetNumGhost(A, A_num_ghost,ierr);
 !   write(nf,*) 'HYPRE_StructMatrixInitialize ierr',ierr

     !
     ! Indicate that the matrix coefficients are ready to be set 
     call HYPRE_StructMatrixInitialize(A,ierr)
  !  write(nf,*) 'HYPRE_StructMatrixInitialize ierr',ierr
     !
     ! Set the stencil values in the interior. Here we set the values
     !	at every node. We will modify the boundary nodes later. 
     !
     invhx2 = 1.0_dp/(H(1)*H(1))
     invhy2 = 1.0_dp/(H(2)*H(2))
     invhz2 = 1.0_dp/(H(3)*H(3))
     !
     !  We have n*m grid points, each with 3 stencil entries 
     !   The order is left-to-right, bottom-to-top 
     ! 
     if(MGOrder == 2) then 
        coeffMatA(1) = 2.0_dp*(invhx2+invhy2+invhz2)
        coeffMatA(2) = -invhx2
        coeffMatA(3) = -invhy2
        coeffMatA(4) = -invhz2
     else if (MGOrder == 4) then 
        coeffMatA(1)  =  8.0_dp/6.0_dp * (invhx2+invhy2+invhz2)
        coeffMatA(2)  = -2.0_dp/6.0_dp  * (invhx2) 
        coeffMatA(3)  = -2.0_dp/6.0_dp  * (invhy2) 
        coeffMatA(4)  = -2.0_dp/6.0_dp  * (invhz2) 
        coeffMatA(5)  = -1.0_dp/12.0_dp  * (invhx2+invhz2) 
        coeffMatA(6)  = -1.0_dp/12.0_dp  * (invhx2+invhz2) 
        coeffMatA(7)  = -1.0_dp/12.0_dp  * (invhy2+invhz2) 
        coeffMatA(8)  = -1.0_dp/12.0_dp  * (invhx2+invhy2)
        coeffMatA(9)  = -1.0_dp/12.0_dp  * (invhy2+invhz2) 
        coeffMatA(10) = -1.0_dp/12.0_dp  * (invhx2+invhy2) 
        coeffMatA(11) =  0.0_dp
        coeffMatA(12) =  0.0_dp
        coeffMatA(13) =  0.0_dp
        coeffMatA(14) =  0.0_dp
    end if
     !
     call HYPRE_StructMatrixSetConstantVa(A, nentries, stencil_indices, coeffMatA(1),ierr)
#ifdef DEBUG_MG
    write(nf,*) '        HYPRE_StructMatrixAssemble   '
#endif
    call HYPRE_StructMatrixAssemble(A,ierr)
!    call HYPRE_StructMatrixPrint(A,0,ierr)
    if(MGOrder == 4) then 
    !---------------------------------------------------------------------------------------------------------------
    ! Set up the rhs Matrix B and update the rhs
    !---------------------------------------------------------------------------------------------------------------
       call HYPRE_StructMatrixCreate(hypre_mpi_comm, grid, stencil, B,ierr)
       call HYPRE_StructMatrixSetConstantEn(B, nentries, stencil_indices,ierr)
       call HYPRE_StructMatrixSetSymmetric(B, 1,ierr)
       call HYPRE_StructMatrixSetNumGhost(B, A_num_ghost,ierr);
       call HYPRE_StructMatrixInitialize(B,ierr)
!       coeffMatB(:)  = 0.0_dp
       coeffMatB(1)  = 6.0_dp/12.0_dp 
       coeffMatB(2)  = 1.0_dp/12.0_dp
       coeffMatB(3)  = 1.0_dp/12.0_dp
       coeffMatB(4)  = 1.0_dp/12.0_dp
       coeffMatB(5)  = 0.0_dp
       coeffMatB(6)  = 0.0_dp
       coeffMatB(7)  = 0.0_dp
       coeffMatB(8)  = 0.0_dp
       coeffMatB(9)  = 0.0_dp
       coeffMatB(10) = 0.0_dp
       coeffMatB(11) = 0.0_dp
       coeffMatB(12) = 0.0_dp
       coeffMatB(13) = 0.0_dp
       coeffMatB(14) = 0.0_dp
       call HYPRE_StructMatrixSetConstantVa(B, nentries, stencil_indices, coeffMatB(1),ierr)
       call HYPRE_StructMatrixAssemble(B,ierr)
    end if
    !---------------------------------------------------------------------------------------------------------------
#ifdef DEBUG_MG
    write(nf,*) "End nossi_hypre_createMatrix"
#endif
#ifdef DEBUG
      call write_debug( '    POS nossi_hypre_createMatrix' )
#endif
    CALL timer ('HYPRE_MATRIX',2)
!
  end subroutine hypre_createMatrix
  !
  ! -------------------------------------------------------------------------
  !   subroutine hypre_solver(rho,energ,V,stress)
  !     The Poisson's solver --  A v = 4 pi rho
  !  input
  !      rho    density function (1D array)  in UNIFORM dist, sequential form  
  !  output
  !      energ  local contribution of the electrostatic energy
  !      V      Electrostatic potential on the grid point that the processor owns
  !      stress stress from electrostatic energy (NOT COMPUTED YET)
  !
  ! -------------------------------------------------------------------------
  !
  subroutine hypre_solver(rhoSiesta,energ,VSiesta,stress)
    !
    use sys,          only : die
    use parallel,     only : Node
    use moreMeshSubs, only : UNIFORM,UNIFORMSPLITX, KEEP
    use moreMeshSubs, only : setMeshDistr, distMeshData
    use mesh,         only : nmeshg,nsm,nsp

    !
    implicit none 
    real(dp), target,    intent(in)  :: rhoSiesta(1)
    real(grid_p),target, intent(out) :: VSiesta(1)
    real(dp),     intent(out)        :: energ
    real(dp),     intent(out)        :: stress(1:3,1:3)
    !
    !  --- Internal definitions ---
    integer               :: ilower(3),iupper(3),ierr,i,j,k,ind
    integer               :: nml(3), nmpl, ntml(3), ntpl
!    real(dp), pointer     :: rho(:)
    real(grid_p), pointer :: V(:)
    real(dp), external    :: ddot
    !
#ifdef DEBUG
      call write_debug( '    PRE hypre_solver' )
#endif
    CALL timer ('HYPRE_SOLVER',1)
#ifdef DEBUG_MG
    write(nf,*) "Begin hypre_solver"
#endif    
!
    if(poissonDist == UNIFORMSPLITX) then
       CALL timer ('HYPRE_DIST',1)
       call setMeshDistr( UNIFORMSPLITX, nsm, nsp, nml, nmpl, ntml, ntpl )
       call distMeshData( UNIFORM, rhoSiesta, UNIFORMSPLITX, rhoMG, KEEP )  !on ne peut pas avec deux rho
       CALL timer ('HYPRE_DIST',2)
       rho => rhoMG
       v   => VMG
    else
       rho => rhoSiesta
       v   => VSiesta
    endif
    V(:) = 0.0_dp    
    call  setRHSandX(rho,V)
    !
    call  hypre_poisson()
    !
!    write(nf,*) 'HYPRE_StructVectorGetBoxValues'
    call setSolutionInPotential(V)
    if(poissonDist == UNIFORMSPLITX) then
       CALL timer ('HYPRE_DIST',1)
       call setMeshDistr( UNIFORM, nsm, nsp, nml, nmpl, ntml, ntpl )
       call distMeshData( UNIFORMSPLITX, V, UNIFORM, VSiesta, KEEP )  
       CALL timer ('HYPRE_DIST',2)
    endif
    !
    !  --- Take into account SIESTA works with 2*v instead of v ---
    !
    call DSCAL(nbSiestaValues,2.0_grid_p,VSiesta,1)
    !
    !  --- Compute the local contribution of the electrostatic energy ---
    !
    energ = ddot(nbSiestaValues,rhoSiesta,1,VSiesta,1)  
    energ = 0.5_dp*energ*h(1)*h(2)*h(3)
   !
   !  --- Compute stress tensor ---
   !   TO DO
   stress(:,:) = 0.0_dp
   !
   nullify(rho,V)
   !
#ifdef DEBUG_MG
    write(nf,*) "End hypre_solver"
#endif   
   CALL timer ('HYPRE_SOLVER',2)
#ifdef DEBUG
      call write_debug( '    POS hypre_solver' ) 
#endif
   !
 end subroutine hypre_solver
 !
 ! -------------------------------------------------------------------------
 !  setSolutionInPotential
 ! -------------------------------------------------------------------------

 subroutine setSolutionInPotential(V)
   !
   real(dp), intent(out)  :: V(*)
   !
   integer :: m,index,i,j,k,ierr,imax,jmax,kmax
   real(dp) :: xx,error,err,skipy
   !
#ifdef DEBUG_MG
    write(nf,*) "Begin setSolutionInPotential"
#endif 
   call HYPRE_StructVectorGetBoxValues(sol, hypreIndex(1,1), hypreIndex(1,2),  hypreBuffer(1), ierr)
   !
   ! Fill inner points of the grid
   !
   m     = 1  
   index =  1 + shiftIndex(3)*NLSiestaPoints(1)*NLSiestaPoints(2)
   skipy = shiftIndex(2)*NLSiestaPoints(1)
   do k = hypreIndex(3,1), hypreIndex(3,2)
      index = index + skipy
      do j =  hypreIndex(2,1), hypreIndex(2,2)
         index = shiftIndex(1) + index
         do i = hypreIndex(1,1), hypreIndex(1,2)
            V(index) = hypreBuffer(m) 
            m = m + 1 ; index = index + 1
         end do
      end do
   end do
   !
#ifdef TEST_MG
   call HYPRE_StructVectorPrint(sol,0,ierr)
   call compareSolution(V)
   stop "Test multigrid ended"
#endif
#ifdef DEBUG_MG
    write(nf,*) "End setSolutionInPotential"
#endif 
   !
 end subroutine setSolutionInPotential
  !
  ! -------------------------------------------------------------------------
  !  Function that approximates the potential V on the boundary of the domain
  ! -------------------------------------------------------------------------
  !
   real function G(x,y,z)
    !  Computes the potential V(x) = integral (rho(y) * 1/|x-y| ) at a given point on the boundary of the box
    !           through a multipole expansion.
    !  See, for example, Cheng, Greegard and Rokhlin JCP 155, 468-498 (1999)
    !
    use nossi_harmonic, only : evalylm

    implicit none
!    real(8)             :: G
    real(8), intent(in) :: x,y,z
    ! ---------------------------------------------------
    real(8) :: x1mc,x2mc,x3mc,invr,contrib,rlp1
    integer :: l,m,ind
    !
#ifdef TEST_MG
     G = solution(x,y,z) 
     return 
#endif
    x1mc = x - centre(1)
    x2mc = y - centre(2)
    x3mc = z - centre(3)
    !
    !  --- Compute electrostatic potential ---
    !
    ! evaluate the spherical harmonic function in point(x,y,z)
    call evalylm(x1mc,x2mc,x3mc,ylm)
    invr = 1.0_dp/sqrt(x1mc*x1mc+x2mc*x2mc+x3mc*x3mc) 
    rlp1 = invr
    ! initialize with l=m=0
    G       = multipoles(1)*ylm(1)!*rlp1 
    contrib = multipoles(2)*ylm(2)+multipoles(3)*ylm(3)+multipoles(4)*ylm(4)
!    x1mc    = multipoles(5)*ylm(5)+multipoles(6)*ylm(6)+multipoles(7)*ylm(7)+multipoles(8)*ylm(8)+multipoles(9)*ylm(9)
    G       = G + contrib*rlp1
 !   rlp1    = rlp1*invr
    ind     = 5
    do  l = 2,P_multipole
       contrib = multipoles(ind)*ylm(ind) ! 0.0_dp
       ind     = ind + 1    ! in peut l'enlever ?
       rlp1    = rlp1*invr
       do m= 2,2*l+1,2         ! -l ... l
          contrib = contrib + multipoles(ind)*ylm(ind)
          contrib = contrib + multipoles(ind+1)*ylm(ind+1)
          ind     = ind + 2
       end do
       G  = G + contrib*rlp1    !  G = G +  contrib/(r**(l+1))
    end do
    G = G*invr
    !
  end function G

  ! -------------------------------------------------------------------------
 !  seRHSandX
 !    set the rhs and x (the unknown) according to
 !    	[A_ii 0; 0 I] [x_i ; x_b] = [rho - A_ib u_0; u_0].
 !      where u_o  the values of x on the boundary
 !            x_b  the values where x is known (b = boundary)
 !            x_i  the values where x is unknown
 !  input
 !    rho the rhs of the linear system
 ! -------------------------------------------------------------------------
 !
 subroutine setRHSandX(rho,V)
   !
    use sys,      only: die
   implicit none
   !
   real(dp), intent(inout)  :: rho(1:nbSiestaValues),V(*)
   !
   real(dp), parameter   :: fourpi=16.0_dp*atan(1.0_dp) 
   real(dp), parameter   :: pi=4.0_dp*atan(1.0_dp) 
   integer               :: ierr,m,i,j,k,n_max
   integer               :: ind,indMG,index,skipy
   integer,save          :: myIter=0
   real(dp), allocatable :: values(:,:)
   real(dp)              :: invhx2,invhy2,invhz2,xp,yp,zp,errMax,coeff
!
#ifdef DEBUG
      call write_debug( '        PRE setRHSandX' )
#endif
#ifdef DEBUG_MG
    write(nf,*) "Begin setRHSandX"
#endif   
   CALL timer ('HYPRE_RHS',1)
   !
   ! Compute multipole expansion to approximate the solution on the boundary
   !
#ifdef MG_ADVANCED
   errMAx =0.0
   errMax = maxval(rho)
   do i = 1,nbSiestaValues
      coeff = max(coeff,rho(i))
   enddo
   print*,myIter,errMax,coeff
   myIter = myIter + 1
#endif
   call eval_multipole_coeffs(rho)
   !
   do i = 1, nbHypreValues
      hypreBuffer(i) =  0.0_dp
   end do
  ! INITIALIZE X AND B TO ZERO
   call HYPRE_StructVectorSetBoxValues(sol, hypreIndex(1,1), hypreIndex(1,2),  hypreBuffer(1), ierr)
    !---------------------------------------------------------------------------------------------------------------
   if(MGOrder == 4) then 
      call HYPRE_StructVectorSetBoxValues(rhs, hypreIndex(1,1), hypreIndex(1,2),  hypreBuffer(1), ierr)
      !       call HYPRE_StructVectorSetBoxValues(Brhs, hypreIndex(1,1), hypreIndex(1,2),  hypreBuffer(1), ierr)
   endif
   !---------------------------------------------------------------------------------------------------------------

! Set rho on the MG grid. 4 pi rho will be the RHS in hypre
!
#ifdef TEST_MG
   call cpu_time(start)
#endif
   m     = 1  
   skipy = shiftIndex(2)*NLSiestaPoints(1)
   index = 1 + shiftIndex(3)*NLSiestaPoints(1)*NLSiestaPoints(2)
   do k =  hypreIndex(3,1), hypreIndex(3,2)
      index = index + skipy 
      do j =  hypreIndex(2,1), hypreIndex(2,2)
         index = shiftIndex(1) + index
         do i = hypreIndex(1,1), hypreIndex(1,2)
#ifdef TEST_MG
            hypreBuffer(m) =  func_rhs(pos_x(shiftIndex(1) + i),pos_y(shiftIndex(2)  + j),pos_z(shiftIndex(3) +k))
#else
            hypreBuffer(m) =  fourpi*rho(index)
#endif
            m = m + 1 ; index = index + 1
         end do
      end do
   end do
   !
   !--------------------------------------------------------------------------------------------
   ! update the rhs
   !--------------------------------------------------------------------------------------------
   !
    if(MGOrder == 4) then 
      ! set in Brhs = 4 pi rho
      call HYPRE_StructVectorSetBoxValues(Brhs,hypreIndex(1,1), hypreIndex(1,2), hypreBuffer(1), ierr)
      call HYPRE_StructVectorAssemble(Brhs,ierr)
      ! set in rhs the boundary conditions
      call setBoundaryInRHSOrder4()
      call HYPRE_StructVectorAssemble(rhs,ierr)
      ! set in rhs :=  B* Brhs + rhs
      call HYPRE_StructMatrixMatvec(1.0_8,B,Brhs,1.0_8,rhs,ierr)
!      call hypre_structvectorprint(rhs,0,ierr)
    else if(MGOrder == 2) then 
      call HYPRE_StructVectorSetBoxValues(rhs,hypreIndex(1,1), hypreIndex(1,2), hypreBuffer(1), ierr)
      call setBoundaryInRHSOrder2()
      call HYPRE_StructVectorAssemble(rhs,ierr)
    end if	   
    !
    !-----------------------------------------------------------------------------------------------
#ifdef DEBUG_MG
    write(nf,*) "   Before  HYPRE_StructVectorAssemble"
#endif   
     call HYPRE_StructVectorAssemble(sol,ierr)
#ifdef DEBUG_MG
    write(nf,*) "End setRHSandX"
#endif   
    CALL timer ('HYPRE_RHS',2)
#ifdef DEBUG
      call write_debug( '        POS setRHSandX' )
#endif
!
    end subroutine setRHSandX

    subroutine setBoundaryInRHSOrder2()
      !
      integer               :: ierr,m,i,j,k,n_max
      integer               :: ilower(3),iupper(3),bc_ilower(3),bc_iupper(3),ind,indMG,index,NL(3)
      real(dp)              :: invhx2,invhy2,invhz2,xp,yp,zp,error,errormax
   !
   !  Put the boundary conditions in b  
   ! then we just need to solve in the interior:
   !                A_ii x_i = b_i - A_ib u_0.
   !
   invhx2 = 1.0/(H(1)*H(1))
   invhy2 = 1.0/(H(2)*H(2))
   invhz2 = 1.0/(H(3)*H(3))
   !
   ! Processors at z = 0  Bottom 
   ! BOTTOM  x= [dx,L-dx], y = [dy,L-dy],  z = dz 
   ! 
     if(bord(1,3)) then
        bc_ilower(1:3) = siestaIndex(1:3,1) 
        bc_iupper(1:2) = siestaIndex(1:2,2)
!        write(nf,*) 'RHS BOTTOM Z = dz ;  get the values of V on Z = 0 '
!        write(nf,*) '  bc_ilower : ', bc_ilower
!        write(nf,*) '  bc_iupper : ', bc_iupper
        if (bord(1,1) ) then          ! x = 0 
           bc_ilower(1) = 2
        end if
        if (bord(1,2) ) then        ! y = 0 
           bc_ilower(2) = 2
        end if
        !
        m = 1
        do  j =  bc_ilower(2), bc_iupper(2)
           do  i =  bc_ilower(1),  bc_iupper(1)
!              write(101,'(2I3,4E14.5)') i,j,pos_x(i),pos_y(j),pos_z(1), G(pos_x(i),pos_y(j),pos_z(1))
              hypreBuffer(m) = G(pos_x(i),pos_y(j),pos_z(1))*invhz2  ! inner plan z=0 
              m = m + 1
           end do
        end do
        bc_ilower(1:2) = hypreIndex(1:2,1) 
        bc_iupper(1:2) = hypreIndex(1:2,2)
        bc_ilower(3) = 1          !   z = dz 
        bc_iupper(3) = 1 
!        write(nf,*) 'RHS BOTTON Z = dz add in vect'
!        write(nf,*) '  bc_ilower : ', bc_ilower
!        write(nf,*) '  bc_iupper : ', bc_iupper
!        write(nf,*) '  m = ',m
        call HYPRE_StructVectorAddToBoxValue(rhs, bc_ilower, bc_iupper,  hypreBuffer(1), ierr)
     end if
!
!
! TOP x= [dx,L-dx], y = [dy,L-dy],  z = L-dz 
!    remember that L-dx and L-dy is already the size of the DFT Box
!
     if(bord(2,3)) then   !     
        bc_ilower(1:2) = siestaIndex(1:2,1) 
        bc_iupper(1:2) = siestaIndex(1:2,2)

        if (bord(1,1) ) then          ! x = 0 
           bc_ilower(1) = 2
        end if
        if (bord(1,2) ) then        ! y = 0 
           bc_ilower(2) = 2
        end if
        !
!!$        write(nf,*) 'RHS TOP Z = L-dz '
!!$        write(nf,*) '  bc_ilower : ', bc_ilower
!!$        write(nf,*) '  bc_iupper : ', bc_iupper
        m = 1
        do  j =  bc_ilower(2), bc_iupper(2)
           do  i =  bc_ilower(1),  bc_iupper(1)
!              write(102,'(2I3,3E14.5)') i,j,pos_x(i),pos_y(j),pos_z(siestaIndex(3,2)+1),G(pos_x(i),pos_y(j),pos_z(siestaIndex(3,2)+1))
              hypreBuffer(m) = G(pos_x(i),pos_y(j),pos_z(siestaIndex(3,2)+1))*invhz2  ! inner plan z=L 
              m = m + 1
           end do
        end do
        bc_ilower(1:2) = hypreIndex(1:2,1) 
        bc_iupper(1:2) = hypreIndex(1:2,2)
        bc_ilower(3)   = hypreIndex(3,2)
        bc_iupper(3)   = hypreIndex(3,2)
!        write(nf,*) 'RHS TOP Z = L-dz add in vect '
!        write(nf,*)' bc_ilower : ', bc_ilower
!        write(nf,*)' bc_iupper : ', bc_iupper
!        write(nf,*) '  m = ',m
        call HYPRE_StructVectorAddToBoxValue(rhs, bc_ilower, bc_iupper,  hypreBuffer(1),ierr)
     end if

!
!  LEFT x = [dx,L-dx],  y = dy, z= [dz,L-dz], */
!
     if (bord(1,2)) then                !  y = 0  
        bc_ilower(1:3) = siestaIndex(1:3,1) 
        bc_iupper(1:3) = siestaIndex(1:3,2)     
        if (bord(1,1)) then             !  x = 0
           bc_ilower(1) = 2
        end if 
        if (bord(1,3)) then             !  z = 0
           bc_ilower(3) = 2 
        end if
        ! Get the values of b on the surface
!!$        write(nf,*) 'RHS LEFT y= dy          X           Z'
!!$        write(nf,*)' bc_ilower : ', bc_ilower(1),bc_ilower(3)
!!$        write(nf,*)' bc_iupper : ', bc_iupper(1),bc_iupper(3)
        m = 1
        do j =  bc_ilower(3), bc_iupper(3)
           do i = bc_ilower(1), bc_iupper(1)
!              write(nf,'(2I3,3E14.5)') i,j,pos_x(i),pos_y(1),pos_z(j)
              hypreBuffer(m) = G(pos_x(i),pos_y(1), pos_z(j))*invhy2
              m = m + 1
           end do
        end do
        bc_ilower(1:3) = hypreIndex(1:3,1) 
        bc_iupper(1:3) = hypreIndex(1:3,2)
        bc_ilower(2) = 1 ; bc_iupper(2) = 1
        
!!$        write(nf,*) 'RHS LEFT y= dy add in vect '
!!$        write(nf,*)' bc_ilower : ', bc_ilower
!!$        write(nf,*)' bc_iupper : ', bc_iupper
        !
        call HYPRE_StructVectorAddToBoxValue(rhs, bc_ilower, bc_iupper,  hypreBuffer(1),ierr)
     end if
!
!  RIGHT x = [dx,L-dx],  y = L-dy, z= [dz,L-dz], */
!
     if (bord(2,2)) then                    !  y = L
        bc_ilower(1:3) = siestaIndex(1:3,1) 
        bc_iupper(1:3) = siestaIndex(1:3,2)     
        if (bord(1,1)) then                 !  x = 0
           bc_ilower(1) = 2
        end if 
        if (bord(1,3)) then                 !  z = 0
           bc_ilower(3) = 2
        end if
!!$        write(nf,*) 'RHS RIGHT y= L-dy       X           Z'
!!$        write(nf,*)' bc_ilower : ', bc_ilower(1), bc_ilower(3)
!!$        write(nf,*)' bc_iupper : ', bc_iupper(1),bc_iupper(3)
        m = 1
        do j =  bc_ilower(3), bc_iupper(3)
           do i = bc_ilower(1), bc_iupper(1)
              !              write(nf,'(2I3,3E14.5)') i,j,pos_x(i),pos_y(siestaIndex(2,2)+1),pos_z(j)
              hypreBuffer(m) = G(pos_x(i),pos_y(siestaIndex(2,2)+1), pos_z(j))*invhy2
              m = m + 1
           end do
        end do
        bc_ilower(1:3) = hypreIndex(1:3,1) 
        bc_iupper(1:3) = hypreIndex(1:3,2)       
        bc_ilower(2)   = hypreIndex(2,2)  
        bc_iupper(2)   = hypreIndex(2,2) 
        !
!!$        write(nf,*) 'RHS RIGHT Y = L-dy add in vect '
!!$        write(nf,*)' bc_ilower : ', bc_ilower
!!$        write(nf,*)' bc_iupper : ', bc_iupper
        !
        call HYPRE_StructVectorAddToBoxValue(rhs, bc_ilower, bc_iupper,  hypreBuffer(1),ierr)
     end if
!
!  BACK x =  dx,  y = [dy,L-dy], zx= [dz,L-dz]
!
     if(bord(1,1)) then                   ! x = 0
        bc_ilower(1:3) = siestaIndex(1:3,1) 
        bc_iupper(1:3) = siestaIndex(1:3,2)
         if (bord(1,3) ) then              ! z = 0 
           bc_ilower(3) = 2
        end if
        if (bord(1,2)) then               ! y = 0
           bc_ilower(2) = 2
        end if
!
        write(nf,*) 'RHS BACK x= dx           Y           Z'
        write(nf,*) '  bc_ilower : ', bc_ilower(2), bc_ilower(3)
        write(nf,*) '  bc_iupper : ', bc_iupper(2), bc_iupper(3)
        m = 1
        do j =  bc_ilower(3), bc_iupper(3)
           do i = bc_ilower(2), bc_iupper(2)
!               write(nf,'(2I3,3E14.5)') i,j,pos_x(1),pos_y(i),pos_z(j)
             hypreBuffer(m) = G(pos_x(1),pos_y(i), pos_z(j))*invhx2
              m = m +1
           end do
        end do
        bc_ilower(2:3) = hypreIndex(2:3,1) 
        bc_iupper(2:3) = hypreIndex(2:3,2)
        bc_ilower(1) = 1 ; bc_iupper(1) = 1
        !
        write(nf,*) 'RHS BACK x= dx'
        write(nf,*)' bc_ilower : ', bc_ilower
        write(nf,*)' bc_iupper : ', bc_iupper
        !
        call HYPRE_StructVectorAddToBoxValue(rhs, bc_ilower, bc_iupper,  hypreBuffer(1),ierr)
     end if
!
!  FRONT x =  L - dx,  y = [dy,L-dy], zx= [dz,L-dz]
!
     if(bord(2,1)) then   !! ????????????????
        bc_ilower(1:3) = siestaIndex(1:3,1) 
        bc_iupper(1:3) = siestaIndex(1:3,2)
         if (bord(1,3) ) then              ! z = 0 
           bc_ilower(3) = 2
        end if
        if (bord(1,2)) then               ! y = 0
           bc_ilower(2) = 2
        end if
        write(nf,*) 'RHS FRONT x= dx          Y           Z'
        write(nf,*) '  bc_ilower : ', bc_ilower(2), bc_ilower(3)
        write(nf,*) '  bc_iupper : ', bc_iupper(2), bc_iupper(3)
        !
        m = 1
        do j =  bc_ilower(3), bc_iupper(3)
           do i = bc_ilower(2), bc_iupper(2)
 !              write(nf,'(2I3,3E14.5)') i,j,pos_x(siestaIndex(1,2)+1),pos_y(i),pos_z(j)
             hypreBuffer(m) = G(pos_x(siestaIndex(1,2)+1),pos_y(i), pos_z(j))*invhx2
              m = m +1
           end do
        end do
        bc_ilower(2:3) = hypreIndex(2:3,1) 
        bc_iupper(2:3) = hypreIndex(2:3,2)
        bc_ilower(1) = hypreIndex(1,2) ; bc_iupper(1) = hypreIndex(1,2) 
        write(nf,*) 'RHS FRONT x= L-dx add in vect '
        write(nf,*)' bc_ilower : ', bc_ilower
        write(nf,*)' bc_iupper : ', bc_iupper
        call HYPRE_StructVectorAddToBoxValue(rhs, bc_ilower, bc_iupper,  hypreBuffer(1),ierr)
     end if

    end subroutine setBoundaryInRHSOrder2

    subroutine setBoundaryInRHSOrder4()
      !
      integer               :: ierr,m,i,j,k,n_max,is,js
      integer               :: ilower(3),iupper(3),bc_ilower(3),bc_iupper(3),ind,indMG,index,NL(3)
      integer               :: ishift,jshift,kshift
      real(dp), allocatable :: values(:,:)
      real(dp)              :: invhx2,invhy2,invhz2,xp,yp,zp,error,errormax,coeff1, coeff2
   !
   !  Put the boundary conditions in b  
   ! then we just need to solve in the interior:
   !                A_ii x_i = b_i - A_ib u_0.
   !
#ifdef DEBUG_MG
    write(nf,*) "       Before  setBoundaryInRHSOrder4"
#endif   
   invhx2 = 1.0_dp/(H(1)*H(1))
   invhy2 = 1.0_dp/(H(2)*H(2))
   invhz2 = 1.0_dp/(H(3)*H(3))
   coeff1 = 1.0_dp/3.0_dp
   coeff2 = 1.0_dp/12.0_dp
   ishift = 1  ! always 1 because we only consider Dirichlet conditions
   jshift = 1  ; kshift =1
   if( BoundaryConditions(2) > 0) jshift = 0 
   if( BoundaryConditions(3) > 0) kshift = 0 
   !
   !
   ! TOP x= [dx,L-dx], y = [dy,L-dy],  z = L-dz 
   !    remember that L-dx and L-dy is already the size of the DFT Box
   !
   if(bord(2,3)) then   !     
      bc_ilower(1:2) = siestaIndex(1:2,1) 
      bc_iupper(1:3) = siestaIndex(1:3,2) + 1 ! bc_iupper+1 car conditions periodiques dans siesta
      bc_ilower(3)   = bc_iupper(3)
      if (.not.(bord(1,1)) ) then          ! x = 0 
         bc_ilower(1) = siestaIndex(1,1) - 1
      end if
      if (.not.(bord(1,2)) ) then          ! y = 0 
         bc_ilower(2) = siestaIndex(2,1) - 1
      end if
    !  reshape(BoundaryZ0, /max(1,siestaIndex(1,1)-1):(siestaIndex(1,2)+1),max(1,siestaIndex(2,1)-1):(siestaIndex(2,2)+1)/ ))

!!$      write(nf,*) 'RHS TOP Z = L-dz         X           Y'
!!$      write(nf,*) '  bc_ilower : ', bc_ilower(1), bc_ilower(2)
!!$      write(nf,*) '  bc_iupper : ', bc_iupper(1), bc_iupper(2)
      do  j =  max(1,bc_ilower(2)), bc_iupper(2)
         do  i =  bc_ilower(1),  bc_iupper(1)
            BoundaryZ0(i,j) = G(pos_x(i),pos_y(j),pos_z(bc_iupper(3)))
         end do
      end do
      !
      bc_ilower(1:2) = hypreIndex(1:2,1) 
      bc_iupper(1:2) = hypreIndex(1:2,2)
      bc_ilower(3)   = hypreIndex(3,2)
      bc_iupper(3)   = hypreIndex(3,2)
!!$      write(nf,*) 'RHS TOP Z = L-dz add in vect '
!!$      write(nf,*)' bc_ilower : ', bc_ilower
!!$      write(nf,*)' bc_iupper : ', bc_iupper
      !
      m = 1
      do  js =  bc_ilower(2) + jshift, bc_iupper(2)  + jshift
         do  is =  bc_ilower(1) + ishift,  bc_iupper(1) + ishift
            hypreBuffer(m) = coeff1*invhz2*BoundaryZ0(is,js) + coeff2*(invhz2+invhx2)*(BoundaryZ0(is+1,js)+BoundaryZ0(is-1,js)) + &
                             coeff2*(invhz2+invhy2)*(BoundaryZ0(is,js+1)+BoundaryZ0(is,js-1))
            m = m + 1
         end do
      end do
      call HYPRE_StructVectorAddToBoxValue(rhs, bc_ilower, bc_iupper,  hypreBuffer(1),ierr)
   end if
   !
   ! Processors at z = 0  Bottom 
   ! BOTTOM  x= [dx,L-dx], y = [dy,L-dy],  z = dz 
   ! 
   if(bord(1,3)) then
!      reshape(BoundaryZ0, /max(1,siestaIndex(1,1)-1):(siestaIndex(1,2)+1),max(1,siestaIndex(2,1)-1):(siestaIndex(2,2)+1)/ ))
      bc_ilower(1:3) = siestaIndex(1:3,1) 
      bc_iupper(1:3) = siestaIndex(1:3,2) + 1 ! /!\ bc_iupper+1 car conditions periodiques
      if (.not.(bord(1,1)) ) then          ! x = 0 
         bc_ilower(1) = siestaIndex(1,1)-1
      end if
      if (.not.(bord(1,2)) ) then          ! y = 0 
          bc_ilower(2) = siestaIndex(2,1)-1
      end if
!!$      write(nf,*) 'RHS BOTTOM Z = dz              X           Y'
!!$      write(nf,*) ' Siesta bc_ilower : ', bc_ilower(1),bc_ilower(2)
!!$      write(nf,*) ' Siesta bc_iupper : ', bc_iupper(1),bc_iupper(2)
      !
!      write(nf,*) 'shape( Boundaryz0) ', shape( Boundaryz0)
      do  j =  bc_ilower(2), bc_iupper(2)
         do  i =  bc_ilower(1),  bc_iupper(1)
            BoundaryZ0(i,j) = G(pos_x(i),pos_y(j),pos_z(1))
         end do
      end do
      !
      bc_ilower(1:2) = hypreIndex(1:2,1) 
      bc_iupper(1:2) = hypreIndex(1:2,2)
      bc_ilower(3) = 1          !   z = dz 
      bc_iupper(3) = 1 
!!$      write(nf,*) 'RHS BOTTOM Z = dz add in vect'
!!$      write(nf,*) '  bc_ilower : ', bc_ilower
!!$      write(nf,*) '  bc_iupper : ', bc_iupper
      !
      m = 1
      do  js =  bc_ilower(2)+jshift, bc_iupper(2)+jshift
         do  is =  bc_ilower(1)+ishift,  bc_iupper(1)+ishift
            hypreBuffer(m) = coeff1*BoundaryZ0(is,js)*invhz2+coeff2*(invhz2+invhx2)*(BoundaryZ0(is+1,js)+BoundaryZ0(is-1,js)) + &
                             coeff2*(invhz2+invhy2)*(BoundaryZ0(is,js+1)+BoundaryZ0(is,js-1))
            m = m + 1
         end do
      end do
      !
      call HYPRE_StructVectorAddToBoxValue(rhs, bc_ilower, bc_iupper,  hypreBuffer(1), ierr)
   end if
!
!  RIGHT x = [dx,L-dx],  y = L-dy, z= [dz,L-dz], */
!
   if (bord(2,2)) then                    !  y = L
      bc_ilower(1:3) = siestaIndex(1:3,1) 
      bc_iupper(1:3) = siestaIndex(1:3,2)+1 ! bc_iupper+1 car conditions periodiques 
      if (.not.(bord(1,1)) ) then          ! x = 0 
          bc_ilower(1) = siestaIndex(1,1) - 1
      end if
      if (.not.(bord(1,3)) ) then          ! z = 0 
         bc_ilower(3) = siestaIndex(3,1) - 1
      end if
!!$      write(nf,*) 'RHS RIGHT y= L-dy       X           Z'
!!$      write(nf,*)' bc_ilower : ', bc_ilower(1), bc_ilower(3)
!!$      write(nf,*)' bc_iupper : ', bc_iupper(1),bc_iupper(3)

      do j =  bc_ilower(3), bc_iupper(3)
         do i = bc_ilower(1), bc_iupper(1)
            !              write(nf,'(2I3,3E14.5)') i,j,pos_x(i),pos_y(siestaIndex(2,2)+1),pos_z(j)
            BoundaryY0(i,j) = G(pos_x(i),pos_y(siestaIndex(2,2)+1),pos_z(j))
         end do
      end do
      bc_ilower(1:3) = hypreIndex(1:3,1) 
      bc_iupper(1:3) = hypreIndex(1:3,2)       
      bc_ilower(2)   = hypreIndex(2,2) 
      bc_iupper(2)   = hypreIndex(2,2)
!!$      write(nf,*) 'RHS RIGHT Y = L-dy add in vect '
!!$      write(nf,*)' bc_ilower : ', bc_ilower
!!$      write(nf,*)' bc_iupper : ', bc_iupper
      m = 1
      do  js =  bc_ilower(3) + kshift, bc_iupper(3) + kshift
         do  is =  bc_ilower(1) + ishift,  bc_iupper(1) + ishift
            hypreBuffer(m) = coeff1*invhy2*BoundaryY0(is,js) + coeff2*(invhy2+invhx2)*(BoundaryY0(is+1,js)+BoundaryY0(is-1,js)) + &
                             coeff2*(invhz2+invhy2)*(BoundaryY0(is,js+1)+BoundaryY0(is,js-1))
            m = m + 1
         end do
      end do
      call HYPRE_StructVectorAddToBoxValue(rhs, bc_ilower, bc_iupper,  hypreBuffer(1),ierr)
   end if
   !
   !  LEFT x = [dx,L-dx],  y = dy, z= [dz,L-dz], */
   !
   if (bord(1,2)) then                !  y = 0  
      bc_ilower(1:3) = siestaIndex(1:3,1) 
      bc_iupper(1:3) = siestaIndex(1:3,2) + 1 ! bc_iupper+1 car conditions periodiques  
      if (.not.(bord(1,1)) ) then          ! x = 0 
          bc_ilower(1) = siestaIndex(1,1) - 1
      end if
      if (.not.(bord(1,3)) ) then          ! z = 0 
          bc_ilower(3) = siestaIndex(3,1) - 1
      end if
!!$      write(nf,*) 'RHS LEFT y= dy          X           Z'
!!$      write(nf,*)' bc_ilower : ', bc_ilower(1), bc_ilower(3)
!!$      write(nf,*)' bc_iupper : ', bc_iupper(1),bc_iupper(3)
     do j =  bc_ilower(3), bc_iupper(3)
         do i = bc_ilower(1), bc_iupper(1)
            BoundaryY0(i,j) = G(pos_x(i),pos_y(1),pos_z(j))
         end do
      end do
        !
      bc_ilower(1:3) = hypreIndex(1:3,1) 
      bc_iupper(1:3) = hypreIndex(1:3,2)
      bc_ilower(2) = 1 
      bc_iupper(2) = 1
!!$      write(nf,*) 'RHS LEFT y= dy add in vect '
!!$      write(nf,*)' bc_ilower : ', bc_ilower
!!$      write(nf,*)' bc_iupper : ', bc_iupper
      !
      m = 1
      do  js =  bc_ilower(3) + kshift, bc_iupper(3) + kshift
         do  is =  bc_ilower(1) + ishift,  bc_iupper(1) + ishift
            hypreBuffer(m) = coeff1*invhy2*BoundaryY0(is,js) + coeff2*(invhy2+invhx2)*(BoundaryY0(is+1,js)+BoundaryY0(is-1,js)) + & 
                             coeff2*(invhz2+invhy2)*(BoundaryY0(is,js+1)+BoundaryY0(is,js-1))
            m = m + 1
         end do
      end do
      call HYPRE_StructVectorAddToBoxValue(rhs, bc_ilower, bc_iupper,  hypreBuffer(1),ierr)
   end if
!
!  BACK x =  dx,  y = [dy,L-dy], zx= [dz,L-dz]
!
   if(bord(1,1)) then                   ! x = 0
      bc_ilower(1:3) = siestaIndex(1:3,1) 
      bc_iupper(1:3) = siestaIndex(1:3,2) + 1 ! bc_iupper+1 car conditions periodiques 
      if (.not.(bord(1,3)) ) then              ! z = 0 
          bc_ilower(3) = siestaIndex(3,1) - 1
      end if
      if (.not.(bord(1,2))) then               ! y = 0
          bc_ilower(2) = siestaIndex(2,1) - 1
      end if
!!$      write(nf,*) 'RHS BACK   X = dx        Y           Z'
!!$      write(nf,*) '  bc_ilower : ', bc_ilower(2),bc_ilower(3)
!!$      write(nf,*) '  bc_iupper : ', bc_iupper(2), bc_iupper(3)
      do j =  bc_ilower(3), bc_iupper(3)
         do i = bc_ilower(2), bc_iupper(2)
            !               write(nf,'(2I3,3E14.5)') i,j,pos_x(1),pos_y(i),pos_z(j)
            BoundaryX0(i,j) = G(pos_x(1),pos_y(i),pos_z(j))
         end do
      end do
      bc_ilower(2:3) = hypreIndex(2:3,1) 
      bc_iupper(2:3) = hypreIndex(2:3,2)
      bc_ilower(1) = 1
      bc_iupper(1) = 1
!!$      write(nf,*) 'RHS BACK  x= dx add in vect '
!!$      write(nf,*)' bc_ilower : ', bc_ilower
!!$      write(nf,*)' bc_iupper : ', bc_iupper
      !
      m = 1
      do  js =  bc_ilower(3) + kshift, bc_iupper(3) + kshift
         do  is =  bc_ilower(2) + jshift,  bc_iupper(2) + jshift
            hypreBuffer(m) = coeff1*invhx2*BoundaryX0(is,js) + coeff2*(invhy2+invhx2)*(BoundaryX0(is+1,js)+BoundaryX0(is-1,js)) + & 
                             coeff2*(invhz2+invhx2)*(BoundaryX0(is,js+1)+BoundaryX0(is,js-1))
            m = m + 1
         end do
      end do
      call HYPRE_StructVectorAddToBoxValue(rhs, bc_ilower, bc_iupper,  hypreBuffer(1),ierr)
   end if
!
!  FRONT x =  L - dx,  y = [dy,L-dy], zx= [dz,L-dz]
!
   if(bord(2,1)) then
      bc_ilower(1:3) = siestaIndex(1:3,1) 
      bc_iupper(1:3) = siestaIndex(1:3,2) + 1 ! bc_iupper+1 car conditions periodiques
      if (.not.(bord(1,3)) ) then              ! z = 0 
          bc_ilower(3) = siestaIndex(3,1) - 1
      end if
      if (.not.(bord(1,2))) then               ! y = 0
          bc_ilower(2) = siestaIndex(2,1) - 1
      end if
!!$      write(nf,*) 'RHS FRONT x= L-dx       Y           Z'
!!$      write(nf,*)' bc_ilower : ', bc_ilower(2),bc_ilower(3)
!!$      write(nf,*)' bc_iupper : ', bc_iupper(2),bc_iupper(3)
      do j =  bc_ilower(3), bc_iupper(3)
         do i = bc_ilower(2), bc_iupper(2)
            !              write(nf,'(2I3,3E14.5)') i,j,pos_x(siestaIndex(1,2)+1),pos_y(i),pos_z(j)
            BoundaryX0(i,j) = G(pos_x(siestaIndex(1,2)+1),pos_y(i),pos_z(j))
         end do
      end do
      bc_ilower(2:3) = hypreIndex(2:3,1) 
      bc_iupper(2:3) = hypreIndex(2:3,2)
      bc_ilower(1)   = hypreIndex(1,2)
      bc_iupper(1)   = hypreIndex(1,2)
!!$      write(nf,*) 'RHS FRONT x= L-dx add in vect '
!!$      write(nf,*)' bc_ilower : ', bc_ilower
!!$      write(nf,*)' bc_iupper : ', bc_iupper
      !
      m = 1
      do  js =  bc_ilower(3) + kshift, bc_iupper(3) + kshift
         do  is =  bc_ilower(2) + jshift,  bc_iupper(2) + jshift
            hypreBuffer(m)=coeff1*BoundaryX0(is,js)*invhx2 + &
              coeff2*(BoundaryX0(is+1,js)+BoundaryX0(is-1,js)) * &
              (invhy2+invhx2) + &
              coeff2*(BoundaryX0(is,js+1)+BoundaryX0(is,js-1))*(invhz2+invhx2)
            m = m + 1
         end do
      end do
      call HYPRE_StructVectorAddToBoxValue(rhs, bc_ilower, bc_iupper,  hypreBuffer(1),ierr)
   end if
#ifdef DEBUG_MG
    write(nf,*) "       End  setBoundaryInRHSOrder4"
#endif  
 end subroutine setBoundaryInRHSOrder4
  !
  !
  ! -------------------------------------------------------------------------
  ! hypre_init_solver  
  !    intialise the PFMG solver in Hypre (a semicoarsening multigrid solver)
  ! -------------------------------------------------------------------------
  !
  subroutine hypre_poisson()
!
   use parallel,    only : Node,IOnode
   use sys,         only : die
#ifdef MPI
      use mpi_siesta
#endif
    implicit none
    integer :: ierr
    integer :: num_iterations
    real(8) :: final_res_norm
#ifdef MPI
      integer                   :: MPIerror
#endif
    !
#ifdef DEBUG_MG
    write(nf,*) "Begin hypre_poisson"
#endif 
#ifdef DEBUG
      call write_debug( '      PRE hypre_poisson' )
#endif
    CALL timer ('HYPRE_POISSON',1)
!
    if(firstSolver ) then   
       !       write(nf,'(A)') "        INIT SOLVER  " 
       call hypre_init_solver()
       firstSolver = .false.
    end if
#ifdef DEBUG_MG
    write(nf,*) "     start HYPRE_StructXXXSolve"
#endif
    if( solver_id == 0 ) then 
       call HYPRE_StructPFMGSolve(solver, A, rhs, sol, ierr)
       call HYPRE_StructPFMGGetFinalRelativ(solver, final_res_norm,ierr)
    else  if( solver_id == 1 ) then 
       call HYPRE_StructPCGSolve(solver, A, rhs, sol, ierr)
       call HYPRE_StructPCGGetFinalRelative(solver, final_res_norm,ierr)
    endif
!
    if( final_res_norm > tolerance ) then
       write(*,'(A,i3,E10.3)') "Convergence is not reached in Poisson solver. Stop",ierr, final_res_norm
       write(nf,*) "     END HYPRE_StructXXXSolve"
       write(*,*) ' Save data '
       call hypre_save_solution()
!   Does we abort ?
       call die("Stop the computation ...") 
    end if
!  Get info
    if(mgVerbose > 0 .and. IOnode ) then 
       if(solver_id == 0 ) then 
          call HYPRE_StructPFMGGetNumIteration(solver, num_iterations,ierr)
       else  if(solver_id == 1 ) then 
          call HYPRE_StructPCGGetNumIterations(solver, num_iterations,ierr)
       end if
       !
       write(*,'(A,I4,A,E10.3)')   ' hypre: iteration number: ',num_iterations, &
            '  Relative Residual Norm: ',final_res_norm
    end if
    !
!
    CALL timer ('HYPRE_POISSON',2)
#ifdef DEBUG_MG
    write(nf,*) "End hypre_poisson"
#endif   
!
#ifdef DEBUG
      call write_debug( '      POS hypre_poisson' )
#endif
  end subroutine hypre_poisson
  ! 
  ! -------------------------------------------------------------------------
  !
  subroutine hypre_init_solver()
    !
    !
    !  Options and setup 
    !
    implicit none
    integer ierr
#ifdef DEBUG_MG
    write(nf,*) "Begin hypre_init_solver"
#endif 
#ifdef DEBUG
    call write_debug( '        PRE hypre_init_solver' )
#endif
    if (solver_id ==0 ) then 
       call HYPRE_StructPFMGCreate(hypre_mpi_comm, solver, ierr)
       call HYPRE_StructPFMGSetMaxIter(solver, maxIterations, ierr)
       call HYPRE_StructPFMGSetTol(solver, tolerance, ierr)
       call HYPRE_StructPFMGSetRelChange(solver, 0, ierr)
       call HYPRE_StructPFMGSetNumPreRelax(solver, num_pre, ierr)
       call HYPRE_StructPFMGSetNumPostRelax(solver, num_post, ierr)
       call HYPRE_StructPFMGSetLogging(solver, 1, ierr)
       if (mgVerbose > 0) call HYPRE_StructPFMGSetPrintLevel(solver, 1, ierr)
       call HYPRE_StructPFMGSetup(solver, A, rhs, sol, ierr)
    else   if (solver_id ==1 ) then 
       call HYPRE_StructPCGCreate(hypre_mpi_comm, solver, ierr)
       call HYPRE_StructPCGSetMaxIter(solver, maxIterations, ierr)
       call HYPRE_StructPCGSetTol(solver, tolerance, ierr)
       call HYPRE_StructPCGSetTwoNorm( solver, 1 ,ierr)
       call  HYPRE_StructPCGSetRelChange( solver, 0 , ierr)
       call HYPRE_StructPCGSetLogging(solver, 1, ierr)
       !------------------------------------------------------------
       !   * The precond_id flags mean :
       !    * 0 - setup a smg preconditioner
       !    * 1 - setup a pfmg preconditioner
       !    * 7 - setup a jacobi preconditioner
       !    * 8 - setup a ds preconditioner
       !    * 9 - dont setup a preconditioner
       !------------------------------------------------------------
       precond =0 
       if (precond_id == 1 ) then
          call HYPRE_StructPFMGCreate(hypre_mpi_comm, precond, ierr) 
          call HYPRE_StructPFMGSetMaxIter(precond, precondIter, ierr) 
          call HYPRE_StructPFMGSetTol(precond,tolerance, ierr)
          call HYPRE_StructPFMGSetZeroGuess(precond,ierr) 
          call HYPRE_StructPFMGSetNumPreRelax(precond, num_pre, ierr) 
          call HYPRE_StructPFMGSetNumPostRelax(precond,num_post, ierr) 
          call HYPRE_StructPFMGSetLogging(precond, 0, ierr) 
       else
          precond_id = 9 ! No preconditioner 
       end if
       call HYPRE_StructPCGSetPrecond(solver,precond_id, precond, ierr)
       if (mgVerbose > 0) call  HYPRE_StructPCGSetPrintLevel(solver, 1, ierr)
       call HYPRE_StructPCGSetup(solver, A, rhs, sol, ierr)
    end if
#ifdef DEBUG_MG
    write(nf,*) "End hypre_init_solver"
#endif  
#ifdef DEBUG
      call write_debug( '        POS hypre_init_solver' )
#endif  
 !
  end subroutine hypre_init_solver
  !
  ! -------------------------------------------------------------------------
  !  hypre_cleanGridAndMatrix 
  !       clean the grid and the matrix in HYPRE. must be called when the 
  !         siesta grid changes.
  ! -------------------------------------------------------------------------
  !  
  subroutine hypre_cleanGridAndMatrix()
    implicit none 
    !
    integer ierr
    !
#ifdef DEBUG
      call write_debug( '    PRE hypre_cleanGridAndMatrix' )
#endif
#ifdef DEBUG_MG
    write(nf,*) "Begin hypre_cleanGridAndMatrix"
#endif 
    CALL timer ('HYPRE_CLEAN',1)
    firstSolver = .true.
    if(firstTime) then
       firstTime = .false.
       CALL timer ('HYPRE_CLEAN',2)
#ifdef DEBUG
      call write_debug( '    POS hypre_cleanGridAndMatrix' )
#endif
       return
    end if
    if(allocated(pos_x) ) then 
!       write(*,*) ' deallocate pos_x'
       deallocate (pos_x,pos_y,pos_z)
    end if
    call HYPRE_StructGridDestroy(grid, ierr)
!    write(*,*) ' HYPRE_StructGridDestroy ', ierr
    call HYPRE_StructMatrixDestroy(A, ierr)
!    write(*,*) ' HYPRE_StructMatrixDestroy ', ierr
    if(MGOrder == 4) then 
      call HYPRE_StructMatrixDestroy(B, ierr)
      deallocate(BoundaryX0,BoundaryY0,BoundaryZ0)
    end if
    call HYPRE_StructVectorDestroy(Brhs, ierr)
    call HYPRE_StructVectorDestroy(rhs, ierr)
!    write(*,*) ' HYPRE_StructVectorDestroy ', ierr

    call HYPRE_StructVectorDestroy(sol, ierr)
!    write(*,*) ' HYPRE_StructVectorDestroy ', ierr
    if (solver_id == 0 ) then    
       call HYPRE_StructPFMGDestroy(solver, ierr)
    else      if (solver_id == 1 ) then    
       call HYPRE_StructPCGDestroy(solver, ierr)
       if (precond_id == 1) call HYPRE_StructPFMGDestroy(precond, ierr)
    else     if (solver_id == 2 ) then    
       call HYPRE_StructGMRESDestroy(solver, ierr)
       if (precond_id == 1) call HYPRE_StructPFMGDestroy(precond, ierr)
    end if
!    write(*,*) ' HYPRE_StructVectorDestroy ', ierr
    if(allocated(hypreBuffer) ) then
!    write(*,*) ' deallocate hypreBuffer '
       deallocate(hypreBuffer)
    end if
!
#ifdef DEBUG_MG
    write(*,*) "Begin hypre_cleanGridAndMatrix"
#endif 
    CALL timer ('HYPRE_CLEAN',2)
!
#ifdef DEBUG
      call write_debug( '    POS hypre_cleanGridAndMatrix' )
#endif
  end subroutine hypre_cleanGridAndMatrix
  !
  ! -------------------------------------------------------------------------
  !  This subroutine uses multipole expansion to compute the
  !  electrostatic potential at the boundaries of the cell that contains
  !  a known localized distribution of charge.
  !  As rho is known on the siesta grid, the loops are done with the associated index.
  ! -------------------------------------------------------------------------
  !
  subroutine eval_multipole_coeffs(rho)
    !
    use parallel,  only: node
    implicit none
!  --- Input arguments ---
    REAL(dp), intent(in) :: rho(*)
!  --- Internal definitions ---
    INTEGER(4) :: i,j,k,ind, l, m, ind1,ind2,indm0
    REAL(dp)   :: dvol
    REAL(dp)   :: x1,x2,x3,r,coeff
#ifdef MPI
    integer                :: ierr
#endif
     !
    !  Element of volume associated to each point
    dvol = H(1)*H(2)*H(3)
    !
    !  Initialize coefficients of multipole expansion  
    multipoles(:)  = 0.0_dp
    !
    ! Compute coefficients of multipole expansion
    ! May be we can win a little if we use the structure of non zero rho
    ! 
    ind   = 1 
    do k = siestaIndex(3,1), siestaIndex(3,2)
       x3 = pos_z(k) - centre(3)
       do j = siestaIndex(2,1), siestaIndex(2,2)
          x2 = pos_y(j) - centre(2)
          do i= siestaIndex(1,1), siestaIndex(1,2)
             x1 = pos_x(i) - centre(1)
!             if(rho(ind) .GT. 0.0_8) then 
                coeff = rho(ind)
                call evalylm(x1,x2,x3,ylm)  ! compute the spherical harmonic at point (x1,x2,x3)
                r     = sqrt(x1*x1+x2*x2+x3*x3)
                !   multipoles += rho(x) |x-centre|^l Y_l^m(theta,phi)
                !   coeff :=  rho(x) r^l and r =  |x-centre|
                do l = 0, P_multipole
                   indm0  =  l*l+l+1
                   !  m = 0
                   multipoles(indm0) = multipoles(indm0) + coeff*ylm(indm0)
                   do m = 1,l   !
                      multipoles(indm0+m) = multipoles(indm0+m) + coeff*ylm(indm0+m)  
                      multipoles(indm0-m) = multipoles(indm0-m) + coeff*ylm(indm0-m)
                   end do
                   coeff = coeff*r   ! multiply by r 
                end do
!             endif
             ind = ind + 1 
          end do
       end do
    end do
    !
    multipoles(:) = dvol*multipoles(:)
#ifdef MPI
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, multipoles(1),num_coeffMultipole,  MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)  
#endif
! 
#ifdef __DEBUG_NOSSI_MG_MULTIPOLES__
    if(firstMultip .AND. node==0 ) then 
       write(*,*) "Spherical multipole coefficients" 
       i=1
       do j = 0, P_multipole 
          write(*,'(A,I3)') ' L: ',j
          do k=-j,j
             write(*,'(I6,I4,E15.5)') i ,k, multipoles(i)
             i = i + 1 
          end do
       end do
       firstMultip = .false.
    end if
#endif
  end subroutine eval_multipole_coeffs
  !
  ! -----------------------------------------------------
  !
  subroutine hypre_checkBox(maxOrbitalRadius)
    use atmfuncs,     only: rcut
    use siesta_geom,  only: xa,isa,ucell
    use parallel,     only: Node,IOnode
    use fdf,          only: fdf_convfac
    !
    real(dp),intent(in) :: maxOrbitalRadius
    integer  :: is
    real(dp) :: CornerLeft(3),CornerRight(3),coeff,length(3)
    logical  :: warn = .false.
    !
    if( IOnode ) then 
       write(*,*) "ATOMS"
       do is = 1, size(xa,2)
          write(*,'(I2,3X,3E13.5)') is, xa(1:3,is)
       end do
    end if
    !
    CornerLeft(:)  = minval(xa,DIM=2) 
    CornerRight(:) = maxval(xa,DIM=2) 
    !
    length(1) = ucell(1,1)
    length(2) = ucell(2,2)
    length(3) = ucell(3,3)
    !
    CornerLeft(:)  = CornerLeft(:)  - maxOrbitalRadius
    CornerRight(:) = CornerRight(:) + maxOrbitalRadius
    if(Node == 0) then      
       coeff = fdf_convfac('bohr', 'ang')
       !
       write(*,'(/,A)') " Check box for non periodic system."
       write(*,'(/,A,E15.8)')  "maxOrbitalRadius (Ang): ",maxOrbitalRadius*coeff
       
       write(*,'(A,3E15.8)')   "Minimal Box (rho support)" 
       write(*,'(A,3E15.8)') "    CornerRight (Ang):      ",CornerRight(:)*coeff
       write(*,'(A,3E15.8,/)')   "    CornerLeft  (Ang):      ",CornerLeft(:)*coeff
       write(*,'(A,3E15.8)')   "Size of rho support (Ang):  ",(CornerRight(:)-CornerLeft(:))*coeff
       write(*,'(A,3E15.8,/)') "CEll (Ang):                 ",Length(:)*coeff
       write(*,'(A,3E15.8)')   "Free distance (Ang):        ",(Length(:)-(CornerRight(:)-CornerLeft(:)))*coeff
       write(*,'(A,3E15.8)')   "     with corner right:     ",(Length(:)-CornerRight(:))*coeff
       write(*,'(A,3E15.8)')   "     with corner left:       ",(CornerLeft(:))*coeff
       write(*,'(/,A)') "--------------------------------------------------"
    end if
    warn = .false.
    if(minval(CornerLeft) <= 0 ) warn = .true.
    CornerRight(:) = CornerRight(:) - length(:)
    if(minval(CornerRight) >= 0 ) warn = .true.
    if(warn .and. IOnode) then 
       write(*,*) "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&"
       write(*,*) "&&&      WARNING  WARNING  WARNING        &&&"
       write(*,*) "&&&                                       &&&"
       write(*,*) "&&& Multipole expantions are not valid    &&&"
       write(*,*) "&&& Increase the  size of your box        &&&"
       write(*,*) "&&&   or check the position of the atoms  &&&"
       write(*,*) "&&&    int he box                         &&&"
       write(*,*) "&&&      WARNING  WARNING  WARNING        &&&"
       write(*,*) "&&&                                       &&&"
       write(*,*) "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&"
    endif
#ifdef MG_FULLGRID
    if  (IOnode) then
       write(*,*)  " **************** FULL MG GRID ***************"       
    endif
#endif
  end subroutine hypre_checkBox
  !
  ! -------------------------------------------------------------------------
  !          End of the nossi_mg_hypre module
  ! -------------------------------------------------------------------------
  !
  subroutine hypre_save_solution(name_in)
    !
    use siesta_geom, only: xa,isa
    use atmfuncs,    only: labelfis
    use parallel,    only: nodes,node
!
    implicit none
    character(*), optional, intent(in) :: name_in

    real(dp), allocatable :: VV(:)
    integer  :: i,j,k, count,ierr,ishifty,ishiftz,NN(3)
    integer  :: iu
    real(dp) :: Ori(3)
    character(200) :: NAME
    !
    call io_assign( iu )
!    write(*,*) 'Save atoms'
    open(iu,file='NOSSI_TRACE_atoms.xyz')
    write(iu,'(I5)') size(xa,2)
    write(iu,'(A)') 'NOSSI_TRACE_atoms'
    
    do i = 1 , size(xa,2)
       write(iu,'(A3,3E15.7)') labelfis(isa(i)),xa(1:3,i)
    end do
   call  io_close(iu)
   !
    allocate (VV(nbHypreValues) )
    ishifty = 0 ; ishiftz = 0
    nn(1:3) =  NGpoints(1:3)
    do i = 1, 3
       if(BoundaryConditions(i) == 0 ) nn(i) = nn(i) -1
    end do
!!$    if(BoundaryConditions(2) > 0) ishifty = 0
!!$    if(BoundaryConditions(3) > 0) ishiftz = 0
    if(BoundaryConditions(2) == 0) ishifty = 1
    if(BoundaryConditions(3) == 0) ishiftz = 1
    call HYPRE_StructVectorGetBoxValues(sol, hypreIndex(1,1), hypreIndex(1,2),  VV, ierr) 
    call io_assign( iu )
    if(present(name_in) ) then 
       write(NAME,'(a,a,I2.2,a)') trim(name_in),'_grid',node,'.vtk'
    else
       write(NAME,'(a,I2.2,a)') 'potential_grid',node,'.vtk'
    endif
    open(iu,file=trim(NAME))
    write(iu,'("# vtk DataFile Version 2.0",/,A,/,A)') 'NOSSI_TRACE_grid ','ASCII'
    write(iu,'(A)') 'DATASET STRUCTURED_POINTS'
    write(iu,'(A,3I7)') 'DIMENSIONS ', NLPoints(1:3)
    write(iu,'(A,3E13.5)') 'ORIGIN     ', pos_x(hypreIndex(1,1)+1), pos_y(hypreIndex(2,1)), pos_z(hypreIndex(3,1))
    write(iu,'(A,3E13.5)') 'SPACING    ',H(1:3)
!!$    write(iu,'(A)') 'DATASET STRUCTURED_GRID'
!!$    write(iu,'(A,3I7)') 'DIMENSIONS ', NLPoints(1:3)
!!$    write(iu,'(A,I10,A)') 'POINTS ', nbHypreValues ,'  double'
!!$    do k =hypreIndex(3,1), hypreIndex(3,2)
!!$       do j=hypreIndex(2,1), hypreIndex(2,2)
!!$          do i =  hypreIndex(1,1), hypreIndex(1,2)
!!$             write(iu,'(3E13.5)')  H(1) * i, H(2) * (j-ishifty), H(3) * (k-ishiftz)
!!$          end do
!!$        end do
!!$    end do
    write(iu,'(A,I10)') 'POINT_DATA ', nbHypreValues 
    write(iu,'(A)') 'SCALARS efield double 1'
    write(iu,'(A)') 'LOOKUP_TABLE default'
    !
    count = 1
    do k =hypreIndex(3,1), hypreIndex(3,2)
       do j=hypreIndex(2,1), hypreIndex(2,2)
          do i =  hypreIndex(1,1), hypreIndex(1,2)
             write(iu,'(E13.5)') VV(count)
             count = count + 1
          end do
       end do
    end do
    CALL  IO_CLOSE(IU)

!!$    if(Nodes ==1 ) then 
!!$       write(NAME,'(a,I2.2,a)') 'grid',node,'.vtr'
!!$    else
!!$       write(NAME,'(a,I2.2,a)') 'grid',node,'.pvtr'
!!$    endif
!!$    open(iu,file=trim(NAME))
!!$    write(iu,'(a)') <VTKFile type=”RectilinearGrid” >
!!$    write(iu,'(a)') </VTKFile>
  deallocate(VV)
    stop 'hkjhkjhkhk'
  end subroutine hypre_save_solution

#ifdef TEST_MG
  function solution(x,y,z)
    !  Computes the potential V(x) = rho * 1/|x-y| a given point on the boundary of the box
    !           through a multipole expansion.
    !
    !  See, for example, "Classical Electrodynamics", J D Jackson, chapter 4 (1998).
    use siesta_geom,  only: xa
    implicit none
   real(dp), parameter   :: twopi=8.0_dp*atan(1.0) 

    real(8), intent(in) :: x,y,z
    real(8) :: solution,xx,yy,zz
!
!    solution =  x*(Lx-x)*y*(Ly-y)*z*(Lz-z)/(lx*ly*lz)
    xx = x / lx; yy = y/Ly ; zz=z/Lz
    solution =  sin(twopi*xx)*sin(twopi*yy)*sin(twopi*zz)
!!$       return
    
     end function SOLUTION
 pure function func_rhs(x,y,z)
    !  Computes the potential V(x) = rho * 1/|x-y| a given point on the boundary of the box
    !           through a multipole expansion.
    !
    !  See, for example, "Classical Electrodynamics", J D Jackson, chapter 4 (1998).
    implicit none
   real(dp), parameter   :: twopi=8.0_dp*atan(1.0) 

    real(8), intent(in) :: x,y,z
    real(8) :: func_rhs,xx,yy,zz
!
!    func_rhs =  2.0_8*( y*(ly-y)*z*(lz-z) + x*(lx-x)* z*(lz-z) + x*(lx-x)*y*(ly-y) )/(lx*ly*lz)
!
    xx = x / lx; yy = y/Ly ; zz=z/Lz
    func_rhs =  twopi*twopi*(1.0_8/ (lx*lx) +1.0_8/(ly*ly) + 1.0_8/(lz*lz))*sin(twopi*xx)*sin(twopi*yy)*sin(twopi*zz)
  end function FUNC_RHS
  !
  ! Tools for Spherical harmonics
  !
  subroutine checkOnPlaneZ0()
    !
    integer :: i,j
    real(dp) :: eval, sol, err,error
    !
    write(*,*) 'Check G on all points'
    !
    ! Check G on all points
    !
    write(*,*) 'Error on plan z = 0 '

    err = 0.0
    do j = 1, NGPoints(2)+1
       do i= 1, NGPoints(1)+1
          !
          eval  = G( pos_x(i), pos_y(j), pos_z(1))
          sol   = solution( pos_x(i), pos_y(j), pos_z(1))
          error = abs(sol- eval)
          err   = max(err,error)
          write(252,'(2I4,7E15.5)')  i,j,pos_x(i), pos_y(j),  sol, eval,error,err
       end do
    end do
    write(*,*) 'Error with G_REAL    ', err

  end subroutine checkOnPlaneZ0

  subroutine compareSolution(V)
   !
   real(dp), intent(inout)  :: V(*)
   !
   integer :: m,index,i,j,k,nl(3),ierr,imax,jmax,kmax
   real(dp) :: xx,error,errorMax
   !
   !
#ifdef DEBUG_MG
    write(nf,*) "Begin compareSolution"
#endif 
   call cpu_time(finish)
   !
   ! set Solution on the boundary
   ! index : i + (j-1)*NL(1) + (k-1)*NL(1)*NL(2)
   NL(:) = siestaIndex(:,2) - siestaIndex(:,1) + 1
!
   if(bord(1,3)) then ! z = 0
      index = 1
      do j =  siestaIndex(2,1), siestaIndex(2,2)
         do i = siestaIndex(1,1), siestaIndex(1,2) 
            V(index) = solution(pos_x(i),pos_y(j),pos_z(1))
            index = index + 1
         end do
      end do
   end if
   if(bord(1,2)) then  ! y = 0
      do k =  siestaIndex(3,1), siestaIndex(3,2)
         do i = siestaIndex(1,1), siestaIndex(1,2) 
            index = i + NL(1)*NL(2)*(k-1)
            V(index) = solution(pos_x(i),pos_y(1),pos_z(k))
         end do
      end do
   end if

   if(bord(1,1)) then ! x = 0
      do k =  siestaIndex(3,1), siestaIndex(3,2)
         do j = siestaIndex(2,1), siestaIndex(2,2) 
            index = 1 + NL(1)*(j-1) + NL(1)*NL(2)*(k-1)
            V(index) = solution(pos_x(1),pos_y(j),pos_z(k))
         end do
      end do
   end if
   !
   errorMax = 0.0
   m = 1  
   index =  1
   do k =  siestaIndex(3,1), siestaIndex(3,2)
      do j =  siestaIndex(2,1), siestaIndex(2,2)
         do i = siestaIndex(1,1), siestaIndex(1,2)
               xx    = SOLUTION(pos_x(i),pos_y(j),pos_z(k))
               error = abs(V(index)- xx)
               if(error > 0.01) then 
                  write(123,'(A,3I4,2x,6E14.6)') "-->",pos_x(i),pos_y(j),pos_z(k),V(index) , xx, error
               end if
               write(122,'(A,3I4,2x,6E14.6)') "-->",i,j,k,pos_x(i),pos_y(j),pos_z(k),V(index) , xx, error
               if ( error > errormax ) then 
                  errormax = error ; imax = i ; jmax =j ; kmax = k
               end if
 !           end if
            m = m + 1 ; index = index + 1
         end do
      end do
   end do
   close(122)
   print*, ' '
   print '("Time:      ",f6.3," seconds.")',finish-start
   print '("Error MAX: ",E11.3)',errormax
!
!   print*,imax,jmax,kmax
!   call hypre_save_solution('test_mg')
#ifdef DEBUG_MG
    write(nf,*) "End compareSolution"
#endif 
   stop 'END TEST'
   !
 end subroutine compareSolution
#endif
end module nossi_mg_hypre

