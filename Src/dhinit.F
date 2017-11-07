        subroutine dhinit(ialr, na, nuotot, maxnh, maxnd, nspin,lasto, 
     &             listh, listhptr, numh, dS, dSmat, dHmat, dT)

C *************************   LINEAR RESPONSE   ************************
C Initializes first order changed overlap matrix elements and
C constructs first order change in non-self-consistent hamiltonian 
C coming from kinetic part.
C Junquera, Ordejon and Pruneda 
C SIESTA implementation L. Riches, March '15
C Cheked and commented S. Illera, April '16
***************************INPUT*************************************
C INTEGER IALR		      : index of displaced atom
C INTEGER NA		      : number of atoms in the supercell
C INTEGER MAXNH		      : Maximum number of orbitals interacting
C INTEGER MAXND		      : Maximum number of nonzero elements of 
C                               each row of density matrix
C INTEGER NSPIN		      : number of different spin polarizations
C INTEGER LASTO		      : position of last orbital of each atom
C INTEGER LISTH(maxnh)        : Nonzero hamiltonian-matrix element  
C                               column indexes for each matrix row
C INTEGER LISTHPTR(nuo)       : Pointer to each row (-1) of the
C                               hamiltonian matrix
C INTEGER NUMH		      : Number of nonzero elements of each row 
C                               of density matrix
C REAL*8 DS		      : change in the overlap matrix elements as 
C				given by SIESTA
C REAL*8 DT		      : change in the kinetic elements 
****************************OUTPUT**********************************
C REAL*8 DSMAT		      : overlap changed matrix elements for moving IALR
C REAL*8 DHMAT		      : perturbed hamiltonian 
*******************************************************************

C Modules-----------------------------------------------------------------
        use precision,     only: dp
        use atomlist,      only: no_l, indxuo, indxua, iaorb
        use alloc,         only: re_alloc, de_alloc
        use parallel,      only: Node, Nodes
        use parallelsubs,  only: GlobalToLocalOrb, WhichNodeOrb
        use listsc_module, only: listsc
        use alloc,         only: re_alloc
        use mesh,          only: iatfold, nmsc, cmesh
        use siesta_geom,   only: xa
#ifdef MPI
        use mpi_siesta
#endif

C -------------------------------------------------------------------------

        implicit none

        integer  :: ialr, na, nuotot, maxnh, maxnd, nspin, lasto(0:na),
     &              listh(maxnh), listhptr(*), numh(*)
        real(dp) :: DS(maxnh,3), DSMAT(maxnh,3),
     &              DHMAT(maxnh,3, nspin), DT(maxnh,3)

C Internal variable types and dimensions
        integer      ::  j, jo, io, iio, ind, ix,
     &                  juo, ind2, k, ko, ispin, i, ia, iua, iu, ju
        integer      :: BNode, MPIerror, maxnhg, maxnumh
        integer, pointer :: numhg(:)=>null(),
     &                      listhptrg(:) => null(),
     &                      listhg(:) => null()
        real(dp),pointer :: DSg(:,:)=>null(),
     &                      DSMATg(:,:)=>null(),
     &                      DHMATg(:,:,:)=>null(),
     &                      DTg(:,:)=>null(),
     &                      DHMAT0(:,:,:)=>null(),
     &                      Daux(:)=>null()

        real(dp)     :: displaat(3), dist(3), qxij, qpoint(3), cqxij, pi
        
      call timer( 'dhinit', 1 )  
C For future implementation ----------------------- 

      qpoint=(/0.0,0.0,0.0/)
      pi = 4*atan(1.0_dp)

C Set to zero output overlap matrix
      dSmat(:,:) = 0.0_dp

C ********************Local to global variables********************
C Create the GLOBAL listhptr and the numh (since they are local)
      nullify(numhg, listhptrg)
C Allocate local memory for global list arrays
      call re_alloc( numhg, 1, nuotot, name='numhg',
     &               routine= 'dhinit' )
      call re_alloc( listhptrg, 1, nuotot, name='listhptrg',
     &               routine= 'dhinit' )

C Build global numh
      do io = 1,nuotot
        call WhichNodeOrb(io,Nodes,BNode)
        if (Node.eq.BNode) then
          call GlobalToLocalOrb(io,Node,Nodes,iio)
          numhg(io) = numh(iio)
        endif
#ifdef MPI
        if (Nodes .gt. 1) then
          call MPI_Bcast(numhg(io),1,MPI_integer,BNode,
     &    MPI_Comm_World,MPIerror)
        endif
#endif
      enddo

#ifdef MPI
      if (Nodes .gt. 1) then
        nullify(DHMAT0)
        call re_alloc(DHMAT0,1,maxnh,1,3,1,nspin,
     &               name='DHMAT0',routine='dhinit')
        DHMAT0(:,:,:)=0.0_dp

        nullify(Daux)
        call re_alloc(Daux, 1, maxnumh, name='Daux',
     &               routine= 'dhinit' )
        Daux(:)=0.0_dp
      endif
#endif

C Build global listhptr
      listhptrg(1) = 0
      do io = 2,nuotot
        listhptrg(io) = listhptrg(io-1) + numhg(io-1)
      enddo

C Global listh
      maxnhg = listhptrg(nuotot) + numhg(nuotot)
      nullify(listhg)
      call re_alloc( listhg, 1, maxnhg, name='listhg',
     &               routine= 'dhinit' )

      do io = 1,nuotot
        call WhichNodeOrb(io,Nodes,BNode)
        if (Node.eq.BNode) then
          call GlobalToLocalOrb(io,Node,Nodes,iio)
          do jo = 1,numhg(io)
            listhg(listhptrg(io)+1:listhptrg(io)+numhg(io)) =
     &        listh(listhptr(iio)+1:listhptr(iio)+numh(iio))
          enddo
        endif
#ifdef MPI
        if (Nodes .gt. 1) then
          call MPI_Bcast(listhg(listhptrg(io)+1),numhg(io),MPI_integer,
     &       BNode,MPI_Comm_World,MPIerror)
        endif
#endif
      enddo

C******************local to global variables******************

#ifdef MPI
      maxnumh = 0
      do io = 1,nuotot
        maxnumh = max(maxnumh,numhg(io))
      enddo
      nullify(Daux)
      call re_alloc( Daux, 1, maxnumh, name='Daux',
     &               routine= 'dhinit' )

#endif

C****************Globalize H and S****************************
C dS and DT (inputs) are local, copy to global variable, operate, and
C for output, change to local again
      nullify( DSg, DTg)
      call re_alloc( DSg, 1, maxnhg,1,3, name='DSg',
     &               routine= 'dhinit' )
      call re_alloc( DTg, 1, maxnhg, 1, 3, name='DTg',
     &               routine= 'dhinit' )

      do io = 1,nuotot
        call WhichNodeOrb(io,Nodes,BNode)
        if (Node.eq.BNode) then
          call GlobalToLocalOrb(io,Node,Nodes,iio)
            do ix = 1, 3
              do jo = 1,numh(iio)
                dTg(listhptrg(io)+jo,ix) =
     &                             DT(listhptr(iio)+jo,ix)
                dSg(listhptrg(io)+jo,ix) =
     &                             DS(listhptr(iio)+jo,ix)
              enddo
            enddo
        endif
#ifdef MPI
        if (Nodes .gt.1 ) then
          do ix = 1, 3
            call MPI_Bcast(dTg(listhptrg(io)+1,ix),numhg(io),
     &        MPI_double_precision,BNode,MPI_Comm_World,MPIerror)
            call MPI_Bcast(dSg(listhptrg(io)+1,ix),numhg(io),
     &        MPI_double_precision,BNode,MPI_Comm_World,MPIerror)
          enddo
        endif
#endif
      enddo

! global copies of the distributed output
      nullify( DSMATg, DHMATg)
      call re_alloc( DSMATg, 1, maxnhg, 1, 3, name='DSMATg',
     &               routine= 'dhinit' )
      call re_alloc( DHMATg,1,maxnhg,1,3,1,nspin,name='DHMATg',
     &               routine= 'dhinit' )

      DSMATg(:,:)=0.0_dp
      DHMATg(:,:,:)=0.0_dp

C*************************************************************

C Loop on orbitals of atom IALR ----------------------------------------
      do io = lasto(IALR-1)+1, lasto(IALR) !orbitals from moving atom
!        call GlobalToLocalOrb(io,Node,Nodes,iio)
        call WhichNodeOrb(io,Nodes,BNode)
        iio=io
        if (Node .eq. BNode) then !local orbital
          do j = 1,numhg(iio) !interacting orbitals
            ind = listhptrg(iio)+j !position to io-j interaction
            jo = listhg(ind)
            juo = indxuo(jo)
            do ix = 1,3 ! calculate mu,nu elements (mu is in the unit cell)
              DSMATg(ind,ix) = DSMATg(ind,ix) - DSg(ind,ix)
              do ispin = 1,nspin
                DHMATg(ind,ix,ispin) = DHMATg(ind,ix,ispin)
     .                               - DTg(ind,ix)
              enddo
            enddo
C         Idetify neighbour index of orbital io relative to orbital jo -
            do k = 1, numhg(juo)
              ind2 = listhptrg(juo)+ k
              ko = listsc(jo,juo,listhg(ind2))
              if ( ko .EQ. iio) then !add contributions nu,mu elements
                ia  = iaorb(listhg(ind2))
                iua = indxua(ia)
                do ix = 1, 3
                  displaat(ix) =
     &                 (iatfold(1,ia)*nmsc(1))*cmesh(ix,1)+
     &                 (iatfold(2,ia)*nmsc(2))*cmesh(ix,2)+
     &                 (iatfold(3,ia)*nmsc(3))*cmesh(ix,3)
                end do
                dist(:) = xa(:,iua) - xa(:,ia) - displaat(:)
                qxij    = qpoint(1) *  dist(1)  +
     &                    qpoint(2) *  dist(2)  +
     &                    qpoint(3) *  dist(3)
                cqxij = cos(qxij)

                do ix = 1,3
                  DSMATg(ind2,ix) = DSMATg(ind2,ix) -
     &                           DSg(ind,ix)
                  DSMATg(ind2,ix) = cqxij * DSMATg(ind2,ix)
                  do ispin = 1,nspin
                    DHMATg(ind2,ix,ispin) = DHMATg(ind2,ix,ispin)
     &                                 - DTg(ind,ix)
                  enddo
                enddo
              endif
            enddo
          enddo
        endif
      enddo

! Transform from global to local the Hamiltonian and the overlap
#ifdef MPI
      if (Nodes .gt. 1) then
        do io = 1, nuotot
          call WhichNodeOrb(io,Nodes,BNode)
          call GlobalToLocalOrb(io,BNode,Nodes,iio)
          do ix = 1,3
            do ispin = 1,nspin
              call MPI_Reduce(DHMATg(listhptrg(io)+1,ix,ispin),
     &          Daux(1),numhg(io),MPI_double_precision,
     &          MPI_sum,BNode,MPI_Comm_World,MPIerror)
              if (Node.eq.BNode) then
                do j = 1,numh(iio)
                  DHMAT0(listhptr(iio)+j,ix,ispin) = Daux(j)
                enddo
              endif
            enddo
            Daux(:)=0.0_dp
            call MPI_Reduce(DSMATg(listhptrg(io)+1,ix),
     &         Daux(1),numhg(io),MPI_double_precision,
     &         MPI_sum,BNode,MPI_Comm_World,MPIerror)
            if (Node.eq.BNode) then
              do j = 1,numh(iio)
                DSMAT(listhptr(iio)+j,ix) = Daux(j)
              enddo
            endif
          enddo
        enddo

! add to the input the here calculated contribution (DHMAT0)
        DHMAT(:,:,:)=DHMAT(:,:,:)+DHMAT0(:,:,:)
      else
! add to the input the here calculated contribution (DHMATg)
        DSMAT(:,:)=DSMATg(:,:)
        DHMAT(:,:,:)=DHMATg(:,:,:)+DHMAT(:,:,:) 
      endif
#else
! add to the input the here calculated contribution (DHMATg)
      DSMAT(:,:)=DSMATg(:,:)
      DHMAT(:,:,:)=DHMATg(:,:,:)+DHMAT(:,:,:) 
#endif

      call de_alloc( DSg, name='DSg', routine= 'dhinit' )
      call de_alloc( DTg, name='DTg', routine= 'dhinit' )
      call de_alloc( DHMATg, name='DSMATg', routine= 'dhinit' )
      call de_alloc( DSMATg, name='DHMATg', routine= 'dhinit' )
      call de_alloc( listhg, name='listhg', routine= 'dhinit' )
      call de_alloc( listhptrg, name='listhptrg', routine= 'dhinit' )
      call de_alloc( numhg, name='numhg', routine= 'dhinit' )

#ifdef MPI
      if (Nodes .gt. 1) then
        call de_alloc(DHMAT0,name='DHMAT0',routine='dhinit')
        call de_alloc(Daux,name='Daux',routine='dhinit')
      endif 
#endif

      call timer( 'dhinit', 2 )

      end subroutine

