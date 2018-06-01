      subroutine ddsmat(dT, dS, d2T, d2S, nspin, no, nuo, na, 
     &                  nua, indxuo, lasto, numh, listh, listhptr, 
     &                  Dscf, Escf, jalr, iaorb, maxnh, dDscf, 
     &                  dEscf, dynmat)

C *************************   LINEAR RESPONSE   ************************
C  Adds contributions to Dynamical Matrix, when atom JALR is moved.  
C  That is, adds the contributions from Kinetic Energy and Overlap and
C  stores these in the corresponding column of DynMat.
C Junquera, Ordejon and Pruneda (Linear response, Jan99)
C Adapted by LR Summer 2015
C    (Energies in Ry.; distances in Bohr)
C DYNMAT terms: terms 7.1, 7.2, 1.2 and 2
C***********************************************************************
      
      use precision,     only: dp
      use listsc_module, only: listsc
      use alloc,         only: re_alloc, de_alloc
      use parallel,      only: Node, Nodes
      use parallelsubs,  only: GlobalToLocalOrb, WhichNodeOrb
#ifdef MPI
      use mpi_siesta
#endif

C Arguments ---------------------------------------------------------
      integer :: na, nua, no, nuo, nspin, maxnh, iaorb(*), 
     &           lasto(0:na), listh(maxnh),listhptr(*),
     &           numh(*), jalr, indxuo(no)

      real(dp) :: dynmat(3,nua,3,nua), dDscf(maxnh,nspin,3),
     &            dEscf(maxnh,nspin,3), dT(maxnh,3), d2T(maxnh,3,3),
     &            dS(maxnh,3), d2S(maxnh,3,3), Escf(maxnh,nspin),
     &            Dscf(maxnh,nspin) 

C Internal variables ------------------------------------------------
      integer :: ind, io, iuo, jo, juo, jn, ko, index(no),
     &           jua, i, ix, jx, ind2, nuotot, BNode, iio,
     &           MPIerror, maxnhg, maxnumh

      integer, pointer :: numhg(:)=>null(),
     &                    listhg(:)=>null(),
     &                    listhptrg(:) => null()

      real(dp), pointer:: d2Haux(:,:,:)=>null(),
     &                    d2Saux(:,:,:)=>null(),
     &                    HSaux(:,:,:)=>null(),
     &                    dDaux(:,:,:)=>null(),
     &                    Daux(:,:)=>null()


      call timer('ddsmat',1)

C ***********************Local to global variables*******************
C Create the GLOBAL listhptr and the numh (since they are local)
      nullify(numhg, listhg, listhptrg)

      nuotot=nuo
C Allocate local memory for global list arrays
      call re_alloc( numhg, 1, nuotot, name='numhg',
     &               routine= 'ddsmat' )
      call re_alloc( listhptrg, 1, nuotot, name='listhptrg',
     &               routine= 'ddsmat' )

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

C Build global listhptr
      listhptrg(1) = 0
      do io = 2,nuotot
        listhptrg(io) = listhptrg(io-1) + numhg(io-1)
      enddo

C Global listh
      maxnhg = listhptrg(nuotot) + numhg(nuotot)
      call re_alloc( listhg, 1, maxnhg, name='listhg',
     &               routine= 'ddsmat' )
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
     &    BNode,MPI_Comm_World,MPIerror)
        endif
#endif
      enddo

C global auxiliar variable: will store dT, d2T, dS and d2S
      maxnumh = 0
      do io = 1,nuotot
        maxnumh = max(maxnumh,numh(io))
      enddo
      nullify( HSaux, d2Haux, d2Saux, dDaux, Daux)
      call re_alloc( HSaux, 1, maxnhg, 1, 3, 1, 2, name='HSaux',
     &               routine= 'ddsmat' ) !1=dT, 2=dS
      call re_alloc( d2Haux, 1, maxnhg, 1, 3, 1, 3, name='d2Haux',
     &               routine= 'ddsmat' ) !d2H
      call re_alloc( d2Saux, 1, maxnhg, 1, 3, 1, 3, name='d2Saux',
     &               routine= 'ddsmat' ) !d2S
      call re_alloc( Daux, 1, maxnhg, 1, 2, name='Daux',
     &               routine= 'ddsmat' )!1=dscf, 2=Escf for each spin
      call re_alloc( dDaux, 1, maxnhg, 1, 3, 1, 2, name='dDaux',
     &               routine= 'ddsmat' )

C Globalize dT and dS into Haux(1,2)
C Globalize d2T and d2S into d2Haux and d2Saux
      HSaux(:,:,:)=0.0_dp
      d2Saux(:,:,:)=0.0_dp
      d2Haux(:,:,:)=0.0_dp

      do io = 1,nuotot
        call WhichNodeOrb(io,Nodes,BNode)
        if (Node.eq.BNode) then
          call GlobalToLocalOrb(io,Node,Nodes,iio)
            do ix = 1, 3
              do jo = 1,numh(iio)
                HSaux(listhptrg(io)+jo,ix,1) =
     &                             DT(listhptr(iio)+jo,ix)
                HSaux(listhptrg(io)+jo,ix,2) =
     &                             DS(listhptr(iio)+jo,ix)
                do jx= 1, 3
                  d2Haux(listhptrg(io)+jo,ix,jx) =
     &                             D2T(listhptr(iio)+jo,ix,jx)
                  d2Saux(listhptrg(io)+jo,ix,jx) =
     &                             D2S(listhptr(iio)+jo,ix,jx)
                enddo
              enddo
            enddo
        endif
#ifdef MPI
        if (Nodes .gt.1 ) then
          do ix = 1, 3
            call MPI_Bcast(HSaux(listhptrg(io)+1,ix,1),numhg(io),
     &      MPI_double_precision,BNode,MPI_Comm_World,MPIerror)
            call MPI_Bcast(HSaux(listhptrg(io)+1,ix,2),numhg(io),
     &      MPI_double_precision,BNode,MPI_Comm_World,MPIerror)
            do jx= 1, 3
              call MPI_Bcast(d2Haux(listhptrg(io)+1,ix,jx),numhg(io),
     &        MPI_double_precision,BNode,MPI_Comm_World,MPIerror)
              call MPI_Bcast(d2Saux(listhptrg(io)+1,ix,jx),numhg(io),
     &        MPI_double_precision,BNode,MPI_Comm_World,MPIerror)
            enddo
          enddo
        endif
#endif
      enddo

      do ispin = 1, nspin
C Globalize Dscf and Escf for each spin component into Daux
C Globalize dDscf and dEscf for each spin component into dDaux
        Daux(:,:)=0.0_dp
        dDaux(:,:,:)=0.0_dp
        do io = 1,nuotot
          call WhichNodeOrb(io,Nodes,BNode)
          if (Node.eq.BNode) then
            call GlobalToLocalOrb(io,Node,Nodes,iio)
            do jo = 1,numh(iio)
              Daux(listhptrg(io)+jo,1) =
     &                           Dscf(listhptr(iio)+jo,ispin)
              Daux(listhptrg(io)+jo,2) =
     &                           Escf(listhptr(iio)+jo,ispin)
              do ix=1,3
                dDaux(listhptrg(io)+jo,ix,1) =
     &                           dDscf(listhptr(iio)+jo,ispin,ix)
                dDaux(listhptrg(io)+jo,ix,2) =
     &                           dEscf(listhptr(iio)+jo,ispin,ix)
              enddo
            enddo
          endif
#ifdef MPI
          if (Nodes .gt.1 ) then
            call MPI_Bcast(Daux(listhptrg(io)+1,1),numhg(io),
     &      MPI_double_precision,BNode,MPI_Comm_World,MPIerror)
            call MPI_Bcast(Daux(listhptrg(io)+1,2),numhg(io),
     &      MPI_double_precision,BNode,MPI_Comm_World,MPIerror)
            do ix=1,3
              call MPI_Bcast(dDaux(listhptrg(io)+1,ix,1),numhg(io),
     &        MPI_double_precision,BNode,MPI_Comm_World,MPIerror)
              call MPI_Bcast(dDaux(listhptrg(io)+1,ix,2),numhg(io),
     &        MPI_double_precision,BNode,MPI_Comm_World,MPIerror)
            enddo
          endif
#endif
        enddo
        do io = 1, nuotot
          iuo = indxuo(io) !EHKEHK
          IALR = iaorb(iuo)

C       Invert the neighb list
          index(:) = 0
          do in = 1,numhg(iuo)
            ind = listhptrg(iuo) + in
            jo = listsc(io,iuo,listhg(ind))
            juo = indxuo(jo)
            do jn = 1,numhg(juo)
              ind2 = listhptrg(juo) + jn
              ko = listsc(jo,juo,listhg(ind2))
              if (ko .eq. io) index(jo) = jn
            enddo              
          enddo
C         ---
C       Loop over neighb orbs
          do j = 1, numhg(iuo)
            ind = listhptrg(iuo) + j
            jo = listhg(ind)
            juo = indxuo(jo)
            jua = iaorb(juo)
C        Identify neighbour indx of orb io relative to orbital jo
            i=index(jo) !EHKEHK
            ind2 = listhptrg(juo) + i !EHKKK22
C        Form changed overlap matrix element --------------------
         !d2S and d2T deviate from Linres1 due to the output of the
         ! interpolation spliu vs splint which is different also 
         ! between the two versions of siesta.
            do ix = 1,3
              do jx = 1,3
                if(IALR.eq.JALR) then
             !terms 7.1 and 7.2
             dynmat(jx,IALR,ix,JALR) = dynmat(jx,IALR,ix,JALR)
     &                            + Daux(ind2,2)*d2Saux(ind2,ix,jx)
     &                            + Daux(ind,2)*d2Saux(ind2,ix,jx)
     &                            - Daux(ind2,1)*d2Haux(ind2,ix,jx)
     &                            - Daux(ind,1)*d2Haux(ind2,ix,jx)
                endif
                if(jua.eq.JALR) then
             ! terms 7.1 and 7.2 
             dynmat(jx,IALR,ix,JALR) = dynmat(jx,IALR,ix,JALR)
     &                             - Daux(ind2,2)*d2Saux(ind2,ix,jx)
     &                             - Daux(ind,2)*d2Saux(ind2,ix,jx)
     &                             + Daux(ind2,1)*d2Haux(ind2,ix,jx)
     &                             + Daux(ind,1)*d2Haux(ind2,ix,jx)
                endif
             ! terms 1.2 and 2
             dynmat(jx,IALR,ix,JALR) = dynmat(jx,IALR,ix,JALR)
     &                         + dDaux(ind2,ix,2)*HSaux(ind2,jx,2)
     &                         + dDaux(ind,ix,2)*HSaux(ind2,jx,2)
     &                         - dDaux(ind,ix,1)*HSaux(ind2,jx,1)
     &                         - dDaux(ind2,ix,1)*HSaux(ind2,jx,1)
              enddo
            enddo
          enddo
        enddo
      enddo
      call de_alloc( Daux, name='Daux', routine= 'ddsmat' )
      call de_alloc( dDaux, name='dDaux', routine= 'ddsmat' )
      call de_alloc( HSaux, name='HSaux', routine= 'ddsmat' )
      call de_alloc( d2Haux, name='d2Haux', routine= 'ddsmat' )
      call de_alloc( d2Saux, name='d2Saux', routine= 'ddsmat' )
      call de_alloc( listhg, name='listhg', routine= 'ddsmat' )
      call de_alloc( listhptrg, name='listhptrg', routine= 'ddsmat' )
      call de_alloc( numhg, name='numhg', routine= 'ddsmat' )


      call timer('ddsmat',2)
C -------------------------------------------------------------------


      end subroutine
