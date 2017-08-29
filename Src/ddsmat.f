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
C Arguments ---------------------------------------------------------
      integer :: na, nua, no, nuo, nspin, maxnh, iaorb(*), 
     &           lasto(0:na), listh(maxnh),listhptr(*),
     &           numh(*), jalr, indxuo(no)

      real(dp) :: dynmat(nua,3,nua,3), dDscf(maxnh,nspin,3),
     &            dEscf(maxnh,nspin,3), dT(maxnh,3), d2T(maxnh,3,3),
     &            dS(maxnh,3), d2S(maxnh,3,3), Escf(maxnh,nspin),
     &            Dscf(maxnh,nspin) 

C Internal variables ------------------------------------------------
      integer :: ind, io, iuo, jo, juo, jn, ko, index(no),
     &           jua, i, ix, jx, ind2

      call timer('ddsmat',1)

      do ispin = 1, nspin
       do io = 1, nuo
        iuo = indxuo(io) !EHKEHK
        IALR = iaorb(iuo)

C       Invert the neighb list
        index(1:nuo) = 0
        do in = 1,numh(iuo)
         ind = listhptr(iuo) + in
         jo = listsc(io,iuo,listh(ind))
         juo = indxuo(jo)
         do jn = 1,numh(juo)
          ind2 = listhptr(juo) + jn
          ko = listsc(jo,juo,listh(ind2))
          if (ko .eq. io) index(jo) = jn
         enddo              
        enddo
C         ---

C       Loop over neighb orbs
        do j = 1, numh(iuo)
         ind = listhptr(iuo) + j
         jo = listh(ind)
         juo = indxuo(jo)
         jua = iaorb(juo)
C        Identify neighbour indx of orb io relative to orbital jo
         i=index(jo) !EHKEHK
         ind2 = listhptr(juo) + i !EHKKK22
C        Form changed overlap matrix element --------------------
         !d2S and d2T deviate from Linres1 due to the output of the
         ! interpolation spliu vs splint which is different also 
         ! between the two versions of siesta.
         do ix = 1,3
          do jx = 1,3
           if(IALR.eq.JALR) then
             !terms 7.1 and 7.2
             dynmat(IALR,jx,JALR,ix) = dynmat(IALR,jx,JALR,ix)
     &                            + escf(ind2,ispin)*d2S(ind2,ix,jx)
     &                            + escf(ind,ispin)*d2S(ind2,ix,jx)
     &                            - dscf(ind2,ispin)*d2t(ind2,ix,jx)
     &                            - dscf(ind,ispin)*d2t(ind2,ix,jx)
           endif
           if(jua.eq.JALR) then
             ! terms 7.1 and 7.2 
             dynmat(IALR,jx,JALR,ix) = dynmat(IALR,jx,JALR,ix)
     &                             - escf(ind2,ispin)*d2S(ind2,ix,jx)
     &                             - escf(ind,ispin)*d2S(ind2,ix,jx)
     &                             + dscf(ind2,ispin)*d2t(ind2,ix,jx)
     &                             + dscf(ind,ispin)*d2t(ind2,ix,jx)
           endif
             ! terms 1.2 and 2
             dynmat(IALR,jx,JALR,ix) = dynmat(IALR,jx,JALR,ix)
     &                             + descf(ind2,ispin,ix)*ds(ind2,jx)
     &                             + descf(ind,ispin,ix)*ds(ind2,jx)
     &                             - ddscf(ind,ispin,ix)*dt(ind2,jx)
     &                             - ddscf(ind2,ispin,ix)*dt(ind2,jx)

           
          enddo
         enddo
        enddo

       enddo
      enddo
      call timer('ddsmat',2)
C -------------------------------------------------------------------


      end subroutine
