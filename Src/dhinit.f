        subroutine dhinit(ialr, na, maxnh, maxnd, nspin,lasto, 
     &             listh, listhptr, numh, dS, dSmat, dHmat, dT,
     &             dDscf,dEscf)

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
C REAL*8 DDSCF		      : perturbed density matrix elements
C REAL*8 DESCF		      : perturbed energy density matrix
*******************************************************************

C Modules-----------------------------------------------------------------
!        use m_iodmlr,      only: read_dmlr
        use precision,     only: dp
        use atomlist,      only: no_l, indxuo
        use alloc,         only: re_alloc, de_alloc
        use parallel,      only : Node, Nodes
        use parallelsubs,  only : GlobalToLocalOrb
        use listsc_module, only: listsc
        use alloc,         only : re_alloc
C -------------------------------------------------------------------------

        implicit none

C Internal variable types and dimensions
        integer      :: nspin_read, maxnd_read, j, jo, 
     &                  io, ialr ,na, lasto(0:na),  ind, ix,
     &                  listh(maxnh), listhptr(*), numh(*), juo,
     &                   ind2, k, ko, ispin
        logical      :: found

        integer      :: maxnh, maxnd, nspin, i, ia
        real(dp)     :: DS(maxnh,3), DSMAT(maxnh,3), 
     &                  DHMAT(maxnh,3, nspin), DT(maxnh,3), 
     &                  dDscf(maxnd,nspin,3), dEscf(maxnd,nspin,3)
      
        
      call timer( 'dHinit', 1 )  
C Read change in the density matrix from file ----------------------- 
!     call read_dmlr ( maxnh, no_l, nspin, numh,
!    &                     listhptr, listh, dDscf, found )

C Initialize changed overlap, DM and energy matrices ----------------
      dSmat = 0.0_dp
      if(.not.found) then
        dDscf(1:maxnd,1:nspin,1:3) = 0.0_dp
        dEscf(1:maxnd,1:nspin,1:3) = 0.0_dp
      endif
C Loop on orbitals of atom IALR ----------------------------------------

      do io = lasto(IALR-1)+1, lasto(IALR) !orbitals from moving atom
        do j = 1,numh(io) !interacting orbitals
          ind = listhptr(io)+j !position to io-j interaction in the matrix
          jo = listh(ind)
          juo = indxuo(jo)
          do ix = 1,3
            DSMAT(ind,ix) = DSMAT(ind,ix) - DS(ind,ix)
            do ispin = 1,nspin
              DHMAT(ind,ix,ispin) = DHMAT(ind,ix,ispin)
     .                               - DT(ind,ix)
            enddo
          enddo
C         Idetify neighbour index of orbital io relative to orbital jo -
          do k = 1, numh(juo)
            ind2 = listhptr(juo)+ k
            ko = listsc(jo,juo,listh(ind2))
            if ( ko .EQ. io) then
              do ix = 1,3
                DSMAT(ind2,ix) = DSMAT(ind2,ix) -
     &                           DS(ind,ix)
                do ispin = 1,nspin
                  DHMAT(ind2,ix,ispin) = DHMAT(ind2,ix,ispin)
     &                                 - DT(ind,ix)
                enddo
              enddo
            endif
          enddo
        enddo
      enddo

      call timer( 'dHinit', 2 )

      end subroutine
