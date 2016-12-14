! ---
! Copyright (C) 1996-2016       The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---

      MODULE m_linresscf
      private
      public :: linresscf

      CONTAINS 

      subroutine linresscf ()
C*************************************************************************
C     Riches 2015: this is the main routine that handles LinRes.      
C     Is is called from siesta_forces once the selfconsistency 
C     for the ground state has been reached.
C
C
C
C
C
C
C
C
C
C

C***************************************************************************


C Modules------------------------------------------------------------------
      use atomlist,      only : lasto, no_s, indxuo, iaorb, indxua,
     &                          iphkb, iphorb, lastkb, datm, no_l, 
     &                          no_u, rmaxv 
      use precision,  only: dp, sp
      use siesta_geom
      use siesta_options,  only: nscf, temp, maxsav, wmix,nkick,
     &                            wmixkick
      use parallel,          only: IOnode
      use m_spin,        only: nspin
      use sparse_matrices
      use m_overfsm,     only: overfsm 
      use listsc_module, only: listsc
      use m_dvnloc,      only: dvnloc
      use m_dhscf,       only: dhscf
      use m_dipol,       only: dipol
      use m_energies
      use files,         only:  filesOut_t,slabel, label_length
      use units,         only: Ang, eV 
      use m_ntm
      use alloc
      use m_gamma
      use m_eo
      use Kpoint_grid,   only: nkpnt, kpoint, kweight
      use m_energies,    only: ef
      use m_pulay,       only: init_pulay_arrays, pulayx
!      use m_iodmlr,      only: write_dmlr
      use linres_matrices, only: resetFirstmatrices
      use fdf
      use m_ddnaefs
      use m_iodynmat
C----------------------------------------------------------------------------
      implicit none 


C Internal variable types and dimensions -----------------------------------      
      integer           :: ialr, iai, iaf, iscf, iiscf,
     &                     maxno, ispin, io, iuo, in, jo,
     &                     ind, juo, ko, jn, i, j, jua, ix, jx,
     &                     istr, ifa, ilr, idyn, iter, ju, iu, unit1,
     &                     init
      
      real(dp)          :: dSmat(maxnh,3),!First Order of Overlap 
     &                     dHmat(maxnh,3,nspin),!Perturbed Hamiltonian
     &                     dDmax,Dmax, 
     &                     dHmat0(maxnh,3,nspin),!non SCF perturbed H 
     &                     dummy_stress(3,3), dummy_fa(1,1), g2max 
      real(dp), pointer :: dDold(:,:,:)
      real(dp)          :: dynmat(na_u,3,na_u,3)

      real(dp)          ::  tollr, eigtollr 

      type(filesOut_t)    :: filesOut
      character(len=label_length+3) :: fname

      logical            :: dummy_use_rhog_in, LRfirst,
     &                      dummy_chargedensonly,mmix,
     &                      readold

C ----------------------------------------------------------------------------

C ----------------------------------------------------------------------------
C Reading linres options 
C ----------------------------------------------------------------------------
      iai = fdf_get("LR.IAI",1)
      iaf = fdf_get("LR.IAF",na_u) ! Move all atoms by default
      tollr = fdf_get("LR.DMTolerance",0.001_dp)
      eigtollr = fdf_physical("LR.EigTolerance",0.001_dp,'Ry')
      readold = fdf_get("LR.readDynmat",.true.)


C Begin to write into output file
      if (IOnode) then
        write(6,'(/,t22,a)') repeat('=',36)
        write(6,'(t32,a,i7)') 'Linres calculation'
        write(6,'(t22,a)') repeat('=',36)

        write(6,'(/,a,i7)')
     .    'Some linres parameters:'
        write(6,'(a,i7)') 'Linres: Initial perturbed atom', iai
        write(6,'(a,i7)') 'Linres: Final perturbed atom', iaf
      endif

C Initialize dDscf, dEscf, dynmat----------------------------------------------
      nullify(dDscf, dEscf)
      call re_alloc(dDscf, 1, maxnh, 1, nspin, 1, 3,
     &                 'dDscf', 'linresscf')
      call re_alloc(dEscf, 1, maxnh, 1, nspin, 1, 3,
     &                 'dEscf', 'linresscf')

      dynmat(:,:,:,:)=0.0_dp

C IS A K CALCULATION????? Could be better....
      if (nkpnt .eq. 1) then
        GAMMA = .TRUE.
      else
        GAMMA = .FALSE.
      endif

C Dummies initialization------------------------------------------------------
      dummy_chargedensonly = .false.
      dummy_use_rhog_in = .false.

      call timer('Linres',1)
C-----------------------------------------------------------------------------

C Read Stored dynamical matrix files-----------------------------------------
      init=iai
      if (readold) then
        call readdynmat(init,iai,dynmat)
      else
        init=iai
      endif
C------------------------------LINRES MAIN LOOP-------------------------------
C Init of Linres calculation. External loop over perturbed atoms (IALR).
C Internal loop scf-loop (ISCF).

      do 200 ialr=init, iaf 

      if (IOnode) then
        write(6,'(/,t22,a)') repeat('=',36)
        write(6,'(t32,a,i7)') 'Linres atom=', ialr
        write(6,'(t22,a)') repeat('=',36)
      endif

      LRfirst=.true. !init flag for first calls


C Doing non-scf elements of the perturbed Hamiltonian and dinamical matrix 
        dHmat0(:,:,:)=0.0_dp
        call dhscf( nspin, no_s, iaorb, iphorb, no_l,
     .            no_u, na_u, na_s, isa, xa, indxua,
     .            ntm, 0, 0, 0, filesOut,
     .            maxnh, numh, listhptr, listh, Dscf, Datm,
     .            maxnh, H, Enaatm, Enascf, Uatm, Uscf, DUscf, DUext,
     .            Exc, Dxc, dipol, dummy_stress, dummy_fa,
     .            dummy_stress, dummy_use_rhog_in,
     .            dummy_chargedensonly, iai, iaf, ialr, lasto,
     .            dynmat, dDscf,dHmat, LRfirst, dHmat0)

C dHmat0 contains the pulay terms of the perturbed hamiltonian
C Dynamat contains all the elements that are non-scf and depend
C on the gradient potential terms

        call init_pulay_arrays(iai)        

        call timer('LRatom', 1)

C       scf loop -----------------------------------------------
        do 100 iscf=1, nscf

          dHmat(1:maxnh,1:3,1:nspin) = 0.0_dp

          if (iscf.eq.1) then
C Calculation of the kinetic terms of the perturbed hamiltonian (added into dHmat)
            call dhinit(IALR, na_s, maxnh, maxnh, nspin, lasto,
     &       listh, listhptr, numh, DS, dSmat, dHmat,dH, dDscf,
     &       dEscf)

C Calculation of the KB terms of the perturbed hamiltonian
C (added into dHmat)
            call dvnloc(scell, na_u, na_s, isa, xa, indxua, Dscf,
     &                 maxnh, maxnh, lasto, lastkb, iphorb,
     &                 iphKB, numh, listhptr, listh, numh,
     &                 listhptr, listh, min(nspin,2), IALR,
     &                 no_s, iaorb, dDscf, dHmat, dynmat, LRfirst)

C Update the non-scf perturbed hamiltonian with the kinetic and KB terms
C stored in dHmat. dHmat0 contains the pulay terms. After this point, dHmat0 
C contains all the non-scf elements and it is constant. It is added every scf 
C step to the perturbed hamiltonian!!
            do ispin=1,nspin
              dHmat0(1:maxnh,1:3,ispin)=dHmat0(1:maxnh,1:3,ispin) + 
     .                                   dHmat(1:maxnh,1:3,ispin)
            enddo
            dHmat(1:maxnh,1:3,1:nspin) = 0.0_dp
          endif
         
 
C Copy non-scf hamiltonian to total Hamiltonian
        do ispin=1,nspin
        dHmat(1:maxnh,1:3,ispin)=dHmat0(1:maxnh,1:3,ispin)
        enddo

C Compute the perturbed potential elements Vxc, Vna and Vh and 
C include them into dHmat. Look that dhscf is called without dHmat0!!!

        call dhscf( nspin, no_s, iaorb, iphorb, no_l,
     .            no_u, na_u, na_s, isa, xa, indxua,
     .            ntm, 0, 0, 0, filesOut,
     .            maxnh, numh, listhptr, listh, Dscf, Datm,
     .            maxnh, H, Enaatm, Enascf, Uatm, Uscf, DUscf, DUext,
     .            Exc, Dxc, dipol, dummy_stress, dummy_fa,
     .            dummy_stress, dummy_use_rhog_in,
     .            dummy_chargedensonly, iai, iaf, ialr, lasto,
     .            dynmat, dDscf, dHmat, LRfirst)

          nullify(dDold)
          call re_alloc(dDold,1,maxnh,1,nspin,1,3,
     &                 'dDold', 'linresscf')

C Copy current density as old one to perform the mixing
          dDold(1:maxnh,1:nspin,1:3) = dDscf(1:maxnh,1:nspin,1:3)
          dDscf(1:maxnh,1:nspin,1:3) = 0.0_dp
          dEscf(1:maxnh,1:nspin,1:3) = 0.0_dp

C Find change in density matrix from perturbed 
C Hamiltonian and Overlap

          call delrho(no_s, na_s, nspin, maxnh, no_l, no_u, GAMMA,
     &               indxuo, xijo, nspin, eo, tollr, eigtollr, nkpnt,
     &               kpoint, kweight, Qo, H, S, dHmat, dSmat, numh,
     &               listh, listhptr, ef, temp, dDscf, dEscf,iscf)

C Perform the density matrix mixing
          dmax=0.0_dp
          do ix=1,3
            call pulayx( iscf , mmix , no_l, maxnh, numh,
     &                  listhptr, nspin,maxsav,wmix,nkick,
     &                  wmixkick,dDscf(:,:,ix),dDold(:,:,ix),dDmax,ix)
            dmax=max(dmax,dDmax)
          enddo

C Print error in the perturbed density
          if (iscf.eq.1) then
            print*, 'siesta Linres: Density tolerance dDTol=', TOLLR
            print*, 'LINRES:',iscf,'dMax=',dMax
          else
            print*, 'LINRES:',iscf,'dMax=',dMax
          endif

          LRfirst=.false.

C Check convergence cryteria
          if (dmax .lt. tollr) goto 999

 100    enddo  ! isc loop
C------------------------------------------------------------------------

 999    continue ! Exit scf loop
        call timer('LRatom', 2)

C Now, is time to calculate the final terms of the dynamical matrix
C------------------------------------------------------------------------

C     Add non-local potential and kinetic part contributions to dynmat---
C     (the terms that depend on the perturbed density and the ones that 
C      do not depend)
C     Note that LRfirst now is false
        call dvnloc(scell, na_u, na_s, isa, xa, indxua, Dscf,
     &                 maxnh, maxnh, lasto, lastkb, iphorb,
     &                 iphKB, numh, listhptr, listh, numh,
     &                 listhptr, listh, min(nspin,2), IALR,
     &                 no_s, iaorb, dDscf, dHmat, dynmat, LRfirst)

C     Add kinetic and ovelap contributions to dynmat -----------------
        call ddsmat(dH, dS, d2H, d2S, nspin, no_s, no_u, na_s, na_u,
     &                indxuo, lasto, numh, listh, listhptr, Dscf, Escf,
     &                ialr, iaorb, maxnh, dDscf, dEscf, dynmat)


c     Add terms that depends on the perturbed density and Vscf 
C     and perturbed Vscf potentials-----------------------------------

        call dhscf( nspin, no_s, iaorb, iphorb, no_l,
     .            no_u, na_u, na_s, isa, xa, indxua,
     .            ntm, 0, 0, 0, filesOut,
     .            maxnh, numh, listhptr, listh, Dscf, Datm,
     .            maxnh, H, Enaatm, Enascf, Uatm, Uscf, DUscf, DUext,
     .            Exc, Dxc, dipol, dummy_stress, dummy_fa,
     .            dummy_stress, dummy_use_rhog_in,
     .            dummy_chargedensonly, iai, iaf, ialr, lasto,
     .            dynmat, dDscf, dHmat)


C Dealloc the non-scf variables (drhoatm, drhoscf0 and dvxc) for next atom
        call resetFirstmatrices()

C Save in a file the dynamical matrix and extra flag of the index ialr of the last atom
        call writedynmat(iai,iaf,ialr,dynmat,.false.)

 200  enddo   ! ialr atoms loop
      call timer('Linres',2)
C-------------------------------------------------------------------------

C There is other contribution to the dynamical matrix (non-scf)
C the Laplacian of the neutral atom potential
      call ddnaefs(na_u, na_s, scell, xa, indxua, rmaxv,
     &            isa, dynmat)

C-------------------------------------------------------------------------
C-----------------Finally, print the dynamical matrix to FC file----------
C-------------------------------------------------------------------------
      call writedynmat(iai,iaf,ialr,dynmat,.true.)

C-------------------------------------------------------------------------
      call de_alloc(dDscf, 'dDscf', 'linresscf')
      call de_alloc(dEscf,'linresscf')

        end subroutine 

      end MODULE m_linresscf

