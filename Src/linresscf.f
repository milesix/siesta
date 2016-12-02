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
c      use m_pulay
!      use m_pulaylr,  only: init_pulaylr_arrays, pulayxlr
!      use m_iodmlr,      only: write_dmlr
      use linres_matrices, only: resetFirstmatrices
      use fdf
!      use m_ddnaefs
C----------------------------------------------------------------------------
      implicit none 


C Internal variable types and dimensions -----------------------------------      
      integer           :: ialr, iai, iaf, iscf, iiscf,
     &                     maxno, ispin, io, iuo, in, jo,
     &                     ind, juo, ko, jn, i, j, jua, ix, jx,
     &                     istr, ifa, ilr, idyn, iter, ju, iu, unit1
      
      real(dp)          :: dSmat(maxnh,3),!First Order of Overlap 
     &                     dHmat(maxnh,3,nspin),!Perturbed Hamiltonian
     &                     dDmax,Dmax, 
     &                     dHmat0(maxnh,3,nspin),!non SCF perturbed H 
     &                     dummy_stress(3,3), dummy_fa(1,1), g2max 
      real(dp), pointer :: dDold(:,:,:)
      real(dp)          :: dynmat(na_u,3,na_u,3)

      real(dp)          ::  tollr, eigtollr !!!!!provisional

      character(len=label_length+3) :: fname

      logical            :: FIRST, dummy_use_rhog_in, 
     &                      dummy_chargedensonly,mmix,check,printdyn!!!!!!

!      external dhinit, delrho, ddsmat, io_assign, io_close 
C ----------------------------------------------------------------------------


CCCC  PROVISIONAL SOLUTION THE READ LINRES OPTIONS
      iai = fdf_get("LR.IAI",1)
      iaf = fdf_get("LR.IAF",1)
      tollr = fdf_get("LR.DMTolerance",0.001_dp)
      eigtollr = fdf_physical("LR.EigTolerance",0.001_dp,'Ry')

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCc

C Initialize dDscf, dEscf, dynmat----------------------------------------------
      nullify(dDscf, dEscf)
      call re_alloc(dDscf, 1, maxnh, 1, nspin, 1, 3,
     &                 'dDscf', 'linresscf')
      call re_alloc(dEscf, 1, maxnh, 1, nspin, 1, 3,
     &                 'dEscf', 'linresscf')

      dynmat(:,:,:,:)=0.0_dp
C ----------------------------------------------------------------------------

C Linres flags
      ilr=1  ! Compute perurbed hamiltonian elemts (1 yes, 0 no)
      idyn=0 ! Compute dynamical terms (1 yes, 0 no)

C IS A K CALCULATION?????








C Dummies initialization------------------------------------------------------
      dummy_chargedensonly = .false.
      dummy_use_rhog_in = .false.

      call timer('Linres atoms loop',1)
C-----------------------------------------------------------------------------

C------------------------------LINRES MAIN LOOP-------------------------------
C Init of Linres calculation. External loop over perturbed atoms (IALR).
C Internal loop scf-loop (ISCF).

      do 200 ialr=iai, iaf 
        first = .true.
        
C init pulay arrays !!!!!!!!!!!!!!!!!!!!!!

        call timer('Linres SCF loop', 1)

C       scf loop -----------------------------------------------
        do 100 iscf=1, nscf

          dHmat(1:maxnh,1:3,1:nspin) = 0.0_dp

          if (first) then
C Calculation of the kinetic terms of the perturbed hamiltonian
            call dhinit(IALR, na_s, maxnh, maxnh, nspin, lasto,
     &       listh, listhptr, numh, DS, dSmat, dHmat,dH, dDscf,
     &       dEscf)

C Calculation of the KB terms of the perturbed hamiltonian
C (added into dHmat)
            call dvnloc(scell, na_u, na_s, isa, xa, indxua, Dscf,
     &                 maxnh, maxnh, lasto, lastkb, iphorb,
     &                 iphKB, numh, listhptr, listh, numh,
     &                 listhptr, listh, min(nspin,2), IALR,
     &                 no_s, iaorb, dDscf, dHmat, dynmat, first)
          endif

          if (first) then
C Up to here, perturbed hamiltonian is Kinetic+KB. These terms are non-scf and 
C are stored in dHmat0 to avoid calculating them in each iscf.
            do ispin=1,nspin
             dHmat0(1:maxnh,1:3,ispin) = dHmat(1:maxnh,1:3,ispin)
            enddo
             dHmat(:,:,:)=0.0_dp
          endif

          if (first) then 
C Non-scf elements of dynamat are calculated and added into dynmat. These elements
C depend on gradients of densities and gradients of orbitals

          ! call dhscf

          endif 

          first = .false.
C Deactivate first calls flag. From now, perturbed hamiltonian elements are going
C to be calculated. First of all, dHmat0 is updated including the pulay terms.
C dHmat is calculated as dHmat=HperSCF+dHmat0 

          ! call dhscf


          nullify(dDold)
          call re_alloc(dDold,1,maxnh,1,nspin,1,3,
     &                 'dDold', 'linresscf')

C Copy current density as old one to perform the mixing
          dDold(1:maxnh,1:nspin,1:3) = dDscf(1:maxnh,1:nspin,1:3)
          dDscf(1:maxnh,1:nspin,1:3) = 0.0_dp
          dEscf(1:maxnh,1:nspin,1:3) = 0.0_dp

C Find change in density matrix from perturbed 
C Hamiltonian and Overlap

          ! call delrho

C Perform the density matrix mixing
          dmax=0.0_dp
          do ix=1,3
!            call pulay
             dmax=max(dmax.dDmax)
          enddo


C Check convergence cryteria
          if (dmax .lt. tollr) goto 999

 100    enddo  ! isc loop

 999  continue 
C End scf loop

 200  enddo   ! ialr atoms loop



      do i=1,maxnh
        print*,'i,1,2,3',i,dhmat(i,1,1),dhmat(i,2,1),dhmat(i,3,1)
      enddo
        stop




      end subroutine 

      end MODULE m_linresscf

