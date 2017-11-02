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

C***************************************************************************


C Modules------------------------------------------------------------------
      use atomlist,      only : lasto, no_s, indxuo, iaorb, indxua,
     &                          iphkb, iphorb, lastkb, datm, no_l, 
     &                          no_u, rmaxv 
      use precision,  only: dp, sp
      use siesta_geom
      use siesta_options,  only: nscf, temp, maxsav, wmix,nkick,
     &                            wmixkick
      use parallel,          only: IOnode, Node
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
      use m_pulay,       only: init_pulay_arrays, pulayx,
     &                         resetPulayArrays
!      use m_iodmlr,      only: write_dmlr
      use linres_matrices, only: resetFirstmatrices
      use fdf
      use m_ddnaefs
      use m_iodynmat
#ifdef MPI
      use m_mpi_utils, only: globalize_max
#endif

C----------------------------------------------------------------------------
      implicit none 


C Internal variable types and dimensions -----------------------------------      
      integer           :: ialr, iai, iaf, iscf,
     &                     maxno, ispin, io, iuo, in, jo,
     &                     ind, juo, ko, jn, i, j, jua, ix, jx,
     &                     ju, iu, unit1,
     &                     init
      
      real(dp)          :: dSmat(maxnh,3),!First Order of Overlap 
     &                     dHmat(maxnh,3,nspin),!Perturbed Hamiltonian
     &                     dDmax,dMax, 
     &                     dHmat0(maxnh,3,nspin),!non SCF perturbed H 
     &                     dummy_stress(3,3), dummy_fa(3,na_u), g2max 
      real(dp), pointer :: dDold(:,:,:)
      real(dp)          :: dynmat(3,na_u,3,na_u)

      real(dp)          ::  tolLR, eigtolLR 

      type(filesOut_t)    :: filesOut

      character(len=label_length+5+5) :: fname
      character(len=6) :: atomdisp

      logical            :: dummy_use_rhog_in, first_LR,
     &                      dummy_chargedensonly,mmix,
     &                      readold, converged, found
#ifdef MPI
      real(dp) :: buffer1
#endif


C ----------------------------------------------------------------------------

C ----------------------------------------------------------------------------
C Reading linres options 
C ----------------------------------------------------------------------------
      iai = fdf_get("MD.FCFirst",1)
      iaf = fdf_get("MD.FCLast",na_u) ! Move all atoms by default
      tolLR = fdf_get("LR.DMTolerance",0.001_dp)
      eigtolLR = fdf_get('LR.EigTolerance',0.001_dp)
      readold = fdf_get("LR.readDynmat",.false.)

C Begin to write into output file
      if (IOnode) then
        write(6,'(/,t22,a)') repeat('=',36)
        write(6,'(t32,a,i7)') 'Linres calculation'
        write(6,'(t22,a)') repeat('=',36)

        write(6,'(a30,i7)') 'Linres: Initial perturbed atom', iai
        write(6,'(a30,i7)') 'Linres: Final perturbed atom', iaf
        write(6,'(a30,f10.5)') 'Linres: Tolerance (DM)', tolLR
      endif

      ! Nullify arrays
      nullify(dDold)

C Initialize ----------------------------------------------
      dynmat(:,:,:,:) = 0.0_dp
      atomdisp = ' '
      
C Dummies initialization------------------------------------------------------
      dummy_chargedensonly = .false.
      dummy_use_rhog_in = .false.

      call timer('Linres',1)
C-----------------------------------------------------------------------------

C Read Stored dynamical matrix files-----------------------------------------
      init = iai
      if (readold) then
        call readdynmat(init,iai,dynmat)
      end if

C Initialize dDscf, dEscf
      nullify(dDscf, dEscf)
      call re_alloc(dDscf, 1, maxnh, 1, nspin, 1, 3,
     &     'dDscf', 'linresscf')
      call re_alloc(dEscf, 1, maxnh, 1, nspin, 1, 3,
     &     'dEscf', 'linresscf')

C------------------------------LINRES MAIN LOOP-------------------------------
C Init of Linres calculation. External loop over perturbed atoms (IALR).
C Internal loop scf-loop (ISCF).

      do ialr = init, iaf 

C Read stored DM in file .LRDMIALR from previous non-converged calculation
        if (readold) then
          if (IOnode) then
            write(6,'(a)')
            write(6,'(a,i7)')'Linres: trying to start from .LRDM file'
            write(atomdisp,'(i0)') ialr
          endif
          fname = trim(slabel)//'.LRDM'//trim(atomdisp)
          inquire( file=fname, exist=found )
          if (found) then
            if (IOnode) then
              write(6,'(a,i7)') 'Linres: reading LRDM file for atom=',
     &                            ialr
            endif
            call read_dmlr( maxnh, no_l, nspin, numh,
     &                     listhptr, listh, dDscf, fname)
              write(6,'(a)')'Linres: Read DM from file...successfully!!'
          else
            write(6,'(a)')'Linres: LRDM file not found'//
     $                    ' or not corresponds to this atom'
          endif
        endif
        call timer('LRatom', 1)

        if (IOnode) then
          write(6,'(/,t22,a)') repeat('=',36)
          write(6,'(t32,a,i7)') 'Linres atom=', ialr
          write(6,'(t22,a)') repeat('=',36)
        endif

        first_LR = .true.       !init flag for first calls

        call init_pulay_arrays(iai)        

C Doing non-scf elements of the perturbed Hamiltonian and dinamical matrix
        dHmat0(:,:,:) = 0.0_dp
        
        call dhscf( nspin, no_s, iaorb, iphorb, no_l,
     .            no_u, na_u, na_s, isa, xa, indxua,
     .            ntm, 0, 0, 0, filesOut,
     .            maxnh, numh, listhptr, listh, Dscf, Datm,
     .            maxnh, H, Enaatm, Enascf, Uatm, Uscf, DUscf, DUext,
     .            Exc, Dxc, dipol, dummy_stress, dummy_fa,
     .            dummy_stress, dummy_use_rhog_in,
     .            dummy_chargedensonly, iai, iaf, ialr, lasto,
     .            dynmat, dDscf, dHmat, first_LR, dHmat0)
        print*,'dynmat 1 dhscf',dynmat
        print*,'dhmat0 dhscf',dhmat0(1:50,1,1)


C dHmat0 contains the pulay terms of the perturbed hamiltonian
C Dynamat contains all the elements that are non-scf and depend
C on the gradient potential terms

C Calculation of the kinetic terms of the perturbed hamiltonian
C (added into dHmat0)
        
      call dhinit(IALR, na_s, no_u, maxnh, maxnh, nspin, lasto,
     &       listh, listhptr, numh, DS, dSmat, dHmat0, dH)
        print*,'dhmat0 dhinit',dhmat0(1:50,1,1)
C Calculation of the KB terms of the perturbed hamiltonian
C (added into dHmat0)

        call dvnloc(scell, na_u, na_s, isa, xa, indxua, Dscf,
     &       maxnh, maxnh, lasto, lastkb, iphorb,
     &       iphKB, numh, listhptr, listh, numh,
     &       listhptr, listh, min(nspin,2), IALR,
     &       no_s, iaorb, dDscf, dHmat0, dynmat, first_LR)
        print*,'dhmat0 dvnloc',dhmat0(1:50,1,1)

C  Linear response SCF loop --------------------------
        iscf = 0
        converged = .false.
        do while ( iscf < nscf )

           ! Conditions of exit:
           !  -- At the top, to catch a non-positive nscf and # of iterations
           !  -- At the bottom, based on convergence
           iscf = iscf + 1

           first_LR = iscf == 1
           
C          Copy non-scf hamiltonian to total Hamiltonian
           dHmat(:,:,:)=0.0_dp
           do ispin = 1 , nspin
              dHmat(:,:,ispin) = dHmat0(:,:,ispin)
           end do

C Compute the perturbed potential elements Vxc, Vna and Vh and 
C include them into dHmat. Look that dhscf is called without dHmat0!!!
           
           call dhscf( nspin, no_s, iaorb, iphorb, no_l,
     .          no_u, na_u, na_s, isa, xa, indxua,
     .          ntm, 0, 0, 0, filesOut,
     .          maxnh, numh, listhptr, listh, Dscf, Datm,
     .          maxnh, H, Enaatm, Enascf, Uatm, Uscf, DUscf, DUext,
     .          Exc, Dxc, dipol, dummy_stress, dummy_fa,
     .          dummy_stress, dummy_use_rhog_in,
     .          dummy_chargedensonly, iai, iaf, ialr, lasto,
     .          dynmat, dDscf, dHmat, first_LR)
        print*,'dhmat iscf',dhmat(1:50,1,1)
           
           call re_alloc(dDold,1,maxnh,1,nspin,1,3,
     &          'dDold', 'linresscf')
           
C Copy current density as old one to perform the mixing
           dDold(1:maxnh,1:nspin,1:3) = dDscf(1:maxnh,1:nspin,1:3)
           dDscf(1:maxnh,1:nspin,1:3) = 0.0_dp
           dEscf(1:maxnh,1:nspin,1:3) = 0.0_dp

C Find change in density matrix from perturbed Hamiltonian and Overlap
           call delrho(no_s, na_s, nspin, maxnh, no_u, no_u, no_l, 
     &          GAMMA, indxuo, xijo, nspin, eo, tolLR, eigtolLR, nkpnt,
     &          kpoint, kweight, Qo, H, S, dHmat, dSmat, numh,
     &          listh, listhptr, ef, temp, dDscf, dEscf, iscf)

C Perform the density matrix mixing
           dMax = 0.0_dp
           do ix = 1 , 3
              call pulayx( iscf , .false. , no_l, maxnh, numh,
     &             listhptr, nspin,maxsav,wmix,nkick,
     &             wmixkick,dDscf(:,:,ix),dDold(:,:,ix),dDmax,ix)
#ifdef MPI
              dMax = max(dMax,dDmax)
              call globalize_max(dMax, buffer1)
              dMax = buffer1
#else
              dMax = max(dMax,dDmax)
#endif
           end do

           converged = dMax < tolLR

           ! Clean-up temporary memory
           call de_alloc(dDold, 'dDold', 'linresscf')

C Print error in the perturbed density
           if (IOnode) then 
             if ( first_LR ) then
               write(6,'(a12,a10)') 'iscf', 'dDmax'
             end if 
             write(*,'(a8,i4,f10.6)')'lr-scf:', iscf, dmax
           endif
          
C Write dDscf to file ------------------------------------------------
C DM file is named as: label.LRDM+'IALR'
           call write_dmlr(maxnh, no_l,nspin,numh,
     &          listhptr,listh,dDscf,ialr)

C     Check convergence cryteria
           if ( converged ) exit

       end do ! isc loop
C------------------------------------------------------------------------

       call timer('LRatom', 2)
       
       call resetPulayArrays()
       
C Adding actions if iscf loop ends but ddscf is NOT converged
C In this case, the FC matrix elements that depend on the dDscf 
C are not calculated
C----------------------------------------------------------------
       if ( .not. converged ) exit

C Now, is time to calculate the final terms of the dynamical matrix
C------------------------------------------------------------------------
        print*,'ddscf conver',ddscf(1:50,1,1)

C     Add non-local potential and kinetic part contributions to dynmat---
C     (the terms that depend on the perturbed density and the ones that 
C     do not depend)
C     Note that first_LR now is false
        print*,'dynmat antes dvnloc final',dynmat
!	dynmat=0.0_dp
!	print*,'dynmat set to 0 antes de dvnloc'
       call dvnloc(scell, na_u, na_s, isa, xa, indxua, Dscf,
     &      maxnh, maxnh, lasto, lastkb, iphorb,
     &      iphKB, numh, listhptr, listh, numh,
     &      listhptr, listh, min(nspin,2), IALR,
     &      no_s, iaorb, dDscf, dHmat, dynmat, first_LR)
        print*,'dynmat dvnloc final',dynmat
!	stop

C     Add kinetic and ovelap contributions to dynmat -----------------
       call ddsmat(dH, dS, d2H, d2S, nspin, no_s, no_u, na_s, na_u,
     &      indxuo, lasto, numh, listh, listhptr, Dscf, Escf,
     &      ialr, iaorb, maxnh, dDscf, dEscf, dynmat)
                      print*,'dynmat ddsmat final',dynmat

c     Add terms that depends on the perturbed density and Vscf 
C     and perturbed Vscf potentials-----------------------------------
!      dynmat=0.0 
       call dhscf( nspin, no_s, iaorb, iphorb, no_l,
     .      no_u, na_u, na_s, isa, xa, indxua,
     .      ntm, 0, 0, 0, filesOut,
     .      maxnh, numh, listhptr, listh, Dscf, Datm,
     .      maxnh, H, Enaatm, Enascf, Uatm, Uscf, DUscf, DUext,
     .      Exc, Dxc, dipol, dummy_stress, dummy_fa,
     .      dummy_stress, dummy_use_rhog_in,
     .      dummy_chargedensonly, iai, iaf, ialr, lasto,
     .      dynmat, dDscf, dHmat)
                      print*,'dynmat dhscf final',dynmat
       
C Dealloc the non-scf variables (drhoatm, drhoscf0 and dvxc) for next atom
       call resetFirstmatrices()

C Save in a file the dynamical matrix and extra flag of the index ialr of the last atom
       call writedynmat(iai,iaf,ialr,dynmat,.false.)
       
      end do                    ! ialr atoms loop
      call timer('Linres',2)
C-------------------------------------------------------------------------

C Dealloc matrices 
      call de_alloc(dDscf, 'dDscf', 'linresscf')
      call de_alloc(dEscf,'dEscf','linresscf')

C There is other contribution to the dynamical matrix (non-scf)
C the Laplacian of the neutral atom potential
      if ( converged ) then
         call ddnaefs(na_u, na_s, scell, xa, indxua, rmaxv,
     &        isa, dynmat)
      end if
C-------------------------------------------------------------------------
C-----------------Finally, print the dynamical matrix to FC file----------
C-------------------------------------------------------------------------

      print *, '----- COMPLETE DYNAMICAL MATRIX -----'
      print *,'Beta:',';','xyz',';','Alpha',';','xyz',';','dynmat'
      do ialr = IAI, IAF
       do ix = 1,3
        do j = 1, na_u
         do jx = 1,3
          print *, j,';',jx,';',ialr,';',ix,';',
     &                   dynmat(jx,j,ix,ialr)
         enddo
        enddo
       enddo
      enddo

      call writedynmat(iai,iaf,ialr,dynmat,.true.)

      end subroutine 

      end MODULE m_linresscf

