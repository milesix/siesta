! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
      Program SIESTA

      use m_siesta_init
      use m_siesta_analysis
      use m_siesta_move
      use m_siesta_end
      use m_siesta_forces
      use m_siesta_tddft
      use m_ts_global_vars, only: onlyS

      USE m_steps, only: inicoor, fincoor
      use precision, only: dp
      use m_forces, only: fa
      use siesta_geom, only: xa, va, na_u
      use siesta_options, only: virtual_dt

      use siesta_options, only: td_elec_dyn, want_virtual_step_md
      use siesta_options, only: virtual_md_translate
      !use m_virtual_step_data, only: substep, BASE_STEP, VIRTUAL_STEP
      use m_virtual_step_data
      use m_virtual_step_procs, only: base_md_step_logic,
     $     virtual_md_step_logic

      use parallel, only: SIESTA_worker ! whether part of Siesta's communicator
      use parallel, only: ionode
#ifdef MPI
      use mpi_siesta, only: mpi_comm_dft ! Includes Solver nodes
#endif
      use m_mpi_utils, only: broadcast

#ifdef TRACING_SOLVEONLY
      use extrae_module
#endif
      use m_io_yaml, only: siesta_write_yaml

      implicit none

      integer :: istep
      logical :: relaxd

      external :: bye

! Notes for PEXSI and ELSI operation (extra procs for Solver)

!
! A subset of nodes carries out non-Solver Siesta operations
! (i.e., setting up H, moving atoms, diagonalizing...).
! These are tagged as "SIESTA_worker" (admittedly, a bad name)
! All nodes are involved in the (PEXSI or ELSI) electronic-structure solver,
! and in the new LocalDOS computation based on selected inversion.
!
! 'siesta_init', 'siesta_forces', and 'siesta_analysis' need to
! be called by all nodes.
!
! In some cases, the result of a computation by "SIESTA_worker" nodes needs
! to be broadcast to guarantee proper control-flow logic (for example,
! the "relaxd" variable from 'siesta_move'.
!
!----------------------------------------------------------------- BEGIN
!      if (ionode) call memory_snapshot("at start of program")

      call siesta_init()
C     V.D.: thtr debug
      if (virtual_md_translate) call equalize_velocities()

#ifdef TRACING_SOLVEONLY
      call extrae_shutdown
#endif

#ifdef SIESTA__PEXSI
      if (ionode) call memory_snapshot("after siesta_init")
#endif

C     Begin of coordinate relaxation iteration
      relaxd = .false.

#if defined (SIESTA__PEXSI) || defined (SIESTA__ELSI)
      ! Broadcast relevant things for program logic
      ! These were set in siesta_options, called only by "SIESTA_workers".
      call broadcast(inicoor,comm=mpi_comm_dft)
      call broadcast(fincoor,comm=mpi_comm_dft)
#endif
      istep  = inicoor

      ! If performing MD with `virtual` step: initialize
      if ( want_virtual_step_md ) then ! do perform two-step MD

         if ( td_elec_dyn )
     $        call die("No TDDFT with thermal transport...")

         call init_virtual_step_md()  ! allocate xa,va,fa aux arrays

         DO WHILE ((istep.le.fincoor) .and. (.not. relaxd))

              istep_vmd = istep
              substep = BASE_STEP

              call siesta_forces( istep )

C     V.D.: thtr debug
              if (virtual_md_translate) call equalize_velocities()

              call base_md_step_logic() ! Perform whatever needs to be done now
              call store_state_before_move() ! Save xa, fa, va

              ! ** Call siesta_move here to record proper step values
              if (SIESTA_worker) call siesta_move( istep, relaxd )

              call store_state_after_move() ! Store future base coords
              ! call store_data_from_base_move()  !-> Stores data we need from the base step on `istep` iteration
              ! TODO: make a separate module for it - already have `ks_flux_data`

              ! To get Jks flux component we need data from the `istep`-iteration base step:
              !  * S
              !  * Dscf
              !  * coeffs
              !  * gradS
              ! as well as the data derived via virtual step:
              !  * Dscf_deriv

              ! We calculate the components on virtual step, but using the data mostly
              ! from the preceding base step - thus we need to store and pass it further.

              call make_virtual_md_move() ! Generate xa for virtual step
              call write_virtual_step_state() ! Append to files with `virtual`-step coords and velocities


              substep = VIRTUAL_STEP

              call siesta_forces( istep )

C     V.D. thtr debug
              if (virtual_md_translate) call check_pcm()

              print*, "[VMD.Logic.Start] ", istep_vmd
              call virtual_md_step_logic()
              print*, "[VMD.Logic.Stop] ", istep_vmd
              ! This should essentially do two things:
              ! 1. compute derivatives e.g. Dscf_deriv
              ! 2. compute Jks-flux components for `istep` iteration,
              !    using the data from preceding base step.
              ! 3. reset specific data storages for the next `Base` step

              call re_store_state_after_move()  ! Restore xa to xa_future, fa to fa_saved

              if (.not. relaxd) then
                 istep = istep + 1
              endif

#if defined (SIESTA__PEXSI) || defined (SIESTA__ELSI)
              call broadcast(relaxd,comm=mpi_comm_dft)
#endif

         ENDDO ! while loop

      else ! Skip `Virtual` step - generic SIESTA cycle

         DO WHILE ((istep.le.fincoor) .AND. (.not. relaxd))

            if ( td_elec_dyn ) then
               call siesta_tddft( istep )
            else
               call siesta_forces( istep )
            end if


         if ( onlyS ) then
             call bye("Saved S only, exiting.")
         end if

        if (SIESTA_worker) call siesta_move( istep, relaxd )

#if defined (SIESTA__PEXSI) || defined (SIESTA__ELSI)
            call broadcast(relaxd,comm=mpi_comm_dft)
#endif
            if (.not. relaxd) then
               istep = istep + 1
            endif
#ifdef SIESTA__PEXSI
            if (ionode) call memory_snapshot("after geometry step")
#endif

         ENDDO

      endif ! One/two-step MD

C     End of coordinate-relaxation loop
      call siesta_analysis( relaxd )
#ifdef SIESTA__PEXSI
      if (ionode) call memory_snapshot("after siesta_analysis")
#endif

      if (ionode) call siesta_write_yaml()

      call siesta_end()
#ifdef SIESTA__PEXSI
      if (ionode) call memory_snapshot("at end of program")
#endif

!-----------------------------------------------------------------------END
      contains

      subroutine init_virtual_step_md ()
        !! Initializes storage arrays for atomic forces, coordinates and
        !! velocities used to save system state in between `base` moves,
        !! as well as to perform `virtual` move.

        ! allocate(xa_before_move(3, na_u), source=xa)
        !NOTE: a bug in gfortran-4.x leads to mandatory typing of the
        ! shape of allocated arrays even when `source` is present.

        allocate(xa_before_move, source=xa)
        allocate(xa_after_move,  source=xa)
        allocate(fa_before_move, source=fa)
        allocate(fa_after_move,  source=fa)
        allocate(va_before_move, source=va)
        allocate(va_after_move,  source=va)

      end subroutine init_virtual_step_md

      subroutine store_state_before_move ()
        !! Store system state (forces, coords, velocities) in corresponding arrays.
        !! Should be called before `[[siesta_move(proc)]]` during `Base` step.
        xa_before_move = xa
        fa_before_move = fa
        va_before_move = va
      end subroutine store_state_before_move

      subroutine store_state_after_move ()
        !! Store system state (forces, coords, velocities) in corresponding arrays.
        !! Should be called after `[[siesta_move(proc)]]` during `Base` step.
        xa_after_move = xa
        fa_after_move = fa
        va_after_move = va
      end subroutine store_state_after_move

      subroutine make_virtual_md_move ()
        !! Obtain state of the system in `Virtual` position
        !! by moving appropriately the coordinates:
        !! x' = x + v*delta_t
        integer :: ia

        do ia=1,na_u
            xa(:,ia) =
     $          xa_before_move(:,ia) + va_before_move(:,ia)*virtual_dt
        end do
      end subroutine make_virtual_md_move


      subroutine write_virtual_step_state ()
        use files,      only : slabel, label_length
        external          io_assign, io_close
        integer :: ia, iu, iv, ix
        character(len=label_length+7) ::   fname

        if (ionode) then
           fname = trim(slabel) // '.vmd.X1'
           call io_assign( iu )
           open( iu, file=fname, form='formatted', status='unknown',
     $          position='append')
           do ia = 1,na_u
              write(iu,'(3x,3f18.12)')  (xa_before_move(ix,ia),ix=1,3)
              write(*,'(a,i3.1,a,3f18.12)') '[XVS] step: ', istep_vmd,
     $             '  X1 ', (xa_before_move(ix,ia),ix=1,3)
           enddo
           call io_close( iu )

           fname = trim(slabel) // '.vmd.X2'
           call io_assign( iu )
           open( iu, file=fname, form='formatted', status='unknown',
     $          position='append')
           do ia = 1,na_u
              write(iu,'(3x,3f18.12)')  (xa(ix,ia),ix=1,3)
              write(*,'(a,i3.1,a,3f18.12)') '[XVS] step: ', istep_vmd,
     $             '  X2 ', (xa(ix,ia),ix=1,3)
           enddo
           call io_close( iu )

           fname = trim(slabel) // '.vmd.VS'
           call io_assign( iu )
           open( iu, file=fname, form='formatted', status='unknown',
     $          position='append')
           do ia = 1,na_u
              write(iu,'(3x,3f18.12)')  (va_before_move(iv,ia),iv=1,3)
              write(*,'(a,i3.1,a,3f18.12)') '[XVS] step: ', istep_vmd,
     $             '  VS ', (va_before_move(iv,ia),iv=1,3)
           enddo
           call io_close( iu )
        end if
      end subroutine write_virtual_step_state


      subroutine re_store_state_after_move ()
        !! In the end of `Virtual` step (should be called there) -
        !! restore forces and positions to continue the program flow.
        xa = xa_after_move
        fa = fa_after_move
        va = va_after_move
      end subroutine re_store_state_after_move

      subroutine check_pcm ()
        !! Debugging routine - checks the value of the system-s
        !! center of mass momentum.
        use atomlist, only : amass, qtot
        real(dp) :: pcm(3), massi, mtot
        integer  :: i

        pcm(:) = 0.0_dp
        mtot   = 0.0_dp

        do i=1,na_u
           massi = amass(i)
           mtot = mtot + massi

           pcm(:) = pcm(:) + massi * va(:,i)
        enddo
        print*, "[pcom] ", pcm(:)
        print*, "[vcom] ", pcm(:)/mtot
        print*, "[chdrift]", pcm(1:3)*qtot/mtot
      end subroutine check_pcm


      subroutine subtract_vcom ()
        use atomlist, only : amass
        real(dp) :: pcm(3), massi, mtot
        integer  :: i, ix

        pcm(:) = 0.0_dp
        mtot   = 0.0_dp

        do i=1,na_u
           massi = amass(i)
           mtot = mtot + massi

           pcm(:) = pcm(:) + massi * va(:,i)
        enddo
        do i = 1,na_u
           do ix = 1,3
              va(ix,i) = va(ix,i) - pcm(ix)/mtot
           enddo
        enddo
      end subroutine subtract_vcom

      subroutine equalize_velocities ()
        !! set all velocities in the system equal to the velocity of the first atom
        use siesta_options, only: virtual_md_translate_vx
        use siesta_options, only: virtual_md_translate_vy
        use siesta_options, only: virtual_md_translate_vz

        integer  :: i, ix

c$$$        va(1,:) = 0.02_dp
        va(1,:) = virtual_md_translate_vx
        va(2,:) = virtual_md_translate_vy
        va(3,:) = virtual_md_translate_vz
c$$$        print*, "[uniform_translation_vel]", va(:,1)

      end subroutine equalize_velocities

      END program siesta
