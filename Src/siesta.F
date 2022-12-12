      Program SIESTA
	

      use m_siesta_init
      use m_siesta_analysis
      use m_siesta_move
      use m_siesta_end
      use m_siesta_forces
      USE m_steps, only: inicoor, fincoor, fincoor_Verlet, 
     .                 inicoor_bath, fincoor_bath,
     .                 inicoor_Verlet,sys_prep,TwoBaths 

      implicit none
      integer :: natoms_1
      integer :: istep
      logical :: relaxd

      character(15) :: typeCall 
!------------------------------------------------------------------------- BEGIN
	call siesta_init()

C DEBUG
!	print *, 'sys prep', sys_prep 
!	print *, 'TwoBaths ', TwoBaths  
C

C     Begin of coordinate relaxation iteration
      if (TwoBaths .eqv. .true.) then
	if (sys_prep .eq. 'Nose') then
	 !NOSE of part1 and part2
	print *, 'Debug: Nose bath'
	print *, 'Tinicial', inicoor_bath
	print *, 'Tfinal', fincoor_bath
	 relaxd = .false.
	 istep  = inicoor_bath
	 if (inicoor_bath .lt. fincoor_bath) then !nos saltamos nose
		typeCall='Nose'
		DO WHILE ((istep.le.fincoor_bath) .and. (.not. relaxd))
			call siesta_forces( istep ) !independent of the relaxation method
			call siesta_move( istep, relaxd, typeCall)
		    if (.not. relaxd) then
	      	    	istep = istep + 1
	      	    endif
	      	ENDDO
	 endif
	elseif (sys_prep .eq. 'Vresc') then
	print *, 'Debug: Velocity rescaling'
	print *, 'Tinicial', inicoor_bath
	print *, 'Tfinal', fincoor_bath
	 !Vresc of part1 and part2
	 relaxd = .false. 
	 istep  = inicoor_bath
	 if (inicoor_bath .lt. fincoor_bath) then !nos saltamos Vresc
		typeCall='Vresc'
		DO WHILE ((istep.le.fincoor_bath) .and. (.not. relaxd))
			call siesta_forces( istep ) !independent of the relaxation method
			call siesta_move( istep, relaxd, typeCall)
		    if (.not. relaxd) then
	      	    	istep = istep + 1
	      	    endif
	      	ENDDO
	 endif
	endif !sys_prep if
	!VERLET relaxation of whole part (no system preparation or ready for relax) 
	print*, 'Debug: Verlet'
	print *, 'Tinicial', inicoor_Verlet
	print *, 'Tfinal', fincoor_Verlet
	istep  = inicoor_Verlet
	relaxd = .false.
	if (inicoor_Verlet .lt. fincoor_Verlet) then !nos saltamos Verlet
	  typeCall='Verlet'
	  DO WHILE ((istep.le.fincoor_Verlet) .AND. (.not. relaxd))
		call siesta_forces( istep ) !independent of the relaxation method
		call siesta_move( istep, relaxd, typeCall)
	    if (.not. relaxd) then
      	    	istep = istep + 1
      	    endif
      	  ENDDO
	endif

      else
      	relaxd = .false.
      	istep  = inicoor
      	DO WHILE ((istep.le.fincoor) .AND. (.not. relaxd))
      	  call siesta_forces( istep )
      	  call siesta_move( istep, relaxd, typeCall )
      	  if (.not. relaxd) then
      	    istep = istep + 1
      	  endif
      	ENDDO

      endif


      ! End of coordinate-relaxation loop 
      call siesta_analysis( relaxd )

      call siesta_end()
!--------------------------------------------------------------------------- END
      END program siesta
