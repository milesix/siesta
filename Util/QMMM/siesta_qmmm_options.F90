
MODULE siesta_qmmm_options

  USE precision
  USE units
  use fdf
  use sys
  use m_qmmm_fdf, only : fdf_block_qmmm 

  implicit none

  PUBLIC

  logical :: default       ! Temporary used to pass default values in fdf reads

  logical :: launch_siesta_flag

  integer :: nfce
  integer :: wricoord
  integer :: mmsteps

  integer :: mneigh_freq

  logical :: center_qm_system
  logical :: siesta_qmmm_usesavexv

  logical, save :: opt_idyn=.false.
  logical, save :: fc_idyn=.false.
  logical, save :: therm_idyn=.false.
  logical, save :: include_qmqm_vdw=.false.

CONTAINS

  subroutine read_siesta_qmmm_options(na_qm,nat)

    integer na_qm, nat, wrifces

    character(len=22) :: dyntype

    integer :: iunit

#ifndef QMMM_BSC
    logical :: leqi
    external ::  leqi
#endif

    ! Write coordinate variable 
    wricoord = fdf_integer('WriteCoordinates',1)
    write(6,'(a,i5,a)') &
         'read: Write coordinates each           = ',wricoord, '  steps'

    ! MMxQM steps
    mmsteps = fdf_integer('MD.MMxQMsteps',1)
    write(6,'(a,i5)') &
         'read: MM x QM steps                    = ',mmsteps

    ! Sets the atoms for whom the forces will be writen
    wrifces = fdf_integer('WriteForcesQMMM',1)
    if(wrifces==0) then
       nfce = 0
    elseif(wrifces==1) then
       nfce = na_qm
    elseif(wrifces==2) then
       nfce = nat
    else
       call die('read: WriteForcesQMMM could only be 0, 1 or 2')
    endif
    write(6,'(a,i5)') &
         'read: Write forces QM-MM               = ',wrifces

    launch_siesta_flag=fdf_boolean('LaunchSiesta',.false.)

    ! Options readed here instead of siesta_options
    default=fdf_boolean('UseSaveData'     , .false.)
    siesta_qmmm_usesavexv=fdf_boolean('MD.UseSaveXV', default)

    mneigh_freq=fdf_integer('NebListFreq',10)

    center_qm_system=fdf_boolean('CenterQmSystem',.false.)

    dyntype= fdf_string('MD.TypeOfRun','cg')

    if (leqi(dyntype,'cg').or.leqi(dyntype,'broyden').or.leqi(dyntype,'fire')) then
       opt_idyn=.true.
    else if (leqi(dyntype,'fc')) then
       fc_idyn=.true.
    else
       therm_idyn=.true.
    endif

    include_qmqm_vdw = fdf_boolean('QMMM.IncludeQmQmVDW',.false.)

    if (include_qmqm_vdw) then
       if (fdf_block_qmmm('MM.Potentials',iunit)) &
            call die('QMMM.IncludeQmQmVDW and block MM.Potentials'//&
            'cannnot appear in the fdf input file at the same time')
    endif

    write(6,102)   
    
2   format(a)
5   format(a,i5,a)
6   format(a,f10.4,a)
100 format(/,'read: ',73(1h*))
101 format('read:                  INPUT ERROR')
102 format('read: ',73(1h*))
103 format('read: ',i4,2x,3f10.5,i3) 

  end subroutine read_siesta_qmmm_options

END MODULE siesta_qmmm_options

