module linkatoms

  use precision, only: dp

  implicit none

  logical, save :: linkatom
  integer, save :: numlink=0, NSpecialLinkAtoms=0
  character(len=2), save :: default_link_symbol='H?'
  integer, save :: default_link_type=0
  integer, dimension(:,:), allocatable, save::linkqm,linkmm
  integer, dimension(:,:,:), allocatable, save::linkmm2
  real(dp), save :: Elink
  real(dp), dimension(:,:,:), allocatable, save :: fa_link
  real(dp), dimension(:), allocatable, save::distl,kch,rch
  real(dp), dimension(:,:), allocatable, save::pclinkmm,Emlink
  character, dimension(:,:), allocatable, save::linkqmtype*4

  integer, dimension(:,:,:,:), allocatable :: parametro
  real(dp), dimension(:,:), allocatable :: link_bondeq
  real(dp), dimension(:,:), allocatable :: link_kbond
  real(dp), dimension(:,:,:), allocatable :: link_kangle

  character(len=5) :: h_link_scaling_type='splam'
  integer, dimension(:), allocatable, save :: num_linkqm, num_linkmm 

  logical, save :: got_numlink=.false.

  ! C. F. Sanz-Navarro (Nov. 2009)
  ! this module deals with the case of having resonance structuresin the QM 
  ! part, which cannot be accurately described by 1 type of link atoms and 
  ! needs the linear combination of different ways (single and double bonds 
  ! with the link atoms) to be able to model the system

  integer :: num_resonances=1
  type resonance_type
     real(dp) :: weight   
     character(len=4), dimension(:), pointer :: atype
     integer, dimension(:), pointer :: isa
     integer, dimension(:), pointer :: iza
     character(len=13) :: path
  end type resonance_type
  type(resonance_type), dimension(:), allocatable, save :: resonance

  real(dp), dimension(:,:,:), allocatable, save :: rlink

  integer, dimension(:), allocatable :: special_linkqm, special_linkmm
  character(len=4), dimension(:,:), allocatable :: special_linkatype
  integer, dimension(:,:), allocatable :: special_isa
  real(dp), dimension(:), allocatable :: special_weight

contains

  subroutine read_special_link_atoms

    use fdf
    use m_qmmm_fdf, only : fdf_block_qmmm    
    use sys

    integer :: i, j
    real(dp) :: total_weight

    character :: acf*22, acf_default*22
    integer :: iscale, iunit                                      
#ifndef QMMM_BSC
    logical :: leqi
#endif

    ! Format of atomic coordinates
    acf_default = 'Ang'
    acf = fdf_string('AtomicCoordinatesFormat',acf_default)

    if (leqi(acf,'NotScaledCartesianBohr') .or. &
         leqi(acf,'Bohr') ) then
       iscale = 0
       write(6,'(a,a)') &
            'read: Atomic-coordinates input format  = ',&
            '    Cartesian coordinates (in Bohr)'
    else if (leqi(acf,'NotScaledCartesianAng') .or.&
         leqi(acf,'Ang') ) then
       iscale = 1
       write(6,'(a,a)')&
            'read: Atomic-coordinates input format  = ',&
            '    Cartesian coordinates (in Ang)'
    else
       write(6,"(/,'read: ',73(1h*))")
       write(6,"('read:                  INPUT ERROR')")
       write(6,'(2a)') 'read: You must use one of the following',&
            'coordinate options:'
       write(6,'(a)') 'read:     - NotScaledCartesianBohr (or Bohr)'
       write(6,'(a)') 'read:     - NotScaledCartesianAng (or Ang) '
       write(6,"('read: ',73(1h*))")
       call die
    endif

    if (fdf_block_qmmm('SpecialLinkAtoms',iunit)) then

       read(iunit,*)NSpecialLinkAtoms,num_resonances

       total_weight=0.0_dp
       read(iunit,*) (special_weight(i),i=1,num_resonances)
       do i=1,num_resonances
          total_weight=total_weight+special_weight(i)
       enddo
       if (abs(1.0_dp-total_weight)>0.0001_dp) then
          write(*,*)'WARNING: The sum of the weights of the ',&
               'resonances give: ',total_weight   
       endif

       do j=1,NSpecialLinkAtoms
          read(iunit,*)(special_isa(i,j),special_linkatype(i,j),i=1,num_resonances),&
               special_linkqm(j),special_linkmm(j)
       enddo

    endif

  end subroutine read_special_link_atoms

  subroutine get_numlink(na_qm,na_mm,rqm,rmm,attype,qmtype,nparm,nbond,bondtype,bondeq,&
                  deb,cell,lattice_type)

    use fdf
    use m_qmmm_fdf, only : fdf_block_qmmm 
    use sys

    integer na_qm, na_mm
    real(dp) :: rqm(3,na_qm),rmm(3,na_mm),rmin(na_qm),r1,dist
    character :: attype(na_mm)*4,qmtype(na_qm)*4,ch4*4,ch1*1
    integer :: nparm,nbond
    character ::  bondtype(nparm)*5
    real(dp) :: bondeq(nparm)
    real(dp) :: cell(3,3), kcell(3,3)
    real(dp) :: dx1, dy1, dz1, dist_v2
    character lattice_type
    logical deb
    real(dp) :: link_rmin
    integer :: min(na_qm)
    integer, dimension(na_qm) :: qm2link
    integer :: iunit
 
    integer :: parameter_set
    real(dp) :: qmmm_range=0.0d0
    character :: tybond*5, chi*4, chj*4

    integer :: i, j, l

    got_numlink=.true.

    if (na_qm==0) return

    call reccel(3,cell,kcell,0)

    qm2link=0

    if (fdf_block_qmmm('SpecialLinkAtoms',iunit)) then

        read(iunit,*)NSpecialLinkAtoms,num_resonances

        if (NSpecialLinkAtoms>0) then

           if (allocated(special_linkqm)) deallocate(special_linkqm)
           allocate(special_linkqm(NSpecialLinkAtoms))
           if (allocated(special_linkmm)) deallocate(special_linkmm)
           allocate(special_linkmm(NSpecialLinkAtoms))
           if (allocated(special_linkatype)) deallocate(special_linkatype)
           allocate(special_linkatype(num_resonances,NSpecialLinkAtoms))
           if (allocated(special_isa)) deallocate(special_isa)
           allocate(special_isa(num_resonances,NSpecialLinkAtoms))
           if (allocated(special_weight)) deallocate(special_weight)
           allocate(special_weight(num_resonances))

           call read_special_link_atoms

           do i=1,NSpecialLinkAtoms
              qm2link(special_linkqm(i))=1
           enddo

        endif
    endif

    numlink=NSpecialLinkAtoms

    ! Count total number of link atoms
    do i=1,na_qm
       if (qm2link(i)==1) cycle
       dx1=rqm(1,i)-rmm(1,1)
       dy1=rqm(2,i)-rmm(2,1)
       dz1=rqm(3,i)-rmm(3,1)
       call pbc_displ_vector(lattice_type,cell,kcell,dx1,dy1,dz1)
       rmin(i)=dist_v2(dx1,dy1,dz1)
       min(i)=1
       do j=2,na_mm
          dx1=rqm(1,i)-rmm(1,j)
          dy1=rqm(2,i)-rmm(2,j)
          dz1=rqm(3,i)-rmm(3,j)
          call pbc_displ_vector(lattice_type,cell,kcell,dx1,dy1,dz1)
          r1=dist_v2(dx1,dy1,dz1)
          if(r1.le.rmin(i)) then
             rmin(i)=r1
             min(i)=j
          endif
       enddo

       if (rmin(i)>2.8/0.529177) cycle  ! Only consider pairs of atoms close enough

       chi=qmtype(i)
       if (chi(1:1)=='H') cycle
       j=min(i)
       chj=attype(j)
       if (chj(1:1)=='H') cycle
       parameter_set=0
       qmmm_range=0.0d0
       do l=1,nbond
          tybond=bondtype(l)
          if(chi(1:2)==tybond(1:2).and.chj(1:2)==tybond(4:5)) then
             qmmm_range=1.20*bondeq(l)/0.529177
             parameter_set=1
             exit
          else if (chi(1:2)==tybond(4:5).and.chj(1:2)==tybond(1:2)) then
             qmmm_range=1.20*bondeq(l)/0.529177
             parameter_set=1
             exit
          endif
       enddo
       if (parameter_set==0) then
          if (check_out_of_range(i,j,qmtype(i),attype(j),rmin(i)*0.529177)) cycle
       endif
       ch4=qmtype(i)
       ch1=ch4(1:1)
       if (ch1=='O' .and. rmin(i)<qmmm_range) then
          if(deb) write(6,*) i,min(i),rmin(i)
          numlink=numlink+1
       else if ((ch1=='C'.or.ch1=='N').and. rmin(i)<qmmm_range) then  ! distance in Bohrs
          if(deb) write(6,*) i,min(i),rmin(i)
          numlink=numlink+1
       endif
    enddo

  end subroutine get_numlink

  ! Reads Link Atoms parameters
  subroutine link1(ng1,na_qm,na_mm,attype,qmtype,nparm,nbond,bondtype,bondeq,rqm,rmm,nspec,deb,&
                           cell,lattice_type)

    use fdf
    use m_fdf_global, only: fdf_global_get
    use m_qmmm_fdf, only : fdf_block_qmmm
    use sys

    integer na_qm,na_mm,i,j,k,l,m,iunit,ng1(na_mm,6),min(na_qm)
    integer :: nspec
    character(len=2) :: atype
    character(len=2), dimension(nspec) :: atsym
    integer, dimension(nspec) :: atnum
    real(dp) :: link_rmin
    character :: exp,attype(na_qm)*4,qmtype(na_qm)*4,ch4*4,ch1*1,ch2*2
    integer :: nparm, nbond
    character ::  bondtype(nparm)*5
    real(dp) :: bondeq(nparm)
    real(dp) :: rqm(3,na_qm),rmm(3,na_mm),rmin(na_qm),r1,dist
    logical deb
    real(dp) :: cell(3,3), kcell(3,3)
    real(dp) :: dx1, dy1, dz1, dist_v2
    character lattice_type
    integer :: default_link_isa, default_link_iza
    integer, dimension(na_qm) :: qm2link

    integer :: parameter_set
    real(dp) :: qmmm_range=0.0d0
    character :: tybond*5, chi*4, chj*4 

    if (na_qm==0) return

    call reccel(3,cell,kcell,0)

    if (.not.got_numlink) call die('link1: Implementation error; call get_numlink first.')

    if (numlink==0) return

    default_link_symbol=fdf_string('DefaultLinkSymbol','H?')

    ! Read the atomic labels
    if ( fdf_block_qmmm('ChemicalSpeciesLabel',iunit) ) then
       do i=1,nspec
          read(iunit,*,err=2,end=2) j,atnum(i),atsym(i)
       enddo
    else 
       call die("link1: You must specify the atomic labels")
    endif

    do k=1,nspec
       atype=atsym(k)
       if (default_link_symbol(1:1)==atype(1:1)) then
          default_link_isa=k
          default_link_iza=atnum(k)
       endif
    enddo

    qm2link=0

    if (fdf_block_qmmm('SpecialLinkAtoms',iunit)) then

        read(iunit,*)NSpecialLinkAtoms,num_resonances

        if (NSpecialLinkAtoms>0) then

           if (allocated(special_linkatype)) deallocate(special_linkatype)
           allocate(special_linkatype(num_resonances,NSpecialLinkAtoms))
           if (allocated(special_isa)) deallocate(special_isa)
           allocate(special_isa(num_resonances,NSpecialLinkAtoms))
           if (allocated(special_weight)) deallocate(special_weight)
           allocate(special_weight(num_resonances))

           call read_special_link_atoms

           do i=1,num_resonances
              resonance(i)%weight=special_weight(i)
           enddo
 
           ! assignation of CQM and CMM for special link atoms
           do i=1,NSpecialLinkAtoms
              linkqm(i,1)=special_linkqm(i)
              qm2link(special_linkqm(i))=1
              linkmm(i,1)=special_linkmm(i)!-na_qm
           enddo

           do i=1,num_resonances
              !  Assignates atomic number (iza)
              do j=1,NSpecialLinkAtoms
                 resonance(i)%atype(j)=special_linkatype(i,j)
                 resonance(i)%isa(j)=special_isa(i,j)
                 do k=1,nspec
                    if(resonance(i)%isa(j)==k) resonance(i)%iza(j)=atnum(k)
                 enddo
                 if(resonance(i)%iza(j)==0) then
                    call die("initiate_resonances: There are "//&
                         " atoms without atomic number")
                 endif
              enddo
           enddo
           deallocate(special_isa)
           deallocate(special_linkatype)

        endif

    endif

    if (numlink>NSpecialLinkAtoms) then
       do j=1,num_resonances
          do i=NSpecialLinkAtoms+1,numlink
             resonance(j)%atype(i)=default_link_symbol
             resonance(j)%isa(i)=default_link_isa
             resonance(j)%iza(i)=default_link_iza
          enddo
       enddo
    endif

    k=NSpecialLinkAtoms
    ! assignation of CQM and CMM 
    do i=1,na_qm
       if (qm2link(i)==1) cycle
       dx1=rqm(1,i)-rmm(1,1)
       dy1=rqm(2,i)-rmm(2,1)
       dz1=rqm(3,i)-rmm(3,1)
       call pbc_displ_vector(lattice_type,cell,kcell,dx1,dy1,dz1)
       rmin(i)=dist_v2(dx1,dy1,dz1)
       min(i)=1
       do j=2,na_mm
          dx1=rqm(1,i)-rmm(1,j)
          dy1=rqm(2,i)-rmm(2,j)
          dz1=rqm(3,i)-rmm(3,j)
          call pbc_displ_vector(lattice_type,cell,kcell,dx1,dy1,dz1)
          r1=dist_v2(dx1,dy1,dz1)
          if(r1<=rmin(i)) then
             rmin(i)=r1
             min(i)=j
          endif
       enddo

       if (rmin(i)>2.8/0.529177) cycle  ! Only consider pairs of atoms close enough

       chi=qmtype(i)
       if (chi(1:1)=='H') cycle
       j=min(i)
       chj=attype(j)
       if (chj(1:1)=='H') cycle
       parameter_set=0
       qmmm_range=0.0d0
       do l=1,nbond
          tybond=bondtype(l)
          if(chi(1:2)==tybond(1:2).and.chj(1:2)==tybond(4:5)) then
             qmmm_range=1.20*bondeq(l)/0.529177
             parameter_set=1
             exit
          else if (chi(1:2)==tybond(4:5).and.chj(1:2)==tybond(1:2)) then
             qmmm_range=1.20*bondeq(l)/0.529177
             parameter_set=1
             exit
          endif
       enddo
       if (parameter_set==0) then
          if (check_out_of_range(i,j,chi,chj,rmin(i)*0.529177)) cycle
       endif
       ch4=qmtype(i)
       ch1=ch4(1:1)
       if (ch1=='O'.and.rmin(i)<qmmm_range) then  ! distance in Bohrs
          rmin(i)=20.0_dp
          if(deb) write(6,*) i,min(i),rmin(i)
          k=k+1
          if (k>numlink) call die('link1: k>numlink???')
          linkqm(k,1)=i
          linkmm(k,1)=min(i)
       else if ((ch1=='C'.or.ch1=='N').and.rmin(i)<qmmm_range) then  ! distance in Bohrs
          rmin(i)=20.0_dp
          if(deb) write(6,*) i,min(i),rmin(i)
          k=k+1
          if (k>numlink) call die('link1: k>numlink???')
          linkqm(k,1)=i
          linkmm(k,1)=min(i)
       endif
    enddo
    if (k/=numlink) then
       print*,'k= ',k,'  numlink= ',numlink
       call die('link1: k/=numlink???')
    endif

    ! assignation of QM 1st neighbors
    rmin=0.0
    min=0
    do i=1,numlink
       m=0
       ! sp oxygen or sp2 carbon/nitrogen
       ch4=qmtype(linkqm(i,1))
       ch1=ch4(1:1)
       if (ch1=='O') then
         m=2
       else if(ch1=='C'.or.ch1=='N') then
         m=3
       endif
       ! sp3 carbon
       ch4=qmtype(linkqm(i,1))
       ch2=ch4(1:2)
       if(ch2=='CT') m=4
       if(m==0) then
          write(6,*) 'Wrong LA QM atom type....Check parameters'
          print*,'QM atom no. ',linkqm(i,1)
          STOP
       endif
       ! loop over CQM neighbors
       do k=2,m
          rmin(i)=10
          min(i)=0
          do j=1,na_qm
             if(j==linkqm(i,1).or.j==linkqm(i,2).or.j==linkqm(i,3)) then
                goto 5
             endif
             dx1=rqm(1,linkqm(i,1))-rqm(1,j)
             dy1=rqm(2,linkqm(i,1))-rqm(2,j)
             dz1=rqm(3,linkqm(i,1))-rqm(3,j)
             call pbc_displ_vector(lattice_type,cell,kcell,dx1,dy1,dz1)
             r1=dist_v2(dx1,dy1,dz1)    
             if(r1.le.rmin(i)) then
                rmin(i)=r1
                min(i)=j
             endif
5         enddo
          if(rmin(i).gt.4) then
             print*,'rmin(i)= ',j,rmin(i)
             write(6,*) 'Wrong LA QM neighbor....Check geometry'
             STOP 
          endif
          linkqm(i,k)=min(i)
       enddo
    enddo

    ! checking CQM neighbors
    do i=1,numlink
       ! sp2 carbon
       ch4=qmtype(linkqm(i,1))
       ch1=ch4(1:1)
       if (ch1=='O') then
          if(linkqm(i,2)==0.and.linkqm(i,3)==0.and.linkqm(i,4)==0) then
             write(6,*) 'Wrong QM neighbor number for a sp oxygen'   
             STOP
          endif
       else if(ch1=='C') then
          if(linkqm(i,2)==0.and.linkqm(i,3)==0.or.linkqm(i,2)==0.and. &
               linkqm(i,4)==0.or.linkqm(i,3)==0.and.linkqm(i,4)==0) then
             write(6,*) 'Wrong QM neighbor number for a sp2 carbon'   
             STOP
          endif
       endif
       ! sp3 carbon
       ch4=qmtype(linkqm(i,1))
       ch2=ch4(1:2)
       if(ch2=='CT') then
          if(linkqm(i,2)==0.or.linkqm(i,3)==0.or.linkqm(i,4)==0) then    
             write(6,*) 'Wrong QM neighbor number for a sp3 carbon' 
             STOP
          endif
       endif
    enddo

    ! asignation of QM link atoms types
    do i=1,numlink
       do j=1,4
          if(linkqm(i,j)==0) linkqmtype(i,j)='XX'
          linkqmtype(i,j)=qmtype(linkqm(i,j))
       enddo
    enddo

    ! assignation of MM 1st neighbors
    do i=1,numlink
       do k=2,4
          linkmm(i,k)=ng1(linkmm(i,1),k-1)
       enddo
    enddo

    ! assignation of MM 2nd neighbors
    do i=1,numlink
       do j=2,4
          if(linkmm(i,j)==0) goto 10
          m=1
          do k=1,4
             if(ng1(linkmm(i,j),k).ne.linkmm(i,1)) then
                linkmm2(i,j,m)=ng1(linkmm(i,j),k)
                m=m+1
             endif
          enddo
10     enddo
    enddo

    num_linkqm=0
    num_linkmm=0
    do i=1,numlink
       do j=1,4
          if (linkqm(i,j)/=0) num_linkqm(i)=num_linkqm(i)+1
          if (linkmm(i,j)/=0) num_linkmm(i)=num_linkmm(i)+1
       enddo
    enddo

    ! write LA params in file
    write(6,'(/,a)') 'siesta-qmmm: Link atom parameters:'
    if (h_link_scaling_type=='splam') then
       write(6,*)'  The link atom will be scaled by the modified version of the SPLAM method '
    else
       write(6,*)'  The link atom will be located at a fix distance from'//&
                     ' the QM atom'
    endif
    write(6,*)
    do i=1,numlink
       write(6,'(a,2x,1I6)') 'linkatom:',na_qm+na_mm+i
       write(6,'(a,2x,4I6)') 'qmatoms: ',(linkqm(i,j),j=1,num_linkqm(i))
       write(6,'(a,2x,4A4)') 'qmtypes: ',(linkqmtype(i,j),j=1,4)
       write(6,'(a,2x,4I6)') 'mmatoms: ',(linkmm(i,j),j=1,num_linkmm(i))
       write(6,'(a,2x,4A4)') 'mmtypes: ',(attype(linkmm(i,j)),j=1,num_linkmm(i))
       !      write(6,'(a,2x,9I6)') 'mmneihg: ',((linkmm2(i,j,k)+na,k=1,3),j=2,4)
    enddo
    write(6,*)

    return

 2     stop 'read: problem reading solute chemical specie label' 

  end subroutine link1

  subroutine get_link_ff_parameters(namber,nparm,nbond,nangle,ndihe,nimp,&
                   bondtype,angletype,dihetype,linkqmtype,attype,&
                         imptype,perdihe,kbond,bondeq,kangle,kdihe)

    use sys

    integer namber,nbond,nangle,ndihe,nimp,nparm
    character   bondtype(nparm)*5,angletype(nparm)*8,&
         dihetype(nparm)*11, imptype(nparm)*11, linkqmtype(numlink,4)*4
    character  attype(namber)*4
    real(dp) perdihe(nparm),perdihe2(nparm)
    real(dp) kbond(nparm),bondeq(nparm),kangle(nparm),kdihe(nparm)
    integer resonance_id, i, j, jm, jq, k, l, m
    character ty1*4,ty2*4,ty3*4,ty4*4,tybond*5,&
         tyangle*8,tydihe*11 
    character(len=4) :: ty1_h, ty2_h, ty3_h, tyl_h
    character(len=5) :: tybond_h
    character(len=8) :: tyangle_h
    integer :: parameter_set, parameter_bond_set
    real(dp) :: qmmm_range

    logical :: wriok

    wriok=.false.
    !wriok=.true.
    parametro=0
    link_bondeq=0.0
    link_kbond=0.0
    link_kangle=0.0

    do resonance_id=1,num_resonances

       ! asignation for E and F        
       do i=1,numlink

          parameter_set=0
          do k=1,nbond          !	bond Cqm -- Cmm parameter 1
             tybond=bondtype(k)
             ty1=tybond(1:2)
             ty2=tybond(4:5)
             if(linkqmtype(i,1)==ty1.and.attype(linkmm(i,1))==ty2) then
                parametro(resonance_id,i,1,1)=k
                parameter_set=1
             elseif(linkqmtype(i,1)==ty2.and.attype(linkmm(i,1))==ty1) then
                parametro(resonance_id,i,1,1)=k
                parameter_set=1
             endif
          enddo
          if (parameter_set==0) then
             print*,'Parameters for bond: ',linkqmtype(i,1),'-',attype(linkmm(i,1))
             call die('No suitable parameter for Cqm-Cmm bond')
          endif

          parameter_set=0  ! 	angles Cqm -- Cqm -- Cmm parameters 2 to 4    
          do j=2,4
             if(linkqm(i,j)==0) cycle
             do k=1,nangle
                tyangle=angletype(k)
                ty1=tyangle(1:2)
                ty2=tyangle(4:5)
                ty3=tyangle(7:8)
                if(linkqmtype(i,j)==ty1.and.linkqmtype(i,1)==ty2.and.&
                     attype(linkmm(i,1))==ty3) then
                   parametro(resonance_id,i,j,1)=k
                   parameter_set=1
                elseif(linkqmtype(i,j)==ty3.and.linkqmtype(i,1)==ty2&
                     .and.attype(linkmm(i,1))==ty1) then
                   parametro(resonance_id,i,j,1)=k
                   parameter_set=1
                endif
             enddo
          enddo
          if (parameter_set==0) then
             print*,'b) Parameters for angle: ',linkqmtype(i,j),'-',linkqmtype(i,1),'-',&
                  attype(linkmm(i,1))
             call die('No suitable parameter for Cqm-Cqm-Cm angle')
          endif

       enddo ! numlink

       do i=1,numlink

          tyl_h=resonance(resonance_id)%atype(i)

          if ('?'==tyl_h(2:2)) then

             parameter_bond_set=0
             bond_loop: do k=1,nbond         !	bond Cqm -- Cmm parameter 1

                parameter_set=0
                tybond_h=bondtype(k)
                ty1_h=tybond_h(1:2)
                ty2_h=tybond_h(4:5)
                
                if(linkqmtype(i,1)==ty1_h.and.tyl_h(1:1)==ty2_h(1:1)) then
                   link_bondeq(resonance_id,i)=bondeq(k)
                   link_kbond(resonance_id,i)=sqrt(kbond(parametro(resonance_id,i,1,1))/kbond(k))
                   parameter_set=1
                   parameter_bond_set=k
                else if (linkqmtype(i,1)==ty2_h.and.tyl_h(1:1)==ty1_h(1:1)) then
                   link_bondeq(resonance_id,i)=bondeq(k)
                   link_kbond(resonance_id,i)=sqrt(kbond(parametro(resonance_id,i,1,1))/kbond(k))
                   parameter_set=1
                   parameter_bond_set=k
                endif
                if (parameter_set==0) cycle bond_loop

                do j=2,4 

                   if(linkqm(i,j)==0) cycle

                   parameter_set=0
                   do l=1,nangle
                      tyangle_h=angletype(l)
                      ty1_h=tyangle_h(1:2)
                      ty2_h=tyangle_h(4:5)
                      ty3_h=tyangle_h(7:8)
                      tyl_h=resonance(resonance_id)%atype(i)
                      if(linkqmtype(i,j)==ty1_h.and.linkqmtype(i,1)==ty2_h.and.&
                           tyl_h(1:1)==ty3_h(1:1)) then
                         link_kangle(resonance_id,i,j-1)=kangle(l)
                         parameter_set=1
                         exit
                      elseif(linkqmtype(i,j)==ty3_h.and.linkqmtype(i,1)==ty2_h.and.&
                           tyl_h(1:1)==ty1_h(1:1)) then
                         link_kangle(resonance_id,i,j-1)=kangle(l)
                         parameter_set=1
                         exit
                      endif
                   enddo
                   if (parameter_set==0)cycle bond_loop

                enddo

                if (parameter_set==1) exit bond_loop

             enddo bond_loop

             if (parameter_set==0) then
                write(6,*)
                write(6,*)'WARNING'
                write(6,*)'In Cqm-Cmm bond: ',linkqmtype(i,1),'-',attype(linkmm(i,1))
                if (parameter_bond_set==0) then
                   write(6,*)'Default parameters set for Cqm-LA bond: ',linkqmtype(i,1),'-',resonance(resonance_id)%atype(i)        
                   call set_default_H_link_bond_parameters(linkqmtype(i,1),link_bondeq(resonance_id,i),&
                     kbond(parametro(resonance_id,i,1,1)),link_kbond(resonance_id,i))
                endif
                write(6,*)'Default parameters are set for angles:'                   
                do j=2,4
                   if(linkqm(i,j)==0) cycle
                   link_kangle(resonance_id,i,j-1)=kangle(parametro(resonance_id,i,j,1))
                   print*,'Angle: ',linkqmtype(i,j),'-',linkqmtype(i,1),'-',&
                        resonance(resonance_id)%atype(i)
                enddo
             endif

          else ! if ('?'/=tyl_h(2))

             parameter_set=0         
             do l=1,nbond
                tybond_h=bondtype(l)
                ty1_h=tybond_h(1:2)
                ty2_h=tybond_h(4:5)
                tyl_h=resonance(resonance_id)%atype(i)
                if(linkqmtype(i,1)==ty1_h.and.tyl_h(1:2)==ty2_h(1:2)) then
                   link_bondeq(resonance_id,i)=bondeq(l)
                   link_kbond(resonance_id,i)=sqrt(kbond(k)/kbond(l))
                   parameter_set=1
                   exit
                else if (linkqmtype(i,1)==ty2_h.and.tyl_h(1:2)==ty1_h(1:2)) then
                   link_bondeq(resonance_id,i)=bondeq(l)
                   link_kbond(resonance_id,i)=sqrt(kbond(k)/kbond(l))
                   parameter_set=1
                   exit
                endif
             enddo
             if (parameter_set==0) then
                print*,'Parameters for bond: ',linkqmtype(i,1),'-',resonance(resonance_id)%atype(i)
                call die('No suitable parameter for Cqm-LA bond')
             endif
       
             do j=2,4

                if(linkqm(i,j)==0) cycle

                parameter_set=0   
                do l=1,nangle  ! 	angles Cqm -- Cqm -- LA parameters 2 to 4 
                   tyangle_h=angletype(l)
                   ty1_h=tyangle_h(1:2)
                   ty2_h=tyangle_h(4:5)
                   ty3_h=tyangle_h(7:8)
                   tyl_h=resonance(resonance_id)%atype(i)
                   if(linkqmtype(i,j)==ty1_h.and.linkqmtype(i,1)==ty2_h.and.&
                        tyl_h(1:2)==ty3_h(1:2)) then
                      link_kangle(resonance_id,i,j-1)=kangle(l)
                      parameter_set=1
                      exit
                   elseif(linkqmtype(i,j)==ty3_h.and.linkqmtype(i,1)==ty2_h.and.&
                        tyl_h(1:2)==ty1_h(1:2)) then
                      link_kangle(resonance_id,i,j-1)=kangle(l)
                      parameter_set=1
                      exit
                   endif
                enddo
                if (parameter_set==0) then
                   print*,'Parameters for angle: ',linkqmtype(i,j),'-',linkqmtype(i,1),'-',&
                        resonance(resonance_id)%atype(i)
                   call die('No suitable parameter for Cqm-Cqm-LA angle')
                endif

             enddo  ! j=2,4

          endif  !  if ('?'==tyl_h(2)) 
 
       enddo ! numlink

       !  H  alyphatihic (1) 340.0 1.090 or aromathic (2) 367.0 1.080

       !	do i=1,numlink
       !
       !        if(linkqmtype(i,1)=='CT'.or.linkqmtype(i,1)=='CF') then
       !	kch(i)=340.0
       !        rch(i)=1.09	
       !	elseif(linkqmtype(i,1)=='CA'.or.linkqmtype(i,1)=='CB'.&
       !     	or.linkqmtype(i,1)=='CC') then
       !	kch(i)=367.0
       !	rch(i)=1.08
       !	else
       !	write(*,*) 'Wrong Hlink type  ',linkqmtype(i,1)
       !	STOP
       !	endif
       !	enddo

       ! 	angles Cqm -- Cmm -- X parameters 5 to 7

       do i=1,numlink
          do j=2,4
             if(linkmm(i,j)==0) cycle
             do k=1,nangle
                tyangle=angletype(k)
                ty1=tyangle(1:2)
                ty2=tyangle(4:5)
                ty3=tyangle(7:8)
                if(linkqmtype(i,1)==ty1.and.attype(linkmm(i,1))==ty2.and.&
                     attype(linkmm(i,j))==ty3) then
                   parametro(resonance_id,i,3+j,1)=k
                elseif(linkqmtype(i,1)==ty3.and.attype(linkmm(i,1))==ty2&
                     .and.attype(linkmm(i,j))==ty1) then
                   parametro(resonance_id,i,3+j,1)=k	
                endif
             enddo
          enddo
       enddo

       !	dihedral X--Cq -- Cmm--Y parameters 5 to 13

       do i=1,numlink
          do jq=2,4
             ! no more neighbours -> out 
             if(linkqm(i,jq)==0) goto 40	

             do jm=2,4
                ! only 2 mm neigh -> out
                if(linkmm(i,jm)==0) goto 39

                m=0
                do k=1,ndihe
                   tydihe=dihetype(k)
                   ty1=tydihe(1:2)
                   ty2=tydihe(4:5)
                   ty3=tydihe(7:8)
                   ty4=tydihe(10:11)
                   ! if X the same constant
                   if(ty1=='X ') then
                      if(linkqmtype(i,1)==ty2.and.&
                           attype(linkmm(i,1))==ty3)  then
                         parametro(resonance_id,i,7+3*(jq-2)+jm-1,1)=k
                      elseif(linkqmtype(i,1)==ty3.and.&
                           attype(linkmm(i,1))==ty2)  then
                         parametro(resonance_id,i,7+3*(jq-2)+jm-1,1)=k
                      endif
                      ! if not X different constants 
                   elseif(ty1.ne.'X ') then
                      if(linkqmtype(i,jq)==ty1.and.linkqmtype(i,1)==&
                           ty2.and.attype(linkmm(i,1))==ty3.and.&
                           attype(linkmm(i,jm))==ty4) then
                         parametro(resonance_id,i,7+3*(jq-2)+jm-1,1)=k

                         if(perdihe2(k).lt.0) then
                            m=m+1
                            parametro(resonance_id,i,7+3*(jq-2)+jm-1,m)=k           
                            parametro(resonance_id,i,7+3*(jq-2)+jm-1,m+1)=k+1   
                         endif
                      elseif(linkqmtype(i,jq)==ty4.and.linkqmtype(i,1)==&
                           ty3.and.attype(linkmm(i,1))==ty2.and.&
                           attype(linkmm(i,jm))==ty1) then
                         parametro(resonance_id,i,7+3*(jq-2)+jm-1,1)=k         

                         if(perdihe2(k).lt.0) then
                            m=m+1
                            parametro(resonance_id,i,7+3*(jq-2)+jm-1,m)=k
                            parametro(resonance_id,i,7+3*(jq-2)+jm-1,m+1)=k+1 
                         endif
                      endif
                   endif
                enddo
39           enddo

40        enddo
       enddo

       !	dihedral Cq--Cmm--Y--Z parameters 14 to 22

       do i=1,numlink
          do jm=2,4
             if(linkmm(i,jm)==0) goto 51
             do j=1,3
                if(linkmm2(i,jm,j)==0) then
                   goto 50
                else
                   m=0
                   do k=1,ndihe
                      tydihe=dihetype(k)
                      ty1=tydihe(1:2)
                      ty2=tydihe(4:5)
                      ty3=tydihe(7:8)
                      ty4=tydihe(10:11)
                      ! if X the same constant     
                      if(ty1=='X ') then
                         if(attype(linkmm(i,1))==ty2.and. &
                              attype(linkmm(i,jm))==ty3)  then
                            parametro(resonance_id,i,13+3*(jm-2)+j,1)=k
                         elseif(attype(linkmm(i,1))==ty3.and. &
                              attype(linkmm(i,jm))==ty2)  then
                            parametro(resonance_id,i,13+3*(jm-2)+j,1)=k
                         endif
                         ! if not X different constants     
                      elseif(ty1.ne.'X ') then

                         if(linkqmtype(i,1)==ty1.and.attype(linkmm(i,1))==&
                              ty2.and.attype(linkmm(i,jm))==ty3.and.&
                              attype(linkmm2(i,jm,j))==ty4) then

                            if(perdihe2(k).lt.0) then
                               parametro(resonance_id,i,13+3*(jm-2)+j,1)=k
                               if(perdihe2(k+1).lt.0) then
                                  parametro(resonance_id,i,13+3*(jm-2)+j,2)=k+1
                                  if(perdihe2(k+2).lt.0) then
                                     parametro(resonance_id,i,13+3*(jm-2)+j,3)=k+2
                                     parametro(resonance_id,i,13+3*(jm-2)+j,4)=k+3
                                  else
                                     parametro(resonance_id,i,13+3*(jm-2)+j,3)=k+2
                                     goto 50
                                  endif
                               else
                                  parametro(resonance_id,i,13+3*(jm-2)+j,2)=k+1
                                  goto 50
                               endif
                            else	  
                               parametro(resonance_id,i,13+3*(jm-2)+j,1)=k
                               goto 50
                            endif

                         elseif(linkqmtype(i,1)==ty4.and.attype(linkmm(i,1))==&
                              ty3.and.attype(linkmm(i,jm))==ty2.and.&
                              attype(linkmm2(i,jm,j))==ty1) then 

                            if(perdihe2(k).lt.0) then
                               parametro(resonance_id,i,13+3*(jm-2)+j,1)=k
                               if(perdihe2(k+1).lt.0) then
                                  parametro(resonance_id,i,13+3*(jm-2)+j,2)=k+1
                                  if(perdihe2(k+2).lt.0) then
                                     parametro(resonance_id,i,13+3*(jm-2)+j,3)=k+2
                                     parametro(resonance_id,i,13+3*(jm-2)+j,4)=k+3
                                  else
                                     parametro(resonance_id,i,13+3*(jm-2)+j,3)=k+2
                                     goto 50
                                  endif
                               else
                                  parametro(resonance_id,i,13+3*(jm-2)+j,2)=k+1
                                  goto 50
                               endif
                            else
                               parametro(resonance_id,i,13+3*(jm-2)+j,1)=k
                               goto 50	      	 
                            endif

                            ! one side or other
                         endif
                         ! with or without X
                      endif
                      ! end dihedrals 
                   enddo
                endif
                ! enddo for neigh 2, neigh 1, numlink 
50           enddo
51        enddo
          ! end numlink
       enddo

       ! end asignation

       !	some printting for debugging	
       if(wriok) then

          write(*,*) 'LINK ATOM VARIABLES'

          ! neigh matrix

          write(*,*)''
          do i=1,numlink
             write(*,*) 'Neigh Matrix: ',i
             write(*,*) 'Cq:',linkqm(i,1),'Vec: ',(linkqm(i,k),k=2,4)
             write(*,*) 'Cm:	  	    ',linkmm(i,1)
             write(*,*) 'Neigh CM:   ',(linkmm(i,k),k=2,4)
             write(*,*) 'Neigh CM 2: ',((linkmm2(i,j,k),k=1,3),j=2,4)	
             write(*,*)''


             ! cq-cm

             write(*,*) 'bond Cq-Cm: ',linkqm(i,1),linkmm(i,1)
             write(*,*) 'types: ',linkqmtype(i,1),attype(linkmm(i,1))
             write(*,*) 'Bondtypenumber(k) ',parametro(resonance_id,i,1,1),&
                  kbond(parametro(resonance_id,i,1,1))
             write(*,*)'************************************'

             ! X -- Cqm -- Cmm parameters 2 to 4

             write(*,*) 'Angle X-Cq-Cm: ',linkqm(i,1),linkmm(i,1)
             write(*,*)''
             write(*,*) ' x: ',linkqm(i,2),linkqm(i,3),linkqm(i,4)
             write(*,*) 'x types: ',linkqmtype(i,2),&
                  linkqmtype(i,3),linkqmtype(i,4)
             write(*,'(a12,i4,f9.3,i4,f9.3,i4,f9.3)')&
                  'Angtype(k): ',&
                  parametro(resonance_id,i,2,1),kangle(parametro(resonance_id,i,2,1)),&
                  parametro(resonance_id,i,3,1),kangle(parametro(resonance_id,i,3,1)),&
                  parametro(resonance_id,i,4,1),kangle(parametro(resonance_id,i,4,1))
             write(*,*)'***********************************'                
!!$
             ! Cqm -- Cmm -- X parameters 5 to 7

             write(*,*) 'Angle Cq-Cm-x: ',linkqm(i,1),linkmm(i,1)
             write(*,*)''
             write(*,*) ' x: ',linkmm(i,2),linkmm(i,3),linkmm(i,4)
             write(*,*) 'x types: ',attype(linkmm(i,2)),&
                  attype(linkmm(i,3)),attype(linkmm(i,4))
             write(*,'(a12,i4,f9.3,i4,f9.3,i4,f9.3)')&
                  'Angtype(k): ',&
                  parametro(resonance_id,i,5,1),kangle(parametro(resonance_id,i,4,1)),&
                  parametro(resonance_id,i,6,1),kangle(parametro(resonance_id,i,6,1)),&
                  parametro(resonance_id,i,7,1),kangle(parametro(resonance_id,i,7,1))
             write(*,*)'***********************************'

             ! X--Cq -- Cmm--Y parameters 8 to 16     

             write(*,*)''
             write(*,*) 'Dihedrals y-Cq-Cm-x: ',linkqm(i,1),linkmm(i,1)
             write(*,*)'' 
             write(*,*) ' x: ',linkqm(i,2),linkqm(i,3),linkqm(i,4)
             write(*,*) 'x types: ',linkqmtype(i,2),&
                  linkqmtype(i,3),linkqmtype(i,4)	
             write(*,*)''
             write(*,*) ' y: ',linkmm(i,2),linkmm(i,3),linkmm(i,4)
             write(*,*) 'y types: ',attype(linkmm(i,2)),&
                  attype(linkmm(i,3)),attype(linkmm(i,4))
             write(*,*)''
             ! x=1
             write(*,'(a20,i2,f6.3,i2,f6.3,i2,f6.3,i2,f6.3,i2,f6.3)')&
                  'dihe x=1,y=1: (mult) ',&
                  (parametro(resonance_id,i,8,k),kdihe(parametro(resonance_id,i,8,k)),k=1,4)
             write(*,'(a20,i2,f6.3,i2,f6.3,i2,f6.3,i2,f6.3,i2,f6.3)')&
                  'dihe x=1,y=2: (mult) ',&
                  (parametro(resonance_id,i,9,k),kdihe(parametro(resonance_id,i,9,k)),k=1,4)
             write(*,'(a20,i2,f6.3,i2,f6.3,i2,f6.3,i2,f6.3,i2,f6.3)')&
                  'dihe x=1,y=3: (mult) ',&
                  (parametro(resonance_id,i,10,k),kdihe(parametro(resonance_id,i,10,k)),k=1,4)
             write(*,*)''
             ! x=2
             write(*,'(a20,i2,f6.3,i2,f6.3,i2,f6.3,i2,f6.3,i2,f6.3)')&
                  'dihe x=2,y=1: (mult) ',&
                  (parametro(resonance_id,i,11,k),kdihe(parametro(resonance_id,i,11,k)),k=1,4)
             write(*,'(a20,i2,f6.3,i2,f6.3,i2,f6.3,i2,f6.3,i2,f6.3)')&
                  'dihe x=2,y=2: (mult) ',&
                  (parametro(resonance_id,i,12,k),kdihe(parametro(resonance_id,i,12,k)),k=1,4)
             write(*,'(a20,i2,f6.3,i2,f6.3,i2,f6.3,i2,f6.3,i2,f6.3)')&
                  'dihe x=2,y=3: (mult) ',&
                  (parametro(resonance_id,i,13,k),kdihe(parametro(resonance_id,i,13,k)),k=1,4)
             write(*,*)''
             ! x=3
             write(*,'(a20,i2,f6.3,i2,f6.3,i2,f6.3,i2,f6.3,i2,f6.3)')&
                  'dihe x=3,y=1: (mult) ',&
                  (parametro(resonance_id,i,14,k),kdihe(parametro(resonance_id,i,14,k)),k=1,4)
             write(*,'(a20,i2,f6.3,i2,f6.3,i2,f6.3,i2,f6.3,i2,f6.3)')&
                  'dihe x=3,y=2: (mult) ',&
                  (parametro(resonance_id,i,15,k),kdihe(parametro(resonance_id,i,15,k)),k=1,4)
             write(*,'(a20,i2,f6.3,i2,f6.3,i2,f6.3,i2,f6.3,i2,f6.3)')&
                  'dihe x=3,y=3: (mult) ',&
                  (parametro(resonance_id,i,16,k),kdihe(parametro(resonance_id,i,16,k)),k=1,4)
             write(*,*)'*********************************************'

             ! Cq--Cmm--Y--Z parameters 17 to 25 

             write(*,*) 'Dihedrals Cq-Cm-y-z: ',linkqm(i,1),linkmm(i,1)
             write(*,*)''
             write(*,*) ' y: ',linkmm(i,2),linkmm(i,3),linkmm(i,4)
             write(*,*) 'y types: ',attype(linkmm(i,2)),&
                  attype(linkmm(i,3)),attype(linkmm(i,4))

             write(*,*)''

             do k=1,3
                if(linkmm2(i,2,k)==0) then
                   goto 60
                else
                   l=16+k
                   write(*,*) ' z:(y1) ',linkmm2(i,2,k)
                   write(*,*) 'z(y1) types: ',attype(linkmm2(i,2,k))
                   write(*,'(a20,i2,f6.3,i2,f6.3,i2,f6.3,i2,f6.3,i2,f6.3)')&
                        'dihedro : (mult)   ',&
                        (parametro(resonance_id,i,l,m),kdihe(parametro(resonance_id,i,l,m)),m=1,4)
                endif
60           enddo

             do k=1,3
                if(linkmm2(i,3,k)==0) then
                   goto 61
                else
                   l=19+k

                   write(*,*) ' z:(y2) ',linkmm2(i,3,k)
                   write(*,*) 'z(y2) types: ',attype(linkmm2(i,3,k))
                   write(*,'(a20,i2,f6.3,i2,f6.3,i2,f6.3,i2,f6.3,i2,f6.3)')&
                        'dihedro : (mult)   ',&
                        (parametro(resonance_id,i,l,m),kdihe(parametro(resonance_id,i,l,m)),m=1,4)
                endif
61           enddo

             do k=1,3
                if(linkmm2(i,4,k)==0) then
                   goto 62
                else
                   l=22+k

                   write(*,*) ' z:(y3) ',linkmm2(i,4,k)
                   write(*,*) 'z(y3) types: ',attype(linkmm2(i,4,k))
                   write(*,'(a20,i2,f6.3,i2,f6.3,i2,f6.3,i2,f6.3,i2,f6.3)')&
                        'dihedro : (mult)   ',&
                        (parametro(resonance_id,i,l,m),kdihe(parametro(resonance_id,i,l,m)),m=1,4)

                endif
62           enddo
             write(*,*) '*********************************************'
             
             ! numlink
          enddo
          write(*,*) 'parametro matrix'
          do i=1,numlink
             do j=1,25
                write(*,*) i,j,parametro(resonance_id,i,j,1:4)
             enddo
          enddo
          ! wriok
       endif

    enddo  ! resonance_id

  end subroutine get_link_ff_parameters

  !*************************************************************
  ! calculates Energy and Forces of HLink s
  ! ramber=rclas
  subroutine link2(ramber,natot,namber,fdummy,stress,ng1,nparm,&
       ndihe,multidihe,multiimp,kbond,bondeq,&
       kangle,angleeq,kdihe,diheeq,kimp,impeq,perdihe,perimp,&
       bondxat,Ener,istp,kch,rch,lattice_type,cell)

    integer natot,namber,nparm,ndihe,ng1(namber,6),bondxat(namber)
    integer i,j,j2,k,l,m,jm,jq,at1,at2,at3,at4,istp,p
    real(dp) Elink_for_h_atom,ramber(3,natot),flink(3,natot),&
         fdummy(3,natot),fce(12),scaling
    real(dp) :: stress(3,3)
    integer multidihe(nparm),multiimp(nparm)
    real(dp) kbond(nparm),bondeq(nparm),kangle(nparm),&
         angleeq(nparm),kdihe(nparm),diheeq(nparm),kimp(nparm),&
         impeq(nparm),perdihe(nparm),perdihe2(nparm),perimp(nparm)
    real(dp)  rij,dist,pi,angle,&
         scal,scalar,dscalar,r12,r32,dr12r32,angulo,&
         dihedro_v2,dihedral,fpp(3),fmod,ip,&
         dversor(3),rhq,rqm,Ener,kch(numlink),rch(numlink),distl(numlink)
    real(dp), dimension(3) :: dr
    real(dp) :: dx, dy, dz, volcel
    real(dp) dx12, dy12, dz12, dx32, dy32, dz32

    integer :: resonance_id
    real(dp) :: angle_v2, scalar_v2
    real :: stress_fact

    ! numerical variables to indicate the parameters
    ! last value eq 0, 1 dihedral, eq 1 more than 1 dihedral
    character lattice_type
    real(dp) cell(3,3), amber_cell(3,3), amber_kcell(3,3)
    real(dp) h_kangle, dist_v2

    real(dp) kfactor

    ! linkmm2 asigantion: linkmm atoms neighbours
    perdihe2=perdihe
    pi=DACOS(-1.d0)

    ! change units 
    ramber(1:3,1:natot)=ramber(1:3,1:natot)*0.529177
    rlink(1:num_resonances,1:3,1:numlink)=rlink(1:num_resonances,1:3,1:numlink)*0.529177

    amber_cell(1:3,1:3)=cell(1:3,1:3)*0.529177
    call reccel(3,amber_cell,amber_kcell,0)

    ! reasignation of perdihe2
    do i=1,ndihe
       if(perdihe2(i).lt.0) then
          perdihe2(i)=-1.0*perdihe2(i)
       endif
    enddo

    stress_fact=1.0_dp/volcel(amber_cell)

    do resonance_id=1,num_resonances

       Ener=0.0

       flink(1:3,1:natot)=0.0

       ! calculation of E and F 

       do i=1,numlink

          Elink_for_h_atom=0.0

          if (h_link_scaling_type/='splam') then

             !       correction to bond Cqm -- X parameter 1 

             at1=linkqm(i,1)
             at2=natot-namber+linkmm(i,1)
             k=parametro(resonance_id,i,1,1)

             dr(1)=ramber(1,at1)-ramber(1,at2)
             dr(2)=ramber(2,at1)-ramber(2,at2)
             dr(3)=ramber(3,at1)-ramber(3,at2)
             call pbc_displ_vector(lattice_type,amber_cell,amber_kcell,dr(1),dr(2),dr(3))

             rij=dist_v2(dr(1),dr(2),dr(3))

             Elink_for_h_atom = Elink_for_h_atom + kbond(k)*(rij-bondeq(k))**2 

             do l=1,3
                dr(l)=(1.0/rij)*dr(l)
                dr(l)=2.0*kbond(k)*(rij-bondeq(k))*dr(l)
                flink(l,at1)=-dr(l)
                flink(l,at2)=dr(l)
                do m=1,3
                   stress(m,l)=stress(m,l)+stress_fact*ramber(m,at1)*dr(l)
                enddo
             enddo

          endif

          ! ***********************************************************************
          !       angle  X -- Cqm -- Cmm parameters 2 to 4

          !	3 angles ( one for each x) j, x index
          dx=0.0
          dy=0.0
          dz=0.0

          do j=2,4

             if(linkqm(i,j)==0) cycle

             at1=natot-namber+linkmm(i,1)
             at2=linkqm(i,1)
             at3=linkqm(i,j)
             k=parametro(resonance_id,i,j,1)

             dx12=ramber(1,at1)-ramber(1,at2)
             dy12=ramber(2,at1)-ramber(2,at2)
             dz12=ramber(3,at1)-ramber(3,at2)
             dx32=ramber(1,at3)-ramber(1,at2)
             dy32=ramber(2,at3)-ramber(2,at2)
             dz32=ramber(3,at3)-ramber(3,at2)
             call pbc_displ_vector(lattice_type,amber_cell,amber_kcell,dx12,dy12,dz12)
             call pbc_displ_vector(lattice_type,amber_cell,amber_kcell,dx32,dy32,dz32)
             angulo=angle_v2(dx12,dy12,dz12,dx32,dy32,dz32)

             h_kangle=kangle(k)-link_kangle(resonance_id,i,j-1)

             Elink_for_h_atom = Elink_for_h_atom + h_kangle*&
                               ((angulo-angleeq(k))*(pi/180))**2

             !	write(*,*) 'contrib cq-cm-x(j) i,j E ',i,j,kangle(k)*&
             !                      ((angulo-angleeq(k))*&
             !                      (pi/180))**2

             scal=scalar_v2(dx12,dy12,dz12,dx32,dy32,dz32)

             r12=dist_v2(dx12,dy12,dz12)
             r32=dist_v2(dx32,dy32,dz32)

             dr12r32=r32*dx12/(r12)
             dx=(dx32*r12*r32-scal*dr12r32)/(r12*r32)**(2.0)
             dx=-1.0/(sqrt(1.0-(scal/(r12*r32))**2.0))*dx
             dx=2.0*h_kangle*(angulo-angleeq(k))*(pi/180)*dx

             dr12r32=r32*dy12/(r12)
             dy=(dy32*r12*r32-scal*dr12r32)/(r12*r32)**(2.0)
             dy=-1.0/(sqrt(1.0-(scal/(r12*r32))**2.0))*dy
             dy=2.0*h_kangle*(angulo-angleeq(k))*(pi/180)*dy

             dr12r32=r32*dz12/(r12)
             dz=(dz32*r12*r32-scal*dr12r32)/(r12*r32)**(2.0)
             dz=-1.0/(sqrt(1.0-(scal/(r12*r32))**2.0))*dz
             dz=2.0*h_kangle*(angulo-angleeq(k))*(pi/180)*dz

             ! atom qm cq force, sum of each particular x
             flink(1,at1)=flink(1,at1)-dx
             flink(2,at1)=flink(2,at1)-dy
             flink(3,at1)=flink(3,at1)-dz

             ! force over x, change at3 by at1 
             dx=0.0
             dy=0.0
             dz=0.0
             dr12r32=r12*dx32/(r32)
             dx=(dx12*r12*r32-scal*dr12r32)/(r12*r32)**(2.0)
             dx=-1.0/(sqrt(1.0-(scal/(r12*r32))**2.0))*dx
             dx=2.0*h_kangle*(angulo-angleeq(k))*(pi/180)*dx

             dr12r32=r12*dy32/(r32)
             dy=(dy12*r12*r32-scal*dr12r32)/(r12*r32)**(2.0)
             dy=-1.0/(sqrt(1.0-(scal/(r12*r32))**2.0))*dy
             dy=2.0*h_kangle*(angulo-angleeq(k))*(pi/180)*dy

             dr12r32=r12*dz32/(r32)
             dz=(dz12*r12*r32-scal*dr12r32)/(r12*r32)**(2.0)
             dz=-1.0/(sqrt(1.0-(scal/(r12*r32))**2.0))*dz
             dz=2.0*h_kangle*(angulo-angleeq(k))*(pi/180)*dz

             ! sum each x  
             flink(1,at3)=flink(1,at3)-dx
             flink(2,at3)=flink(2,at3)-dy
             flink(3,at3)=flink(3,at3)-dz

             ! middle atom force (change in derivate) 
             dx=0.0
             dy=0.0
             dz=0.0

             dscalar=-dx12-dx32
             dr12r32=(r32*(-dx12)/r12)+(r12*(-dx32)/r32)
             dx=(dscalar*r12*r32-scal*dr12r32)/(r12*r32)**(2.0)
             dx=-1.0/(sqrt(1.0-(scal/(r12*r32))**2.0))*dx
             dx=2.0*h_kangle*(angulo-angleeq(k))*(pi/180)*dx

             dscalar=-dy12-dy32
             dr12r32=(r32*(-dy12)/r12)+(r12*(-dy32)/r32)
             dy=(dscalar*r12*r32-scal*dr12r32)/(r12*r32)**(2.0)
             dy=-1.0/(sqrt(1.0-(scal/(r12*r32))**2.0))*dy
             dy=2.0*h_kangle*(angulo-angleeq(k))*(pi/180)*dy

             dscalar=-dz12-dz32
             dr12r32=(r32*(-dz12)/r12)+(r12*(-dz32)/r32)
             dz=(dscalar*r12*r32-scal*dr12r32)/(r12*r32)**(2.0)
             dz=-1.0/(sqrt(1.0-(scal/(r12*r32))**2.0))*dz
             dz=2.0*h_kangle*(angulo-angleeq(k))*(pi/180)*dz

             ! central atom force 
             flink(1,at2)=flink(1,at2)-dx
             flink(2,at2)=flink(2,at2)-dy
             flink(3,at2)=flink(3,at2)-dz

             ! enddo for each possible x 
          enddo


          !	write(*,*) 'contrib cq-cm i E F',i,&
          !        kbond(k)*(1.0-kbond(k)/kch(i))&
          !          *(rij-bondeq(k))**2,&
          !         sqrt(flink(1,at1)**2+flink(2,at1)**2+flink(3,at1)**2)

          ! ***********************************************************************
          !       angle  Cqm -- Cmm -- X parameters 5 to 7

          !	3 angles ( one for each x) j, x index
          dx=0.0
          dy=0.0
          dz=0.0

          do j=2,4

             if(linkmm(i,j)==0) cycle

             at1=linkqm(i,1)
             at2=natot-namber+linkmm(i,1)
             at3=natot-namber+linkmm(i,j)
             k=parametro(resonance_id,i,3+j,1)

             dx12=ramber(1,at1)-ramber(1,at2)
             dy12=ramber(2,at1)-ramber(2,at2)
             dz12=ramber(3,at1)-ramber(3,at2)
             dx32=ramber(1,at3)-ramber(1,at2)
             dy32=ramber(2,at3)-ramber(2,at2)
             dz32=ramber(3,at3)-ramber(3,at2)
             call pbc_displ_vector(lattice_type,amber_cell,amber_kcell,dx12,dy12,dz12)
             call pbc_displ_vector(lattice_type,amber_cell,amber_kcell,dx32,dy32,dz32)
             angulo=angle_v2(dx12,dy12,dz12,dx32,dy32,dz32)

             Elink_for_h_atom = Elink_for_h_atom + kangle(k)*((angulo-angleeq(k))*&
                  (pi/180))**2

             !	write(*,*) 'contrib cq-cm-x(j) i,j E ',i,j,kangle(k)*&
             !                      ((angulo-angleeq(k))*&
             !                      (pi/180))**2

             scal=scalar_v2(dx12,dy12,dz12,dx32,dy32,dz32)

             r12=dist_v2(dx12,dy12,dz12)
             r32=dist_v2(dx32,dy32,dz32)

             dr12r32=r32*dx12/(r12)
             dx=(dx32*r12*r32-scal*dr12r32)/(r12*r32)**(2.0)
             dx=-1.0/(sqrt(1.0-(scal/(r12*r32))**2.0))*dx
             dx=2.0*kangle(k)*(angulo-angleeq(k))*(pi/180)*dx

             dr12r32=r32*dy12/(r12)
             dy=(dy32*r12*r32-scal*dr12r32)/(r12*r32)**(2.0)
             dy=-1.0/(sqrt(1.0-(scal/(r12*r32))**2.0))*dy
             dy=2.0*kangle(k)*(angulo-angleeq(k))*(pi/180)*dy

             dr12r32=r32*dz12/(r12)
             dz=(dz32*r12*r32-scal*dr12r32)/(r12*r32)**(2.0)
             dz=-1.0/(sqrt(1.0-(scal/(r12*r32))**2.0))*dz
             dz=2.0*kangle(k)*(angulo-angleeq(k))*(pi/180)*dz

             ! atom qm cq force, sum of each particular x
             flink(1,at1)=flink(1,at1)-dx
             flink(2,at1)=flink(2,at1)-dy
             flink(3,at1)=flink(3,at1)-dz

             ! force over x, change at3 by at1 
             dx=0.0
             dy=0.0
             dz=0.0
             dr12r32=r32*dx32/(r12)
             dx=(dx12*r12*r32-scal*dr12r32)/(r12*r32)**(2.0)
             dx=-1.0/(sqrt(1.0-(scal/(r12*r32))**2.0))*dx
             dx=2.0*kangle(k)*(angulo-angleeq(k))*(pi/180)*dx

             dr12r32=r32*dy32/(r12)
             dy=(dy12*r12*r32-scal*dr12r32)/(r12*r32)**(2.0)
             dy=-1.0/(sqrt(1.0-(scal/(r12*r32))**2.0))*dy
             dy=2.0*kangle(k)*(angulo-angleeq(k))*(pi/180)*dy

             dr12r32=r32*dz32/(r12)
             dz=(dz12*r12*r32-scal*dr12r32)/(r12*r32)**(2.0)
             dz=-1.0/(sqrt(1.0-(scal/(r12*r32))**2.0))*dz
             dz=2.0*kangle(k)*(angulo-angleeq(k))*(pi/180)*dz

             ! sum each x  
             flink(1,at3)=flink(1,at3)-dx
             flink(2,at3)=flink(2,at3)-dy
             flink(3,at3)=flink(3,at3)-dz

             ! middle atom force (change in derivate) 
             dx=0.0
             dy=0.0
             dz=0.0

             dscalar=-dx12-dx32
             dr12r32=(r32*(-dx12)/r12)+(r12*(-dx32)/r32)
             dx=(dscalar*r12*r32-scal*dr12r32)/(r12*r32)**(2.0)
             dx=-1.0/(sqrt(1.0-(scal/(r12*r32))**2.0))*dx
             dx=2.0*kangle(k)*(angulo-angleeq(k))*(pi/180)*dx

             dscalar=-dy12-dy32
             dr12r32=(r32*(-dy12)/r12)+(r12*(-dy32)/r32)
             dy=(dscalar*r12*r32-scal*dr12r32)/(r12*r32)**(2.0)
             dy=-1.0/(sqrt(1.0-(scal/(r12*r32))**2.0))*dy
             dy=2.0*kangle(k)*(angulo-angleeq(k))*(pi/180)*dy

             dscalar=-dz12-dz32
             dr12r32=(r32*(-dz12)/r12)+(r12*(-dz32)/r32)
             dz=(dscalar*r12*r32-scal*dr12r32)/(r12*r32)**(2.0)
             dz=-1.0/(sqrt(1.0-(scal/(r12*r32))**2.0))*dz
             dz=2.0*kangle(k)*(angulo-angleeq(k))*(pi/180)*dz

             ! central atom force 
             flink(1,at2)=flink(1,at2)-dx
             flink(2,at2)=flink(2,at2)-dy
             flink(3,at2)=flink(3,at2)-dz

             ! enddo for each possible x 
          enddo

          !********************************************************
          !	dihedrals  X--Cq -- Cmm--Y parameters 7 to 16
          !	for each 3 X 
          !	for each 3 Y
          do jq=2,4
             do jm=2,4

                if(linkmm(i,jm)==0) goto 71
                at1=linkqm(i,jq)
                at2=linkqm(i,1)   
                at3=natot-namber+linkmm(i,1)
                at4=natot-namber+linkmm(i,jm)      

                dihedral=dihedro_v2(ramber(1,at1),ramber(2,at1),ramber(3,at1),&
                     ramber(1,at2),ramber(2,at2),ramber(3,at2),&
                     ramber(1,at3),ramber(2,at3),ramber(3,at3),&
                     ramber(1,at4),ramber(2,at4),ramber(3,at4),amber_cell,amber_kcell,&
                     lattice_type)

                ! who many parameters for dihedrals 7 to 16
                do m=1,4		
                   k=parametro(resonance_id,i,7+3*(jq-2)+jm-1,m)
                   if(k.ne.0) then	

                      Elink_for_h_atom =Elink_for_h_atom+(kdihe(k)/multidihe(k))*&
                           (1+dCOS((pi/180)*(perdihe2(k)*dihedral-diheeq(k))))


                      !	write(*,*) 'contrib x-cq-cm-y ',i,jq,jm,(kdihe(k)/multidihe(k))*&
                      !       (1+dCOS((pi/180)*(perdihe2(k)*dihedral-diheeq(k))))

                      ! force for each 4 atoms 
                      do j=1,4
                         call diheforce_v2(namber,ramber,at1,at2,at3,at4,j,&
                              kdihe(k),diheeq(k),perdihe2(k),multidihe(k),&
                              fce,amber_cell,amber_kcell,lattice_type)

                         ! force assignation 
                         flink(1,at1)=flink(1,at1)-fce(1)
                         flink(2,at1)=flink(2,at1)-fce(2)
                         flink(3,at1)=flink(3,at1)-fce(3) 
                         flink(1,at2)=flink(1,at2)-fce(4)
                         flink(2,at2)=flink(2,at2)-fce(5)
                         flink(3,at2)=flink(3,at2)-fce(6)
                         flink(1,at3)=flink(1,at3)-fce(7)
                         flink(2,at3)=flink(2,at3)-fce(8)
                         flink(3,at3)=flink(3,at3)-fce(9)
                         flink(1,at4)=flink(1,at4)-fce(10)
                         flink(2,at4)=flink(2,at4)-fce(11)
                         flink(3,at4)=flink(3,at4)-fce(12)	

                      enddo
                   else
                      goto 70
                   endif
                   ! enddo more than 1 parameter 
70              enddo

                ! enddo X and Y 
71           enddo
          enddo

          ! 	dihedrals Cq -- Cmm--Y--Z  parameters 17 to 25
          ! 	for each 3 X
          ! 	for each 3 Y
          
          do jm=2,4
             if(linkmm(i,jm)==0) goto 81
             do jq=1,3 
                ! is other Z atom? 
                if(linkmm2(i,jm,jq)==0) then
                   goto 80
                else
                   l=16+3*(jm-2)+jq
                   at1=linkqm(i,1)
                   at2=natot-namber+linkmm(i,1)
                   at3=natot-namber+linkmm(i,jm)
                   at4=natot-namber+linkmm2(i,jm,jq) 

                   dihedral=dihedro_v2(ramber(1,at1),ramber(2,at1),ramber(3,at1),&
                        ramber(1,at2),ramber(2,at2),ramber(3,at2),ramber(1,at3),&
                        ramber(2,at3),ramber(3,at3),ramber(1,at4),ramber(2,at4),&
                        ramber(3,at4),amber_cell,amber_kcell,lattice_type)


                   ! how many parameters for dihedrals 17 to 25 
                   ! more than 1 paramter
                   do m=1,4
                      k=parametro(resonance_id,i,l,m)
                      if(k.ne.0) then

                         Elink_for_h_atom =Elink_for_h_atom+(kdihe(k)/multidihe(k))*&
                              (1+dCOS((pi/180)*(perdihe2(k)*dihedral-diheeq(k))))

                         ! force over 4 atoms
                         do j=1,4
                            call diheforce_v2(namber,ramber,at1,at2,at3,at4,j,&
                                 kdihe(k),diheeq(k),perdihe2(k),multidihe(k),fce,&
                                    amber_cell,amber_kcell,lattice_type)

                            ! force assignation 
                            flink(1,at1)=flink(1,at1)-fce(1)
                            flink(2,at1)=flink(2,at1)-fce(2)
                            flink(3,at1)=flink(3,at1)-fce(3)
                            flink(1,at2)=flink(1,at2)-fce(4)
                            flink(2,at2)=flink(2,at2)-fce(5)
                            flink(3,at2)=flink(3,at2)-fce(6)
                            flink(1,at3)=flink(1,at3)-fce(7)
                            flink(2,at3)=flink(2,at3)-fce(8)
                            flink(3,at3)=flink(3,at3)-fce(9)
                            flink(1,at4)=flink(1,at4)-fce(10)
                            flink(2,at4)=flink(2,at4)-fce(11)
                            flink(3,at4)=flink(3,at4)-fce(12)
                         enddo

                      else
                         goto 75
                      endif
                      ! enddo more than 1 parameter 
75                 enddo

                endif
                ! enddo to Y and Z
80           enddo
81        enddo

          ! ***************************************************
          ! For each LA Force over HL scaling and sum where corresponds 
          ! HL force division in parallel fpar and perpendicular fpp
          ! unitary versor dversor(3) at3=LA
          !scaling=0.1
          at1=linkqm(i,1)
          at2=natot-namber+linkmm(i,1)

          dversor(1)=ramber(1,at1)-rlink(resonance_id,1,i)
          dversor(2)=ramber(2,at1)-rlink(resonance_id,2,i)
          dversor(3)=ramber(3,at1)-rlink(resonance_id,3,i)

          call pbc_displ_vector(lattice_type,amber_cell,amber_kcell,dversor(1),&
               dversor(2),dversor(3))

          rhq=dist_v2(dversor(1),dversor(2),dversor(3))

          dversor(1:3)=dversor(1:3)/rhq

          ip=dversor(1)*fa_link(resonance_id,1,i)+&
                dversor(2)*fa_link(resonance_id,2,i)+&
                  dversor(3)*fa_link(resonance_id,3,i)

          if (h_link_scaling_type=='splam') then
             flink(1:3,at2)=flink(1:3,at2)+link_kbond(resonance_id,i)*ip*dversor(1:3)
             flink(1:3,at1)=flink(1:3,at1)+(1-link_kbond(resonance_id,i))*ip*dversor(1:3)
          else
             flink(1:3,at1)=flink(1:3,at1)+ip*dversor(1:3)
          endif
          fpp(1:3)=fa_link(resonance_id,1:3,i)-ip*dversor(1:3)          

          dx=ramber(1,at1)-ramber(1,at2)
          dy=ramber(2,at1)-ramber(2,at2)
          dz=ramber(3,at1)-ramber(3,at2)
          call pbc_displ_vector(lattice_type,amber_cell,amber_kcell,dx,dy,dz)
          rqm=dist_v2(dx,dy,dz)

          flink(1,at2)=flink(1,at2)+fpp(1)*rhq/rqm
          flink(2,at2)=flink(2,at2)+fpp(2)*rhq/rqm 
          flink(3,at2)=flink(3,at2)+fpp(3)*rhq/rqm 

          flink(1,at1)=flink(1,at1)+fpp(1)*(1.0-rhq/rqm)
          flink(2,at1)=flink(2,at1)+fpp(2)*(1.0-rhq/rqm)
          flink(3,at1)=flink(3,at1)+fpp(3)*(1.0-rhq/rqm)

          Ener=Ener+resonance(resonance_id)%weight*Elink_for_h_atom

          !      enddo numlink
       enddo
       
       fdummy(1:3,1:natot)=fdummy(1:3,1:natot)+&
                 resonance(resonance_id)%weight*flink(1:3,1:natot)   

    enddo

    !	write(*,*) 'Total LA energy: ',Ener 

    ! change units 
    ramber(1:3,1:natot)=ramber(1:3,1:natot)/0.529177
    rlink(1:num_resonances,1:3,1:numlink)=rlink(1:num_resonances,1:3,1:numlink)/0.529177

  end subroutine link2

  !*****************************************************************
  ! calculates HL position 

  subroutine link3(rclas,na_qm,natot,nparm,bondeq,cell,lattice_type)

    integer na_qm, natot,nparm

    integer resonance_id, i, at1, at2, k
    real(dp) :: rclas(3,natot), bondeq(nparm)
    real(dp) :: dx, dy, dz, dist, dist_v2, rlink2qm
    character lattice_type
    real(dp) cell(3,3), kcell(3,3)

    call reccel(3,cell,kcell,0)

    ! dist HL-CQM  
    do resonance_id=1,num_resonances
       do i=1,numlink
          at1=linkqm(i,1)
          at2=na_qm+linkmm(i,1)
          dx=rclas(1,at2)-rclas(1,at1)
          dy=rclas(2,at2)-rclas(2,at1)
          dz=rclas(3,at2)-rclas(3,at1)
          call pbc_displ_vector(lattice_type,cell,kcell,dx,dy,dz)
          dist=dist_v2(dx,dy,dz)
          dx=dx/dist
          dy=dy/dist
          dz=dz/dist
          k=parametro(resonance_id,i,1,1)
          if (h_link_scaling_type=='splam') then
              rlink2qm=link_bondeq(resonance_id,i)/0.529177+&
                  link_kbond(resonance_id,i)*(dist-bondeq(k)/0.529177)
          else
              rlink2qm=link_bondeq(resonance_id,i)/0.529177
          endif
          rlink(resonance_id,1,i)=rclas(1,at1)+rlink2qm*dx
          rlink(resonance_id,2,i)=rclas(2,at1)+rlink2qm*dy
          rlink(resonance_id,3,i)=rclas(3,at1)+rlink2qm*dz
       enddo
    enddo

  end subroutine link3

  ! *************************************************************************

  logical function check_out_of_range(atomi,atomj,chi,chj,rmin)

    use sys

    integer, intent(in) :: atomi, atomj
    character(len=4), intent(in) :: chi, chj
    real(dp), intent(in) :: rmin

    real(dp), parameter :: bond_range=1.2_dp
    character :: ch1, ch2

    ! List of typical single bond between species
    ! Bond lengths taken from 
    ! http://www.science.uwaterloo.ca/~cchieh/cact/c120/bondel.html

    check_out_of_range=.false.


    ch1=chi(1:1)  ; ch2=chj(1:1)

    if (ch1=='C'.and.ch2=='C') then
       if (rmin>1.54_dp*bond_range) check_out_of_range=.true.
    else if (ch1=='C'.and.ch2=='H'.or.ch1=='H'.and.ch2=='C') then
       if (rmin>1.09_dp*bond_range) check_out_of_range=.true.
    else if (ch1=='C'.and.ch2=='O'.or.ch1=='O'.and.ch2=='C') then
       if (rmin>1.43_dp*bond_range) check_out_of_range=.true.
    else if (ch1=='C'.and.ch2=='N'.or.ch1=='N'.and.ch2=='C') then
       if (rmin>1.47_dp*bond_range) check_out_of_range=.true.
    else if (ch1=='N'.and.ch2=='N') then
       if (rmin>1.45_dp*bond_range) check_out_of_range=.true.
    else if (ch1=='N'.and.ch2=='H'.or.ch1=='H'.and.ch2=='N') then
       if (rmin>1.01_dp*bond_range) check_out_of_range=.true.
    else if (ch1=='N'.and.ch2=='O'.or.ch1=='O'.and.ch2=='N') then
       !          if (rmin>1.47_dp*bond_range) check_out_of_range=.true.
    else if (ch1=='H'.and.ch2=='H') then
       if (rmin>0.74_dp*bond_range) check_out_of_range=.true.
    else if (ch1=='O'.and.ch2=='O') then
       if (rmin>1.48_dp*bond_range) check_out_of_range=.true.
    else if (ch1=='C'.and.ch2=='S'.or.ch1=='S'.and.ch2=='C') then
       if (rmin>1.82_dp*bond_range) check_out_of_range=.true.
    else if (ch1=='C'.and.ch2=='F'.or.ch1=='F'.and.ch2=='C') then
       if (rmin>1.35_dp*bond_range) check_out_of_range=.true.
    endif

    if (.not.check_out_of_range) then
       print*,'Pair distance(A): ',rmin
       print*,'QM atom: ',atomi
       print*,'MM atom: ',atomj
       print*,'Parameters for bond equilibrium distance: ',chi,'-',chj
       call die('No suitable parameter for Cqm-LA bond. Please set a value in the '//&
            'bonds block of amber.parm This value is needed to decide if the atom '//&
            'are closed enough to make a bond')
    endif

  end function check_out_of_range

  subroutine set_default_H_link_bond_parameters(qmtype,H_bond_length,MM_kbond,H_kbond)

    use sys

     character(len=4), intent(in) :: qmtype
     real(dp), intent(out) :: H_bond_length
     real(dp), intent(in) :: MM_kbond
     real(dp), intent(out) :: H_kbond

     if (qmtype(1:1)=='C') then
        H_bond_length=1.09_dp
        H_kbond=sqrt(MM_kbond/413.0_dp)
     else if (qmtype(1:1)=='N') then
        H_bond_length=1.01_dp
        H_kbond=sqrt(MM_kbond/391.0_dp)
     else if (qmtype(1:1)=='O') then
        H_bond_length=0.96_dp
        H_kbond=sqrt(MM_kbond/366.0_dp)
     else
        call die('No default parameters fof H link bond')
     endif    

  end subroutine set_default_H_link_bond_parameters

end module linkatoms
