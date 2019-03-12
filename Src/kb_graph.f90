module kb_graph

  public :: kb_graph_generate

  integer, allocatable, public :: numkb(:), listkbptr(:), listkb(:)

CONTAINS
  
  subroutine kb_graph_generate( scell, na_u, na, isa, xa )

    use precision, only: dp
    use atmfuncs, only : rcut, nofis, nkbfis, floating
    use atm_types, only : nspecies, species, species_info
    use sorting
    use neighbour, only: jna=>jan, xij, r2ij, maxna => maxnna
    use neighbour, only: mneighb

    use alloc, only : re_alloc, de_alloc

    integer,  intent(in) :: na_u        ! num of atoms in unit cell
    integer,  intent(in) :: na          ! num of atoms in supercell
    integer,  intent(in) :: isa(na)
    real(dp), intent(in) :: scell(3,3)  ! Aux super-cell vectors
    real(dp), intent(in) :: xa(3,na)
    

    real(dp), allocatable :: rkbmax(:), rorbmax(:)
    integer, dimension(:),  pointer :: index => null()

    integer  :: is, io, ikb, isel, ia, kna, ka, ks, ind, nna
    real(dp) :: rmaxo, rmaxkb, rmax, rci, rik, rck
    ! tolerance for comparing vector-coordinates   -- probably superfluous
    real(dp), parameter        :: tol = 1.0d-8   
    
    ! Find maximum radius of orbs and KB projectors of each specie
    allocate(rkbmax(nspecies), rorbmax(nspecies))
    do is = 1 , nspecies
       rorbmax(is) = 0.0_dp
       do io = 1 , nofis(is)
          rorbmax(is) = max(rorbmax(is),rcut(is,io))
       end do
       rkbmax(is) = 0.0_dp
       do ikb = 1 , nkbfis(is)
          rkbmax(is) = max(rkbmax(is),rcut(is,-ikb))
       end do
    end do
      
    ! Find maximum range of basis orbitals and KB projectors
    rmaxo    = maxval(rorbmax (1:nspecies))
    rmaxkb   = maxval(rkbmax  (1:nspecies))


    !    rmax = 2._dp * (rmaxo+rmaxkb)
    rmax = (rmaxo+rmaxkb)
    isel = 0
    ! Initialize internal data structures in neighb
    call mneighb( scell, rmax, na, xa, 0, isel, nna )

    if (allocated(numkb)) deallocate(numkb)
    if (allocated(listkbptr)) deallocate(listkbptr)
    if (allocated(listkb)) deallocate(listkb)
    
    allocate(numkb(na_u),listkbptr(na_u))
    do ia = 1 , na_u

       ! initialize number of columns
       numkb(ia) = 0

       ! Find neighbour atoms within maximum range
       call mneighb( scell, rmax, na, xa, ia, isel, nna )
                                ! in case neighbor arrays have expanded
       call re_alloc(index,1,maxna,name="index",routine="kb_graph")

       ! Order neighbours in a well defined way
       call ordvec( tol, 3, nna, xij, index )
       call iorder( jna, 1, nna, index )
       call order ( r2ij, 1, nna, index )

       is  = isa(ia)
       rci = rorbmax(is)

       do kna = 1,nna
          ka = jna(kna)
          ks = isa(ka)
          if (floating(ks)) CYCLE
          rik = sqrt( r2ij(kna) )
          ! It is only necessary to check with
          ! the *largest* projector
          rck = rkbmax(ks)
          if ( rci + rck > rik ) then
             numkb(ia) = numkb(ia) + 1
          endif
       end do
    end do
    ! Count number of non-zeroes and allocate column index
    n_nzskb = sum(numkb(1:na_u))
    ! Create pointers
    listkbptr(1) = 0
    do ia = 2 , na_u
       listkbptr(ia) = listkbptr(ia-1) + numkb(ia-1)
    end do
    allocate(listkb(n_nzskb))

    ! Now fill-in
    do ia = 1 , na_u

       ! initialize number of columns
       numkb(ia) = 0

       ! Find neighbour atoms within maximum range
       call mneighb( scell, rmax, na, xa, ia, isel, nna )
                                ! in case neighbor arrays have expanded
       call re_alloc(index,1,maxna,name="index",routine="kb_graph")

       ! Order neighbours in a well defined way
       call ordvec( tol, 3, nna, xij, index )
       call iorder( jna, 1, nna, index )
       call order ( r2ij, 1, nna, index )

       is  = isa(ia)
       rci = rorbmax(is)

       do kna = 1,nna
          ka = jna(kna)
          ks = isa(ka)
          if (floating(ks)) CYCLE
          rik = sqrt( r2ij(kna) )
          ! It is only necessary to check with
          ! the *largest* projector 
          rck = rkbmax(ks)
          if ( rci + rck > rik ) then
             numkb(ia)   = numkb(ia) + 1
             ind         = listkbptr(ia) + numkb(ia)
             listkb(ind) = ka  ! supercell index
          endif
       end do
    end do
  end subroutine kb_graph_generate
  
end module kb_graph
