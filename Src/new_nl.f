      module new_nl_m
      implicit none

      public :: new_nl
      
      CONTAINS
      subroutine new_nl( scell, nua, na, isa, xa, indxua,
     .                   maxnh, maxnd, lasto, lastkb, iphorb, 
     .                   iphKB, numd, listdptr, listd, numh, 
     .                   listhptr, listh, nspin, Dscf, Enl, 
     .                   fa, stress, H , matrix_elements_only)
C *********************************************************************
C Calculates non-local (NL) pseudopotential contribution to total 
C energy, atomic forces, stress and hamiltonian matrix elements.
C Energies in Ry. Lengths in Bohr.
C     Written by A. Garcia (March 2019), based on
C     code by J.Soler and P.Ordejon, June 1997.
C **************************** INPUT **********************************
C real*8  scell(3,3)       : Supercell vectors SCELL(IXYZ,IVECT)
C integer nua              : Number of atoms in unit cell
C integer na               : Number of atoms in supercell
C integer isa(na)          : Species index of each atom
C real*8  xa(3,na)         : Atomic positions in cartesian coordinates
C integer indxua(na)       : Index of equivalent atom in unit cell
C integer maxnh            : First dimension of H and listh
C integer maxnd            : Maximum number of elements of the
C                            density matrix
C integer lasto(0:na)      : Position of last orbital of each atom
C integer lastkb(0:na)     : Position of last KB projector of each atom
C integer iphorb(no)       : Orbital index of each orbital in its atom,
C                            where no=lasto(na)
C integer iphKB(nokb)      : Index of each KB projector in its atom,
C                            where nokb=lastkb(na)
C integer numd(nuo)        : Number of nonzero elements of each row of the
C                            density matrix
C integer listdptr(nuo)    : Pointer to the start of each row (-1) of the
C                            density matrix
C integer listd(maxnd)     : Nonzero hamiltonian-matrix element column 
C                            indexes for each matrix row
C integer numh(nuo)        : Number of nonzero elements of each row of the
C                            hamiltonian matrix
C integer listhptr(nuo)    : Pointer to the start of each row (-1) of the
C                            hamiltonian matrix
C integer listh(maxnh)     : Nonzero hamiltonian-matrix element column 
C                            indexes for each matrix row
C integer nspin            : Number of spinor components (really min(nspin,2))
C                            If computing only matrix elements, it
C                            can be set to 1.
C logical matrix_elements_only:
C integer Dscf(maxnd,nspin): Density matrix. Not touched if computing
C                            only matrix elements.
C ******************* INPUT and OUTPUT *********************************
C real*8 fa(3,na)          : NL forces (added to input fa)
C real*8 stress(3,3)       : NL stress (added to input stress)
C real*8 H(maxnh,nspin)    : NL Hamiltonian (added to input H)
C **************************** OUTPUT *********************************
C real*8 Enl               : NL energy
C *********************************************************************
C
C  Modules
C
      use precision,     only : dp
      use parallel,      only : Node, Nodes
      use parallelsubs,  only : GetNodeOrbs, LocalToGlobalOrb
      use atmfuncs,      only : rcut, epskb, orb_gindex, kbproj_gindex
      use m_new_matel,   only : new_matel
      use kb_graph,      only : numkb, listkbptr, listkb
      use kb_graph,      only : kb_graph_generate
      use kb_graph,      only : kb_graph_generate, kb_graph_print
      use intrinsic_missing, only: VNORM
      use atomlist,      only : iaorb

      integer, intent(in) ::
     .   maxnh, na, maxnd, nspin, nua

      integer, intent(in)  ::
     .  indxua(na), iphKB(*), iphorb(*), isa(na),  
     .  lasto(0:na), lastkb(0:na), listd(maxnd), listh(maxnh),
     .  numd(*), numh(*), listdptr(*), listhptr(*)

      real(dp), intent(in) :: scell(3,3), Dscf(maxnd,nspin),
     .                        xa(3,na)
      real(dp), intent(inout) :: fa(3,nua), stress(3,3)
      real(dp), intent(inout) :: H(maxnh,nspin)
      real(dp), intent(out)   :: Enl
      logical, intent(in)     :: matrix_elements_only

      real(dp) ::   volcel
      external ::   timer, volcel

      integer
     .  ia, ik, ikb, ind, indkb, ja, js,
     .  io, iio, ioa, is, joa, ispin, ix, ig, jg, kg,
     .  j, jo, jx, ka, ko, koa, ks, kua,
     .  no, nuo, nuotot

      real(dp) :: rci, rcj, rck, dik, djk
      real(dp) :: Di, xki(3), xkj(3)
      
      real(dp)
     .  Cijk, epsk, fik, rki, rmax, rmaxkb, rmaxo, 
     .  Sik, Sjk, volume, grSik(3), grSjk(3)

C ......................

C Start time counter
      call timer( 'new_nl', 1 )

C Find unit cell volume
      volume = volcel( scell ) * nua / na

      no = lasto(na)
      nuotot = lasto(nua)
      call GetNodeOrbs(nuotot,Node,Nodes,nuo)

      call timer( 'kb_graph', 1 )
      call kb_graph_generate( scell, nua, na, isa, xa )
      call kb_graph_print( nua, na, isa, xa )
      call timer( 'kb_graph', 2 )

      Enl = 0.0d0

      do io = 1, nuo
           ! Find out in which atom we are
         call LocalToGlobalOrb(io,Node,Nodes,iio)
         ia = iaorb(iio)
         is = isa(ia)
         ioa = iphorb(iio)
         rci = rcut(is,ioa)
         ig = orb_gindex(is,ioa)
         Di = 0.0_dp
         
         do j = 1, numh(io)
            ind = listhptr(io) + j
            jo = listh(ind)
            ja = iaorb(jo)
            js = isa(ja)
            joa = iphorb(jo)
            rcj = rcut(js,joa)
            jg = orb_gindex(js,joa)
            if (.not. matrix_elements_only) then
               do ispin = 1, nspin ! Both spins add up... (nspin=spin%spinor here)
                  Di = Di + Dscf(ind,ispin)
               enddo
            endif

            do ik = 1, numkb(ia)
               indkb = listkbptr(ia) + ik
               ka = listkb(indkb)
               kua = indxua(ka) ! Used only if forces and energies are comp.
               ks = isa(ka)
               ! xki is from k to i
               xki(:) = xa(:,ia) - xa(:,ka)
               xkj(:) = xa(:,ja) - xa(:,ka)
               dik = VNORM( xki )
               djk = VNORM( xkj )

               do ikb = lastkb(ka-1)+1,lastkb(ka)
                  koa = iphKB(ikb)  ! koa is negative
                  rck = rcut(ks,koa)
                  !
                  if ( (rck + rci) < dik ) CYCLE
                  if ( (rck + rcj) < djk ) CYCLE
                  !
                  epsk = epskb(ks,koa)
                  kg = kbproj_gindex(ks,koa)
                  call new_MATEL( 'S', kg, ig, xki, Sik, grSik)
                  call new_MATEL( 'S', kg, jg, xkj, Sjk, grSjk)
                  !
                  do ispin = 1,nspin ! Only diagonal parts
                     H(ind,ispin) = H(ind,ispin) + epsk * Sik*Sjk
                  enddo

                  if (.not. matrix_elements_only) then
                     Cijk = Di * epsk    
                     Enl = Enl + Cijk * Sik * Sjk
                     do ix = 1,3
                        fik = 2.d0 * Cijk * Sjk * grSik(ix)  ! Check this
                        fa(ix,ia)  = fa(ix,ia)  - fik
                        fa(ix,kua) = fa(ix,kua) + fik
                        do jx = 1,3
                           stress(jx,ix) = stress(jx,ix) +
     &                          xki(jx) * fik / volume       ! Check this
                        enddo
                     enddo
                  endif
                  
               enddo            ! ikb
            enddo               ! ik: atoms with KB close to ia
         enddo                  ! jo
      enddo                     ! io


      call timer( 'new_nl', 2 )
      
      end subroutine new_nl
      end module new_nl_m
