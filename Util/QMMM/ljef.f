c This subroutine calculates the solute-solvent LJ energy and forces
c----------------------------------------------------------------------------
      subroutine ljef(na_qm,nac,natot,r,Em,Rm,f,stress,Es,
     .                         rcorteqmmm,cell,lattice_type)

      use precision, only : dp
      use qmmm_neighbour

      implicit real(dp) (a-h,o-z)
      real(dp)
     .  r(3,natot), f(3,natot),Em(natot), Rm(natot), rr(3),flj(3)
      real(dp) stress(3,3), cell(3,3)
      real(dp) rcorteqmmm, rcorteqmmm2
      character lattice_type

c variables de la lista de vecinos
      real(dp), dimension(3) :: drij
      integer n_pointer

      real(dp) :: factor, dfactor, arg, earg
      real(dp) :: d_grimme=20.0_dp

      Elj=0.0D0
      flj=0.0D0

      rcorteqmmm2=(rcorteqmmm/0.529177d0)**2

c----------------------------------------------------------------------------

      stress_fact=1.0_dp/volcel(cell)

      n_pointer=1   
c     loop over QM atoms without considering the link atoms
      do  j1=1,na_qm          

C     n_pointer: points the first neb atom of i in the neb list.
         do k=n_pointer,qmmm_veclistxat(j1)
            
            j2=qmmm_veclist(k)

            if (lattice_type=='D') then
               drij(1)=r(1,j1)-r(1,j2)+qmmm_nr(1,k)*cell(1,1)
               drij(2)=r(2,j1)-r(2,j2)+qmmm_nr(2,k)*cell(2,2)
               drij(3)=r(3,j1)-r(3,j2)+qmmm_nr(3,k)*cell(3,3)
            else
               drij(1)=r(1,j1)-r(1,j2)
               drij(2)=r(2,j1)-r(2,j2)
               drij(3)=r(3,j1)-r(3,j2)
               do l=1,3
                  do m=1,3
                     drij(l)=drij(l)+qmmm_nr(m,k)*cell(l,m)
                  enddo
               enddo
            endif

            dd2=drij(1)*drij(1) + drij(2)*drij(2) + 
     .           drij(3)*drij(3)

            if (dd2 .le. rcorteqmmm2) then
               dd=sqrt(dd2)
c     Energy and forces from LJ term
               epsilon=sqrt(Em(j2)*Em(j1))
               if (epsilon==0.0d0) cycle
               sigma=0.50D0*(Rm(j2)+Rm(j1))
               t1=sigma**6
               B=4.0D0*epsilon*t1
               A=B*t1

               if (j2<=na_qm) then
                  earg=exp(-d_grimme*(0.5D0*dd/sigma-1.0d0))
                  factor=1.0d0/(1.0+arg)
                  dfactor=0.5D0*d_grimme/sigma*arg/factor**2
                  arg=A/dd**12 -B/dd**6
                  Elj = Elj+ factor*arg 
                  fej = (-12.0D0*A/dd**14 + 6.0D0*B/dd**8)*factor+
     .                 arg*dfactor
               else
                  Elj = Elj+  A/dd**12 -B/dd**6
                  fej = -12.0D0*A/dd**14 + 6.0D0*B/dd**8
               endif

               flj(1)=-2.0*fej*drij(1)
               flj(2)=-2.0*fej*drij(2)
               flj(3)=-2.0*fej*drij(3)

               do i=1,3
                  f(i,j1)=f(i,j1) + flj(i)
                  f(i,j2)=f(i,j2) - flj(i) 
                  do j=1,3
                     stress(j,i)=stress(j,i)+stress_fact*
     .                    drij(j)*flj(i)
                  enddo
               enddo

            endif

         enddo                  !! at st

         n_pointer = qmmm_veclistxat(j1) + 1

      enddo                     !! at sv

      Es = 2.0*Elj
      
      return
      end

