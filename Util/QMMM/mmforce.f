c This subroutine computes forces exerted by electrons and nucleus
c on the point charges
c----------------------------------------------------------------------------
      subroutine mmforce(na_qm,na_mm,natot,nesp,ntpl,ntm,r,
     .                   f,stress,pc,ucell,rcorteqmmm,rcorteqm,dvol,Rho,
     .                        lattice_type)

      use precision, only: sp, dp, grid_p
      use qmmm_neighbour

      implicit real(dp) (a-h,o-z)
      integer na_qm,na_mm,natot,nesp
      integer ntpl,ntm(3),n,l
      real(dp) rcorte,vna,grvna(3),fel(3),rr(3)
      real(dp) r(3,natot),ucell(3,3),f(3,natot),dvol,pc(na_mm)
      real(dp) rcorteqmmm, rcorteqmmm2, rcorteqm, rcorteqm2
      real(grid_p) Rho(ntpl)
      real(dp) kcell(3,3), stress(3,3)
      real(dp) :: stress_fact
      character lattice_type

c variables de la lista de vecinos
      real(dp), dimension(3) :: drij
 
      call reccel(3,ucell,kcell,0)

      if (num_mmvecs==0) return

      rcorteqmmm2=(rcorteqmmm/0.529177d0)**2
      rcorteqm2=(rcorteqm/0.529177d0)**2

      stress_fact=1.0/volcel(ucell)

! For now, we put a range of orbital equals to rmax0 (bohrs)

c     --------------------------------------------------------------------------

c     loop over the mesh points 
      do iz=0,ntm(3)-1
         do iy=0,ntm(2)-1
            do ix=0,ntm(1)-1

               imesh = 1 + ix + ntm(1)*iy + ntm(1)*ntm(2)*iz

C     Only for grid points where there is density of charge
               if (abs(Rho(imesh)) .gt. 0.0d0) then 

                  xm = ucell(1,1) * ix/ntm(1) + ucell(1,2) * 
     .                 iy/ntm(2) + ucell(1,3) * iz/ntm(3)
                  ym = ucell(2,1) * ix/ntm(1) + ucell(2,2) * 
     .                 iy/ntm(2) + ucell(2,3) * iz/ntm(3)
                  zm = ucell(3,1) * ix/ntm(1) + ucell(3,2) * 
     .                 iy/ntm(2) + ucell(3,3) * iz/ntm(3)

!     loop over MM atoms
                  do n=1,num_mmvecs
                     
                     js=grid_veclist(n)

                     if (lattice_type=='D') then
                        drij(1)=xm-r(1,js)+grid_nr(1,n)*ucell(1,1)
                        drij(2)=ym-r(2,js)+grid_nr(2,n)*ucell(2,2)
                        drij(3)=zm-r(3,js)+grid_nr(3,n)*ucell(3,3)
                     else
                        drij(1)=xm-r(1,js)
                        drij(2)=ym-r(2,js)
                        drij(3)=zm-r(3,js)
                        do l=1,3
                           do m=1,3
                              drij(l)=drij(l)+grid_nr(m,n)*ucell(l,m)
                           enddo
                        enddo
                     endif
                     
                     d2 = drij(1)*drij(1) + drij(2)*drij(2) + 
     .                    drij(3)*drij(3)

C     Forces exerted on point charges by electrons 
                        if (d2.gt.rcorteqm2) then
                           d=sqrt(d2)
                           De = -2.*dvol*Rho(imesh)*pc(js-na_qm)/d**3
                        else
                           De = 0.0 
                        endif
                        
                        do m=1,3
                           f(m,js) = f(m,js) - De*drij(m)
                           do l=1,3
                              stress(l,m)=stress(l,m)-
     .                             stress_fact*r(l,js)*De*drij(m)
                           enddo
                        enddo

                  enddo         !! at sv

               endif            !! abs(Rho(imesh)) .gt. 0.0d0
               
            enddo               !!grid
         enddo
      enddo

      return
      end

c This subroutine computes forces exerted by electrons and nucleus
c on the point charges using the Ewald method
c----------------------------------------------------------------------------
      subroutine mmforce_ewald(na_qm,na_mm,natot,nesp,ntpl,ntm,r,
     .            f,stress,pc,ucell,rcorteqmmm,rcorteqm,
     .            ewald_alpha,kewald_cutoff,dvol,Rho,lattice_type)

      use sys, only : die
      use precision, only: sp, dp, grid_p
      use qmmm_neighbour

      implicit real(dp) (a-h,o-z)
      integer na_qm,na_mm,natot,nesp
      integer ntpl,ntm(3),n,l
      real(dp) rcorteqmmm, rcorteqm, rcorteqmmm2, rcorteqm2
      real(dp) vna,grvna(3),fel(3),rr(3)
      real(dp) r(3,natot),f(3,natot),dvol,pc(na_mm)
      real(dp) ucell(3,3), kcell(3,3), stress(3,3)
      real(dp) krecip(3)
      real(grid_p) Rho(ntpl)

      real(dp) ewald_alpha
      real(dp) qm_ewald_alpha, qm_sqrt_ewald_alpha
      real(dp) const1, const2, const3, const6, const7, kmod2
      integer n1, n2, n3, n1max, n2max, n3max
      real(dp) kewald_cutoff, qm_kewald_cutoff, kcut2
      real(dp) stress_fact
      integer ewald_nmax
      parameter (ewald_nmax=20)
      real(dp) S_qm_real(-ewald_nmax:ewald_nmax,
     .                 -ewald_nmax:ewald_nmax,-ewald_nmax:ewald_nmax)
      real(dp) S_qm_imag(-ewald_nmax:ewald_nmax,
     .                 -ewald_nmax:ewald_nmax,-ewald_nmax:ewald_nmax)
      real(dp) S_mm_real(-ewald_nmax:ewald_nmax,
     .                 -ewald_nmax:ewald_nmax,-ewald_nmax:ewald_nmax)
      real(dp) S_mm_imag(-ewald_nmax:ewald_nmax,
     .                 -ewald_nmax:ewald_nmax,-ewald_nmax:ewald_nmax)
      real(dp) kronij
      real(dp) twopi, kr, lattice_volume
      character lattice_type

c variables de la lista de vecinos
      real(dp), dimension(3) :: drij

      rcorteqmmm2=(rcorteqmmm/0.529177)**2
      rcorteqm2=(rcorteqm/0.529177d0)**2

      stress_fact=1.0/volcel(ucell)

! For now, we put a range of orbital equals to rmax0 (bohrs)

      qm_kewald_cutoff=kewald_cutoff*0.529177
      qm_ewald_alpha=ewald_alpha*(0.529177)**2
      qm_sqrt_ewald_alpha=sqrt(qm_ewald_alpha)
      twopi=2.0d0*dacos(-1.0d0)

      lattice_volume=volcel(ucell) 

      const1=qm_sqrt_ewald_alpha/sqrt(acos(-1.0d0))
      const2=4.0*acos(-1.0d0)/lattice_volume

      call reccel(3,ucell,kcell,0)

      n1max=INT(qm_kewald_cutoff/(twopi*
     .                  sqrt(scalar_v2(kcell(1,1),kcell(1,2),kcell(1,3),
     .                        kcell(1,1),kcell(1,2),kcell(1,3)))))
      if(n1max>ewald_nmax) call die('mmforce_ewald: n1max>ewald_nmax')
      n2max=INT(qm_kewald_cutoff/(twopi*
     .                  sqrt(scalar_v2(kcell(2,1),kcell(2,2),kcell(2,3),
     .                        kcell(2,1),kcell(2,2),kcell(2,3)))))
      if(n2max>ewald_nmax) call die('mmforce_ewald: n2max>ewald_nmax')
      n3max=INT(qm_kewald_cutoff/(twopi*
     .                  sqrt(scalar_v2(kcell(3,1),kcell(3,2),kcell(3,3),
     .                        kcell(3,1),kcell(3,2),kcell(3,3)))))
      if(n3max>ewald_nmax) call die('mmforce_ewald: n3max>ewald_nmax')

      kcut2=qm_kewald_cutoff*qm_kewald_cutoff

      if (num_mmvecs==0) return

!----------------------------------------------------------------------------
C     REAL PART OF EWALD SUM
!----------------------------------------------------------------------------  

c     loop over the mesh points 
      do iz=0,ntm(3)-1
         do iy=0,ntm(2)-1
            do ix=0,ntm(1)-1

               imesh = 1 + ix + ntm(1)*iy + ntm(1)*ntm(2)*iz

C     Only for grid points where there is density of charge
               if (abs(Rho(imesh)) .gt. 0.0d0) then 

                  xm = ucell(1,1) * ix/ntm(1) + ucell(1,2) * 
     .                 iy/ntm(2)   + ucell(1,3) * iz/ntm(3)
                  ym = ucell(2,1) * ix/ntm(1) + ucell(2,2) * 
     .                 iy/ntm(2)   + ucell(2,3) * iz/ntm(3)
                  zm = ucell(3,1) * ix/ntm(1) + ucell(3,2) * 
     .                 iy/ntm(2)   + ucell(3,3) * iz/ntm(3)

!     loop over MM atoms
                  do n=1,num_mmvecs
                     
                     js=grid_veclist(n)

                     if (lattice_type=='D') then
                        drij(1)=xm-r(1,js)+grid_nr(1,n)*ucell(1,1)
                        drij(2)=ym-r(2,js)+grid_nr(2,n)*ucell(2,2)
                        drij(3)=zm-r(3,js)+grid_nr(3,n)*ucell(3,3)
                     else
                        drij(1)=xm-r(1,js)
                        drij(2)=ym-r(2,js)
                        drij(3)=zm-r(3,js)
                        do l=1,3
                           do m=1,3
                              drij(l)=drij(l)+grid_nr(m,n)
     .                             *ucell(l,m)
                           enddo
                        enddo
                     endif

                     d2 = drij(1)*drij(1) + drij(2)*drij(2) 
     .                    + drij(3)*drij(3)  

C     Forces exerted on point charges by electrons 
                     if (d2.gt.rcorteqm2) then
                        d=sqrt(d2)
                        De = -2.*dvol*Rho(imesh)*pc(js-na_qm)/d2*
     .                       (erfc(qm_sqrt_ewald_alpha*d)
     .                       /d+2.0*const1*exp(-qm_ewald_alpha*d2))          
                     else
                        De = 0.0 
                     endif

                     do m=1,3
                        f(m,js) = f(m,js) - De*drij(m)
                        do l=1,3
                           stress(l,m)=stress(l,m)-
     .                          stress_fact*r(l,js)*De*drij(m)
                        enddo
                     enddo

                  enddo         !! at sv

               endif            !! abs(Rho(imesh)) .gt. 0.0d0

            enddo               !!grid
         enddo
      enddo

!----------------------------------------------------------------------------
!     RECIPROCAL-SPACE PART OF EWALD SUM OF THE POTENTIAL ENERGY & FORCES
!----------------------------------------------------------------------------
!     Calculate structure factors for QM atoms
      S_qm_real=0.0d0
      S_qm_imag=0.0d0
c     loop over the mesh points 
      do iz=0,ntm(3)-1
         do iy=0,ntm(2)-1
            do ix=0,ntm(1)-1
               imesh = 1 + ix + ntm(1)*iy + ntm(1)*ntm(2)*iz 
C     Only for grid points where there is density of charge
               if (abs(Rho(imesh)) .gt. 0.0d0) then 
                  xm = ucell(1,1) * ix/ntm(1) + ucell(1,2) * 
     .                 iy/ntm(2)  + ucell(1,3) * iz/ntm(3)     
                  ym = ucell(2,1) * ix/ntm(1) + ucell(2,2) * 
     .                 iy/ntm(2)  + ucell(2,3) * iz/ntm(3)
                  zm = ucell(3,1) * ix/ntm(1) + ucell(3,2) * 
     .                 iy/ntm(2)  + ucell(3,3) * iz/ntm(3)
                  if (lattice_type.eq.'D') then
!     Reciprocal-space sum
                     do n1=-n1max,n1max
                        do  n2=-n2max,n2max
                           do  n3=-n3max,n3max
                              if (.not.(n1.eq.0.and.n2.eq.0.and.n3.eq.0) 
     .                             ) then 
                                 krecip(1)=n1*twopi*kcell(1,1)
                                 krecip(2)=n2*twopi*kcell(2,2)
                                 krecip(3)=n3*twopi*kcell(3,3)
                                 kmod2=krecip(1)*krecip(1)
     .                                   +krecip(2)*krecip(2)
     .                                   +krecip(3)*krecip(3)
                                 if (kmod2<=kcut2) then
                                    kr=krecip(1)*xm+krecip(2)*ym
     .                                     +krecip(3)*zm
                                    S_qm_real(n1,n2,n3)=
     .                                   S_qm_real(n1,n2,n3)+
     .                                   2.*dvol*Rho(imesh)*cos(kr)
                                    S_qm_imag(n1,n2,n3)=
     .                                   S_qm_imag(n1,n2,n3)+
     .                                   2.*dvol*Rho(imesh)*sin(kr)
                                 endif
                              endif
                           enddo
                        enddo
                     enddo
                  else
!     Reciprocal-space sum
                     do n1=-n1max,n1max
                        do  n2=-n2max,n2max
                           do  n3=-n3max,n3max
                              if (.not.(n1.eq.0.and.n2.eq.0.and.n3.eq.0)
     .                             ) then 
                                 krecip(1)=twopi*(n1*kcell(1,1)
     .                                +n2*kcell(2,1)+n3*kcell(3,1))
                                 krecip(2)=twopi*(n1*kcell(1,2)
     .                                +n2*kcell(2,2)+n3*kcell(3,2))
                                 krecip(3)=twopi*(n1*kcell(1,3)
     .                                +n2*kcell(2,3)+n3*kcell(3,3))
                                 kmod2=krecip(1)*krecip(1)
     .                                   +krecip(2)*krecip(2)
     .                                   +krecip(3)*krecip(3)
                                 if (kmod2<=kcut2) then
                                    kr=krecip(1)*xm+krecip(2)*ym
     .                                   +krecip(3)*zm
                                    S_qm_real(n1,n2,n3)=
     .                                   S_qm_real(n1,n2,n3)+
     .                                   2.*dvol*Rho(imesh)*cos(kr)
                                    S_qm_imag(n1,n2,n3)=
     .                                   S_qm_imag(n1,n2,n3)+
     .                                   2.*dvol*Rho(imesh)*sin(kr)
                                 endif
                              endif
                           enddo
                        enddo
                     enddo
                  endif
               endif            !! abs(Rho(imesh)) .gt. 0.0d0
            enddo
         enddo
      enddo

C     Calculate structure factors for all MM atoms:
      S_mm_real=0.0d0
      S_mm_imag=0.0d0
      if (lattice_type.eq.'D') then
!     Reciprocal-space sum
         do n1=-n1max,n1max
            do  n2=-n2max,n2max
               do  n3=-n3max,n3max
                  if (.not.(n1.eq.0.and.n2.eq.0.and.n3.eq.0)) then
!     loop over MM atoms
                     krecip(1)=n1*twopi*kcell(1,1)
                     krecip(2)=n2*twopi*kcell(2,2)
                     krecip(3)=n3*twopi*kcell(3,3)
                     kmod2=krecip(1)*krecip(1)+krecip(2)*krecip(2)
     .                    +krecip(3)*krecip(3)
                     if (kmod2<=kcut2) then
                        do i=1,na_mm
                           kr=krecip(1)*r(1,i)
     .                          +krecip(2)*r(2,i)
     .                          +krecip(3)*r(3,i)
                           S_mm_real(n1,n2,n3)=S_mm_real(n1,n2,n3)+
     .                          pc(i)*cos(kr)
                           S_mm_imag(n1,n2,n3)=S_mm_imag(n1,n2,n3)+
     .                          pc(i)*sin(kr)
                        enddo
                     endif
                  endif
               enddo
            enddo
         enddo
      else
!     Reciprocal-space sum
         do n1=-n1max,n1max
            do  n2=-n2max,n2max
               do  n3=-n3max,n3max
                  if (.not.(n1.eq.0.and.n2.eq.0.and.n3.eq.0)) then
!     loop over MM atoms
                     krecip(1)=twopi*(n1*kcell(1,1)+
     .                    n2*kcell(2,1)+n3*kcell(3,1))
                     krecip(2)=twopi*(n1*kcell(1,2)+
     .                    n2*kcell(2,2)+n3*kcell(3,2))
                     krecip(3)=twopi*(n1*kcell(1,3)+
     .                    n2*kcell(2,3)+n3*kcell(3,3))
                     kmod2=krecip(1)*krecip(1)+krecip(2)*krecip(2)
     .                    +krecip(3)*krecip(3)
                     if (kmod2<=kcut2) then
                        do i=1,na_mm
                           kr=krecip(1)*r(1,i)
     .                          +krecip(2)*r(2,i)
     .                          +krecip(3)*r(3,i)
                           S_mm_real(n1,n2,n3)=S_mm_real(n1,n2,n3)+
     .                          pc(i)*cos(kr)
                           S_mm_imag(n1,n2,n3)=S_mm_imag(n1,n2,n3)+
     .                          pc(i)*sin(kr)
                        enddo
                     endif
                  endif
               enddo
            enddo
         enddo
      endif

      const4=1.0_dp/lattice_volume
      const6=1.0_dp/(4.0_dp*ewald_alpha)
C     Calculate the reciprocal part of the force on classical atoms
C     due to QM grid points
      if (lattice_type=='D') then
         do n1=-n1max,n1max
            do  n2=-n2max,n2max
               do  n3=-n3max,n3max
                  if (.not.(n1.eq.0.and.n2.eq.0.and.n3.eq.0)) then
c     loop over MM atoms
                     krecip(1)=n1*twopi*kcell(1,1)
                     krecip(2)=n2*twopi*kcell(2,2)
                     krecip(3)=n3*twopi*kcell(3,3)
                     kmod2=krecip(1)*krecip(1)+krecip(2)*krecip(2)
     .                    +krecip(3)*krecip(3)
                     if (kmod2<=kcut2) then
                        const3=const2/kmod2*
     .                          exp(-kmod2/(4.0*qm_ewald_alpha))
                        do js=na_qm+1,natot
                           xq=r(1,js)
                           yq=r(2,js)
                           zq=r(3,js)
                           kr=krecip(1)*xq+krecip(2)*yq
     .                          +krecip(3)*zq
                           De=const3*pc(js)*
     .                          (cos(kr)*S_qm_imag(n1,n2,n3)
     $                          -sin(kr)*S_qm_real(n1,n2,n3))
                           do k=1,3
                              f(k,js) = f(k,js) + De*krecip(k)
                           enddo
                        enddo
                        De2=const3*const4*
     .                       (S_mm_imag(n1,n2,n3)*S_qm_imag(n1,n2,n3)
     $                       +S_mm_real(n1,n2,n3)*S_qm_real(n1,n2,n3))
                        const7=2.0_dp*(1.0_dp+kmod2*const6)/kmod2
                        do k=1,3
                           do l=1,3
                              kronij=real(int(((l+k)-abs(l-k))/
     $                             ((l+k)+abs(l-k))),kind=dp)
                              stress(l,k)=stress(l,k)+
     $                      De2*(kronij-const7*krecip(l)*krecip(k))
                           enddo
                        enddo
                     endif
                  endif
               enddo
            enddo
         enddo
      else
C     Calculate the reciprocal part of the force on classical atoms
!     due to QM grid points
         do n1=-n1max,n1max
            do  n2=-n2max,n2max
               do  n3=-n3max,n3max
                  if (.not.(n1.eq.0.and.n2.eq.0.and.n3.eq.0)) then
c     loop over MM atoms
                     krecip(1)=twopi*(n1*kcell(1,1)+n2*kcell(2,1)+
     .                    n3*kcell(3,1))
                     krecip(2)=twopi*(n1*kcell(1,2)+n2*kcell(2,2)+
     .                    n3*kcell(3,2))
                     krecip(3)=twopi*(n1*kcell(1,3)+n2*kcell(2,3)+
     .                    n3*kcell(3,3))
                     kmod2=krecip(1)*krecip(1)+krecip(2)*krecip(2)
     .                    +krecip(3)*krecip(3)
                     if (kmod2<=kcut2) then
                        const3=const2/kmod2*
     .                          exp(-kmod2/(4.0*qm_ewald_alpha))
c     loop over MM atoms
                        do js=na_qm+1,natot
                           xq=r(1,js)
                           yq=r(2,js)
                           zq=r(3,js)
                           kr=krecip(1)*xq+krecip(2)*yq
     .                          +krecip(3)*zq
                           De=const3*
     .                          (cos(kr)*S_qm_imag(n1,n2,n3)
     $                          -sin(kr)*S_qm_real(n1,n2,n3))
                           do k=1,3
                              f(k,js) = f(k,js) + De*krecip(k)
                           enddo
                        enddo
                        De2=const3*const4*
     .                       (S_mm_imag(n1,n2,n3)*S_qm_imag(n1,n2,n3)
     $                       +S_mm_real(n1,n2,n3)*S_qm_real(n1,n2,n3))
                        const7=2.0_dp*(1.0_dp+kmod2*const6)/kmod2
                        do k=1,3
                           do l=1,3
                              kronij=real(int(((l+k)-abs(l-k))/
     $                             ((l+k)+abs(l-k))),kind=dp)
                              stress(l,k)=stress(l,k)+
     $                         De2*(kronij-const7*krecip(l)*krecip(k))
                           enddo
                        enddo
                     endif
                  endif
               enddo
            enddo
         enddo
      endif

      return
      end

