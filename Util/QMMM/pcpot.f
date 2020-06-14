! This subroutine calculates the potential due to solvent
! at each point of the mesh and keeps it in Vqm        
!----------------------------------------------------------------------------
      subroutine pcpot(na_qm,na_mm,natot,ntpl,ntm,list,r,
     .  qmcell,pc,rcorteqm,Rho,Vqm,
     .  lattice_type)

      use precision, only: sp, dp, grid_p
      use qmmm_neighbour

      implicit real(dp) (a-h,o-z)
      integer na_qm,na_mm,natot
      integer ntpl,ntm(3)
      real(dp) r(3,natot),qmcell(3,3),pc(na_mm)
      real(dp) rcorteqm, rcorteqm2
      real(grid_p) Rho(ntpl), Vqm(ntpl)
      character lattice_type

      real(dp), dimension(3) :: drij

      real(dp) :: rcorteqm_A

      Vqm=0.0

      rcorteqm_A=rcorteqm/0.529177d0
      rcorteqm2=rcorteqm_A**2     

      if (num_mmvecs==0) return

!----------------------------------------------------------------------------

!     loop over the mesh points
      do iz=0,ntm(3)-1
         do iy=0,ntm(2)-1
            do ix=0,ntm(1)-1

               imesh = 1 + ix + ntm(1)*iy + ntm(1)*ntm(2)*iz

C     Only for grid points where there is density of charge
               if (abs(Rho(imesh)) .gt. 0.0d0) then 

                  xm = qmcell(1,1) * ix/ntm(1) + qmcell(1,2) *
     .                 iy/ntm(2) + qmcell(1,3) * iz/ntm(3)
                  ym = qmcell(2,1) * ix/ntm(1) + qmcell(2,2) *
     .                 iy/ntm(2) + qmcell(2,3) * iz/ntm(3)
                  zm = qmcell(3,1) * ix/ntm(1) + qmcell(3,2) *
     .                 iy/ntm(2) + qmcell(3,3) * iz/ntm(3)

!     loop over MM atoms
                  do n=1,num_mmvecs
                     
                     js=grid_veclist(n)

                     if (lattice_type=='D') then
                        drij(1)=xm-r(1,js)+grid_nr(1,n)*qmcell(1,1)
                        drij(2)=ym-r(2,js)+grid_nr(2,n)*qmcell(2,2)
                        drij(3)=zm-r(3,js)+grid_nr(3,n)*qmcell(3,3)
                     else
                        drij(1)=xm-r(1,js)
                        drij(2)=ym-r(2,js)
                        drij(3)=zm-r(3,js)
                        do l=1,3
                           do m=1,3
                              drij(l)=drij(l)+grid_nr(m,n)*qmcell(l,m)
                           enddo
                        enddo
                     endif

                     d2 = drij(1)*drij(1) + drij(2)*drij(2) + 
     .                    drij(3)*drij(3)

!     calculation of the external potential due to point charges
                     if (d2.gt.rcorteqm2) then
                        d=sqrt(d2)
                        Vqm(imesh)=Vqm(imesh)+pc(js-na_qm)/d 
                     else
                        Vqm(imesh)=Vqm(imesh)+pc(js-na_qm)/rcorteqm_A
                     endif

                  enddo         !!sv atoms

               endif            !! abs(Rho(imesh)) .gt. 0.0d0

            enddo               !!grid
         enddo
      enddo


!     change units
      Vqm=-2.*Vqm

      return
      end

! This subroutine calculates the potential due to solvent
! at each point of the mesh and keeps it in Vqm  using the Ewald
! method
!----------------------------------------------------------------------------
      subroutine pcpot_ewald(na_qm,na_mm,natot,ntpl,ntm,list,
     .      r,qmcell,pc,rcorteqm,ewald_alpha,
     .      kewald_cutoff,Rho,Vqm,lattice_type)

      use precision, only: sp, dp, grid_p
      use qmmm_neighbour

      implicit real(dp) (a-h,o-z)
      integer na_qm,na_mm,natot
      integer ntpl,ntm(3)
      real(dp) r(3,natot),qmcell(3,3),pc(na_mm)
      real(dp) rcorteqm, rcorteqm2
      real(grid_p) Rho(ntpl), Vqm(ntpl)
      character lattice_type

      real(dp) ewald_alpha, pi, sqrt_pi
      real(dp) qm_ewald_alpha, qm_sqrt_ewald_alpha
      real(dp) const2, kmod2
      integer n1, n2, n3, n1max, n2max, n3max
      integer ewald_nmax
      parameter (ewald_nmax=20)
      real(dp) S_real(-ewald_nmax:ewald_nmax,
     .                 -ewald_nmax:ewald_nmax,-ewald_nmax:ewald_nmax)
      real(dp) S_imag(-ewald_nmax:ewald_nmax,
     .                 -ewald_nmax:ewald_nmax,-ewald_nmax:ewald_nmax)
      real(dp) kr, twopi, lattice_volume
      real(dp) krecip(3)
      real(dp) qm_kewald_cutoff, kewald_cutoff, kcut2
      real(dp) kcell(3,3)
      real(dp) scalar_v2

      real(dp), dimension(3) :: drij 

      real(dp) :: rcorteqm_A

      Vqm=0.0
      rcorteqm_A=rcorteqm/0.529177d0
      rcorteqm2=rcorteqm_A**2
      qm_kewald_cutoff=kewald_cutoff*0.529177d0
      qm_ewald_alpha=ewald_alpha*(0.529177d0)**2
      qm_sqrt_ewald_alpha=sqrt(qm_ewald_alpha)

      sqrt_pi=sqrt(acos(-1.0d0))
      twopi=2.0d0*acos(-1.0d0)

      lattice_volume=volcel(qmcell)

! Several constants in the Ewald expresions
      const2=4.0*acos(-1.0d0)/lattice_volume

      call reccel(3,qmcell,kcell,0)

      n1max=INT(qm_kewald_cutoff/(twopi*
     .                  sqrt(scalar_v2(kcell(1,1),kcell(1,2),kcell(1,3),
     .                        kcell(1,1),kcell(1,2),kcell(1,3)))))
      n2max=INT(qm_kewald_cutoff/(twopi*
     .                  sqrt(scalar_v2(kcell(2,1),kcell(2,2),kcell(2,3),
     .                        kcell(2,1),kcell(2,2),kcell(2,3)))))
      n3max=INT(qm_kewald_cutoff/(twopi*
     .                  sqrt(scalar_v2(kcell(3,1),kcell(3,2),kcell(3,3),
     .                        kcell(3,1),kcell(3,2),kcell(3,3)))))
      kcut2=qm_kewald_cutoff*qm_kewald_cutoff       

      if (num_mmvecs==0) return
!----------------------------------------------------------------------------
C     REAL PART OF EWALD SUM
!----------------------------------------------------------------------------
!     loop over the mesh points
      do iz=0,ntm(3)-1
         do iy=0,ntm(2)-1
            do ix=0,ntm(1)-1

               imesh = 1 + ix + ntm(1)*iy + ntm(1)*ntm(2)*iz

C     Only for grid points where there is density of charge
               if (abs(Rho(imesh)) .gt. 0.0d0) then 

                  xm = qmcell(1,1) * ix/ntm(1) + qmcell(1,2) * 
     .                 iy/ntm(2)  + qmcell(1,3) * iz/ntm(3)
                  ym = qmcell(2,1) * ix/ntm(1) + qmcell(2,2) * 
     .                 iy/ntm(2)  + qmcell(2,3) * iz/ntm(3)
                  zm = qmcell(3,1) * ix/ntm(1) + qmcell(3,2) * 
     .                 iy/ntm(2)  + qmcell(3,3) * iz/ntm(3)

!     loop over MM atoms
                  do n=1,num_mmvecs
                     
                     js=grid_veclist(n)

                     if (lattice_type=='D') then
                        drij(1)=xm-r(1,js)+grid_nr(1,n)*qmcell(1,1)
                        drij(2)=ym-r(2,js)+grid_nr(2,n)*qmcell(2,2)
                        drij(3)=zm-r(3,js)+grid_nr(3,n)*qmcell(3,3)
                     else
                        drij(1)=xm-r(1,js)
                        drij(2)=ym-r(2,js)
                        drij(3)=zm-r(3,js)
                        do l=1,3
                           do m=1,3
                              drij(l)=drij(l)+grid_nr(m,n)*
     .                             qmcell(l,m)
                           enddo
                        enddo
                     endif

                     d2 = drij(1)*drij(1) + drij(2)*drij(2) + 
     .                    drij(3)*drij(3)       

!     Calculation of the external potential due to point charges + gaussian 
!     distributions

!     Real-space sum
                     if (d2.gt.rcorteqm2) then
                        d=sqrt(d2)
                        Vqm(imesh)=Vqm(imesh)+pc(js-na_qm)*
     .                       erfc(qm_sqrt_ewald_alpha*d)/d 
                     else
                        Vqm(imesh)=Vqm(imesh)+pc(js-na_qm)*
     .                       erfc(qm_sqrt_ewald_alpha*rcorteqm_A)/
     .                             rcorteqm_A   
                     endif

                  enddo         !!sv atoms

               endif            !! abs(Rho(imesh)) .gt. 0.0d0

            enddo               !!grid
         enddo
      enddo

!----------------------------------------------------------------------------
!         RECIPROCAL-SPACE PART OF EWALD SUM
!----------------------------------------------------------------------------
! Calculate structure factors for classical atoms
      S_real=0.0d0
      S_imag=0.0d0
      if (lattice_type.eq.'D') then
         do n1=-n1max,n1max
            do  n2=-n2max,n2max
               do  n3=-n3max,n3max
                  if (.not.(n1.eq.0.and.n2.eq.0.and.n3.eq.0)) then
                     krecip(1)=n1*twopi*kcell(1,1)
                     krecip(2)=n2*twopi*kcell(2,2)
                     krecip(3)=n3*twopi*kcell(3,3)
                     kmod2=krecip(1)*krecip(1)+krecip(2)*krecip(2)+
     .                         krecip(3)*krecip(3)
                     if (kmod2<=kcut2) then
!     loop over MM atoms
                        do js=na_qm+1,natot
                           xq=r(1,js)
                           yq=r(2,js)
                           zq=r(3,js)
                           kr=krecip(1)*xq+krecip(2)*yq+krecip(3)*zq
                           S_real(n1,n2,n3)=S_real(n1,n2,n3)+
     .                          pc(js-na_qm)*cos(kr)
                           S_imag(n1,n2,n3)=S_imag(n1,n2,n3)+
     .                          pc(js-na_qm)*sin(kr)
                        enddo
                     endif
                  endif
               enddo
            enddo
         enddo
      else
         do n1=-n1max,n1max
            do  n2=-n2max,n2max
               do  n3=-n3max,n3max
                  if (.not.(n1.eq.0.and.n2.eq.0.and.n3.eq.0)) then
                     krecip(1)=twopi*(n1*kcell(1,1)+n2*kcell(2,1)+
     .                    n3*kcell(3,1))
                     krecip(2)=twopi*(n1*kcell(1,2)+n2*kcell(2,2)+
     .                    n3*kcell(3,2))
                     krecip(3)=twopi*(n1*kcell(1,3)+n2*kcell(2,3)+
     .                    n3*kcell(3,3))
                     kmod2=krecip(1)*krecip(1)+krecip(2)*krecip(2)+
     .                         krecip(3)*krecip(3)
                     if (kmod2<=kcut2) then
!     loop over MM atoms
                        do js=na_qm+1,natot
!     Consider only ListQMMM solvent atoms
                           xq=r(1,js)
                           yq=r(2,js)
                           zq=r(3,js)
                           kr=krecip(1)*xq+krecip(2)*yq+krecip(3)*zq
                           S_real(n1,n2,n3)=S_real(n1,n2,n3)+
     .                          pc(js-na_qm)*cos(kr)
                           S_imag(n1,n2,n3)=S_imag(n1,n2,n3)+
     .                          pc(js-na_qm)*sin(kr)
                        enddo
                     endif
                  endif
               enddo
            enddo
         enddo
      endif

C     Calculate the reciprocal part of the potential on QM grid points
c     loop over the mesh points 
      do iz=0,ntm(3)-1
         do iy=0,ntm(2)-1
            do ix=0,ntm(1)-1
               imesh = 1 + ix + ntm(1)*iy + ntm(1)*ntm(2)*iz 
C     Only for grid points where there is density of charge
               if (abs(Rho(imesh)) .gt. 0.0d0) then
                  xm = qmcell(1,1) * ix/ntm(1) + qmcell(1,2) * iy/ntm(2)
     .                 + qmcell(1,3) * iz/ntm(3)
                  ym = qmcell(2,1) * ix/ntm(1) + qmcell(2,2) * iy/ntm(2)
     .                 + qmcell(2,3) * iz/ntm(3)
                  zm = qmcell(3,1) * ix/ntm(1) + qmcell(3,2) * iy/ntm(2)
     .                 + qmcell(3,3) * iz/ntm(3) 
                  if (lattice_type.eq.'D') then
                     do n1=-n1max,n1max
                        do  n2=-n2max,n2max
                           do  n3=-n3max,n3max
                              if (.not.(n1.eq.0.and.n2.eq.0.and.n3.eq.0) 
     $                             )then
                                 krecip(1)=n1*twopi*kcell(1,1)
                                 krecip(2)=n2*twopi*kcell(2,2)
                                 krecip(3)=n3*twopi*kcell(3,3)
                                 kmod2=krecip(1)*krecip(1)
     .                                +krecip(2)*krecip(2)
     .                                +krecip(3)*krecip(3)
                                 if (kmod2<=kcut2) then
                                    kr=krecip(1)*xm+krecip(2)*ym
     .                                     +krecip(3)*zm
                                    Vqm(imesh)=Vqm(imesh) + const2/kmod2
     $                                *exp(-kmod2/(4.0*qm_ewald_alpha))
     $                                *(cos(kr)*S_real(n1,n2,n3)
     $                                     +sin(kr)*S_imag(n1,n2,n3))
                                endif
                              endif
                           enddo
                        enddo
                     enddo
                  else
                     do n1=-n1max,n1max
                        do  n2=-n2max,n2max
                           do  n3=-n3max,n3max
                              if (.not.(n1.eq.0.and.n2.eq.0.and.n3.eq.0) 
     $                             )then
                                 krecip(1)=twopi*(n1*kcell(1,1)+
     .                                n2*kcell(2,1)+n3*kcell(3,1))
                                 krecip(2)=twopi*(n1*kcell(1,2)+
     .                                n2*kcell(2,2)+n3*kcell(3,2))
                                 krecip(3)=twopi*(n1*kcell(1,3)+
     .                                n2*kcell(2,3)+n3*kcell(3,3))
                                 kmod2=krecip(1)*krecip(1)
     .                                +krecip(2)*krecip(2)
     .                                +krecip(3)*krecip(3)
                                 if (kmod2<=kcut2) then
                                    kr=krecip(1)*xq+krecip(2)*yq
     .                                     +krecip(3)*zq
                                    Vqm(imesh)=Vqm(imesh)+const2/kmod2*
     $                                   exp(-kmod2/(4.0*ewald_alpha))*
     $                                   (cos(kr)*S_real(n1,n2,n3)
     $                                   +sin(kr)*S_imag(n1,n2,n3))
                                 endif
                              endif
                           enddo
                        enddo
                     enddo
                  endif         !! lattice_type.eq.'D'
               endif            !! abs(Rho(imesh)) .gt. 0.0d0
            enddo
         enddo
      enddo

! change units
       Vqm=-2.*Vqm

       return
       end

