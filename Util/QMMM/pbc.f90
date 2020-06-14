
      subroutine reinserting_atoms_in_box(lattice_type,natot,na_qm,nac,&
                                                     noaa,rclas,cell)

      use precision, only: dp
      implicit none

      character lattice_type
      integer natot, na_qm, nac
      real(dp) rclas(3,natot)
      real(dp) cell(3,3), kcell(3,3) 
      character noaa(nac)*4

      real(dp) twopi, rnk
      integer i, j, k, l, nk
      integer shift

      shift=0
      i=0
      if (lattice_type=='D') then

         do 
            i=i+1+shift
            if (i>natot) exit
            shift=0
            if (i.gt.na_qm) then
               if (noaa(i-na_qm).eq.'HOH') then
                  shift=2
               else
                  shift=0
               endif
            endif
            if (rclas(1,i)<0.0) then
               do j=0,shift
                  rclas(1,i+j)=rclas(1,i+j)+cell(1,1)
               enddo
            else if (rclas(1,i)>cell(1,1)) then
               do j=0,shift
                  rclas(1,i+j)=rclas(1,i+j)-cell(1,1)
               enddo
            endif
            if (rclas(2,i)<0.0) then
               do j=0,shift
                  rclas(2,i+j)=rclas(2,i+j)+cell(2,2)
               enddo
            else if (rclas(2,i)>cell(2,2)) then
               do j=0,shift
                  rclas(2,i+j)=rclas(2,i+j)-cell(2,2)
               enddo
            endif
            if (rclas(3,i)<0.0) then
               do j=0,shift
                  rclas(3,i+j)=rclas(3,i+j)+cell(3,3)
               enddo
            else if (rclas(3,i)>cell(3,3)) then
               do j=0,shift
                  rclas(3,i+j)=rclas(3,i+j)-cell(3,3)
               enddo
            endif
         enddo
   
      else

         call reccel(3,cell,kcell,0)
         do
            i=i+1+shift
            if (i>natot) exit
            shift=0
            if (i.gt.na_qm) then
               if (noaa(i-na_qm).eq.'HOH') then
                  shift=2
               else
                  shift=0
               endif
            endif
            do k=1,3
               rnk=0.0
!     Here it is assumed that kcell is the matrix, whose rows are 
!     the reciprocal vectors without 2*pi factor
               do l=1,3
                  rnk=rnk+rclas(k,i)*kcell(l,k)
               enddo
               if (rnk>=0.0_dp) then
                  nk=INT(rnk)
               else
                  nk=NINT(rnk-0.5_dp)
               endif
               if (nk.ne.0) then
                  do j=0,shift
                     do l=1,3
                        rclas(k,i+j)=rclas(k,i+j)-cell(l,k)*nk
                     enddo
                  enddo
               endif
            enddo   
         enddo

      endif

      end subroutine reinserting_atoms_in_box

      subroutine pbc_displ_vector(lattice_type,cell,kcell,dx,dy,dz)

      use precision, only: dp
      implicit none

      character lattice_type
      real(dp) dx, dy, dz
      real(dp) dr(3)
      real(dp) cell(3,3), kcell(3,3) 

      integer i, k, l, nk

      dr(1)=dx
      dr(2)=dy
      dr(3)=dz

      if (lattice_type=='D') then

         do k=1,3
            dr(k) = dr(k) - ANINT(dr(k)/cell(k,k))*cell(k,k)
         enddo

      else

         do k=1,3
!     Here it is assumed that kcell is matrix of the reciprocal vectors
!     without 2*pi factor
            nk=ANINT(dr(1)*kcell(k,1)+dr(2)*kcell(k,2)+dr(3)*kcell(k,3))
            if (nk.ne.0) then
               do l=1,3
                  dr(k)=dr(k)-cell(l,k)*nk
               enddo
            endif
         enddo   

      endif

      dx=dr(1)
      dy=dr(2)
      dz=dr(3)

      end subroutine pbc_displ_vector

      subroutine get_pbc_vectors(lattice_type,cell,kcell,drij,nr)

      use precision, only: dp
      implicit none

      character lattice_type
      real(dp) cell(3,3), kcell(3,3)
      real(dp) drij(3)
      integer nr(3)

      integer k

      if (lattice_type=='D') then
         do k=1,3
            nr(k) = -ANINT(drij(k)/cell(k,k))
         enddo
      else
         do k=1,3
!     Here it is assumed that kcell is matrix of the reciprocal vectors
!     without 2*pi factor
            nr(k)=-ANINT(drij(1)*kcell(k,1)+ &
                 drij(2)*kcell(k,2)+drij(3)*kcell(k,3))
         enddo   
      endif

      end subroutine get_pbc_vectors

      character function get_lattice_type(cell)

       use precision, only: dp
       implicit none

! RETURN 'D' IF THE LATTICE VECTORS LIE ALONG THE X-, Y-
! AND Z-AXIS RESPECTIVELY. IN OTHER WORDS, THE CELL MATRIX IS DIAGONAL.
! RETURN 'G' OTHERWISE.

        real(dp) cell(3,3)
        integer i, j

        get_lattice_type='D'

        do i=1,3
           do j=1,3
              if (i.ne.j) then
                 if (cell(i,j).ne.0.0) get_lattice_type='G' 
              endif
           enddo
        enddo

      end function get_lattice_type

      subroutine reccel( N, A, B, IOPT )

        use precision, only: dp
        use sys, only: die

        implicit none

!  CALCULATES RECIPROCAL LATTICE VECTORS B.
!  THEIR PRODUCT WITH DIRECT LATTICE VECTORS A IS 1 (IF IOPT=0) OR
!  2*PI (IF IOPT=1). N IS THE SPACE DIMENSION.
!  WRITTEN BY J.M.SOLER.

      integer :: n, iopt
      real(dp):: A(N,N),B(N,N), c, ci

      integer :: i

      C=1.D0
      IF (IOPT.EQ.1) C=2.D0*ACOS(-1.D0)

      IF (N .EQ. 1) THEN
        B(1,1) = C / A(1,1)
      ELSEIF (N .EQ. 2) THEN
        C = C / (A(1,1)*A(2,2) - A(1,2)*A(2,1))
        B(1,1) =  A(2,2)*C
        B(1,2) = (-A(2,1))*C
        B(2,1) = (-A(1,2))*C
        B(2,2) =  A(1,1)*C
      ELSEIF (N .EQ. 3) THEN
        B(1,1)=A(2,2)*A(3,3)-A(3,2)*A(2,3)
        B(2,1)=A(3,2)*A(1,3)-A(1,2)*A(3,3)
        B(3,1)=A(1,2)*A(2,3)-A(2,2)*A(1,3)
        B(1,2)=A(2,3)*A(3,1)-A(3,3)*A(2,1)
        B(2,2)=A(3,3)*A(1,1)-A(1,3)*A(3,1)
        B(3,2)=A(1,3)*A(2,1)-A(2,3)*A(1,1)
        B(1,3)=A(2,1)*A(3,2)-A(3,1)*A(2,2)
        B(2,3)=A(3,1)*A(1,2)-A(1,1)*A(3,2)
        B(3,3)=A(1,1)*A(2,2)-A(2,1)*A(1,2)
        DO 20 I=1,3
          CI=C/(A(1,I)*B(1,I)+A(2,I)*B(2,I)+A(3,I)*B(3,I))
          B(1,I)=B(1,I)*CI
          B(2,I)=B(2,I)*CI
          B(3,I)=B(3,I)*CI
  20    CONTINUE
      ELSE
         call die('RECCEL: NOT PREPARED FOR N>3')
      ENDIF

      end subroutine reccel


