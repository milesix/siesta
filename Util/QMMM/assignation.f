c subroutine that asignates masses and specie
      subroutine assign(na_qm,na_mm,atname,iza,izs,masst)

      use precision
      use periodic_table
      use fdf
      use sys
      use m_qmmm_fdf, only : fdf_block_qmmm 

      implicit          none
      integer           i,j,k, na_qm, na_mm, iu
      integer           iza(na_qm), izs(na_qm+na_mm)
      real(dp)          masst(na_qm+na_mm), mass
      character         atname(na_mm)*4,atn4*4,atn1*1,atn2*2
      character*2       sym(na_qm+na_mm), name 

C     Assigantes solute atomic masses (masst) and atomic symbol (sym)
       masst=0.0
       do i=1,na_qm
          masst(i) = ATMASS(iza(i))
          sym(i) = SYMBOL(iza(i))
          if(masst(i).eq.0.0) then
             call die("assign: There are solute atoms without mass")
          endif
       enddo

C Assigantes solvent atomic masses (masst) and atomic symbol (sym)
        izs = 0
        izs(1:na_qm)=iza(1:na_qm)
        k=na_qm+1
        do i=1,na_mm
        atn4 = atname(i)
        atn1 = atn4(1:1)
        atn2 = atn4(1:2)
        if (atn1.eq.'H')  izs(k)=1
        if (atn1.eq.'C')  izs(k)=6
        if (atn1.eq.'N')  izs(k)=7
        if (atn1.eq.'O')  izs(k)=8
        if (atn1.eq.'F')  izs(k)=9
        if (atn2.eq.'Na') izs(k)=11
        if (atn2.eq.'Mg') izs(k)=12
        if (atn1.eq.'P')  izs(k)=15
        if (atn1.eq.'S')  izs(k)=16
        if (atn2.eq.'Cl') izs(k)=17
        if (atn1.eq.'K')  izs(k)=19
        if (atn2.eq.'Ca') izs(k)=20
        if (atn2.eq.'Mn') izs(k)=25
        if (atn2.eq.'FE') izs(k)=26
        if (atn2.eq.'Cu') izs(k)=29
        if (atn2.eq.'Zn') izs(k)=30
        if (atn2.eq.'Br') izs(k)=35
        if (atn1.eq.'I')  izs(k)=53
        if (atn2.eq.'Pt') izs(k)=78
       if(izs(k).eq.0) then
       call die('assign: There are solvent atoms without atomic number')
       endif
        masst(k) = ATMASS (izs(k))
        sym(k)= SYMBOL(izs(k))
       if(masst(k).eq.0.0) then
       call die('assign: There are solvent atoms without mass')
       endif
        k=k+1
        enddo

c Read atomic masses of different atoms from fdf block 

      if ( fdf_block_qmmm('NewMasses', iu) ) then
 5     continue
       read(iu,'(A2)',advance='no',err=10,end=10) name 
       if(name.eq.'%e'.or.name.eq.'%E') return 
       read(iu,*,err=10,end=10) mass 
       write(6,"(/,a, A2, a, f6.2)")
     . 'assign: Read atomic mass of:  ',name,'as',mass

c assignates new masses
        do i=1,na_qm+na_mm
        if(name.eq.sym(i)) masst(i) = mass
        enddo

       goto 5
 6     continue
      endif

      return
 10   stop 'assign: Problem reading from "NewMasses" block'
      end 

