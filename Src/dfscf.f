      subroutine dfscf( NO, endC,  listC, C, iOrb,
     .                  NX, gradCi, NSPmax, NSP,
     .                  NP, endCt, listCt, CtoCt,
     .                  NDmax, numDs, listDs, Dscf, Datm,
     .                  VolCel, Vscf, Vatm, NVmax, numVs, listVs,
     .                  fa, MaxDi, Di )

C ********************************************************************
C Adds the SCF contribution of a single orbital to atomic forces or
C stress.
C a single orbital.
C Written by P.Ordejon and J.M.Soler. May'95.
C *********************** INPUT **************************************
C integer NO              : Number of basis orbitals
C integer endC(NO)        : Accumulated umber of nonzero elements in
C                           each row of C
C integer listC(*)        : List of nonzero elements in each row of C
C real*4  C(NSPmax,*)     : Values of nonzero elements in C
C integer iOrb            : Orbital whose force is to be found
C integer NX              : Space dimension of forces or stress
C                           NX=3 for forces, NX=9 for stress.
C real*8  gradCi(NX,NSPmax,*): Gradient of orbital iOrb at mesh points,
C                             or tensor x(i)*gradCi(j)/VolCel for stress
C integer NSPmax          : First dimension of array C
C integer NSP             : Number of sub-points of each mesh point
C integer NP              : Number of columns in C (num. of mesh points)
C integer endCt(NP)       : Accumulated number of nonzero elements in
C                           each column of C
C integer listCt(*)       : List of nonzero elements in each column of C
C integer CtoCt(*)        : Index of matrix C where the nonzero
C                           elements in each column are stored
C integer NDmax           : First dimension of listD and Dscf, and
C                           maximum number of nonzero elements in
C                           any row of Dscf
C integer numDs(NO)       : Number of nonzero elemts in each row of Dscf
C integer listDs(NVmax,NO): List of nonzero elements in each row of Dscf
C real*8  Dscf(NVmax,NO)  : Value of nonzero elemens in each row of Dscf
C real*8  Datm(NO)        : Occupations of basis orbitals in free atom
C real*8  VolCel          : Unit cell volume
C real*4  Vscf(NP)        : Value of SCF potential at the mesh points
C real*4  Vatm(NP)        : Value of Harris potential (Hartree potential
C                           of sum of atomic desities) at mesh points
C integer NVmax           : First dimension of listV and Vs, and maxim.
C                           number of nonzero elements in any row of Vs
C integer numVs(NO)       : Number of nonzero elements in each row of Vs
C integer listVs(NVmax,NO): List of nonzero elements in each row of Vs
C integer MaxDi           : Total size of auxiliary array Di
C ******************INPUT AND OUTPUT **********************************
C real*8  fa(NX)          : Atomic force or stress, to which the SCF
C                           contribution is added
C ********************* AUXILIARY *************************************
C real*8  Di(MaxDi)       : Array that must be initialized to zero
C                           before call for first orbital and which
C                           MUST NOT BE MODIFIED before call for last
C                           orbital. MaxDi must be at least NO
C *********************************************************************

      implicit none

      integer
     .   NO, NP, NDmax, NVmax, NSP, NSPmax, NX, MaxDi, iOrb,
     .   endC(0:NO),  listC(*),
     .   endCt(0:NP), listCt(*), CtoCt(*),
     .   numDs(NO), listDs(NDmax,NO),
     .   numVs(NO), listVs(NVmax,NO)

      real
     .   C(NSPmax,*), Vscf(NSPmax,NP), Vatm(NSPmax,NP)

      double precision
     .   gradCi(NX,NSPmax,*),
     .   Dscf(NDmax,NO), Datm(NO), VolCel,
     .   fa(NX), Di(NO)

      integer i, in, ini, ip, isp, ix, j, jn, kn, MaxNSP, MaxX
      parameter ( MaxNSP = 8, MaxX = 9 )
      double precision DiiCiV, DijCj, dVol, VgrCi(MaxX,MaxNSP)
      
      call timer('dfscf',1)

C     Check array sizes
      if (NSP .gt. MaxNSP) stop 'DFSCF: dimension MaxNSP too small'
      if (NO  .gt. MaxDi)  stop 'DFSCF: dimension MaxDi too small'
      if (NX  .gt. MaxX)   stop 'DFSCF: dimension MaxX too small'

C     A shorter name for orbital index
      i = iOrb
      dVol = VolCel / (NP*NSP)

C     Copy full row i of density matrix to Di(j)
      do 70 in = 1, numDs(i)
        j = listDs(in,i)
        Di(j) = Dscf(in,i)
   70 continue

C     Loop on mesh points of orbital i
      do 130 in = 1+endC(i-1), endC(i)
        ip = listC(in)

C       Loop over sub-points and add 2*Datm_ii*<Ci|Vatm|gradCi>
        ini = in - endC(i-1)
        do 90 isp = 1, NSP
          DiiCiV = 2.d0 * dVol * Datm(i) * C(isp,in) * Vatm(isp,ip)
          do 80 ix = 1,NX
            fa(ix) = fa(ix) - DiiCiV * gradCi(ix,isp,ini)
            VgrCi(ix,isp) = Vscf(isp,ip) * gradCi(ix,isp,ini)
   80     continue
   90   continue

C       Loop on other non-zero orbitals in mesh point
        do 110 kn = 1+endCt(ip-1), endCt(ip)
          j = listCt(kn)
          jn = CtoCt(kn)

C         Loop over sub-points and add 2*Dscf_ij*<Cj|Vscf|gradCi>
          do 105 isp = 1, NSP
            DijCj = 2.d0 * dVol * Di(j) * C(isp,jn)
            do 100 ix = 1,NX
              fa(ix) = fa(ix) + DijCj * VgrCi(ix,isp)
  100       continue
  105     continue

  110   continue
  130 continue

C     Restore Di for next orbital i
      do 180 in = 1, numDs(i)
        j = listDs(in,i)
        Di(j) = 0.d0
  180 continue
  
      call timer('dfscf',2)
      end

