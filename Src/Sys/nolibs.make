# 
#
SIESTA_ARCH=nolibs
#
# Serial compilation without the need of any installed libraries.
# You still need to change the name of the compiler, etc.
#
FC=f90
FC_ASIS=$(FC)
#
FFLAGS= -O
FFLAGS_DEBUG= -g -O0
RANLIB=echo "do we need to ranlib this? : "
LIBS=  
SYS=nag
FPPFLAGS=
MPILIB=
COMP_LIBS=linalg.a
MPI_INCLUDE=/usr/local/include
#
.F.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS)  $(FPPFLAGS) $<
.f.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS)   $<
.F90.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS)  $(FPPFLAGS) $<
.f90.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS)   $<
#
