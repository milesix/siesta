# 
#
SIESTA_ARCH=darwin
#
FC=f90
FC_ASIS=$(FC)
#
FFLAGS=  -O
FFLAGS_DEBUG= -g
RANLIB=echo "... do we need to ranlib this? :"
LDFLAGS=
COMP_LIBS=linalg.a
#
NETCDF_LIBS=
NETCDF_INTERFACE=
FPPFLAGS_CDF=
#
MPI_INTERFACE=
MPI_INCLUDE=
FPPFLAGS_MPI=
#
LIBS= -lblas -lU77 $(NETCDF_LIBS)
SYS=bsd
FPPFLAGS= $(FPPFLAGS_CDF) $(FPPFLAGS_MPI)
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
