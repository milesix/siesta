# 
# Copyright (C) 1996-2016	The SIESTA group
#  This file is distributed under the terms of the
#  GNU General Public License: see COPYING in the top directory
#  or http://www.gnu.org/copyleft/gpl.txt.
# See Docs/Contributors.txt for a list of contributors.
#
# Makefile for Vibra package

.SUFFIXES:
.SUFFIXES: .f .F .o .a .f90 .F90

.PHONY: default dep clean

PROGS:= fcbuild vibra
default: $(PROGS)

override WITH_MPI=
override WITH_NETCDF=
override WITH_FLOOK=

MAIN_OBJDIR=.
TOPDIR=.

ARCH_MAKE=$(MAIN_OBJDIR)/arch.make
include $(ARCH_MAKE)
include $(MAIN_OBJDIR)/check_for_build_mk.mk

VPATH=$(TOPDIR)/Util/Vibra/Src:$(TOPDIR)/Src

FC_DEFAULT:=$(FC)
FC_SERIAL?=$(FC_DEFAULT)
FC:=$(FC_SERIAL)         # Make it non-recursive


# Define these if desired
#
### FPPFLAGS += -DSANKEY_DIAG
### FPPFLAGS += -DVIBRA_DEBUG

# Objects from Siesta

SIESTA_SRCS = precision.F m_io.f io.f reclat.f 
SIESTA_OBJS = $(call change_extensions, $(SIESTA_SRCS))

COMMON_SRCS = recoor.f local_sys.f
COMMON_OBJS = recoor.o local_sys.o

BUILD_SRCS = fcbuild.f 
BUILD_OBJS = fcbuild.o $(COMMON_OBJS)

VIBRA_SRCS = vibra.F hermdp.F klines.f outbands.f
VIBRA_OBJS = vibra.o hermdp.o klines.o outbands.o $(COMMON_OBJS)

INCFLAGS+= $(FDF_INCFLAGS) $(UNITS_INCFLAGS)


fcbuild: $(SIESTA_OBJS) $(BUILD_OBJS) $(UNITS)
	$(FC) -o fcbuild \
	       $(LDFLAGS) $(BUILD_OBJS) $(SIESTA_OBJS) $(FDF_LIBS) $(UNITS)

vibra:  $(COMP_LIBS) $(SIESTA_OBJS) $(VIBRA_OBJS) $(UNITS)
	$(FC) -o vibra \
	       $(LDFLAGS) $(VIBRA_OBJS) $(SIESTA_OBJS) $(FDF_LIBS) $(UNITS) $(COMP_LIBS) $(LIBS)

clean: 
	@echo "==> Cleaning object and executable files"
	rm -f fcbuild vibra *.o *.mod

install: $(PROGS)
	cp -p $(PROGS) $(SIESTA_INSTALL_DIRECTORY)/bin


# Dependencies
ABS_TOP_SRCS:= $(addprefix $(TOPDIR)/Src/,$(SIESTA_SRCS))
ABS_LOCAL_SRCS:= $(addprefix $(TOPDIR)/Util/Vibra/Src/,$(COMMON_SRCS) $(VIBRA_SRCS) $(BUILD_SRCS))

.PHONY: dep
dep:
	sfmakedepend --depend=obj --modext=o $(ABS_TOP_SRCS) $(ABS_LOCAL_SRCS) || true


# DO NOT DELETE THIS LINE - used by make depend
io.o: m_io.o
m_io.o: local_sys.o
fcbuild.o: precision.o
vibra.o: precision.o
sys.o: local_sys.o
