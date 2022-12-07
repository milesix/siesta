# ---
# Copyright (C) 1996-2016	The SIESTA group
#  This file is distributed under the terms of the
#  GNU General Public License: see COPYING in the top directory
#  or http://www.gnu.org/copyleft/gpl.txt .
# See Docs/Contributors.txt for a list of contributors.
# ---
.SUFFIXES:
.SUFFIXES: .f .f90 .F .F90 .o
#
TOPDIR=.
MAIN_OBJDIR=.

VPATH=$(TOPDIR)/Util/JobList/Src:$(TOPDIR)/Src

PROGS= countJobs runJobs getResults horizontal
default: $(PROGS)

include $(MAIN_OBJDIR)/arch.make
include $(MAIN_OBJDIR)/check_for_build_mk.mk

FC:=$(FC_SERIAL)         # Make it non-recursive

#
COUNTJOBS_OBJS=jobList.o countJobs.o posix_calls.o
RUNJOBS_OBJS=jobList.o runJobs.o posix_calls.o
GETRESULTS_OBJS=jobList.o getResults.o posix_calls.o
HORIZONTAL_OBJS=jobList.o horizontal.o posix_calls.o
#
dep: 
	sfmakedepend --depend=obj --modext=o *.f *.f90 *.F *.F90
#
countJobs: $(COUNTJOBS_OBJS)
	$(FC) -o $@ $(COUNTJOBS_OBJS)
#
runJobs: $(RUNJOBS_OBJS)
	$(FC) -o $@ $(RUNJOBS_OBJS)
#
getResults: $(GETRESULTS_OBJS)
	$(FC) -o $@ $(GETRESULTS_OBJS)
#
horizontal: $(HORIZONTAL_OBJS)
	$(FC) -o $@ $(HORIZONTAL_OBJS)
#
clean:
	rm -f *.o *.mod $(PROGS)
#
install: $(PROGS)
	cp -p $(PROGS) $(SIESTA_INSTALL_DIRECTORY)/bin

# DO NOT DELETE THIS LINE - used by make depend
countJobs.o: jobList.o
getResults.o: jobList.o
jobList.o: posix_calls.o
runJobs.o: jobList.o
