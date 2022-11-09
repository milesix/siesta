#
# test makefile template
# for creating different tests with different
# list elements
#
# This makefile implements a simplistic way of
# running a list of elements.
#
# The list of elements must be provided in the
# variable:
#   LIST
# while the fdf-flag to be runned should be named
#   FDF_LIST


TOPDIR=.
MAIN_OBJDIR=.

MPI=mpirun -np 2
SIESTA=$(MAIN_OBJDIR)/Src/siesta

# Make compatibility layer for old test-runs
ifeq ($(strip $(firstword $(SIESTA))),mpirun)
MPI=
endif
ifeq ($(strip $(firstword $(SIESTA))),mpiexec)
MPI=
endif

#----------------------------------------------------------------------------
REFERENCE_DIR?=$(TOPDIR)/Tests/Reference
REFERENCE_CHECKER?=../cmp_digest.sh

label=work


.PHONY: completed

# Create list of jobs
_JOBS = $(addprefix completed_$(label)_,$(LIST))

# Create completed jobs target
completed: $(_JOBS)

.PHONY: $(_JOBS)

$(_JOBS): completed_$(label)_%:
	@echo ">>>> Running $(name)_$* test..."
	@if [ -d $(label)_$* ] ; then rm -rf $(label)_$* ; fi; mkdir $(label)_$*
	@if [ -n "$(EXTRAFILES)" ] ; then cp -f $(EXTRAFILES) $(label)_$* ; fi
	@for ps in `cat $(name).pseudos` ; do \
          echo "    ==> Copying pseudopotential file $$ps ..." ;\
          ln ../Pseudos/$$ps $(label)_$*/$$ps ;\
         done
	@echo "    ==> Running SIESTA as $(MPI) $(SIESTA) -L $(name)_$* -fdf XML.Write -fdf $(FDF_LIST):$* ../$(name).fdf"
	@(cd $(label)_$* ; $(MPI) $(SIESTA) -L $(name)_$* -fdf XML.Write -fdf $(FDF_LIST):$* ../$(name).fdf 2>&1 > $(name)_$*.out ) \
          && touch completed_$(label)_$*
	@if [ -f completed_$(label)_$* ] ; then cp $(label)_$*/$(name)_$*.out $(label)_$*/$(name)_$*.xml .;\
           echo "    ===> SIESTA finished successfully";\
         else \
           echo " **** Test $(name)_$* did not complete successfully";\
         fi

check: completed check-only

# Create list of jobs
_CHECK_JOBS = $(addprefix check-only-,$(LIST))

check-only: $(_CHECK_JOBS)
$(_CHECK_JOBS): check-only-%:
	@echo "    ==> Running check for system $(name)_$*"
	@REFERENCE_DIR=$(REFERENCE_DIR) sh $(REFERENCE_CHECKER) $(name)_$*.out

_CLEAN_JOBS = $(addprefix clean-,$(LIST))

clean: $(_CLEAN_JOBS)
$(_CLEAN_JOBS): clean-%:
	@echo ">>>> Cleaning $(name)_$* [label=$(label)] test..."
	rm -rf $(label)_$* completed_$(label)_$* $(name)_$*.out $(name)_$*.xml

