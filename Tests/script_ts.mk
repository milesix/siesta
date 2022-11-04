#
# Single-test makefile template for script usage
#
# We do not have popd and pushd in "sh" scripts, thus force the SHELL to be bash
SHELL=/bin/bash

TOPDIR=.
MAIN_OBJDIR=.

MPI=mpirun -np 2
SIESTA=$(MAIN_OBJDIR)/Src/siesta

# Make compatibility layer for old test-runs
ifeq ($(strip $(firstword $(TS))),"mpirun")
MPI=
endif
ifeq ($(strip $(firstword $(TS))),"mpiexec")
MPI=
endif

label=work

.PHONY: completed
completed: completed_$(label)

completed_$(label):
	@echo ">>>> Running $(name) test..."
	@if [ -d $(label) ] ; then rm -rf $(label) ; fi; mkdir $(label)
	@if [ -f script.sh ] ; then cp -f script.sh $(label) ; fi
	@echo "    ==> Running script with TranSIESTA as $(TS)"
	@(cd $(label) ; $(SHELL) script.sh "$(MPI) $(TS)" )
	@if [ -f completed ] ; then \
           echo "    ===> Script finished successfully";\
         else \
           echo " **** Test $(name) did not complete successfully";\
         fi

clean:
	@echo ">>>> Cleaning $(name) [label=$(label)] test..."
	rm -rf $(label) completed_$(label) $(name).out $(name).xml
