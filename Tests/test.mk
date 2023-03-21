#
# Single-test makefile template
#
# You can edit the SIESTA macro here, or pass it on the command line

TOPDIR=/Users/ag/code/GITLAB/Siesta/garalb
MAIN_OBJDIR=/Users/ag/code/GITLAB/Siesta/garalb/,pdos

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
completed: completed_$(label)

completed_$(label):
	@echo ">>>> Running $(name) test..."
	@if [ -d $(label) ] ; then rm -rf $(label) ; fi; mkdir $(label)
	@if [ -n "$(EXTRAFILES)" ] ; then cp -f $(EXTRAFILES) $(label) ; fi
	@echo "    ==> Running SIESTA as SIESTA_PS_PATH=../../Pseudos $(MPI) $(SIESTA) -fdf XML.Write ../$(name).fdf "
	@(cd $(label) ; SIESTA_PS_PATH=../../Pseudos $(MPI) $(SIESTA) -fdf XML.Write ../$(name).fdf 2>&1 > $(name).out ) \
          && touch completed_$(label)
	@if [ -f completed_$(label) ] ; then cp $(label)/$(name).out $(label)/$(name).xml .;\
           echo "    ===> SIESTA finished successfully";\
         else \
           echo " **** Test $(name) did not complete successfully";\
         fi

check: completed check-only

check-only:
	@echo "    ==> Running check for system $(name)"
	@REFERENCE_DIR=$(REFERENCE_DIR) sh $(REFERENCE_CHECKER) $(name).out

clean:
	@echo ">>>> Cleaning $(name) [label=$(label)] test..."
	rm -rf $(label) completed_$(label) $(name).out $(name).xml
