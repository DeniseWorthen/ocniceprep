# GNU Makefile template for user ESMF application

################################################################################
################################################################################
## This Makefile must be able to find the "esmf.mk" Makefile fragment in the  ##
## 'include' line below. Following the ESMF User's Guide, a complete ESMF     ##
## installation should ensure that a single environment variable "ESMFMKFILE" ##
## is made available on the system. This variable should point to the         ##
## "esmf.mk" file.                                                            ##
##                                                                            ##
## This example Makefile uses the "ESMFMKFILE" environment variable.          ##
##                                                                            ##
## If you notice that this Makefile cannot find variable ESMFMKFILE then      ##
## please contact the person responsible for the ESMF installation on your    ##
## system.                                                                    ##
## As a work-around you can simply hardcode the path to "esmf.mk" in the      ##
## include line below. However, doing so will render this Makefile a lot less ##
## flexible and non-portable.                                                 ##
################################################################################

ifneq ($(origin ESMFMKFILE), environment)
$(error Environment variable ESMFMKFILE was not set.)
endif

include $(ESMFMKFILE)

#localFopt = -C -fPIC

#localFopt = -check bounds -fpe0 -ftrapuv -fPIC -warn unused -check uninit
#localFopt = -check bounds -fpe0 -ftrapuv -fPIC -check uninit
localFopt = -check bounds -fpe0 -ftrapuv -fPIC -init=snan,arrays -check uninit

#localFopt = -Wunused
################################################################################
################################################################################

.SUFFIXES: .f90 .F90 .c .C

%.o : %.f90
	$(ESMF_F90COMPILER) -c $(ESMF_F90COMPILEOPTS) $(ESMF_F90COMPILEPATHS) $(ESMF_F90COMPILEFREENOCPP) $<

%.o : %.F90
	$(ESMF_F90COMPILER) -c $(localFopt) $(ESMF_F90COMPILEOPTS) $(ESMF_F90COMPILEPATHS) $(ESMF_F90COMPILEFREECPP) $(ESMF_F90COMPILECPPFLAGS) -DESMF_VERSION_MAJOR=$(ESMF_VERSION_MAJOR) $<

%.o : %.c
	$(ESMF_CXXCOMPILER) -c $(ESMF_CXXCOMPILEOPTS) $(ESMF_CXXCOMPILEPATHSLOCAL) $(ESMF_CXXCOMPILEPATHS) $(ESMF_CXXCOMPILECPPFLAGS) $<

%.o : %.C
	$(ESMF_CXXCOMPILER) -c $(ESMF_CXXCOMPILEOPTS) $(ESMF_CXXCOMPILEPATHSLOCAL) $(ESMF_CXXCOMPILEPATHS) $(ESMF_CXXCOMPILECPPFLAGS) $<


# -----------------------------------------------------------------------------
ictest: ocniceprep.o utils_mod.o utils_esmf_mod.o init_mod.o arrays_mod.o restarts_mod.o
	$(ESMF_F90LINKER) $(ESMF_F90LINKOPTS) $(ESMF_F90LINKPATHS) $(ESMF_F90LINKRPATHS) -o $@ $^ $(ESMF_F90ESMFLINKLIBS)

# module dependencies:
ocniceprep.o: utils_mod.o utils_esmf_mod.o init_mod.o arrays_mod.o restarts_mod.o
restarts_mod.o: utils_mod.o init_mod.o arrays_mod.o
utils_mod.o: init_mod.o
utils_esmf_mod.o: init_mod.o arrays_mod.o
init_mod.o:
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
.PHONY: dust clean distclean info edit
dust:
	rm -f PET*.ESMF_LogFile *.nc *.stdout
clean:
	rm -f ictest *.o *.mod
distclean: dust clean

info:
	@echo ==================================================================
	@echo ESMFMKFILE=$(ESMFMKFILE)
	@echo ==================================================================
	@cat $(ESMFMKFILE)
	@echo ==================================================================

edit:
	nedit rhtest.F90

run:
	mpirun -np 4 ./rhtest
