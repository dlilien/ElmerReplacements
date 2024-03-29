#
# Makefile
# dlilien, 2016-05-03 13:54
#
#
# vim:ft=make
#
#
SHELL=/bin/bash
FCB=mpifort -fPIC
FC=elmerf90
FCNS=./elmerf90-nosh
CC=mpicc
LIB_SOURCES=borstad_damage.f90 DJDmu_Adjoint_lilien.f90 MATC_Replacements.f90 DAViscosityInversion.f90 g2di.f90 Cost_Functions.f90 DummySolver.f90 USF_Sliding.f90 USF_Contact.f90 USF_GetFrictionHeating.f90 lilien_sliding.f90 USF_Zs.f90 AIFlowFricHeat.f90 sidewall_drag.f90 MF.f90
SOLVER_SOURCES=melt_solver.f90 manual_grounding.f90 GroundedSolverSSA.f90 FabricSolve.f90
MODULE_SOURCES=read_routines.f90 MeltFunctions.f90
MODULE=$(MODULE_SOURCES:.f90=.o)
LIB_OBJECTS=$(LIB_SOURCES:.f90=.o)
SOLVER_OBJECTS=$(SOLVER_SOURCES:.f90=.o)
SOLVER_SINGLES=$(SOLVER_SOURCES:.f90=.so)
LIB=lilien_lib.so
SOLVER=lilien_solver.so

FORT_SOURCES=MshGlacierDEM.f90
FORT_OBJECTS=$(FORT_SOURCES:.f90=)

C_SOURCES=ExtrudeMesh.c
C_OBJECTS=$(C_SOURCES:.c=)

all: compile AdjointSSASolvers

compile: $(MODULE) $(LIB_OBJECTS) $(SOLVER_SINGLES) $(LIB) $(C_OBJECTS) $(SOLVER)
	
AdjointSSASolvers:
	$(MAKE) -C src/AdjointSSA compile

$(LIB): $(MODULE) $(LIB_OBJECTS)
	$(FC) $(LIB_OBJECTS) $(SOLVER_OBJECTS) $(MODULE_SOURCES) -I . -o $@

$(SOLVER): $(SOLVER_OBJECTS)
	$(FC) $(SOLVER_OBJECTS) -o $@

$(MODULE): %.o: %.f90
	@ if [[ `hostname` =~ g-* ]] ; then $(FC) -c $< -o $@; else $(FC) -c $< -o $@ -lm ; fi

$(LIB_OBJECTS): %.o: %.f90
	$(FC) -c $< -o $@

$(SOLVER_OBJECTS): %.o: %.f90
	$(FC) -c $< -o $@

$(SOLVER_SINGLES): %.so: %.o
	$(FC) $< -o $@

$(FORT_OBJECTS): %: %.f90
	./elmerf90-nosh $< -o $@

$(C_OBJECTS): %: %.c
	@ if [[ `hostname` =~ pfe* ]] ; then icc $< -o $@; else $(CC) $< -o $@ -lm ; fi

test: testMelt testRead
	./testMelt
	./testRead

testMelt: $(LIB) testMelt.f90
	$(FCNS) -I./ MeltFunctions.f90 -o testMelt testMelt.f90

testRead: $(LIB) testRead.f90
	$(FCNS) --I./ read_routines.f90 -o testRead testRead.f90

clean:
	-rm -f *.o *.so MshGlacierDEM ExtrudeMesh *.mod
	$(MAKE) -C src/AdjointSSA clean
