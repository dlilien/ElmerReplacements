#
# Makefile
# dlilien, 2016-05-03 13:54
#
#
# vim:ft=make
#
#
SHELL=/bin/bash
FCB=ifort -fPIC
FC=elmerf90
FCNS=./elmerf90-nosh
CC=/usr/bin/gcc
LIB_SOURCES=borstad_damage.f90 DJDmu_Adjoint_lilien.f90 MATC_Replacements.f90 DAViscosityInversion.f90 g2di.f90 Cost_Functions.f90 DummySolver.f90 GroundedSolverSSA.f90 USF_Sliding.f90 USF_Contact.f90 USF_GetFrictionHeating.f90 lilien_sliding.f90
SOLVER_SOURCES=melt_solver.f90 manual_grounding.f90
MODULE_SOURCES=read_routines.f90 MeltFunctions.f90
MODULE=$(MODULE_SOURCES:.f90=.o)
LIB_OBJECTS=$(LIB_SOURCES:.f90=.o)
SOLVER_OBJECTS=$(SOLVER_SOURCES:.f90=.o)
LIB=lilien_lib.so
SOLVER=lilien_solver.so

FORT_SOURCES=MshGlacierDEM.f90
FORT_OBJECTS=$(FORT_SOURCES:.f90=)

C_SOURCES=ExtrudeMesh.c
C_OBJECTS=$(C_SOURCES:.c=)

all: compile AdjointSSASolvers

compile: $(MODULE) $(LIB_OBJECTS) $(LIB) $(C_OBJECTS) $(SOLVER)
	
AdjointSSASolvers:
	$(MAKE) -C src/AdjointSSA compile

$(LIB): $(MODULE) $(LIB_OBJECTS)
	$(FC) $(LIB_OBJECTS) $(MODULE) -I . -o $@

$(SOLVER): $(SOLVER_OBJECTS)
	$(FC) $(SOLVER_OBJECTS) -o $@

$(MODULE): %.o: %.f90
	@ if [[ `hostname` =~ pfe* ]] ; then $(FC) -assume byterecl -c $< -o $@; else $(FC) -c $< -o $@ -lm ; fi

$(LIB_OBJECTS): %.o: %.f90
	$(FC) -c $< -o $@

$(SOLVER_OBJECTS): %.o: %.f90
	$(FC) -c $< -o $@

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
	$(FCNS) -assume byterecl --I./ read_routines.f90 -o testRead testRead.f90

clean:
	-rm -f *.o *.so MshGlacierDEM ExtrudeMesh
	$(MAKE) -C src/AdjointSSA clean
