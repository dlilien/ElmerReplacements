#
# Makefile
# dlilien, 2016-05-03 13:54
#
#

FC=elmerf90
LIB_SOURCES=borstad_damage.f90 DJDmu_Adjoint_lilien.f90 MATC_Replacements.f90 lilien_sliding.f90 DAViscosityInversion.f90 g2di.f90
LIB_OBJECTS=$(LIB_SOURCES:.f90=.o)
LIB=lilien_lib.so

FORT_SOURCES=MshGlacierDEM.f90
FORT_OBJECTS=$(FORT_SOURCES:.f90=)

C_SOURCES=ExtrudeMesh.c
C_OBJECTS=$(C_SOURCES:.c=)

all: $(LIB_OBJECTS) $(LIB)

$(LIB): $(LIB_OBJECTS)
	$(FC) $(LIB_OBJECTS) -o $@

$(LIB_OBJECTS): %.o: %.f90
	$(FC) -c $< -o $@

$(FORT_OBJECTS): %: %.f90
	elmerf90-nosh $< -o $@

$(C_OBJECTS): %: %.c
	$(CC) $< -o $@

#MshGlacierDEM: MshGlacierDEM.f90
#	elmerf90-nosh -o MshGlacierDEM MshGlacierDEM.f90

#ExtrudeMesh: ExtrudeMesh.c
#    if [[ `hostname` =~ pfe* ]]; then
#        module load gcc
#    fi
#	gcc -o ExtrudeMesh ExtrudeMesh.c

clean:
	rm -f *.o *.so MshGlacierDEM ExtrudeMesh


# vim:ft=make
#