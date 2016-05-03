#
# Makefile
# dlilien, 2016-05-03 13:54
#

all: DJDmu_Adjoint_lilien.so MATC_Replacements.so lilien_sliding.so DAViscosityInversion.so MshGlacierDEM ExtrudeMesh
	elmerf90 DJDmu_Adjoint_lilien.so MATC_Replacements.so lilien_sliding.so -o lilien_lib.so

DJDmu_Adjoint_lilien.so: DJDmu_Adjoint_lilien.f90
	elmerf90 -o DJDmu_Adjoint_lilien.so DJDmu_Adjoint_lilien.f90

MATC_Replacements.so: MATC_Replacements.f90
	elmerf90 -o MATC_Replacements.so MATC_Replacements.f90

lilien_sliding.so: lilien_sliding.f90
	elmerf90 -o lilien_sliding.so lilien_sliding.f90

DAViscosityInversion.so: DAViscosityInversion.F90
	elmerf90 -o DAViscosityInversion.so DAViscosityInversion.F90

MshGlacierDEM: MshGlacierDEM.f90
	elmerf90-nosh -o MshGlacierDEM MshGlacierDEM.f90

ExtrudeMesh: ExtrudeMesh.c
	$(CC) -o ExtrudeMesh ExtrudeMesh.c




# vim:ft=make
#
