##
##
## Introduction
## ============
##
## Makefile for flink/src directory.
##
## Usage
## =====
##
## 'make' or 'make clean'
##
## Author
## ======
##
## This building system is designed, created, and maintained by
##
## Li Huang // email: lihuang.dmft@gmail.com
##
## History
## =======
##
## 09/18/2014 by li huang (created)
## 08/06/2021 by li huang (last modified)
##
##

.SUFFIXES: .f90

# Fortran compiler, linker, and archiver
#-------------------------------------------------------------------------
F90    = mpiifort
LINKER = $(F90)
ARCHIVER = ar -ruv

# Fortran preprocessor options
#-------------------------------------------------------------------------
MPI    = -DMPI
OMP    = #-qopenmp
FPP    = -fpp
CPP    = $(FPP) $(MPI) $(OMP)

# Machine tuning options
#-------------------------------------------------------------------------
CHECK  = -nogen-interfaces -warn all #-check all -traceback -g
MTUNE  = -O3 -xHost

# Flags for fortran compiler and linker
#-------------------------------------------------------------------------
FFLAGS = -c $(CPP) $(CHECK) $(MTUNE)
LFLAGS = $(OMP)

# External linear algebra library
#-------------------------------------------------------------------------
LIBS   = -L${MKLROOT}/lib -lmkl_core -lmkl_sequential -lmkl_rt

mods = m_constants.o m_linkedlist.o m_mpi.o m_parser.o m_sparse.o m_spring.o m_stack.o m_tetra.o
subs = s_error.o s_fourier.o s_function.o s_integrator.o s_matrix.o s_spline.o s_util.o s_vector.o
objects = $(mods) $(subs)

default: all

all: lib

lib: $(objects)
	$(ARCHIVER) libflink.a $(objects)

.f90.o:
	$(F90) $(FFLAGS) $*.f90

clean:
	rm -f *.mod
	rm -f *.o
	rm -f libflink.a

clean-dat:
	rm -f *.dat
	rm -f *.out

clean-all: clean clean-dat
