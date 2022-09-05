#!/bin/bash

#============Parsing install.var=======
#============Compilers===============
CXX=`grep "CXX" install.vars | grep -v "MPI_CXX" | awk '{sub(/^[ \t]+/,"");print $3}'`        #C++ compiler
if [ -z CXX ];then
  echo "ERROR: CXX compiler is not specified" 
  exit
fi

CC=`grep "CC" install.vars | grep -v "MPI_CC" | awk '{sub(/^[ \t]+/,"");print $3}'`           #C complier
if [ -z CC ];then
  echo "ERROR: CC compiler is not specified" 
  exit
fi

FC=`grep "FC" install.vars | grep -v "MPI_FC" | awk '{sub(/^[ \t]+/,"");print $3}'`           #Fortran compiler
if [ -z FC ];then
  echo "ERROR: FC compiler is not specified" 
  exit
fi

MPI_FC=`grep "MPI_FC" install.vars | awk '{sub(/^[ \t]+/,"");print $3}'`                      #MPI version fortran compiler
if [ -z MPI_FC ];then
  echo "ERROR: MPI_FC compiler is not specified" 
  exit
fi

MPI_CXX=`grep "MPI_CXX" install.vars | awk '{sub(/^[ \t]+/,"");print $3}'`                    #MPI version C++ compiler
if [ -z MPI_CXX ];then
  echo "ERROR: MPI_CXX compiler is not specified" 
  exit
fi

MPI_CC=`grep "MPI_CC" install.vars | awk '{sub(/^[ \t]+/,"");print $3}'`                      #MPI version C compiler
if [ -z MPI_CC ];then
  echo "ERROR: MPI_CC compiler is not specified" 
  exit
fi

#============DFT softwares path=========
FHIaims_lib_path=`grep "FHIaims_lib_path" install.vars | awk '{sub(/^[ \t]+/,"");print $3}'`
ABACUS_lib_path=`grep "ABACUS_lib_path" install.vars | awk '{sub(/^[ \t]+/,"");print $3}'`

#========Starting compilation========
cd ../
root_dir=$(dirname $(readlink -f "$0"))
cd build

#====================================
#    PART 1: PACS CTHYB
#====================================
test -d impurities || mkdir impurities
cd impurities 
test -d PACS || mkdir PACS
cd PACS

if [ ! -f libpacs.a ]
then
cat > Makefile << EOF
FC         = $MPI_FC
FC_FLAGS   = -O3 -xHost -nogen-interface # -ipo -traceback -implicitnone -warn all -check bounds   
#FC        = mpif90
#FC_FLAGS  = -Wall -Wextra -pedantic -Wconversion -fbacktrace -fbounds-check -ffree-line-length-none
LIB        =  -L\${MKLROOT}/lib/intel64 -lmkl_rt -lpthread -lm -ldl

VPATH=$root_dir/src/impurities/PACS

.suffixes : .mod .o .f90

%.o:  %.f90
	\$(FC) \$(FC_FLAGS) -c \$<

SRC = GlobalVariables.f90 \\
      ctqmc_math.f90 MPI_mod.f90 \\
      Segment_Util.f90 Segment_MonteCarlo.f90 Segment_Phys.f90 \\
      Segment_LibraryMode.f90  Input_Parameters.f90 main.f90 

OBJ = \$(SRC:.f90=.o)

intel: \$(OBJ)
	ar -rc libpacs.a \$(OBJ)
clean:
	rm -f *.o *.mod  *.a
EOF
  make
  if [ $? -ne 0 ]; then
    echo "Errors in building PACS-CTHYB"
    exit
  fi
fi
cd $root_dir/build

#====================================
#    PART 2: iQIST CTHYB
#====================================
cd impurities 

test -d iQIST || mkdir iQIST
cd iQIST

test -d atomic || mkdir atomic
test -d ct_hyb1 || mkdir ct_hyb1
test -d ct_hyb2 || mkdir ct_hyb2
test -d dependencies || mkdir dependencies

cd dependencies
test -d Flink || mkdir Flink
cd Flink

if [ ! -f libflink.a ]
then
cat > Makefile << EOF
.SUFFIXES: .f90

# Fortran compiler, linker, and archiver
#-------------------------------------------------------------------------
F90    = $MPI_FC
LINKER = \$(F90)
ARCHIVER = ar -ruv

# Fortran preprocessor options
#-------------------------------------------------------------------------
MPI    = -DMPI
OMP    = #-qopenmp
FPP    = -fpp
CPP    = \$(FPP) \$(MPI) \$(OMP)

# Machine tuning options
#-------------------------------------------------------------------------
CHECK  = -nogen-interfaces -warn all #-check all -traceback -g
MTUNE  = -O3 -xHost

# Flags for fortran compiler and linker
#-------------------------------------------------------------------------
FFLAGS = \$(CPP) \$(CHECK) \$(MTUNE)
LFLAGS = \$(OMP)

# External linear algebra library
#-------------------------------------------------------------------------
LIBS   = -L\${MKLROOT}/lib -lmkl_core -lmkl_sequential -lmkl_rt

VPATH=${root_dir}/src/impurities/iQIST/dependencies/Flink/src

mods = m_constants.o m_linkedlist.o m_mpi.o m_parser.o m_sparse.o m_spring.o m_stack.o m_tetra.o
subs = s_error.o s_fourier.o s_function.o s_integrator.o s_matrix.o s_spline.o s_util.o s_vector.o
objects = \$(mods) \$(subs)

default: all

all: lib

lib: \$(objects)
	\$(ARCHIVER) libflink.a \$(objects)

.f90.o:
	\$(F90) \$(FFLAGS) -c \$<

clean:
	rm -f *.mod
	rm -f *.o
	rm -f libflink.a

clean-dat:
	rm -f *.dat
	rm -f *.out

clean-all: clean clean-dat
EOF
  make
  if [ $? -ne 0 ]; then
    echo "Errors in building Flink library"
    exit
  fi
fi
cd $root_dir/build

cd $root_dir/build/impurities/iQIST/ct_hyb1
if [ ! -f libnarcissus.a ]
then
cat > Makefile << EOF
.SUFFIXES: .f90

# Fortran compiler, linker, and archiver
#-------------------------------------------------------------------------
F90    = $MPI_FC
LINKER = \$(F90)
ARCHIVER = ar -ruv

# Fortran preprocessor options
#-------------------------------------------------------------------------
MPI    = -DMPI
OMP    = #-qopenmp
FPP    = -fpp
CPP    = \$(FPP) \$(MPI) \$(OMP)

# Machine tuning options
#-------------------------------------------------------------------------
CHECK  = -nogen-interfaces -warn all #-check all -traceback -g
MTUNE  = -O3 -xHost

# Flags for fortran compiler and linker
#-------------------------------------------------------------------------
FFLAGS = \$(CPP) \$(CHECK) \$(MTUNE)
LFLAGS = \$(OMP)

# External linear algebra library
#-------------------------------------------------------------------------
LIBS   = -L\${MKLROOT}/lib -lmkl_core -lmkl_sequential -lmkl_rt
FLINK  = $root_dir/build/impurities/iQIST/dependencies/Flink

VPATH=${root_dir}/src/impurities/iQIST/src/ct_hyb1

modc = ctqmc_control.o ctqmc_context.o
dmft = ctqmc_dmft.o
core = ctqmc_solver.o
lev1 = ctqmc_flavor.o ctqmc_hybmat.o ctqmc_update.o
lev2 = ctqmc_record.o ctqmc_status.o ctqmc_stream.o ctqmc_util.o
lev3 = ctqmc_dump.o ctqmc_print.o
main = ctqmc_main.o
mlib = libflink.a

objects = \$(modc) \$(dmft) \$(core) \$(lev1) \$(lev2) \$(lev3) \$(main) #\$(mlib)

default: all

all: lib

flink: flink_lib flink_mod

flink_lib:
	cp \$(FLINK)/libflink.a .

flink_mod:
	cp \$(FLINK)/constants.mod .
	cp \$(FLINK)/mmpi.mod .
	cp \$(FLINK)/spring.mod .
	cp \$(FLINK)/stack.mod .
	cp \$(FLINK)/parser.mod .
	cp \$(FLINK)/linkedlist.mod .

lib: flink \$(objects)
	\$(ARCHIVER) libnarcissus.a \$(objects)

.f90.o:
	\$(F90) \$(FFLAGS) -c \$<

clean:
	rm -f *.mod
	rm -f *.o
	rm -f ctqmc
	rm -f libflink.a
	rm -f libnarcissus.a

clean-dat:
	rm -f *.dat
	rm -f *.out

clean-all: clean clean-dat
EOF
  make all

  if [ $? -ne 0 ]; then
    echo "Errors in building iQIST-CTHYB"
    exit
  fi
fi
cd $root_dir/build

#====================================
#    PART 3: Rutgers CTHYB
#====================================
cd impurities

test -d Rutgers || mkdir Rutgers
cd Rutgers

test -d dependencies || mkdir dependencies
cd dependencies

cp $root_dir/src/impurities/Rutgers/dependencies/gsl-2.6.tar.gz .

if [ ! -d gsl-2.6 ];then
  tar -zxvf gsl-2.6.tar.gz
fi

rm gsl-2.6.tar.gz

cd gsl-2.6
if [ ! -f libgsl.la ];then
  ./configure CC=$MPI_CC --prefix=$root_dir/build/impurities/Rutgers/dependencies/gsl-2.6

  if [ $? -eq 0 ]
  then
    make -j
    if [ $? -eq 0 ]
    then
      make install
    else 
      echo "Errors in building gsl-2.6"
      exit
    fi
  else
    echo "Errors in building gsl-2.6"
    exit
  fi
fi
# Check the lib folder name for good
if [ ! -d lib ] && [ -d lib64 ]; then
  ln -s lib64 lib
fi

cd $root_dir/build

cd impurities/Rutgers
if [ ! -f librutgers.a ]
then 
  cp $root_dir/src/impurities/Rutgers/src/* .

cat > Makefile << EOF
PC++ = $MPI_CXX
F77 = $FC

GSLINC = -I$root_dir/build/impurities/Rutgers/dependencies/gsl-2.6/include
GSLLIB =  -L$root_dir/build/impurities/Rutgers/dependencies/gsl-2.6/lib -l:libgslcblas.a -l:libgsl.a

LLIBS  =  -mkl
PLIBS = \$(LLIBS) \$(GSLLIB)

PFLAGS   = -D_MPI -DMPICH_IGNORE_CXX_SEEK -O3 #-restrict -ipo -no-prec-div
CFLAGS = \$(PFLAGS) \$(GSLINC) -D_TIME #-DAS -D_TIME #-fast -xAVX # -D_TIME #-D_LOGGING #-DAS 
FFLAGS = -O2  -free -no-prec-div -pc80 -qopenmp 

GHEADERS = assert.h complex.h random.h sblas.h sfunction.h smesh.h sutil.h zeroin.h
QHEADERS = common.h  inout.h intervals.h local.h matrixm.h mpi.h bcast.h number.h operators.h state.h stateim.h segment.h svdfunc.h tanmesh.h

ctqmc : ctqmc.o SMatrix1.o 
	ar -rc librutgers.a ctqmc.o SMatrix1.o

all :  ctqmc ctqmcf

ctqmcf : ctqmcf.o SMatrix1.o
	ar -rc librutgersf.a ctqmcf.o SMatrix1.o

ctqmcf.o : ctqmc.cc
	\$(PC++) \$(CFLAGS) -DAS -c -o ctqmcf.o \$<

SMatrix1.o : SMatrix1.cc sfunction.h
	\$(PC++) -c \$(CFLAGS) SMatrix1.cc 

ctqmc.o : ctqmc.cc \$(GHEADERS) \$(QHEADERS) 
	\$(PC++) \$(CFLAGS) -c ctqmc.cc

clean :
	rm -f ctqmc.o ctqmcf.o ctqmcf ctqmc SMatrix1.o librutgers*

.SUFFIXES : .cc
.cc.o:
	\$(PC++) \$(CFLAGS) -c \$<

.SUFFIXES : .f
.f.o:
	\$(F77) \$(FFLAGS) -c \$<
EOF

  make all
  if [ $? -ne 0 ]; then
    echo "Errors in building Rutgers-CTHYB"
    exit
  fi
fi
cd $root_dir/build

#====================================
#    PART 4: DFT+DMFT
#====================================
# if [ ! -f ../bin/DFTDMFT ]
if [ ! -f DFTDMFT ]
then
  if [ -z $ABACUS_lib_path ];then
    MACRO_ABACUS=
  else
    MACRO_ABACUS=-D__ABACUS
  fi

  if [ -z $FHIaims_lib_path ];then    #FHI-aims has not been built
cat > Makefile <<EOF
CPLUSPLUS_MPI = $MPI_CXX
OPTIONS = -g -O3 -std=c++14
MACRO = $MACRO_ABACUS
INCLUDES = -I\${MKLROOT}/include
LIB_MKL = -L\${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core
LIB_PACS = -L${root_dir}/build/impurities/PACS -l:libpacs.a
LIB_iQIST = -L${root_dir}/build/impurities/iQIST/ct_hyb1 -l:libnarcissus.a \\
-L${root_dir}/build/impurities/iQIST/dependencies/Flink -l:libflink.a

LIB_RUTGERS = -L${root_dir}/build/impurities/Rutgers -l:librutgers.a \\
-L${root_dir}/build/impurities/Rutgers/dependencies/gsl-2.6/lib -l:libgslcblas.a -l:libgsl.a

LIBRARY = \$(LIB_PACS) \$(LIB_MKL) \$(LIB_iQIST) \$(LIB_RUTGERS)

VPATH=../src \\
:../src/para \\
:../src/solver

#============================
#     Objects
#============================
OBJS=main.o \\
parameters.o \\
mpi_environment.o \\
global_variables.o \\
utilities.o \\
correlated_atoms.o \\
KS_bands.o \\
tetrahedron.o \\
KS_eigenvectors.o \\
overlap_matrix.o \\
overlap_matrix_aims.o \\
overlap_matrix_abacus.o \\
input.o \\
timer.o \\
solver.o \\
projector.o \\
chemical_potential.o \\
Hilbert_space.o \\
self_energy.o \\
self_energy_imaginary_aixs.o \\
self_energy_real_aixs.o \\
double_counting.o \\
coulomb_tensor.o \\
Kanamori_parameterization.o \\
Anderson_impurity.o \\
PACS_cthyb.o \\
alps_cthyb.o \\
alps_cthyb_segment.o \\
rutgers_cthyb.o \\
iQIST_narcissus.o \\
math_zone.o \\
spectrum.o \\
charge_scf.o \\
charge_scf_aims.o

#====================
#   Target
#====================
all:\${OBJS}
	\${CPLUSPLUS_MPI} \${OPTIONS} \${OBJS} \${LIBRARY} -o ../bin/DFTDMFT

.PHONY:clean
clean:
	rm *.o
	rm -f ../bin/DFTDMFT
#	rm -f ../bin/maxent

#==========================
#       rules
#==========================
.cpp.o:
	\${CPLUSPLUS_MPI} \${OPTIONS} \${INCLUDES} -c \${MACRO} \$< -o \$@
EOF
  else #FHI-aims has been built
  cat > Makefile <<EOF
CPLUSPLUS_MPI = $MPI_CXX
OPTIONS = -g -O3 -std=c++14
MACRO = -D__FHIaims $MACRO_ABACUS
INCLUDES = -I\${MKLROOT}/include -I${FHIaims_lib_path}/include

LIB_MKL = -L\${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lmkl_blacs_intelmpi_lp64 \\
-lmkl_scalapack_lp64 -lifport -lifcoremt

LIB_PACS = -L${root_dir}/build/impurities/PACS -l:libpacs.a
LIB_iQIST = -L${root_dir}/build/impurities/iQIST/ct_hyb1 -l:libnarcissus.a \\
-L${root_dir}/build/impurities/iQIST/dependencies/Flink -l:libflink.a

LIB_RUTGERS = -L${root_dir}/build/impurities/Rutgers -l:librutgers.a \\
-L${root_dir}/build/impurities/Rutgers/dependencies/gsl-2.6/lib -l:libgslcblas.a -l:libgsl.a

LIB_aims = -L${FHIaims_lib_path}/lib -laims -lelsi -lelpa  -lOMM -lMatrixSwitch -lNTPoly -lfortjson

LIBRARY = \$(LIB_MKL) \$(LIB_PACS) \$(LIB_iQIST) \$(LIB_RUTGERS) \$(LIB_aims)

VPATH=../src \\
:../src/para \\
:../src/solver

#============================
#     Objects
#============================
OBJS=main.o \\
parameters.o \\
mpi_environment.o \\
global_variables.o \\
utilities.o \\
correlated_atoms.o \\
KS_bands.o \\
tetrahedron.o \\
KS_eigenvectors.o \\
overlap_matrix.o \\
overlap_matrix_aims.o \\
overlap_matrix_abacus.o \\
input.o \\
timer.o \\
solver.o \\
projector.o \\
chemical_potential.o \\
Hilbert_space.o \\
self_energy.o \\
self_energy_imaginary_aixs.o \\
self_energy_real_aixs.o \\
double_counting.o \\
coulomb_tensor.o \\
Kanamori_parameterization.o \\
Anderson_impurity.o \\
PACS_cthyb.o \\
alps_cthyb.o \\
alps_cthyb_segment.o \\
rutgers_cthyb.o \\
iQIST_narcissus.o \\
math_zone.o \\
spectrum.o \\
charge_scf.o \\
charge_scf_aims.o

#====================
#   Target
#====================
all:\${OBJS}
	\${CPLUSPLUS_MPI} \${OPTIONS} \${OBJS} \${LIBRARY} -o ../bin/DFTDMFT

.PHONY:clean
clean:
	rm *.o
	rm -f ../bin/DFTDMFT
#	rm -f ../bin/maxent

#==========================
#       rules
#==========================
.cpp.o:
	\${CPLUSPLUS_MPI} \${OPTIONS} \${INCLUDES} -c \${MACRO} \$< -o \$@
EOF
  fi
  make -j
fi

cd $root_dir/build

#====================================
#    PART 4: maxent
#====================================
test -d maxent || mkdir maxent
cd maxent 
if [ ! -f ../../bin/maxent ]; then

cat > Makefile << EOF
FC         = $MPI_FC
LIB        =  -L\${MKLROOT}/lib/intel64 -lmkl_rt -lpthread -lm -ldl

VPATH=../../src/maxent

.suffixes : .mod .o .f90

%.o:  %.f90
	\$(FC) -c \$<

SRC = maxentropy.f90

OBJ = \$(SRC:.f90=.o)

intel: \$(OBJ)
	\${FC} \${LIB} \$(OBJ) -o ../../bin/maxent
clean:
	rm -f *.o *.mod
EOF
  make

  if [ $? -ne 0 ]; then
    echo "Errors in building maxent"
    exit
  fi

fi
cd $root_dir/build


cd $root_dir/bin
#========Gw_AC.py=======
sed -i "/^maxent_exe=/c"maxent_exe=\"$root_dir/bin/maxent\""" Gw_AC.py
if [ ! -x Gw_AC.py ];then
  chmod +x Gw_AC.py
fi

#========Sigma_AC.py=======
sed -i "/^maxent_exe=/c"maxent_exe=\"$root_dir/bin/maxent\""" Sigma_AC.py
if [ ! -x Sigma_AC.py ];then
  chmod +x Sigma_AC.py
fi
cd $root_dir/build
