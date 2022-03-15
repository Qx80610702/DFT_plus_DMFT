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
#The installing path of FHI-aims if FHI-aims has been built  
FHIaims_install_dir=`grep "FHIaims_install_dir" install.vars | awk '{sub(/^[ \t]+/,"");print $3}'`
FHIaims_exe=`grep "FHIaims_exe" install.vars | awk '{sub(/^[ \t]+/,"");print $3}'`

#The installing path of ABACUS if ABACUS has been built  
ABACUS_install_dir=`grep "ABACUS_install_dir" install.vars | awk '{sub(/^[ \t]+/,"");print $3}'`
ABACUS_exe=`grep "ABACUS_exe" install.vars | awk '{sub(/^[ \t]+/,"");print $3}'`

#========Starting compilation========
cd ../
root_dir=$(dirname $(readlink -f "$0"))

#====================================
#    PART 1: libraries
#====================================

#====================================
#      gsl-2.6
#====================================
if [ ! -d $root_dir/build/libraries/gsl-2.6 ]
then
  cd $root_dir/libraries/
  test -d gsl-2.6 && rm -rf gsl-2.6
  tar -zxvf gsl-2.6.tar.gz
  cd gsl-2.6
  
  ./configure CC=mpiicc --prefix=$root_dir/build/libraries/gsl-2.6

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

  cd $root_dir
fi

#====================================
#    PART 2: impurity_solver
#====================================
test -d $root_dir/build/impurity_solver || mkdir $root_dir/build/impurity_solver

#========  PACS_CTHYB ========
if [ ! -f $root_dir/build/impurity_solver/PACS/pacs.cthyb ]
then
  test -d $root_dir/build/impurity_solver/PACS || mkdir $root_dir/build/impurity_solver/PACS
  cd $root_dir/src/impurity_solver
  test -d PACS_cthyb && rm -rf PACS_cthyb

  tar -zxvf PACS_cthyb.tar.gz
  cd PACS_cthyb

cat > Makefile << EOF
FC         = $MPI_FC
FC_FLAGS   = -O3 -xHost -nogen-interface # -ipo -traceback -implicitnone -warn all -check bounds   
#FC        = mpif90
#FC_FLAGS  = -Wall -Wextra -pedantic -Wconversion -fbacktrace -fbounds-check -ffree-line-length-none
LIB        =  -L\${MKLROOT}/lib/intel64 -lmkl_rt -lpthread -lm -ldl

.suffixes : .mod .o .f90

%.o:  %.f90
	\$(FC) \$(FC_FLAGS) -c \$<

SRC = GlobalVariables.f90 \\
      ctqmc_math.f90 MPI_mod.f90 \\
      Segment_Util.f90 Segment_MonteCarlo.f90 Segment_Phys.f90 \\
      Segment_LibraryMode.f90  Input_Parameters.f90 main.f90 

OBJ = \$(SRC:.f90=.o)

intel: \$(OBJ)
	\$(FC) \$(OBJ) \$(FC_FLAGS) -o ./pacs.cthyb \$(LIB)  
clean:
	rm -f *.o *.mod *~ ./pacs.cthyb
EOF

  make

  if [ $? -eq 0 ]; then
    mv pacs.cthyb $root_dir/build/impurity_solver/PACS/
  else
    echo "Errors in building PACS-CTHYB"
    exit
  fi
  cd $root_dir
else
  cd $root_dir/src/impurity_solver/PACS_cthyb
  # make clean
  # make
  # if [ $? -eq 0 ]; then
  #   cp CTHYB $root_dir/build/impurity_solver/CTHYB-LG/
  # else
  #   echo "Errors in building LG-CTHYB"
  #   exit
  # fi
  cd $root_dir
fi

#========  iQIST ========
if [ ! -f $root_dir/build/impurity_solver/iQIST/cthyb.narcissus ]
then
  test -d $root_dir/build/impurity_solver/iQIST || mkdir $root_dir/build/impurity_solver/iQIST
  cd $root_dir/src/impurity_solver
  test -d iQIST && rm -rf iQIST

  tar -zxvf iQIST.tar.gz
  cd iQIST

  #compile Flink
  cd dependencies/Flink/build
  sed -i "/^F90    =/cF90    = $MPI_FC" make.inc
  make
  if [ $? -ne 0 ]; then
    echo "Errors in compiling Flink"
    exit
  fi
  cd ../../../

  #compile iQIST
  cd build
  sed -i "/^F90    =/cF90    = $MPI_FC" make.inc
  sed -i "/^FLINK  =/cFLINK  = $root_dir/src/impurity_solver/iQIST/dependencies/Flink/src" make.inc
  make all
  if [ $? -eq 0 ]; then
    cd ../src/ct_hyb1
    cp cthyb.narcissus $root_dir/build/impurity_solver/iQIST/
  else
    echo "Errors in building iQIST"
    exit
  fi
  cd $root_dir
else
  cd $root_dir/src/impurity_solver/iQIST
  # make clean
  # make
  # if [ $? -eq 0 ]; then
  #   cp CTHYB $root_dir/build/impurity_solver/CTHYB-LG/
  # else
  #   echo "Errors in building LG-CTHYB"
  #   exit
  # fi
  cd $root_dir
fi

#=========  Rutgers_CTHYB  ========= 
if [ ! -f $root_dir/build/impurity_solver/Rutgers/ctqmc -o  ! -f $root_dir/build/impurity_solver/Rutgers/ctqmcf ];then
  test -d $root_dir/build/impurity_solver/Rutgers || mkdir $root_dir/build/impurity_solver/Rutgers
  cd $root_dir/src/impurity_solver
  test -d Rutgers && rm -rf Rutgers
  tar -zxvf Rutgers.tar.gz && cd Rutgers

cat > Makefile << EOF
PC++ = $MPI_CXX
F77 = $FC

GSLINC = -I$root_dir/build/libraries/gsl-2.6/include
GSLLIB =  -L$root_dir/build/libraries/gsl-2.6/lib -l:libgslcblas.a -l:libgsl.a

LLIBS  =  -mkl
PLIBS = \$(LLIBS) \$(GSLLIB)

PFLAGS   = -D_MPI -DMPICH_IGNORE_CXX_SEEK -O3 #-restrict -ipo -no-prec-div
CFLAGS = \$(PFLAGS) \$(GSLINC) -D_TIME #-DAS -D_TIME #-fast -xAVX # -D_TIME #-D_LOGGING #-DAS 
FFLAGS = -O2  -free -no-prec-div -pc80 -qopenmp 

GHEADERS = assert.h complex.h random.h sblas.h sfunction.h smesh.h sutil.h zeroin.h
QHEADERS = common.h  inout.h intervals.h local.h matrixm.h mpi.h bcast.h number.h operators.h state.h stateim.h segment.h svdfunc.h tanmesh.h

ctqmc : ctqmc.o SMatrix1.o 
	\$(PC++) \$(CFLAGS) -o \$@ ctqmc.o SMatrix1.o \$(PLIBS)

all : ctqmc ctqmcf

ctqmcf : ctqmcf.o SMatrix1.o
	\$(PC++) \$(CFLAGS) -o \$@ ctqmcf.o SMatrix1.o \$(PLIBS)

ctqmcf.o : ctqmc.cc
	\$(PC++) \$(CFLAGS) -DAS -c -o ctqmcf.o \$<

SMatrix1.o : SMatrix1.cc sfunction.h
	\$(PC++) -c \$(CFLAGS) SMatrix1.cc 

ctqmc.o : ctqmc.cc \$(GHEADERS) \$(QHEADERS) 
	\$(PC++) \$(CFLAGS) -c ctqmc.cc

clean :
	rm -f ctqmc.o ctqmcf.o ctqmcf ctqmc SMatrix1.o

.SUFFIXES : .cc
.cc.o:
	\$(PC++) \$(CFLAGS) -c \$<

.SUFFIXES : .f
.f.o:
	\$(F77) \$(FFLAGS) -c \$<
EOF

  make all
  if [ $? -eq 0 ]; then
    mv ctqmc $root_dir/build/impurity_solver/Rutgers/ctqmc
    mv ctqmcf $root_dir/build/impurity_solver/Rutgers/ctqmcf
  else
    echo "Errors in building Rutgers-CTHYB"
    exit
  fi
  cd $root_dir
fi

#====================================
#    PART 3: projecting_embedning
#====================================
if [ ! -f $root_dir/build/projection_embeding ]
then
  cd $root_dir/src/projecting_embending/build/
  rm ./*

  if [ -z $ABACUS_install_dir ];then
    MACRO_ABACUS=
  else
    MACRO_ABACUS=-D__ABACUS
  fi

  if [ -z $FHIaims_install_dir ];then    
cat > Makefile <<EOF
CPLUSPLUS_MPI = $MPI_CXX
OPTIONS = -g -qopenmp -O3 -std=c++14
MACRO = $MACRO_ABACUS
INCLUDES = -I\${MKLROOT}/include
LIBRARY = -L\${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core

VPATH=../src \\
:../src/para \\
:../src/solver

#============================
#     Objects
#============================
OBJS=main.o \\
parameters.o \\
mpi_environment.o \\
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
	\${CPLUSPLUS_MPI} \${OPTIONS} \${OBJS} \${LIBRARY} -o ../bin/projection_embeding

.PHONY:clean
clean:
	rm *.o ../bin/projection_embeding

#==========================
#       rules
#==========================
.cpp.o:
	\${CPLUSPLUS_MPI} \${OPTIONS} \${INCLUDES} -c \${MACRO} \$< -o \$@
EOF
  else   #FHI-aims has been built
  cat > Makefile <<EOF
CPLUSPLUS_MPI = $MPI_CXX
OPTIONS = -g -qopenmp -O3 -std=c++14
MACRO = -D__FHIaims $MACRO_ABACUS
INCLUDES = -I\${MKLROOT}/include -I${FHIaims_install_dir}/include
LIBRARY = -L$FHIaims_install_dir/lib -lelsi -lelpa  -lOMM -lMatrixSwitch -lNTPoly -lfortjson \\
-L\${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lmkl_blacs_intelmpi_lp64 \\
-lmkl_scalapack_lp64 -lifport -lifcoremt

VPATH=../src \\
:../src/para \\
:../src/solver

#============================
#     Objects
#============================
OBJS=main.o \\
parameters.o \\
mpi_environment.o \\
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
	\${CPLUSPLUS_MPI} \${OPTIONS} \${OBJS} \${LIBRARY} -o ../bin/projection_embeding

.PHONY:clean
clean:
	rm *.o ../bin/projection_embeding

#==========================
#       rules
#==========================
.cpp.o:
	\${CPLUSPLUS_MPI} \${OPTIONS} \${INCLUDES} -c \${MACRO} \$< -o \$@
EOF
  fi

  make -j

  if [ $? -eq 0 ]
  then
    mv ../bin/projection_embeding $root_dir/build/
  else
    echo "Errors in building projecting_embeding"
    exit
  fi
  cd $root_dir
else
  cd $root_dir/src/projecting_embending/build/
  # make clean
  make -j
  if [ $? -ne 0 ]
  then
    echo "Errors in building projecting_embeding"
    exit
  fi
  mv ../bin/projection_embeding $root_dir/build/
  cd $root_dir
fi

#====================================
#    PART 4: maxent
#====================================
if [ ! -f $root_dir/build/maxent/maxent ]; then
  test -d $root_dir/build/maxent || mkdir $root_dir/build/maxent
  cd $root_dir/src/maxent

  ifort -L\${MKLROOT}/lib/intel64 -lmkl_rt -lpthread -lm -ldl maxentropy.f90 -o maxent
  mv maxent $root_dir/build/maxent/

  cd $root_dir
fi

#====================================
#    PART 5: job script
#====================================
cd $root_dir/bin/
test -f run_dmft && rm run_dmft
cat > run_dmft <<EOF1
#!/bin/bash

#===========Execution path====================
EXE_DMFT=$root_dir/build/projection_embeding
EXE_ALPS_CTHYB=$root_dir/build/impurity_solver/ALPS-CTHYB/bin/hybmat
EXE_ALPS_CTHYB_SEGMENT=$root_dir/build/impurity_solver/ALPS-CTHYB-SEGMENT/bin/alps_cthyb
EXE_PACS_CTHYB=$root_dir/build/impurity_solver/PACS/pacs.cthyb
EXE_RUTGERS_CTHYB=$root_dir/build/impurity_solver/Rutgers/ctqmc
EXE_iQIST_CTHYB1=$root_dir/build/impurity_solver/iQIST/cthyb.narcissus
EXE_FHIaims=$FHIaims_exe
EXE_ABACUS=$ABACUS_exe

#===========================================
#            PART 1
#Determine the number of process and threads
#===========================================
nodes=1
num_threads=1

if [ \$# -eq 0 ];then
  echo "Warrning: Arguments nodes and num_threads are not given, and you will run the job with 1 process and 1 thread!!!"
fi

while [ \$# != 0 ]
do
  case \$1 in
  "-nodes")
    shift
    nodes=\$1
  ;;
  "-num_threads")
    shift
    num_threads=\$1
  ;;
  *)
    echo "ERROR: option \$1 does not exist"
  esac

  shift;
done

nprocess=\`echo "\$nodes * \$num_threads" | bc\`

#===================================================
#      PART 2: parsing file DMFT.in
#===================================================
#==============DFT solver type=================
dft_solver=\`grep -i "dft_solver" DMFT.in | awk '{sub(/^[ \\t]+/,"");print \$0}' | awk '{print \$1, \$2}' | grep -v "#" | awk '{print \$2}'\`
dft_solver_lower_case=\`echo \$dft_solver | tr A-Z a-z\`
case \$dft_solver_lower_case in
  "aims")
    dft_solver_type=1
  ;;
  "abacus")
    dft_solver_type=2
  ;;
  *)
    echo "ERROR: unsupported DFT solver \$dft_solver" 
    exit
esac
#echo "dft_solver_type: \$dft_solver_type"

#==============impurity solver type=================
impurity_solver=\`grep -i "impurity_solver" DMFT.in | awk '{sub(/^[ \\t]+/,"");print \$0}' | awk '{print \$1, \$2}' | grep -v "#" | awk '{print \$2}'\`
impurity_solver_lower_case=\`echo \$impurity_solver | tr A-Z a-z\`
case \$impurity_solver_lower_case in
  "alps-cthyb")
    impurity_solver_type=1
  ;;
  "alps-cthyb-segment")
    impurity_solver_type=2
  ;;
  "pacs")
    impurity_solver_type=3
  ;;
  "rutgers-cthyb")
    impurity_solver_type=4
  ;;
  "iqist")
    impurity_solver_type=5
  ;;
  "") #default value
    impurity_solver_type=3
  ;;
  *)
    echo "ERROR: unsupported impurity solver \$impurity_solver" 
    exit
esac
#echo "impurity_solver_type: \$impurity_solver_type"

#===============magnetism===================
magnetism=\`grep -i "magnetism" DMFT.in | awk '{sub(/^[ \\t]+/,"");print \$0}' | awk '{print \$1, \$2}' | grep -v "#" | awk '{print \$2}' | tr A-Z a-z\`

#=============Maximum charge step=================
max_charge_step=\`grep -i "max_charge_step" DMFT.in | awk '{sub(/^[ \\t]+/,"");print \$0}' | awk '{print \$1, \$2}' | grep -v "#" | awk '{print \$2}'\`
if [ -z \$max_charge_step ];then
  max_charge_step=1
fi

#=============Maximum DMFT step=================
max_DMFT_step=\`grep -i "max_DMFT_step" DMFT.in | awk '{sub(/^[ \\t]+/,"");print \$0}' | awk '{print \$1, \$2}' | grep -v "#" | awk '{print \$2}'\`
if [ -z \$max_DMFT_step ];then
  max_DMFT_step=5
fi

#========Current charge and DMFT step
start_charge_step=\`grep -i "restart" DMFT.in | awk '{sub(/^[ \\t]+/,"");print \$0}' | awk '{print \$1, \$2, \$3}' | grep -v "#" | awk '{print \$2}'\`
if [ -z \$start_charge_step ];then
  start_charge_step=1
fi

start_dmft_step=\`grep -i "restart" DMFT.in | awk '{sub(/^[ \\t]+/,"");print \$0}' | awk '{print \$1, \$2, \$3}' | grep -v "#" | awk '{print \$3}'\`
if [ -z \$start_dmft_step ];then
  start_dmft_step=1
fi

# if [ -d dmft_solving ]
# then
#   cd impurity_solving
#   current_step=1
#   for step_dir in \`ls\`
#   do
#     match=\`echo \$step_dir | grep "step"\`
#     if [ -n \$match ];then
#       step_tmp=\`echo \$step_dir | awk -F 'step' '{print \$2}'\`
#       if [ \$step_tmp -gt \$current_step ]
#       then
#         current_step=\$step_tmp
#       fi
#     fi
#   done
#   cd ..
# else
#   current_step=1
# fi

#================================================
#           PART 3: DFT+DMFT iteration
#================================================
for((char_step=\$start_charge_step;char_step<=\$max_charge_step;char_step=char_step+1))
do
  for((dmft_step=\$start_dmft_step;dmft_step<=\$max_DMFT_step;dmft_step=dmft_step+1))
  do
    #==================================================
    #     Run projecting and embeding
    #==================================================
    mpirun -n \$nodes -env OMP_NUM_THREADS=\$num_threads \$EXE_DMFT -charge.step \$char_step -dmft.step \$dmft_step -eva.density 0

    if [ \$? -ne 0 ]; then
      echo "Errors occured in running projecting_embeding"
      exit
    fi

    #=========judge whether convergency is reached==============
    flag_conver=0
    if [ \$dmft_step -gt 1 ]; then
      convergency=\`grep "DMFT self-consistency in DMFT loop" DMFT_running.log | awk 'END{print \$0}' | awk -F ": "   '{print \$2}'\`
      if [ "\$convergency" = "true" ]; then
        flag_conver=1
      else
        flag_conver=0
      fi
    fi
    if [ \$flag_conver -eq 1 ]; then
      echo "DMFT loop self-consistency has been reached. DMFT loop stopped" >> DMFT_running.log 
      break
    fi

    #==================================================
    #       impurity solving
    #==================================================
    echo "Impurity solver starts working......" >> DMFT_running.log
    echo "  Impurities       starting time              ending time     " >> DMFT_running.log
   #echo "impurity0    2021.03.30--15:44:28      2021.03.30--15:44:28"

    cd dmft/charge_step\$char_step/dmft_step\$dmft_step
    for dir_imp in \`ls\`
    do
      echo -e "  \$dir_imp    "\`date  "+%Y.%m.%d--%H:%M:%S"\`"\c" >> ../../../DMFT_running.log
      cd \$dir_imp

      if [ \$impurity_solver_type -eq 1 ]; then
        #================================================
        #             ALPS-CTHYB
        #================================================
        #seen_num=\`echo "\$nprocess + 10" | bc\`
        #sed -i "1a seed=\$seen_num" input.ini
        #==============================
        #    run impurity solver
        #==============================
        mpirun \$EXE_ALPS_CTHYB input.ini >> ./ALPS_CTHYB.log

        if [ \$? -ne 0 ]
        then
          echo "Errors occured in running impurity solver solving \$dir_imp !!!"
          exit
        fi

        python $root_dir/utilities/read_data.py >> ./ALPS_CTHYB.log

      elif [ \$impurity_solver_type -eq 2 ]; then
        #================================================
        #             ALPS-CTHYB-SEGMENT
        #================================================
        #========== parsing hyb.param ================
        FLAVORS=\`grep "FLAVORS=" hyb.param | awk -F '=' '{print \$2}'\`
        N_TAU=\`grep "N_TAU=" hyb.param | awk -F '=' '{print \$2}'\`
        BETA=\`grep "BETA=" hyb.param | awk -F '=' '{print \$2}'\`
        THERMALIZATION=\`grep "THERMALIZATION=" hyb.param | awk -F '=' '{print \$2}'\`
        N_MEAS=\`grep "N_MEAS=" hyb.param | awk -F '=' '{print \$2}'\`
        SWEEPS=\`grep "SWEEPS=" hyb.param | awk -F '=' '{print \$2}'\`
        DELTA=\`grep "DELTA=" hyb.param | awk -F '=' '{print \$2}'\`
        MU_VECTOR=\`grep "MU_VECTOR=" hyb.param | awk -F '=' '{print \$2}'\`
        U_MATRIX=\`grep "U_MATRIX=" hyb.param | awk -F '=' '{print \$2}'\`
        N_LEGENDRE=\`grep "N_LEGENDRE=" hyb.param | awk -F '=' '{print \$2}'\`
        N_HISTOGRAM_ORDERS=\`grep "N_HISTOGRAM_ORDERS=" hyb.param | awk -F '=' '{print \$2}'\`
        NMATSUBARA=\`grep "NMATSUBARA=" hyb.param | awk -F '=' '{print \$2}'\`

        #========== run alps_cthyb ================
        mpirun \$EXE_ALPS_CTHYB_SEGMENT \\
        -FLAVORS=\$FLAVORS -SEED=1000 -BETA=\$BETA \\
        -cthyb.DELTA=\$DELTA  -NMATSUBARA=\$NMATSUBARA -cthyb.MEASURE_sector_statistics=1 \\
        -MU_VECTOR=\$MU_VECTOR -U_MATRIX=\$U_MATRIX \\
        -cthyb.N_HISTOGRAM_ORDERS=\$N_HISTOGRAM_ORDERS -N=\$N_TAU \\
        -cthyb.N_MEAS=\$N_MEAS -cthyb.SWEEPS=\$SWEEPS -cthyb.THERMALIZATION=\$THERMALIZATION \\
        -cthyb.TEXT_OUTPUT=1 -MAX_TIME=2592000 \\
        -cthyb.N_LEGENDRE=\$N_LEGENDRE -VERBOSE=1 1>alps_cthyb.log 2>alps_cthyb.error

        if [ \$? -ne 0 ];then
          echo "Errors occured in running impurity solver solving \$dir_imp !!!"
          exit
        else
          if [ "\$magnetism"="para" -o "\$magnetism"="none" ]; then
            python $root_dir/utilities/alps_cthyb_segment_average.py
          fi
        fi
      elif [ \$impurity_solver_type -eq 3 ]; then
        mpirun \$EXE_PACS_CTHYB 1>PACS_cthyb.log 2>PACS_cthyb.error

        if [ \$? -ne 0 ];then
          echo "Errors occured in running impurity solver solving \$dir_imp !!!"
          exit
        fi

      elif [ \$impurity_solver_type -eq 4 ]; then
        mpirun \$EXE_RUTGERS_CTHYB 1>Rutgers_cthyb.log 2>Rutgers_cthyb.error

        if [ \$? -ne 0 ];then
          echo "Errors occured in running impurity solver solving \$dir_imp !!!"
          exit
        fi
      elif [ \$impurity_solver_type -eq 5 ]; then
        mpirun \$EXE_iQIST_CTHYB1 1>iQIST_cthyb.log 2>iQIST_cthyb.error

        if [ \$? -ne 0 ];then
          echo "Errors occured in running impurity solver solving \$dir_imp !!!"
          exit
        fi

      fi
      cd ..
      echo "      "\`date  "+%Y.%m.%d--%H:%M:%S"\` >> ../../../DMFT_running.log
    done #dir_imp

    cd ../../../
    #if [ \$i -eq 1 ]; then break; fi

  done #DMFT self-consistency loop

  start_dmft_step=1

  #============Charge update===============
  mpirun -n \$nodes -env OMP_NUM_THREADS=\$num_threads \$EXE_DMFT -charge.step \$char_step -dmft.step \$dmft_step -eva.density 1
  if [ \$? -ne 0 ];then
    echo "Errors occured in updating charge density!!!"
    exit
  fi

  cd dft
  rm -r outputs_to_DMFT

  if [ \$dft_solver_type -eq 1 ];then
    mpirun \$EXE_FHIaims 1>./job.log 2>./job.error
  elif [ \$dft_solver_type -eq 2 ];then
    mpirun \$EXE_ABACUS 1>./job.log 2>./job.error
  else
    echo "Errors: unspported DFT solver!!!"
    exit
  fi

  if [ \$? -ne 0 ];then
    echo "Errors occured in running dft solver!!!"
    exit
  fi
  cd ..
done #DFT+DMFT charge self-consistent loop 
EOF1

test -f cal_spectrum && rm cal_spectrum
cat > cal_spectrum << EOF2
#!/bin/bash
nodes=1
num_threads=1

if [ \$# -eq 0 ];then
  echo "Warrning: Arguments nodes and num_threads are not given, and you will run the job with 1 process and 1 thread!!!"
fi

while [ \$# != 0 ]
do
  case \$1 in
  "-nodes")
    shift
    nodes=\$1
  ;;
  "-num_threads")
    shift
    num_threads=\$1
  ;;
  *)
    echo "ERROR: option \$1 does not exist"
  esac

  shift;
done

mpirun -n \$nodes -env OMP_NUM_THREADS=\$num_threads $root_dir/build/projection_embeding -eva.spectrum 1
EOF2

chmod +x run_dmft
chmod +x cal_spectrum

#========Gw_AC.py=======
sed -i "/^maxent_exe=/c"maxent_exe=\"$root_dir/build/maxent/maxent\""" Gw_AC.py
if [ ! -x Gw_AC.py ];then
  chmod +x Gw_AC.py
fi

#========Sigma_AC.py=======
sed -i "/^maxent_exe=/c"maxent_exe=\"$root_dir/build/maxent/maxent\""" Sigma_AC.py
if [ ! -x Sigma_AC.py ];then
  chmod +x Sigma_AC.py
fi

#echo -e "\n" | cat >> ~/.bashrc
#echo "#################DFT+DMFT#################" | cat >> ~/.bashrc
#echo "export PATH=$root_dir/bin:\$PATH" | cat >> ~/.bashrc