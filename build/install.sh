#!/bin/bash

#============Compilers===============
CXX=icpc         #C++ compiler
CC=icc           #C complier
FC=ifort         #Fortran compiler
MPI_FC=mpiifort  #MPI version fortran compiler
MPI_CXX=mpiicpc  #MPI version C++ compiler
MPI_CC=mpiicc    #MPI version C compiler

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

#========  LG_CTHYB ========
if [ ! -f $root_dir/build/impurity_solver/CTHYB-LG/CTHYB ]
then
  test -d $root_dir/build/impurity_solver/CTHYB-LG || mkdir $root_dir/build/impurity_solver/CTHYB-LG
  cd $root_dir/src/impurity_solver
  test -d CTHYB-LG || mkdir CTHYB-LG
  cd CTHYB-LG

  cp ../CTHYB-LG.tar.gz .
  tar -zxvf CTHYB-LG.tar.gz
  rm CTHYB-LG.tar.gz

cat > Makefile << EOF
FC        = $MPI_FC
FC_FLAGS  = -O3 -xHost -nogen-interface # -ipo -traceback -implicitnone -warn all -check bounds  
#LIB   = -lmkl_intel_thread -lmkl_intel_lp64 -lmkl_core -liomp5 -lpthread
LIB   =  -L\${MKLROOT}/lib/intel64 -lmkl_rt -lpthread -lm -ldl

.suffixes : .mod .o .f90

%.o:  %.f90
	\$(FC) \$(FC_FLAGS) -c \$<

SRC = GlobalVariables.f90 ctqmc_math.f90 MPI_mod.f90 Segment_Util.f90 Segment_MonteCarlo.f90 Segment_Phys.f90 Segment_Main.f90  

OBJ = \$(SRC:.f90=.o)

intel: \$(OBJ)
	\$(FC) \$(OBJ) \$(FC_FLAGS) -o CTHYB \$(LIB)  
clean:
	rm -f *.o *.mod *~ CTHYB
EOF

  make
  if [ $? -eq 0 ]; then
    cp CTHYB $root_dir/build/impurity_solver/CTHYB-LG/
  else
    echo "Errors in building LG-CTHYB"
    exit
  fi
  cd $root_dir
else
  cd $root_dir/src/impurity_solver/CTHYB-LG
  make clean
  make
  if [ $? -eq 0 ]; then
    cp CTHYB $root_dir/build/impurity_solver/CTHYB-LG/
  else
    echo "Errors in building LG-CTHYB"
    exit
  fi
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
    cp ctqmc $root_dir/build/impurity_solver/Rutgers/ctqmc
    cp ctqmcf $root_dir/build/impurity_solver/Rutgers/ctqmcf
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

cat > Makefile <<EOF
CPLUSPLUS_MPI = $MPI_CXX
OPTIONS = -g -qopenmp -O3 -std=c++14
INCLUDES = -I\${MKLROOT}/include
LIBRARY = -L\${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core

VPATH=../src \\
:../src/para \\
:../src/debug \\
:../src/solver \\
:../src/post_processing \\

#============================
#     Objects
#============================
OBJ=main.o \\
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
LG_cthyb.o \\
alps_cthyb.o \\
alps_cthyb_segment.o \\
rutgers_cthyb.o \\
math_zone.o \\
post_processing.o \\
spectral_function.o 

#====================
#   Target
#====================
all:\${OBJ}
	\${CPLUSPLUS_MPI} \${OPTIONS} \${LIBRARY} \${OBJ} -o ../bin/projection_embeding

.PHONY:clean
clean:
	rm *.o $root_dir/build/projection_embeding

#==========================
#       rules
#==========================
.cpp.o:
	\${CPLUSPLUS_MPI} \${OPTIONS} \${INCLUDES} -c \$< -o \$@
EOF

  make -j

  if [ $? -eq 0 ]
  then
    cp ../bin/projection_embeding $root_dir/build/
  else
    echo "Errors in building projecting_embeding"
    exit
  fi
  cd $root_dir
else
  cd $root_dir/src/projecting_embending/build/
  make clean
  make -j
  if [ $? -ne 0 ]
  then
    echo "Errors in building projecting_embeding"
    exit
  fi
  cp ../bin/projection_embeding $root_dir/build/
  cd $root_dir
fi

#====================================
#    PART 4: maxent
#====================================
if [ ! -f $root_dir/build/analy_con/maxent ]; then
  test -d $root_dir/build/analy_continuation || mkdir $root_dir/build/analy_continuation
  cd $root_dir/src/analytic_continuation/maxent

  ifort -L\${MKLROOT}/lib/intel64 -lmkl_rt -lpthread -lm -ldl maxentropy.f90 -o maxent
  cp maxent $root_dir/build/analy_continuation/maxent

  cd $root_dir
fi

#====================================
#    PART 5: skrams
#====================================
if [ ! -f $root_dir/build/analy_continuation/skrams ];then
  test -d $root_dir/build/analy_continuation || mkdir $root_dir/build/analy_continuation
  cd $root_dir/src/analytic_continuation/skrams
  test -f skrams && rm skrams
  
  echo "Start building skrams ......"
  $CXX -O2 -funroll-all-loops -DNO_ARG_CHECK skrams.cc -o skrams
  echo "skrams building finish"

  cp skrams $root_dir/build/analy_continuation/skrams
  cd $root_dir
fi

#====================================
#    PART 6: job script
#====================================
cd $root_dir/bin/
test -f run_dmft && rm run_dmft
cat > run_dmft <<EOF1
#!/bin/bash

EXE_DMFT=$root_dir/build/projection_embeding
EXE_ALPS_CTHYB=$root_dir/build/impurity_solver/ALPS-CTHYB/bin/hybmat
EXE_ALPS_CTHYB_SEGMENT=$root_dir/build/impurity_solver/ALPS-CTHYB-SEGMENT/bin/alps_cthyb
EXE_LG_CTHYB=$root_dir/build/impurity_solver/CTHYB-LG/CTHYB
EXE_RUTGERS_CTHYB=$root_dir/build/impurity_solver/Rutgers/ctqmc

#===========================================
#            PART 1
#Determine the number of process and threads
#===========================================
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
    exit
  esac

  shift;
done

nprocess=\`echo "\$nodes * \$num_threads" | bc\`

#===================================================
#      PART 2: parsing file DMFT.in
#===================================================

#==============impurity solver type=================
impurity_solver=\`grep -i "impurity_solver" DMFT.in | grep -v "#" | awk '{print \$2}'\`
impurity_solver_lower_case=\`echo \$impurity_solver | tr A-Z a-z\`
case \$impurity_solver_lower_case in
  "alps_cthyb")
    impurity_solver_type=1
  ;;
  "alps_cthyb_segment")
    impurity_solver_type=2
  ;;
  "lg_cthyb")
    impurity_solver_type=3
  ;;
  "rutgers_cthyb")
    impurity_solver_type=4
  ;;
  *)
    echo "ERROR: unsupported impurity solver"
    exit
esac
#echo "impurity_solver_type: \$impurity_solver_type"

#===============magnetism===================
magnetism=\`grep -i "magnetism" DMFT.in | grep -v "#" | awk '{print \$2}' | tr A-Z a-z\`

#=============DMFT iteration step=================
DMFT_step=\`grep -i "DMFT_step" DMFT.in | grep -v "#" | awk '{print \$2}'\`

if [ -d impurity_solving ]
then
  cd impurity_solving
  current_step=1
  for step_dir in \`ls\`
  do
    step_tmp=\`echo \$step_dir | awk -F 'p' '{print \$2}'\`
    if [ \$step_tmp -gt \$current_step ]
    then
      current_step=\$step_tmp
    fi
  done
  cd ..
else
  current_step=1
fi

#================================================
#           PART 3: DMFT iteration
#================================================
for((i=\$current_step;i<=\$DMFT_step;i=i+1))
do
  #==================================================
  #     Run projecting and embeding
  #==================================================
  mpirun -n \$nodes -env OMP_NUM_THREADS=\$num_threads \$EXE_DMFT \\
  -dmft.step=\$i -eva.sigma_only=0

  if [ \$? -ne 0 ]; then
    echo "Errors occured in running projecting_embeding"
    exit
  fi

  #=========judge whether convergency is reached==============
  flag_conver=0
  if [ \$i -gt 1 ]; then
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
# echo "impurity0    2021.03.30--15:44:28      2021.03.30--15:44:28"

  cd impurity_solving/step\$i
  for dir_imp in \`ls\`
  do
    echo -e "  \$dir_imp    "\`date  "+%Y.%m.%d--%H:%M:%S"\`"\c" >> ../../DMFT_running.log
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
      mpirun \$EXE_LG_CTHYB 1>CTHYB.log 2>CTHYB.error

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
    fi
    cd ..
    echo "      "\`date  "+%Y.%m.%d--%H:%M:%S"\` >> ../../DMFT_running.log
  done #dir_imp

  cd ../../
  #if [ \$i -eq 1 ]; then break; fi

done #DMFT self-consistency loop
EOF1

chmod +x run_dmft

#========Gw_AC.py=======
sed -i "/^maxent_exe=/c"maxent_exe=\"$root_dir/build/analy_continuation/maxent\""" Gw_AC.py
if [ ! -x Gw_AC.py ];then
  chmod +x Gw_AC.py
fi

#========Sigma_AC.py=======
sed -i "/^skrams_exe=/c"skrams_exe=\"$root_dir/build/analy_continuation/skrams\""" Sigma_AC.py
sed -i "/^maxent_exe=/c"maxent_exe=\"$root_dir/build/analy_continuation/maxent\""" Sigma_AC.py
if [ ! -x Sigma_AC.py ];then
  chmod +x Sigma_AC.py
fi

#echo -e "\n" | cat >> ~/.bashrc
#echo "#################DFT+DMFT#################" | cat >> ~/.bashrc
#echo "export PATH=$root_dir/bin:\$PATH" | cat >> ~/.bashrc