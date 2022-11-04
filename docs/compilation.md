
## Prerequisites

To install this package, the minimal requirements are listed as following:
- `C++ compiler:` Intel >= 18.0
- `Fortran compiler`
- `System environment:` gcc>=6.3.0, python3, scipy, numpy, mpi4py
- `Libraries:` Intel MKL

## Complilation

### Building DFT codes
If you want to do DFT+DMFT calculations based on the ABACUS, you need to build the ABACUS at first (Refer to [Building ABACUS](#building-abacus)). If you want to do DFT+DMFT calculations based on the FHI-aims, you need to build the FHI-aims at first (Refer to [Building FHI-aims](#building-fhi-aims)). We suggest you to use the same compliers and libraries when building the DFT codes and this DFT+DMFT package.

#### Building ABACUS

This is optional. It's required if you want to do DFT+DMFT calcualtions with ABACUS.

- `Download:`ABACUS is an open-source package, and you can download the latest version (a safe suggestion) from the [github](https://github.com/abacusmodeling/abacus-develop).

- `Compilation:` Compile the ABACUS following the introduction in its manual.


#### Building FHI-aims

This is optional. It's required if you want to do DFT+DMFT calcualtions with FHI-aims.

- `Download:`FHI-aims is a full-potential all-electrons LCNAO DFT package, and you can get the latest version (a safe suggestion) from its [offical website](https://fhi-aims.org/).

- `Build the FHI-aims executable:`Compile the FHI-aims as introduced in its manual. 

- `Build the FHI-aims shared library:`
    i) copy the file libaims.intel.cmake to the builindg directory from the aims_root_dir/cmake/toolchains/
    ii) Type the following command in the teminate
    ```bash
    cmake -C libaims.intel.cmake aims_root_dir -DCMAKE_INSTALL_PREFIX=path_to_installing_aims_shared_lib
    ```
    iii) Type the following command in the teminate
    ```bash
    make -j
    ```
    iv) Type the following command in the teminate
    ```bash 
    make install
    ```
    v) Type the following command in the teminate
    ```bash
    cp libaims.so path_to_installing_aims_shared_lib/lib/
    ```
    vi) Add the follwing path 
    ```bash 
    export LD_LIBRARY_PATH=path_to_installing_aims_shared_lib/lib:$LD_LIBRARY_PATH
    ```
    to the file .bashrc

### Building DFT+DMFT codes
i) Go to the directory `build/` under the root directory of the DFT+DMFT code

ii) Adjust the follwing variables according to the environment of your system in the file install.vars
 
  ```bash 
  #============Compilers===============
  CXX     =   icpc          #C++ compiler
  CC      =   icc           #C complier
  FC      =   ifort         #Fortran compiler
  MPI_FC  =   mpiifort      #MPI version fortran compiler
  MPI_CXX =   mpiicpc       #MPI version C++ compiler
  MPI_CC  =   mpiicc        #MPI version C compiler

  #============DFT library path=========
  FHIaims_lib_path   =  path_to_installing_aims_shared_lib  #Required if DFT calculations were carried out by FHI-aims
  ABACUS_lib_path    =  path_of_ABACUS                      #Required if DFT calculations were carried out by ABACUS
  ```
  
iii) Type the following command in the teminate to start compilation
  ```bash 
  bash install.sh
   ```
The compilation will take several minutes. If the compilation runs successfully, there will be two binary executables, i.e., DFTDMFT and maxent, and three executable python scripts, i.e., Gw_AC.py, kpath_generate.py and Sigma_AC.py, in the directory `bin/` under the root directory of the DFT+DMFT code

iv) Add the following paht 
  ```bash 
  export PATH=root_DFT_plus_DMFT/bin:$PATH
  ```
to the file .bashrc