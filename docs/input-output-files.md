## Input files
Two kinds of input files, i.e., DFT and DMFT files, are required if one wants to do DFT+DMFT calculations. The DFT files include all input and output files of DFT, and they must be grouped together into the folder named `dft/` under the root directory where the DFT+DMFT calculations are carried out. The DMFT files include all parameters controlling the DFT+DMFT calculations.

### How to prepare DFT files
  -i) Create a folder named dft under the root working directory. 
  
  -ii) Do a scf DFT calculation in the `dft/` directory as introduced in the manual of FHI-aims and ABACUS. In the DFT calculations, several extra parameters are required in the DFT input files, e.g., the INPUT file in ABACUS and control.in file in FHI-aims. Thoese parameters are listed in the manual of FHI-aims and ABACUS with explanations. There will be a folder named outputs_to_DMFT, which contains all data need to be passed to DMFT by DFT, after the DFT calculation finished.

### How to prepare DMFT files
  - `DMFT.in:` This file controls how the DMFT calulation will be carried out. An example is presented as follows
    ```bash
    DFT_solver   AIMS
    calculation scf
    temperature   1160
    impurity_solver rutgers-cthyb
    magnetism     para
    max_charge_step 10
    max_DMFT_step     1
    max_DFT_step    10
    charge_mix_param  0.05
    delta_sigma 0.01
    delta_rho  1.0e-4
    MC_step      50000000
    energy_window -2.0 1.0
    local_symmetry 1
    ```
    All available parameters are explictly explained in [here](list_of_parameters.md). Note that all parameters and their values are case-insensitive.


## Output files
After the DFT+DMFT calculations finished, there will be following output files

  -`DMFT_running.log`: This is the log file of the DFT+DMFT calculations. It contains information about the DFT+DMFT calculations.

  -`dmft/`:This folder contains all outputs given by impurity solvers.
