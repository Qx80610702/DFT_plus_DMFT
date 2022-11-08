## List of all available parameters
  [dft_solver](#dft_solver) | [calculation](#calculation) | [temperature](#temperature) | [magnetism](#magnetism) | [projection_window](#projection_window) | [DOS_window](#dos_window) | [impurity_solver](#impurity_solver) | [max_charge_step](#max_charge_step) | [max_dmft_step](#max_dmft_step) | [max_dft_step](#max_dft_step) | [mc_step](#mc_step) | [charge_mix_param](#charge_mix_param) | [delta_sigma](#delta_sigma) | [delta_rho](#delta_rho) | [local_symmetry](#local_symmetry) | [restart](#restart)

[back to top](#list-of-all-available-parameters)

### dft_solver
- **Type**: String
- **Description**: Specify the DFT package used in DFT+DMFT calculation.
  - *aims*: the DFT calculations are carried out by FHI-aims
  - *abacus*: the DFT calculations are carried out by ABACUS
- **Default**: aims.

[back to top](#list-of-all-available-parameters)

### calculation
- **Type**: String
- **Description**: Specify the DFT+DMFT calculation type.
  - *scf*: self-consistency iteration
  - *spectra*: calculating the spectra
- **Default**: scf.

[back to top](#list-of-all-available-parameters)

### temperature
- **Type**: real
- **Description**: Specify the temperatue in the unit of Kelvin.
- **Default**: 300.

### magnetism
- **Type**: String
- **Description**: Specify the magnetism of the system.
  - *para*: paramagnetic
  - *none*: none-magnetic
  - *fm*: ferromagnetic
  - *afm*: antiferromagnetic
- **Default**: para.
- **!!!NOTE**: At present, only paramagnetism and none-magnetism are supported.

[back to top](#list-of-all-available-parameters)

### projection_window
- **Type**: Two reals
- **Description**: Specify the energy window for correlated subset of Kohn-Sham bands by two real values, where the first and second values are the lower and upper limit of the window respectively. The unit is electronvolt and the Fermi energy level is zero.
- **Default**: -5.0 5.0.

[back to top](#list-of-all-available-parameters)

### dos_window
- **Type**: Two reals
- **Description**: Specify the energy window for DOS claculation by two real values, where the first and second values are the lower and upper limit of the window respectively. The unit is electronvolt and the Fermi energy level is zero.
- **Default**: same as projection_window.

[back to top](#list-of-all-available-parameters)

### impurity_solver
- **Type**: String
- **Description**: Specify the impurity solver used for DMFT calculations.
  - *rutgers-cthyb*: A CT-HYB impurity solver. Refer to [Here](http://hauleweb.rutgers.edu/tutorials/Tutorial0.html)
  - *pacs*: A CT-HYB impurity solver
  - *iqist*: A CT-HYB impurity solver. Refer to [Here](https://github.com/huangli712/iQIST)
- **Default**: rutgers-cthyb.

[back to top](#list-of-all-available-parameters)

### max_charge_step
- **Type**: Integer
- **Description**: Specify the maximum charge loop. If you want to do one-shot DFT+DMFT, it should be set to be 1.
- **Default**: 1.

[back to top](#list-of-all-available-parameters)

### max_dmft_step
- **Type**: Integer
- **Description**: Specify the maximum DMFT loop under each chagre step. 
- **Default**: 1.

[back to top](#list-of-all-available-parameters)

### max_dft_step
- **Type**: Integer
- **Description**: Specify the maximum DFT loop under each chagre step. This parameter works only when max_charge_step is greater than 1.  
- **Default**: 1.

[back to top](#list-of-all-available-parameters)

### mc_step
- **Type**: Integer
- **Description**: The number of MC step for solving the impurity problem.  
- **Default**: 5000000.

[back to top](#list-of-all-available-parameters)

### charge_mix_param
- **Type**: Real
- **Description**: Specify the mixing parameter for charge desity mixing in Pulay-DIIS method.  
- **Default**: 0.05.

[back to top](#list-of-all-available-parameters)

### delta_sigma
- **Type**: Real
- **Description**: The convergrncy criterion for self-energy. Unit electronvolt.  
- **Default**: 0.1.

[back to top](#list-of-all-available-parameters)

### delta_rho
- **Type**: Real
- **Description**: The convergrncy criterion for charge density.
- **Default**: 1.0e-3.

[back to top](#list-of-all-available-parameters)

### local_symmetry
- **Type**: Integer
- **Description**: Specify the local symmetry of the strongly correlated electrons.
  - *0*: no symmtry 
  - *1*: cubic symmetry
  - *2*: only t2g orbitals
  - *3*: only eg orbitals 
- **Default**: 0.

[back to top](#list-of-all-available-parameters)

### restart
- **Type**: Four integers
- **Description**: Specify where to restart the DFT+DMFT calculations by four integer values, where the first, second, third and fourth values are the current charge step, current DMFT step, last charge step and last DMFT step respectively.
- **Default**: 1 1 1 1.

[back to top](#list-of-all-available-parameters)