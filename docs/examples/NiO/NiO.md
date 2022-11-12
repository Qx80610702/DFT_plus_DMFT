### Introduction

NiO is a classic example that clearly demonstrate the failure of traditional band theory, e.g., LDA and GGA. LDA and GGA predict NiO to be metallic, whereas it's a wide-gap insulator experimentally. Even though a gap appears when the spin-symmetry-broken antiferromagnetic state is considered, the calculated gap is still much smaller than the experimental values. By merging the Hubbard model and ab initio DFT,  DFT+DMFT successfully opens the gaps in this prototypical transitional metal oxides.

### DFT+DMFT calculation 

Now we will explain how to do DFT+DMFT calculations on NiO step by step.

- `Step 1:` Do a self-soncistent DFT calculation.

  - `FHI-aims:` If you carry out DFT calculation with FHI-aims, you can use the following `geometry.in`
  ```bash
  atom_frac    0.000000000         0.000000000         0.000000000    Ni
  atom_frac    0.500000000         0.500000000         0.500000000    O
  #
  #
  lattice_vector  2.0890000000         2.0890000000         0.0000000000
  lattice_vector  0.0000000000         2.0890000000         2.0890000000
  lattice_vector  2.0890000000         0.0000000000         2.0890000000
  ```

  For the file `control.in`, you can use the follwing file.
  ```bash
  #--------------------------------
  #  Physical model settings
  #--------------------------------

      xc               pw-lda
      charge           0.
      spin             none
  #   relativistic     none
  #   override_relativity .true.
      relativistic     atomic_zora scalar 1E-12

  #   fixed_spin_moment 1 
  #   default_initial_moment 0.01

  #--------------------------------
  #  SCF convergence settings
  #--------------------------------

      occupation_type  	gaussian 	0.01
      mixer            	pulay
      n_max_pulay             		10
      charge_mix_param        		0.4
      # ini_linear_mixing     		10
      # ini_linear_mix_param  		0.05
       preconditioner kerker 		1.5
      # precondition_max_l    		0
      # preconditioner turnoff charge 	1e-4
      # preconditioner turnoff sum_ev 	1e-1

      sc_accuracy_rho  	1E-4
      sc_accuracy_eev  	1E-2
      sc_accuracy_etot 	1E-5
      sc_iter_limit    	80

  #    sc_accuracy_forces 1E-4


  #----------------------------------
  #  For periodic boundary conditions
  #----------------------------------

      k_grid 11 11 11
      k_offset 0. 0. 0.

  #----------------------------------
  # DFT_plus_DMFT
  #----------------------------------
    DFT_plus_DMFT  .true.

  ####################################################################### #########
  #
  #  FHI-aims code project
  # Volker Blum, Fritz Haber Institute Berlin, 2010
  #
  #  Suggested "tight" defaults for Ni atom (to be pasted into control. in file)
  #
  #  Revised June 14, 2011.
  #     p and d functions of tier 2 now enabled by default, as the atom   may change
  #     its occupation. Similar to Co.
  #
  ####################################################################### #########
    species        Ni
  #     global species definitions
      nucleus             28
      mass                58.6934
  #
      l_hartree           6
  #
      cut_pot             4.0          2.0  1.0
      basis_dep_cutoff    1e-4
  #
      radial_base         52 7.0
      radial_multiplier   2
      angular_grids       specified
        division   0.2935   50
        division   0.6132  110
        division   0.9287  194
        division   1.1299  302
        division   1.3700  434
  #      division   1.5675  590
  #      division   1.7612  770
  #      division   1.9438  974
  #      division   2.5441 1202
        outer_grid  434
  ####################################################################### #########
  #
  #  Definition of "minimal" basis
  #
  ####################################################################### #########
  #     valence basis states
      valence      4  s   2.
      valence      3  p   6.
      valence      3  d   8.
  #     ion occupancy
      ion_occ      4  s   1.
      ion_occ      3  p   6.
      ion_occ      3  d   7.
  ####################################################################### #########
  #
  #  Suggested additional basis functions. For production calculations, 
  #  uncomment them one after another (the most important basis   functions are
  #  listed first).
  #
  #  Constructed for dimers: 1.65 A, 2.00 A, 2.50 A, 3.00 A, 4.00 A
  #
  ####################################################################### #########
  #  "First tier" - improvements: -123.08 meV to -11.61 meV 
       hydro 3 p 6
       hydro 4 f 9
       hydro 5 g 12.4
       hydro 3 d 5.2
       ionic 4 s auto
  #  "Second tier" - improvements: -6.71 meV to -1.07 meV
       ionic 4 p auto
       hydro 4 d 6
  #     hydro 6 h 18
  #     hydro 4 f 9.4
  #     hydro 4 f 16.4
  #     hydro 1 s 0.75
  #  "Third tier" - improvements: -0.57 meV to -0.07 meV
  #     hydro 4 p 18.4
  #     hydro 4 d 8
  #     hydro 5 g 13.2
  #     hydro 5 f 8.4
  #     hydro 6 h 16.8
  #     hydro 4 s 4.4
  #  Further functions: improvements -0.07 meV and below
  #     hydro 5 f 16.8
  #     hydro 4 p 10

        DMFT 3 d 8.0 1.0

  ####################################################################### #########
  #
  #  FHI-aims code project
  # Volker Blum, Fritz Haber Institute Berlin, 2009
  #
  #  Suggested "tight" defaults for O atom (to be pasted into control.in  file)
  #
  ####################################################################### #########
    species        O
  #     global species definitions
      nucleus             8
      mass                15.9994
  #
      l_hartree           6
  #
      cut_pot             4.0  2.0  1.0
      basis_dep_cutoff    1e-4
  #
      radial_base         36 7.0
      radial_multiplier   2
       angular_grids specified
        division   0.1817   50
        division   0.3417  110
        division   0.4949  194
        division   0.6251  302
        division   0.8014  434
  #      division   0.8507  590
  #      division   0.8762  770
  #      division   0.9023  974
  #      division   1.2339 1202
  #      outer_grid 974
        outer_grid  434
  ####################################################################### #########
  #
  #  Definition of "minimal" basis
  #
  ####################################################################### #########
  #     valence basis states
      valence      2  s   2.
      valence      2  p   4.
  #     ion occupancy
      ion_occ      2  s   1.
      ion_occ      2  p   3.
  ####################################################################### #########
  #
  #  Suggested additional basis functions. For production calculations, 
  #  uncomment them one after another (the most important basis   functions are
  #  listed first).
  #
  #  Constructed for dimers: 1.0 A, 1.208 A, 1.5 A, 2.0 A, 3.0 A
  #
  ####################################################################### #########
  #  "First tier" - improvements: -699.05 meV to -159.38 meV
       hydro 2 p 1.8
       hydro 3 d 7.6
       hydro 3 s 6.4
  #  "Second tier" - improvements: -49.91 meV to -5.39 meV
       hydro 4 f 11.6
       hydro 3 p 6.2
       hydro 3 d 5.6
       hydro 5 g 17.6
       hydro 1 s 0.75
  #  "Third tier" - improvements: -2.83 meV to -0.50 meV
  #     ionic 2 p auto
  #     hydro 4 f 10.8
  #     hydro 4 d 4.7
  #     hydro 2 s 6.8
  #  "Fourth tier" - improvements: -0.40 meV to -0.12 meV
  #     hydro 3 p 5
  #     hydro 3 s 3.3
  #     hydro 5 g 15.6
  #     hydro 4 f 17.6
  #     hydro 4 d 14
  # Further basis functions - -0.08 meV and below
  #     hydro 3 s 2.1
  #     hydro 4 d 11.6
  #     hydro 3 p 16
  #     hydro 2 s 17.2
  ####################################################################### #########
  #
  # For methods that use the localized form of the "resolution of   identity" for
  # the two-electron Coulomb operator (RI_method LVL), particularly   Hartree-Fock and
  # hybrid density functional calculations, the highest accuracy can be   obtained by
  # uncommenting the line beginning with "for_aux"  below, thus adding  an extra g radial
  # function to the construction of the product basis set for the   expansion.
  # See Ref. New J. Phys. 17, 093020 (2015) for more information,   particularly Figs. 1 and 6.
  #
  ####################################################################### #########
  #
  # for_aux hydro 5 g 6.0
  ```

  As clearly explained in the case of [SrVO$_3$](../SrVO3/SrVO3.md), the U and J parameters of Ni 3$d$ electrons are 8.0 and 1.0eV repectively.

  - `ABACUS:` If you carry out DFT calculation with ABACUS. The follwing lines
  ```bash 
  dft_plus_dmft    1
  orbital_corr    2 -1
  hubbard_u    8.0  0.0
  hund_j       1.0  0.0
  ```
  should be added to the `INPUT` file. After the DFT iterartion converged, a folder named `outputs_to_DMFT` will be created. Unfortunately, ABACUS can not output the structure symmetry, a file named `symmetry.dat` must be added to the directory `outputs_to_DMFT`. For the primitive cell of NiO, i.e., the `geometry.in` in FHI-aims, the `symmetry.dat` likes as 
  ```bash 
  1
  1
  0    0
  ```
  The format of this file could be refered to case of [SrVO$_3$](../SrVO3/SrVO3.md).

  - `Step 2:` Create a working folder, which will works as the root directory of the DFT+DMFT. Then copy the full folder of the above DFT scf calculation to the `dft` folder.

  - `Step 3:` Create a file `DMFT.in`
  ```bash 
  DFT_solver   AIMS
  calculation  spectra
  temperature   1160
  impurity_solver rutgers-cthyb
  magnetism     para
  max_charge_step  1
  max_DMFT_step     10
  max_DFT_step 1
  delta_sigma  0.01
  delta_rho  1.0e-4
  MC_step      100000000
  projection_window -2.0 1.0
  dos_window   -20.0 20.0
  local_symmetry 1
  ```
  The meaning of all the parameters please refer to [here](../../list_of_parameters.md).

- `Step 4:` Start DFT+DMFT calculation through follwong command
  ```bash
  mpirun -n 48 DFTDMFT
  ```
  This step is very computationally cost. You can monitor the job through the log fie `DMFT_running.log`.

  **!!!Note:** The number of cores must be equal to the job of the DFT scf calculation if you want to do charge self-consistent DFT+DMFT calculation.

- `Step 5:` Calculate the 3$d$ eletrons spectral function

  - `maxent_params.dat:` Prepare the follwing file `maxent_params.dat` for analytical continuation.
  ```bash 
  params={
        'Ntau'      : 400,     # Number of time points
        'bwdth'     : 0.1,     # smoothing width
        'Nw'        : 200,     # number of frequency points on real axis
        'Dw'        : 0.1,     # energy increment
        'idg'       : 1,       # error scheme: idg=1 -> sigma=deltag ; idg=0 -> sigma=deltag*G(tau)
        'deltag'    : 0.001,   # error
        'Asteps'    : 8000,    # anealing steps
        'alpha0'    : 2000,    # starting alpha
        'iflat'     : 1,       # iflat=0 : constant model, iflat=1 : gaussian of width gwidth, iflat=2 : input using file model$
        'Nr'        : 0,       # number of smoothing runs
        'eta'       : 0.01,     # Broadening in Kramer-Kronig tranform
    }
  ```

  - `Analytical continuation:` Using the MEM to do analytical continuation of the impurity Green's function through the follwing command line
  ```bash 
  mpirun -n 10 Gw_AC.py
  ```
  The Gw_AC.py are orbital-paralleled. For NiO, there are 10 3$d$ orbitals (including spin freedom), so you can run this script with 10 cores. After the analytical continuation finished, a folder named `Aw_loc` will be created. In this folder, the `spectrum.dat` contians the spectrum of the $3d$ electrons of NiO. Plot it and you will get the follwing spectrum

<p align="center">
   <img src="./NiO-Aw.png">
</p>

- `Step 6:` Calculate the k-resolved spectral function


