### Introduction

SrVO$_3$ has a simple cubic perovskite structure without magnetism. The V$^{4+}$ cation with only one 3$d$ electron is located in the center of the octahedron formed by its six surrounding axial O$^{2−}$ ligand anions. In the presence of an octahedral crystal field, the five degenerate d orbitals split into two subsets: three folddegenerate $t_{2g}$ orbitals (i.e., $d_{xy}$, $d_{yz}$, and $d_{xz}$), and twofold degenerate $e_g$ orbitals (i.e., $d_{z^2}$ and $d_{x^2−y^2}$). The single 3$d$ electron of V$^{4+}$ occupies the lower-energy $t_{2g}$ orbitals, and the higher-energy $e_g$ orbitals are empty. As a common practice, we considered only the three degenerate $t_{2g}$ orbitals in our DFT+DMFT calculation.

### DFT+DMFT calculation 

Now we will explain how to do DFT+DMFT calculations on SrVO$_3$ step by step.

- `Step 1:` Do a self-soncistent DFT calculation.

  - `FHI-aims:` If you carry out DFT calculation with FHI-aims, you can use the following `geometry.in`
  ```bash 
  atom_frac    0.000000000         0.000000000         0.000000000    V
  atom_frac    0.500000000         0.500000000         0.500000000    Sr 
  atom_frac    0.500000000         0.000000000         0.000000000    O
  atom_frac    0.000000000         0.500000000         0.000000000    O
  atom_frac    0.000000000         0.000000000         0.500000000    O
  #
  #
  lattice_vector  3.8408999443         0.0000000000         0.0000000000
  lattice_vector  0.0000000000         3.8408999443         0.0000000000
  lattice_vector  0.0000000000         0.0000000000         3.8408999443

  ```
  
  For the file `control.in`, you can use the follwing file.
  ```bash 
  #--------------------------------
  #  Physical model settings
  #--------------------------------

      xc               pbe
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

      occupation_type  	gaussian 	0.05
      mixer            	pulay
      n_max_pulay             		10
      charge_mix_param    0.05
      # ini_linear_mixing     		10
      # ini_linear_mix_param  		0.05
       preconditioner kerker 		1.5
      # precondition_max_l    		0
      # preconditioner turnoff charge 	1e-4
      # preconditioner turnoff sum_ev 	1e-1

      sc_accuracy_rho  	1E-4
      sc_accuracy_eev  	1E-2
      sc_accuracy_etot 	1E-5

  #    sc_accuracy_forces 1E-4
      override_illconditioning .true.

  #----------------------------------
  #  For periodic boundary conditions
  #----------------------------------

      k_grid 11 11 11
      k_offset 0. 0. 0.

  #----------------------------------
  # DFT_plus_DMFT
  #----------------------------------
  DFT_plus_DMFT  .true.
    sc_iter_limit    80

  ########################################################################  ########
  #
  #  FHI-aims code project
  # Volker Blum, Fritz Haber Institute Berlin, 2009
  #
  #  Suggested "tight" defaults for V atom (to be pasted into control.in  file)
  #
  ########################################################################  ########
    species        V
  #     global species definitions
      nucleus             23
      mass                50.9415
  #
      l_hartree           6
  #
      cut_pot             4.0          2.0  1.0
      basis_dep_cutoff    1e-4
  #
      radial_base         49 7.0
      radial_multiplier   2
      angular_grids       specified
        division   0.2753   50
        division   0.6242  110
        division   0.9885  194
        division   1.1666  302
        division   1.3189  434
  #      division   1.5211  590
  #      division   1.6850  770
  #      division   1.8688  974
  #      division   3.0666 1202
  #      outer_grid  974
        outer_grid  434
  ########################################################################  ########
  #
  #  Definition of "minimal" basis
  #
  ########################################################################  ########
  #     valence basis states
      valence      4  s   2.
      valence      3  p   6.
      valence      3  d   3.
  #     ion occupancy
      ion_occ      4  s   1.
      ion_occ      3  p   6.
      ion_occ      3  d   2.
  ########################################################################  ########
  #
  #  Suggested additional basis functions. For production calculations, 
  #  uncomment them one after another (the most important basis functions   are
  #  listed first).
  #
  #  Constructed for dimers: 1.45 A, 1.65 A, 2.25 A, 3.00 A, 4.00 A
  #
  ########################################################################  ########
  #  "First tier" - improvements: -573.19 meV to -17.48 meV 
       hydro 4 f 9
       hydro 3 d 3
       ionic 4 p auto
       hydro 5 g 12.8
       ionic 4 s auto
  #  "Second tier" - improvements: -21.58 meV to -1.18 meV
  #     hydro 3 d 5.4
  #     hydro 5 f 11.2
  #     hydro 6 h 18.4
  #     hydro 4 d 7
  #     hydro 4 f 11.2
  #     hydro 4 p 5.6
  #     hydro 5 g 14
  #     hydro 1 s 0.6
  #  "Third tier" - improvements: -0.56 meV to -0.32 meV
  #     hydro 3 d 8.8
  #     hydro 4 p 7.8
  #     hydro 6 h 18.8
  #     hydro 4 f 24.8  
  #     hydro 4 s 4.0   
  #  "Fourth tier" - improvements: -0.30 meV to -0.09 meV
  #     hydro 5 p 12
  #     hydro 5 g 15.2
  #     hydro 5 f 8
  #     hydro 5 p 6.4
  #     hydro 4 d 5.2
  #     hydro 5 s 7.8
  #  Further functions - impr. -0.09 meV and below
  #     hydro 3 s 12
  #     hydro 6 h 20
  #     hydro 5 g 7
    DMFT 3 d 4.00 0.65

  ########################################################################  ########
  #
  #  FHI-aims code project
  # Volker Blum, Fritz Haber Institute Berlin, 2009
  #
  #  Suggested "tight" defaults for Sr atom (to be pasted into control.in   file)
  #
  #  2016/03/28 Added "g" function from tier 2 to default basis set   (noticeable
  #             improvement in Delta Project)
  #
  #  2018/01/30 Increased the cut_pot radius for "tight" to 6.0 AA to be  consistent
  #             with Rb, Cs, Ba. This is a VERY large value and will make   especially
  #             hybrid functional calculations very costly. However, our  data for
  #             neutral Sr containing systems support this value (it's  just a big atom).
  #             Consider reducing the cut_pot value for ionic bonding   situations,
  #             which could be entirely justified (see intermediate   settings).
  #
  ########################################################################  ########
    species          Sr
  #     global species definitions
      nucleus        38
      mass           87.62
  #
      l_hartree      6
  #
      cut_pot        6.0  2.0  1.0
      basis_dep_cutoff    1e-4
  #
      radial_base    57  7.0
      radial_multiplier  2
      angular_grids specified
        division   0.6981  110
        division   0.9394  194
        division   1.1230  302
        division   1.2482  434
  #      division   1.3391  590
  #      division   1.4365  770
  #      division   7.0005  974
  #      outer_grid  974
        outer_grid  434
  ########################################################################  ########
  #
  #  Definition of "minimal" basis
  #
  ########################################################################  ########
  #     valence basis states
      valence      5  s   2.
      valence      4  p   6.
      valence      3  d  10.
  #     ion occupancy
      ion_occ      5  s   1.
      ion_occ      4  p   6.
      ion_occ      3  d  10.
  ########################################################################  ########
  #
  #  Suggested additional basis functions. For production calculations, 
  #  uncomment them one after another (the most important basis functions   are
  #  listed first).
  #
  #  Constructed for dimers: 2.75, 3.50, 4.40, 5.00 A
  #
  ########################################################################  ########
  #  "First tier" - improvements: -289.57 meV to -14.02 meV
       ionic 4 d auto
       ionic 5 p auto
       hydro 4 f 5.6
       ionic 5 s auto
  #  "Second tier" - improvements: -4.95 meV to -0.45 meV
       hydro 5 g 7.4
  #     hydro 4 d 4.4
  #     hydro 3 p 3.3
  #     hydro 6 h 10.4
  #     hydro 5 s 4.9
  #     hydro 5 f 13.2
  #  "Third tier" - improvements: -0.38 meV to -0.11 meV
  #     hydro 6 p 4.8
  #     hydro 5 f 6
  #     hydro 2 p 1.2
  #     hydro 1 s 0.55
  #     hydro 5 d 3.6   
  #  "Fourth tier" - improvements: -0.12 meV and lower.
  #     hydro 5 p 5.2
  #     hydro 4 f 14.8
  #     hydro 5 g 7.6
  #     hydro 4 p 4.5
  #     hydro 5 d 5.4
  #     hydro 6 s 6.8   

  ########################################################################  ########
  #
  #  FHI-aims code project
  # Volker Blum, Fritz Haber Institute Berlin, 2009
  #
  #  Suggested "tight" defaults for O atom (to be pasted into control.in  file)
  #
  ########################################################################  ########
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
  ########################################################################  ########
  #
  #  Definition of "minimal" basis
  #
  ########################################################################  ########
  #     valence basis states
      valence      2  s   2.
      valence      2  p   4.
  #     ion occupancy
      ion_occ      2  s   1.
      ion_occ      2  p   3.
  ########################################################################  ########
  #
  #  Suggested additional basis functions. For production calculations, 
  #  uncomment them one after another (the most important basis functions   are
  #  listed first).
  #
  #  Constructed for dimers: 1.0 A, 1.208 A, 1.5 A, 2.0 A, 3.0 A
  #
  ########################################################################  ########
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
  ########################################################################  ########
  #
  # For methods that use the localized form of the "resolution of   identity" for
  # the two-electron Coulomb operator (RI_method LVL), particularly   Hartree-Fock and
  # hybrid density functional calculations, the highest accuracy can be   obtained by
  # uncommenting the line beginning with "for_aux"  below, thus adding an   extra g radial
  # function to the construction of the product basis set for the   expansion.
  # See Ref. New J. Phys. 17, 093020 (2015) for more information,   particularly Figs. 1 and 6.
  #
  ########################################################################  
  ########
  #
  # for_aux hydro 5 g 6.0
  ```

  The parameter `DFT_plus_DMFT` in `control.in` should be switched on. The line `DMFT 3 d 4.00 0.65` is composed by the paramter `DMFT` followed by four values, i.e., the main quantum number of the correlated orbitals, the angular momentum of the correlated orbitals, the U and J values. This line should be added to the end of the species, which you want to employ DMFT correction. 

  To get a good DFT band structure, a dense $k$-mesh, e.g., $11\times 11\times 11$, is recomended. After the DFT iterartion converged, a folder named `outputs_to_DMFT` will be created.

  - `ABACUS:` If you carry out DFT calculation with ABACUS. The follwing lines
  ```bash 
  dft_plus_dmft    1
  orbital_corr    2 -1 -1
  hubbard_u    4.00  0.0 0.0
  hund_j       0.65  0.0 0.0
  ```
  should be added to the `INPUT` file. After the DFT iterartion converged, a folder named `outputs_to_DMFT` will be created. Unfortunately, ABACUS can not output the structure symmetry, a file named `symmetry.dat` must be added to the directory `outputs_to_DMFT`. For the conventional cell of SrVO$_3$, i.e., the `geometry.in` in FHI-aims, the `symmetry.dat` likes as 
  ```bash 
  1
  1
  0    0
  ```
  The first line specifies the total numuber of correlated atoms in the cell. The second line specifies the total numuber of non-equivalent correlated atoms in the cell. Starting from the third line, each line contians two integers, which are the numbering of the correlated atoms in the cell and the numbering of their equivalent correlated atoms.

- `Step 2:` Create a working folder, which will works as the root directory of the DFT+DMFT. Then copy the full folder of the above DFT scf calculation to the `dft` folder.

- `Step 3:` Create a file `DMFT.in`
  ```bash 
  DFT_solver   AIMS
  calculation scf
  temperature   300
  impurity_solver rutgers-cthyb
  magnetism     para
  max_charge_step 1
  max_DMFT_step    15
  max_DFT_step    1
  charge_mix_param  0.05
  delta_sigma 0.01
  delta_rho  1.0e-4
  MC_step      50000000
  projection_window -2.0 2.0
  dos_window  -20.0 20.0
  local_symmetry 2
  ```
  The meaning of the parameters please refer to [here](../../list_of_parameters.md).

- `Step 4:` Start DFT+DMFT calculation through follwong command
  ```bash
  mpirun -n 48 DFTDMFT
  ```
  This step is very computationally cost. You can monitor the job through the log fie `DMFT_running.log`.

  **!!!Note:** The number of cores must be equal to the job of the DFT scf calculation if you want to do charge self-consistent DFT+DMFT calculation.

- `Step 5:` Calculate the spectral function

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
  mpirun -n 6 Gw_AC.py
  ```
  The Gw_AC.py are orbital-paralleled. For SrVO$_3$, there are 6 $t_{2g}$ orbitals (including spin freedom), so you can run this script with 6 cores. After the analytical continuation finished, a folder named `Aw_loc` will be created. In this folder, the `spectrum.dat` contians the spectrum of the $3d$ electrons of SrVO$_3$. Plot it and you will get the follwing spectrum

<p align="center">
   <img src="./SrVO3-Aw.png">
</p>


[Back to hand-on examples](../hands-on-examples.md)