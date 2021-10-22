#pragma once

#include "../para/correlated_atoms.h"
#include "../para/input.h"
#include "../para/KS_bands.h"

#include <complex>
#include <vector>

namespace DMFT
{
  class self_energy_real_aixs
  {
    public:
    self_energy_real_aixs(){;}
    ~self_energy_real_aixs(){;}
    
    public:
    void read_AC_sigma(DFT_output::KS_bands& band,
              DMFT::input_info& in, 
              DFT_output::atoms_info& atom);

    void read_maxent_parms();

    public: //interface
    int& nomega(){return n_omega;} 
    std::vector<double>& frequency(){return real_frequency;}

    inline std::vector<std::vector<std::vector<std::vector<std::complex<double>>>>>& 
            sigma_new_access(){return local_sigma_new;}

    inline std::vector<std::vector<std::vector<std::vector<std::complex<double>>>>>& 
            sigma_save_access(){return local_sigma_save;}

    inline std::vector<std::vector<std::vector<std::vector<std::complex<double>>>>>& 
            correlated_sigma_access(){return correlated_sigma_omega_n;}

    inline std::vector<std::vector<std::vector<std::complex<double>>>>&
            lattice_sigma_access(){return lattice_sigma;}

    private:
    int n_omega;   //Num of real frequency points

    //The real frquence
    std::vector<double> real_frequency;

    //local_sigma_new[ineq][ispin][omega_n][m_index]
    std::vector<std::vector<std::vector<std::vector<std::complex<double>>>>> local_sigma_new;
    
    //the local self energy of last DMFT step which is read from file
    //local_sigma_save[ineq][ispin][omega_n][m_index]
    std::vector<std::vector<std::vector<std::vector<std::complex<double>>>>> local_sigma_save;

    //the self energy where double counting term is subtracted
    //correlated_sigma_omega_n[ineq][ispin][omega_n][m_index]
    std::vector<std::vector<std::vector<std::vector<std::complex<double>>>>> correlated_sigma_omega_n;

    //the lattice self energy 
    //lattice_sigma[ispin][omega_n][band_index]
    std::vector<std::vector<std::vector<std::complex<double>>>> lattice_sigma;

  };
}
