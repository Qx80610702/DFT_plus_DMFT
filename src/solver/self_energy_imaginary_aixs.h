#pragma once

#include "../para/correlated_atoms.h"
#include "projector.h"

#include <complex>
#include <vector>
#include <map>

namespace DMFT
{
  class self_energy_imaginary_aixs
  {
    public:
    self_energy_imaginary_aixs(){};
    ~self_energy_imaginary_aixs(){};

    bool read_sigma_save(const bool SOC, const int nspin, DFT_output::atoms_info& atom);

    //===============================================
    //              Test
    //===============================================
    void out();    //test whether reading worked corrected

    //===============================================
    //       interface
    //===============================================
    inline std::vector<std::vector<std::vector<std::vector<std::complex<double>>>>>& 
            sigma_new_access(){return local_sigma_new;}

    inline std::vector<std::vector<std::vector<std::vector<std::complex<double>>>>>& 
            sigma_save_access(){return local_sigma_save;}

    inline std::vector<std::vector<std::vector<std::vector<std::complex<double>>>>>& 
            correlated_sigma_access(){return correlated_sigma_omega_n;}

    inline std::vector<std::vector<std::vector<std::complex<double>>>>&
            lattice_sigma_access(){return lattice_sigma;}

    std::vector<double>& Matsubara_freq(){return imag_frequency;}

    int& nomega(){return n_omega;} 
    double& inverse_T(){return beta;}

    public: //static menber function
    static void second_order_self_energy();

    private: 
    //The Matsubara frquence
    std::vector<double> imag_frequency;

    //the local self energy of current DMFT step which is calculated by Dyson equation
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

    std::map<int, bool> self_energy_convergence;    //whether self-energy converged

    bool flag_initital_guess;

    double beta;
    int n_omega;
  };
}
