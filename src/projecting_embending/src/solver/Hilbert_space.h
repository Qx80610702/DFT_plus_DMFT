#pragma once

#include "../para/KS_bands.h"
#include "../para/correlated_atoms.h"
#include "../para/input.h"

#include <complex>
#include <vector>

namespace DFT_plus_DMFT
{
  class Hilbert_space
  {
    public:
    Hilbert_space(){};
    ~Hilbert_space(){};

    //dectect which KS bands need self energy correction
    void KS_bands_window(
          const int charge_step,
          DFT_output::KS_bands& band, 
          DFT_output::atoms_info& atom, 
          DMFT::input_info& in );

    void down_folding_first_charge_step(
          DMFT::input_info& in,
          DFT_output::KS_bands& band );

    void read_corr_bands_windows(
          DFT_output::KS_bands& band );
    
    // interfaces
    std::vector<int>& Wbands(){return n_wbands;}
    std::vector<std::vector<int>>& wbands2ibands(){return wbands_ibands;}
    std::vector<std::vector<int>>& ibands2wbands(){return ibands_wbands;}
    std::vector<std::vector<bool>>& correction_flag(){return sigma_correction;}
    double valence(){return n_valence;}
    const std::vector<double>& Ener_window(){return energy_up_down;}
    const std::vector<std::vector<std::vector<double>>>& eigen_val(){return this->eigen_values;}

    private:

    int nspin;   
    double n_valence;               //Number of valence electrons  
    int n_KS_bands;
    int n_k_points;

    //The energy window of KS bands
    std::vector<double> energy_up_down;       // energy_up_down[0(dw) or 1(up)]

    //whether KS band[is][n_KS_bands] is in the window or not
    std::vector<std::vector<bool>> sigma_correction;

    //number of KS bands in the window; n_KS_bands_corrected[ispin]
    std::vector<int> n_wbands;
    std::vector<std::vector<std::vector<double>>> eigen_values;      //eigen_values[i_spin][ik][i_wbands]

    //mapping between the band indices of self energy corrected bands and the band indices of KS bands
    //ibands_wbands[is][n_KS_bands]
    //wbands_ibands[is][n_wbands]
    std::vector<std::vector<int>> wbands_ibands;
    std::vector<std::vector<int>> ibands_wbands;

  };
}
