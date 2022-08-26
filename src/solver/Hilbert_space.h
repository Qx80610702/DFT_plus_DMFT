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
          const int type,
          DFT_output::KS_bands& band, 
          DFT_output::atoms_info& atom, 
          DMFT::input_info& in );
    
    // interfaces
    std::vector<int>& Wbands(){return n_wbands;}
    std::vector<int>& DOSWbands(){return DOS_n_wbands;}
    std::vector<std::vector<int>>& wbands2ibands(){return wbands_ibands;}
    std::vector<std::vector<int>>& ibands2wbands(){return ibands_wbands;}
    std::vector<std::vector<int>>& correction_flag(){return sigma_correction;}
    double valence(){return n_valence;}
    const std::vector<double>& Ener_window(){return energy_up_down;}
    const std::vector<std::vector<std::vector<double>>>& eigen_val(){return this->eigen_values;}
    const std::vector<std::vector<std::vector<double>>>& DOS_eigen_val(){return this->DOS_eigen_values;}
    const std::vector<std::vector<int>>& corr_bands2DOS_bands(){return corr_bands_DOS_bands;}

    private:
    void deter_bands_window_first_charge_step(
          DMFT::input_info& in,
          DFT_output::KS_bands& band );

    void read_bands_windows(
          DFT_output::KS_bands& band );

    private:

    int nspin;   
    double n_valence;               //Number of valence electrons  
    int n_KS_bands;
    int n_k_points;

    //===========Projection================
    //The energy window of KS bands
    std::vector<double> energy_up_down;       // energy_up_down[0(dw) or 1(up)]

    //whether KS band[is][n_KS_bands] is in the window or not; 1 for in and 0 for not
    std::vector<std::vector<int>> sigma_correction;

    //number of KS bands in the window; n_KS_bands_corrected[ispin]
    std::vector<int> n_wbands;
    std::vector<std::vector<std::vector<double>>> eigen_values;      //eigen_values[i_spin][ik][i_wbands]

    //mapping between the band indices of self energy corrected bands and the band indices of KS bands
    //ibands_wbands[is][n_KS_bands]
    //wbands_ibands[is][n_wbands]
    std::vector<std::vector<int>> wbands_ibands;
    std::vector<std::vector<int>> ibands_wbands;

    //===========DOS=================
    //The energy window of KS bands
    std::vector<double> DOS_energy_up_down;       // energy_up_down[0(dw) or 1(up)]
    //whether KS band[is][n_KS_bands] is in the DOS window or not; 1 for in and 0 for not
    std::vector<std::vector<int>> DOS_band;
    //number of KS bands in the DOS_window; n_KS_bands_corrected[ispin]
    std::vector<int> DOS_n_wbands;
    std::vector<std::vector<std::vector<double>>> DOS_eigen_values;      //DOS_eigen_values[i_spin][ik][i_wbands]
    //mapping between the band indices of DOS bands and the band indices of KS bands
    std::vector<std::vector<int>> DOS_wbands_ibands;         //DOS_wbands_ibands[is][DOS_n_wbands]
    std::vector<std::vector<int>> DOS_ibands_wbands;         //DOS_ibands_wbands[is][n_KS_bands]

    std::vector<std::vector<int>> corr_bands_DOS_bands;      //corr_bands_DOS_bands[is][n_wbands] 

  };
}
