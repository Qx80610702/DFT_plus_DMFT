//=================================================================
//    Computer Physics Communications 215 (2017) 128–136,
//    https://doi.org/10.1016/j.cpc.2017.01.003
//    The energy unit used is eV
//=================================================================
#pragma once

#include "../para/KS_bands.h"
#include "../para/correlated_atoms.h"
#include "../para/input.h"

#include <vector>
#include <complex>

namespace DMFT
{
  class ALPS_CTHYB
  {
    public:
    ALPS_CTHYB(){};
    ~ALPS_CTHYB(){};

    void read_last_step(const int istep, 
          DFT_output::KS_bands& band,
          DMFT::input_info& in, 
          DFT_output::atoms_info& atom,
          std::vector<std::vector<std::
          vector<std::vector<
          std::complex<double>>>>>& Gf_qmc,
          std::vector<std::vector<std::
          vector<std::vector<
          std::complex<double>>>>>& Weiss,
          std::vector<std::vector<std::
          vector<std::vector<
          std::complex<double>>>>>& Gf_save);

    void output(const int istep, const double mu, DMFT::input_info& in, 
          DFT_output::atoms_info& atom, DFT_output::KS_bands& band,
          std::vector<std::vector<std::vector<std::complex<double>>>>& Eimp,
          std::vector<std::vector<std::vector<std::vector<std::complex<double>>>>>& Gf_in,
          std::vector<std::vector<std::vector<std::vector<std::complex<double>>>>>& Weiss,
          std::vector<std::vector<std::vector<std::vector<std::complex<double>>>>>& hyb_omega);

    void out_sigma_last_step(
          const int istep, DFT_output::KS_bands& band,
          DMFT::input_info& in, DFT_output::atoms_info& atom,
          const  std::vector<std::vector<std::vector<
          std::vector<std::complex<double>>>>>& sigma);

    //Interfaces
    std::vector<std::vector<std::vector<std::vector<std::complex<double>>>>>&
        hybridization_func_tau(){return hyb_tau;}

    private:
    //hyb_tau[ineq][is][i_tau][m_index];
    std::vector<std::vector<std::vector<std::vector<std::complex<double>>>>> hyb_tau;


  };
}
