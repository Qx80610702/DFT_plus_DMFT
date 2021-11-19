//================================================================================
//    Phys. Rev. B, 2007,215:155113. https://doi.org/10.1103/PhysRevB.75.155113
//    Phys. Rev. B, 2010,81:195107. https://doi.org/10.1103/PhysRevB.81.195107
//    The energy unit used is eV
//=================================================================================
#pragma once

#include "../para/KS_bands.h"
#include "../para/correlated_atoms.h"
#include "../para/input.h"
#include "coulomb_tensor.h"

#include <vector>
#include <complex>

namespace DMFT
{
  class Rutgers_CTHYB
  {
    public:
    Rutgers_CTHYB(){};
    ~Rutgers_CTHYB(){};

    void read_last_step(const int istep, 
          DFT_output::KS_bands& band,
          DMFT::input_info& in, 
          DFT_output::atoms_info& atom,
          std::vector<std::vector<std::
          vector<std::vector<
          std::complex<double>>>>>& Gw_qmc,
          std::vector<std::vector<std::
          vector<std::vector<
          std::complex<double>>>>>& Gw_save,
          std::vector<std::vector<std::
          vector<std::vector<
          std::complex<double>>>>>& Sw );

    void output(const int istep, const double mu, DMFT::input_info& in, 
          DFT_output::atoms_info& atom, DFT_output::KS_bands& band, const std::vector<double>& freq,
          std::vector<std::vector<std::vector<std::complex<double>>>>& Eimp,
          std::vector<std::vector<std::vector<std::vector<std::complex<double>>>>>& Gf_in,
          std::vector<std::vector<std::vector<std::vector<std::complex<double>>>>>& Weiss,
          std::vector<std::vector<std::vector<std::vector<std::complex<double>>>>>& hyb_omega,
          DMFT::coulomb_tensor& Umat );

    void write_params(DMFT::input_info& in, const std::string file, const double mu);
    void write_delta_omega(DMFT::input_info& in, const std::string file,
          const int nspin, const int nomega, 
          const int m_tot, const std::vector<double>& freq,
          const std::vector<std::vector<
          std::vector<std::complex<double>>>>& hyb_omega);

    void write_cix(const int nspin, const int m_tot, const std::string file, 
        const std::vector<std::vector<std::complex<double>>>& Esplit,
        const std::vector<std::vector<std::vector<std::vector<double>>>>& Utensor);

    int site_occ(int n);
    double site_Sz(int n, const int nbath);
    int F_dagger_state(int n, const int ibath, const int nbath);
    double state_energy(const int n, const int m_tot, const int nspin, 
            const std::vector<std::vector<std::vector<std::vector<double>>>>& Utensor,
            const std::vector<std::vector<std::complex<double>>>& Eimp);

    double F_dagger_sign(int n, const int ibath);

    //========================interface=========================
    std::vector<std::vector<std::vector<std::complex<double>>>>&
          hopping_term(){return hopping_matrix;}

    std::vector<std::vector<std::vector<std::vector<std::complex<double>>>>>&
        hybridization_func_tau(){return hyb_tau;}

    private:
    //hopping_term[ineq][is][m_index]
    std::vector<std::vector<std::vector<std::complex<double>>>> hopping_matrix;

    //hyb_tau[ineq][is][i_tau][m_index];
    std::vector<std::vector<std::vector<std::vector<std::complex<double>>>>> hyb_tau;

  };
}
