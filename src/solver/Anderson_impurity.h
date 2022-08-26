#pragma once

#include "self_energy.h"
#include "../para/KS_bands.h"
#include "../para/correlated_atoms.h"
#include "projector.h"
#include "../para/input.h"
#include "alps_cthyb.h"
#include "alps_cthyb_segment.h"
#include "PACS_cthyb.h"
#include "rutgers_cthyb.h"
#include "iQIST_narcissus.h"
#include "Hilbert_space.h"
#include "coulomb_tensor.h"

#include <vector>
#include <complex>

namespace DMFT
{
  class impurity
  {
    public:
    impurity(){};
    ~impurity(){};

    self_energy sigma;
    PACS_CTHYB pacs;
    Rutgers_CTHYB Rutgers;
    ALPS_CTHYB ALPS_hyb;
    ALPS_CTHYB_SEGMENT ALPS_hyb_segment;
    IQIST_NARCISSUS  iQIST_narcissus;

    void evaluate_impurity_level(
          DFT_plus_DMFT::Hilbert_space& space,
          DFT_output::KS_bands& band, 
          DFT_plus_DMFT::projector& proj,
          DFT_output::atoms_info& atom);
    
    //solving Anderson impurity model on imaginary axis
    void evaluate_Weiss_hybridization_imag(
          DFT_output::KS_bands& band, 
          DFT_plus_DMFT::projector& proj,
          DFT_output::atoms_info& atom,
          DFT_plus_DMFT::Hilbert_space& space,
          const double mu,
          const int nomega,
          const int mag);

    void evaluate_delta_Weiss_tau(
          const int impurity_solver,
          DFT_output::atoms_info& atom,
          const int nspin,
          const double beta,
          const int ntau,
          const int nomega);
    
    void evaluate_local_occupation(
          DMFT::input_info& in, 
          DFT_output::KS_bands& band, 
          DFT_plus_DMFT::projector& proj,
          DFT_output::atoms_info& atom,
          DFT_plus_DMFT::Hilbert_space& space,
          const int nomega,
          const int mag);

    void read_self_energy(
          const bool restart,
          const int char_step,
          const int DMFT_step,
          const int impurity_solver, 
          DFT_output::KS_bands& band,
          DMFT::input_info& in, 
          DFT_output::atoms_info& atom );
    
    bool scf_condition(const int impurity_solver, 
                      DFT_output::KS_bands& band,
                      DFT_output::atoms_info& atom,
                      DMFT::input_info& in);

    void update_self_energy(DFT_output::KS_bands& band,
                      DFT_output::atoms_info& atom,
                      DMFT::input_info& in);

    void out(const int char_step,
            const int DMFT_step,
            const int impurity_solver, 
            const double mu, DFT_output::KS_bands& band,
            DMFT::input_info& in, DFT_output::atoms_info& atom,
            DMFT::coulomb_tensor& Umat);

    // interfaces
    std::vector<double>& delta_scf(){return scf_delta;}

    std::vector<std::vector<std::vector<std::vector<std::complex<double>>>>>&
        Matsubara_Gf(){return Green_fun_omega;}

    std::vector<std::vector<std::vector<std::vector<std::complex<double>>>>>&
        Matsubara_Gf_save(){return Green_fun_omega_save;}

    std::vector<std::vector<std::vector<std::vector<std::complex<double>>>>>&
        Weiss(){return Weiss_omega;}

    std::vector<std::vector<std::vector<std::vector<std::complex<double>>>>>& 
        hybdization_tau(const int impurity_solver);

    public: //static member function

    static void Fourier_trans_omega_tau(
            const double beta,
            const int ntau,
            const int nomega,
            const double* freq,
            const double* tau_mesh,
            const std::complex<double>* fw,
            std::complex<double>* ftau);

    private:
    //Current step, impurity_level[ineq][is][m_index]
    std::vector<std::vector<std::vector<std::complex<double>>>> impurity_level;

    //hyb_omega[ineq][is][iomega][m_index]
    std::vector<std::vector<std::vector<std::vector<std::complex<double>>>>> hyb_omega;

    //Weiss_omega[ineq][is][iomega][m_index]
    std::vector<std::vector<std::vector<std::vector<std::complex<double>>>>> Weiss_omega;

    //Weiss_tau[ineq][is][itau][m_index]
    std::vector<std::vector<std::vector<std::vector<std::complex<double>>>>> Weiss_tau;

    //Green_fun_omega[ineq][is][iomega][m_index]
    std::vector<std::vector<std::vector<std::vector<std::complex<double>>>>> Green_fun_omega;

    //Green_fun_omega_save[ineq][is][iomega][m_index]
    std::vector<std::vector<std::vector<std::vector<std::complex<double>>>>> Green_fun_omega_save;

    std::vector<double> scf_delta;

  };
}
