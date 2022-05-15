
#ifdef __FHIaims
#include "charge_scf_aims.h"
#endif

#include "../para/KS_bands.h"
#include "../para/correlated_atoms.h"
#include "self_energy.h"
#include "../para/input.h"
#include "projector.h"
#include "Hilbert_space.h"
#include "chemical_potential.h"

#include <vector>
#include <complex>
#include <string>
#include <deque>

namespace DFT_plus_DMFT
{
  class Charge_SCF
  {
    public:
    Charge_SCF(){;}
    ~Charge_SCF(){;}

    public:
    #ifdef __FHIaims
    Charge_SCF_aims char_scf_aims;
    #endif

    auto& char_ref();

    public:
    void init(const int DFT_solver,
        const double mixing_parameter,
        const int nks, const int n_spin );

    void update_charge_density_matrix(
        const int axis_flag,
        DFT_plus_DMFT::chemical_potential& Mu,
        DFT_output::KS_bands& band, 
        DFT_output::atoms_info& atom, 
        DFT_plus_DMFT::projector& proj,
        DMFT::self_energy& sigma,
        DMFT::input_info& in,
        DFT_plus_DMFT::Hilbert_space& space );

    void eva_new_char_dens(
        const int axis_flag,
        DFT_plus_DMFT::chemical_potential& Mu,
        DFT_output::KS_bands& band, 
        DFT_output::atoms_info& atom, 
        DFT_plus_DMFT::projector& proj,
        DMFT::self_energy& sigma,
        DMFT::input_info& in,
        DFT_plus_DMFT::Hilbert_space& space );

    void eva_fik_DMFT_imag_axis(
        const double mu,
        DFT_output::KS_bands& band, 
        DFT_output::atoms_info& atom, 
        DFT_plus_DMFT::projector& proj,
        DMFT::self_energy& sigma,
        DMFT::input_info& in,
        DFT_plus_DMFT::Hilbert_space& space );

    void eva_fik_DMFT_real_axis(
        const double mu,
        DFT_output::KS_bands& band, 
        DFT_output::atoms_info& atom, 
        DFT_plus_DMFT::projector& proj,
        DMFT::self_energy& sigma,
        DMFT::input_info& in,
        DFT_plus_DMFT::Hilbert_space& space );

    void eva_k_densmat(
        const int dft_solver,
        const int nspin,
        const int ik,
        const int i_k_point,
        DFT_plus_DMFT::Hilbert_space& space,
        std::vector<std::vector<
        std::complex<double>>>& eigenvector,
        std::vector<std::vector<
        std::complex<double>>>& dense_cmplx );

    void read_charge(
      const bool initial_charge,
      const bool DMFT_charge );

    void read_charge_density_matrix(const int nks);

    void charge_mixing(const int mix_step, double& charge_change);

    void update_data(const int mix_step);

    void update_alpha(const int mix_step, std::vector<double>& alpha);

    void update_density(
      const int mix_step, 
      std::vector<double>& alpha,
      double& charge_change );

    void output_charge_density_matrix(const int nks);

    void prepare_nscf_dft();

    private:
    int flag_DFT_solver;
    double mixing_beta;
    int nkpoints;
    int nspin;

    //DM_mat_pulay[istep][ik][ispin][nbasis*nbasis]
    std::deque<std::vector<std::vector<std::vector<std::complex<double>>>>> DM_mat_pulay;

    std::vector<std::vector<std::vector<double>>> fik_DMFT;                     //fik_DMFT[is][ik][iband]

    std::vector<std::vector<std::vector<std::complex<double>>>> dens_mat;  //dens_mat_last[ik][ispin][nbasis*nbasis]
    std::vector<std::vector<std::vector<std::complex<double>>>> dens_mat_last;  //dens_mat_last[ik][ispin][nbasis*nbasis]
    
    // std::vector<std::vector<std::vector<double>>> fik_test;       //fik_test[is][ik][iband]

    // Rstep: iteration step for Rrho
    // dRstep: iteration step for dRrho

  };
}