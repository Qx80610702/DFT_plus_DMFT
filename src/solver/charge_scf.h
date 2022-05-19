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

// Pulay DIIS method. 
//Refs.: P. Pulay. Chem. Phys. Lett. 73, 393(1980); G. Kress, J. Furthmuller. Phys. Rev. B 54, 11169(1996)
//Mixing step 8.
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
        const int mixing_step,
        const double delta_rho,
        const int nks, const int n_spin );

    void update_charge_density_matrix(
        const int axis_flag,
        DFT_plus_DMFT::chemical_potential& Mu,
        DFT_output::KS_bands& band, 
        DFT_output::atoms_info& atom, 
        DFT_plus_DMFT::projector& proj,
        DMFT::self_energy& sigma,
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

    void eva_k_densmat(
        const int dft_solver,
        const int ik,
        const int i_k_point,
        DFT_plus_DMFT::Hilbert_space& space,
        std::vector<std::vector<
        std::complex<double>>>& eigenvector,
        std::vector<std::vector<
        std::complex<double>>>& dense_cmplx );

    void read_charge_density(
      const bool initial_charge,
      const bool DMFT_charge );

    void read_initial_charge_density_matrix();

    bool charge_mixing(const int mix_step);

    void update_data(const int mix_step);

    void update_alpha(const int mix_step, std::vector<double>& alpha);

    void mixing_density_matrix(
      const int mix_step, 
      const std::vector<double>& alpha);

    void output_mixed_charge_density_matrix();

    void output_DMFT_charge_density_matrix();

    void prepare_nscf_dft();

    private:
    int flag_DFT_solver;
    double sc_delta_rho;
    double mixing_beta;
    int max_mixing_step;
    int nkpoints;
    int nspin;

    //Opt_DM_mat[istep][ik][ispin][nbasis*nbasis]
    std::deque<std::vector<std::vector<std::vector<std::complex<double>>>>> Opt_DM_mat;

    //Res_DM_mat[istep][ik][ispin][nbasis*nbasis]
    std::deque<std::vector<std::vector<std::vector<std::complex<double>>>>> Res_DM_mat;

    std::vector<std::vector<std::vector<double>>> fik_DMFT;                     //fik_DMFT[is][ik][iband]

    std::vector<std::vector<std::vector<std::complex<double>>>> dens_mat_out;   //dens_mat_ouy[ik][ispin][nbasis*nbasis]

  };
}