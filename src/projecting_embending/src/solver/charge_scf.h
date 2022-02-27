
#ifdef __FHIaims
#include "charge_scf_aims.h"
#endif

#include "../para/KS_bands.h"
#include "../para/correlated_atoms.h"
#include "self_energy.h"
#include "../para/input.h"
#include "projector.h"
#include "Hilbert_space.h"

#include <vector>
#include <complex>
#include <string>

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

    public:
    void update_char_dens(
        const int axis_flag,
        const double mu,
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

    void output_char_dense(
        const int dft_solver,
        const int nks );

    void prepare_nscf_dft(const int dft_solver);

    private:
    std::vector<std::vector<std::vector<std::complex<double>>>> dens_mat;       //dens_mat[ik][ispin][nbasis*nbasis]
    std::vector<std::vector<std::vector<double>>> fik_DMFT;       //fik_DMFT[is][ik][iband]

  };
}