#ifdef __FHIaims

#include <complex>
#include <vector>

namespace DFT_plus_DMFT
{
  class Charge_SCF_aims
  {
    public:
    Charge_SCF_aims(){;}
    ~Charge_SCF_aims(){;}

    public:
    void output_charge_density(
          const int nks, 
          std::vector<std::vector<std::vector<
          std::complex<double>>>>& dens_mat_cmplx);

    void prepare_nscf_dft();

  };
}

#endif