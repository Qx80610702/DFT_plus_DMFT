#ifdef __FHIaims

#include <string>
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
    void output_charge_density(std::string file, std::vector<std::complex<double>>& dens_cmplx);

  };
}

#endif