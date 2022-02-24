
#ifdef __FHIaims
#include "charge_scf_aims.h"
#endif

#include <vector>
#include <complex>

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
    void update_char_dens();
    void eva_char_dens();

    private:
    std::vector<std::vector<double>> dens_mat;       //dens_mat[ispin][nbasis*nbasis]

  };
}