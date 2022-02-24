#ifdef __FHIaims

#include "charge_update_aims.h"
#include "../debug.h"

#include <elsi.h>

namespace DFT_plus_DMFT
{
  void Charge_update_aims::output_charge_density(
    std::string file, std::vector<std::complex<double>>& dens_cmplx)
  {
    debug::codestamp("Charge_update_aims::output_charge_density");

    elsi_rw_handle rwh;

    c_elsi_init_rw(&rwh,1,0,0,0.0);

    c_elsi_write_mat_complex(rwh, &file[0], reinterpret_cast<double _Complex*>(dens_cmplx.data()));

    c_elsi_finalize_rw(rwh);

    return;
  }
}

#endif
