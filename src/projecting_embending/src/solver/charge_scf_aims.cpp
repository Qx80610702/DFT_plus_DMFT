#ifdef __FHIaims

#include "charge_scf_aims.h"
#include "../debug.h"
#include "../mpi_environment.h"

#include <cmath>
#include <sstream>
#include <string>
#include <iomanip>
#include <string>

#include <elsi.h>

namespace DFT_plus_DMFT
{
  void Charge_SCF_aims::output_charge_density(
        const int nks, 
        std::vector<std::vector<std::vector<
        std::complex<double>>>>& dens_mat_cmplx)
  {
    debug::codestamp("Charge_SCF_aims::output_charge_density");
    const int nbasis = (int) std::sqrt(dens_mat_cmplx[0][0].size());

    std::vector<int> k_map;
    int task_nks=0;
    for(int ik=0; ik<nks; ik++){
      if(ik%mpi_ntasks() != mpi_rank()) continue;  //k_points are divided acording to the process id
      task_nks +=1;
      k_map.push_back(ik);
    }

    elsi_rw_handle rwh;
    c_elsi_init_rw(&rwh,1,0,nbasis,0.0);

    for(int ik=0; ik<task_nks; ik++){
      for(int ispin=0; ispin<dens_mat_cmplx[0].size(); ispin++){
        std::stringstream ss;
        ss << "../DFT/D_spin_" 
           << std::setfill('0') << std::setw(2) << ispin+1
           << "_kpt_"
           << std::setfill('0') << std::setw(6) << k_map[ik]+1
           << ".csc";

        std::string file = ss.str();

        c_elsi_write_mat_complex(rwh, &file[0], reinterpret_cast<double _Complex*>(dens_mat_cmplx[ik][ispin].data()));
      }
    }

    c_elsi_finalize_rw(rwh);

    return;
  }
}

#endif
