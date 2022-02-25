#include "charge_scf.h"
#include "../debug.h"
#include "../mpi_environment.h"

#include <omp.h>
#include <mpi.h>

namespace DFT_plus_DMFT
{
  void Charge_SCF::update_char_dens(
        const int axis_flag,
        DFT_output::KS_bands& band, 
        DFT_output::atoms_info& atom, 
        DFT_plus_DMFT::projector& proj,
        DMFT::self_energy& sigma,
        DMFT::input_info& in,
        DFT_plus_DMFT::Hilbert_space& space)
  {
    debug::codestamp("Charge_SCF::update_char_dens");

    const int nks=band.nk();
    const int nspin=band.nspins();
    const std::vector<std::vector<int>>& wb2ib = space.wbands2ibands();

    int task_nks=0;
    for(int ik=0; ik<band.nk(); ik++) 
    {
      if(ik%mpi_ntasks() != mpi_rank()) continue;  //k_points are divided acording to process id
      task_nks +=1;
    }

    //Allocation
    this->fik_DMFT.resize(nspin);
    for(int is=0; is<nspin; is++){
      this->fik_DMFT[is].resize(task_nks);
      for(int ik=0; ik<task_nks; ik++){
        this->fik_DMFT[is][ik].resize(wb2ib[is][wb2ib[is].back()]+1, 0.0);
      }
    }

    this->eva_fik_DMFT( axis_flag, band, 
      atom, proj, sigma, in, space );


    return;
  }

  void Charge_SCF::eva_fik_DMFT(
        const int axis_flag,
        DFT_output::KS_bands& band, 
        DFT_output::atoms_info& atom, 
        DFT_plus_DMFT::projector& proj,
        DMFT::self_energy& sigma,
        DMFT::input_info& in,
        DFT_plus_DMFT::Hilbert_space& space)
  {
    debug::codestamp("Charge_SCF::eva_fik_DMFT");

    const double beta = *(double*)in.parameter("beta");
    const int nomega = *(int*)in.parameter("n_omega");
    const int nks=band.nk();
    const int nspin=band.nspins();
    std::vector<double> k_weight = band.kweight();
    const int magnetism = *(int*)in.parameter("magnetism");
    const std::vector<int>& wbands=space.Wbands();
    const int n_valence = space.valence();
    const std::vector<std::vector<int>>& wb2ib = space.wbands2ibands();
    std::vector<std::vector<bool>>& corb_flag = space.correction_flag();

    std::vector<int> k_map;
    int task_nks=0;
    for(int ik=0; ik<band.nk(); ik++){
      if(ik%mpi_ntasks() != mpi_rank()) continue;  //k_points are divided acording to the process id
      task_nks +=1;
      k_map.push_back(ik);
    }

    if(nspin==1 && !band.soc()) 
      for(int ik=0; ik<k_weight.size(); ik++)
        k_weight[ik] *= 2.0;

    //Uncorrelated bands
    for(int is=0; is<nspin; is++)
      for(int ik=0; ik<task_nks; ik++)
        for(int iband=0; iband<=wb2ib[is].back(); iband++)
          if(!corb_flag[is][iband]) this->fik_DMFT[is][ k_map[ik] ][iband] = k_weight[ k_map[ik] ];

    //Correlated bands

    return;
  }

  void Charge_SCF::eva_char_dens()
  {
    debug::codestamp("Charge_SCF::eva_char_dens");


    return;
  }
  
}