#include "Hilbert_space.h"
#include "../mpi_environment.h"
#include "../debug.h"
#include "../constants.h"

#include <iostream>
#include <iomanip>
#include <cstdlib>   //Use exit function

namespace DFT_plus_DMFT
{
  void Hilbert_space::KS_bands_window(
        DFT_output::KS_bands& band, 
        DFT_output::atoms_info& atom,
        DMFT::input_info& in)
  {
    debug::codestamp("Hilbert_space::KS_bands_window");

    this->nspin = band.nspins();
    this->n_KS_bands = band.nband();
    const std::vector<std::vector<std::vector<double>>>& KS_eigenvalues = band.enk();
    const auto& Ewindow = in.window();

    this->energy_up_down.resize(2);
    this->energy_up_down[1] = band.Fermi_level()+Ewindow[1]/Hartree_to_eV;
    this->energy_up_down[0] = band.Fermi_level()+Ewindow[0]/Hartree_to_eV;
    
    // std::cout << "Up energy   : " << this->energy_up_down[1] << " Hartree\n";
    // std::cout << "Down energy : " << this->energy_up_down[0] << " Hartree\n";

    this->n_valence = band.n_electrons();

    this->n_wbands.resize(this->nspin);
    this->wbands_ibands.resize(this->nspin);
    this->ibands_wbands.resize(this->nspin);

    this->sigma_correction.resize(this->nspin);
    for(int is=0; is<this->nspin; is++)
      this->sigma_correction.at(is).resize(this->n_KS_bands,false);

    this->eigen_values.resize(this->nspin);
    
    for(int is=0; is<this->nspin; is++)
    {
      int upper_band_index;
      int lower_band_index;
      int iband;

      bool get_lower_band = false;
      iband=0;
      while(!get_lower_band)
      {
        if(iband>=this->n_KS_bands)
        {
          std::cout << "Error in finding lower band" <<std::endl;
          std::exit(EXIT_FAILURE);
        }

        for(int ik=0; ik<band.nk(); ik++)
        {
          if(KS_eigenvalues[is][ik][iband]>=this->energy_up_down[0])
          {
            get_lower_band=true;
            lower_band_index=iband;
            break;
          }
        }
        iband++;
      }

      bool get_upper_band = false;
      iband=this->n_KS_bands-1;
      while(!get_upper_band)
      {
        if(iband<0)
        {
          std::cout << "Error in finding upper band" <<std::endl;
          std::exit(EXIT_FAILURE);
        }

        for(int ik=0; ik<band.nk(); ik++)
        {
          if(KS_eigenvalues[is][ik][iband]<=this->energy_up_down[1])
          {
            get_upper_band=true;
            upper_band_index=iband;
            break;
          }
        }  
        iband--;
      }

      this->n_wbands.at(is) = upper_band_index-lower_band_index+1;
      this->wbands_ibands.at(is).resize(upper_band_index-lower_band_index+1);
      this->ibands_wbands.at(is).resize(this->n_KS_bands,-1);
      
      if(mpi_rank()==0)
      {
        std::cout << "\n==============Downfolding==============\n";
        std::cout << "Energy window: " 
                  << std::setw(11) << std::fixed << std::setprecision(6)
                  << this->energy_up_down[0]*Hartree_to_eV
                  << "  ~ "
                  << std::setw(11) << std::fixed << std::setprecision(6)
                  << this->energy_up_down[1]*Hartree_to_eV << '\n';

        if(this->nspin==1)
        {
          std::cout << "Bands number of spin up:   " << lower_band_index << " ~ " << upper_band_index << std::endl;
          std::cout << "Bands number of spin down: " << lower_band_index << " ~ " << upper_band_index << std::endl;
        }
        else if(this->nspin==2)
        {
          if(is==0) std::cout << "Bands number of spin up:   " << lower_band_index << " ~ " << upper_band_index << std::endl;
          else std::cout << "Bands number of spin down: " << lower_band_index << " ~ " << upper_band_index << std::endl;
        }
      }

      for(iband=lower_band_index; iband<=upper_band_index; iband++)
      {
        this->sigma_correction.at(is).at(iband) = true;
        this->wbands_ibands.at(is).at(iband-lower_band_index) = iband;
        this->ibands_wbands.at(is).at(iband) = iband-lower_band_index;
      }

      if(this->nspin==1) this->n_valence -= 2*lower_band_index;
      else this->n_valence -= lower_band_index;

      this->eigen_values[is].resize(band.nk());
      for(int ik=0; ik<band.nk(); ik++)
      {
        this->eigen_values[is][ik].resize(this->n_wbands[is]);
        for(iband=0; iband<this->n_wbands[is]; iband++)
          this->eigen_values[is][ik][iband] = 
          KS_eigenvalues[is][ik][this->wbands_ibands[is][iband]];
      }

    }//is

    return;
  }

}
