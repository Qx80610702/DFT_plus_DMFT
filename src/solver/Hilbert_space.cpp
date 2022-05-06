#include "Hilbert_space.h"
#include "../mpi_environment.h"
#include "../debug.h"
#include "../constants.h"
#include "../global_variables.h"

#include <fstream>
#include <iostream>
#include <iomanip>
#include <cstdlib>   //Use exit function

namespace DFT_plus_DMFT
{
  void Hilbert_space::KS_bands_window(
        const int charge_step,
        DFT_output::KS_bands& band, 
        DFT_output::atoms_info& atom,
        DMFT::input_info& in)
  {
    debug::codestamp("Hilbert_space::KS_bands_window");

    this->nspin = band.nspins();
    this->n_KS_bands = band.nband();
    this->n_k_points = band.nk();
    const std::vector<std::vector<std::vector<double>>>& KS_eigenvalues = band.enk();

    // GlobalV::ofs_running << "Up energy   : " << this->energy_up_down[1] << " Hartree\n";
    // GlobalV::ofs_running << "Down energy : " << this->energy_up_down[0] << " Hartree\n";

    this->n_valence = band.n_electrons();

    this->n_wbands.resize(this->nspin);
    this->wbands_ibands.resize(this->nspin);
    this->ibands_wbands.resize(this->nspin);

    this->sigma_correction.resize(this->nspin);
    for(int is=0; is<this->nspin; is++)
      this->sigma_correction.at(is).resize(this->n_KS_bands,false);

    this->eigen_values.resize(this->nspin);
    
    if(charge_step==1) this->down_folding_first_charge_step(in, band);
    else this->read_corr_bands_windows(band);

    for(int is=0; is<this->nspin; is++){
      this->eigen_values[is].resize(band.nk());
      for(int ik=0; ik<band.nk(); ik++){
        this->eigen_values[is][ik].resize(this->n_wbands[is]);
        for(int iband=0; iband<this->n_wbands[is]; iband++)
          this->eigen_values[is][ik][iband] = 
          KS_eigenvalues[is][ik][this->wbands_ibands[is][iband]];
      }
    }//is

    return;
  }

  void Hilbert_space::down_folding_first_charge_step(
                      DMFT::input_info& in,
                      DFT_output::KS_bands& band )
  {
    debug::codestamp("Hilbert_space::down_folding_first_charge_step");

    const std::vector<std::vector<std::vector<double>>>& KS_eigenvalues = band.enk();
    const auto& Ewindow = in.window();

    this->energy_up_down.resize(2);
    this->energy_up_down[1] = band.Fermi_level()+Ewindow[1]/GlobalC::Hartree_to_eV;
    this->energy_up_down[0] = band.Fermi_level()+Ewindow[0]/GlobalC::Hartree_to_eV;

    double max_en, min_en;

    for(int is=0; is<this->nspin; is++){
      int upper_band_index;
      int lower_band_index;
      int iband;

      bool get_lower_band = false;
      iband=0;
      while(!get_lower_band)
      {
        if(iband>=this->n_KS_bands)
        {
          GlobalV::ofs_error << "Error in finding lower band" <<std::endl;
          std::exit(EXIT_FAILURE);
        }

        for(int ik=0; ik<this->n_k_points; ik++)
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
          GlobalV::ofs_error << "Error in finding upper band" <<std::endl;
          std::exit(EXIT_FAILURE);
        }

        for(int ik=0; ik<this->n_k_points; ik++)
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

      //Find maximum and minimum energy of bands in the window
      if(is==0){
        min_en = KS_eigenvalues[0][0][ lower_band_index ];
        max_en = KS_eigenvalues[0][0][ upper_band_index ];
      }

      for(int ik=0; ik<this->n_k_points; ik++){
        if( min_en > KS_eigenvalues[is][ik][ lower_band_index ] ) min_en = KS_eigenvalues[is][ik][ lower_band_index ];
        if( max_en < KS_eigenvalues[is][ik][ upper_band_index ] ) max_en = KS_eigenvalues[is][ik][ upper_band_index ];
      }

      
      GlobalV::ofs_running << "\n==============Correlated Kohn-Sham bands subset==============\n";
      if(this->nspin==1)
      {
        GlobalV::ofs_running << "Bands number of spin up:   " << lower_band_index << " ~ " << upper_band_index << std::endl;
        GlobalV::ofs_running << "Bands number of spin down: " << lower_band_index << " ~ " << upper_band_index << std::endl;
      }
      else if(this->nspin==2)
      {
        if(is==0) GlobalV::ofs_running << "Bands number of spin up:   " << lower_band_index << " ~ " << upper_band_index << std::endl;
        else GlobalV::ofs_running << "Bands number of spin down: " << lower_band_index << " ~ " << upper_band_index << std::endl;
      }
      

      this->n_wbands.at(is) = upper_band_index-lower_band_index+1;
      this->wbands_ibands.at(is).resize(upper_band_index-lower_band_index+1);
      this->ibands_wbands.at(is).resize(this->n_KS_bands,-1);

      for(iband=lower_band_index; iband<=upper_band_index; iband++)
      {
        this->sigma_correction.at(is).at(iband) = true;
        this->wbands_ibands.at(is).at(iband-lower_band_index) = iband;
        this->ibands_wbands.at(is).at(iband) = iband-lower_band_index;
      }

      if(this->nspin==1) this->n_valence -= 2*lower_band_index;
      else this->n_valence -= lower_band_index;
    }//is

    this->energy_up_down[1] = max_en;
    this->energy_up_down[0] = min_en;
    
    return;
  }

  void Hilbert_space::read_corr_bands_windows(
                      DFT_output::KS_bands& band )
  {
    debug::codestamp("Hilbert_space::read_corr_bands_windows");

    const std::vector<std::vector<std::vector<double>>>& KS_eigenvalues = band.enk();

    std::vector<int> upper_band_index(this->nspin);
    std::vector<int> lower_band_index(this->nspin);

    char word[200], word_low[200];

    std::ifstream ifs("DMFT_running.log", std::ios::in);
    if(!ifs){
      GlobalV::ofs_error << "Fail to oepn DMFT_running.log" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    ifs.seekg(0);      //set the position at the beginning of the file

    while(ifs.good())
    {
      ifs.getline(word, 200);   
      if(ifs.eof()) break;
      
      DMFT::input_info::strtolower(word, word_low);

      std::string line(word_low);

      if(!line.empty()){
        line.erase(0, line.find_first_not_of(" "));
        line.erase(line.find_last_not_of(" ") + 1);
      }

      if(!line.empty()){
        size_t pos_up=line.find("bands number of spin up:");
        size_t pos_dw=line.find("bands number of spin down:");

        if(pos_up!=std::string::npos || pos_dw!=std::string::npos){
          std::string index_str = line.substr(line.find_first_of(':')+1);
          index_str.erase(0, index_str.find_first_not_of(' '));
          index_str.erase(index_str.find_last_not_of(' ') + 1);

          std::string index_dw_str = index_str.substr(0,index_str.find_first_of('~'));
          index_dw_str.erase(0, index_dw_str.find_first_not_of(' '));
          index_dw_str.erase(index_dw_str.find_last_not_of(' ') + 1);

          std::string index_up_str = index_str.substr(index_str.find_first_of('~')+1);
          index_up_str.erase(0, index_up_str.find_first_not_of(' '));
          index_up_str.erase(index_up_str.find_last_not_of(' ') + 1);

          int index_dw = atoi(index_dw_str.c_str());
          int index_up = atoi(index_up_str.c_str());

          if(this->nspin==1 && pos_up!=std::string::npos){
            upper_band_index[0] = index_up;
            lower_band_index[0] = index_dw;
          }
          else if(this->nspin==2){
            if(pos_up!=std::string::npos){
              upper_band_index[0] = index_up;
              lower_band_index[0] = index_dw;
            }
            else if(pos_dw!=std::string::npos){
              upper_band_index[1] = index_up;
              lower_band_index[1] = index_dw;
            }
          }         
        }
      }
    }
    ifs.close();

    //Find maximum and minimum energy of bands in the window
    double edw = KS_eigenvalues[0][0][ lower_band_index[0] ];
    double eup = KS_eigenvalues[0][0][ upper_band_index[0] ];
    for(int is=0; is<this->nspin; is++){
      for(int ik=0; ik<this->n_k_points; ik++){
        if( edw > KS_eigenvalues[is][ik][ lower_band_index[0] ] ) edw = KS_eigenvalues[is][ik][ lower_band_index[0] ];
        if( eup < KS_eigenvalues[is][ik][ upper_band_index[0] ] ) eup = KS_eigenvalues[is][ik][ upper_band_index[0] ];
      }
    }

    this->energy_up_down.resize(2);
    this->energy_up_down[1] = eup;
    this->energy_up_down[0] = edw;

    for(int is=0; is<this->nspin; is++){
      this->n_wbands.at(is) = upper_band_index[is]-lower_band_index[is]+1;
      this->wbands_ibands.at(is).resize(upper_band_index[is]-lower_band_index[is]+1);
      this->ibands_wbands.at(is).resize(this->n_KS_bands,-1);

      for(int iband=lower_band_index[is]; iband<=upper_band_index[is]; iband++){
        this->sigma_correction.at(is).at(iband) = true;
        this->wbands_ibands.at(is).at(iband-lower_band_index[is]) = iband;
        this->ibands_wbands.at(is).at(iband) = iband-lower_band_index[is];
      }

      if(this->nspin==1) this->n_valence -= 2*lower_band_index[is];
      else this->n_valence -= lower_band_index[is];
    }//is

    for(int is=0; is<this->nspin; is++){
      GlobalV::ofs_running << "\n==============Correlated Kohn-Sham bands subset==============\n";
      if(this->nspin==1){
        GlobalV::ofs_running << "Bands number of spin up:   " << lower_band_index[0] << " ~ " << upper_band_index[0] << std::endl;
        GlobalV::ofs_running << "Bands number of spin down: " << lower_band_index[0] << " ~ " << upper_band_index[0] << std::endl;
      }
      else if(this->nspin==2){
        if(is==0) GlobalV::ofs_running << "Bands number of spin up:   " << lower_band_index[0] << " ~ " << upper_band_index[0] << std::endl;
        else GlobalV::ofs_running << "Bands number of spin down: " << lower_band_index[1] << " ~ " << upper_band_index[1] << std::endl;
      }  
    }

    return;
  }

}