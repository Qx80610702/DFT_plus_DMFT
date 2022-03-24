#include "KS_bands.h"
#include "../debug.h"
#include "../mpi_environment.h"

#include <fstream>
#include <iostream>
#include <cstdlib>   //Use exit function
#include <iomanip>   //use setw and setprecison function

namespace DFT_output
{
  bool KS_bands::read()
  {
    debug::codestamp("KS_bands::read");
    // if(mpi_rank()==0) std::cout << "Reading bands information ......" << std::endl;

    int val;
    int spin_index,bands_index,k_index;

    std::ifstream if_bands("dft/outputs_to_DMFT/bands.dat", std::ios::in);

    if (!if_bands)  
	  {
	  	std::cout << "Fail to oepn dft/outputs_to_DMFT/bands.dat" << std::endl;
      std::exit(EXIT_FAILURE);
    }

    if_bands.seekg(0);      //set the position at the beginning of the file
    if_bands >> this->flag_SOC;
    if_bands.ignore(150,'\n');

    if(flag_SOC) // SOC
    {

    }
    else //non_SOC
    {
      if_bands >> this->tot_elec_num;
      if_bands >> this->nspin;
      if_bands >> this->nbands;
      if_bands >> this->kpoints; 
      if_bands >> this->Efermi;
      if_bands.ignore(150, '\n');

      this->eigen_values.resize(this->nspin);
      this->DFT_occ_numbers.resize(this->nspin);
      for(int ispin=0; ispin<nspin; ispin++)
      {
        this->eigen_values[ispin].resize(this->kpoints);
        this->DFT_occ_numbers[ispin].resize(this->kpoints);
        for(int ik=0; ik<kpoints; ik++)
        {
          this->eigen_values[ispin][ik].resize(this->nbands);
          this->DFT_occ_numbers[ispin][ik].resize(this->nbands);
        }
      }

      while(if_bands.good())
      {
        if_bands >> spin_index;
        if(if_bands.eof()) break;
        if_bands >> bands_index;
        if_bands >> k_index;

        if_bands >> this->eigen_values[spin_index][k_index][bands_index];
        if_bands >> this->DFT_occ_numbers[spin_index][k_index][bands_index];

        if_bands.ignore(150, '\n');
        if(if_bands.eof()) break;  //Check whether end of file is reached 
      }
      if_bands.close();

    }

    std::ifstream ifs("dft/outputs_to_DMFT/k_weight.dat", std::ios::in);

    if (!ifs)  
	  {
	  	std::cout << "Fail to oepn dft/outputs_to_DMFT/k_weight.dat" << std::endl;
      std::exit(EXIT_FAILURE);
    }

    ifs.seekg(0);      //set the position at the beginning of the file
    ifs >> val;
    this->k_weight.resize(val);
    ifs.ignore(150,'\n');

    while(ifs.good())
    {
      ifs >> val;
      if(ifs.eof()) break;
      ifs >> this->k_weight[val];
      ifs.ignore(150, '\n');

      if(ifs.eof()) break;  //Check whether end of file is reached 
    }
    ifs.close();

    return true;
  }

  void KS_bands::out()
  {
    debug::codestamp("KS_bands::out");

    std::string file;

    //===============================
    //  bands.dat
    //===============================
    file="bands.dat";
    std::ofstream ofs(file.c_str(),std::ios::out);
    if(!ofs)
    {
      std::cout << "Fail to oepn " << file.c_str() << std::endl;
      std::exit(EXIT_FAILURE);
    }

    ofs << this->flag_SOC << '\n';

    ofs << this->nspin << std::setw(6) << this->nbands << std::setw(6) 
        << std::setw(6) << this->kpoints 
        << std::setw(10) << std::setprecision(6) << this->Efermi << '\n';

    for(int ispin=0; ispin<this->nspin; ispin++)
    {
      for(int ib=0; ib<this->nbands; ib++)
      {
        for(int ik=0; ik<this->kpoints; ik++)
        {
          ofs << ispin << std::setw(6) << ib << std::setw(6) << ik
          << std::setw(18) << std::fixed << std::setprecision(9) << this->eigen_values[ispin][ik][ib] << std::endl; 
        }
      }
    }

    ofs.close();
  }
  
}
