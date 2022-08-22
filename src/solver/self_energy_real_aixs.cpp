#include "self_energy_real_aixs.h"
#include "../debug.h"
#include "../constants.h"
#include "../global_variables.h"

#include <fstream>
#include <iostream>
#include <string>
#include <cstring>
#include <cstdlib>
#include <sstream>

namespace DMFT
{
  void self_energy_real_aixs::read_AC_sigma(
        DFT_output::KS_bands& band,
        DMFT::input_info& in, 
        DFT_output::atoms_info& atom )
  {
    debug::codestamp("self_energy_real_aixs::read_AC_sigma");

    const int ineq_num = atom.inequ_atoms();
    const std::vector<int>& norb_sub = atom.iatom_norb();
    const int nspin = band.nspins();
    const int magnetism = *(int*)in.parameter("magnetism");

    const std::complex<double> zeroc(0.0,0.0);
    double omega;

    this->read_maxent_parms();

    //Allocation
    if(this->local_sigma_new.empty())
    {
      this->local_sigma_new.resize(ineq_num);
      for(int ineq=0; ineq<ineq_num; ineq++)
      {
        const int iatom = atom.ineq_iatom(ineq);
        const int m_tot = norb_sub[iatom];

        this->local_sigma_new[ineq].resize(nspin);
        for(int is=0; is<nspin; is++)
        {
          this->local_sigma_new[ineq][is].resize(this->n_omega);
          for(int iomega=0; iomega<this->n_omega; iomega++)
          {
            this->local_sigma_new[ineq][is][iomega].resize(m_tot*m_tot,zeroc);
          }//is
        }//i_omega
      }//ineq
    }

    if(this->real_frequency.empty())
      this->real_frequency.resize(this->n_omega,0.0);

    for(int ineq=0; ineq<ineq_num; ineq++)
    {
      const int iatom = atom.ineq_iatom(ineq);
      const int m_tot=norb_sub[iatom];

      std::stringstream dir_ss;
      dir_ss << "./self-energy/impurity" << ineq;
      std::string file= dir_ss.str() + "/Sigma_omega.dat";

      std::vector<std::vector<double>> Gw_real(2);
      std::vector<std::vector<double>> Gw_im(2);
      for(int is=0; is<2; is++){
        Gw_real[is].resize(m_tot);
        Gw_im[is].resize(m_tot);
      }

      std::ifstream ifs(file.c_str(), std::ios::in);
      if (!ifs){
	    	std::cerr << "Fail to oepn " << file << std::endl;
        std::exit(EXIT_FAILURE);}
      ifs.seekg(0);    //set the position at the beginning of the file

      int count=0;
      while(ifs.good())
      {
        ifs >> omega;
        if(ifs.eof()) break; //Check whether end of file is reached
        this->real_frequency[count] = omega/GLC::Hartree_to_eV;

        for(int is=0; is<2; is++)
          for(int m=0; m<m_tot; m++)
          {
            ifs >> Gw_real[is][m];
            ifs >> Gw_im[is][m];
          }
        ifs.ignore(150,'\n');

        for(int is=0; is<nspin; is++)
        {
          for(int m=0; m<m_tot; m++)
          {
            if(magnetism==3 || magnetism==4) //none magnetic or paramamagnetic
            {
              if(nspin==1) 
              {
                Gw_real[0][m] = (Gw_real[0][m]+Gw_real[1][m])/2.0;
                Gw_im[0][m] = (Gw_im[0][m]+Gw_im[1][m])/2.0;
              }
              else
              {
                Gw_real[0][m] = (Gw_real[0][m]+Gw_real[1][m])/2.0;
                Gw_real[1][m] = Gw_real[0][m];
                Gw_im[0][m] = (Gw_im[0][m]+Gw_im[1][m])/2.0;
                Gw_im[1][m] = Gw_im[0][m];
              }
            }
          }
          DFT_output::atoms_info::symmetry_operation_vector<double>(
              atom.local_sym(), atom.L(ineq), m_tot, &Gw_real[is][0] );

          DFT_output::atoms_info::symmetry_operation_vector<double>(
              atom.local_sym(), atom.L(ineq), m_tot, &Gw_im[is][0] );
        }

        for(int is=0; is<nspin; is++)
          for(int m=0; m<m_tot; m++)
            this->local_sigma_new[ineq][is][count][m*m_tot+m] = 
              std::complex<double>(Gw_real[is][m], Gw_im[is][m])/GLC::Hartree_to_eV;

        count++;
        if(ifs.eof()) break;  //Check whether end of file is reached       
      }
      ifs.close();

      if(count!=this->n_omega)
      {
        std::cerr << "The number of frequency points of Sigma_omega.dat is not equal to nomega" << std::endl;
        std::exit(EXIT_FAILURE);
      }
    }

    return;
  }

  void self_energy_real_aixs::read_maxent_parms()
  {
    debug::codestamp("self_energy_real_aixs::read_maxent_parms");

    char word[100], word_low[100];

    std::ifstream ifs("maxent_params.dat", std::ios::in);
    if (!ifs)  
	  {
	  	std::cerr << "Error: fail to oepn maxent_params.dat!!!" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    ifs.seekg(0);   //set the position at the beginning of the file

    while(ifs.good())
    {
      ifs.getline(word, 100);   
      if(ifs.eof()) break;

      std::string line(word);

      if(!line.empty())
      {
        line.erase(0, line.find_first_not_of(" "));
        line.erase(line.find_last_not_of(" ") + 1);
      }

      if(!line.empty() && (line[0]!='}' || line[0]!='#'))
      {
        
        if(line.find_first_of('#') != std::string::npos) line.erase(line.find_first_of('#'));

        size_t pos=line.find("Nw");
        if(pos==std::string::npos) continue;

        line.erase(0,line.find_first_of(':')+1);
        line.erase(line.find_last_of(','));

        line.erase(0, line.find_first_not_of(" "));
        line.erase(line.find_last_not_of(" ") + 1); 

        this->n_omega = 2*atoi(line.c_str()) + 1;
        break;
      }
    }
    ifs.close();
    return;
  }

}
