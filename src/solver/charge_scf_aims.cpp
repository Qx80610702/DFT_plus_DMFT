#ifdef __FHIaims

#include "charge_scf_aims.h"
#include "../debug.h"
#include "../mpi_environment.h"
#include "../para/input.h"
#include "../global_variables.h"

#include <cmath>
#include <sstream>
#include <string>
#include <iomanip>
#include <string>
#include <cstring>
#include <fstream>
#include <cstdlib>   //Use exit function

#include <elsi.h>

namespace DFT_plus_DMFT
{
  void Charge_SCF_aims::output_charge_density(
        const int nks, 
        std::vector<std::vector<std::vector<
        std::complex<double>>>>& dens_mat_cmplx)
  {
    debug::codestamp("Charge_SCF_aims::output_charge_density");
    const int nbasis = (int)std::sqrt(dens_mat_cmplx[0][0].size());

    std::vector<int> k_map;
    int task_nks=0;
    for(int ik=0; ik<nks; ik++){
      if(ik%mpi_ntasks() != mpi_rank()) continue;  //k_points are divided acording to the process id
      task_nks +=1;
      k_map.push_back(ik);
    }

    elsi_rw_handle rwh;
    c_elsi_init_rw(&rwh,1,0,nbasis,0.0);

    std::vector<std::complex<double>> dense_mat_tmp(nbasis*nbasis);

    for(int ik=0; ik<task_nks; ik++){
      for(int ispin=0; ispin<dens_mat_cmplx[0].size(); ispin++){
        std::stringstream ss;
        ss << "dft/D_spin_" 
           << std::setfill('0') << std::setw(2) << ispin+1
           << "_kpt_"
           << std::setfill('0') << std::setw(6) << k_map[ik]+1
           << ".csc";

        std::string file = ss.str();
        
        // std::stringstream rmfl;
        // rmfl << "test -f " << file << " && rm " << file;
        // system(rmfl.str().c_str());

        //Fortran order : column major
        for(int ibasis1=0; ibasis1<nbasis; ibasis1++)
          for(int ibasis2=0; ibasis2<nbasis; ibasis2++)
            dense_mat_tmp[ibasis1*nbasis + ibasis2] = 
              std::conj(dens_mat_cmplx[ik][ispin][ibasis1*nbasis + ibasis2]);

        c_elsi_write_mat_complex(rwh, &file[0], reinterpret_cast<double _Complex*>(dense_mat_tmp.data()));
        
      }
    }

    c_elsi_finalize_rw(rwh);

    return;
  }

    void Charge_SCF_aims::read_charge_density(
        const int istep, const bool DMFT_charge)
  {
    debug::codestamp("Charge_SCF_aims::read_charge_density");
    
    std::stringstream rho_file, part_file;

    part_file << "dft/outputs_to_DMFT/charge_density/partition_tab" << mpi_rank() << ".dat";

    if(DMFT_charge)
      rho_file << "dft/outputs_to_DMFT/charge_density/rho" << mpi_rank() << ".read";
    else
      rho_file << "dft/outputs_to_DMFT/charge_density/rho" << mpi_rank() << ".out";

    //partition table
    if(this->partition_tab.empty()){
      std::ifstream ifsp(part_file.str().c_str(), std::ios::in);
      if (!ifsp){
      	GlobalV::ofs_error << "Fail to oepn " << part_file.str().c_str() << std::endl;
        std::exit(EXIT_FAILURE);
      }
      ifsp.seekg(0); //set the position at the beginning of the file

      double value;
      while(ifsp.good())
      {
        ifsp >> value;
        if(ifsp.eof()) break;
        this->partition_tab.push_back(value);
        ifsp.ignore(150,'\n');

        if(ifsp.eof()) break;   //Check whether end of file is reached 
      }
      ifsp.close();
    }

    //charge density rho
    int nspin_tmp, ngrids_tmp, ispin;
    double value_tmp;
    std::vector<std::vector<double>> rho_tmp;
    std::ifstream ifsr(rho_file.str().c_str(), std::ios::in);
    if (!ifsr){
    	GlobalV::ofs_error << "Fail to oepn " << rho_file.str().c_str() << std::endl;
      std::exit(EXIT_FAILURE);
    }
    ifsr.seekg(0); //set the position at the beginning of the file
    
    ifsr >> nspin_tmp;
    ifsr >> ngrids_tmp;
    rho_tmp.resize(nspin_tmp);

    while(ifsr.good())
    {
      ifsr >> ispin;
      if(ifsr.eof()) break;
      ifsr >> value_tmp;
      rho_tmp[ispin].push_back(value_tmp);
      ifsr.ignore(150,'\n');

      if(ifsr.eof()) break;   //Check whether end of file is reached 
    }
    ifsr.close();

    if(istep<=8)
    {
      this->rho.push_back(rho_tmp);
      
      if(istep==1) this->istep2index.push_back(1);
      else{
        for(auto& iter : this->istep2index) iter++;
        this->istep2index.push_back(1);
      }
    }
    else{
      for(int index = 0; index<8; index++){
        if(this->istep2index[index]==8){
          this->istep2index[index] = 1;

          for(int ispin=0; ispin<rho_tmp.size(); ispin++)
            for(int igrid=0; igrid<rho_tmp[ispin].size(); igrid++)
              this->rho[index][ispin][igrid] = rho_tmp[ispin][igrid];
        }
        else this->istep2index[index] += 1;
      }
    }

    return;
  }

  void Charge_SCF_aims::read_charge_density_matrix(
        const int nks, 
        std::vector<std::vector<std::vector<
        std::complex<double>>>>& dens_mat_cmplx)
  {
    debug::codestamp("Charge_SCF_aims::read_charge_density_matrix");
    const int nbasis = (int)std::sqrt(dens_mat_cmplx[0][0].size());

    std::vector<int> k_map;
    int task_nks=0;
    for(int ik=0; ik<nks; ik++){
      if(ik%mpi_ntasks() != mpi_rank()) continue;  //k_points are divided acording to the process id
      task_nks +=1;
      k_map.push_back(ik);
    }

    elsi_rw_handle rwh;
    c_elsi_init_rw(&rwh,0,0,nbasis,0.0);

    for(int ik=0; ik<task_nks; ik++){
      for(int ispin=0; ispin<dens_mat_cmplx[0].size(); ispin++){
        std::stringstream ss;
        ss << "dft/D_spin_" 
           << std::setfill('0') << std::setw(2) << ispin+1
           << "_kpt_"
           << std::setfill('0') << std::setw(6) << k_map[ik]+1
           << ".csc";

        std::string file = ss.str();
        
        // std::stringstream rmfl;
        // rmfl << "test -f " << file << " && rm " << file;
        // system(rmfl.str().c_str());

        c_elsi_read_mat_complex(rwh, &file[0], reinterpret_cast<double _Complex*>(dens_mat_cmplx[ik][ispin].data()));

        //C++ order : row major
        for(int ibasis1=0; ibasis1<nbasis; ibasis1++)
          for(int ibasis2=0; ibasis2<nbasis; ibasis2++)
            dens_mat_cmplx[ik][ispin][ibasis1*nbasis + ibasis2] = 
              std::conj(dens_mat_cmplx[ik][ispin][ibasis1*nbasis + ibasis2]);
      }
    }

    c_elsi_finalize_rw(rwh);

    return;
  }

  void Charge_SCF_aims::prepare_nscf_dft()
  {
    debug::codestamp("Charge_SCF_aims::prepare_nscf_dft");

    if(mpi_rank() != 0) return; 

    //==================================
    //      Update control.in
    //==================================
    std::vector<std::string> lines;
    char word[200];

    //=========Read control.in==========
    std::ifstream ifs("dft/control.in", std::ios::in);
    if(!ifs){
      GlobalV::ofs_error << "Fail to oepn file control.in" << std::endl;
      std::exit(EXIT_FAILURE);
    }

    ifs.seekg(0);   //set the position at the beginning of the file

    while(ifs.good())
    {
      ifs.getline(word, 200);   
      if(ifs.eof()) break;

      std::string line(word);
      lines.push_back(line);

      if(ifs.eof()) break;
    }
    ifs.close();

    //============Cpoy control.in as control.in-original==================
    std::ofstream ofs_copy("dft/control.in-original", std::ios::out);
    for(int i=0; i<lines.size(); i++){
      std::string line_str = lines[i];
      ofs_copy << line_str << std::endl;
    }
    ofs_copy.close();


    //============Renew control.in==================
    std::ofstream ofs("dft/control.in", std::ios::out);
    char word_low[200];

    for(int i=0; i<lines.size(); i++){
      std::string line_str = lines[i];

      if( line_str.empty() ) ofs << std::endl;
      else{
        DMFT::input_info::strtolower(&line_str[0], word_low);
        std::string line(word_low);

        line.erase(0, line.find_first_not_of(" "));
        line.erase(line.find_last_not_of(" ") + 1);

        if(!line.empty() && line[0]!='#'){
          std::string param = line.substr( 0, line.find_first_of(' ') );
          
          if(std::strcmp(param.c_str(), "dft_plus_dmft")==0){
            ofs << line_str << std::endl;
            ofs << "  mixer               linear" << std::endl;
            ofs << "  charge_mix_param     1.0" << std::endl;
            ofs << "  sc_iter_limit         0"  << std::endl;
            ofs << "  elsi_restart read 1" << std::endl;
          }
          else{
            if(std::strcmp(param.c_str(), "charge_mix_param")!=0 &&
              std::strcmp(param.c_str(), "elsi_restart") !=0 &&
              std::strcmp(param.c_str(), "sc_iter_limit") !=0 &&
              std::strcmp(param.c_str(), "mixer") !=0 )
              ofs << line_str << std::endl;
          } 
        }
        else{
          ofs << line_str << std::endl;
        }
      }
    }
    ofs.close();

    return;
  }

}

#endif
