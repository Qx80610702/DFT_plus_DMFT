#ifdef __FHIaims

#include "charge_scf_aims.h"
#include "../debug.h"
#include "../mpi_environment.h"
#include "../para/input.h"
#include "../global_variables.h"
#include "math_zone.h"

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
  void Charge_SCF_aims::output_charge_density_matrix(
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
        const bool initial_charge,
        const bool DMFT_charge )
  {
    debug::codestamp("Charge_SCF_aims::read_charge_density");

    if(initial_charge && DMFT_charge){
      GLV::ofs_error << "There some logical errors in read the charge density" <<std::endl;
      std::exit(EXIT_FAILURE);
    }
    
    std::stringstream rho_file, part_file;

    part_file << "dft/outputs_to_DMFT/charge_density/partition_tab" 
              << std::setfill('0') << std::setw(6) 
              << mpi_rank() << ".dat";

    if(DMFT_charge)
      rho_file << "dft/outputs_to_DMFT/charge_density/rho" 
               << std::setfill('0') << std::setw(6) 
               << mpi_rank() << ".read";
    else
      rho_file << "dft/outputs_to_DMFT/charge_density/rho" 
               << std::setfill('0') << std::setw(6) 
               << mpi_rank() << ".out";

    //partition table
    if(this->partition_tab.empty()){
      std::ifstream ifsp(part_file.str().c_str(), std::ios::in);
      if (!ifsp){
      	GLV::ofs_error << "Fail to oepn " << part_file.str().c_str() << std::endl;
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
    	GLV::ofs_error << "Fail to oepn " << rho_file.str().c_str() << std::endl;
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

    if(DMFT_charge){
      if(this->rho_out.empty()) this->rho_out = rho_tmp;
      else{
        for(int ispin=0; ispin<rho_tmp.size(); ispin++)
          for(int igrid=0; igrid<rho_tmp[ispin].size(); igrid++)
            this->rho_out[ispin][igrid] = rho_tmp[ispin][igrid];
      }
    }
    else{
      if(initial_charge){
        if(!this->Opt_rho.empty()) this->Opt_rho.clear();

        this->Opt_rho.push_back(rho_tmp);
      }
      else{
        if(this->rho_out.empty()) this->rho_out = rho_tmp;
        else{
          for(int ispin=0; ispin<rho_tmp.size(); ispin++)
            for(int igrid=0; igrid<rho_tmp[ispin].size(); igrid++)
              this->rho_out[ispin][igrid] = rho_tmp[ispin][igrid];
        }
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
    c_elsi_init_rw(&rwh, 0, 0, nbasis, 0.0);

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

  void Charge_SCF_aims::update_data(
        const int mix_step,
        const int max_mixing_step)
  {
    debug::codestamp("Charge_SCF_aims::update_data");

    //Residual vector
    std::vector<std::vector<double>> vector_tmp;
    vector_tmp = this->Opt_rho.back();

    for(int is=0; is<vector_tmp.size(); is++)
        for(int igrid=0; igrid<vector_tmp[is].size(); igrid++)
          vector_tmp[is][igrid] = this->rho_out[is][igrid] - vector_tmp[is][igrid];

    if(this->Rrho.size()<max_mixing_step) 
      this->Rrho.push_back(vector_tmp);   //Rrho.size()<=max_mixing_step
    else{
      this->Rrho.pop_front();
      this->Rrho.push_back(vector_tmp);
    }

    //dRrho
    if(mix_step>1){
      for(int is=0; is<this->Rrho[0].size(); is++)
        for(int igrid=0; igrid<this->Rrho[0][is].size(); igrid++)
          vector_tmp[is][igrid] = this->Rrho.back()[is][igrid] - 
                        this->Rrho[this->Rrho.size()-2][is][igrid];

      if(this->dRrho.size() < max_mixing_step-1) 
        this->dRrho.push_back(vector_tmp); //dRrho.size()<=7
      else{
        this->dRrho.pop_front();
        this->dRrho.push_back(vector_tmp);
      }
    }

    return;
  }

  void Charge_SCF_aims::update_alpha(
      const int mix_step,
      std::vector<double>& alpha)
  {
    debug::codestamp("Charge_SCF_aims::update_alpha");

    if(mix_step==1) return;

    //Calculate Abar
    std::vector<double> Abar(this->dRrho.size()*this->dRrho.size(), 0.0);

    for(int is=0; is<this->dRrho.back().size();is++)
      for(int i=0; i<this->dRrho.size(); i++)
        for(int j=0; j<this->dRrho.size(); j++){
          for(int igrid=0; igrid<this->dRrho[i][is].size(); igrid++)
            Abar[i*this->dRrho.size()+j] += std::pow(this->partition_tab[igrid],2)*
              this->dRrho[i][is][igrid]*this->dRrho[j][is][igrid];
        }

    std::vector<int> ipiv(this->dRrho.size());
    general_real_matrix_inverse(&Abar[0], this->dRrho.size(), &ipiv[0]);

    //dRR
    std::vector<double> dRR(this->dRrho.size(), 0.0);
    for(int is=0; is<this->dRrho.back().size();is++)
      for(int i=0; i<this->dRrho.size(); i++)
        for(int igrid=0; igrid<this->dRrho[i][is].size(); igrid++)
          dRR[i] += std::pow(this->partition_tab[igrid],2)*
           this->dRrho[i][is][igrid]*this->Rrho.back()[is][igrid];

    //alpha
    alpha.resize(this->dRrho.size(), 0.0);
    for(int i=0; i<this->dRrho.size(); i++)
      for(int j=0; j<this->dRrho.size(); j++)
        alpha[i] -= Abar[j*this->dRrho.size()+i]*dRR[j];

    return;
  }

  void Charge_SCF_aims::mixing_density(
    const int mix_step,
    const double mixing_beta,
    const int max_mixing_step,
    std::vector<double>& alpha,
    double& charge_change)
  {
    debug::codestamp("Charge_SCF_aims::update_density");

    //change of charge density
    charge_change = 0.0;
    for(int is=0; is<this->Rrho[0].size(); is++)
      for(int igrid=0; igrid<this->Rrho[0][is].size(); igrid++)
        charge_change += this->partition_tab[igrid]*
            std::pow(this->Rrho.back()[is][igrid], 2);

    charge_change = std::sqrt(charge_change);

    if(mix_step==1){//plain mixing
      std::vector<std::vector<double>> rho_tmp;
      rho_tmp = this->Opt_rho.back();

      for(int is=0; is<this->Opt_rho.back().size(); is++)
        for(int igrid=0; igrid<this->Opt_rho.back()[is].size(); igrid++){
          rho_tmp[is][igrid] = this->rho_out[is][igrid]*mixing_beta +
              (1.0-mixing_beta)*rho_tmp[is][igrid];
        }

      this->Opt_rho.push_back(rho_tmp);
    }
    else{//Pulay mixing
      std::vector<std::vector<double>> rho_tmp;

      //First part: \rho^{opt}
      rho_tmp = this->Opt_rho.back();

      for(int istep=0; istep<alpha.size(); istep++){
        for(int is=0; is<this->Opt_rho[istep].size(); is++)
          for(int igrid=0; igrid<this->Opt_rho[istep][is].size(); igrid++)
            rho_tmp[is][igrid] += alpha[istep]*
              ( this->Opt_rho[istep+1][is][igrid] - 
              this->Opt_rho[istep][is][igrid] );
      }

      //Second part: \Rrho^{opt}
      for(int is=0; is<this->Rrho.back().size(); is++)
          for(int igrid=0; igrid<this->Rrho.back()[is].size(); igrid++)
            rho_tmp[is][igrid] += mixing_beta*this->Rrho.back()[is][igrid];

      for(int istep=0; istep<alpha.size(); istep++){
        for(int is=0; is<this->Rrho[istep].size(); is++)
          for(int igrid=0; igrid<this->Rrho[istep][is].size(); igrid++)
            rho_tmp[is][igrid] += mixing_beta*alpha[istep]*
              ( this->Rrho[istep+1][is][igrid] - 
              this->Rrho[istep][is][igrid] );
      }

      if(this->Opt_rho.size()<max_mixing_step) 
        this->Opt_rho.push_back(rho_tmp);   //this->Opt_rho.seize()<=max_mixing_step
      else{
        this->Opt_rho.pop_front();
        this->Opt_rho.push_back(rho_tmp);
      }
    }

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
      GLV::ofs_error << "Fail to oepn file control.in" << std::endl;
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
    std::ifstream if_test("dft/control.in-original", std::ios::in);
    if(!if_test){ //file dft/control.in-original not exist
      std::ofstream ofs_copy("dft/control.in-original", std::ios::out);
      for(int i=0; i<lines.size(); i++){
        std::string line_str = lines[i];
        ofs_copy << line_str << std::endl;
      }
      ofs_copy.close();

      //============Renew control.in==================
      std::ofstream ofs("dft/control.in", std::ios::out);

      for(int i=0; i<lines.size(); i++){
        std::string line_str = lines[i];

        if( line_str.empty() ) ofs << std::endl;
        else{
          std::string line = DMFT::input_info::strtolower(&line_str[0]);

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
    }

    return;
  }

}

#endif
