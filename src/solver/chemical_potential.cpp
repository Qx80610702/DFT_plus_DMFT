#include "chemical_potential.h"
#include "../mpi_environment.h"
#include "../debug.h"
#include "../constants.h"
#include "math_zone.h"
#include "../global_variables.h"
#include "../timer.h"

#include <omp.h>
#include <mpi.h>
#include <memory>
#include <iostream>
#include <cmath>
#include <cstdlib> 
#include <iomanip>
#include <fstream>
#include <string>
#include <cstdlib>

namespace DFT_plus_DMFT
{
  void chemical_potential::update_chemical_potential(
        const int axis_flag,
        DFT_output::KS_bands& band, 
        DFT_output::atoms_info& atom, 
        DFT_plus_DMFT::projector& proj,
        DMFT::self_energy& sigma,
        DMFT::input_info& in,
        DFT_plus_DMFT::Hilbert_space& space )
  {
    debug::codestamp("chemical_potential::update_chemical_potential");

    double time;
    double seconds;
    timer::timestamp(time);

    switch(axis_flag)
    {
    case 0:
      this->evaluate_mu_bisection_imag(
          band, atom, proj, sigma, in, space);
      break;
    case 1:
      // ;
      break;
    default:
      GLV::ofs_error << "Error parameter of axis_flag" << std::endl;
      std::exit(EXIT_FAILURE);
    }

    timer::get_time(time, seconds);
    GLV::ofs_running << "Time consuming for evaluate chemical potention: " 
                << (int)seconds << "s" << std::endl;

    return;
  }

  void chemical_potential::evaluate_mu_bisection_imag( 
          DFT_output::KS_bands& band, 
          DFT_output::atoms_info& atom, 
          DFT_plus_DMFT::projector& proj, 
          DMFT::self_energy& sigma,
          DMFT::input_info& in,
          DFT_plus_DMFT::Hilbert_space& space)
  {
    debug::codestamp("chemical_potential::evaluate_mu_bisection");

    const double beta = *(double*)in.parameter("beta");
    const int nomega = *(int*)in.parameter("n_omega");
    const int nks=band.nk();
    const int nspin=band.nspins();
    const int magnetism = *(int*)in.parameter("magnetism");
    const std::vector<int>& wbands=space.Wbands();
    const int n_valence = space.valence();
    const std::vector<double>& window = space.Ener_window();
    const std::complex<double> zero(0.0,0.0);

    int task_nks=0;
    for(int ik=0; ik<band.nk(); ik++) 
    {
      if(ik%mpi_ntasks() != mpi_rank()) continue;  //k_points are divided acording to process id
      task_nks +=1;
    }

    if(this->fik_wind.empty()){
      this->fik_wind.resize(nspin);
      for(int is=0; is<nspin; is++){
        this->fik_wind[is].resize(task_nks);
        for(int ik=0; ik<task_nks; ik++){
          this->fik_wind[is][ik].resize(wbands[is], 0.0);
        }
      }
    }
    else{
      for(auto& iter1 : this->fik_wind)
        for(auto& iter2 : iter1)
          for(auto& iter3 : iter2)
            iter3 = 0.0;
    }

    double max_U=0.0;
    for(int ineq=0; ineq<atom.inequ_atoms(); ineq++)
      max_U = max_U > atom.Uval(atom.ineq_iatom(ineq)) ? max_U : atom.Uval(atom.ineq_iatom(ineq));

    bool converged=false;
    double elw = window[0], eup = window[1];
    double mid_n_elec;
    
    for(int istep=0; istep<50; istep++) //50 may be a very safe value
    {
      this->sigma_corrected_mu = (elw + eup) / 2.0;
      
      mid_n_elec=this->evaluate_electrons_number_imag(
              band, space, sigma, atom, proj, beta, 
              magnetism, this->sigma_corrected_mu );

      if((mid_n_elec - n_valence)>=-1.0e-9 &&  (mid_n_elec - n_valence)<= 1.0e-9)
        converged = true;
      else if((mid_n_elec - n_valence) < -1.0e-9)
        elw = this->sigma_corrected_mu;
      else
        eup = this->sigma_corrected_mu;

// GLV::ofs_error << "mu        " << std::setw(14) << std::fixed << std::setprecision(9) << this->sigma_corrected_mu << '\n';
// GLV::ofs_error << "electrons " << std::setw(12) << std::fixed << std::setprecision(6) << mid_n_elec << '\n';
// GLV::ofs_error << "step " << istep << '\n';

      if(converged) break;
    }

    if(!converged){
      GLV::ofs_error << "Error in calculating chemical potential!!!!" << std::endl;
      std::exit(EXIT_FAILURE);
    }

    GLV::ofs_running << "\nChemical potential: " << std::setw(15) << std::fixed << std::setprecision(9) 
                << this->sigma_corrected_mu*GLC::Hartree_to_eV << " eV" << std::endl;

    return;
  }

  double chemical_potential::evaluate_electrons_number_imag(
            DFT_output::KS_bands& band,
            DFT_plus_DMFT::Hilbert_space& space,
            DMFT::self_energy& sigma,
            DFT_output::atoms_info& atom,
            DFT_plus_DMFT::projector& proj,
            const double beta, const int mag, const double mu)
  {
    debug::codestamp("chemical_potential::evaluate_electrons_number_imag");

    const int nks = band.nk();
    const int nstates = band.nband();
    const int nspin = band.nspins();
    std::vector<double> weight = band.kweight();
    const std::vector<int>& wbands = space.Wbands();
    const std::vector<double>& freq=sigma.sigma_imag.Matsubara_freq();
    const int nomega=freq.size();
    const std::complex<double> im(0.0,1.0), one(1.0,0.0);

    if(nspin==1 && !band.soc()) 
    {
      for(int ik=0; ik<weight.size(); ik++)
      {
        weight[ik] *= 2.0;
      }
    }

    double nele_sum=0.0;
    double sum_tmp=0.0;

    int ik_count = 0;
    for(int ik=0; ik<nks; ik++)
    {
      if(ik%mpi_ntasks() != mpi_rank()) continue;  //k_points are divided acording to process id

      sigma.evalute_lattice_sigma(
          0, mag, nspin, wbands, atom, 
          proj.proj_access(ik_count) );

      const std::vector<std::vector<std::vector<std::complex<double>>>>&
            latt_sigma = sigma.lattice_sigma(0);

      for(int is=0; is<nspin; is++)
      {
        const auto& epsilon=space.eigen_val()[is][ik];

        // for(int iband=0; iband<wbands[is]; iband++)
        // {
        //   const double sigma_N = latt_sigma[ik_count][is][nomega-1][iband].real();

        //   for(int iomega=0; iomega<nomega; iomega++)
        //     sum_tmp += 2.0/beta*weight[ik]*( one/(im*freq[iomega]+mu-epsilon[iband]
        //               -latt_sigma[ik_count][is][iomega][iband]) -
        //               one/(im*freq[iomega]+mu-epsilon[iband]-sigma_N) ).real();

        //   sum_tmp += weight[ik]/(1.0+std::exp(beta*(epsilon[iband]+sigma_N-mu)));

        // }//iband

        std::vector<std::vector<std::complex<double>>> KS_Gw(nomega);
        for(int iomega=0; iomega<nomega; iomega++)
          KS_Gw[iomega].resize(wbands[is]*wbands[is]);

        const int mkl_threads = mkl_get_max_threads();
        mkl_set_num_threads(1);  //set the number of threads of MKL library function to 1
        #pragma omp parallel
        {
          std::unique_ptr<int[]> ipiv(new int [wbands[is]]);

          #pragma omp for
          for(int iomega=0; iomega<nomega; iomega++)
          {
            for(int iband1=0; iband1<wbands[is]; iband1++)
            {
              for(int iband2=0; iband2<wbands[is]; iband2++)
              {           
                if(iband1==iband2)
                  KS_Gw[iomega][iband1*wbands[is]+iband2] = im*freq[iomega] + mu 
                    -epsilon[iband1]-latt_sigma[is][iomega][iband1*wbands[is]+iband2];
                else
                  KS_Gw[iomega][iband1*wbands[is]+iband2] = 
                    -latt_sigma[is][iomega][iband1*wbands[is]+iband2];
              }
            }
            general_complex_matrix_inverse(&KS_Gw[iomega][0], wbands[is], &ipiv[0]);
          }//iomega
        }
        mkl_set_num_threads(mkl_threads);

        for(int iband=0; iband<wbands[is]; iband++)
        {
          double fik_tmp = 0.0;
          const double sigma_oo = latt_sigma[is][nomega-1][iband*wbands[is]+iband].real();

          for(int iomega=0; iomega<nomega; iomega++)
            fik_tmp += 2.0/beta*( KS_Gw[iomega][iband*wbands[is]+iband] -
                      one/(im*freq[iomega] + mu - epsilon[iband]-sigma_oo) ).real();

          fik_tmp += 1.0/(1.0+std::exp(beta*(epsilon[iband]+sigma_oo-mu)));

          this->fik_wind[is][ik_count][iband] = fik_tmp;

          sum_tmp += fik_tmp*weight[ik];
        }//iband
      }//is
      ik_count++;
    }//ik

    MPI_Allreduce(&sum_tmp, &nele_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    
    return nele_sum;
  }

  void chemical_potential::evaluate_mu_bisection_imag_DFT(
        DFT_output::KS_bands& band, 
        DFT_output::atoms_info& atom,
        DMFT::input_info& in,
        DFT_plus_DMFT::Hilbert_space& space)
  {
    debug::codestamp("chemical_potential::evaluate_mu_bisection_DFT");

    const double beta = *(double*)in.parameter("beta");
    const int nks=band.nk();
    const int nspin=band.nspins();
    const std::vector<int>& wbands=space.Wbands();
    const int n_valence = space.valence();
    const std::vector<double>& window = space.Ener_window();
    const std::complex<double> zero(0.0,0.0);

    bool converged=false;
    double elw = window[0], eup = window[1];
    for(int istep=0; istep<40; istep++)//40 may be a very safe value
    {
      this->DFT_mu = (elw + eup) / 2.0;
      
      double mid_n_elec=this->evaluate_electrons_number_imag_DFT(
                        band, space, beta, this->DFT_mu );

      if((mid_n_elec - n_valence)>=-1.0e-6 &&  (mid_n_elec - n_valence)<= 1.0e-6)
        converged = true;
      else if((mid_n_elec - n_valence) < -1.0e-6)
        elw = this->DFT_mu;
      else
        eup = this->DFT_mu;

      if(converged) break;
    }

    return;
  }

  double chemical_potential::evaluate_electrons_number_imag_DFT(
            DFT_output::KS_bands& band, 
            DFT_plus_DMFT::Hilbert_space& space,
            const double beta, const double mu)
  {
    debug::codestamp("chemical_potential::evaluate_electrons_number_imag_DFT");

    const int nks=band.nk();
    const int nspin=band.nspins();
    std::vector<double> weight = band.kweight();
    const std::vector<int>& wbands = space.Wbands();

    if(nspin==1 && !band.soc()) 
      for(int ik=0; ik<weight.size(); ik++)
        weight[ik] *= 2.0;

    double nele_sum=0.0;
    double sum_tmp=0.0;

    for(int is=0; is<nspin; is++)
    {
      int ik_count=-1;
      for(int ik=0; ik<nks; ik++)
      {
        if(ik%mpi_ntasks() != mpi_rank()) continue;
        ik_count++;

        const auto& epsilon=space.eigen_val()[is][ik];

        for(int iband=0; iband<wbands[is]; iband++)
          sum_tmp += weight[ik]/(1.0+std::exp(beta*(epsilon[iband]-mu)));

      }//ik
    }//is

    MPI_Allreduce(&sum_tmp, &nele_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    
    return nele_sum;
  }

  void chemical_potential::read_chemical_potential()
  {
    debug::codestamp("chemical_potential::read_chemical_potential");

    char word[100], word_low[100];

    std::ifstream ifs("DMFT_running.log", std::ios::in);
    if (!ifs)  {
      GLV::ofs_error << "Error: fail to oepnDMFT_running.log!!!" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    ifs.seekg(0);   //set the position at the beginning of the file

    while(ifs.good())
    {
      ifs.getline(word, 100);   
      if(ifs.eof()) break;
      
      std::string line = DMFT::input_info::strtolower(word);

      if(!line.empty())
      {
        line.erase(0, line.find_first_not_of(" "));
        line.erase(line.find_last_not_of(" ") + 1);
      }

      if(!line.empty())
      {
        size_t pos=line.find("chemical potential:");
        if(pos==std::string::npos) continue;

        line.erase(0,line.find_first_of(':')+1);
        line.erase(line.rfind("ev"));

        line.erase(0, line.find_first_not_of(" "));
        line.erase(line.find_last_not_of(" ") + 1);

        this->sigma_corrected_mu = atof(line.c_str())/GLC::Hartree_to_eV;
      }
    }
    ifs.close();

    return;
  }

}
