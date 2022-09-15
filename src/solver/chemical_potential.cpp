#include "chemical_potential.h"
#include "../mpi_environment.h"
#include "../debug.h"
#include "../constants.h"
#include "math_zone.h"
#include "../global_variables.h"
#include "../timer.h"

#include <mkl.h>
#include <mpi.h>
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
      std::cerr << "Error parameter of axis_flag" << std::endl;
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
    
    //==================================================
    // Calculate the eigen values of the Hamiltonian 
    // H_{ij}(\bfk,i\omega_{\infty}) = \epsilon_{i\bfk}\delta_{ij} + \bar{\Sigma}_{ij}(\mathbf{k},i\omega_{\infty})
    //==================================================
    if(this->epsilon_infty.empty()){
      this->epsilon_infty.resize(nspin);
      for(int is=0; is<nspin; is++){
        this->epsilon_infty[is].resize(task_nks);
        for(int ik_count=0; ik_count<task_nks; ik_count++)
          this->epsilon_infty[is][ik_count].resize(wbands[is], 0.0);
      }
    }

    for(int is=0; is<nspin; is++)
    {
      std::vector<std::complex<double>> Hij(wbands[is]*wbands[is],zero);

      int ik_count = 0;
      for(int ik=0; ik<band.nk(); ik++){
        if(ik%mpi_ntasks() != mpi_rank()) continue;  //k_points are divided acording to process id

        sigma.evalute_lattice_sigma_infty(
          0, magnetism, is, wbands, atom, 
          proj.proj_access(ik_count), Hij);

        const auto& epsilon = space.eigen_val()[is][ik];
        
        for(int iband=0; iband<wbands[is]; iband++)
          Hij[iband*wbands[is]+iband] += epsilon[iband];

        // lapack_int LAPACKE_zheev( int matrix_layout, char jobz, char uplo, 
        // lapack_int n, lapack_complex_double* a, lapack_int lda, double* w );
 
        lapack_int info_zheev = LAPACKE_zheev(LAPACK_ROW_MAJOR, 'N', 'U', wbands[is], 
                                              reinterpret_cast<MKL_Complex16*>(Hij.data()), 
                                              wbands[is], &this->epsilon_infty[is][ik_count][0] );
        if(info_zheev != 0){
          std::cerr << "Failed to compute the eigenvalues of the H_ij(\\bfk,i\\omega_{\\infty})" << std::endl;
          std::exit(EXIT_FAILURE);
        }
        ik_count++;
      }//ik
    }//is

    //==================================================
    // Update chemical potential
    //==================================================
    double max_U=0.0;
    for(int ineq=0; ineq<atom.inequ_atoms(); ineq++)
      max_U = max_U > atom.Uval(atom.ineq_iatom(ineq)) ? max_U : atom.Uval(atom.ineq_iatom(ineq));

    bool converged = false;
    double elw = window[0]-2.0*max_U, eup = window[1]+2.0*max_U;
    double mid_n_elec;
    
    for(int istep=0; istep<50; istep++) //50 may be a very safe value
    {
      this->sigma_corrected_mu = (elw + eup) / 2.0;
      
      mid_n_elec = this->evaluate_electrons_number_imag(
              band, space, sigma, atom, proj, beta, 
              magnetism, this->sigma_corrected_mu );

      if((mid_n_elec - n_valence)>=-1.0e-9 &&  (mid_n_elec - n_valence)<= 1.0e-9)
        converged = true;
      else if((mid_n_elec - n_valence) < -1.0e-9)
        elw = this->sigma_corrected_mu;
      else
        eup = this->sigma_corrected_mu;

// std::cerr << "mu        " << std::setw(14) << std::fixed << std::setprecision(9) << this->sigma_corrected_mu << '\n';
// std::cerr << "electrons " << std::setw(12) << std::fixed << std::setprecision(6) << mid_n_elec << '\n';
// std::cerr << "step " << istep << '\n';

      if(converged) break;
    }

    if(!converged){
      std::cerr << "Error in calculating chemical potential!!!" << std::endl;
      std::exit(EXIT_FAILURE);
    }

    GLV::ofs_running << "Chemical potential: " << std::setw(15) << std::fixed << std::setprecision(9) 
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
    const std::complex<double> im(0.0,1.0), one(1.0,0.0), zero(0.0,0.0);

    if(nspin==1 && !band.soc()) 
    {
      for(int ik=0; ik<weight.size(); ik++)
      {
        weight[ik] *= 2.0;
      }
    }

    double nele_sum=0.0;
    double sum_tmp=0.0;
    for(int is=0; is<nspin; is++)
    {
      std::vector<std::vector<std::complex<double>>> lattice_Sigma(nomega);
      for(int iomega=0; iomega<nomega; iomega++)
        lattice_Sigma[iomega].resize(wbands[is]*wbands[is], zero);

      std::vector<std::complex<double>> lattice_Sigma_infty(wbands[is]*wbands[is], zero);

      int ik_count = 0;
      for(int ik=0; ik<nks; ik++)
      {
        if(ik%mpi_ntasks() != mpi_rank()) continue;  //k_points are divided acording to process id

        sigma.evalute_lattice_sigma(
            0, mag, is, wbands, atom, 
            proj.proj_access(ik_count), 
            lattice_Sigma);

        sigma.evalute_lattice_sigma_infty(
            0, mag, is, wbands, atom, 
            proj.proj_access(ik_count), 
            lattice_Sigma_infty);

        const auto& epsilon = space.eigen_val()[is][ik];
        std::vector<std::vector<std::complex<double>>> KS_Gw(nomega);
        std::vector<std::vector<std::complex<double>>> KS_Gw_infty(nomega);
        for(int iomega=0; iomega<nomega; iomega++){
          KS_Gw[iomega].resize(wbands[is]*wbands[is]);
          KS_Gw_infty[iomega].resize(wbands[is]*wbands[is]);
        }

        int* ipiv = new int [wbands[is]];
        for(int iomega=0; iomega<nomega; iomega++)
        {
          for(int iband1=0; iband1<wbands[is]; iband1++)
          {
            for(int iband2=0; iband2<wbands[is]; iband2++)
            {           
              if(iband1==iband2){
                KS_Gw[iomega][iband1*wbands[is]+iband2] = im*freq[iomega] + mu 
                  -epsilon[iband1]-lattice_Sigma[iomega][iband1*wbands[is]+iband2];
                
                KS_Gw_infty[iomega][iband1*wbands[is]+iband2] = im*freq[iomega] + mu 
                  -epsilon[iband1]-lattice_Sigma_infty[iband1*wbands[is]+iband2];
              }
              else{
                KS_Gw[iomega][iband1*wbands[is]+iband2] = 
                  -lattice_Sigma[iomega][iband1*wbands[is]+iband2];
                
                KS_Gw_infty[iomega][iband1*wbands[is]+iband2] = 
                  -lattice_Sigma_infty[iband1*wbands[is]+iband2];
              }
            }
          }
          general_complex_matrix_inverse(&KS_Gw[iomega][0], wbands[is], &ipiv[0]);

          general_complex_matrix_inverse(&KS_Gw_infty[iomega][0], wbands[is], &ipiv[0]);
        }//iomega
        delete [] ipiv;

        for(int iband=0; iband<wbands[is]; iband++)
        {
          double fik_tmp = 0.0;

          for(int iomega=0; iomega<nomega; iomega++)
            fik_tmp += 2.0/beta*( KS_Gw[iomega][iband*wbands[is]+iband] -
                        KS_Gw_infty[iomega][iband*wbands[is]+iband] ).real();

          fik_tmp += 1.0/(1.0+std::exp(beta*(this->epsilon_infty[is][ik_count][iband]-mu)));

          this->fik_wind[is][ik_count][iband] = fik_tmp;
          
          sum_tmp += fik_tmp*weight[ik];
        }//iband

        ik_count++;
      }//ik
    }//is

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
    if (!ifs){
      std::cerr << "Error: fail to oepnDMFT_running.log!!!" << std::endl;
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
