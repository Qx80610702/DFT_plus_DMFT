#include "charge_scf.h"
#include "../debug.h"
#include "../mpi_environment.h"
#include "math_zone.h"

#include <omp.h>
#include <mpi.h>
#include <memory>

#include <mkl.h>

//test
// #include <fstream>
// #include <iomanip>
// #include <sstream>
// #include <string>
// #include <cmath>

namespace DFT_plus_DMFT
{
  void Charge_SCF::update_char_dens(
        const int axis_flag,
        const double mu,
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
        this->fik_DMFT[is][ik].resize(wb2ib[is].back()+1, 0.0);
      }
    }

    if(axis_flag == 0) //imaginary axis
      this->eva_fik_DMFT_imag_axis(
            mu, band,atom, proj, 
            sigma, in, space );
    else  //real axis
      this->eva_fik_DMFT_real_axis(
            mu, band, atom, proj, 
            sigma, in, space );

    this->eva_char_dens(in);


    return;
  }

  void Charge_SCF::eva_fik_DMFT_imag_axis(
        const double mu,
        DFT_output::KS_bands& band, 
        DFT_output::atoms_info& atom, 
        DFT_plus_DMFT::projector& proj,
        DMFT::self_energy& sigma,
        DMFT::input_info& in,
        DFT_plus_DMFT::Hilbert_space& space)
  {
    debug::codestamp("Charge_SCF::eva_fik_DMFT_imag_axis");

    const double beta = *(double*)in.parameter("beta");
    const int nomega = *(int*)in.parameter("n_omega");
    const int nks=band.nk();
    const int nspin=band.nspins();
    std::vector<double> k_weight = band.kweight();
    const int magnetism = *(int*)in.parameter("magnetism");
    const std::vector<int>& wbands=space.Wbands();
    const int n_valence = space.valence();
    const std::vector<double>& freq=sigma.sigma_imag.Matsubara_freq();
    const std::vector<std::vector<int>>& wb2ib = space.wbands2ibands();
    std::vector<std::vector<bool>>& corb_flag = space.correction_flag();
    const std::complex<double> im(0.0,1.0), one(1.0,0.0);

    std::vector<int> k_map;
    int task_nks=0;
    for(int ik=0; ik<nks; ik++){
      if(ik%mpi_ntasks() != mpi_rank()) continue;  //k_points are divided acording to the process id
      task_nks +=1;
      k_map.push_back(ik);
    }

    //Uncorrelated bands
    for(int is=0; is<nspin; is++)
      for(int ik=0; ik<task_nks; ik++)
        for(int iband=0; iband<=wb2ib[is].back(); iband++)
          if(!corb_flag[is][iband]) this->fik_DMFT[is][ik][iband] = 1.0;

    //Correlated bands
    for(int ik=0; ik<task_nks; ik++)
    {
      if(ik%mpi_ntasks() != mpi_rank()) continue;  //k_points are divided acording to process id

      sigma.evalute_lattice_sigma(
          0, magnetism, 
          nspin, wbands, atom, 
          proj.proj_access(ik) );

      const std::vector<std::vector<std::vector<std::complex<double>>>>&
            latt_sigma = sigma.lattice_sigma(0);
  
      for(int is=0; is<nspin; is++)
      {
        const auto& epsilon=space.eigen_val()[is][ k_map[ik] ];

        std::vector<std::vector<std::complex<double>>> KS_Gw(nomega);
        for(int iomega=0; iomega<nomega; iomega++)
          KS_Gw[iomega].resize(wbands[is]*wbands[is]);

        const int mkl_threads = mkl_get_max_threads();
        mkl_set_num_threads(1);  //set the number of threads of MKL library function to 1

        #pragma omp parallel
        {
          std::unique_ptr<int[]> ipiv(new int [wbands[is]]);
          int info_trf, info_tri;

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
            general_complex_matrix_inverse(&KS_Gw[iomega][0], wbands[is], &ipiv[0], info_trf, info_tri);
          }//iomega
        }
        mkl_set_num_threads(mkl_threads);

        for(int iband=0; iband<wbands[is]; iband++){
          double sum_tmp = 0.0;
          const double sigma_oo = latt_sigma[is][nomega-1][iband*wbands[is]+iband].real();

          for(int iomega=0; iomega<nomega; iomega++)
            sum_tmp += 2.0/beta*( KS_Gw[iomega][iband*wbands[is]+iband] -
                      one/(im*freq[iomega] + mu - epsilon[iband]-sigma_oo) ).real();

          sum_tmp += 1.0/(1.0+std::exp(beta*(epsilon[iband]+sigma_oo-mu)));

          this->fik_DMFT[is][ik][ wb2ib[is][iband] ] = sum_tmp;

        }//iband 
      }//is
    }//ik

    //TEST
    // std::ofstream ofs("fik.dat", std::ios::out);
    // for(int is=0; is<nspin; is++)
    //   for(int ik=0; ik<task_nks; ik++){
    //     for(int iband=0; iband<=wb2ib[is].back(); iband++)
    //       ofs << std::setw(6) << std::fixed << std::setprecision(3) << this->fik_DMFT[is][ik][iband];
    //     ofs << std::endl;
    //   }

    if(nspin==1 && !band.soc()) 
      for(int ik=0; ik<k_weight.size(); ik++)
        k_weight[ik] *= 2.0;

    for(int is=0; is<nspin; is++)
      for(int ik=0; ik<task_nks; ik++)
        for(int iband=0; iband<=wb2ib[is].back(); iband++)
          this->fik_DMFT[is][ik][iband] *= k_weight[ik];

    return;
  }

    void Charge_SCF::eva_fik_DMFT_real_axis(
        const double mu,
        DFT_output::KS_bands& band, 
        DFT_output::atoms_info& atom, 
        DFT_plus_DMFT::projector& proj,
        DMFT::self_energy& sigma,
        DMFT::input_info& in,
        DFT_plus_DMFT::Hilbert_space& space)
  {
    debug::codestamp("Charge_SCF::eva_fik_DMFT_real_axis");

    return;
  }

  void Charge_SCF::eva_char_dens(
        DMFT::input_info& in )
  {
    debug::codestamp("Charge_SCF::eva_char_dens");

    std::string file;
    std::vector<std::complex<double>> dens_cmplx;

    this->output_char_dense(
      *(int*)in.parameter("dft_solver"),
      file, dens_cmplx);

    return;
  }

  void Charge_SCF::output_char_dense(
      const int dft_solver,
      std::string file, 
      std::vector<std::complex<double>>& dens_cmplx)
  {
    debug::codestamp("Charge_SCF::output_char_dense");

    switch(dft_solver)
    {
      case 1: //aims
        #ifdef __FHIaims
        this->char_scf_aims.output_charge_density(file, dens_cmplx);
        #else
        std::cout << "FHI-aims has not been installed!!!  ";
        std::cout << "Suggestion:Install FHI-aims and then re-compile the codes." << std::endl;
        std::exit(EXIT_FAILURE);
        #endif   
        break;
      case 2: //ABACUS
        #ifdef __ABACUS
        // this->char_scf_aims.output_charge_density(file, dens_cmplx);
        std::cout << "Charge sel-consistent DMFT does not support ABACUS at present!!!  ";
        std::exit(EXIT_FAILURE);
        #else
        std::cout << "ABACUS has not been installed!!!  ";
        std::cout << "Suggestion:Install ABACUS and then re-compile the codes." << std::endl;
        std::exit(EXIT_FAILURE);
        #endif
        break;
      default:
        std::cout << "Not supported DFT_solver" << std::endl;
        std::exit(EXIT_FAILURE);
    }

    return;
  }
  
}