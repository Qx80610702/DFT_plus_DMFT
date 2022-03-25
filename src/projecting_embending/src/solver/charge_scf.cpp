#include "charge_scf.h"
#include "../debug.h"
#include "../mpi_environment.h"
#include "math_zone.h"
#include "../para/KS_eigenvectors.h"

#include <omp.h>
#include <mpi.h>
#include <memory>
#include <fstream>
#include <iomanip>

#include <mkl.h>

namespace DFT_plus_DMFT
{
  void Charge_SCF::update_char_dens(
      const int axis_flag,
      DFT_plus_DMFT::chemical_potential& Mu,
      DFT_output::KS_bands& band, 
      DFT_output::atoms_info& atom, 
      DFT_plus_DMFT::projector& proj,
      DMFT::self_energy& sigma,
      DMFT::input_info& in,
      DFT_plus_DMFT::Hilbert_space& space)
  {
    debug::codestamp("Charge_SCF::update_char_dens");

    this->eva_new_char_dens(
          axis_flag, Mu, band, atom, 
          proj, sigma, in, space );

    // this->output_char_dense(
    //       *(int*)in.parameter("dft_solver"), band.nk() );

    this->read_char_dense(
        *(int*)in.parameter("dft_solver"), band.nk() );

    // //TEST
    // double delta_tmp=0.0, delta=0.0;
    // for(int ik=0; ik<this->dens_mat.size(); ik++){
    //   for(int is=0; is<this->dens_mat[ik].size(); is++){
    //     for(int i=0; i<this->dens_mat[ik][is].size(); i++){
    //       delta_tmp += std::sqrt(std::norm( this->dens_mat[ik][is][i] 
    //             - this->dens_mat_last[ik][is][i] ));
    //     }
    //   }
    // }

    // MPI_Allreduce(&delta_tmp, &delta, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    // std::cout << "The difference between the input and output rho: " 
    //           << std::setw(15) << std::setprecision(9) << delta << std::endl;

    this->mix_char_dense(*(double*)in.parameter("charge_mix_beta"));

    return;
  }

  void Charge_SCF::eva_new_char_dens(
        const int axis_flag,
        DFT_plus_DMFT::chemical_potential& Mu,
        DFT_output::KS_bands& band, 
        DFT_output::atoms_info& atom, 
        DFT_plus_DMFT::projector& proj,
        DMFT::self_energy& sigma,
        DMFT::input_info& in,
        DFT_plus_DMFT::Hilbert_space& space)
  {
    debug::codestamp("Charge_SCF::eva_new_char_dens");

    const int nks=band.nk();
    const int nspin=band.nspins();
    std::vector<double> k_weight = band.kweight();
    const std::vector<std::vector<int>>& wb2ib = space.wbands2ibands();
    const std::vector<int>& wbands = space.Wbands();
    const int n_valence = space.valence();
    std::vector<std::vector<std::vector<double>>> fik_wind = Mu.fik();
    const int dft_solver = *(int*)in.parameter("dft_solver");

    std::vector<int> k_map;
    int task_nks=0;
    for(int ik=0; ik<nks; ik++){
      if(ik%mpi_ntasks() != mpi_rank()) continue;  //k_points are divided acording to the process id
      task_nks += 1;
      k_map.push_back(ik);
    }

    //Allocation
    this->dens_mat.resize(task_nks);
    for(int ik=0; ik<task_nks; ik++)
        this->dens_mat[ik].resize(nspin);

    this->fik_DMFT.resize(nspin);
    for(int is=0; is<nspin; is++){
      this->fik_DMFT[is].resize(task_nks);
      for(int ik=0; ik<task_nks; ik++){
        this->fik_DMFT[is][ik].resize(wb2ib[is].back()+1, 1.0);
      }
    }

    // this->eva_fik_DMFT_imag_axis(
    //   Mu.mu_corrected(), band,
    //   atom, proj, sigma, in, space);

    //fik
    for(int is=0; is<nspin; is++)
      for(int ik=0; ik<task_nks; ik++)
        for(int iband=0; iband<wbands[is]; iband++)
          this->fik_DMFT[is][ik][ wb2ib[is][iband] ] = fik_wind[is][ik][iband];

    //fik*k_weight
    for(int ik=0; ik<task_nks; ik++){
      for(int is=0; is<nspin; is++){
        for(int iband=0; iband<=wb2ib[is].back(); iband++){
          if(this->fik_DMFT[is][ik][iband] < -1.0e-2){
            std::cout << "Error in fik : negative occupation number!!!" << std::endl;
            std::cout << "spin:" << std::setw(2) << is
                      << "; ik:" << std::setw(4) << k_map[ik]
                      << "; iband:" << std::setw(4) << iband
                      << "; fik:" << std::setw(15) << std::fixed << std::setprecision(12) 
                      << this->fik_DMFT[is][ik][iband] << std::endl;
            std::exit(EXIT_FAILURE);
          }
          else if(std::fabs(this->fik_DMFT[is][ik][iband]) < 1.0e-2){
             this->fik_DMFT[is][ik][iband] = 1.0e-12;
          }
        }
      }
    }

    if(nspin==1 && !band.soc())
      for(int ik=0; ik<k_weight.size(); ik++)
        k_weight[ik] *= 2.0;

    double val_sum = 0.0, val_tmp=0.0;
    for(int ik=0; ik<task_nks; ik++){
      for(int is=0; is<nspin; is++){
        for(int iband=0; iband<wbands[is]; iband++){
          val_tmp += this->fik_DMFT[is][ik][ wb2ib[is][iband] ]*k_weight[ k_map[ik] ];
        }
      }
    }
    MPI_Allreduce(&val_tmp, &val_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    for(int ik=0; ik<task_nks; ik++){
      for(int is=0; is<nspin; is++){
        for(int iband=0; iband<wbands[is]; iband++){
          this->fik_DMFT[is][ik][ wb2ib[is][iband] ] *= (n_valence/val_sum);
        }
      }
    }

    //=======TEST=======
    std::ofstream ofs("fik-test.dat", std::ios::out);
    for(int is=0; is<nspin; is++)
      for(int ik=0; ik<task_nks; ik++){
        ofs << std::setw(4) << k_map[ik];
        for(int iband=0; iband<=wb2ib[is].back(); iband++)
          ofs << std::setw(15) << std::fixed << std::setprecision(12) << this->fik_DMFT[is][ik][iband];
        ofs << std::endl;
      }
    ofs.close();

    for(int is=0; is<nspin; is++)
      for(int ik=0; ik<task_nks; ik++)
        for(int iband=0; iband<=wb2ib[is].back(); iband++)
          this->fik_DMFT[is][ik][iband] *= k_weight[ k_map[ik] ]; 
    // Test
    double val_fik = 0.0, core_fik=0.0;
    double val_fik_tmp = 0.0, core_fik_tmp=0.0;
    for(int ik=0; ik<task_nks; ik++){
      for(int is=0; is<nspin; is++){
        for(int iband=0; iband<wbands[is]; iband++)
              val_fik_tmp += this->fik_DMFT[is][ik][ wb2ib[is][iband] ];

        for(int iband=0; iband<wb2ib[is][0]; iband++){
          core_fik_tmp += this->fik_DMFT[is][ik][iband];
        }
      }
    }
    MPI_Allreduce(&val_fik_tmp, &val_fik, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&core_fik_tmp, &core_fik, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    std::cout << "Sum of fik of valence bands: " << val_fik << std::endl;
    std::cout << "Sum of fik of core bands: " << core_fik << std::endl; 
    
    int max_threads=omp_get_max_threads();
    int threads_num = (max_threads>k_map.size()? k_map.size() : max_threads);
    const int mkl_threads = mkl_get_max_threads();
    mkl_set_num_threads(1);      //set the number of threads of MKL library function to 1
const std::complex<double> zero(0.0,0.0), one(1.0,0.0);
double Nele=0.0;
    #pragma omp parallel num_threads(threads_num)
    {
      std::vector<std::vector<std::complex<double>>> eigenvector;
      // std::vector<std::vector<std::complex<double>>> dens_cmplx;

      DFT_output::overlap_matrix ovlp(dft_solver, "dft/outputs_to_DMFT/overlap_matrix");
      const int n_basis = ovlp.ovlp_aims.nbasis();
      std::vector<std::complex<double>> mat_tmp(n_basis*n_basis);

      #pragma omp for
      for(int ik=0; ik<task_nks; ik++){
        this->eva_k_densmat( 
          dft_solver, nspin, 
          ik, k_map[ik], space, 
          eigenvector, this->dens_mat[ik] );

        //========Test==========
        ovlp.evaluate_ovlp_k(k_map[ik], atom);        //caculate overlap matrix in k-space
        const std::vector<std::complex<double>>& ovlp_mat = ovlp.ovlp_aims.ovlp_mat_work();

        for(int is=0; is<nspin; is++){
          double tmp=0.0;
          cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                      n_basis, n_basis, n_basis,
                      &one,
                      &this->dens_mat[ik][is][0], n_basis,
                      &ovlp_mat[0], n_basis,
                      &zero,
                      &mat_tmp[0], n_basis);

          for(int ibasis=0; ibasis<n_basis; ibasis++)
            tmp += mat_tmp[ibasis*n_basis+ibasis].real();
          
          #pragma omp atomic
            Nele += tmp;
        } 

      }
    }
    mkl_set_num_threads(mkl_threads);

    double ntmp;
    MPI_Allreduce(&Nele, &ntmp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    std::cout << "Number of total electrons : " << ntmp << std::endl;

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
    for(int ik=0; ik<task_nks; ik++){
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
    //     ofs << std::setw(4) << k_map[ik];
    //     for(int iband=0; iband<=wb2ib[is].back(); iband++)
    //       ofs << std::setw(15) << std::fixed << std::setprecision(12) << this->fik_DMFT[is][ik][iband];
    //     ofs << std::endl;
    //   }
    // ofs.close();

    // this->fik_test.resize(nspin);
    // for(int is=0; is<nspin; is++){
    //   this->fik_test[is].resize(nks);
    //   for(int ik=0; ik<nks; ik++)
    //     this->fik_test[is][ik].resize(wb2ib[is].back()+1, 0.0);
    // }

    // int i_k_point, i_state;

    // std::ifstream ifs("dft/occu.dat", std::ios::in);

    // if (!ifs){
	  // 	std::cout << "Fail to oepn dft/occu.dat" << std::endl;
    //   std::exit(EXIT_FAILURE);
    // }
    // ifs.seekg(0);      //set the position at the beginning of the file

    // while(ifs.good()){
    //   ifs >> i_k_point;
    //   if(ifs.eof()) break;

    //   for(int i_state=0; i_state<wb2ib[0].back()+1; i_state++)
    //     ifs >> this->fik_test[0][i_k_point][i_state];

    //   ifs.ignore(150, '\n'); 

    //   if(ifs.eof()) break;  //Check whether end of file is reached 
    // }
    // ifs.close();

    // //TEST
    // std::ofstream ofs("fik.dat", std::ios::out);
    // for(int is=0; is<nspin; is++)
    //   for(int ik=0; ik<nks; ik++){
    //     ofs << std::setw(4) << ik;
    //     for(int iband=0; iband<=wb2ib[is].back(); iband++)
    //       ofs << std::setw(15) << std::fixed << std::setprecision(12) << this->fik_test[is][ik][iband];
    //     ofs << std::endl;
    //   }

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

  void Charge_SCF::eva_k_densmat(
        const int dft_solver,
        const int nspin,
        const int ik,
        const int i_k_point,
        DFT_plus_DMFT::Hilbert_space& space,
        std::vector<std::vector<
        std::complex<double>>>& eigenvector,
        std::vector<std::vector<
        std::complex<double>>>& dense_cmplx )
  {
    debug::codestamp("Charge_SCF::eva_k_densmat");

    wave_function wfc;

    wfc.read_DMFT_occ_subset(
        "dft/outputs_to_DMFT/KS_eigenvector/", 
        i_k_point, space, eigenvector );       //read KS-eigenvector

    const int nbasis = wfc.basis_n();

    for(int ispin=0; ispin<nspin; ispin++)
      for(int iband=0; iband<this->fik_DMFT[ispin][ik].size(); iband++)
        for(int ibasis=0; ibasis<nbasis; ibasis++){
          if(this->fik_DMFT[ispin][ik][iband] < -1.0e-6){
            std::cout << "Error in fik : negative occupation number!!!" << std::endl;
            std::exit(EXIT_FAILURE);
          }
          else if(std::fabs(this->fik_DMFT[ispin][ik][iband]) < 1.0e-6){
             eigenvector[ispin][ibasis*this->fik_DMFT[ispin][ik].size()+iband] *= 0.0;
          }
          else{
            eigenvector[ispin][ibasis*this->fik_DMFT[ispin][ik].size()+iband] *= std::sqrt(this->fik_DMFT[ispin][ik][iband]);
          }
        }
    
    // void cblas_zherk (const CBLAS_LAYOUT Layout, const CBLAS_UPLO uplo, const
    //    CBLAS_TRANSPOSE trans, const MKL_INT n, const MKL_INT k, const double alpha, const void
    //    *a, const MKL_INT lda, const double beta, void *c, const MKL_INT ldc);

    for(int ispin=0; ispin<nspin; ispin++){
      if(dense_cmplx[ispin].empty())
        dense_cmplx[ispin].resize(nbasis*nbasis);

      cblas_zherk(CblasRowMajor, CblasUpper, CblasNoTrans, 
                  nbasis, this->fik_DMFT[ispin][ik].size(), 
                  1.0, &eigenvector[ispin][0], this->fik_DMFT[ispin][ik].size(), 
                  0.0, &dense_cmplx[ispin][0], nbasis );

      for(int ibasis1=0; ibasis1<nbasis; ibasis1++)
        for(int ibasis2=0; ibasis2<ibasis1; ibasis2++)
          dense_cmplx[ispin][ibasis1*nbasis + ibasis2] = std::conj(dense_cmplx[ispin][ibasis2*nbasis + ibasis1]);
    }

    return;
  }

  void Charge_SCF::read_char_dense(
      const int dft_solver,
      const int nks )
  {
    debug::codestamp("Charge_SCF::read_char_dense");

    const std::complex<double> zero(0.0,0.0);
    
    //Allocation
    this->dens_mat_last.resize(this->dens_mat.size());
    for(int ik=0; ik<this->dens_mat.size(); ik++){
      this->dens_mat_last[ik].resize( this->dens_mat[ik].size() );
      for(int is=0; is<this->dens_mat[ik].size(); is++){
        this->dens_mat_last[ik][is].resize( this->dens_mat[ik][is].size(), zero);
      }
    }

    switch(dft_solver)
    {
      case 1: //aims
        #ifdef __FHIaims
        this->char_scf_aims.read_charge_density(nks, this->dens_mat_last);
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

  void Charge_SCF::mix_char_dense(
    const double mix_beta )
  {
    debug::codestamp("Charge_SCF::mix_char_dense");

    for(int ik=0; ik<this->dens_mat.size(); ik++){
      for(int is=0; is<this->dens_mat[ik].size(); is++){
        for(int i=0; i<this->dens_mat[ik][is].size(); i++){
          this->dens_mat[ik][is][i] = 
            (1.0-mix_beta)*this->dens_mat[ik][is][i] +
            mix_beta*this->dens_mat_last[ik][is][i];
        }
      }
    }

    return;
  }

  void Charge_SCF::output_char_dense(
      const int dft_solver,
      const int nks )
  {
    debug::codestamp("Charge_SCF::output_char_dense");

    switch(dft_solver)
    {
      case 1: //aims
        #ifdef __FHIaims
        this->char_scf_aims.output_charge_density(nks, this->dens_mat);
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

  void Charge_SCF::prepare_nscf_dft(
      const int dft_solver, 
      const int max_DFT_step )
  {
    debug::codestamp("Charge_SCF::prepare_nscf_dft");

    switch(dft_solver)
    {
      case 1: //aims
        #ifdef __FHIaims
        this->char_scf_aims.prepare_nscf_dft(max_DFT_step);
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