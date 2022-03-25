#include "Anderson_impurity.h"
#include "../mpi_environment.h"
#include "math_zone.h"
#include "../debug.h"
#include "../timer.h"
#include "../constants.h"

#include <mpi.h>
#include <omp.h>
#include <memory>
#include <cmath>
#include <iostream>
#include <cstdlib>
#include <iomanip>

//test
// #include <fstream>
// #include <sstream>
namespace DMFT
{
  void impurity::evaluate_impurity_level(
        DFT_plus_DMFT::Hilbert_space& space,
        DFT_output::KS_bands& band, 
        DFT_plus_DMFT::projector& proj,
        DFT_output::atoms_info& atom)
  {
    debug::codestamp("impurity::evaluate_impurity_level");

    const int nspin = band.nspins();
    const int nks = band.nk();
    const std::vector<int>& wbands = space.Wbands();
    const auto& epsilon = space.eigen_val();
    const int nprocs = mpi_ntasks();
    const int myid = mpi_rank();
    const int ineq_num = atom.inequ_atoms();
    const std::vector<int>& norb_sub = atom.iatom_norb();
    const std::vector<std::vector<int>>& orb_index = atom.Im2iorb();
    const int norb = atom.norb();
    const std::vector<double>& fk = band.kweight();
    const auto& TB_Hk = proj.Hloc();

    const std::complex<double> zero(0.0,0.0);

    //Allocation
    if(this->impurity_level.empty())
    {
      this->impurity_level.resize(ineq_num);
      for(int ineq=0; ineq<ineq_num; ineq++)
      {
        const int iatom = atom.ineq_iatom(ineq);
        const int m_tot=norb_sub[iatom];
  
        this->impurity_level[ineq].resize(nspin);
        for(int is=0; is<nspin; is++)
          this->impurity_level[ineq][is].resize(m_tot*m_tot, zero);
      }
    }

    //Evaluation
    for(int ineq=0; ineq<ineq_num; ineq++)
    {
      const int iatom = atom.ineq_iatom(ineq);
      const int m_tot=norb_sub[iatom];

      std::unique_ptr<std::complex<double>[]> tmp(new std::complex<double> [m_tot*m_tot]);

      for(int is=0; is<nspin; is++)
      {
        for(int m_index=0; m_index<m_tot*m_tot; m_index++) tmp[m_index] = zero;

        const std::vector<std::complex<double>>& dc_term = this->sigma.dc.Vdc()[ineq][is];

        int ik_count=-1;
        for(int ik=0; ik<nks; ik++)
        {
          if(ik%nprocs != myid) continue;
          ik_count++;

          for(int m=0; m<m_tot; m++)
            tmp[m*m_tot+m] += fk[ik]*TB_Hk[ik_count][iatom][is][m];                                           
        }//ik

        MPI_Allreduce(&tmp[0], &this->impurity_level[ineq][is][0], m_tot*m_tot,
                      MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD);

        for(int m=0; m<m_tot; m++)
          this->impurity_level[ineq][is][m*m_tot+m] -= dc_term[m*m_tot+m];

        //=============================================
        //   Local symmetry operation
        //=============================================
        const int local_symmetry = atom.local_sym();
        const int corr_L=atom.L(ineq);
        DFT_output::atoms_info::symmetry_operation_matrix<std::complex<double>>(
              local_symmetry, corr_L, m_tot, &this->impurity_level[ineq][is][0]);

      }//is
    }//ineq

    if(mpi_rank()==0)
    {
      std::cout << "\n================Impurity level(eV)===============\n";
      for(int ineq=0; ineq<atom.inequ_atoms(); ineq++)
      {
        const int iatom=atom.ineq_iatom(ineq);
        const int m_tot=norb_sub[iatom];

        for(int is=0; is<nspin; is++)
        {
          std::cout << "===========impurity:" << ineq << " spin:" << is << "========\n";

          for(int m=0; m<m_tot; m++)
            std::cout << std::setw(12) << std::fixed << std::setprecision(6) << 
              Hartree_to_eV*this->impurity_level[ineq][is][m*m_tot+m].real();

          std::cout << std::endl;
        }
      }
    }

    return;
  }

  void impurity::evaluate_Weiss_hybridization_imag(
        DFT_output::KS_bands& band, 
        DFT_plus_DMFT::projector& proj,
        DFT_output::atoms_info& atom,
        DFT_plus_DMFT::Hilbert_space& space,
        const double mu,
        const int nomega,
        const int mag )
  {
    debug::codestamp("impurity::evaluate_hybridization_function");

    const int nspin = band.nspins();
    const int nks = band.nk();
    const std::vector<int> wbands = space.Wbands();
    const int nprocs = mpi_ntasks();
    const int myid = mpi_rank();
    const int ineq_num = atom.inequ_atoms();
    const std::vector<int>& norb_sub = atom.iatom_norb();
    const std::vector<std::vector<int>>& orb_index = atom.Im2iorb();
    const int norb = atom.norb();
    const std::vector<double>& freq = this->sigma.sigma_imag.Matsubara_freq();
    const std::vector<double>& fk = band.kweight();

    const std::complex<double> zero(0.0,0.0), im(0.0,1.0), one(1.0,0.0);

    //==========================================================================
    //                        PART 1:Local Green function
    //==========================================================================
    //Allocation and set zero
    if(this->Green_fun_omega.empty())
    {
      this->Green_fun_omega.resize(ineq_num);
      for(int ineq=0; ineq<ineq_num; ineq++)
      {
        const int iatom = atom.ineq_iatom(ineq);
        const int m_tot=norb_sub[iatom];

        this->Green_fun_omega[ineq].resize(nspin);
        for(int is=0; is<nspin; is++)
        {
          this->Green_fun_omega[ineq][is].resize(nomega);
          for(int iomega=0; iomega<nomega; iomega++)
          {
            this->Green_fun_omega[ineq][is][iomega].resize(m_tot*m_tot,zero);
          }
        }
      }
    }
    else
    {
      for(int ineq=0; ineq<ineq_num; ineq++)
      {
        const int iatom = atom.ineq_iatom(ineq);
        const int m_tot=norb_sub[iatom];

        auto& Gfa = this->Green_fun_omega[ineq];
        for(int is=0; is<nspin; is++)
        {
          auto& Gfb = Gfa[is];
          for(int iomega=0; iomega<nomega; iomega++)
          {
            auto& Gfc = Gfb[iomega];
            for(int m_index=0; m_index<m_tot*m_tot; m_index++)
            {
              Gfc[m_index] = zero;
            }//m_index
          }//iomge
        }//is
      }//ineq
    }
    
    //evaluation
    int ik_count=0;
    for(int ik=0; ik<nks; ik++)
    {
      if(ik%nprocs != myid) continue;

      this->sigma.evalute_lattice_sigma(
            0, mag, nspin, wbands, atom,
            proj.proj_access(ik_count) );

      const std::vector<std::vector<std::vector<std::complex<double>>>>&
            latt_sigma = this->sigma.lattice_sigma(0);

      for(int is=0; is<nspin; is++)
      {
        const std::vector<double>& epsilon = space.eigen_val()[is][ik];    
        
        const int mkl_threads = mkl_get_max_threads();
        mkl_set_num_threads(1);  //set the number of threads of MKL library function to 1

        #pragma omp parallel
        {
          std::unique_ptr<std::complex<double>[]> Gf_latt(new std::complex<double> [wbands[is]*wbands[is]]);
          std::unique_ptr<int[]> ipiv(new int [wbands[is]]);
          
          #pragma omp for
          for(int iomega=0; iomega<nomega; iomega++)
          {
            const std::vector<std::complex<double>>& latt_sigmb = latt_sigma[is][iomega];

            for(int iband1=0; iband1<wbands[is]; iband1++)
            {
              for(int iband2=0; iband2<wbands[is]; iband2++)
              {
                int band_index = iband1*wbands[is]+iband2;

                if(iband1==iband2)
                    Gf_latt[band_index] = im*freq[iomega] + mu
                          -epsilon[iband1] - latt_sigmb[band_index];
                else
                    Gf_latt[band_index] = -latt_sigmb[band_index];

              }//iband2
            }//iband1

            general_complex_matrix_inverse(&Gf_latt[0], wbands[is], &ipiv[0]);

            for(int ineq=0; ineq<ineq_num; ineq++)
            {
              const int iatom = atom.ineq_iatom(ineq);
              const int m_tot=norb_sub[iatom];

              auto& Gf = this->Green_fun_omega[ineq][is];

              const std::vector<std::complex<double>>& projector = proj.proj_access(ik_count)[iatom][is];

              for(int m1=0; m1<m_tot; m1++)
              {
                for(int m2=0; m2<m_tot; m2++)
                {
                  const int m_index = m1*m_tot+m2;
                  for(int iband1=0; iband1<wbands[is]; iband1++)
                  {
                    int index1 = iband1*m_tot + m1;
                    for(int iband2=0; iband2<wbands[is]; iband2++)
                    {
                      int index2 = iband2*m_tot + m2;
                      Gf[iomega][m_index] += std::conj(projector[index1])*Gf_latt[iband1*wbands[is]+iband2]
                                              *projector[index2]*fk[ik];
                    }//iband2
                  }//iband1
                }//m2
              }//m1
            }//ineq

          }//iomega
        }//omp parallel
        mkl_set_num_threads(mkl_threads);
      }//is
      ik_count++;
    }//ik
    
    for(int ineq=0; ineq<ineq_num; ineq++)
    {
      const int iatom = atom.ineq_iatom(ineq);
      const int m_tot=norb_sub[iatom];
      
      //reduce operation  
      for(int is=0; is<nspin; is++)
      {
        for(int iomega=0; iomega<nomega; iomega++)
        {
          std::vector<std::complex<double>> tmp = this->Green_fun_omega[ineq][is][iomega];

          MPI_Allreduce( &tmp[0], &this->Green_fun_omega[ineq][is][iomega][0],
                        m_tot*m_tot, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD );
        }
      }

      //=============================================
      //   Local symmetry operation
      //=============================================
      const int local_symmetry = atom.local_sym();
      const int corr_L=atom.L(ineq);
      for(int is=0; is<nspin; is++)
        #pragma omp parallel for
        for(int iomega=0; iomega<nomega; iomega++)
          DFT_output::atoms_info::symmetry_operation_matrix<std::complex<double>>(
            local_symmetry, corr_L, m_tot, &this->Green_fun_omega[ineq][is][iomega][0] );

// tests
// if(mpi_rank()==0)
// {
// std::ofstream ofs_gf("Gf_omega.dat", std::ios::out);
// for(int iomega=0; iomega<nomega; iomega++)
// {
//   for(int is=0; is<2; is++)
//   {
//     for(int m1=0; m1<m_tot; m1++)
//        for(int m2=0; m2<m_tot; m2++)
//         if(nspin==1)
//           ofs_gf << std::setw(5) << iomega
//            << std::setw(2) << is << std::setw(3) << m1 << std::setw(3) << m2
//            << std::setw(15) << std::fixed << std::setprecision(9)
//            << this->Green_fun_omega[ineq][0][iomega][m1*m_tot+m2].real()/Hartree_to_eV
//            << std::setw(15) << std::fixed << std::setprecision(9) 
//            << this->Green_fun_omega[ineq][0][iomega][m1*m_tot+m2].imag()/Hartree_to_eV
//            << std::setw(15) << std::fixed << std::setprecision(9) 
//            << std::sqrt(std::pow(this->Green_fun_omega[ineq][0][iomega][m1*m_tot+m2].imag()/Hartree_to_eV,2) +
//               std::pow(this->Green_fun_omega[ineq][0][iomega][m1*m_tot+m2].real()/Hartree_to_eV,2))<< '\n';
//         else
//           ofs_gf << std::setw(5) << iomega 
//            << std::setw(2) << is << std::setw(3) << m1 << std::setw(3) << m2
//            << std::setw(15) << std::fixed << std::setprecision(9)
//            << this->Green_fun_omega[ineq][is][iomega][m1*m_tot+m2].real()/Hartree_to_eV
//            << std::setw(15) << std::fixed << std::setprecision(9) 
//            << this->Green_fun_omega[ineq][is][iomega][m1*m_tot+m2].imag()/Hartree_to_eV 
//            << std::setw(15) << std::fixed << std::setprecision(9) 
//            << std::sqrt(std::pow(this->Green_fun_omega[ineq][is][iomega][m1*m_tot+m2].imag()/Hartree_to_eV,2) +
//               std::pow(this->Green_fun_omega[ineq][is][iomega][m1*m_tot+m2].real()/Hartree_to_eV,2))<< '\n';  
//   }
// }
// ofs_gf.close();
// }
    }//ineq

    //==========================================================================
    //             PART 2:Hybridization function and Weiss Green's function
    //==========================================================================
    //Allocation
    if(this->hyb_omega.empty())
    {
      this->hyb_omega.resize(ineq_num);
      for(int ineq=0; ineq<ineq_num; ineq++)
      {
        const int iatom = atom.ineq_iatom(ineq);
        const int m_tot=norb_sub[iatom];
  
        this->hyb_omega[ineq].resize(nspin);
        for(int is=0; is<nspin; is++)
        {
          this->hyb_omega[ineq][is].resize(nomega);
          for(int iomega=0; iomega<nomega; iomega++)
          {
            this->hyb_omega[ineq][is][iomega].resize(m_tot*m_tot, zero);
          }
        }
      }
    }

    if(this->Weiss_omega.empty())
    {
      this->Weiss_omega.resize(ineq_num);
      for(int ineq=0; ineq<ineq_num; ineq++)
      {
        const int iatom = atom.ineq_iatom(ineq);
        const int m_tot=norb_sub[iatom];
  
        this->Weiss_omega[ineq].resize(nspin);
        for(int is=0; is<nspin; is++)
        {
          this->Weiss_omega[ineq][is].resize(nomega);
          for(int iomega=0; iomega<nomega; iomega++)
            this->Weiss_omega[ineq][is][iomega].resize(m_tot*m_tot, zero);
        }
      }
    }

    // evaluation
    for(int ineq=0; ineq<ineq_num; ineq++)
    {
      const int iatom = atom.ineq_iatom(ineq);
      const int m_tot=norb_sub[iatom];

      std::vector<std::vector<std::vector<std::complex<double>>>>&
            hyb = this->hyb_omega[ineq];

      std::vector<std::vector<std::vector<std::complex<double>>>>&
            G0 = this->Weiss_omega[ineq];

      const std::vector<std::vector<std::vector<std::complex<double>>>>&
            local_sigma = this->sigma.sigma_new(0)[ineq];

      #pragma omp parallel for
      for(int iomega=0; iomega<nomega; iomega++)
      {
        for(int is=0; is<nspin; is++)
        {
          for(int m=0; m<m_tot; m++)
          {
            G0[is][iomega][m*m_tot+m] = one/this->Green_fun_omega[ineq][is][iomega][m*m_tot+m] + local_sigma[is][iomega][m*m_tot+m];
          
            hyb[is][iomega][m*m_tot+m] = im*freq[iomega] + mu - this->impurity_level[ineq][is][m*m_tot+m] 
                                      - G0[is][iomega][m*m_tot+m];

            G0[is][iomega][m*m_tot+m] = one/G0[is][iomega][m*m_tot+m];                             
          }//m
        }//is
      }//iomega

//tests
// std::ofstream ofs_hyb("Hyb_omega.dat", std::ios::out);
// for(int iomega=0; iomega<nomega; iomega++)
// {
//   for(int is=0; is<nspin; is++)
//   {
//     for(int m1=0; m1<m_tot; m1++)
//       for(int m2=0; m2<m_tot; m2++)
//        ofs_hyb << std::setw(5) << iomega << std::setw(2) << is
//        << std::setw(3) << m1 << std::setw(3) << m2
//        << std::setw(15) << std::fixed << std::setprecision(9) 
//        << hyb[is][iomega][m1*m_tot+m2].real()
//        << std::setw(15) << std::fixed << std::setprecision(9) 
//        << hyb[is][iomega][m1*m_tot+m2].imag() << '\n';
//   }
// }
// ofs_hyb.close();

    }//ineq

    return;
  }

  void impurity::evaluate_delta_Weiss_tau(
                  const int impurity_solver,
                  DFT_output::atoms_info& atom,
                  const int nspin,
                  const double beta,
                  const int ntau,
                  const int nomega)
  {
    debug::codestamp("impurity::evaluate_delta_Weiss_tau");

    const int ineq_num = atom.inequ_atoms();
    const std::vector<int>& norb_sub = atom.iatom_norb();
    const std::complex<double> zero(0.0,0.0);
    const std::complex<double> im(0.0,1.0);

    std::vector<std::vector<std::vector<std::vector<std::complex<double>>>>>&
        hyb = this->hybdization_tau(impurity_solver);
    
    // std::vector<std::vector<std::vector<std::vector<std::complex<double>>>>> hyb;

    //Allocation and set zero
    if(hyb.empty())
    {
      hyb.resize(ineq_num);
      for(int ineq=0; ineq<ineq_num; ineq++)
      {
        const int iatom = atom.ineq_iatom(ineq);
        const int m_tot=norb_sub[iatom];

        hyb[ineq].resize(nspin);
        for(int is=0; is<nspin; is++)
        {
          hyb[ineq][is].resize(ntau+1);
          for(int itau=0; itau<ntau+1; itau++)
            hyb[ineq][is][itau].resize(m_tot*m_tot,zero);
        }
      }
    }
    else
    {
      for(int ineq=0; ineq<ineq_num; ineq++)
      {
        const int iatom = atom.ineq_iatom(ineq);
        const int m_tot=norb_sub[iatom];

        auto& hyba = hyb[ineq];
        for(int is=0; is<nspin; is++)
        {
          auto& hybb = hyba[is];
          for(int itau=0; itau<ntau+1; itau++)
          {
            auto& hybc = hybb[itau];
            for(int m_index=0; m_index<m_tot*m_tot; m_index++)
              hybc[m_index] = zero;
          }//itau
        }//is
      }//ineq
    }

    if(this->Weiss_tau.empty())
    {
      this->Weiss_tau.resize(ineq_num);
      for(int ineq=0; ineq<ineq_num; ineq++)
      {
        const int iatom = atom.ineq_iatom(ineq);
        const int m_tot=norb_sub[iatom];

        this->Weiss_tau[ineq].resize(nspin);
        for(int is=0; is<nspin; is++)
        {
          this->Weiss_tau[ineq][is].resize(ntau+1);
          for(int itau=0; itau<ntau+1; itau++)
            this->Weiss_tau[ineq][is][itau].resize(m_tot*m_tot,zero);
        }
      }
    }
    else
    {
      for(int ineq=0; ineq<ineq_num; ineq++)
      {
        const int iatom = atom.ineq_iatom(ineq);
        const int m_tot=norb_sub[iatom];

        auto& Weissa = this->Weiss_tau[ineq];
        for(int is=0; is<nspin; is++)
        {
          auto& Weissb = Weissa[is];
          for(int itau=0; itau<ntau+1; itau++)
          {
            auto& Weissc = Weissb[itau];
            for(int m_index=0; m_index<m_tot*m_tot; m_index++)
              Weissc[m_index] = zero;
          }//itqu
        }//is
      }//ineq
    }

    std::unique_ptr<double[]> tau(new double [ntau+1]);
    for(int itau=0; itau<ntau+1; itau++)
      tau[itau] = (beta/ntau)*itau;
    tau[0] = 1.0e-3*beta/ntau;
    tau[ntau] = beta - 1.0e-3*beta/ntau;

    std::vector<std::complex<double>> hyb_tau_tmp(ntau+1);
    std::vector<std::complex<double>> hyb_omega_tmp(nomega);

    //Fourier transform
    const int threads_num = omp_get_max_threads();
    for(int ineq=0; ineq<ineq_num; ineq++)
    {
      const int iatom = atom.ineq_iatom(ineq);
      const int m_tot=norb_sub[iatom];
      auto& hyb_tau1 = hyb[ineq];

      for(int is=0; is<nspin; is++)
      {
        for(int m=0; m<m_tot; m++)
        {
          //Hybridization function
          for(int iomega=0; iomega<nomega; iomega++)
            hyb_omega_tmp[iomega] = this->hyb_omega[ineq][is][iomega][m*m_tot+m];

          Fourier_trans_omega_tau( beta, ntau+1, nomega,
            &this->sigma.sigma_imag.Matsubara_freq()[0],
            &tau[0], &hyb_omega_tmp[0], &hyb_tau_tmp[0],
            threads_num );
       
          for(int itau=0; itau<ntau+1; itau++)
            hyb_tau1[is][itau][m*m_tot+m] = hyb_tau_tmp[itau];

          //Weiss Green's function
          for(int iomega=0; iomega<nomega; iomega++)
            hyb_omega_tmp[iomega] = this->Weiss_omega[ineq][is][iomega][m*m_tot+m];

          Fourier_trans_omega_tau( beta, ntau+1, nomega,
            &this->sigma.sigma_imag.Matsubara_freq()[0],
            &tau[0], &hyb_omega_tmp[0], &hyb_tau_tmp[0],
            threads_num );
          
          for(int itau=0; itau<ntau+1; itau++)
            this->Weiss_tau[ineq][is][itau][m*m_tot+m] = hyb_tau_tmp[itau];
        }//m
      }//is
//===========================================
//           write delta_tau.txt
//===========================================
// if(mpi_rank()==0)
// {
//   std::ofstream ofs("Delta_tau.dat", std::ios::out);
//   for(int itau=0; itau<ntau+1; itau++)
//   {
//     ofs << std::setw(5) << itau;
  
//     for(int m=0; m<m_tot; m++)
//     {
//       if(nspin==2)
//       {
//         ofs << std::setw(22) << std::fixed << std::setprecision(15) << 
//         std::pow(Hartree_to_eV,2)*hyb[ineq][0][itau][m*m_tot+m].real()
//         << std::setw(22) << std::fixed << std::setprecision(15)
//         << std::pow(Hartree_to_eV,2)*hyb[ineq][1][itau][m*m_tot+m].real();
//       }
//       else
//       {
//         ofs << std::setw(22) << std::fixed << std::setprecision(15) << 
//         std::pow(Hartree_to_eV,2)*hyb[ineq][0][itau][m*m_tot+m].real()
//         << std::setw(22) << std::fixed << std::setprecision(15)
//         << std::pow(Hartree_to_eV,2)*hyb[ineq][0][itau][m*m_tot+m].real();
//       }
//     }//m_index
//     ofs << '\n';
//   }//itau
//   ofs.close();
// }
    }//ineq

    return;
  }

  void impurity::out(const int char_step,
                    const int DMFT_step,
                    const int impurity_solver,
                    const double mu, 
                    DFT_output::KS_bands& band,
                    DMFT::input_info& in, 
                    DFT_output::atoms_info& atom,
                    DMFT::coulomb_tensor& Umat )
  {
    switch(impurity_solver)
    {
      case 1:
        this->ALPS_hyb.output(char_step, DMFT_step, 
               mu, in, atom, band, this->impurity_level, 
               this->sigma.sigma_new(0), this->Weiss_omega, this->hyb_omega);
        break;
      case 2:
        this->ALPS_hyb_segment.output(char_step, DMFT_step,
               mu, in, atom, band, this->impurity_level, 
               this->sigma.sigma_new(0), this->Weiss_omega);
        break;
      case 3:
        this->pacs.output(char_step, DMFT_step, 
               mu, in, atom, band, this->sigma.sigma_imag.Matsubara_freq(), 
               this->impurity_level, this->sigma.sigma_new(0),
               this->Weiss_omega, this->hyb_omega );
        break;
      case 4:
        this->Rutgers.output(char_step, DMFT_step, 
              mu, in, atom, band, this->sigma.sigma_imag.Matsubara_freq(), 
              this->impurity_level, this->sigma.sigma_new(0),
              this->Weiss_omega, this->hyb_omega, Umat );
        break;
      case 5: 
        this->iQIST_narcissus.output(char_step, DMFT_step, 
              mu, in, atom, band, this->sigma.sigma_imag.Matsubara_freq(), 
              this->impurity_level, this->sigma.sigma_new(0),
              this->Weiss_omega, this->hyb_omega, Umat );
        break;
      default:
        std::cout << "Not supported impurity solver" << std::endl;
        std::exit(EXIT_FAILURE);
    }
    return;
  }

  void impurity::read_last_step(
        const int char_step,
        const int DMFT_step,
        const int impurity_solver, 
        DFT_output::KS_bands& band,
        DMFT::input_info& in, 
        DFT_output::atoms_info& atom)
  {
    debug::codestamp("impurity::read_last_step");

    // if(mpi_rank()==0){
    //   std::cout << "Reading the last loop: charge step " << char_step << "  DMFT step " << DMFT_step << std::endl;
    // }

    const int ineq_num = atom.inequ_atoms();
    const int nspin = band.nspins();
    const std::vector<int>& norb_sub = atom.iatom_norb();
    const int nomega = *(int*)in.parameter("n_omega");
    const std::complex<double> zero(0.0,0.0);

    // this->Green_fun_omega.resize(ineq_num);
    // this->hyb_omega.resize(ineq_num);
    // this->Weiss_omega.resize(ineq_num);
    // this->Green_fun_omega_save.resize(ineq_num);
    // for(int ineq=0; ineq<ineq_num; ineq++)
    // {
    //   const int iatom = atom.ineq_iatom(ineq);
    //   const int m_tot=norb_sub[iatom];

    //   this->Green_fun_omega[ineq].resize(nspin);
    //   this->hyb_omega[ineq].resize(nspin);
    //   this->Weiss_omega[ineq].resize(nspin);
    //   this->Green_fun_omega_save[ineq].resize(nspin);

    //   for(int is=0; is<nspin; is++)
    //   {
    //     this->Green_fun_omega[ineq][is].resize(nomega);
    //     this->hyb_omega[ineq][is].resize(nomega);
    //     this->Weiss_omega[ineq][is].resize(nomega);
    //     this->Green_fun_omega_save[ineq][is].resize(nomega);
    //     for(int iomega=0; iomega<nomega; iomega++)
    //     {
    //       this->Green_fun_omega[ineq][is][iomega].resize(m_tot*m_tot,zero);
    //       this->hyb_omega[ineq][is][iomega].resize(m_tot*m_tot,zero);
    //       this->Weiss_omega[ineq][is][iomega].resize(m_tot*m_tot,zero);
    //       this->Green_fun_omega_save[ineq][is][iomega].resize(m_tot*m_tot,zero);
    //     }
    //   }//is
    // }//ineq

    //Self-energy allocation
    std::vector<std::vector<std::vector<std::vector<std::complex<double>>>>>&
        sigma_new = this->sigma.sigma_imag.sigma_new_access();

    if(sigma_new.empty())
    {
      sigma_new.resize(ineq_num);
      for(int ineq=0; ineq<ineq_num; ineq++)
      {
        const int iatom = atom.ineq_iatom(ineq);
        const int m_tot = norb_sub[iatom];

        sigma_new[ineq].resize(nspin);
        for(int is=0; is<nspin; is++)
        {
          sigma_new[ineq][is].resize(nomega);
          for(int iomega=0; iomega<nomega; iomega++)
          {
            sigma_new[ineq][is][iomega].resize(m_tot*m_tot,zero);
          }//is
        }//i_omega
      }//ineq
    }

    std::vector<std::vector<std::vector<std::vector<std::complex<double>>>>>&
        sigma_save = this->sigma.sigma_imag.sigma_save_access();

    if(sigma_save.empty())
    {
      sigma_save.resize(ineq_num);
      for(int ineq=0; ineq<ineq_num; ineq++)
      {
        const int iatom = atom.ineq_iatom(ineq);
        const int m_tot = norb_sub[iatom];

        sigma_save[ineq].resize(nspin);
        for(int is=0; is<nspin; is++)
        {
          sigma_save[ineq][is].resize(nomega);
          for(int iomega=0; iomega<nomega; iomega++)
          {
            sigma_save[ineq][is][iomega].resize(m_tot*m_tot,zero);
          }//is
        }//i_omega
      }//ineq
    }

    switch(impurity_solver)
    {
      case 1:
        // this->ALPS_hyb.read_last_step(
        //         char_step, DMFT_step, 
        //         band, in, atom,
        //         this->Green_fun_omega, 
        //         this->Weiss_omega, 
        //         this->Green_fun_omega_save );
        break;
      case 2:
        // this->ALPS_hyb_segment.read_last_step(istep, band, in, atom, this->hyb_omega,
        //                     this->Green_fun_omega, this->Green_fun_omega_save);
        break;
      case 3:
        this->pacs.read_last_step(
                char_step, DMFT_step,
                band, in, atom,
                sigma_new, sigma_save);
        break;
      case 4:
        this->Rutgers.read_last_step(
                char_step, DMFT_step,
                band, in, atom, 
                sigma_new, sigma_save );
        break;
      case 5:
        this->iQIST_narcissus.read_last_step(
                char_step, DMFT_step,
                band, in, atom, 
                sigma_new, sigma_save );
        break;
      default:
        std::cout << "Not supported impurity solver" << std::endl;
        std::exit(EXIT_FAILURE);
    }

// if(mpi_rank()==0)
// {
//   std::ofstream ofs("G_omga.dat", std::ios::out);
//   for(int ineq=0; ineq<ineq_num; ineq++)
//   {
//     const int iatom = atom.ineq_iatom(ineq);
//     const int m_tot=norb_sub[iatom];
    
//     for(int iomega=0; iomega<nomega; iomega++)
//     {
//       for(int is=0; is<2; is++)
//       {
//         for(int m_index=0; m_index<m_tot*m_tot; m_index++)
//         {
//           int m1 = m_index/m_tot;
//             int m2 = m_index%m_tot;

//             if(nspin==1)
//               ofs << std::left << std::setw(5) << iomega << std::setw(2) << is 
//                   << std::setw(3) << m1 << std::setw(3) << m2
//                   << std::setw(20) << std::fixed << std::setprecision(12) 
//                   << this->Green_fun_omega[ineq][0][iomega][m_index].real()/Hartree_to_eV << " "
//                   << std::setw(20) << std::fixed << std::setprecision(12) 
//                   << this->Green_fun_omega[ineq][0][iomega][m_index].imag()/Hartree_to_eV << '\n';
//             else
//               ofs << std::left << std::setw(5) << iomega << std::setw(2) << is 
//                   << std::setw(3) << m1 << std::setw(3) << m2
//                   << std::setw(20) << std::fixed << std::setprecision(12) 
//                   << this->Green_fun_omega[ineq][is][iomega][m_index].real()/Hartree_to_eV << " "
//                   << std::setw(20) << std::fixed << std::setprecision(12) 
//                   << this->Green_fun_omega[ineq][is][iomega][m_index].imag()/Hartree_to_eV << '\n';
//         }
//       }
//     }//is
//   }//ineq
//   ofs.close();
// }

    return;
  }

  bool impurity::scf_condition(
        const int flag_axis, 
        DFT_output::KS_bands& band,
        DFT_output::atoms_info& atom,
        DMFT::input_info& in )
  {
    debug::codestamp("impurity::scf_condition");

    bool convergency=true;

    const int ineq_num = atom.inequ_atoms();
    const std::vector<int>& norb_sub = atom.iatom_norb();
    const int nspin = band.nspins();
    const int nomega = *(int*)in.parameter("n_omega");
    const double sc_criteria = *(double*)in.parameter("delta_sigma");

    const std::vector<std::vector<std::vector<std::vector<std::complex<double>>>>>&
        sigma_new = this->sigma.sigma_new(flag_axis);

    const std::vector<std::vector<std::vector<std::vector<std::complex<double>>>>>&
        sigma_save = this->sigma.sigma_save(flag_axis);

    this->scf_delta.resize(ineq_num,0.0);
    std::unique_ptr<bool[]> flag_conver(new bool [ineq_num]);

    for(int ineq=0; ineq<ineq_num; ineq++)
    {
      const int iatom = atom.ineq_iatom(ineq);
      const int m_tot=norb_sub[iatom];

      for(int is=0; is<nspin; is++)
      {
        for(int iomega=0; iomega<nomega; iomega++)
        {
          for(int m=0; m<m_tot; m++)
          {
            double tmp = std::sqrt(std::norm( sigma_new[ineq][is][iomega][m*m_tot+m] 
                        - sigma_save[ineq][is][iomega][m*m_tot+m] ));
            if(tmp>this->scf_delta[ineq]) this->scf_delta[ineq]=tmp;
          }//m_index
        }//iomega
      }//is
      this->scf_delta[ineq] *= Hartree_to_eV;
      // this->scf_delta[ineq] = this->scf_delta[ineq]/(nspin*nomega*m_tot);

      if(this->scf_delta[ineq]>sc_criteria)
        flag_conver[ineq] = false;
      else
        flag_conver[ineq] = true;

    }//ineq

    for(int ineq=0; ineq<ineq_num; ineq++)
    {
      if(!flag_conver[ineq])
      {
        convergency=false;
        break;
      }
    }

    return convergency;
  }

  void impurity::update_self_energy( 
              DFT_output::KS_bands& band,
              DFT_output::atoms_info& atom,
              DMFT::input_info& in)
  {
    debug::codestamp("impurity::update_self_energy");

    const double beta = *(double*)in.parameter("beta");
    const std::complex<double> im = std::complex<double>(0.0,1.0);

    const int ineq_num = atom.inequ_atoms();
    const std::vector<int>& norb_sub = atom.iatom_norb();
    const int nspin = band.nspins();
    const int nomega = *(int*)in.parameter("n_omega"); 

    std::vector<std::vector<std::vector<std::vector<std::complex<double>>>>>&
        sigma_new = this->sigma.sigma_imag.sigma_new_access();
    
    std::vector<std::vector<std::vector<std::vector<std::complex<double>>>>>&
        sigma_correlatd = this->sigma.sigma_imag.correlated_sigma_access();

    std::vector<double>& freq=this->sigma.sigma_imag.Matsubara_freq();

    //evalute self-energy by Dyson equation
    for(int ineq=0; ineq<ineq_num; ineq++)
    {
      const int iatom = atom.ineq_iatom(ineq);
      const int m_tot = norb_sub[iatom];

      std::vector<std::vector<std::vector<std::complex<double>>>>
          Gf_omega = this->Green_fun_omega[ineq];

      std::vector<std::vector<std::vector<std::complex<double>>>>
          Weiss = this->Weiss_omega[ineq];

      std::vector<std::vector<std::vector<std::complex<double>>>>&
          sigma_newa =sigma_new[ineq];

      const int mkl_threads = mkl_get_max_threads();
      mkl_set_num_threads(1);      //set the number of threads of MKL library function to 1

      #pragma omp parallel
      {
        std::unique_ptr<int[]> ipiv(new int [m_tot]);

        #pragma omp for
        for(int iomega=0; iomega<nomega; iomega++)   
        {          
          for(int is=0; is<nspin; is++)
          {         
            general_complex_matrix_inverse(&Gf_omega[is][iomega][0], m_tot, &ipiv[0]);

            general_complex_matrix_inverse(&Weiss[is][iomega][0], m_tot, &ipiv[0]);
            
            for(int m_index=0; m_index<m_tot*m_tot; m_index++)
              sigma_newa[is][iomega][m_index] = Weiss[is][iomega][m_index] - Gf_omega[is][iomega][m_index];

          }//is
        }//iomega     
      } 
      mkl_set_num_threads(mkl_threads); 
    }//ineq

    return;
  }

  void impurity::evaluate_local_occupation(
                DMFT::input_info& in, 
                DFT_output::KS_bands& band, 
                DFT_plus_DMFT::projector& proj,
                DFT_output::atoms_info& atom,
                DFT_plus_DMFT::Hilbert_space& space,
                const double mu,
                const int nomega,
                const int mag )
  {
    debug::codestamp("impurity::evaluate_local_occupation");

    const int ineq_num = atom.inequ_atoms();
    const int nspin = band.nspins();
    const int nks = band.nk();
    const std::vector<int>& wbands = space.Wbands();
    const int nprocs = mpi_ntasks();
    const int myid = mpi_rank();
    const double beta= *(double*)in.parameter("beta");
    const std::vector<int>& sub_norb = atom.iatom_norb();
    const std::vector<std::vector<int>>& orb_index = atom.Im2iorb();
    const int norb = atom.norb();
    const std::vector<double>& fk = band.kweight();
    std::vector<std::vector<std::vector<double>>> dft_occ_num = band.dft_occ();
    const std::vector<std::vector<int>>& wb2ib = space.wbands2ibands();

    const std::complex<double> zero(0.0,0.0), im(0.0,1.0), one(1.0,0.0);
    
    std::vector<double> freq(nomega);
    for(int iomega=0; iomega<nomega; iomega++)
      freq[iomega] = (2*iomega+1)*PI/beta;

    //Allocation and set zero
    std::vector<std::vector<std::vector<std::vector<std::complex<double>>>>>
          local_gf(ineq_num);
    for(int ineq=0; ineq<ineq_num; ineq++)
    {
      const int iatom = atom.ineq_iatom(ineq);
      const int m_tot = sub_norb[iatom];

      local_gf[ineq].resize(nspin);
      for(int is=0; is<nspin; is++)
      {
        local_gf[ineq][is].resize(nomega);
        for(int iomega=0; iomega<nomega; iomega++)
          local_gf[ineq][is][iomega].resize(m_tot*m_tot,zero);
      }
    }

    std::vector<std::vector<std::vector<double>>>& local_occ = atom.occ_num_m();
    for(auto& iter1 : local_occ)
      for(auto& iter2 : iter1)
        for(auto& iter3 : iter2)
          iter3 = 0.0;

    std::vector<std::vector<std::vector<double>>> local_occ_tmp = local_occ;

    std::vector<std::vector<double>>& occ_num = atom.occ_num_ref();
    for(auto& iter1 : occ_num)
      for(auto& iter2 : iter1)
        iter2=0.0;

    if(mpi_rank()==0) std::cout << "\n===========Local occupation number========" << std::endl;

    //evaluation
    if(nspin==1 && !band.soc()){
      for(auto& iter1 : dft_occ_num)
        for(auto& iter2 : iter1)
          for(auto& iter3 : iter2)
            iter3 /= 2.0;
    }

    for(int ineq=0; ineq<ineq_num; ineq++)
    {
      const int iatom = atom.ineq_iatom(ineq);
      const int m_tot = sub_norb[iatom];

      int ik_count=0;
      for(int ik=0; ik<nks; ik++)
      {
        if(ik%nprocs != myid) continue;

        for(int is=0; is<nspin; is++)
        {
          const std::vector<double>& epsilon = space.eigen_val()[is][ik];     
          const std::vector<std::complex<double>>& projector = proj.proj_access(ik_count)[iatom][is];
      
            for(int m=0; m<m_tot; m++)
              for(int iband=0; iband<wbands[is]; iband++)
                local_occ_tmp[iatom][is][m] += 
                  ( std::conj(projector[iband*m_tot + m])*
                    dft_occ_num[is][ik][ wb2ib[is][iband] ]*fk[ik]*
                    projector[iband*m_tot + m] ).real();
        }//is
        ik_count++;
      }//ik

      //reduce operation
      for(int is=0; is<nspin; is++){
        MPI_Allreduce( &local_occ_tmp[iatom][is][0], 
                       &local_occ[iatom][is][0], m_tot, 
                       MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
      }

      if(nspin==1) local_occ[iatom][1] = local_occ[iatom][0];

      //=============================================
      //   Local symmetry operation
      //=============================================
      const int local_symmetry = atom.local_sym();
      const int corr_L=atom.L(ineq);
      for(int is=0; is<2; is++)
        DFT_output::atoms_info::symmetry_operation_vector<double>(
          local_symmetry, corr_L, m_tot, &local_occ[iatom][is][0] );
      //=============================================

      for(int is=0; is<2; is++)
        for(int m=0; m<m_tot; m++)
          occ_num[iatom][is] += local_occ[iatom][is][m];
      
      for(int iatom1=0; iatom1<atom.total_atoms(); iatom1++)
      {
        if(iatom1!=iatom && atom.equ_atom(iatom1)==ineq)
        {
          int symm=0;       //symmetry

          if(mag==1 && atom.magnetic(iatom)*atom.magnetic(iatom1)==-1) symm=1;  //anti_symmetry;AFM

          for(int ispin=0; ispin<2; ispin++)
            for(int m=0; m<m_tot; m++)
              local_occ[iatom1][ispin][m] = local_occ[iatom][(ispin+symm)%2][m];

          for(int is=0; is<2; is++)
            for(int m=0; m<m_tot; m++)
              occ_num[iatom1][is] += local_occ[iatom1][is][m];
        }
      }

      if(mpi_rank()==0)
      {
        for(int is=0; is<nspin; is++)
        {
          std::cout << "===========impurity:" << ineq << " spin:" << is << "========\n";
          for(int m=0; m<m_tot; m++)
            if(nspin==1)
              std::cout << std::setw(12) << std::fixed << std::setprecision(6) << local_occ[iatom][0][m];
            else
              std::cout << std::setw(12) << std::fixed << std::setprecision(6) << local_occ[iatom][is][m]; 
          std::cout << std::endl;
        }
      }

    }//ineq
    
    //evaluation
    /*
    for(int ineq=0; ineq<ineq_num; ineq++)
    {
      const int iatom = atom.ineq_iatom(ineq);
      const int m_tot = sub_norb[iatom];

      int ik_count=0;
      for(int ik=0; ik<nks; ik++)
      {
        if(ik%nprocs != myid) continue;

        for(int is=0; is<nspin; is++)
        {
          const std::vector<double>& epsilon = space.eigen_val()[is][ik];     
          const std::vector<std::complex<double>>& projector = proj.proj_access(ik_count)[iatom][is];
      
            for(int m=0; m<m_tot; m++)
              for(int iband=0; iband<wbands[is]; iband++)
                local_occ_tmp[iatom][is][m] += 
                  ( std::conj(projector[iband*m_tot + m])/
                    (1.0+std::exp(beta*(epsilon[iband]-mu)))*
                    projector[iband*m_tot + m]*fk[ik] ).real();
        }//is
        ik_count++;
      }//ik

      //reduce operation
      for(int is=0; is<nspin; is++){
        MPI_Allreduce( &local_occ_tmp[iatom][is][0], 
                       &local_occ[iatom][is][0], m_tot, 
                       MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
      }

      if(nspin==1) local_occ[iatom][1] = local_occ[iatom][0];

      //=============================================
      //   Local symmetry operation
      //=============================================
      const int local_symmetry = atom.local_sym();
      const int corr_L=atom.L(ineq);
      for(int is=0; is<2; is++)
        DFT_output::atoms_info::symmetry_operation_vector<double>(
          local_symmetry, corr_L, m_tot, &local_occ[iatom][is][0] );
      //=============================================

      for(int is=0; is<2; is++)
        for(int m=0; m<m_tot; m++)
          occ_num[iatom][is] += local_occ[iatom][is][m];
      
      for(int iatom1=0; iatom1<atom.total_atoms(); iatom1++)
      {
        if(iatom1!=iatom && atom.equ_atom(iatom1)==ineq)
        {
          int symm=0;       //symmetry

          if(mag==1 && atom.magnetic(iatom)*atom.magnetic(iatom1)==-1) symm=1;  //anti_symmetry;AFM

          for(int ispin=0; ispin<2; ispin++)
            for(int m=0; m<m_tot; m++)
              local_occ[iatom1][ispin][m] = local_occ[iatom][(ispin+symm)%2][m];

          for(int is=0; is<2; is++)
            for(int m=0; m<m_tot; m++)
              occ_num[iatom1][is] += local_occ[iatom1][is][m];
        }
      }

      if(mpi_rank()==0)
      {
        for(int is=0; is<nspin; is++)
        {
          std::cout << "===========impurity:" << ineq << " spin:" << is << "========\n";
          for(int m=0; m<m_tot; m++)
            if(nspin==1)
              std::cout << std::setw(12) << std::fixed << std::setprecision(6) << local_occ[iatom][0][m];
            else
              std::cout << std::setw(12) << std::fixed << std::setprecision(6) << local_occ[iatom][is][m]; 
          std::cout << std::endl;
        }
      }

    }//ineq */

    /*
    for(int ineq=0; ineq<ineq_num; ineq++)
    {
      const int iatom = atom.ineq_iatom(ineq);
      const int m_tot = sub_norb[iatom];

      int ik_count=0;
      for(int ik=0; ik<nks; ik++)
      {
        if(ik%nprocs != myid) continue;

        for(int is=0; is<nspin; is++)
        {
          auto& Gf = local_gf[ineq][is];
          const std::vector<double>& epsilon = space.eigen_val()[is][ik];     
          const std::vector<std::complex<double>>& projector = proj.proj_access(ik_count)[iatom][is];
      
          #pragma omp parallel for
          for(int iomega=0; iomega<nomega; iomega++)
            for(int m=0; m<m_tot; m++)
              for(int iband=0; iband<wbands[is]; iband++)
                Gf[iomega][m*m_tot+m] += std::conj(projector[iband*m_tot + m])/
                (im*freq[iomega] + mu - epsilon[iband])*projector[iband*m_tot + m]*fk[ik];
        }//is
        ik_count++;
      }//ik

      //reduce operation
      for(int is=0; is<nspin; is++)
      {
        for(int iomega=0; iomega<nomega; iomega++)
        {
          std::vector<std::complex<double>> tmp = local_gf[ineq][is][iomega];

          MPI_Allreduce( &tmp[0], &local_gf[ineq][is][iomega][0], m_tot*m_tot, 
                         MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD );
        }
      }

      //evaluate local occupation
     
      for(int m=0; m<m_tot; m++)
      {
        for(int is=0; is<nspin; is++)
        {
          const auto& gf = local_gf[ineq][is];
          double real, imag;
          double c1, c2;
          c1 = -gf[nomega-1][m*m_tot+m].imag()*freq[nomega-1];
          c2 = -gf[nomega-1][m*m_tot+m].real()*std::pow(freq[nomega-1],2);

          double sum = 0.0;
          for(int iomega=0; iomega<nomega; iomega++)
          {
            const double omegan = freq[iomega];

            std::complex<double> val_tmp = gf[iomega][m*m_tot+m] + im*c1/omegan
                                          + c2/std::pow(omegan,2);
            
            sum -= 2.0*val_tmp.real()/beta;
          }//iomega
          sum += (-0.5*c1 + c2*beta/4.0);
          local_occ[iatom][is][m] -= sum;
        }//is
        if(nspin==1) local_occ[iatom][1][m] = local_occ[iatom][0][m];
      }//m

      //=============================================
      //   Local symmetry operation
      //=============================================
      const int local_symmetry = atom.local_sym();
      const int corr_L=atom.L(ineq);
      for(int is=0; is<2; is++)
        DFT_output::atoms_info::symmetry_operation_vector<double>(
          local_symmetry, corr_L, m_tot, &local_occ[iatom][is][0] );
      //=============================================

      for(int is=0; is<2; is++)
        for(int m=0; m<m_tot; m++)
          occ_num[iatom][is] += local_occ[iatom][is][m];
      
      for(int iatom1=0; iatom1<atom.total_atoms(); iatom1++)
      {
        if(iatom1!=iatom && atom.equ_atom(iatom1)==ineq)
        {
          int symm=0;       //symmetry

          if(mag==1 && atom.magnetic(iatom)*atom.magnetic(iatom1)==-1) symm=1;  //anti_symmetry;AFM

          for(int ispin=0; ispin<2; ispin++)
            for(int m=0; m<m_tot; m++)
              local_occ[iatom1][ispin][m] = local_occ[iatom][(ispin+symm)%2][m];

          for(int is=0; is<2; is++)
            for(int m=0; m<m_tot; m++)
              occ_num[iatom1][is] += local_occ[iatom1][is][m];
        }
      }

      if(mpi_rank()==0)
      {
        for(int is=0; is<nspin; is++)
        {
          std::cout << "===========impurity:" << ineq << " spin:" << is << "========\n";
          for(int m=0; m<m_tot; m++)
            if(nspin==1)
              std::cout << std::setw(12) << std::fixed << std::setprecision(6) << local_occ[iatom][0][m];
            else
              std::cout << std::setw(12) << std::fixed << std::setprecision(6) << local_occ[iatom][is][m]; 
          std::cout << std::endl;
        }
      }

    }//ineq */

    return;
  }

  std::vector<std::vector<std::vector<std::vector<std::complex<double>>>>>& 
      impurity::hybdization_tau(const int impurity_solver)
  {
    switch(impurity_solver)
    {
      case 1: //ALPS_CTHYB
        return this->ALPS_hyb.hybridization_func_tau();
        break;
      case 2: //
        return this->ALPS_hyb_segment.hybridization_func_tau();
        break;
      case 3: //pacs
        return this->pacs.hybridization_func_tau();
        break;
      case 4: //Rutgers
        return this->Rutgers.hybridization_func_tau();
        break;
      case 5:  //iQIST
        return this->iQIST_narcissus.hybridization_func_tau();
        break;
      default:
        std::cout << "Not supported impurity solver" << std::endl;
        std::exit(EXIT_FAILURE);
    }
  }

  void impurity::Fourier_trans_omega_tau(
            const double beta,
            const int ntau,
            const int nomega,
            const double* freq,
            const double* tau_mesh,
            const std::complex<double>* fw,
            std::complex<double>* ftau,
            const int threads_num)
  {
    debug::codestamp("impurity::Fourier_trans_omega_tau");

    const std::complex<double> zero(0.0,0.0), im(0.0,1.0);
    for(int itau=0; itau<ntau; itau++)
      ftau[itau] = zero;

    //  ================Part 1 : polynomial regression================
    //  Calculate the imaginary-time hybridization function Delta(tau) using the discrete
    //  Matsubara frequency summation augmented with high-frequency tail:
    //  f(iw_n) = C1/iw_n + C2/(iw_n)^2 + c3/(iw_n)^3
    //  f(tau) = -C1/2 + C2/4*(-beta + 2tau) + C3/4(beta tau - tau^2)
    //  NOTE: C1, C2, C3 are averaged over the last 20 Matsubara frequency points
    
    std::vector<std::complex<double>> fw_tail(20);
    std::vector<double> freq_tail(20);
    for(int iomega=0; iomega<20; iomega++)
    {
      fw_tail[iomega] = fw[nomega-20+iomega];
      freq_tail[iomega] = freq[nomega-20+iomega];
    }

    double C1, C2, C3;
    polynomial_regression(&fw_tail[0], &freq_tail[0], 20, C1, C2, C3);

    #pragma omp parallel for num_threads(threads_num)
    for(int itau=0; itau<ntau; itau++)
    {
      const double tau_val = tau_mesh[itau];                
      for(int iomega=0; iomega<nomega; iomega++)
      {
        const double omegan = freq[iomega];
        std::complex<double> val_tmp = fw[iomega] + im*C1/omegan
                              + C2/std::pow(omegan,2) - im*C3/std::pow(omegan,3);

        ftau[itau] += ( val_tmp*(std::cos(omegan*tau_val)
                        -im*std::sin(omegan*tau_val)) 
                        + std::conj(val_tmp)*(std::cos(omegan*tau_val)
                        + im*std::sin(omegan*tau_val))
                        )/beta;
      }//iomega

      ftau[itau] += -0.5*C1 + C2*(2.0*tau_val-beta)/4.0 + 
                    C3*(beta*tau_val-tau_val*tau_val)/4.0;
                              
//===========TESTS===================
// const double pi=3.14159265358979323846;
// std::ofstream ofs_hwr("H0w_real.dat",std::ios::out);
// for(int iomega=0; iomega<nomega+200; iomega++)
// {
//   double omegan = (2*iomega+1)*pi/beta;
//   ofs_hwr << std::setw(5) << iomega << std::setw(8) << std::fixed << std::setprecision(3) << omegan;
//   for(int is=0; is<nspin; is++)
//   {
//     for(int m=0; m<m_tot; m++)
//     {
//       if(iomega<nomega) ofs_hwr << std::setw(15) << std::fixed << std::setprecision(9) << hyb_omega[ineq][is][iomega][m*m_tot+m].real();
//       else  ofs_hwr << std::setw(15) << std::fixed << std::setprecision(9) << -C2[is][m]/std::pow(omegan,2);
//     }//m
//   }//is
//   ofs_hwr << '\n';
// }//iomega
// ofs_hwr.close();

// std::ofstream ofs_hwi("H0w_imag.dat",std::ios::out);
// for(int iomega=0; iomega<nomega+200; iomega++)
// {
//   double omegan = (2*iomega+1)*pi/beta;
//   ofs_hwi << std::setw(5) << iomega << std::setw(8) << std::fixed << std::setprecision(3) << omegan;
//   for(int is=0; is<nspin; is++)
//   {
//     for(int m=0; m<m_tot; m++)
//     {
//       if(iomega<nomega) ofs_hwi << std::setw(15) << std::fixed << std::setprecision(9) << hyb_omega[ineq][is][iomega][m*m_tot+m].imag();
//       else  ofs_hwi << std::setw(15) << std::fixed << std::setprecision(9) << -C1[is][m]/omegan;
//     }//m
//   }//is
//   ofs_hwi << '\n';
// }//iomega
// ofs_hwi.close();

// std::ofstream ofs_ht("H0t.dat",std::ios::out);
// ofs_ht << "ntau " << ntau << '\n';
// for(int itau=0; itau<ntau+1; itau++)
// {
//   ofs_ht << std::setw(5) << itau << std::setw(12) << std::fixed << std::setprecision(3) << (beta/ntau)*itau;
//   for(int is=0; is<nspin; is++)
//   {
//     for(int m=0; m<m_tot; m++)
//     {
//       ofs_ht << std::setw(15) << std::fixed << std::setprecision(9) << hyb_tau[ineq][is][itau][m*m_tot+m].real();
//     }//m
//   }//is
//   ofs_ht << '\n';
// }//iomega
// ofs_ht.close();
    }//itau

    return;
  }

}
