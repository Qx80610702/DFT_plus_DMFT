#include "projector.h"
#include "../mpi_environment.h"
#include "../para/overlap_matrix.h"
#include "math_zone.h"
#include "../constants.h"
#include "../global_variables.h"
#include "../timer.h"
#include "../debug.h"

#include <mpi.h>
#include <omp.h> 
#include <memory>

#define MKL_Complex16 std::complex<double>
#include <mkl.h>

//test
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <cmath>

namespace DFT_plus_DMFT
{
  void projector::evaluate_projector(
        const int DFT_solver,
        DFT_output::KS_bands& band,
        DFT_plus_DMFT::Hilbert_space& space,
        DFT_output::atoms_info& atom )
  {
    debug::codestamp("projector::evaluate_projector");

    double time;
    double seconds;
    timer::timestamp(time);

    const std::vector<int>& wbands = space.Wbands();
    const std::vector<int>& norb_sub = atom.iatom_norb();
    const std::vector<std::vector<int>>& orb_index = atom.Im2iorb();
    const int natom=atom.total_atoms();
    const int nkpoints = band.nk();
    const int nspin = band.nspins();
    const auto& epsilon = space.eigen_val();
    const int norb=atom.norb();
    const std::complex<double> zero(0.0,0.0), one(1.0,0.0);

    int wbands_max=wbands[0];
    if(nspin==2)
      wbands_max = wbands[0] > wbands[1] ? wbands[0] : wbands[1];

    std::vector<int> k_map;
    for(int ik=0; ik<nkpoints; ik++)
      if(ik%mpi_ntasks() == mpi_rank()) k_map.push_back(ik);

    if(this->projector_mat.empty()){
      this->projector_mat.resize(k_map.size());
      for(int ik=0; ik<k_map.size(); ik++)
      {
        this->projector_mat[ik].resize(natom);

        for(int iatom=0; iatom<natom; iatom++)
        {
          const int m_tot=norb_sub[iatom];

          this->projector_mat[ik][iatom].resize(nspin);
          for(int is=0; is<nspin; is++)
            this->projector_mat[ik][iatom][is].resize(wbands[is]*m_tot);
        }
      }
    }

    if(this->TB_Hk.empty()){
      this->TB_Hk.resize(k_map.size());
      for(int ik=0; ik<k_map.size(); ik++)
      {
        this->TB_Hk[ik].resize(natom);

        for(int iatom=0; iatom<natom; iatom++)
        {
          const int m_tot=norb_sub[iatom];

          this->TB_Hk[ik][iatom].resize(nspin);
          for(int is=0; is<nspin; is++)
            this->TB_Hk[ik][iatom][is].resize(m_tot);
        }
      }
    }

    if(this->unitary_trans.empty()){
      this->unitary_trans.resize(k_map.size());
      for(int ik=0; ik<k_map.size(); ik++)
      {
        this->unitary_trans[ik].resize(natom);

        for(int iatom=0; iatom<natom; iatom++)
        {
          const int m_tot=norb_sub[iatom];

          this->unitary_trans[ik][iatom].resize(nspin);
          for(int is=0; is<nspin; is++)
            this->unitary_trans[ik][iatom][is].resize(m_tot*m_tot, zero);
        }
      }
    }
    else{
      for(auto& iter1 : this->unitary_trans)
        for(auto& iter2 : iter1)
          for(auto& iter3 : iter2)
           for(auto& iter4 : iter3)
            iter4 = zero;
    }

    DFT_output::overlap_matrix ovlp(DFT_solver, "dft/outputs_to_DMFT/overlap_matrix");
    wave_function wfc;

    std::vector<std::complex<double>> Lowdin(norb*norb);
    std::vector<std::complex<double>> proj_mat_tmp(wbands_max*norb);
    std::vector<std::complex<double>> mat_tmp(wbands_max*norb);
    std::vector<std::complex<double>> norm_orth(norb*norb);
    std::vector<std::complex<double>> sqrt_inver(norb*norb);
    std::vector<std::complex<double>> work_mat(norb*norb);
    std::vector<double> eigen_val(norb);
    int info_zheev;

    std::vector<std::vector<std::complex<double>>> eigenvec;

    for(int ik=0; ik<k_map.size(); ik++)
    {
      const int i_k_point = k_map[ik];

      ovlp.evaluate_ovlp_k(i_k_point, atom);        //caculate overlap matrix in k-space

      wfc.read_corr_subset("dft/outputs_to_DMFT/KS_eigenvector/", i_k_point, space, eigenvec);       //read KS-eigenvector

      this->n_basis=wfc.basis_n();
      
      const std::vector<std::complex<double>>& ovlp_mat = ovlp.overlap();

      std::vector<std::complex<double>>& ovlp_localorb = ovlp.overlap_localorb();
      Hermitian_matrix_sqrt_inver(
        &ovlp_localorb[0], &eigen_val[0], 
        &work_mat[0], &Lowdin[0], norb, info_zheev );

      for(int is=0; is<nspin; is++)
      {
        // void cblas_zgemm (const CBLAS_LAYOUT Layout, const CBLAS_TRANSPOSE transa, const
        // CBLAS_TRANSPOSE transb, const MKL_INT m, const MKL_INT n, const MKL_INT k, const void
        // *alpha, const void *a, const MKL_INT lda, const void *b, const MKL_INT ldb, const void
        // *beta, void *c, const MKL_INT ldc);

        //===============test the orthonormality of wavefunctions=====
        /*
        std::vector<std::complex<double>> wfc_norm(wbands[is]*wbands[is]);
        std::vector<std::complex<double>> wfc_ovlp_tmp(wbands[is]*this->n_basis);
        const std::vector<std::complex<double>>& ovlp_all = ovlp.ovlp_aims.ovlp_mat_work();
        cblas_zgemm(CblasRowMajor, CblasConjTrans, CblasNoTrans,
                    wbands[is], this->n_basis, this->n_basis,
                    &one,
                    &eigenvec[is][0], wbands[is],
                    &ovlp_all[0], this->n_basis,
                    &zero,
                    &wfc_ovlp_tmp[0], this->n_basis );

        cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                    wbands[is], wbands[is], this->n_basis,
                    &one,
                    &wfc_ovlp_tmp[0], this->n_basis,
                    &eigenvec[is][0], wbands[is],
                    &zero,
                    &wfc_norm[0], wbands[is] );
        
        std::stringstream ss;
        ss << "wfc-norm/wfc_norm_is" << is << "_ik" << ik << ".dat";
        std::ofstream ofs(ss.str().c_str(), std::ios::out);

        for(int m1=0; m1<wbands[is]; m1++)
        {
          for(int m2=0; m2<wbands[is]; m2++)
          {
            ofs << std::setw(15) << std::fixed << std::setprecision(9) << wfc_norm[m1*wbands[is]+m2].real()
                << std::setw(15) << std::fixed << std::setprecision(9) << wfc_norm[m1*wbands[is]+m2].imag(); 
          }
          ofs << std::endl;
        }//m
        ofs.close(); */


        cblas_zgemm(CblasRowMajor, CblasConjTrans, CblasNoTrans,
                    wbands[is], norb, this->n_basis,
                    &one,
                    &eigenvec[is][0], wbands[is],
                    &ovlp_mat[0], norb,
                    &zero,
                    &proj_mat_tmp[0], norb);

        cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                    wbands[is], norb, norb,
                    &one,
                    &proj_mat_tmp[0], norb,
                    &Lowdin[0], norb,
                    &zero,
                    &mat_tmp[0], norb);

        cblas_zgemm(CblasRowMajor, CblasConjTrans, CblasNoTrans,
                    norb, norb, wbands[is],
                    &one,
                    &mat_tmp[0], norb,
                    &mat_tmp[0], norb,
                    &zero,
                    &norm_orth[0], norb );

        Hermitian_matrix_sqrt_inver(
                    &norm_orth[0], &eigen_val[0], &work_mat[0],
                    &sqrt_inver[0], norb, info_zheev );

        cblas_zgemm( CblasRowMajor, CblasNoTrans, CblasNoTrans,
                    wbands[is], norb, norb,
                    &one,
                    &mat_tmp[0], norb,
                    &sqrt_inver[0], norb,
                    &zero,
                    &proj_mat_tmp[0], norb);

        for(int iatom=0; iatom<natom; iatom++)
        {
          const int m_tot = norb_sub[iatom];
          for(int m=0; m<m_tot; m++)
            for(int iband=0; iband<wbands[is]; iband++)
              this->projector_mat[ik][iatom][is][iband*m_tot+m] = 
                    proj_mat_tmp[iband*norb+orb_index[iatom][m]];
        
          //Local Haimiltonian
          for(int m1=0; m1<m_tot; m1++)
            for(int m2=0; m2<m_tot; m2++)
              for(int iband=0; iband<wbands[is]; iband++)
                this->unitary_trans[ik][iatom][is][m1*m_tot+m2] += 
                  std::conj(proj_mat_tmp[iband*norb+orb_index[iatom][m1]])*
                  epsilon[is][i_k_point][iband]*
                  proj_mat_tmp[iband*norb+orb_index[iatom][m2]];

          for(int m=0; m<m_tot; m++)
            this->TB_Hk[ik][iatom][is][m] = this->unitary_trans[ik][iatom][is][m*m_tot+m];
        }

      }//is
    }//ik

    timer::get_time(time, seconds);
    GLV::ofs_running << "Time consuming for contruction of the projector: " 
                << (int)seconds << "s" << std::endl;

    //=========TEST orthonormality========================
    /*
    for(int is=0; is<nspin; is++)
    {
      std::vector<std::complex<double>> norm(norb*norb);
      int ik_count=-1;
      for(int ik=0; ik<nkpoints; ik++)
      {
        if(ik%mpi_ntasks() != mpi_rank()) continue;  //k_points are divided acording to process id
        ik_count++;
        for(int iatom1=0; iatom1<natom; iatom1++)
        {
          const int m_tot1 = norb_sub[iatom1];
          for(int m1=0; m1<m_tot1; m1++)
          {
            const int iorb1 = orb_index[iatom1][m1];
            for(int iatom2=0; iatom2<natom; iatom2++)
            {
              const int m_tot2 = norb_sub[iatom2];
              for(int m2=0; m2<m_tot2; m2++)
              {
                const int iorb2 = orb_index[iatom2][m2];
                norm[iorb1*norb+iorb2] = zero;
                for(int iband=0; iband<wbands[is]; iband++)
                  norm[iorb1*norb+iorb2] += std::conj(this->projector_mat[ik_count][iatom1][is][iband*m_tot1+m1])*
                                       this->projector_mat[ik_count][iatom2][is][iband*m_tot2+m2];
              }
            }//iatom2
          }//m1
        }//iatom1

        std::stringstream ss;
        ss << "projector_orthonormality/norm_ik" << ik << ".dat";
        std::ofstream ofs(ss.str().c_str(), std::ios::out);

        for(int iatom1=0; iatom1<natom; iatom1++)
        {
          const int m_tot1 = norb_sub[iatom1];
          for(int m1=0; m1<m_tot1; m1++)
          {
            const int iorb1 = orb_index[iatom1][m1];
            for(int iatom2=0; iatom2<natom; iatom2++)
            {
              const int m_tot2 = norb_sub[iatom2];
              for(int m2=0; m2<m_tot2; m2++)
              {
                const int iorb2 = orb_index[iatom2][m2];
                ofs << is << std::setw(5) << ik << std::setw(3) << iorb1 << std::setw(3) << iorb2 <<
                std::setw(15) << std::fixed << std::setprecision(9) << norm[iorb1*norb+iorb2].real()
                << std::setw(15) << std::fixed << std::setprecision(9) << norm[iorb1*norb+iorb2].imag() << std::endl;
              }
            }//iatom2
          }//m1
        }//iatom1
        ofs.close();
      }//ik
    }//is 
    */

    //=============================================
    //   unitary tranform of local orbital to
    //   get diagonal Green function, hybridization
    //   function and etc.
    //   Tight binding Hamiltonian.
    //=============================================
    if(atom.local_sym()==0 || atom.local_sym()==1 || atom.local_sym()==2 || atom.local_sym()==3) return;
    auto projector_mat_tmp = this->projector_mat;

    for(int ineq=0; ineq<atom.inequ_atoms(); ineq++)
    {
      const int iatom = atom.ineq_iatom(ineq);
      const int m_tot = norb_sub[iatom];
      int info_zheev;
      std::vector<double> eigen_val(m_tot);

      for(int ik=0; ik<k_map.size(); ik++)
      {
        for(int is=0; is<nspin; is++)
        {
          auto& unitary_trans_isk = this->unitary_trans[ik][iatom][is];

          info_zheev=LAPACKE_zheev(LAPACK_ROW_MAJOR, 'V', 'U', m_tot,
                        &unitary_trans_isk[0], m_tot, &eigen_val[0] );

          for(int m=0; m<m_tot; m++) this->TB_Hk[ik][iatom][is][m] = std::complex<double>(eigen_val[m], 0.0);

          cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                      wbands[is], m_tot, m_tot,
                      &one,
                      &projector_mat_tmp[ik][iatom][is][0], m_tot,
                      &unitary_trans_isk[0], m_tot,
                      &zero,
                      &this->projector_mat[ik][iatom][is][0], m_tot );

          for(int iatom1=0; iatom1<natom; iatom1++)
            if(iatom1!=iatom && atom.equ_atom(iatom1)==ineq)
              cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                      wbands[is], m_tot, m_tot,
                      &one,
                      &projector_mat_tmp[ik][iatom1][is][0], m_tot,
                      &unitary_trans_isk[0], m_tot,
                      &zero,
                      &this->projector_mat[ik][iatom1][is][0], m_tot );
            
        }//is
      }//ik  
    }//iatom

//=========TEST orthonormality========================
// for(int is=0; is<nspin; is++)
// {
//   std::vector<std::complex<double>> norm(norb*norb);
//   int ik_count=-1;
//   for(int ik=0; ik<nkpoints; ik++)
//   {
//     if(ik%mpi_ntasks() != mpi_rank()) continue;  //k_points are divided acording to process id
//     ik_count++;
//     for(int iatom1=0; iatom1<natom; iatom1++)
//     {
//       const int m_tot1 = norb_sub[iatom1];
//       for(int m1=0; m1<m_tot1; m1++)
//       {
//         const int iorb1 = orb_index[iatom1][m1];
//         for(int iatom2=0; iatom2<natom; iatom2++)
//         {
//           const int m_tot2 = norb_sub[iatom2];
//           for(int m2=0; m2<m_tot2; m2++)
//           {
//             const int iorb2 = orb_index[iatom2][m2];
//             norm[iorb1*norb+iorb2] = zero;
//             for(int iband=0; iband<wbands[is]; iband++)
//               norm[iorb1*norb+iorb2] += std::conj(this->projector_mat[ik_count][iatom1][is][iband*m_tot1+m1])*
//                                    this->projector_mat[ik_count][iatom2][is][iband*m_tot2+m2];
//           }
//         }//iatom2
//       }//m1
//     }//iatom1

//     std::stringstream ss;
//     ss << "projector_orthonormality/norm_ik" << ik << ".dat";
//     std::ofstream ofs(ss.str().c_str(), std::ios::out);

//     for(int iatom1=0; iatom1<natom; iatom1++)
//     {
//       const int m_tot1 = norb_sub[iatom1];
//       for(int m1=0; m1<m_tot1; m1++)
//       {
//         const int iorb1 = orb_index[iatom1][m1];
//         for(int iatom2=0; iatom2<natom; iatom2++)
//         {
//           const int m_tot2 = norb_sub[iatom2];
//           for(int m2=0; m2<m_tot2; m2++)
//           {
//             const int iorb2 = orb_index[iatom2][m2];
//             ofs << is << std::setw(5) << ik << std::setw(3) << iorb1 << std::setw(3) << iorb2 <<
//             std::setw(15) << std::fixed << std::setprecision(9) << norm[iorb1*norb+iorb2].real()
//             << std::setw(15) << std::fixed << std::setprecision(9) << norm[iorb1*norb+iorb2].imag() << std::endl;
//           }
//         }//iatom2
//       }//m1
//     }//iatom1
//     ofs.close();
//   }//ik
// }//is


// //=========TEST completeness========================
// for(int is=0; is<nspin; is++)
// {
//   std::vector<std::complex<double>> norm(wbands[is]*wbands[is]);
//   int ik_count=-1;
//   for(int ik=0; ik<nkpoints; ik++)
//   {
//     if(ik%mpi_ntasks() != mpi_rank()) continue;  //k_points are divided acording to process id
//     ik_count++;

//     for(int iband1=0; iband1<wbands[is]; iband1++)
//     {
//       for(int iband2=0; iband2<wbands[is]; iband2++)
//       {
//         norm[iband1*wbands[is]+iband2] = zero;
//         for(int iorb=0; iorb<norb; iorb++)
//           norm[iband1*wbands[is]+iband2] += std::conj(this->projector_mat[ik_count][iband2*norb+iorb])*
//                                this->projector_mat[ik_count][is][iband1*norb+iorb];
//       }
//     }//m

//     std::stringstream ss;
//     ss << "projector_completeness/norm_ik" << ik << ".dat";
//     std::ofstream ofs(ss.str().c_str(), std::ios::out);

//     for(int iband1=0; iband1<wbands[is]; iband1++)
//     {
//       for(int iband2=0; iband2<wbands[is]; iband2++)
//       {
//         ofs << is << std::setw(5) << ik << std::setw(3) << iband1 << std::setw(3) << iband2 <<
//         std::setw(15) << std::fixed << std::setprecision(9) << norm[iband1*wbands[is]+iband2].real()
//         << std::setw(15) << std::fixed << std::setprecision(9) << norm[iband1*wbands[is]+iband2].imag() << std::endl;
//       }
//     }
//     ofs.close();
//   }//ik
// }//is

// DFT_output::overlap_matrix ovlp1(DFT_solver,"../DFT/outputs/overlap_matrix");
// wave_function wfc1;

// int ik_count_tmp=-1;
// for(int ik=0; ik<nkpoints; ik++)
// {
//   if(ik%mpi.ntasks() != mpi.rank()) continue;  //k_points are divided acording to process id
//   ik_count_tmp++;
  
//   ovlp1.evaluate_ovlp_k(ik, atom);   //caculate overlap matrix in k-space
//   this->evalute_projector_k(ovlp1, orth, wfc1, atom, ik, ik_count_tmp);   
// }

    return;
  }

//   void projector::evalute_projector_k(
//         DFT_output::overlap_matrix& ovlp,
//         DFT_output::KS_eigenvectors& wfc,
//         DFT_output::atoms_info& atom,
//         const int ik,
//         const int ik_count)
//   {
//     debug::codestamp("DFT_plus_DMFT::evalute_projector_k");

//     const std::vector<std::vector<std::complex<double>>>& 
//           eigenvec=wfc.wave_c();

// const int nbands=wfc.nband();
// const int NBANDS_square=nbands*nbands;
// const int nspin_tmp=wfc.nspin();
// const int this->n_basis = wfc.basis_n();
// const auto& ovlp_mat_all=ovlp.ovlp_aims.ovlp_mat_work();
// const std::complex<double> zero(0.0,0.0);
// std::vector<std::vector<std::complex<double>>> norm_tmp;
// norm_tmp.resize(nspin_tmp);
// for(int is=0; is<nspin_tmp; is++)
// {
//   const std::complex<double> alpha(1.0,0.0);
//   const std::complex<double> beta(0.0,0.0);
//   norm_tmp[is].resize(NBANDS_square);

//   std::unique_ptr<std::complex<double>[]> product_mat(new std::complex<double> [nbands*this->n_basis]);

//   cblas_zgemm(CblasRowMajor, CblasConjTrans, CblasNoTrans,
//             nbands, this->n_basis, this->n_basis,
//             &alpha,
//             &eigenvec[is][0], nbands,
//             &ovlp_mat_all[0], this->n_basis,
//             &beta,
//             &product_mat[0], this->n_basis);
  
//   cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
//             nbands, nbands, this->n_basis,
//             &alpha,
//             &product_mat[0], this->n_basis,
//             &eigenvec[is][0], nbands,
//             &beta,
//             &norm_tmp[is][0], nbands);
  
// }
//   std::stringstream ss;
//   ss << "wavefuction_orthomormality/norm_ik" << ik << ".dat";
//   std::ofstream ofs(ss.str().c_str(), std::ios::out);
//   for(int is=0; is<nspin_tmp; is++)
//   {
//     for(int band_index=0; band_index<NBANDS_square; band_index++)
//     {
//       const int iband1=band_index/nbands;
//       const int iband2=band_index%nbands;
//       ofs << std::setw(3) << is << std::setw(5) << iband1 << std::setw(5) << iband2
//       << std::setw(12) << std::fixed << std::setprecision(6) << norm_tmp[is][band_index].real() 
//       << std::setw(12) << std::fixed << std::setprecision(6) << norm_tmp[is][band_index].imag() << '\n';
//     }
//   }
//   ofs.close();

//     return;
//   }


}
