#pragma once

#include "correlated_atoms.h"

#include <string>
#include <complex>
#include <vector>

namespace DFT_output
{
  namespace aims
  {
    class overlap_matrix
    {
      public:
      overlap_matrix(){};
      ~overlap_matrix(){};

      bool read(const std::string dir);

      void evaluate_ovlp_k(const int ik, atoms_info& at_info);

      void out();    //test whether reading worked corrected

      //=================
      //  interfaces
      //=================
      int ncells();
      int kpoints();
      std::complex<double> kphase(int icell, int ik);
      int index_ham(int rc, int icell, int i_basis);
      int index_col_ham(int index);
      double ovlp_mat_real(int index);
      int nbasis(){return n_basis;}
      std::vector<std::complex<double>>& ovlp_mat(){return this->ovlp_matrix_k;}
      std::vector<std::complex<double>>& local_ovlp(){return this->ovlp_localorb_k;}
      std::vector<std::complex<double>>& ovlp_mat_work(){return this->ovlp_matrix_work;}

      private:
      int n_cells_in_hamiltonian;
      int n_k_points; 
      int n_basis;
      long size_colum_index_hamiltonian;
      long size_overlap_matrix;

      std::vector<std::vector<std::complex<double>>> k_phase;     //k_phase[icell][ik]
      std::vector<std::vector<std::vector<int>>> index_hamiltonian;  //index_hamiltonian[1 or 2][icell][i_basis]
      std::vector<int> colum_index_hamiltonian; 
      std::vector<double> ovlp_matrix_real; 

      //ovlp_matrix_k[nbasis*n_DNFT_orb];
      std::vector<std::complex<double>> ovlp_matrix_k;

      //ovlp_localorb_k[n_DNFT_orb*n_DNFT_orb];
      std::vector<std::complex<double>> ovlp_localorb_k;

      //dimension:nabsis*nabsis; ovlp_matrix_k[nbasis*nbasis];
      std::vector<std::complex<double>> ovlp_matrix_work;
    };
  }
}
