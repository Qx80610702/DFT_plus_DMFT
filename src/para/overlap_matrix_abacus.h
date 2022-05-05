#pragma once

#include "correlated_atoms.h"

#include <vector>
#include <string>
#include <complex>

namespace DFT_output
{
  namespace abacus
  {
    class overlap_matrix
    {
      public:
      overlap_matrix(){;}
      overlap_matrix(const std::string dir);
      ~overlap_matrix(){;}

      public:
      void read(const std::string dir);
      void evaluate_ovlp_k(const int ik, atoms_info& atom);

      //Interface
      std::vector<std::complex<double>>& ovlp_mat(){return this->ovlp_matrix_k;}
      std::vector<std::complex<double>>& local_ovlp(){return this->ovlp_localorb_k;}


      private:
      void get_last_number_in_line(std::string& line, int& val);

      private:
      std::vector<std::vector<double>> kvector;   //kvector[nkpoints][3]

      std::vector<std::vector<double>> R_vector;  //R_vector[ncells][3]
      int ncells;
      int nbasis;

      std::vector<std::vector<double>> ovlp_csr;   //ovlp_csr[ncells][N]
      std::vector<std::vector<int>> col_indices;   //col_indices[ncells][N]
      std::vector<std::vector<int>> row_indptr;    //row_indptr[ncells][nbasis+1]

      //ovlp_matrix_k[nbasis*n_DNFT_orb];
      std::vector<std::complex<double>> ovlp_matrix_k;

      //ovlp_localorb_k[n_DNFT_orb*n_DNFT_orb];
      std::vector<std::complex<double>> ovlp_localorb_k;

      //dimension:nabsis*nabsis; ovlp_matrix_k[nbasis*nbasis];
      // std::vector<std::complex<double>> ovlp_matrix_work;
    };
  }
}