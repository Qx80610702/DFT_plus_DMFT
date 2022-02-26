#pragma once

#include <string>
#include <complex>
#include <vector>
#include "../solver/Hilbert_space.h"

namespace DFT_output
{
  class KS_eigenvectors
  {
    public:
    KS_eigenvectors(){};
    ~KS_eigenvectors(){};

    bool read_corr_subset(const std::string dir, const int ik,
              DFT_plus_DMFT::Hilbert_space& space,
              std::vector<std::vector<std::complex<double>>>& eigenvector );

    bool read_DMFT_occ_subset(const std::string dir, const int ik,
              DFT_plus_DMFT::Hilbert_space& space,
              std::vector<std::vector<std::complex<double>>>& eigenvector );

    void evalute_k_wave_c_mat(
      std::complex<double>* wave_c_mat, const int is,
      const int iband1, const int iband2,
      std::vector<std::vector<std::complex<double>>>& eigenvector);

    //=================
    //  interfaces
    //=================
    bool soc(){return flag_SOC;}
    int nspin(){return nspins;}
    int nband(){return nbands;}
    int basis_n(){return nbasis;}

    // std::vector<std::vector<std::complex<double>>>&
    // wave_c(){return eigenvector;}
    

    private:
    int nbands;                //Total number of bands of each spin
    int nbasis;                //Number of basis
    int i_kpoint;              //the kpoint index of the eigenvector
    bool flag_SOC;             //true for SOC, false for non_SOC

    //=============================
    //   non_SOC case
    //=============================
    int nspins;                   //Total number of spin


    //eigen_values[ispin][nbasis*nbands];
    // std::vector<std::vector<std::complex<double>>> eigenvector; 

  };
  
}
typedef DFT_output::KS_eigenvectors wave_function;
