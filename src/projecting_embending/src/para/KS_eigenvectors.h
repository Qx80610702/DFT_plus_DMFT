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

    // void initialize();
    bool read(const std::string dir, const int ik,
              DFT_plus_DMFT::Hilbert_space& space);

    void evalute_k_wave_c_mat(
      std::complex<double>* wave_c_mat, const int is,
      const int iband1, const int iband2);

    //=================
    //  interfaces
    //=================
    bool soc(){return flag_SOC;}
    int nspin(){return nspins;}
    int nband(){return nbands;}
    inline int basis_n(){return nbasis;}

    std::vector<std::vector<std::complex<double>>>&
    wave_c(){return eigenvector;}
    

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
    std::vector<std::vector<std::complex<double>>> eigenvector; 

  };
  
}
typedef DFT_output::KS_eigenvectors wave_function;
