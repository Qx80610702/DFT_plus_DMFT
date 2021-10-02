#include "KS_eigenvectors.h"
#include "../debug/debug.h"

#define MKL_Complex16 std::complex<double>
#include <mkl.h>

#include <fstream>
#include <iostream>
#include <sstream>
#include <cstdlib>   //Use exit function
#include <iomanip>   //use setw and setprecison function

namespace DFT_output
{
  bool KS_eigenvectors::read(
      const std::string dir, const int ik,
      DFT_plus_DMFT::Hilbert_space& space)
  {
    debug::codestamp("KS_eigenvectors::read");

    const std::vector<std::vector<bool>>& correction=space.correction_flag();
    const std::vector<int>& wbands=space.Wbands();
    const std::vector<std::vector<int>>& index=space.ibands2wbands();

    std::stringstream ss;

    int ispin, iband, ibasis;
    double real, imag;

    this->i_kpoint = ik;

    ss << dir << "eigenvector" << ik << ".dat";

    std::ifstream ifs(ss.str().c_str(), std::ios::in);

    if (!ifs)
	  {
	  	std::cout << "Fail to oepn " << ss.str().c_str() << std::endl;
      std::exit(EXIT_FAILURE);
    }

    ifs.seekg(0);      //set the position at the beginning of the file
    ifs >> this->flag_SOC;
    ifs.ignore(150,'\n');

    if(flag_SOC) // SOC
    {

    }
    else //non_SOC
    {
      ifs >> this->nspins;
      ifs >> this->nbands;
      ifs >> this->nbasis;
      ifs.ignore(150, '\n');

      if(this->eigenvector.empty())
      {
        this->eigenvector.resize(this->nspins);
        for(ispin=0; ispin<nspins; ispin++)
          this->eigenvector[ispin].resize(this->nbasis*wbands[ispin]);
      }

      while(ifs.good())
      {
        ifs >> ispin;
        if(ifs.eof()) break;
        ifs >> iband;
        ifs >> ibasis;

        ifs >> real;
        ifs >> imag;

        if(correction[ispin][iband])
          this->eigenvector[ispin][ibasis*wbands[ispin]+index[ispin][iband]] = std::complex<double>(real,imag);

        ifs.ignore(150, '\n'); 

        if(ifs.eof()) break;  //Check whether end of file is reached 
      }
    } 
    ifs.close();

    return true;
  }
  
  void KS_eigenvectors::evalute_k_wave_c_mat(
    std::complex<double>* wave_c_mat, const int is,
    const int iband1, const int iband2)
  {
    // void cblas_zgemm (const CBLAS_LAYOUT Layout, const CBLAS_TRANSPOSE transa, const
          // CBLAS_TRANSPOSE transb, const MKL_INT m, const MKL_INT n, const MKL_INT k, const void
          // *alpha, const void *a, const MKL_INT lda, const void *b, const MKL_INT ldb, const void
          // *beta, void *c, const MKL_INT ldc);
    // debug::codestamp("KS_eigenvectors::evalute_k_wave_c_mat");

    std::complex<double> alpha(1.0,0.0);
    std::complex<double> beta(0.0,0.0);
    
    cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasConjTrans,
                this->nbasis, this->nbasis, 1, 
                &alpha,
                &this->eigenvector[is][iband1], 1,
                &this->eigenvector[is][iband2], 1,
                &beta,
                wave_c_mat, this->nbasis);

    return;
  }

}
