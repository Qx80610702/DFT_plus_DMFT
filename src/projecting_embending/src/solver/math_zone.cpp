#include "math_zone.h"

void polynomial_regression(std::complex<double>* fw, double* freq, 
                    const int nomega, double& C1, double& C2, double& C3)
{
  //=========================polynomial regression========================
  //         Adapted from the CT-HYB solver of G.L. 
  //  Calculate the coefficients of high-frequency tail of 
  //  imaginary-frequency hybridization function and Green's function etc.
  //  f(iw_n) = C1/iw_n + C2/(iw_n)^2 + c3/(iw_n)^3
  //  f(tau) = -C1/2 + C2/4*(-beta + 2tau) + C3/4(beta tau - tau^2)
  //  NOTE: C1, C2, C3 are averaged over the last 20 Matsubara frequency points

  const std::complex<double> zero(0.0,0.0), im(0.0,1.0);

  std::complex<double> C[3][3];
  for(int i=0; i<3; i++)
    for(int j=0; j<3; j++)
      C[i][j] = zero;

  std::complex<double> D1(0.0,0.0), D2(0.0,0.0), D3(0.0,0.0);

  for(int iomega=0; iomega<nomega; iomega++)
  {
    C[0][0] += -1.0/std::pow(freq[iomega],2);
    C[0][1] += im/std::pow(freq[iomega],3);
    C[0][2] += 1.0/std::pow(freq[iomega],4);
    C[1][2] += -im/std::pow(freq[iomega],5);
    C[2][2] += -1.0/std::pow(freq[iomega],6);
    
    D1 += fw[iomega]/(im*freq[iomega]);
    D2 += -fw[iomega]/std::pow(freq[iomega],2);
    D3 += im*fw[iomega]/std::pow(freq[iomega],3);
  }
  
  C[1][0] = C[0][1];
  C[1][1] = C[0][2];
  C[2][0] = C[0][2];
  C[2][1] = C[1][2];

  int info_tri, info_trf;
  int ipiv[3];

  const int mkl_threads = mkl_get_max_threads();
  mkl_set_num_threads(1);  //set the number of threads of MKL library function to 1
  general_complex_matrix_inverse(&C[0][0], 3, ipiv);
  mkl_set_num_threads(mkl_threads);

  C1 = (C[0][0]*D1 + C[0][1]*D2 + C[0][2]*D3).real();
  C2 = (C[1][0]*D1 + C[1][1]*D2 + C[1][2]*D3).real();
  C3 = (C[2][0]*D1 + C[2][1]*D2 + C[2][2]*D3).real();

  return;
}