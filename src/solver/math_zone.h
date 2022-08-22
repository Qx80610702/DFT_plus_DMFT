#pragma once

#include <complex>

void general_complex_matrix_inverse(
      std::complex<double>* a, 
      const int dim, int* ipiv );

void general_real_matrix_inverse(
      double* a, const int dim, int* ipiv );

void Hermitian_matrix_sqrt_inver(
      std::complex<double>* a, 
      double* eigen_val,
      std::complex<double>* work_mat, 
      std::complex<double>* b,
      const int dim, 
      int& info_zheev );

void polynomial_regression(
      std::complex<double>* fw, 
      double* freq, 
      const int nomega, 
      double& C1, 
      double& C2, 
      double& C3 );


//=================================
// Simpson integral got from ABACUS
// ================================
void Simpson_Integral(
      const int mesh,
      const double * const func,
      const double * const rab,
      double &asum );

void Simpson_Integral(
      const int mesh,
      const double * const func,
      const double dr,
      double &asum );

void Simpson_Integral_0toall(
      const int mesh,
      const double * const func,
      const double * const rab,
      double * const asum );

void Simpson_Integral_alltoinf(
      const int mesh,
      const double * const func,
      const double * const rab,
      double * const asum );
