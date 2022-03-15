#pragma once

#include <complex>
#define MKL_Complex16 std::complex<double>
#include <mkl.h>

#include <assert.h>
#include <cmath>
#include <iostream>
#include <cstdlib>   //Use exit function

inline void general_complex_matrix_inverse(
              std::complex<double>* a, 
              const int dim, int* ipiv )
{
  lapack_int info=LAPACKE_zgetrf(LAPACK_ROW_MAJOR, dim, dim, a, dim, ipiv);
  if(info != 0 ){
    std::cout << "Error in LAPACKE_zgetrf!!!" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  info=LAPACKE_zgetri(LAPACK_ROW_MAJOR, dim, a, dim, ipiv);
  if(info != 0 ){
    std::cout << "Error in LAPACKE_zgetri!!!" << std::endl;
    std::exit(EXIT_FAILURE);
  }

}

inline void Hermitian_matrix_sqrt_inver(
              std::complex<double>* a,
              double* eigen_val,
              std::complex<double>* work_mat,
              std::complex<double>* b,
              const int dim,
              int& info_zheev)
{
  // lapack_int LAPACKE_zheev( int matrix_layout, char jobz, char uplo, lapack_int n,
  //     lapack_complex_double* a, lapack_int lda, double* w );

  info_zheev=LAPACKE_zheev(LAPACK_ROW_MAJOR, 'V', 'U', dim, a, dim, eigen_val);

  std::complex<double> alpha(1.0,0.0), beta(0.0,0.0);

  for(int index=0; index<dim*dim; index++) 
    b[index] = std::complex<double>(0.0,0.0);

  for(int i=0; i<dim; i++){
    if(std::fabs(eigen_val[i])<1.0e-12){
      std::cout << "Error in Hermitian_matrix_sqrt_inver : zero eigenvalues!!!" << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }

  for(int i=0; i<dim; i++)
    b[i*dim+i] = std::complex<double>(1.0/std::sqrt(eigen_val[i]),0.0);

  cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
              dim, dim, dim,
              &alpha,
              a, dim,
              b, dim,
              &beta,
              work_mat, dim);

  cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasConjTrans,
              dim, dim, dim,
              &alpha,
              work_mat, dim,
              a, dim,
              &beta,
              b, dim);

}

void polynomial_regression(std::complex<double>* fw, double* freq, 
            const int nomega, double& C1, double& C2, double& C3);

