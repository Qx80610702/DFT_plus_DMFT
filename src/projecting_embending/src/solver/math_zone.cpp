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

void Simpson_Integral(
      const int mesh,
      const double* const func,
      const double* const rab,
      double &asum )
{
  /*     simpson's rule integration. On input:
  !      mesh = mhe number of grid points (should be odd)
  !      func(i)= function to be integrated
  !      rab(i) = r(i) * dr(i)/di * di
  !      For the logarithmic grid not including r=0 :
  !      r(i) = r_0*exp((i-1)*dx) ==> rab(i)=r(i)*dx
  !      For the logarithmic grid including r=0 :
  !      r(i) = a(exp((i-1)*dx)-1) ==> rab(i)=(r(i)+a)*dx
  !      Output in asum = \sum_i c_i f(i)*rab(i) = \int_0^\infty f(r) dr
  !      where c_i are alternativaly 2/3, 4/3 except c_1 = c_mesh = 1/3
  */
  //  simpson's rule integrator for function stored on the
  //  radial logarithmic mesh
  //  routine assumes that mesh is an odd number so run check

  if(mesh%2!=1){
    std::cout << "The number of mesh must be odd in Simpson integral" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  asum = 0.00;
  const size_t end = mesh-2;
  for( size_t i=1; i!=end; i+=2 )
  {
    const double f1 = func[i]*rab[i];
    asum += f1 + f1 + func[i+1]*rab[i+1];
  }
  const double f1 = func[mesh-2]*rab[mesh-2];
  asum += f1+f1;
  asum += asum;
  asum += func[0]*rab[0] + func[mesh-1]*rab[mesh-1];
  asum /= 3.0;
  return;
}// end subroutine simpson

void Simpson_Integral(
    const int mesh,
    const double* const func,
    const double dr,
    double &asum )
{
  /*     simpson's rule integration. On input:
  !      mesh = mhe number of grid points (should be odd)
  !      func(i)= function to be integrated
  !      rab(i) = r(i) * dr(i)/di * di
  !      For the logarithmic grid not including r=0 :
  !      r(i) = r_0*exp((i-1)*dx) ==> rab(i)=r(i)*dx
  !      For the logarithmic grid including r=0 :
  !      r(i) = a(exp((i-1)*dx)-1) ==> rab(i)=(r(i)+a)*dx
  !      Output in asum = \sum_i c_i f(i)*rab(i) = \int_0^\infty f(r) dr
  !      where c_i are alternativaly 2/3, 4/3 except c_1 = c_mesh = 1/3
  */
  //  simpson's rule integrator for function stored on the
  //  radial logarithmic mesh
  //  routine assumes that mesh is an odd number so run check

  if(mesh%2!=1){
    std::cout << "The number of mesh must be odd in Simpson integral" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  
  asum = 0.00;
  const size_t end = mesh-2;
  for(size_t i=1; i!=end; i+=2 )
  {
    const double f1 = func[i];
    asum += f1 + f1 + func[i+1];
  }
  const double f1 = func[mesh-2];
  asum += f1+f1;
  asum += asum;
  asum += func[0] + func[mesh-1];
  asum *= dr/3.0;
  return;
}// end subroutine simpson


void Simpson_Integral_0toall(
    const int mesh,
    const double * const func,
    const double * const rab,
    double * const asum )
{
  const double r2=1.00/2.00, r3=1.00/3.00;
  asum[0] = 0.00;
  double f3 = func [0] * rab [0];
  for( int i=1; i<mesh; i+=2)
  {
    const double f1 = f3;
    const double f2 = func[i] * rab[i] ;
    f3 = func[i+1] * rab[i+1] ;
    asum[i] = asum[i-1] + r2*( f1 + f2);
    if(i+1<mesh)
    {
        asum[i+1] = asum[i-1] + r3*( f1 + 4.00*f2 + f3 );
    }
  }
  return;
}

void Simpson_Integral_alltoinf(
    const int mesh,
    const double* const func,
    const double* const rab,
    double* const asum )
{
  Simpson_Integral_0toall( mesh, func, rab, asum );

  const double asum_all = asum[mesh-1];
  for (int i = 0;i < mesh; ++i)
  {
    asum[i] = asum_all - asum[i];
  }
  return;
}
