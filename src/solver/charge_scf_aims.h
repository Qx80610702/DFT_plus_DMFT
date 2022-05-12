#ifdef __FHIaims

#include <complex>
#include <vector>

namespace DFT_plus_DMFT
{
  class Charge_SCF_aims
  {
    public:
    Charge_SCF_aims(){;}
    ~Charge_SCF_aims(){;}

    public:
    void output_charge_density(
          const int nks, 
          std::vector<std::vector<std::vector<
          std::complex<double>>>>& dens_mat_cmplx);

    void read_charge_density(
          const int istep,
          const bool DMFT_charge);

    void read_charge_density_matrix(
          const int nks, 
          std::vector<std::vector<std::vector<
          std::complex<double>>>>& dens_mat_cmplx);

    void prepare_nscf_dft();

    private:
    std::vector<std::vector<std::vector<double>>> rho;     //rho[istep][ispin][igrid]   
    std::vector<double> partition_tab;                    //partition_tab][igrid]

    std::vector<int> istep2index;
  };
}

#endif