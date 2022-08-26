#ifdef __FHIaims

#include <complex>
#include <vector>
#include <deque>

namespace DFT_plus_DMFT
{
  class Charge_SCF_aims
  {
    public:
    Charge_SCF_aims(){;}
    ~Charge_SCF_aims(){;}

    public:
    void output_charge_density_matrix(
          const int nks,
          std::vector<std::vector<std::vector<
          std::complex<double>>>>& dens_mat_cmplx);

    void read_charge_density(
          const bool initial_charge,
          const bool DMFT_charge );

    void read_charge_density_matrix(
          const int nks,
          std::vector<std::vector<std::vector<
          std::complex<double>>>>& dens_mat_cmplx);

    void update_data(const int mix_step, const int max_mixing_step);

    void update_alpha(const int mix_step, std::vector<double>& alpha);

    void mixing_density(
          const int mix_step, 
          const double mixing_param, 
          const int max_mixing_step,
          std::vector<double>& alpha,
          double& charge_change);

    void prepare_nscf_dft();

    private:
    std::deque<std::vector<std::vector<double>>> Rrho;        //Rrho[istep][ispin][igrid]; Rrho(istep)= rho_in(istep) - rho_out(istep)
    std::deque<std::vector<std::vector<double>>> dRrho;       //dRrho[istep][ispin][igrid]; dRrho(istep)= Rrho(istep+1) - Rrho(istep)
    std::deque<std::vector<std::vector<double>>> Opt_rho;     //Opt_rho[istep][ispin][igrid];
    std::vector<std::vector<double>> rho_out;                 //Output Kohn-Sham or DFT+DMFT rho: rho_out[ispin][igrid]   
    std::vector<double> partition_tab;                        //partition_tab][igrid]
    
  };
}

#endif