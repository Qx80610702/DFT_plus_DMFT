#include "solver.h"
#include "../mpi_environment.h"
#include "../debug.h"
#include "../timer.h"
#include "projector.h"
#include "../constants.h"

#include <iostream>
#include <map>
#include <complex>
#include <vector>
#include <iomanip>

namespace DFT_plus_DMFT
{
  solver::solver(argument_lists& args):
  DMFT_iteration_step(args.global_step),
  flag_eva_sigma_only(args.sigma_only),
  flag_eva_spectrum(args.cal_spectrum)
  {;}

  void solver::solve()
  {
    debug::codestamp("solver::sovle");

    if(this->flag_eva_spectrum) this->cal_spectrum_func();
    else this->DMFT_solve();
    
    return;
  }

  void solver::DMFT_solve()
  {
    debug::codestamp("solver::DMFT_sovle");

    double time;
    double seconds;
    timer::timestamp(time);

    this->reading_inputs();

    const int impurity_solver = *(int*)pars.in.parameter("impurity_solver");
    if(impurity_solver==1 || 
       impurity_solver==2 || 
       impurity_solver==3 ||
       impurity_solver==4 ||
       impurity_solver==5 )
       this->flag_axis = 0;       //imaginary axis
    else this->flag_axis = 1;     //real axis

    this->flag_convergency = this->scf_update();
    if(this->flag_eva_sigma_only || this->flag_convergency) return;  //Only calculated the self enegy of last step

    if(mpi_rank()==0)
    {
      std::cout << "<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>\n";
      std::cout << "<><><><><><>  Current DMFT iteration step:" 
                << std::setw(4) << this->DMFT_iteration_step << "  <><><><><><>\n";
      std::cout << "<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>" << std::endl;  
    }

    this->space.KS_bands_window(
          this->pars.bands, 
          this->pars.atom,
          this->pars.in );

    this->Mu.evaluate_mu_bisection_imag_DFT(
          this->pars.bands, this->pars.atom, 
          this->pars.in, this->space);
    
    this->proj.elaluate_projector(
          *(int*)pars.in.parameter("DFT_solver"),
          pars.bands, this->space, this->pars.atom );
    
    this->imp.evaluate_local_occupation( 
          this->pars.in, this->pars.bands, 
          this->proj, this->pars.atom,
          this->space, this->Mu.mu_DFT(),
          *(int*)this->pars.in.parameter("n_omega"),
          *(int*)this->pars.in.parameter("magnetism") );

    this->imp.sigma.dc.cal_double_counting( 
          *(int*)this->pars.in.parameter("double_counting"), 
          this->pars.bands.soc(), this->pars.bands.nspins(), 
          this->pars.atom, this->pars.in);

    if(this->DMFT_iteration_step==1)
      this->imp.sigma.initial_guess( 
            this->flag_axis, 
            this->pars.bands.soc(), 
            this->pars.bands.nspins(), 
            this->pars.atom );
    
    this->imp.sigma.subtract_double_counting(this->flag_axis);

    this->Umat.update_coulomb_tensor(1, this->pars.atom);

    this->Mu.update_chemical_potential(
             this->flag_axis, this->pars.bands, 
             this->pars.atom, this->proj, 
             this->imp.sigma, this->pars.in, this->space );

    this->update_Anderson_impurities();

    //prepare the file needed by impurity solver
    if(mpi_rank()==0) this->output_to_impurity_solver(); 

    timer::get_time(time, seconds);
    if(mpi_rank()==0)
      std::cout << "\nProjecting and embending time consuming (seconds): " << seconds << "\n" << std::endl;

    return;
  }

  void solver::cal_spectrum_func()
  {
    debug::codestamp("solver::cal_spectrum_func");

    double time;
    double seconds;
    timer::timestamp(time);

    if(mpi_rank()==0)
    {
      std::cout << "<><><><><><><><><><><><><><><><><><><><><><><><><><><><><>\n";
      std::cout << "<><><><><><><><><> Calculating spectrum <><><><><><><><><>\n";
      std::cout << "<><><><><><><><><><><><><><><><><><><><><><><><><><><><><>" << std::endl;  
    }

    this->flag_axis = 1;   //real axis

    this->reading_inputs();

    this->imp.sigma.sigma_real.read_AC_sigma(
          this->pars.bands, 
          this->pars.in, 
          this->pars.atom );

    this->space.KS_bands_window(
          this->pars.bands, 
          this->pars.atom,
          this->pars.in );
    
    this->Mu.evaluate_mu_bisection_imag_DFT(
          this->pars.bands, this->pars.atom, 
          this->pars.in, this->space );
    
    this->proj.elaluate_projector(
          *(int*)pars.in.parameter("DFT_solver"),
          pars.bands, this->space, this->pars.atom );
    
    this->imp.evaluate_local_occupation( 
          this->pars.in, this->pars.bands, 
          this->proj, this->pars.atom,
          this->space, this->Mu.mu_DFT(),
          *(int*)this->pars.in.parameter("n_omega"),
          *(int*)this->pars.in.parameter("magnetism") );

    this->imp.sigma.dc.cal_double_counting( 
          *(int*)this->pars.in.parameter("double_counting"), 
          this->pars.bands.soc(), this->pars.bands.nspins(), 
          this->pars.atom, this->pars.in);
    
    this->imp.sigma.subtract_double_counting(this->flag_axis);
    
    this->Mu.read_chemical_potential();

    this->Aw.eva_spectrum(
          this->Mu.mu_corrected(),
          this->pars.bands, this->pars.atom, 
          this->proj, this->imp.sigma, 
          this->pars.in, this->space );

    if(mpi_rank()==0) this->Aw.out_spectrum();

    timer::get_time(time, seconds);
    if(mpi_rank()==0)
      std::cout << "\nSpectrum evaluation time consuming (seconds): " << seconds << "\n" << std::endl;

    return;
  }

  void solver::reading_inputs()
  {
    debug::codestamp("solver::reading_inputs");
    
    this->pars.in.read();

    this->pars.atom.read(this->pars.in);

    this->pars.bands.read();

    // this->pars.tetra.read();
  }

  bool solver::scf_update()
  {
    debug::codestamp("solver::scf_update");

    const double beta = *(double*)pars.in.parameter("beta");
    const int nomega = *(int*)pars.in.parameter("n_omega");

    bool convergency=false;

    if(this->flag_axis==0)
    {
      this->imp.sigma.sigma_imag.nomega() = nomega;
      this->imp.sigma.sigma_imag.inverse_T() = beta;
      this->imp.sigma.sigma_imag.Matsubara_freq().resize(nomega);
      auto& freq=this->imp.sigma.sigma_imag.Matsubara_freq();
      for(int iomega=0; iomega<nomega; iomega++)
          freq[iomega] = (2*iomega+1)*PI/beta;
    }

    if(this->DMFT_iteration_step>1)
    {
      this->imp.read_last_step( this->DMFT_iteration_step,
            *(int*)pars.in.parameter("impurity_solver"),
            this->pars.bands, this->pars.in, this->pars.atom );

      convergency = this->imp.scf_condition(
            *(int*)pars.in.parameter("impurity_solver"),
            this->pars.bands, this->pars.atom, this->pars.in );

      if(mpi_rank()==0 && !this->flag_eva_sigma_only)
      {
        if(convergency)
          std::cout << "\nDMFT self-consistency in DMFT loop " << DMFT_iteration_step-1 
                    << " : true\n";
        else
          std::cout << "\nDMFT self-consistency in DMFT loop " << DMFT_iteration_step-1 
                    << " : false\n";
          
        std::cout << "    impuritys            delta_scf\n";
        for(int ineq=0; ineq<this->pars.atom.inequ_atoms(); ineq++)
        {
          std::cout << "    impurity" << ineq << "          " 
                    << std::setiosflags(std::ios::scientific) 
                    << this->imp.delta_scf()[ineq] << std::endl;
        }
      }

      if(*(int*)pars.in.parameter("impurity_solver")==1) //ALPS-CTHYB do not output self-enrgy
      {
        this->imp.update_self_energy( this->pars.bands, this->pars.atom, this->pars.in );

        if(mpi_rank()==0)
          this->imp.ALPS_hyb.out_sigma_last_step( 
                this->DMFT_iteration_step, this->pars.bands, 
                this->pars.in, this->pars.atom,
                this->imp.sigma.sigma_imag.sigma_new_access());
      }
    }

    return convergency;
  }

  void solver::update_Anderson_impurities()
  {
    debug::codestamp("solver::update_Anderson_impurities"); 

    this->imp.evaluate_impurity_level(
            this->space, this->pars.bands, 
            this->proj, this->pars.atom );

    if(this->flag_axis==0)
    {
      this->imp.evaluate_Weiss_hybridization_imag(
            this->pars.bands, this->proj, 
            this->pars.atom, this->space, this->Mu.mu_corrected(),
            *(int*)this->pars.in.parameter("n_omega"),
            *(int*)this->pars.in.parameter("magnetism"));

      // if(*(int*)this->pars.in.parameter("impurity_solver")!=4) 
      this->imp.evaluate_delta_Weiss_tau(
            *(int*)this->pars.in.parameter("impurity_solver"),
            this->pars.atom, this->pars.bands.nspins(),
            *(double*)this->pars.in.parameter("beta"),
            *(int*)this->pars.in.parameter("n_tau"),
            *(int*)this->pars.in.parameter("n_omega"));
    }

// timer::get_time(time,seconds);
// std::cout << "evaluate_hybridization_function_imag time consuming " << seconds << '\n';

    return;
  }

  void solver::output_to_impurity_solver()
  {
    this->imp.out(this->DMFT_iteration_step,
                  *(int*)this->pars.in.parameter("impurity_solver"),
                  this->Mu.mu_corrected(), this->pars.bands, 
                  this->pars.in, this->pars.atom, this->Umat);

    this->Umat.out_coulomb_tensor(1, this->DMFT_iteration_step,
                  *(int*)this->pars.in.parameter("impurity_solver"),
                  this->pars.atom, this->pars.bands);
  }

}
