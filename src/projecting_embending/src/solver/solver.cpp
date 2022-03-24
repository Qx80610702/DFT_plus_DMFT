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
  current_charge_step(args.current_charge_step),
  current_DMFT_step(args.current_DMFT_step),
  last_charge_step(args.last_charge_step),
  last_DMFT_step(args.last_DMFT_step),
  flag_eva_sigma_only(args.sigma_only),
  flag_eva_spectrum(args.cal_spectrum),
  flag_update_density(args.update_density )
  {;}

  void solver::solve()
  {
    debug::codestamp("solver::sovle");

    double time;
    double seconds;
    timer::timestamp(time);
    
    this->reading_inputs();
    
    //Check whwther dmft iteration convergency achieveed in last dmft step
    if( !this->flag_eva_spectrum ){
      this->flag_convergency = this->self_energy_scf_update();
      if(this->flag_eva_sigma_only) return;   //Only calculated the self enegy of last step
      if(!this->flag_update_density && this->flag_convergency && this->current_DMFT_step>1) return;
    }

    if(mpi_rank()==0){
      if(this->flag_eva_spectrum){
        std::cout << "\n<><><><><><><><><><><><><><><><><><><><><><><><><><><><><>\n";
        std::cout << "<><><><><><><><><> Calculating spectrum <><><><><><><><><>\n";
        std::cout << "<><><><><><><><><><><><><><><><><><><><><><><><><><><><><>" << std::endl; 
      }
      else if(this->flag_update_density){
        std::cout << "<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>\n";
        std::cout << "<><><><><><><><><>  Update charge density <><><><><><><><><>\n";
        std::cout << "<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>" << std::endl;
      }
      else{
          std::cout << "\n<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>\n";
          std::cout << "<><><><><><>  Current charge step" 
                    << std::setw(4) << this->current_charge_step 
                    << " DMFT step" << std::setw(4) << this->current_DMFT_step 
                    << " <><><><>\n"; 
          std::cout << "<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>" << std::endl;
      }
    }
    
    this->space.KS_bands_window(
          this->current_charge_step,
          this->pars.bands, 
          this->pars.atom,
          this->pars.in );

    // this->Mu.evaluate_mu_bisection_imag_DFT(
    //       this->pars.bands, this->pars.atom, 
    //       this->pars.in, this->space);
    
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
          this->pars.atom, this->pars.in,
          *(bool*)this->pars.in.parameter("hyb_func"),
          *(double*)this->pars.in.parameter("hyf_xc_alpha") );

    if(this->flag_eva_spectrum) this->cal_spectrum_func();
    else if(this->flag_update_density) this->charge_solve();
    else this->DMFT_solve();

    timer::get_time(time, seconds);
    if(mpi_rank()==0){
      if(this->flag_eva_spectrum)
        std::cout << "\nSpectrum evaluation time consuming (seconds): " << seconds  << std::endl;
      else if(this->flag_update_density)
        std::cout << "\nUpdating charge density time consuming (seconds): " << seconds  << std::endl;
      else
        std::cout << "\nProjecting and embending time consuming (seconds): " << seconds  << std::endl;
    }
    
    return;
  }

  void solver::DMFT_solve()
  {
    debug::codestamp("solver::DMFT_sovle");

    const int impurity_solver = *(int*)pars.in.parameter("impurity_solver");
    if(impurity_solver==1 || 
       impurity_solver==2 || 
       impurity_solver==3 ||
       impurity_solver==4 ||
       impurity_solver==5 )
       this->flag_axis = 0;       //imaginary axis
    else this->flag_axis = 1;     //real axis

    // this->flag_convergency = this->scf_update();
    // if(this->flag_eva_sigma_only || this->flag_convergency) return;  //Only calculated the self enegy of last step

    if(this->current_charge_step == 1 && this->current_DMFT_step==1)
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

    return;
  }

  void solver::cal_spectrum_func()
  {
    debug::codestamp("solver::cal_spectrum_func");

    this->flag_axis = 1;   //real axis

    this->imp.sigma.sigma_real.read_AC_sigma(
          this->pars.bands, 
          this->pars.in, 
          this->pars.atom );
    
    this->imp.sigma.subtract_double_counting(this->flag_axis);
    
    this->Mu.read_chemical_potential();

    // this->Aw.eva_spectrum(
    //       this->Mu.mu_corrected(),
    //       this->pars.bands, this->pars.atom, 
    //       this->proj, this->imp.sigma, 
    //       this->pars.in, this->space );

    this->Aw.eva_spectrum_normalization(
          this->Mu.mu_corrected(),
          this->pars.bands, this->pars.atom, 
          this->proj, this->imp.sigma, 
          this->pars.in, this->space );

    if(mpi_rank()==0) this->Aw.out_spectrum();

    return;
  }
  
  void solver::charge_solve()
  {
    debug::codestamp("solver::charge_solve");

    const int impurity_solver = *(int*)pars.in.parameter("impurity_solver");
    if(impurity_solver==1 || 
       impurity_solver==2 || 
       impurity_solver==3 ||
       impurity_solver==4 ||
       impurity_solver==5 )
       this->flag_axis = 0;       //imaginary axis
    else this->flag_axis = 1;     //real axis

    const double beta = *(double*)pars.in.parameter("beta");
    const int nomega = *(int*)pars.in.parameter("n_omega");

    if(this->flag_axis==0){       //imaginary axis
      this->imp.sigma.sigma_imag.nomega() = nomega;
      this->imp.sigma.sigma_imag.inverse_T() = beta;
      this->imp.sigma.sigma_imag.Matsubara_freq().resize(nomega);
      auto& freq=this->imp.sigma.sigma_imag.Matsubara_freq();
      for(int iomega=0; iomega<nomega; iomega++)
          freq[iomega] = (2*iomega+1)*PI/beta;
    }

    this->imp.read_last_step( 
            this->last_charge_step,
            this->last_DMFT_step,
            *(int*)pars.in.parameter("impurity_solver"),
            this->pars.bands, this->pars.in, this->pars.atom );

    this->imp.sigma.subtract_double_counting(this->flag_axis);

    this->Mu.update_chemical_potential(
            this->flag_axis, this->pars.bands, 
            this->pars.atom, this->proj, 
            this->imp.sigma, this->pars.in, this->space );

    this->Char_scf.update_char_dens(
            this->flag_axis, this->Mu,
            this->pars.bands, this->pars.atom, 
            this->proj, this->imp.sigma, 
            this->pars.in, this->space );

    this->Char_scf.output_char_dense(
            *(int*)pars.in.parameter("dft_solver"),
            this->pars.bands.nk() );

    this->Char_scf.prepare_nscf_dft(
            *(int*)pars.in.parameter("dft_solver"),
            *(int*)pars.in.parameter("max_dft_step") );
    
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

  bool solver::self_energy_scf_update()
  {
    debug::codestamp("solver::self_energy_scf_update");

    const double beta = *(double*)pars.in.parameter("beta");
    const int nomega = *(int*)pars.in.parameter("n_omega");

    bool convergency=false;

    if(this->flag_axis==0){
      this->imp.sigma.sigma_imag.nomega() = nomega;
      this->imp.sigma.sigma_imag.inverse_T() = beta;
      this->imp.sigma.sigma_imag.Matsubara_freq().resize(nomega);
      auto& freq=this->imp.sigma.sigma_imag.Matsubara_freq();
      for(int iomega=0; iomega<nomega; iomega++)
          freq[iomega] = (2*iomega+1)*PI/beta;
    }

    if(this->current_charge_step>1 || this->current_DMFT_step>1)
    {
      this->imp.read_last_step(
            this->last_charge_step,
            this->last_DMFT_step,
            *(int*)pars.in.parameter("impurity_solver"),
            this->pars.bands, this->pars.in, this->pars.atom );

      convergency = this->imp.scf_condition(
            this->flag_axis,this->pars.bands, 
            this->pars.atom, this->pars.in );

      if(mpi_rank()==0 && !this->flag_eva_sigma_only)
      {
        if(convergency)
          std::cout << "\nSelf-consistency of self-energy in charge step" 
                    << std::setw(4) << this->last_charge_step 
                    << "  dmft step" << std::setw(4) << this->last_DMFT_step 
                    << " : true\n";
        else
          std::cout << "\nSelf-consistency of self-energy in charge step" 
                    << std::setw(4) << this->last_charge_step 
                    << "  dmft step" << std::setw(4) << this->last_DMFT_step 
                    << " : false\n";
          
        std::cout << "    impuritys            Delta_Sigma\n";
        for(int ineq=0; ineq<this->pars.atom.inequ_atoms(); ineq++){
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
                this->last_DMFT_step, this->pars.bands, 
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

      if( *(int*)this->pars.in.parameter("impurity_solver")==1 ||
          *(int*)this->pars.in.parameter("impurity_solver")==2 ) 
        this->imp.evaluate_delta_Weiss_tau(
              *(int*)this->pars.in.parameter("impurity_solver"),
              this->pars.atom, this->pars.bands.nspins(),
              *(double*)this->pars.in.parameter("beta"),
              *(int*)this->pars.in.parameter("n_tau"),
              *(int*)this->pars.in.parameter("n_omega") );
    }

// timer::get_time(time,seconds);
// std::cout << "evaluate_hybridization_function_imag time consuming " << seconds << '\n';

    return;
  }

  void solver::output_to_impurity_solver()
  {
    this->imp.out(this->current_charge_step,
                  this->current_DMFT_step,
                  *(int*)this->pars.in.parameter("impurity_solver"),
                  this->Mu.mu_corrected(), this->pars.bands, 
                  this->pars.in, this->pars.atom, this->Umat );

    this->Umat.out_coulomb_tensor(
                  1, this->current_charge_step, 
                  this->current_DMFT_step,
                  *(int*)this->pars.in.parameter("impurity_solver"),
                  this->pars.atom, this->pars.bands );
  }

}
