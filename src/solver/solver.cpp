#include "solver.h"
#include "../mpi_environment.h"
#include "../debug.h"
#include "../timer.h"
#include "projector.h"
#include "../constants.h"
#include "../global_variables.h"

#include <mpi.h>
#include <iostream>
#include <map>
#include <complex>
#include <vector>
#include <iomanip>
#include <unistd.h>

extern "C" void aims_(int* mpi_comm, int* unit, bool* mpi_switch);

namespace DFT_plus_DMFT
{
  void solver::solve()
  {
    debug::codestamp("solver::sovle");
    
    this->set_ios(GlobalV::ofs_running, GlobalV::ofs_error);

    GlobalV::ofs_running << "Welcome to DFT+DMFT calculation" << std::endl;
    GlobalV::ofs_running << "Number of processes of the job: " << mpi_ntasks() << "\n" << std::endl;

    this->reading_inputs();

    this->working_init();

    switch(this->calculation_type)
    {
    case 0: //DFT+DMFT scf
      return this->DFT_DMFT_scf();
      break;
    case 1: //specctra
      return this->cal_spectrum_func();
      break;
    default:
      GlobalV::ofs_error << "Unknown calculation type" << std::endl;
      std::exit(EXIT_FAILURE);
    }

    return;
  }

  void solver::DFT_DMFT_scf()
  {
    debug::codestamp("solver::DFT_DMFT_scf");
    
    double time, seconds;
    int days, hours, minutes;
    std::string start_date, end_date;

    timer::timestamp(time);
    timer::get_date_time(start_date);

    int loop_count = 1, charge_loop_count = 1;
    bool charge_convergency = false;
    bool sigma_convergency = false;
    int last_char_step = *(int*)this->pars.in.parameter("last_charge_step");
    int last_DMFT_step = *(int*)this->pars.in.parameter("last_dmft_step");

    //===============================================
    //           Charge loop
    //===============================================
    for(int char_step = *(int*)this->pars.in.parameter("start_charge_step"); 
        char_step <= *(int*)this->pars.in.parameter("max_charge_step") && 
        (!charge_convergency || !sigma_convergency); char_step++)
    {
      GlobalV::ofs_running << "\n================ Begin DFT+DMFT scf loop ================" << std::endl;
      GlobalV::ofs_running << "<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>\n";
      GlobalV::ofs_running << "<><><><><><><><><> Charge step " << char_step << " <><><><><><><><><><><><><>\n"; 
      GlobalV::ofs_running << "<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>" << std::endl;

      this->space.KS_bands_window(
          char_step,
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
            this->pars.atom, this->pars.in,
            *(bool*)this->pars.in.parameter("hyb_func"),
            *(double*)this->pars.in.parameter("hyf_xc_alpha") );

      //=========================================
      //            DMFT loop
      //=========================================
      for(int dmft_step = *(int*)this->pars.in.parameter("start_dmft_step");
          dmft_step <= *(int*)this->pars.in.parameter("max_dmft_step") &&
          !sigma_convergency; dmft_step++)
      {
        GlobalV::ofs_running << "\n<><><><><><><><><><><>  DMFT step "  << dmft_step << " <><><><><><><><><><><>\n";

        if(dmft_step>1) last_char_step = char_step;

        if(char_step==1 && dmft_step==1)
          this->imp.sigma.initial_guess( 
                this->flag_axis,
                this->pars.bands.soc(), 
                this->pars.bands.nspins(), 
                this->pars.atom );
        else
          if(*(bool*)this->pars.in.parameter("restart") && loop_count==1)
            this->imp.read_self_energy(
              true, last_char_step, last_DMFT_step,
              *(int*)pars.in.parameter("impurity_solver"),
              this->pars.bands, this->pars.in, this->pars.atom );

        this->projecting_embeding(char_step, dmft_step);

        this->impurity_solving(
            *(int*)this->pars.in.parameter("impurity_solver"),
              char_step, dmft_step, this->pars.atom );

        sigma_convergency = this->sigma_scf_update(char_step, dmft_step );

        last_DMFT_step = dmft_step;

        loop_count++;
      }//dmft loop

      if( *(int*)this->pars.in.parameter("max_charge_step") > 1 && 
          char_step < *(int*)this->pars.in.parameter("max_charge_step") )
      {
        this->DMFT_charge_updating();

        //=========================================
        //            DFT loop
        //=========================================
        for(int dft_step=1; dft_step <= *(int*)this->pars.in.parameter("max_dft_step"); dft_step++)
        {
          this->Char_scf.mix_char_dense(
              *(double*)this->pars.in.parameter("charge_mix_beta") );

          this->Char_scf.output_char_dense(
            *(int*)pars.in.parameter("dft_solver"),
            this->pars.bands.nk() );

          if(charge_loop_count==1){
            this->Char_scf.prepare_nscf_dft(
                  *(int*)pars.in.parameter("dft_solver"),
                  *(int*)pars.in.parameter("max_dft_step") );
          }

          this->run_nscf_dft(
              *(int*)this->pars.in.parameter("dft_solver") );
        }//DFT loop
      }

      charge_loop_count++;
    }//charge loop

    timer::get_time(time, seconds, minutes, hours, days);
    timer::get_date_time(end_date);

    GlobalV::ofs_running << "\nCongratulations!!! The DFT+DMFT calculation finished" << std::endl;
    GlobalV::ofs_running << "========================== The full time records ==========================" << std::endl;
    GlobalV::ofs_running << "       starting time              ending time              time-consumption" << std::endl;
    GlobalV::ofs_running << "      " << start_date
                         << "       " << end_date 
                         << "           " 
                         << days << "d "
                         << hours << "h " 
                         << minutes << "m "
                         << (int)seconds << "s" << std::endl;

    return;
  }
  
  void solver::projecting_embeding(
        const int charge_step, 
        const int DMFT_step)
  {
    debug::codestamp("solver::projecting_embeding");

    this->imp.sigma.subtract_double_counting(this->flag_axis);

    this->Umat.update_coulomb_tensor(1, this->pars.atom);

    this->Mu.update_chemical_potential(
             this->flag_axis, this->pars.bands, 
             this->pars.atom, this->proj, 
             this->imp.sigma, this->pars.in, this->space );

    this->update_Anderson_impurities();

    //prepare the files needed by the impurity solvers
    if(mpi_rank()==0) this->output_to_impurity_solver(charge_step, DMFT_step); 
    MPI_Barrier(MPI_COMM_WORLD);  //Blocks until all processes reach here

    /*
    this->space.KS_bands_window(
          this->current_charge_step,
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
          this->pars.atom, this->pars.in,
          *(bool*)this->pars.in.parameter("hyb_func"),
          *(double*)this->pars.in.parameter("hyf_xc_alpha") );

    if(this->flag_eva_spectrum) this->cal_spectrum_func();
    else if(this->flag_update_density) this->charge_solve();
    else this->DMFT_solve();
    */
    
    return;
  }

  void solver::impurity_solving(
        const int impurity_solver,
        const int char_step,
        const int DMFT_step,
        DFT_output::atoms_info& atom )
  {
    debug::codestamp("solver::impurity_solving");

    GlobalV::ofs_running << "\nImpurity solver starts working......" << std::endl;
    GlobalV::ofs_running << "  Impurities       starting time              ending time              time-consumption" << std::endl;

    switch(impurity_solver)
    {
      case 1:
        // this->ALPS_hyb.solving(char_step, DMFT_step, 
        //        mu, in, atom, band, this->impurity_level, 
        //        this->sigma.sigma_new(0), this->Weiss_omega, this->hyb_omega);
        break;
      case 2:
        // this->ALPS_hyb_segment.output(char_step, DMFT_step,
        //        mu, in, atom, band, this->impurity_level, 
        //        this->sigma.sigma_new(0), this->Weiss_omega);
        break;
      case 3:
        this->imp.pacs.impurities_solving(char_step, DMFT_step, atom);
        break;
      case 4:
        this->imp.Rutgers.impurities_solving(char_step, DMFT_step, atom);
        break;
      case 5: 
        this->imp.iQIST_narcissus.impurities_solving(char_step, DMFT_step, atom);
        break;
      default:
        GlobalV::ofs_error << "Not supported impurity solver" << std::endl;
        std::exit(EXIT_FAILURE);
    }

    return;
  }

  void solver::working_init()
  {
    debug::codestamp("solver::working_init");

    this->calculation_type = *(int*)this->pars.in.parameter("calculation");

    if(this->calculation_type == 1) this->flag_axis = 1;   //spectra
    else if(this->calculation_type == 0){                  //DFT+DMFT scf
      const int impurity_solver = *(int*)this->pars.in.parameter("impurity_solver");
      const double beta = *(double*)pars.in.parameter("beta");
      const int nomega = *(int*)pars.in.parameter("n_omega");

      if(impurity_solver==1 || 
         impurity_solver==2 || 
         impurity_solver==3 ||
         impurity_solver==4 ||
         impurity_solver==5 )
         this->flag_axis = 0;       //imaginary axis
      else this->flag_axis = 1;     //real axis

      if(this->flag_axis==0){
        this->imp.sigma.sigma_imag.nomega() = nomega;
        this->imp.sigma.sigma_imag.inverse_T() = beta;
        this->imp.sigma.sigma_imag.Matsubara_freq().resize(nomega);
        auto& freq=this->imp.sigma.sigma_imag.Matsubara_freq();
        for(int iomega=0; iomega<nomega; iomega++)
            freq[iomega] = (2*iomega+1)*GlobalC::PI/beta;
      }
    }

    return;
  }

  void solver::cal_spectrum_func()
  {
    debug::codestamp("solver::cal_spectrum_func");

    GlobalV::ofs_running << "\n================ Begin calculating DFT+DMFT spectra ================" << std::endl;

    this->space.KS_bands_window(
          2, this->pars.bands, 
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
          this->pars.atom, this->pars.in,
          *(bool*)this->pars.in.parameter("hyb_func"),
          *(double*)this->pars.in.parameter("hyf_xc_alpha") );

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

    GlobalV::ofs_running << "DFT+DMFT spectra calculation stop successfuly!!!" << std::endl;

    return;
  }
  
  void solver::DMFT_charge_updating()
  {
    debug::codestamp("solver::DMFT_charge_updating");

    /*
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
          freq[iomega] = (2*iomega+1)*GlobalC::PI/beta;
    }

    this->imp.read_last_step( 
            this->last_charge_step,
            this->last_DMFT_step,
            *(int*)pars.in.parameter("impurity_solver"),
            this->pars.bands, this->pars.in, this->pars.atom );
    */

    GlobalV::ofs_running << "\nStarting update DMFT corrected charge density" << std::endl;

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

  bool solver::sigma_scf_update(
        const int charge_step, 
        const int DMFT_step)
  {
    debug::codestamp("solver::sigma_scf_update");

    bool convergency = false;

    this->imp.read_self_energy(
          false, charge_step, DMFT_step,
          *(int*)pars.in.parameter("impurity_solver"),
          this->pars.bands, this->pars.in, this->pars.atom );

    convergency = this->imp.scf_condition(
          this->flag_axis, this->pars.bands, 
          this->pars.atom, this->pars.in );

    GlobalV::ofs_running << "\nSelf-consistency of self-energy in charge step " 
              << charge_step 
              << " DMFT step " << DMFT_step 
              << " : false\n";
      
    GlobalV::ofs_running << "    impuritys            Delta_Sigma (eV)\n";
    for(int ineq=0; ineq<this->pars.atom.inequ_atoms(); ineq++){
      GlobalV::ofs_running << "    impurity" << ineq << "              " 
                << std::setprecision(3) << std::setiosflags(std::ios::scientific)
                << this->imp.delta_scf()[ineq] << std::endl;
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
// GlobalV::ofs_running << "evaluate_hybridization_function_imag time consuming " << seconds << '\n';

    return;
  }

  void solver::output_to_impurity_solver(
        const int charge_step, 
        const int DMFT_step )
  {
    this->imp.out(charge_step, DMFT_step,
            *(int*)this->pars.in.parameter("impurity_solver"),
            this->Mu.mu_corrected(), this->pars.bands, 
            this->pars.in, this->pars.atom, this->Umat );

    this->Umat.out_coulomb_tensor(
            1, charge_step, DMFT_step,
            *(int*)this->pars.in.parameter("impurity_solver"),
            this->pars.atom, this->pars.bands );
  }

  void solver::set_ios(
    std::ofstream& ofs_running, 
    std::ofstream& ofs_error )
  {
    if(mpi_rank()==0){
      ofs_running.open("DMFT_running.log", std::ios_base::app);
      ofs_error.open("DMFT_running.error", std::ios_base::out);
    }
  }

  void solver::run_nscf_dft(
      const int dft_solver )
  {
    debug::codestamp("solver::run_nscf_dft");

    double time, seconds;
    int minutes;
    timer::timestamp(time);

    GlobalV::ofs_running << "\nStart non-self-consistent DFT..." << std::endl;

    int ierr = chdir("./dft");
    if(ierr != 0){
      std::cout << "Process " << mpi_rank() << " fails to enter to the directory dft" << std::endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);  //Blocks until all processes reach here

    switch(dft_solver)
    {
      case 1: //FHI-aims
        #ifdef __FHIaims
          int unit = 6;
          bool use_mpi = true;
          int mpi_comm_global = MPI_COMM_WORLD;

          aims_(&mpi_comm_global, &unit, &use_mpi);    //run FHI-aims
        #else
          GlobalV::ofs_error << "FHI-aims has not been installed!!!  ";
          GlobalV::ofs_error << "Suggestion:Install FHI-aims and then re-compile the codes." << std::endl;
          std::exit(EXIT_FAILURE);
        #endif   
        break;
      case 2: //ABACUS
        #ifdef __ABACUS
          // this->char_scf_aims.output_charge_density(file, dens_cmplx);
          GlobalV::ofs_error << "Charge sel-consistent DMFT does not support ABACUS at present!!!  ";
          std::exit(EXIT_FAILURE);
        #else
          GlobalV::ofs_error << "ABACUS has not been installed!!!  ";
          GlobalV::ofs_error << "Suggestion:Install ABACUS and then re-compile the codes." << std::endl;
          std::exit(EXIT_FAILURE);
        #endif
        break;
      default:
        GlobalV::ofs_error << "Not supported DFT_solver" << std::endl;
        std::exit(EXIT_FAILURE);
    }

    ierr = chdir("../");
    if(ierr != 0){
      std::cout << "Process " << mpi_rank() << " fails to return to root directory " << std::endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);  //Blocks until all processes reach here

    timer::get_time(time, seconds, minutes);
    GlobalV::ofs_running << "End non-self-consistent DFT. The time consumption of this DFT step : " 
                         << minutes << "m "
                         << (int)seconds << "s" << std::endl;

    return;
  }

}

