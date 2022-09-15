#include "solver.h"
#include "../mpi_environment.h"
#include "../debug.h"
#include "../timer.h"
#include "projector.h"
#include "../constants.h"
#include "../global_variables.h"
#include "../utilities.h"

#include <mpi.h>
#include <iostream>
#include <complex>
#include <vector>
#include <iomanip>
#include <unistd.h>
#include <stdio.h>

extern "C" void aims_(int* mpi_comm, int* unit, bool* mpi_switch);

namespace DFT_plus_DMFT
{
  void solver::solve()
  {
    debug::codestamp("solver::sovle");
    
    this->set_ios(GLV::ofs_running);

    GLV::ofs_running << "Welcome to DFT+DMFT calculation" << std::endl;
    GLV::ofs_running << "Number of processes of the job: " << mpi_ntasks() << "\n" << std::endl;

    this->pars.in.read();

    this->pars.atom.read(this->pars.in);

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
      std::cerr << "Unknown calculation type" << std::endl;
      std::exit(EXIT_FAILURE);
    }

    this->unset_ios(GLV::ofs_running);

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

    int loop_count = 1, dmft_loop_count = 1, charge_loop_count = 1;
    int mix_step = 1;
    bool density_convergency = false;
    bool sigma_convergency = false;
    int last_char_step = *(int*)this->pars.in.parameter("last_charge_step");
    int last_DMFT_step = *(int*)this->pars.in.parameter("last_dmft_step");

    GLV::ofs_running << "================ Begin DFT+DMFT scf loop ================" << std::endl;
    //===============================================
    //           Charge loop
    //===============================================
    for(int char_step = *(int*)this->pars.in.parameter("start_charge_step"); 
        char_step <= *(int*)this->pars.in.parameter("max_charge_step"); char_step++)
    {
      GLV::ofs_running << "<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>\n";
      GLV::ofs_running << "<><><><><><><><><> Charge step " << char_step << " <><><><><><><><><><><><><>\n"; 
      GLV::ofs_running << "<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>" << std::endl;

      this->pars.bands.read();

      if(loop_count==1){
        if(char_step==1)
          this->space.KS_bands_window(
              0, this->pars.bands, 
              this->pars.atom,
              this->pars.in );
        else
          this->space.KS_bands_window(
              1, this->pars.bands, 
              this->pars.atom,
              this->pars.in );
      }

      // this->Mu.evaluate_mu_bisection_imag_DFT(
      //       this->pars.bands, this->pars.atom, 
      //       this->pars.in, this->space);

      this->proj.evaluate_projector(
            *(int*)pars.in.parameter("DFT_solver"),
            pars.bands, this->space, this->pars.atom );

      this->imp.evaluate_local_occupation( 
            this->pars.in, this->pars.bands, 
            this->proj, this->pars.atom, this->space,
            *(int*)this->pars.in.parameter("n_omega"),
            *(int*)this->pars.in.parameter("magnetism") );

      this->imp.sigma.dc.cal_double_counting( 
            *(int*)this->pars.in.parameter("double_counting"), 
            this->pars.bands.soc(), this->pars.bands.nspins(), 
            this->pars.atom, this->pars.in,
            *(bool*)this->pars.in.parameter("hyb_func"),
            *(double*)this->pars.in.parameter("hyf_xc_alpha") );

      GLV::ofs_running << "\n<><><><><><><><><> DMFT loop <><><><><><><><><>" << std::endl;
      //=============================================================
      //            DMFT loop
      //=============================================================
      for(int dmft_step = (charge_loop_count>1? 1 : *(int*)this->pars.in.parameter("start_dmft_step"));
          dmft_step <= *(int*)this->pars.in.parameter("max_dmft_step"); dmft_step++)
      {
        if(dmft_step>1 && sigma_convergency) break;   //Run at least one dmft step under each charge step

        GLV::ofs_running << "================  DMFT step "  << dmft_step << " ================" << std::endl;

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

        sigma_convergency = this->sigma_scf_update(char_step, dmft_step);

        last_DMFT_step = dmft_step;

        loop_count++;
        dmft_loop_count++;
      }//dmft loop
      GLV::ofs_running << "================ End of DMFT loop ================\n" << std::endl;

      if( *(int*)this->pars.in.parameter("max_charge_step") > 1 && 
          char_step < *(int*)this->pars.in.parameter("max_charge_step") )
      {
        //==============================================================
        //            DMFT charge density updating
        //==============================================================
        GLV::ofs_running << "================ Starting charge updating... ================" << std::endl;
        if(mix_step==1){
          this->Char_scf.init(*(int*)this->pars.in.parameter("dft_solver"),
                *(double*)this->pars.in.parameter("charge_mix_param"),
                *(int*)this->pars.in.parameter("mixing_step"),
                *(double*)this->pars.in.parameter("delta_rho"),
                this->pars.bands.nk(), this->pars.bands.nspins(), proj.nbasis() );
          
          this->Char_scf.prepare_nscf_dft();

          //Reading initial input charge density 
          this->Char_scf.read_charge_density(true, false);

          //Calculate and store initial input DFT charge density matrix
          this->Char_scf.evaluate_DFT_charge_density_matrix(this->pars.bands);
          this->Char_scf.store_initial_DFT_charge_density_matrix();
        }

        this->DMFT_charge_updating();

        density_convergency = this->Char_scf.charge_mixing(mix_step);
        mix_step++;

        if(density_convergency && sigma_convergency){
          GLV::ofs_running << "\nThe DFT+DMFT calculation has converged!!!" << std::endl;
          break;
        }

        GLV::ofs_running << "================ End of charge updating ================\n" << std::endl;
        GLV::ofs_running << "<><><><><><><><><> DFT loop <><><><><><><><><>" << std::endl;
        //=========================================
        //            DFT loop
        //=========================================
        for(int dft_step=1; dft_step <= *(int*)this->pars.in.parameter("max_dft_step"); dft_step++)
        {
          GLV::ofs_running << "================  DFT step "  << dft_step << " ================" << std::endl;
          this->run_nscf_dft(*(int*)pars.in.parameter("dft_solver")); //Run at least one dft step under each charge step
          loop_count++;

          if(dft_step>1 && density_convergency) break;

          if(dft_step < *(int*)this->pars.in.parameter("max_dft_step")){

            density_convergency = this->DFT_loop_charge_mixing(mix_step);
            mix_step++;
              
            GLV::ofs_running << std::endl;
          }
        }//DFT loop
        GLV::ofs_running << "<><><><><><><><><> End of DFT loop <><><><><><><><><>" << std::endl;
      }

      last_char_step = char_step;
      GLV::ofs_running << "<><><><><><><><><> End charge step " << char_step << " <><><><><><><><><><><><><>\n" << std::endl;
      charge_loop_count++;
    }//charge loop

    timer::get_time(time, seconds, minutes, hours, days);
    timer::get_date_time(end_date);

    GLV::ofs_running << "Congratulations!!! The DFT+DMFT calculation finished" << std::endl;
    GLV::ofs_running << "========================== The full time records ==========================" << std::endl;
    GLV::ofs_running << "       starting time              ending time              time-consumption" << std::endl;
    GLV::ofs_running << "      " << start_date
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
    
    this->proj.evaluate_projector(
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

    GLV::ofs_running << "\nImpurity solver starts working......" << std::endl;
    GLV::ofs_running << "  Impurities       starting time              ending time              time-consumption" << std::endl;

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
        std::cerr << "Not supported impurity solver" << std::endl;
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
            freq[iomega] = (2*iomega+1)*GLC::PI/beta;
      }
    }

    return;
  }

  void solver::cal_spectrum_func()
  {
    debug::codestamp("solver::cal_spectrum_func");

    GLV::ofs_running << "\n================ Begin calculating DFT+DMFT spectra ================" << std::endl;

    this->pars.bands.read();

    this->space.KS_bands_window(
          1, this->pars.bands,
          this->pars.atom,
          this->pars.in );

    // this->Mu.evaluate_mu_bisection_imag_DFT(
    //       this->pars.bands, this->pars.atom, 
    //       this->pars.in, this->space);

    this->proj.evaluate_projector(
          *(int*)pars.in.parameter("DFT_solver"),
          pars.bands, this->space, this->pars.atom );

    this->imp.evaluate_local_occupation( 
          this->pars.in, this->pars.bands, 
          this->proj, this->pars.atom, this->space,
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

    GLV::ofs_running << "DFT+DMFT spectra calculation stop successfuly!!!" << std::endl;

    return;
  }
  
  void solver::DMFT_charge_updating()
  {
    debug::codestamp("solver::DMFT_charge_updating");

    GLV::ofs_running << "Start updating DMFT corrected charge density..." << std::endl;

    MPI_Barrier(MPI_COMM_WORLD);  //Blocks until all processes reach here
    int ierr = chdir("./dft");
    if(ierr != 0){
      std::cerr << "Process " << mpi_rank() << " fails to enter to the directory dft" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    if(mpi_rank()==0){
      if(access("outputs_to_DMFT",0)==0){
        // mv_dir("outputs_to_DMFT", ".outputs_to_DMFT");
        system("cp -a outputs_to_DMFT .outputs_to_DMFT");
      }
    }
    ierr = chdir("../");
    if(ierr != 0){
      std::cerr << "Process " << mpi_rank() << " fails to enter to the root directory of the DFT+DMFT calculation" << std::endl;
      std::exit(EXIT_FAILURE);
    }

    MPI_Barrier(MPI_COMM_WORLD);   //Blocks until all processes reach here

    this->imp.sigma.subtract_double_counting(this->flag_axis);

    this->Mu.update_chemical_potential(
            this->flag_axis, this->pars.bands, 
            this->pars.atom, this->proj, 
            this->imp.sigma, this->pars.in, this->space );

    this->Char_scf.update_charge_density_matrix(
            this->flag_axis, this->Mu,
            this->pars.bands, this->pars.atom, 
            this->proj, this->imp.sigma, this->space );

    this->Char_scf.output_DMFT_charge_density_matrix();

    //Call DFT solver to compute the charge density corresponding to 
    //the DMFT corrected density matrix
    // GLV::ofs_running << "Calculating DMFT corrected charge density..." << std::endl;
    this->run_nscf_dft(*(int*)pars.in.parameter("dft_solver"));
  
    //Read the charge density corresponding to 
    //the DMFT corrected density matrix
    this->Char_scf.read_charge_density(false, true);

    MPI_Barrier(MPI_COMM_WORLD);  //Blocks until all processes reach here
    ierr = chdir("./dft");
    if(ierr != 0){
      std::cerr << "Process " << mpi_rank() << " fails to enter to the directory dft" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    if(mpi_rank()==0){
      rm_dir("outputs_to_DMFT");
      mv_dir(".outputs_to_DMFT", "outputs_to_DMFT");
    }
    ierr = chdir("../");
    if(ierr != 0){
      std::cerr << "Process " << mpi_rank() << " fails to enter to the root directory of the DFT+DMFT calculation" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    MPI_Barrier(MPI_COMM_WORLD);  //Blocks until all processes reach here

    GLV::ofs_running << "End updating DMFT corrected charge density\n" << std::endl;
   
    return;
  }

  bool solver::DFT_loop_charge_mixing(
        const int mix_step )
  {
    debug::codestamp("solver::DFT_loop_charge_mixing");
    bool density_convergency = false;

    GLV::ofs_running << "\nStart updating charge density..." << std::endl;

    this->pars.bands.read();

    // this->Mu.evaluate_mu_bisection_imag_DFT(
    //         this->pars.bands, this->pars.atom, 
    //         this->pars.in, this->space);

    this->proj.evaluate_projector(
          *(int*)pars.in.parameter("DFT_solver"),
          pars.bands, this->space, this->pars.atom );

    this->imp.evaluate_local_occupation( 
          this->pars.in, this->pars.bands, 
          this->proj, this->pars.atom, this->space,
          *(int*)this->pars.in.parameter("n_omega"),
          *(int*)this->pars.in.parameter("magnetism") );

    this->imp.sigma.dc.cal_double_counting( 
          *(int*)this->pars.in.parameter("double_counting"), 
          this->pars.bands.soc(), this->pars.bands.nspins(), 
          this->pars.atom, this->pars.in,
          *(bool*)this->pars.in.parameter("hyb_func"),
          *(double*)this->pars.in.parameter("hyf_xc_alpha") );

    this->DMFT_charge_updating();

    GLV::ofs_running << "End updating charge density" << std::endl;

    density_convergency = this->Char_scf.charge_mixing(mix_step);

    return density_convergency;
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

    GLV::ofs_running << "Self-consistency of self-energy in charge step " 
              << charge_step 
              << " DMFT step " << DMFT_step;
    if(convergency)
      GLV::ofs_running << ": true" << std::endl;
    else
      GLV::ofs_running << ": false" << std::endl;

    GLV::ofs_running.setf(std::ios::scientific, std::ios_base::floatfield);
    GLV::ofs_running.setf(std::ios::dec, std::ios_base::basefield);  
    GLV::ofs_running << "    impuritys            Delta_Sigma (eV)" << std::endl;
    for(int ineq=0; ineq<this->pars.atom.inequ_atoms(); ineq++){
      GLV::ofs_running << "    impurity" << ineq << "              "
                // << std::setprecision(3) << std::setiosflags(std::ios::scientific)
                << std::setprecision(3) << this->imp.delta_scf()[ineq] 
                << std::endl;
    }
    GLV::ofs_running << std::endl;
    GLV::ofs_running.unsetf(std::ios::scientific);

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
// GLV::ofs_running << "evaluate_hybridization_function_imag time consuming " << seconds << '\n';

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
    std::ofstream& ofs_running )
  {
    if(mpi_rank()==0){
      ofs_running.open("DMFT_running.log", std::ios::app);
      // ofs_error.open("DMFT_running.error", std::ios::out);
    }

    //debug
    bool debug = false;
    if(debug){
      std::stringstream ss1;
      ss1 << "debug_proc" << mpi_rank() << ".log";
      GLV::ofs_debug.open(ss1.str(), std::ios_base::app);
    }
    
    // std::stringstream ss2;
    // ss2 << "DMFT_running_proc" << mpi_rank() << ".error";
    // ofs_error.open(ss2.str(), std::ios_base::out);

  }

  void solver::unset_ios(
    std::ofstream& ofs_running)
  {
    if(mpi_rank()==0){
      ofs_running.close();
      // ofs_error.close();
    }

    // ofs_running.close();
    // ofs_error.close();
  }

  void solver::run_nscf_dft(
      const int dft_solver )
  {
    debug::codestamp("solver::run_nscf_dft");

    double time, seconds;
    int minutes;
    timer::timestamp(time);

    GLV::ofs_running << "Start non-self-consistent DFT..." << std::endl;

    int ierr = chdir("./dft");
    if(ierr != 0){
      std::cerr << "Process " << mpi_rank() << " fails to enter to the directory dft" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    if(mpi_rank()==0){
      if(access("outputs_to_DMFT",0)==0){
        int info = rm_dir("outputs_to_DMFT");
        if(info==-1){
          std::cerr << "Fail to remove the directory outputs_to_DMFT!!!" << std::endl;
          std::exit(EXIT_FAILURE);
        }
      }
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
          
          /*
          if(mpi_rank()==0){
            std::string dft_exe = "mpirun " + *(std::string*)this->pars.in.parameter("dft_solver_exe") + " 1>job.log 2>job.error";
            // std::string dft_exe = *(std::string*)this->pars.in.parameter("dft_solver_exe") + " 1>job.log 2>job.error";
            std::ofstream ofs("task.sh", std::ios::out);
            ofs << "#!/bin/sh\n" << std:: endl;
            ofs << dft_exe << std::endl;
            ofs.close();

            system("bash task.sh");
          }
          */
         
        #else
          std::cerr << "FHI-aims has not been installed!!!  ";
          std::cerr << "Suggestion:Install FHI-aims and then re-compile the codes." << std::endl;
          std::exit(EXIT_FAILURE);
        #endif   
        break;
      case 2: //ABACUS
        #ifdef __ABACUS
          // this->char_scf_aims.output_charge_density(file, dens_cmplx);
          std::cerr << "Charge sel-consistent DMFT does not support ABACUS at present!!!  ";
          std::exit(EXIT_FAILURE);
        #else
          std::cerr << "ABACUS has not been installed!!!  ";
          std::cerr << "Suggestion:Install ABACUS and then re-compile the codes." << std::endl;
          std::exit(EXIT_FAILURE);
        #endif
        break;
      default:
        std::cerr << "Not supported DFT_solver" << std::endl;
        std::exit(EXIT_FAILURE);
    }
    
    MPI_Barrier(MPI_COMM_WORLD);  //Blocks until all processes reach here

    ierr = chdir("../");
    if(ierr != 0){
      std::cerr << "Process " << mpi_rank() << " fails to return to root directory " << std::endl;
      std::exit(EXIT_FAILURE);
    }

    timer::get_time(time, seconds, minutes);
    GLV::ofs_running << "End non-self-consistent DFT.\nThe time consumption of this DFT step: " 
                         << minutes << "m "
                         << (int)seconds << "s" << std::endl;

    return;
  }

}

