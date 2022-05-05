#pragma once

#include "../para/parameters.h"
#include "projector.h"
#include "chemical_potential.h"
#include "Hilbert_space.h"
#include "Anderson_impurity.h"
#include "coulomb_tensor.h"
#include "spectrum.h"
#include "charge_scf.h"

//================Arguments lists(command line)==================
// -dmft.step : current dmft step, default value : 1;
// -eva.sigma_only : whether calculate self-energy only 
//  (only for impurity solver ALPS-CTHYB), default value : false
// ==============================================================
struct argument_lists
{
  int current_charge_step;
  int current_DMFT_step;
  int last_charge_step;
  int last_DMFT_step;
  bool sigma_only;
  bool cal_spectrum;
  bool update_density;
};

namespace DFT_plus_DMFT
{
  class solver
  {
    public:
    solver(){};
    ~solver(){};

    param pars;
    projector proj;
    chemical_potential Mu;
    Hilbert_space space;
    DMFT::impurity imp;
    DMFT::coulomb_tensor Umat;
    spectrum Aw;
    Charge_SCF Char_scf;
    
    void solve();

    void projecting_embeding(
        const int charge_step, 
        const int DMFT_step);

    void impurity_solving(
        const int impurity_solver,
        const int char_step,
        const int DMFT_step,
        DFT_output::atoms_info& atom );
    
    bool sigma_scf_update(
        const int charge_step, 
        const int DMFT_step);

    void init_working_axis();

    void cal_spectrum_func();

    void charge_solve();
    
    void update_Anderson_impurities();

    void output_to_impurity_solver(
        const int charge_step, 
        const int DMFT_step );

    void reading_inputs();

    public:  //static member function
    static void set_ios(std::ofstream& ofs_running, std::ofstream& ofs_error);

    //==================================
    //     interfaces
    //==================================
    // int& curr_dmft_step(){return current_DMFT_step;}
    // int& curr_char_step(){return current_charge_step;}
    // bool convergency(){return this->flag_convergency;}

    private:
    int flag_axis;                     //Do calculations on imaginary (0) or real (1) axis
    // int current_DMFT_step;             //The step number of current DMFT iteration
    // int current_charge_step;           //The step number of current charge scf iteration
    // int last_DMFT_step;                //The step number of current DMFT iteration
    // int last_charge_step;              //The step number of current charge scf iteration
    // bool flag_convergency;             //Whether DMFT iteration convergrncy is achieved
    // bool flag_eva_sigma_only;          //
    bool flag_eva_spectrum;            //calculate spectrum
    // bool flag_update_density;          //update charge density

  };
}

// typedef DFT_plus_DMFT::solver DFT_DMFT_solver;