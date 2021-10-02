#pragma once

#include "../para/parameters.h"
#include "projector.h"
#include "chemical_potential.h"
#include "Hilbert_space.h"
#include "Anderson_impurity.h"
#include "coulomb_tensor.h"
#include "../post_processing/post_processing.h"

//================Arguments lists(command line)==================
// -dmft.step : current dmft step, default value : 1;
// -eva.sigma_only : whether calculate self-energy only 
//  (only for impurity solver ALPS-CTHYB), default value : false
// ==============================================================
struct argument_lists
{
  int global_step;
  bool sigma_only;
};

namespace DFT_plus_DMFT
{
  class solver
  {
    public:
    solver(argument_lists& args);
    ~solver(){};

    param pars;
    projector proj;
    chemical_potential Mu;
    Hilbert_space space;
    DMFT::impurity imp;
    DMFT::coulomb_tensor Umat;
    post_processing post;
    
    bool solve();

    bool scf_update();
    
    void update_Anderson_impurities();

    void output_to_impurity_solver();

    void post_processing();

    void reading_inputs();

    //==================================
    //     interfaces
    //==================================
    int& istep(){return DMFT_iteration_step;}
    bool convergency(){return this->flag_convergency;}

    private:
    int DMFT_iteration_step;           //The step number of current DMFT iteration
    bool flag_convergency;             //Whether DMFT iteration convergrncy is achieved
    bool flag_eva_sigma_only;         //

  };
}
typedef DFT_plus_DMFT::solver DFT_DMFT_solver;