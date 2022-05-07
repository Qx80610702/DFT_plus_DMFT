#pragma once

#include "../para/parameters.h"
#include "projector.h"
#include "chemical_potential.h"
#include "Hilbert_space.h"
#include "Anderson_impurity.h"
#include "coulomb_tensor.h"
#include "spectrum.h"
#include "charge_scf.h"

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

    void working_init();

    void DFT_DMFT_scf();

    void cal_spectrum_func();

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

    void DMFT_charge_updating();
    
    void update_Anderson_impurities();

    void output_to_impurity_solver(
        const int charge_step, 
        const int DMFT_step );

    void reading_inputs();

    void run_nscf_dft(const int dft_solver);

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
    int calculation_type;              //DFT+DMFT scf:0; spectra:1; default:0

  };
}

// typedef DFT_plus_DMFT::solver DFT_DMFT_solver;