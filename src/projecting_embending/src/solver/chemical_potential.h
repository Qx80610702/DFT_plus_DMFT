#pragma once

#include "../para/KS_bands.h"
#include "../para/correlated_atoms.h"
#include "self_energy.h"
#include "../para/input.h"
#include "projector.h"
#include "Hilbert_space.h"

#include <complex>
#include <vector>

namespace DFT_plus_DMFT
{
  class chemical_potential
  {
    public:
    chemical_potential(){};
    ~chemical_potential(){};
    
    void update_chemical_potential(
         const int impurity_solver,
         DFT_output::KS_bands& band, 
         DFT_output::atoms_info& atom, 
         DFT_plus_DMFT::projector& proj,
         DMFT::self_energy& sigma,
         DMFT::input_info& in,
         DFT_plus_DMFT::Hilbert_space& space);

    void evaluate_mu_bisection_imag(
         DFT_output::KS_bands& band, 
         DFT_output::atoms_info& atom, 
         DFT_plus_DMFT::projector& proj,
         DMFT::self_energy& sigma,
         DMFT::input_info& in, 
         DFT_plus_DMFT::Hilbert_space& space);
    
    double evaluate_electrons_number_imag(
         DFT_output::KS_bands& band,
         DFT_plus_DMFT::Hilbert_space& space,
         DMFT::self_energy& sigma,
         DFT_output::atoms_info& atom,
         DFT_plus_DMFT::projector& proj,
         const double beta, const int mag, const double mu);

    void evaluate_mu_bisection_imag_DFT(
         DFT_output::KS_bands& band, 
         DFT_output::atoms_info& atom, 
         DMFT::input_info& in, 
         DFT_plus_DMFT::Hilbert_space& space);
    
    double evaluate_electrons_number_imag_DFT(
         DFT_output::KS_bands& band,
         DFT_plus_DMFT::Hilbert_space& space,
         const double beta, const double mu);
    
    //interfaces
    double mu_corrected(){return sigma_corrected_mu;}
    double mu_DFT(){return DFT_mu;}

    private:

    double DFT_mu;
    double sigma_corrected_mu;
    double local_ele_num;

  };
}
