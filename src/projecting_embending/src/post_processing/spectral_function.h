#pragma once

#include "../para/input.h"
#include "../para/correlated_atoms.h"
#include "../para/KS_bands.h"

#include <vector>
#include <complex>

class spectral_function
{
  public:
  spectral_function(){;}
  ~spectral_function(){;}

  void evaluate_local_spectrum(
        DMFT::input_info& in,
        DFT_output::atoms_info& atom,
        DFT_output::KS_bands& band,
        std::vector<std::vector<std::vector<
        std::vector<std::complex<double>>>>>& Gf,
        std::vector<std::vector<std::vector<
        std::vector<std::complex<double>>>>>& Gf_save,
        std::vector<double>& freq     
        );
  
};

