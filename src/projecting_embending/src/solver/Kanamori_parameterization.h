//======================================================================
//   Acta Cryst.(2014). c70, 137-159. DOI:10.1107/s2053229613032312
//======================================================================
#pragma once

#include "../para/correlated_atoms.h"
#include "../para/KS_bands.h"

#include <vector>

namespace DMFT
{
  class Kanamori_parameterization
  {
    public:
    Kanamori_parameterization(){};
    ~Kanamori_parameterization(){};

    void evaluate_coulomb_tensor(DFT_output::atoms_info& atom,
              std::vector<std::vector<std::vector<std::vector<
              std::vector<double>>>>>& U_matrix);

    void out_ALPS_CTHYB(const int istep, DFT_output::atoms_info& atom,
                          DFT_output::KS_bands& band);

    void out_ALPS_CTHYB_SEGMENT(const int istep, DFT_output::atoms_info& atom,
                          DFT_output::KS_bands& band);

    private:
    
  };
}
typedef DMFT::Kanamori_parameterization Kanamori_para;
