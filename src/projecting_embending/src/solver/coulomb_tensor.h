#pragma once

#include "Kanamori_parameterization.h"
#include "../para/correlated_atoms.h"
#include "../para/KS_bands.h"

#include <vector>

namespace DMFT
{
  class coulomb_tensor
  {
    public:
    coulomb_tensor(){};
    ~coulomb_tensor(){};

    Kanamori_para Kana_coulomb;

    void update_coulomb_tensor(const int type, DFT_output::atoms_info& atom);

    void out_coulomb_tensor(const int type, const int istep, 
                            const int impurity_solver, 
                            DFT_output::atoms_info& atom,
                            DFT_output::KS_bands& band);

    std::vector<std::vector<std::vector<std::vector<
    std::vector<double>>>>>& Coulomb_matrix(){return U_matrix;}

    private:

    std::vector<std::vector<std::vector<std::vector<
    std::vector<double>>>>> U_matrix; //U_matrix[ineq][m1][m2][m3][m4]

  };
}
