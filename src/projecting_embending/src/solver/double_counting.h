#pragma once

#include "../para/correlated_atoms.h"
#include "../para/input.h"

#include <complex>
#include <vector>

namespace DMFT
{
  class double_counting
  {
    public:
    double_counting(){};
    ~double_counting(){};

    void cal_double_counting(const int type, 
            const bool SOC, const int nspin, 
            DFT_output::atoms_info& atom,
            DMFT::input_info& in );

    //============================================
    //        interface
    //============================================
    inline std::vector<std::vector<std::vector<std::complex<double>>>>&
    Vdc(){return V_dc;}

    private:
    std::vector<std::vector<std::vector<std::complex<double>>>> V_dc;   //V_dc[ineq][ispin][m_index]

  };
}
