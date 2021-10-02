#include "coulomb_tensor.h"

#include "../debug/debug.h"

namespace DMFT
{
  void coulomb_tensor::update_coulomb_tensor(
    const int type, DFT_output::atoms_info& atom)
  {
    switch(type)
    {
      case 1:
        this->Kana_coulomb.evaluate_coulomb_tensor(atom, this->U_matrix);
        break;
      default:
        break;
    }
  }

  void coulomb_tensor::out_coulomb_tensor(
            const int type, const int istep, 
            const int impurity_solver, 
            DFT_output::atoms_info& atom,
            DFT_output::KS_bands& band)
  {
    switch(type)
    {
      case 1:
        if(impurity_solver==1)
          this->Kana_coulomb.out_ALPS_CTHYB(istep, atom, band);
        else if(impurity_solver==2)
          this->Kana_coulomb.out_ALPS_CTHYB_SEGMENT(istep, atom, band);
        break;
      default:
        break;
    }
  }
  
}
