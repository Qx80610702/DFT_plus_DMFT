#pragma once

#include <string>

#include "input.h"
#include "KS_bands.h"
#include "correlated_atoms.h"
#include "tetrahedron.h"

class parameters
{
  public:
  parameters(){;}
  ~parameters(){;}

  DMFT::input_info in;
  DFT_output::atoms_info atom;
  DFT_output::KS_bands bands;
  DFT_output::tretrahedron tetra;

  void out();

};
  
typedef parameters param;
