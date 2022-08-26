#pragma once

#include <fstream>
#include <iostream>

namespace DFT_plus_DMFT
{
  namespace GlobalV
  {
    extern std::ofstream ofs_running;
    // extern std::ofstream ofs_error;
    extern std::ofstream ofs_debug;
  }
}

namespace GLV = DFT_plus_DMFT::GlobalV;