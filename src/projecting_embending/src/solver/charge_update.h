
#ifdef __FHIaims
#include "charge_update_aims.h"
#endif

namespace DFT_plus_DMFT
{
  class Charge_update
  {
    public:
    Charge_update(){;}
    ~Charge_update(){;}

    public:
    #ifdef __FHIaims
    Charge_update_aims charge_aims;
    #endif

    public:
    void update_char_dens();
    void eva_char_dens();

  };
}