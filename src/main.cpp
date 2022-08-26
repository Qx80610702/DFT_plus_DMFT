#include "./solver/solver.h"

#include <mpi.h>

int main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);

  DFT_plus_DMFT::solver DFT_DMFT_solver;
  DFT_DMFT_solver.solve();

  MPI_Finalize();

  return 0;
}


