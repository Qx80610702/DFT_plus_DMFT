#include "mpi_environment.h"
#include <mpi.h>

int mpi_ntasks()
{
  int nprocs;
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  return nprocs;
}

int mpi_rank()
{
  int myid;
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  return myid;
}
