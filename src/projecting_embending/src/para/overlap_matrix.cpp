#include "overlap_matrix.h"
#include "../debug.h"

#include <iostream>
#include <cstdlib>

namespace DFT_output
{
  overlap_matrix::overlap_matrix(const int DFT_solver, const std::string dir)
  : flag_DFT_solver(DFT_solver)
  {
    switch(DFT_solver)
    {
      case 1: //aims
        this->ovlp_aims.read(dir);
        break;
      case 2: //ABACUS
        this->ovlp_abacus.read(dir);
        break;
      default:
        std::cout << "Not supported DFT_solver" << std::endl;
        std::exit(EXIT_FAILURE);
    }
  }


  void overlap_matrix::evaluate_ovlp_k(const int ik, atoms_info& at_info)
  {
    switch(this->flag_DFT_solver)
    {
      case 1: //aims
        this->ovlp_aims.evaluate_ovlp_k(ik, at_info);
        break;
      case 2: //ABACUS
        this->ovlp_abacus.evaluate_ovlp_k(ik, at_info);
        break;
      default:
        std::cout << "Not supported DFT_solver" << std::endl;
        std::exit(EXIT_FAILURE);
    }  

  }

  std::vector<std::complex<double>>& overlap_matrix::overlap()
  {
    switch(this->flag_DFT_solver)
    {
      case 1: //aims
        return this->ovlp_aims.ovlp_mat();
        break;
      case 2: //ABACUS
        return this->ovlp_abacus.ovlp_mat();
        break;
      default:
        std::cout << "Not supported DFT_solver" << std::endl;
        std::exit(EXIT_FAILURE);
    } 
  }

  std::vector<std::complex<double>>& overlap_matrix::overlap_localorb()
  {
    switch(this->flag_DFT_solver)
    {
      case 1: //aims
        return this->ovlp_aims.local_ovlp();
        break;
      case 2: //ABACUS
        return this->ovlp_abacus.local_ovlp();
        break;
      default:
        std::cout << "Not supported DFT_solver" << std::endl;
        std::exit(EXIT_FAILURE);
    } 
  }

  void overlap_matrix::out()
  {
    switch(this->flag_DFT_solver)
    {
      case 1: //aims
        this->ovlp_aims.out();
        break;
      case 2: //ABACUS

        break;
      default:
        std::cout << "Not supported DFT_solver" << std::endl;
        std::exit(EXIT_FAILURE);
    }
  }


}
