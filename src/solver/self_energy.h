#pragma once

#include "../para/correlated_atoms.h"
#include "self_energy_real_aixs.h"
#include "self_energy_imaginary_aixs.h"
#include "double_counting.h"
#include "projector.h"

#include <vector>

#include <complex>

namespace DMFT
{
  class self_energy
  {
    public:
    self_energy();
    ~self_energy();

    self_energy_real_aixs sigma_real;
    self_energy_imaginary_aixs sigma_imag;
    double_counting dc;

    bool read_sigma_save(const int impurity_solver, const bool SOC,
                         const int nspin, DFT_output::atoms_info& atom);
    
    void initial_guess(const int axis_flag, const bool SOC, 
                       const int nspin, DFT_output::atoms_info& atom);

    void subtract_double_counting(const int axis_flag);  //0: imaginary axis, 1:real axis

    void evalute_lattice_sigma(
        const int axis_flag, const int mag, 
        const int ispin, const std::vector<int>& wbands, 
        DFT_output::atoms_info& atom,
        const std::vector<std::vector<std::vector<
        std::complex<double>>>>&  projector,
        std::vector<std::vector<std::complex<double>>>& Simga);

    void evalute_lattice_sigma_infty(
        const int axis_flag, const int mag, 
        const int nspin, const std::vector<int>& wbands, 
        DFT_output::atoms_info& atom,
        const std::vector<std::vector<std::vector<
        std::complex<double>>>>&  projector,
        std::vector<std::complex<double>>& Simga_infty);

    // void embeding_self_energy(const int impurity_solver, )

    void out();    //test whether reading worked corrected

    // =================================================
    //             interface
    // =================================================
    int nomega(const int axis_flag); //0: imaginary axis, 1:real axis

    std::vector<std::vector<std::vector<std::vector<std::complex<double>>>>>&
              sigma_new(const int axis_flag); //0: imaginary axis, 1:real axis

    std::vector<std::vector<std::vector<std::vector<std::complex<double>>>>>&
              sigma_save(const int axis_flag); //0: imaginary axis, 1:real axis

    std::vector<std::vector<std::vector<std::vector<std::complex<double>>>>>&
              correlated_sigma(const int axis_flag); //0: imaginary axis, 1:real axis
  
  };
}
