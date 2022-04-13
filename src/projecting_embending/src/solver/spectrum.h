#pragma once

#include "../para/KS_bands.h"
#include "../para/correlated_atoms.h"
#include "self_energy.h"
#include "../para/input.h"
#include "projector.h"
#include "Hilbert_space.h"

#include <vector>

class spectrum
{
  public:
  spectrum(){;}
  ~spectrum(){;}

  public:
  void eva_spectrum(
        const double mu,
        DFT_output::KS_bands& band, 
        DFT_output::atoms_info& atom, 
        DFT_plus_DMFT::projector& proj,
        DMFT::self_energy& sigma,
        DMFT::input_info& in,
        DFT_plus_DMFT::Hilbert_space& space );

  void eva_spectrum_normalization(
        const double mu,
        DFT_output::KS_bands& band, 
        DFT_output::atoms_info& atom, 
        DFT_plus_DMFT::projector& proj,
        DMFT::self_energy& sigma,
        DMFT::input_info& in,
        DFT_plus_DMFT::Hilbert_space& space );

  void out_spectrum();

  private:
  //Spectrum function; Awk[is][iomega][ik]
  std::vector<std::vector<std::vector<double>>> Awk;

  //Spectrum function; DOS[is][iomega]
  std::vector<std::vector<double>> DOS;

  //Spectrum function; Aw_loc[ineq][is][iomega][m] 
  std::vector<std::vector<std::vector<std::vector<double>>>> Aw_loc;

  //Real frequency  
  std::vector<double> freq;

  //Spectrum function; Aw_iband_ik[is][ik][iband][iomega]
  std::vector<std::vector<std::vector<std::vector<double>>>> Aw_ik_iband;

};