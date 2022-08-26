#pragma once

#include "../para/overlap_matrix.h"
#include "../para/KS_eigenvectors.h"
#include "../para/correlated_atoms.h"
#include "../para/KS_bands.h"
#include "Hilbert_space.h"

#include <complex>
#include <vector>


//TEST
#include "../para/KS_bands.h"
namespace DFT_plus_DMFT
{
  class projector
  {
    public:
    projector(){};
    ~projector(){};

    void evaluate_projector(
        const int DFT_solver,
        DFT_output::KS_bands& band,
        DFT_plus_DMFT::Hilbert_space& space,
        DFT_output::atoms_info& atom);

    void evalute_projector_k(
          DFT_output::overlap_matrix& ovlp,
          DFT_output::KS_eigenvectors& wfc,
          DFT_output::atoms_info& atom,
          const int ik,
          const int ik_count);


    //===========================
    //   interface
    //===========================
    std::vector<std::vector<std::vector<std::complex<double>>>>& 
          proj_access(const int ik){return projector_mat[ik];}
    
    std::vector<std::vector<std::vector<std::vector<std::complex<double>>>>>&
          Hloc(){return TB_Hk;}

    const int& nbasis(){return n_basis;} 
    
    private:
    int n_basis;

    //projector_mat[ik][iatom][nspin][nbands*m_tot];
    std::vector<std::vector<std::vector<std::vector<std::complex<double>>>>> projector_mat;
  
    //unitary_trans[ik][iatom][nspin][m_tot*m_tot];
    std::vector<std::vector<std::vector<std::vector<std::complex<double>>>>> unitary_trans;
    //Tight binding Haimiltonian:TB_Hk[ik][iatom][nspin][m_tot];
    std::vector<std::vector<std::vector<std::vector<std::complex<double>>>>> TB_Hk;
    
  };
}
