#include "self_energy.h"

#include "../debug.h"
#include "../timer.h"
#include "../constants.h"
#include "../global_variables.h"

#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <cstdlib>

#define MKL_Complex16 std::complex<double>
#include <mkl.h>

namespace DMFT
{
  self_energy::self_energy(){}

  self_energy::~self_energy(){}

  bool self_energy::read_sigma_save(const int impurity_solver, const bool SOC,
                    const int nspin, DFT_output::atoms_info& atom)
  {
    if(impurity_solver==1)
    {
      this->sigma_imag.read_sigma_save(SOC, nspin, atom);
    }
    return true;
  }

  void self_energy::initial_guess(const int axis_flag, const bool SOC,
            const int nspin, DFT_output::atoms_info& atom)
  {
    debug::codestamp("self_energy::initial_guess");

    const std::vector<int>& norb_sub = atom.iatom_norb();
    const std::complex<double> zero(0.0,0.0);
    
    std::vector<std::vector<std::vector<std::vector<std::complex<double>>>>>&
            Sigma_new = this->sigma_new(axis_flag);
    
    std::vector<std::vector<std::vector<std::vector<std::complex<double>>>>>&
            Sigma_save = this->sigma_save(axis_flag);

    Sigma_new.resize(atom.inequ_atoms());
    for(int ineq=0; ineq<atom.inequ_atoms(); ineq++)
    {
      const int iatom = atom.ineq_iatom(ineq);
      const int m_tot=norb_sub[iatom];
      int nspin_tmp=nspin;
      // if(SOC) nspin_tmp=4;

      std::vector<std::vector<std::complex<double>>>& 
          Vdca = this->dc.Vdc()[ineq];

      const int n_omega=this->nomega(axis_flag);
        
      Sigma_new[ineq].resize(nspin_tmp);
 
      for(int is=0; is<nspin_tmp; is++)
      {
        Sigma_new[ineq][is].resize(n_omega);

        for(int i_omega=0; i_omega<n_omega; i_omega++)
        {
          Sigma_new[ineq][is][i_omega].resize(m_tot*m_tot, zero);

          for(int m_index=0; m_index<m_tot*m_tot; m_index++)
              Sigma_new[ineq][is][i_omega][m_index] = Vdca[is][m_index];
        }//is
      }//i_omega
    }//iatom

    Sigma_save = Sigma_new;

    return;
  }

  void self_energy::evalute_lattice_sigma(
        const int axis_flag, const int mag, 
        const int ispin, const std::vector<int>& wbands,
        DFT_output::atoms_info& atom, 
        const std::vector<std::vector<std::vector<
        std::complex<double>>>>&  projector,
        std::vector<std::vector<std::complex<double>>>& Simga )
  {
    debug::codestamp("self_energy::evalute_lattice_sigma");

    const int natom = atom.total_atoms();
    const std::vector<int>& norb_sub = atom.iatom_norb();
    const int norb = atom.norb();
    const int omega_num = this->nomega(axis_flag);
    const std::complex<double> one(1.0,0.0), zero(0.0,0.0);

    if(omega_num != Simga.size()){
      std::cerr << "Fatal error in calculating lattice self-energy!!!" << std::endl;
      std::exit(EXIT_FAILURE);
    }

    // int nbands = wbands[ispin];
    // if(wbands.size()==2) nbands = wbands[0] > wbands[1] ? wbands[0] : wbands[1];

    const std::vector<std::vector<std::vector<std::vector<std::complex<double>>>>>&
          dleta_sigma = this->correlated_sigma(axis_flag);

    for(std::vector<std::complex<double>>& iter1 : Simga)
      for(std::complex<double>& iter2 : iter1)
        iter2 = zero;

    for(int ineq=0; ineq<atom.inequ_atoms(); ineq++)
    {
      const int iatom = atom.ineq_iatom(ineq);
      const int m_tot = norb_sub[iatom];
      
      std::vector<std::complex<double>> mat_tmp(wbands[ispin]*m_tot);
      for(int iomega=0; iomega<omega_num; iomega++)
      {
        for(int iband=0; iband<wbands[ispin]; iband++)
          for(int m=0; m<m_tot; m++)
            mat_tmp[iband*m_tot+m] = projector[iatom][ispin][iband*m_tot+m]*
                                    dleta_sigma[ineq][ispin][iomega][m*m_tot+m];

        // void cblas_zgemm (const CBLAS_LAYOUT Layout, const CBLAS_TRANSPOSE transa, const
        //      CBLAS_TRANSPOSE transb, const MKL_INT m, const MKL_INT n, const MKL_INT k, const void
        //      *alpha, const void *a, const MKL_INT lda, const void *b, const MKL_INT ldb, const void
        //      *beta, void *c, const MKL_INT ldc);

        cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasConjTrans,
                    wbands[ispin], wbands[ispin], m_tot,
                    &one,
                    &mat_tmp[0], m_tot,
                    &projector[iatom][ispin][0], m_tot,
                    &one,
                    &Simga[iomega][0], wbands[ispin]);
     
        for(int iatom1=0; iatom1<natom; iatom1++)
        {
          if(iatom1==iatom) continue;
          else if(iatom1!=iatom && atom.equ_atom(iatom1)==ineq)
          {
            int symm=0;            //symmetry
            if(mag==1)             //AFM
            {
              if(atom.magnetic(iatom)*atom.magnetic(iatom1)==-1) symm=1;  //anti_symmetry
            }

            for(int iband=0; iband<wbands[ispin]; iband++)
              for(int m=0; m<m_tot; m++)
                mat_tmp[iband*m_tot+m] = projector[iatom1][ispin][iband*m_tot+m]*
                            dleta_sigma[ineq][(ispin+symm)%2][iomega][m*m_tot+m];

            cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasConjTrans,
                        wbands[ispin], wbands[ispin], m_tot,
                        &one,
                        &mat_tmp[0], m_tot,
                        &projector[iatom1][ispin][0], m_tot,
                        &one,
                        &Simga[iomega][0], wbands[ispin]);

          }
        }//iatom1
      }//iomega
    }//ineq

    return;
  }

  void self_energy::evalute_lattice_sigma_infty(
        const int axis_flag, const int mag, 
        const int ispin, const std::vector<int>& wbands, 
        DFT_output::atoms_info& atom,
        const std::vector<std::vector<std::vector<
        std::complex<double>>>>&  projector,
        std::vector<std::complex<double>>& Simga_infty)
  {
    debug::codestamp("self_energy::evalute_lattice_sigma_infty");

    const int natom = atom.total_atoms();
    const std::vector<int>& norb_sub = atom.iatom_norb();
    const int norb = atom.norb();
    const int omega_num = this->nomega(axis_flag);
    const std::complex<double> one(1.0,0.0), zero(0.0,0.0);

    // int nbands = wbands[ispin];
    // if(wbands.size()==2) nbands = wbands[0] > wbands[1] ? wbands[0] : wbands[1];

    const std::vector<std::vector<std::vector<std::vector<std::complex<double>>>>>&
          dleta_sigma = this->correlated_sigma(axis_flag);

    for(std::complex<double>& iter : Simga_infty)
      iter = zero;

    for(int ineq=0; ineq<atom.inequ_atoms(); ineq++)
    {
      const int iatom = atom.ineq_iatom(ineq);
      const int m_tot = norb_sub[iatom];

      std::vector<std::complex<double>> mat_tmp(wbands[ispin]*m_tot);

      for(int iband=0; iband<wbands[ispin]; iband++)
        for(int m=0; m<m_tot; m++)
          mat_tmp[iband*m_tot+m] = projector[iatom][ispin][iband*m_tot+m]*
                                  dleta_sigma[ineq][ispin][omega_num-1][m*m_tot+m].real();

      // void cblas_zgemm (const CBLAS_LAYOUT Layout, const CBLAS_TRANSPOSE transa, const
      //      CBLAS_TRANSPOSE transb, const MKL_INT m, const MKL_INT n, const MKL_INT k, const void
      //      *alpha, const void *a, const MKL_INT lda, const void *b, const MKL_INT ldb, const void
      //      *beta, void *c, const MKL_INT ldc);

      cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasConjTrans,
                  wbands[ispin], wbands[ispin], m_tot,
                  &one,
                  &mat_tmp[0], m_tot,
                  &projector[iatom][ispin][0], m_tot,
                  &one,
                  &Simga_infty[0], wbands[ispin]);
   
      for(int iatom1=0; iatom1<natom; iatom1++)
      {
        if(iatom1==iatom) continue;
        else if(iatom1!=iatom && atom.equ_atom(iatom1)==ineq)
        {
          int symm=0;            //symmetry
          if(mag==1)             //AFM
          {
            if(atom.magnetic(iatom)*atom.magnetic(iatom1)==-1) symm=1;  //anti_symmetry
          }

          for(int iband=0; iband<wbands[ispin]; iband++)
            for(int m=0; m<m_tot; m++)
              mat_tmp[iband*m_tot+m] = projector[iatom1][ispin][iband*m_tot+m]*
                          dleta_sigma[ineq][(ispin+symm)%2][omega_num-1][m*m_tot+m];
          
          cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasConjTrans,
                      wbands[ispin], wbands[ispin], m_tot,
                      &one,
                      &mat_tmp[0], m_tot,
                      &projector[iatom1][ispin][0], m_tot,
                      &one,
                      &Simga_infty[0], wbands[ispin]);
        }
      }//iatom1
    }//ineq

    return;
  }

  void self_energy::subtract_double_counting(const int axis_flag)
  {
    debug::codestamp("self_energy::subtract_double_counting");

    std::vector<std::vector<std::vector<std::vector<std::complex<double>>>>>&
            sigma_correlatd = this->correlated_sigma(axis_flag);
    
    const std::vector<std::vector<std::vector<std::vector<std::complex<double>>>>>&
            sigma_new = this->sigma_new(axis_flag);

    //Allocation
    if(sigma_correlatd.empty()) sigma_correlatd = sigma_new;
    // {
    //   sigma_correlatd.resize(this->dc.Vdc().size());
    //   for(int ineq=0; ineq<this->dc.Vdc().size(); ineq++)
    //   {
    //     sigma_correlatd[ineq].resize(this->dc.Vdc()[ineq].size());    
    //     for(int is=0; is<this->dc.Vdc().at(0).size(); is++)
    //     {
    //       sigma_correlatd[ineq][is].resize(this->nomega(axis_flag));
    //       for(int iomega=0; iomega<this->nomega(axis_flag); iomega++)
    //         sigma_correlatd[ineq][is][iomega].resize(this->dc.Vdc()[ineq][is].size());
    //     }//i_omega
    //   }//ineq
    // }

    for(int ineq=0; ineq<this->dc.Vdc().size(); ineq++)
    {
      const int nspin = this->dc.Vdc().at(0).size();
      const int omega_number = this->nomega(axis_flag);

      for(int ispin=0; ispin<nspin; ispin++)
      {
        const int m_index_number = this->dc.Vdc().at(ineq).at(ispin).size();
        for(int m_index=0; m_index<m_index_number; m_index++)
        {
          std::complex<double> tmp=this->dc.Vdc().at(ineq).at(ispin).at(m_index);
          for(int iomega=0; iomega<omega_number; iomega++){
            sigma_correlatd[ineq][ispin][iomega][m_index] = sigma_new[ineq][ispin][iomega][m_index] - tmp;
          }
        }//m_index
      }//ispin
    }//ineq

    return;  
  }

  int self_energy::nomega(const int axis_flag)
  {
    switch(axis_flag)
    {
    case 0:
      return this->sigma_imag.nomega();
      break;
    case 1:
      return this->sigma_real.nomega();
      break;
    default:
      std::cerr << "Error parameter of axis_flag" << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }

  std::vector<std::vector<std::vector<std::vector<std::complex<double>>>>>&
              self_energy::sigma_new(const int axis_flag)
  {
    switch(axis_flag)
    {
    case 0:
      return this->sigma_imag.sigma_new_access();
      break;
    case 1:
      return this->sigma_real.sigma_new_access();
      break;
    default:
      std::cerr << "Error parameter of axis_flag" << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }

  std::vector<std::vector<std::vector<std::vector<std::complex<double>>>>>&
            self_energy::sigma_save(const int axis_flag)
  {
    switch(axis_flag)
    {
    case 0:
      return this->sigma_imag.sigma_save_access();
      break;
    case 1:
      return this->sigma_real.sigma_save_access();
      break;
    default:
      std::cerr << "Error parameter of axis_flag" << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }

  std::vector<std::vector<std::vector<std::vector<std::complex<double>>>>>&
              self_energy::correlated_sigma(const int axis_flag)
  {
    switch(axis_flag)
    {
    case 0:
      return this->sigma_imag.correlated_sigma_access();
      break;
    case 1:
      return this->sigma_real.correlated_sigma_access();
      break;
    default:
      std::cerr << "Error parameter of axis_flag" << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }
  
}