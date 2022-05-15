#include "self_energy.h"

#include "../debug.h"
#include "../timer.h"
#include "../constants.h"
#include "../global_variables.h"

#include <omp.h>
#include <memory>
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
        const int nspin, const std::vector<int>& wbands,
        DFT_output::atoms_info& atom, 
        const std::vector<std::vector<std::vector<
        std::complex<double>>>>&  projector )
  {
    debug::codestamp("self_energy::evalute_lattice_sigma");

    const int natom = atom.total_atoms();
    const std::vector<int>& norb_sub = atom.iatom_norb();
    const int norb = atom.norb();
    const int omega_num = this->nomega(axis_flag);
    const std::complex<double> one(1.0,0.0), zero(0.0,0.0);

    int nbands = wbands[0];
    if(nspin==2) nbands = wbands[0] > wbands[1] ? wbands[0] : wbands[1];

    std::vector<std::vector<std::vector<std::complex<double>>>>&
          latt_sigma=this->lattice_sigma(axis_flag);

    const std::vector<std::vector<std::vector<std::vector<std::complex<double>>>>>&
          dleta_sigma = this->correlated_sigma(axis_flag);

    if(latt_sigma.empty())
    {
      latt_sigma.resize(nspin);
      for(int ispin=0; ispin<nspin; ispin++)
      {
        latt_sigma[ispin].resize(omega_num);
        for(int iomega=0; iomega<omega_num; iomega++)
        {
          latt_sigma[ispin][iomega].resize(wbands[ispin]*wbands[ispin], zero);
        }
      }
    }
    else
    {
      #pragma omp parallel for    //expensive
      for(int index=0; index<omega_num*nspin; index++) //omega and spin
      {
        int ispin = index/omega_num;
        int iomega = index%omega_num;

        auto& la=latt_sigma[ispin][iomega];
        for(int band_index=0; band_index<wbands[ispin]*wbands[ispin]; band_index++)
          la[band_index] = zero;
      }
    }

    for(int ineq=0; ineq<atom.inequ_atoms(); ineq++)
    {
      const int iatom = atom.ineq_iatom(ineq);
      const int m_tot = norb_sub[iatom];

      const int mkl_threads = mkl_get_max_threads();
      mkl_set_num_threads(1);      //set the number of threads of MKL library function to 1   
      #pragma omp parallel
      {        
        std::vector<std::complex<double>> mat_tmp(nbands*m_tot);

        #pragma omp for
        for(int index=0; index<omega_num*nspin; index++)
        {
          int ispin = index/omega_num;
          int iomega = index%omega_num;

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
                      &latt_sigma[ispin][iomega][0], wbands[ispin]);
     
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
                          &latt_sigma[ispin][iomega][0], wbands[ispin]);

            }
          }//iatom1
        }//index
      }
      mkl_set_num_threads(mkl_threads);
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
    if(sigma_correlatd.empty())
    {
      sigma_correlatd.resize(this->dc.Vdc().size());
      for(int ineq=0; ineq<this->dc.Vdc().size(); ineq++)
      {
        sigma_correlatd[ineq].resize(this->dc.Vdc()[ineq].size());    
        for(int is=0; is<this->dc.Vdc().at(0).size(); is++)
        {
          sigma_correlatd[ineq][is].resize(this->nomega(axis_flag));
          for(int iomega=0; iomega<this->nomega(axis_flag); iomega++)
            sigma_correlatd[ineq][is][iomega].resize(this->dc.Vdc()[ineq][is].size());
        }//i_omega
      }//ineq
    }

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
          // sigma_correlatd[ineq][ispin][iomega][m_index] = std::complex<double>(0.0,0.0);
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
      GLV::ofs_error << "Error parameter of axis_flag" << std::endl;
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
      GLV::ofs_error << "Error parameter of axis_flag" << std::endl;
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
      GLV::ofs_error << "Error parameter of axis_flag" << std::endl;
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
      GLV::ofs_error << "Error parameter of axis_flag" << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }

  std::vector<std::vector<std::vector<std::complex<double>>>>&
        self_energy::lattice_sigma(const int axis_flag)
  {
    switch(axis_flag)
    {
    case 0:
      return this->sigma_imag.lattice_sigma_access();
      break;
    case 1:
      return this->sigma_real.lattice_sigma_access();
      break;
    default:
      GLV::ofs_error << "Error parameter of axis_flag" << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }

  void self_energy::evalute_lattice_sigma_test(
        const int impurity_solver, const int mag, 
        const int nspin, const std::vector<int>& wbands,
        DFT_output::atoms_info& atom, 
        const std::vector<std::vector<std::vector<
        std::complex<double>>>>&  projector,
        std::vector<std::vector<std::vector<
        std::vector<std::complex<double>>>>>& sigma_loc,
        std::vector<std::vector<std::vector<
        std::complex<double>>>>& lattice_sigma_tmp)
  {
    debug::codestamp("self_energy::evalute_lattice_sigma");

    const int natom=atom.total_atoms();
    const std::vector<int>& norb_sub = atom.iatom_norb();
    const int omega_num=this->nomega(impurity_solver);
    const std::complex<double> zero = std::complex<double>(0.0,0.0);
    const std::complex<double> alpha(1.0,0.0), beta(0.0,0.0);

    std::vector<std::vector<std::vector<std::complex<double>>>>&
          latt_sigma=this->lattice_sigma(impurity_solver);
    
    //======================================================
    //               Allocation and setting zeros
    //======================================================
    if(latt_sigma.empty())
    {
      latt_sigma.resize(nspin);
      for(int ispin=0; ispin<nspin; ispin++)
      {
        latt_sigma[ispin].resize(omega_num);
        for(int i_omega=0; i_omega<omega_num; i_omega++)
        {
          latt_sigma[ispin][i_omega].resize(wbands[ispin]*wbands[ispin], zero);
        }
      }
    }
    else
    {
      #pragma omp parallel for    //expensive
      for(int index=0; index<omega_num*nspin; index++) //omega and spin
      {
        int ispin = index/omega_num;
        int iomega = index%omega_num;

        auto& la=latt_sigma[ispin][iomega];
        for(int band_index=0; band_index<wbands[ispin]*wbands[ispin]; band_index++)
          la[band_index] = zero;
      }
    }

    for(int ineq=0; ineq<atom.inequ_atoms(); ineq++)
    {
      const int iatom = atom.ineq_iatom(ineq);
      const int m_tot=norb_sub[iatom];

      const int mkl_threads = mkl_get_max_threads();
      mkl_set_num_threads(1);      //set the number of threads of MKL library function to 1

      #pragma omp parallel
      {
        int nbands = wbands[0];
        if(nspin==2) nbands = wbands[0] > wbands[1] ? wbands[0] : wbands[1];
        std::unique_ptr<std::complex<double>[]> mat_tmp(new std::complex<double> [nbands*m_tot]);

        #pragma omp for
        for(int index=0; index<omega_num*nspin; index++)
        {
          int ispin = index/omega_num;
          int iomega = index%omega_num;

          // void cblas_zgemm (const CBLAS_LAYOUT Layout, const CBLAS_TRANSPOSE transa, const
          //      CBLAS_TRANSPOSE transb, const MKL_INT m, const MKL_INT n, const MKL_INT k, const void
          //      *alpha, const void *a, const MKL_INT lda, const void *b, const MKL_INT ldb, const void
          //      *beta, void *c, const MKL_INT ldc);

          cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                  wbands[ispin], m_tot, m_tot,
                  &alpha,
                  &projector[ineq][ispin][0], m_tot,
                  &sigma_loc[ineq][ispin][iomega][0], m_tot,
                  &beta,
                  &mat_tmp[0], m_tot);
          
          cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasConjTrans,
                  wbands[ispin], wbands[ispin], m_tot,
                  &alpha,
                  &mat_tmp[0], m_tot,
                  &projector[ineq][ispin][0], m_tot,
                  &beta,
                  &lattice_sigma_tmp[ispin][iomega][0], wbands[ispin]);
        }
      } 

      mkl_set_num_threads(mkl_threads);

      for(int iatom1=0; iatom1<atom.total_atoms(); iatom1++)
      {
        if(iatom1==iatom)
        {
          #pragma omp parallel for    //expensive
          for(int index=0; index<omega_num*nspin; index++)
          {
            int ispin = index/omega_num;
            int iomega = index%omega_num;

            auto& la=latt_sigma[ispin][iomega];
            auto& la1=lattice_sigma_tmp[ispin][iomega]; 

            for(int band_index=0; band_index<wbands[ispin]*wbands[ispin]; band_index++)
              la[band_index] += la1[band_index];
          }
        }
        else if(iatom1!=iatom && atom.equ_atom(iatom1)==ineq)
        {
          int symm=0;            //symmetry
          if(mag==1)             //AFM
          {
            if(atom.magnetic(iatom)*atom.magnetic(iatom1)==-1) symm=1;  //anti_symmetry
          }

          #pragma omp parallel for    //expensive
          for(int index=0; index<omega_num*nspin; index++)
          {
            int ispin = index/omega_num;
            int iomega = index%omega_num;

            auto& la=latt_sigma[ispin][iomega];
            auto& la1=lattice_sigma_tmp[(ispin+symm)%2][iomega]; 
               
            for(int band_index=0; band_index<wbands[ispin]*wbands[ispin]; band_index++)
              la[band_index] += la1[band_index];
          }
        }
      }//iatom

    }//ineq

    return;
  }

}
