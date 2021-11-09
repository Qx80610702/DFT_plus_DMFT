#include "spectrum.h"
#include "../debug.h"
#include "../mpi_environment.h"
#include "math_zone.h"
#include "../constants.h"

#define MKL_Complex16 std::complex<double>
#include <mkl.h>

#include <mpi.h>
#include <omp.h>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <cstring>

void spectrum::eva_spectrum(
      const double mu,
      DFT_output::KS_bands& band, 
      DFT_output::atoms_info& atom, 
      DFT_plus_DMFT::projector& proj,
      DMFT::self_energy& sigma,
      DMFT::input_info& in,
      DFT_plus_DMFT::Hilbert_space& space )
{
  debug::codestamp("spectrum::eva_spectrum");

  const int nomega = sigma.sigma_real.nomega();
  const int nks=band.nk();
  const int nspin=band.nspins();
  const int mag = *(int*)in.parameter("magnetism");
  const std::vector<int>& wbands=space.Wbands();
  const int ineq_num = atom.inequ_atoms();
  const std::vector<int>& norb_sub = atom.iatom_norb();
  const std::vector<double>& fk = band.kweight();

  std::vector<int> k_map;
  for(int ik=0; ik<nks; ik++)
    if(ik%mpi_ntasks() == mpi_rank()) k_map.push_back(ik);

  //===========Allocation===============
  this->freq = sigma.sigma_real.frequency();
  this->Awk.resize(nspin);
  for(int is=0; is<nspin; is++)
  {
    this->Awk[is].resize(nomega);
    for(int iomega=0; iomega<nomega; iomega++)
      this->Awk[is][iomega].resize(nks,0.0);
  }

  this->DOS.resize(nspin);
  for(int is=0; is<nspin; is++)
    this->DOS[is].resize(nomega, 0.0);

  this->Aw_loc.resize(ineq_num);
  for(int ineq=0; ineq<ineq_num; ineq++)
  {
    const int iatom = atom.ineq_iatom(ineq);
    const int m_tot=norb_sub[iatom];

    this->Aw_loc[ineq].resize(nspin);
    for(int is=0; is<nspin; is++)
    {
      this->Aw_loc[ineq][is].resize(nomega);
      for(int iomega=0; iomega<nomega; iomega++)
        this->Aw_loc[ineq][is][iomega].resize(m_tot, 0.0);
    }
  }

  //===========evaluation===============
  for(int ik=0; ik<k_map.size(); ik++)
  {
    const int i_k_point = k_map[ik];

    sigma.evalute_lattice_sigma(
        1, mag, nspin, wbands, atom, 
        proj.proj_access(ik) );

    const std::vector<std::vector<std::vector<std::complex<double>>>>&
          latt_sigma = sigma.lattice_sigma(1);
 
    for(int is=0; is<nspin; is++)
    {
      const auto& epsilon=space.eigen_val()[is][i_k_point];

      std::vector<std::vector<std::complex<double>>> KS_Gw(nomega);
      for(int iomega=0; iomega<nomega; iomega++)
        KS_Gw[iomega].resize(wbands[is]*wbands[is]);
 
      const int mkl_threads = mkl_get_max_threads();
      mkl_set_num_threads(1);  //set the number of threads of MKL library function to 1
      #pragma omp parallel
      {
        int* ipiv = new int [wbands[is]];
        int info_trf, info_tri;

        #pragma omp for
        for(int iomega=0; iomega<nomega; iomega++)
        {        
          for(int iband1=0; iband1<wbands[is]; iband1++)
          {
            for(int iband2=0; iband2<wbands[is]; iband2++)
            {           
              if(iband1==iband2)
                KS_Gw[iomega][iband1*wbands[is]+iband2] = this->freq[iomega] + mu 
                  -epsilon[iband1]-latt_sigma[is][iomega][iband1*wbands[is]+iband2];
              else
                KS_Gw[iomega][iband1*wbands[is]+iband2] = 
                  -latt_sigma[is][iomega][iband1*wbands[is]+iband2];
            }
          }

          general_complex_matrix_inverse(&KS_Gw[iomega][0], wbands[is], &ipiv[0], info_trf, info_tri);
          
          for(int iband=0; iband<wbands[is]; iband++)
            this->Awk[is][iomega][i_k_point] -= KS_Gw[iomega][iband*wbands[is]+iband].imag()/PI;

          for(int ineq=0; ineq<ineq_num; ineq++)
          {
            const int iatom = atom.ineq_iatom(ineq);
            const int m_tot = norb_sub[iatom];

            const std::vector<std::complex<double>>& projector = proj.proj_access(ik)[iatom][is];

            for(int m=0; m<m_tot; m++)
            {
              const int m_index = m*m_tot+m;
              for(int iband1=0; iband1<wbands[is]; iband1++)
              {
                int index1 = iband1*m_tot + m;
                for(int iband2=0; iband2<wbands[is]; iband2++)
                {
                  int index2 = iband2*m_tot + m;
                  this->Aw_loc[ineq][is][iomega][m] -= 
                      ( std::conj(projector[index1])*
                      KS_Gw[iomega][iband1*wbands[is]+iband2]
                      *projector[index2]*fk[i_k_point] ).imag()/PI;
                }//iband2
              }//iband1
            }//m
          }//ineq
        }//iomega

        delete [] ipiv;
      }
      mkl_set_num_threads(mkl_threads);
    }//is
  }//ik

  for(int is=0; is<nspin; is++)
  {
    for(int iomega=0; iomega<nomega; iomega++)
    {
      std::vector<double> array = this->Awk[is][iomega];
      MPI_Allreduce(&array[0], &this->Awk[is][iomega][0], nks,
                    MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
    }
  }
  
  for(int is=0; is<nspin; is++)
    for(int iomega=0; iomega<nomega; iomega++)
      for(int ik=0; ik<nks; ik++)
        this->DOS[is][iomega] += this->Awk[is][iomega][ik]*fk[ik];

  for(int ineq=0; ineq<ineq_num; ineq++)
  {
    const int iatom = atom.ineq_iatom(ineq);
    const int m_tot=norb_sub[iatom];
    const int local_symmetry = atom.local_sym();
    const int corr_L=atom.L(ineq);

    for(int is=0; is<nspin; is++)
    {
      for(int iomega=0; iomega<nomega; iomega++)
      {
        std::vector<double> array = this->Aw_loc[ineq][is][iomega];
        
        MPI_Allreduce(&array[0], &this->Aw_loc[ineq][is][iomega][0],
                      m_tot, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
        
        DFT_output::atoms_info::symmetry_operation_vector<double>(
                      local_symmetry, corr_L, m_tot, 
                      &this->Aw_loc[ineq][is][iomega][0] );
      }
    }
  }

  return;
}

void spectrum::out_spectrum()
{
  debug::codestamp("spectrum::out_spectrum");

  for(int is=0; is<this->Awk.size(); is++)
  {
    std::stringstream ss;
    ss << "Awk_spin" << is << ".dat";
    std::ofstream ofs(ss.str().c_str(), std::ios::out);

    for(std::vector<std::vector<double>>::iterator iter1=this->Awk[is].begin(); iter1!=this->Awk[is].end(); iter1++)
    {
      for(double iter2 : *iter1)
        ofs << std::fixed << std::setprecision(6) << iter2 << " ";
      ofs << std::endl;
    }
    ofs.close();
  }

  for(int is=0; is<this->Awk.size(); is++)
  {
    std::stringstream ss;
    ss << "DOS_spin" << is << ".dat";
    std::ofstream ofs(ss.str().c_str(), std::ios::out);

    for(int iomega=0; iomega<this->DOS[is].size(); iomega++)
    {
      double omega=this->freq[iomega]*Hartree_to_eV;
      
      ofs << std::fixed << std::setprecision(6) << omega << " "
          << std::fixed << std::setprecision(6) << this->DOS[is][iomega] << std::endl ;
    }
    
    ofs.close();
  }

  for(int ineq=0; ineq<this->Aw_loc.size(); ineq++)
  {
    for(int is=0; is<this->Aw_loc[ineq].size(); is++)
    {
      std::stringstream ss;
      ss << "Aw_imp" << ineq <<"_spin" << is << ".dat";
      std::ofstream ofs(ss.str().c_str(), std::ios::out);

      for(int iomega=0; iomega<this->Aw_loc[ineq][is].size(); iomega++)
      {
        double omega=this->freq[iomega]*Hartree_to_eV;
        
        ofs << std::fixed << std::setprecision(6) << omega << " ";
        for(double iter : this->Aw_loc[ineq][is][iomega])
          ofs << std::fixed << std::setprecision(6) << iter << " ";

        ofs << std::endl;
      }
      ofs.close();
    }
  }

  return;
}
