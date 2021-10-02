#include "double_counting.h"
#include "../mpi_environment.h"
#include "../debug/debug.h"
#include "../constants.h"

#include <iomanip>
#include <cstdlib>

namespace DMFT
{
  void double_counting::cal_double_counting(
      const int type, const bool SOC, const int nspin, 
      DFT_output::atoms_info& atom,
      DMFT::input_info& in )
  {
    debug::codestamp("double_counting::cal_double_counting");

    const std::vector<int>& norb_sub = atom.iatom_norb();
    const std::complex<double> zero(0.0,0.0);

    this->V_dc.resize(atom.inequ_atoms());

    if(type==1)
    {
      if(mpi_rank()==0) std::cout << "\n=============Double counting(FLL)============" << std::endl;

      for(int ineq=0; ineq<atom.inequ_atoms(); ineq++)
      {
        const int iatom=atom.ineq_iatom(ineq);
        const int m_tot=norb_sub[iatom];

        int nspin_tmp=nspin;
        if(SOC) nspin_tmp=4;

        const double Norb = 2*m_tot;
        const double averU = ( atom.Uval(iatom) + (atom.Uval(iatom)-2.0*atom.Jval(iatom))
                              *(Norb-2) )/(Norb-1.0);

        const double averJ = atom.Jval(iatom);

        this->V_dc.at(ineq).resize(nspin_tmp);
        for(int is=0; is<nspin_tmp; is++)
        {
          this->V_dc.at(ineq).at(is).resize(m_tot*m_tot,zero);
          const double tmp = averU*(atom.occ_num(iatom,0)+atom.occ_num(iatom,1)-0.5)
                            -averJ*(atom.occ_num(iatom,is)-0.5);

          for(int m=0; m<m_tot; m++)
            this->V_dc.at(ineq).at(is).at(m*m_tot+m) = std::complex<double>(tmp,0.0); //*atom.occ_num_m()[iatom][is][m];

          if(mpi_rank()==0)
          {
            std::cout << "===========impurity:" << ineq << " spin:" << is << "========\n";
            
            for(int m=0; m<m_tot; m++)
              std::cout << std::setw(12) << std::fixed << std::setprecision(6) << 
                Hartree_to_eV*this->V_dc[ineq][is][m*m_tot+m].real();

            std::cout << std::endl;
          }
        }
      }//ineq

      // for(int ineq=0; ineq<atom.inequ_atoms(); ineq++)
      // {
      //   const int iatom=atom.ineq_iatom(ineq);
      //   const int angular_L=atom.L(iatom);
      //   const int m_tot = 2*angular_L+1;

      //   int nspin_tmp=nspin;
      //   if(SOC) nspin_tmp=4;

      //   std::vector<double> Nelec(2,0.0);
      //   for(int is=0; is<2; is++)
      //   {
      //     for(int m=0; m<m_tot; m++)
      //     {
      //       if(m==2 || m==4) continue;
      //       Nelec[is] += atom.occ_num_m()[iatom][is][m];
      //     }
      //   }

      //   this->V_dc.at(ineq).resize(nspin_tmp);
      //   for(int is=0; is<nspin_tmp; is++)
      //   {
      //     this->V_dc.at(ineq).at(is).resize(m_tot*m_tot,zero);
      //     const double tmp = atom.Uval(iatom)*(Nelec[0]+Nelec[1]-0.5)
      //                       -atom.Jval(iatom)*(Nelec[is]-0.5) ;
      //     for(int m=0; m<m_tot; m++)
      //     {
      //       if(m==2 || m==4) continue;
      //       this->V_dc.at(ineq).at(is).at(m*m_tot+m) = std::complex<double>(tmp,0.0)*atom.occ_num_m()[iatom][is][m];
      //     }
      //   }
      // }//ineq
    }//type==1

    return;
  }

}
