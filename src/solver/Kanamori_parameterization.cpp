#include "Kanamori_parameterization.h"
#include "../constants.h"
#include "../debug.h"

#include <fstream>
#include <iomanip>
#include <string>
#include <sstream>

namespace DMFT
{
  void Kanamori_parameterization::evaluate_coulomb_tensor(
            DFT_output::atoms_info& atom,
            std::vector<std::vector<std::vector<std::vector<
            std::vector<double>>>>>& U_matrix)
  {
    debug::codestamp("Kanamori_parameterization::evaluate_coulomb_tensor");

    const std::vector<int>& norb_sub = atom.iatom_norb();

    //Allocation
    U_matrix.resize(atom.inequ_atoms());
    for(int ineq=0; ineq<atom.inequ_atoms(); ineq++)
    {
      const int iatom = atom.ineq_iatom(ineq);
      const int m_tot = norb_sub[iatom];

      U_matrix[ineq].resize(2*m_tot);
      for(int m1=0; m1<2*m_tot; m1++)
      {
        U_matrix[ineq][m1].resize(2*m_tot);
        for(int m2=0; m2<2*m_tot; m2++)
        {
          U_matrix[ineq][m1][m2].resize(2*m_tot);
          for(int m3=0; m3<2*m_tot; m3++)
          {
            U_matrix[ineq][m1][m2][m3].resize(2*m_tot,0.0);
          }
        }
      }
    }

    //evaluation
    for(int ineq=0; ineq<atom.inequ_atoms(); ineq++)
    {
      const int iatom = atom.ineq_iatom(ineq);
      const int m_tot=norb_sub[iatom];
      const double U = atom.Uval(iatom);
      const double J = atom.Jval(iatom);

      int count1 = 0;
      for(int iorb1=0; iorb1<2*m_tot; iorb1++)
      {
        int is1 = iorb1/m_tot;
        int m1 = iorb1<m_tot ? iorb1 : (iorb1-m_tot);

        for(int iorb2=0; iorb2<2*m_tot; iorb2++)
        {
          int is2 = iorb2/m_tot;
          int m2 = iorb2<m_tot ? iorb2 : (iorb2-m_tot);

          if(is1==is2)
          {
            if(m1==m2) U_matrix[ineq][iorb1][iorb2][iorb2][iorb1] = 0.0;
            else U_matrix[ineq][iorb1][iorb2][iorb2][iorb1] = U-3.0*J;
          }
          else
          {
            if(m1==m2) U_matrix[ineq][iorb1][iorb2][iorb2][iorb1] = U;
            else U_matrix[ineq][iorb1][iorb2][iorb2][iorb1] = U-2.0*J;
          }
        }
      }
    }

    return;
  }

  void Kanamori_parameterization::out_ALPS_CTHYB(
                      const int char_step,
                      const int DMFT_step, 
                      DFT_output::atoms_info& atom,
                      DFT_output::KS_bands& band)
  {
    debug::codestamp("Kanamori_parameterization::out_ALPS_CT_HYB");

    const int ineq_num = atom.inequ_atoms();
    const std::vector<int>& norb_sub = atom.iatom_norb();
    const int nspin = band.nspins();

    std::string dir_impurity_solving = "dmft";

    std::stringstream char_dir_ss;
    char_dir_ss << "/charge_step" << char_step;
    std::string char_step_dir= char_dir_ss.str();

    std::stringstream step_dir_ss;
    step_dir_ss << "/dmft_step" << DMFT_step;
    std::string step_dir= step_dir_ss.str();

    for(int ineq=0; ineq<ineq_num; ineq++)
    {
      const int iatom = atom.ineq_iatom(ineq);
      const int m_tot=norb_sub[iatom];

      std::stringstream site_dir_ss;
      site_dir_ss << "/impurity" << ineq;
      std::string site_dir= site_dir_ss.str();

      std::stringstream current_dir_ss;
      current_dir_ss << dir_impurity_solving << char_step_dir << step_dir << site_dir;
      std::string current_dir = current_dir_ss.str();

      std::string Uijkl = current_dir+"/Uijkl.txt";
      std::ofstream ofs(Uijkl.c_str(), std::ios::out);

      ofs << 2*m_tot*(2*m_tot-1) << '\n';

      const double zero = 0.0;
      int count=0;
      int count1=0;
      for(int m1=0; m1<m_tot; m1++)
      {
        for(int is1=0; is1<2; is1++)
        {
          int count2=0;
          for(int m2=0; m2<m_tot; m2++)
          {
            for(int is2=0; is2<2; is2++)
            {
              if(is1!=is2)
                if(m1==m2)
                {
                  ofs << std::left << std::setw(6) << count 
                  << std::setw(3) << count1 << std::setw(3) << count2
                  << std::setw(3) << count2 << std::setw(3) << count1
                  << std::setw(20) << std::fixed << std::setprecision(12) 
                  << GlobalC::Hartree_to_eV*atom.Uval(iatom)
                  << std::setw(20) << std::fixed << std::setprecision(12) << zero << '\n';
                  count++;
                }
                else
                {
                  ofs << std::left << std::setw(6) << count 
                    << std::setw(3) << count1 << std::setw(3) << count2
                    << std::setw(3) << count2 << std::setw(3) << count1
                    << std::setw(20) << std::fixed << std::setprecision(12) 
                    << GlobalC::Hartree_to_eV*(atom.Uval(iatom)-2*atom.Jval(iatom))
                    << std::setw(20) << std::fixed << std::setprecision(12) << zero << '\n';
                  count++;
                }
              else
                if(m1!=m2)
                {
                  ofs << std::left << std::setw(6) << count 
                    << std::setw(3) << count1 << std::setw(3) << count2
                    << std::setw(3) << count2 << std::setw(3) << count1
                    << std::setw(20) << std::fixed << std::setprecision(12) 
                    << GlobalC::Hartree_to_eV*(atom.Uval(iatom)-3*atom.Jval(iatom))
                    << std::setw(20) << std::fixed << std::setprecision(12) << zero << '\n';
                    count++;
                }
              count2++;
            }//is2
          }//m2
          count1++;
        }//is1
      }//m1
      ofs.close();

    }//ineq

    return;
  }

  void Kanamori_parameterization::out_ALPS_CTHYB_SEGMENT(
            const int char_step,
            const int DMFT_step, 
            DFT_output::atoms_info& atom,
            DFT_output::KS_bands& band )
  {
    debug::codestamp("Kanamori_parameterization::out_ALPS_CTHYB_SEGMENT");

    const int ineq_num = atom.inequ_atoms();
    const std::vector<int>& norb_sub = atom.iatom_norb();
    const int nspin = band.nspins();

    std::string dir_impurity_solving = "dmft";

    std::stringstream char_dir_ss;
    char_dir_ss << "/charge_step" << char_step;
    std::string char_step_dir= char_dir_ss.str();

    std::stringstream step_dir_ss;
    step_dir_ss << "/dmft_step" << DMFT_step;
    std::string step_dir= step_dir_ss.str();

    for(int ineq=0; ineq<ineq_num; ineq++)
    {
      const int iatom = atom.ineq_iatom(ineq);
      const int m_tot=norb_sub[iatom];

      std::stringstream site_dir_ss;
      site_dir_ss << "/impurity" << ineq;
      std::string site_dir= site_dir_ss.str();

      std::stringstream current_dir_ss;
      current_dir_ss << dir_impurity_solving << char_step_dir << step_dir << site_dir;
      std::string current_dir = current_dir_ss.str();

      std::string Uij = current_dir+"/Uij.dat";
      std::ofstream ofs(Uij.c_str(), std::ios::out);

      const double zero = 0.0;
      for(int m1=0; m1<m_tot; m1++)
      {
        for(int is1=0; is1<2; is1++)
        {
          for(int m2=0; m2<m_tot; m2++)
          {
            for(int is2=0; is2<2; is2++)
            {
              if(is1==is2)
              {
                if(m1==m2)
                  ofs << std::left << std::setw(10) << std::fixed 
                      << std::setprecision(6) << zero;
                else
                  ofs << std::left << std::setw(10) << std::fixed 
                      << std::setprecision(6) << GlobalC::Hartree_to_eV*
                        (atom.Uval(iatom) - 3*atom.Jval(iatom));
              }//is==is2
              else
              {
                if(m1==m2)
                  ofs << std::left << std::setw(10) << std::fixed 
                      << std::setprecision(6) << GlobalC::Hartree_to_eV*atom.Uval(iatom);
                else
                  ofs << std::left << std::setw(10) << std::fixed 
                      << std::setprecision(6) << GlobalC::Hartree_to_eV*
                          (atom.Uval(iatom) - 2*atom.Jval(iatom));
              }//is1!=is2
            }//is2
          }//m2
          ofs << '\n';
        }//is1
      }//m1
      ofs.close();

    }//ineq

    return;
  }

}
