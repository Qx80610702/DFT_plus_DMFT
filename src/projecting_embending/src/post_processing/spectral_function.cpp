#include "../debug/debug.h"
#include "spectral_function.h"
#include "../constants.h"

#include <string>
#include <sstream>
#include <fstream>
#include <iomanip>

void spectral_function::evaluate_local_spectrum(
        DMFT::input_info& in,
        DFT_output::atoms_info& atom,
        DFT_output::KS_bands& band,
        std::vector<std::vector<std::vector<
        std::vector<std::complex<double>>>>>& Gf,
        std::vector<std::vector<std::vector<
        std::vector<std::complex<double>>>>>& Gf_save,
        std::vector<double>& freq )
{
  debug::codestamp("spectral_function::evaluate_local_spectrum");

  const int ineq_num = atom.inequ_atoms();
  const bool soc = band.soc();
  const int nspin = band.nspins();
  const int nomega = freq.size();

  //Create directory post_processing 
  std::string post_processing = "post_processing";
  std::stringstream make_dir1;
  make_dir1 << "test -d " << post_processing << " || mkdir " << post_processing;
  system(make_dir1.str().c_str());

  //Create directory spectrum
  std::string spectrum_dir="/spectrum";
  std::stringstream make_dir2;
  make_dir2 << "test -d " << post_processing << spectrum_dir
          << " || mkdir " << post_processing << spectrum_dir;
  system(make_dir2.str().c_str());

  for(int ineq=0; ineq<ineq_num; ineq++)
  {
    const int iatom = atom.ineq_iatom(ineq);
    const int angular_L = atom.L(iatom);
    const int m_tot = 2*angular_L+1;

    auto& Gfa = Gf[ineq];

    std::stringstream site_dir_ss;
    site_dir_ss << "/impurity" << ineq;
    std::string site_dir= site_dir_ss.str();
    std::stringstream make_dir3;
    make_dir3 << "test -d " << post_processing << spectrum_dir << site_dir
            << " || mkdir " << post_processing << spectrum_dir << site_dir;
    system(make_dir3.str().c_str());

    for(int is=0; is<nspin; is++)
    {
      auto& Gfb = Gfa[is];

      std::stringstream spin_dir_ss;
      spin_dir_ss << "/spin" << is;
      std::string spin_dir= spin_dir_ss.str();

      std::stringstream make_dir4;
      make_dir4 << "test -d " << post_processing << spectrum_dir << site_dir << spin_dir 
              << " || mkdir " << post_processing << spectrum_dir << site_dir << spin_dir;
      system(make_dir4.str().c_str());

      for(int m_index=0; m_index<m_tot; m_index++)
      {
        std::stringstream m_dir_ss;
        m_dir_ss << "/m" << m_index;
        std::string m_dir= m_dir_ss.str();
        std::stringstream make_dir5;
        make_dir5 << "test -d " << post_processing << spectrum_dir << site_dir << spin_dir << m_dir
                << " || mkdir " << post_processing << spectrum_dir << site_dir << spin_dir << m_dir;
        system(make_dir5.str().c_str());

        std::stringstream current_dir_ss;
        current_dir_ss << post_processing << spectrum_dir << site_dir << spin_dir << m_dir;
        std::string current_dir = current_dir_ss.str();

        //Green functioin
        std::string Gf_file = current_dir+"/G_omega.in";
        std::ofstream ofs_gf(Gf_file.c_str(), std::ios::out);
        for(int iomega=0; iomega<nomega; iomega++)
        {
          ofs_gf << std::left << std::setw(20) << std::fixed << std::setprecision(12) << freq[iomega]*Hartree_to_eV
          << std::setw(18) << std::fixed << std::setprecision(12) << Gfb[iomega][m_index*m_tot+m_index].real()/Hartree_to_eV
          << " 1.0e-5 "
          << std::setw(18) << std::fixed << std::setprecision(12) << Gfb[iomega][m_index*m_tot+m_index].imag()/Hartree_to_eV
          << " 1.0e-5 \n";
        }
        ofs_gf.close();

        //in.param
        std::string param_file = current_dir+"/in.param";
        std::ofstream ofs_param(param_file.c_str(), std::ios::out);

        ofs_param << std::left << "BETA=" << std::fixed << std::setprecision(9)
                  << (*(double*)in.parameter("beta"))/Hartree_to_eV << '\n';
        ofs_param << std::left << "NDAT=" << 2*nomega << '\n';
        ofs_param << std::left << "NFREQ=" << nomega << '\n';
        ofs_param << std::left << "DATASPACE=frequency" << '\n';
        ofs_param << std::left << "KERNEL=fermionic" << '\n';
        ofs_param << std::left << "PARTICLE_HOLE_SYMMETRY=false" << '\n';
        ofs_param << std::left << "DATA=\"G_omega.in\"" << '\n';
        ofs_param << std::left << "NORM=.75643" << '\n';

        ofs_param.close();
      }//m
    }//is

    //Input Green's function
    // auto& Gf_savea = Gf_save[ineq];
    // for(int is=0; is<nspin; is++)
    // {
    //   auto& Gf_saveb = Gf_savea[is];

    //   std::stringstream spin_dir_ss;
    //   spin_dir_ss << "/save_spin" << is;
    //   std::string spin_dir= spin_dir_ss.str();

    //   std::stringstream make_dir4;
    //   make_dir4 << "test -d " << post_processing << spectrum_dir << site_dir << spin_dir 
    //           << " || mkdir " << post_processing << spectrum_dir << site_dir << spin_dir;
    //   system(make_dir4.str().c_str());

    //   for(int m_index=0; m_index<m_tot; m_index++)
    //   {
    //     std::stringstream m_dir_ss;
    //     m_dir_ss << "/m" << m_index;
    //     std::string m_dir= m_dir_ss.str();
    //     std::stringstream make_dir5;
    //     make_dir5 << "test -d " << post_processing << spectrum_dir << site_dir << spin_dir << m_dir
    //             << " || mkdir " << post_processing << spectrum_dir << site_dir << spin_dir << m_dir;
    //     system(make_dir5.str().c_str());

    //     std::stringstream current_dir_ss;
    //     current_dir_ss << post_processing << spectrum_dir << site_dir << spin_dir << m_dir;
    //     std::string current_dir = current_dir_ss.str();

    //     //Green functioin
    //     std::string Gf_file = current_dir+"/G_omega.in";
    //     std::ofstream ofs_gf(Gf_file.c_str(), std::ios::out);
    //     for(int iomega=0; iomega<nomega; iomega++)
    //     {
    //       ofs_gf << std::left << std::setw(20) << std::fixed << std::setprecision(12) << freq[iomega]*Hartree_to_eV
    //       << std::setw(18) << std::fixed << std::setprecision(12) << Gf_saveb[iomega][m_index*m_tot+m_index].real()/Hartree_to_eV
    //       << " 1.0e-5 "
    //       << std::setw(18) << std::fixed << std::setprecision(12) << Gf_saveb[iomega][m_index*m_tot+m_index].imag()/Hartree_to_eV
    //       << " 1.0e-5 \n";
    //     }
    //     ofs_gf.close();

    //     //in.param
    //     std::string param_file = current_dir+"/in.param";
    //     std::ofstream ofs_param(param_file.c_str(), std::ios::out);

    //     ofs_param << std::left << "BETA=" << std::fixed << std::setprecision(9)
    //               << (*(double*)in.parameter("beta"))/Hartree_to_eV << '\n';
    //     ofs_param << std::left << "NDAT=" << 2*nomega << '\n';
    //     ofs_param << std::left << "NFREQ=" << nomega << '\n';
    //     ofs_param << std::left << "DATASPACE=frequency" << '\n';
    //     ofs_param << std::left << "KERNEL=fermionic" << '\n';
    //     ofs_param << std::left << "PARTICLE_HOLE_SYMMETRY=false" << '\n';
    //     ofs_param << std::left << "DATA=\"G_omega.in\"" << '\n';
    //     ofs_param << std::left << "NORM=.75643" << '\n';

    //     ofs_param.close();
    //   }//m
    // }//is

  }//ineq

  return;
}