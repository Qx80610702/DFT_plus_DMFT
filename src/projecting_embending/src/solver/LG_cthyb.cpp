#include "LG_cthyb.h"

#include "../debug.h"
#include "../timer.h"
#include "math_zone.h"
#include "../constants.h"
#include "math_zone.h"

#include <omp.h>
#include <memory>
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>
#include <cstdlib>   //Use exit function

namespace DMFT
{
  void LG_CTHYB::output(const int istep, 
        const double mu, DMFT::input_info& in, 
        DFT_output::atoms_info& atom, DFT_output::KS_bands& band,
        std::vector<std::vector<std::vector<std::complex<double>>>>& Eimp,
        std::vector<std::vector<std::vector<std::vector<std::complex<double>>>>>& Gf_in,
        std::vector<std::vector<std::vector<std::vector<std::complex<double>>>>>& Weiss,
        std::vector<std::vector<std::vector<std::vector<std::complex<double>>>>>& hyb_omega)
  {
    //======================================================
    //             energy unit: eV
    //======================================================
    debug::codestamp("LG_CTHYB::output");

    const int ineq_num = atom.inequ_atoms();
    const std::vector<int>& norb_sub = atom.iatom_norb();
    const bool soc = band.soc();
    const int nspin = band.nspins();
    const int ntau = *(int*)in.parameter("n_tau");
    const int nomega = *(int*)in.parameter("n_omega");
    double zero=0.0;

    //Create directory impurity 
    std::string dir_impurity_solving = "impurity_solving";
    std::stringstream make_dir1;
    make_dir1 << "test -d " << dir_impurity_solving << " || mkdir " << dir_impurity_solving;
    system(make_dir1.str().c_str());

    //Create directory step+num
    std::stringstream step_dir_ss;
    step_dir_ss << "/step" << istep;
    std::string step_dir= step_dir_ss.str();
    std::stringstream make_dir2;
    make_dir2 << "test -d " << dir_impurity_solving << step_dir
            << " || mkdir " << dir_impurity_solving << step_dir;
    system(make_dir2.str().c_str());

    for(int ineq=0; ineq<ineq_num; ineq++)
    {
      const int iatom = atom.ineq_iatom(ineq);
      const int m_tot = norb_sub[iatom];
      
      const std::vector<std::vector<std::complex<double>>>&
            hoppinga = Eimp[ineq];
      
      const std::vector<std::vector<std::vector<std::complex<double>>>>&
            Gf_ina = Gf_in[ineq];

      const std::vector<std::vector<std::vector<std::complex<double>>>>& 
            hyb = hyb_omega[ineq];
      
      const std::vector<std::vector<std::vector<std::complex<double>>>>& 
            G0 = Weiss[ineq];

      const std::vector<std::vector<std::vector<std::complex<double>>>>& 
            Gw = Gf_in[ineq];

      const std::vector<std::vector<std::vector<std::complex<double>>>>& 
            hybt = this->hyb_tau[ineq];

      std::stringstream site_dir_ss;
      site_dir_ss << "/impurity" << ineq;
      std::string site_dir= site_dir_ss.str();
      std::stringstream make_dir3;
      make_dir3 << "test -d " << dir_impurity_solving << step_dir << site_dir
              << " || mkdir " << dir_impurity_solving << step_dir << site_dir;
      system(make_dir3.str().c_str());

      std::stringstream current_dir_ss;
      current_dir_ss << dir_impurity_solving << step_dir << site_dir;
      std::string current_dir = current_dir_ss.str();

      //===========================================
      //           write  control.in
      //===========================================
      std::string input_file = current_dir+"/input.dat";
      std::ofstream ofs_input(input_file.c_str(), std::ios::out);

      ofs_input << "&Model_Parameter\n";
      ofs_input << "       nSpin = 2\n";
      ofs_input << "      nOrbit = " << m_tot << '\n';
      ofs_input << "        Beta = " << std::fixed << std::setprecision(9) 
                << *(double*)in.parameter("beta")/Hartree_to_eV  << '\n';
      ofs_input << "          xU = " << atom.Uval(iatom)*Hartree_to_eV  << '\n';
      ofs_input << "          xJ = " << atom.Jval(iatom)*Hartree_to_eV  << '\n';
      ofs_input << "Paramagnetic = .True.\n";
      ofs_input << "      nOmega = " << nomega << '\n';
      ofs_input << "        nTau = " << ntau << '\n';
      ofs_input << "  nLegenPoly = 50\n";
      ofs_input << "/\n\n";

      ofs_input << "&MC_Parameter\n";
      ofs_input << "   nDMFT = 1\n";
      ofs_input << "   nWarm = 100000\n";
      ofs_input << "   nMeasure = " << *(long long*)in.parameter("mc_step")/10 << '\n';
      ofs_input << "    nBin = 10\n";
      ofs_input << "/\n";

      ofs_input.close();

      //===========================================
      //           write  delta_omega.dat
      //===========================================
      std::string delta_file = current_dir+"/delta_omega.dat";
      std::ofstream ofs_delta(delta_file.c_str(), std::ios::out);

      for(int iomega=0; iomega<nomega; iomega++)
      {
        // ofs_delta << std::setw(22) << std::fixed << std::setprecision(15)
        //   << Hartree_to_eV*freq[iomega];
        for(int is=0; is<2; is++)
        {
          for(int m=0; m<m_tot; m++)
          {
            if(nspin==1)
              ofs_delta << std::setw(22) << std::fixed << std::setprecision(15) 
              << Hartree_to_eV*hyb[0][iomega][m*m_tot+m].real() 
              << std::setw(22) << std::fixed << std::setprecision(15) 
              << Hartree_to_eV*hyb[0][iomega][m*m_tot+m].imag();
            else
              ofs_delta << std::setw(22) << std::fixed << std::setprecision(15) 
              << Hartree_to_eV*hyb[is][iomega][m*m_tot+m].real() 
              << std::setw(22) << std::fixed << std::setprecision(15) 
              << Hartree_to_eV*hyb[is][iomega][m*m_tot+m].imag();       
          }//m
        }//is
        ofs_delta << '\n';
      }//iomega
      ofs_delta.close();

      //===========================================
      //           write  delta_t.dat
      //===========================================
      std::string deltat_file = current_dir+"/delta_t.dat";
      std::ofstream ofs_deltat(deltat_file.c_str(), std::ios::out);

      for(int itau=0; itau<ntau+1; itau++)
      {
        for(int is=0; is<2; is++)
        {
          for(int m=0; m<m_tot; m++)
          {
            if(nspin==1)
              ofs_deltat << std::left << std::setw(5) << itau
              << std::setw(4) << m+1 << std::setw(4) << is+1
              << std::setw(22) << std::fixed << std::setprecision(15) 
              << std::pow(Hartree_to_eV,2)*hybt[0][itau][m*m_tot+m].real() << '\n';
            else
              ofs_deltat << std::left << std::setw(5) << itau
              << std::setw(4) << m+1 << std::setw(4) << is+1
              << std::setw(22) << std::fixed << std::setprecision(15) 
              << std::pow(Hartree_to_eV,2)*hybt[is][itau][m*m_tot+m].real() << '\n';      
          }//m
        }//is
      }//iomega
      ofs_deltat.close();

      //===========================================
      //           write  mu_vector.dat
      //===========================================
      std::string hopping_file = current_dir+"/mu_vector.dat";
      std::ofstream ofs_hopping(hopping_file.c_str(), std::ios::out);

      int count1=0;
      for(int m=0; m<m_tot; m++)
      {
        for(int is=0; is<2; is++)
        {
          if(nspin==1)
            ofs_hopping << std::setw(4) << m+1 << std::setw(4) << is+1
              << std::setw(22) << std::fixed << std::setprecision(15) 
              << Hartree_to_eV*(mu - hoppinga[0][m*m_tot+m].real()) << '\n';
          else
            ofs_hopping << std::setw(4) << m+1 << std::setw(4) << is+1
              << std::setw(22) << std::fixed << std::setprecision(15) 
              << Hartree_to_eV*(mu - hoppinga[is][m*m_tot+m].real()) << '\n';
        }//is
      }//m
      ofs_hopping.close();

      //===========================================
      //           write  Weiss.dat
      //===========================================
      std::string G0_file = current_dir+"/Weiss.dat";
      std::ofstream ofs_G0(G0_file.c_str(), std::ios::out);

      for(int iomega=0; iomega<nomega; iomega++)
      {
        for(int is=0; is<2; is++)
        {
          // if(m==2 || m==4) continue;
          for(int m=0; m<m_tot; m++)
          {
            if(nspin==1)
              ofs_G0 << std::left << std::setw(5) << iomega+1 
              << std::setw(4) << m+1 << std::setw(4) << is+1
              << std::setw(22) << std::fixed << std::setprecision(15) 
              << G0[0][iomega][m*m_tot+m].real()/Hartree_to_eV << " "
              << std::setw(22) << std::fixed << std::setprecision(15) 
              << G0[0][iomega][m*m_tot+m].imag()/Hartree_to_eV << '\n';
            else
              ofs_G0 << std::left << std::setw(5) << iomega+1
              << std::setw(4) << m+1 << std::setw(4) << is+1
              << std::setw(22) << std::fixed << std::setprecision(15) 
              << G0[is][iomega][m*m_tot+m].real()/Hartree_to_eV << " "
              << std::setw(22) << std::fixed << std::setprecision(15) 
              << G0[is][iomega][m*m_tot+m].imag()/Hartree_to_eV << '\n';       
          }//m
        }//is
      }//iomega
      ofs_G0.close();

      //==================================================
      //   write  Gf.in; the input Green function 
      //   of current step (in Matrsubara frequency), 
      //   which will be read by last step to judge whether
      //   the self-consistency is achieved
      //==================================================
      std::string Gf_file = current_dir+"/Gf.in";
      std::ofstream ofs_gf(Gf_file.c_str(), std::ios::out);

      for(int iomega=0; iomega<nomega; iomega++)
      {
        for(int is=0; is<2; is++)
        {
          for(int m=0; m<m_tot; m++)
          {
            if(nspin==1)
              ofs_gf << std::left << std::setw(5) << iomega+1 
              << std::setw(4) << m+1 << std::setw(4) << is+1
              << std::setw(22) << std::fixed << std::setprecision(15) 
              << Gw[0][iomega][m*m_tot+m].real()/Hartree_to_eV << " "
              << std::setw(22) << std::fixed << std::setprecision(15) 
              << Gw[0][iomega][m*m_tot+m].imag()/Hartree_to_eV << '\n';
            else
              ofs_gf << std::left << std::setw(5) << iomega+1
              << std::setw(4) << m+1 << std::setw(4) << is+1
              << std::setw(22) << std::fixed << std::setprecision(15) 
              << Gw[is][iomega][m*m_tot+m].real()/Hartree_to_eV << " "
              << std::setw(22) << std::fixed << std::setprecision(15) 
              << Gw[is][iomega][m*m_tot+m].imag()/Hartree_to_eV << '\n';       
          }//m
        }//is
      }//iomega
      ofs_gf.close();

    }//ineq

    return;
  }

  void LG_CTHYB::read_last_step(
          const int istep, 
          DFT_output::KS_bands& band,
          DMFT::input_info& in, 
          DFT_output::atoms_info& atom,
          std::vector<std::vector<std::
          vector<std::vector<
          std::complex<double>>>>>& Gf_qmc,
          std::vector<std::vector<std::
          vector<std::vector<
          std::complex<double>>>>>& Weiss,
          std::vector<std::vector<std::
          vector<std::vector<
          std::complex<double>>>>>& Gf_save)
  {
    debug::codestamp("ALPS_CTHYB::read_last_step");

    const int ineq_num = atom.inequ_atoms();
    const std::vector<int>& norb_sub = atom.iatom_norb();
    const bool soc = band.soc();
    const int nspin = band.nspins();

    int iomega, m1, m2;
    double real, imag;

    //directory impurity 
    std::string dir_impurity_solving = "impurity_solving";

    //directory step+num
    std::stringstream step_dir_ss;
    step_dir_ss << "/step" << istep-1;
    std::string step_dir= step_dir_ss.str();

    for(int ineq=0; ineq<ineq_num; ineq++)
    {
      const int iatom = atom.ineq_iatom(ineq);
      const int m_tot = norb_sub[iatom];

      std::stringstream site_dir_ss;
      site_dir_ss << "/impurity" << ineq;
      std::string site_dir= site_dir_ss.str();

      std::vector<std::vector<std::vector<std::complex<double>>>>&
            Gf_savea = Gf_save[ineq];

      std::vector<std::vector<std::vector<std::complex<double>>>>&
            Gf_qmca = Gf_qmc[ineq];

      std::vector<std::vector<std::vector<std::complex<double>>>>&
            Weissa =Weiss[ineq];
                  
      for(int is=0; is<nspin; is++)
      {
        std::stringstream spin_dir_ss;
        spin_dir_ss << "/spin" << is;
        std::string spin_dir= spin_dir_ss.str();

        std::stringstream current_dir_ss;
        current_dir_ss << dir_impurity_solving << step_dir << site_dir << spin_dir;
        std::string current_dir = current_dir_ss.str();

        //=================================================
        //    Read interacting Green function of last step
        //=================================================
        std::vector<std::vector<std::complex<double>>>& 
            Gf_qmcb = Gf_qmca[is];

        std::string Gf_file = current_dir+"/G_omega.dat";
        std::ifstream ifs_gf(Gf_file.c_str(), std::ios::in);

        if (!ifs_gf)  
	      {
	      	std::cout << "Fail to oepn " << Gf_file.c_str() << std::endl;
          std::exit(EXIT_FAILURE);
        }

        ifs_gf.seekg(0);    //set the position at the beginning of the file
        ifs_gf >> iomega;
        ifs_gf.ignore(150, '\n');

        while(ifs_gf.good())
        {
          ifs_gf >> iomega;
          if(ifs_gf.eof()) break; //Check whether end of file is reached 
          ifs_gf >> m1;
          ifs_gf >> m2;
          ifs_gf >> real;
          ifs_gf >> imag;
          ifs_gf.ignore(150,'\n');

          Gf_qmcb[iomega][m1*m_tot+m2] = Hartree_to_eV*std::complex<double>(real,imag);

          if(ifs_gf.eof()) break;//Check whether end of file is reached       
        }
        ifs_gf.close();

        //=============================================
        //    Read input Green function of last step
        //=============================================
        std::vector<std::vector<std::complex<double>>>&
            Gf_saveb = Gf_savea[is];

        std::string Gf_save_file = current_dir+"/Gf.in";
        std::ifstream ifs_gfsave(Gf_save_file.c_str(), std::ios::in);

        if (!ifs_gfsave)  
	      {
	      	std::cout << "Fail to oepn " << Gf_save_file.c_str() << std::endl;
          std::exit(EXIT_FAILURE);
        }

        ifs_gfsave.seekg(0);    //set the position at the beginning of the file
        while(ifs_gfsave.good())
        {
          ifs_gfsave >> iomega;
          if(ifs_gfsave.eof()) break; //Check whether end of file is reached 
          ifs_gfsave >> m1;
          ifs_gfsave >> m2;
          ifs_gfsave >> real;
          ifs_gfsave >> imag;
          ifs_gfsave.ignore(150,'\n');

          Gf_saveb[iomega][m1*m_tot+m2] = std::complex<double>(real,imag);

          if(ifs_gfsave.eof()) break; //Check whether end of file is reached
        }
        ifs_gfsave.close();

        //===================================================
        //    Read input hybridization function of last step
        //===================================================
        std::string hyb_save_file = current_dir+"/hybridization.in";
        std::ifstream ifs_hybsave(hyb_save_file.c_str(), std::ios::in);

        if (!ifs_hybsave)  
	      {
	      	std::cout << "Fail to oepn " << hyb_save_file.c_str() << std::endl;
          std::exit(EXIT_FAILURE);
        }

        ifs_hybsave.seekg(0);    //set the position at the beginning of the file
        while(ifs_hybsave.good())
        {
          ifs_hybsave >> iomega;
          if(ifs_hybsave.eof()) break; //Check whether end of file is reached 
          ifs_hybsave >> m1;
          ifs_hybsave >> m2;
          ifs_hybsave >> real;
          ifs_hybsave >> imag;
          ifs_hybsave.ignore(150,'\n');

          // hyb_saveb[iomega][m1*m_tot+m2] = std::complex<double>(real,imag);

          if(ifs_hybsave.eof()) break; //Check whether end of file is reached       
        }
        ifs_hybsave.close();

      }//is
    }//ineq

    return;
  }

  void LG_CTHYB::out_sigma_last_step(
        const int istep, DFT_output::KS_bands& band,
        DMFT::input_info& in, DFT_output::atoms_info& atom,
        const  std::vector<std::vector<std::vector<
        std::vector<std::complex<double>>>>>& sigma)
  {
    debug::codestamp("LG_CTHYB::out_sigma_last_step");

    const int ineq_num = atom.inequ_atoms();
    const std::vector<int>& norb_sub = atom.iatom_norb();
    const bool soc = band.soc();
    const int nspin = band.nspins();
    const int nomega = *(int*)in.parameter("n_omega");

    std::string dir_impurity_solving = "impurity_solving";

    std::stringstream step_dir_ss;
    step_dir_ss << "/step" << istep-1;
    std::string step_dir= step_dir_ss.str();

    for(int ineq=0; ineq<ineq_num; ineq++)
    {
      const int iatom = atom.ineq_iatom(ineq);
      const int m_tot=norb_sub[iatom];

      const std::vector<std::vector<std::vector<std::complex<double>>>>& 
            sigma_a = sigma[ineq];

      std::stringstream site_dir_ss;
      site_dir_ss << "/impurity" << ineq;
      std::string site_dir= site_dir_ss.str();

      std::stringstream current_dir_ss;
      current_dir_ss << dir_impurity_solving << step_dir << site_dir;
      std::string current_dir = current_dir_ss.str();

      std::string sigma_file = current_dir+"/sigma_omega.dat";
      std::ofstream ofs_sig(sigma_file.c_str(), std::ios::out);

      for(int iomega=0; iomega<nomega; iomega++)
      {
        for(int is=0; is<2; is++)
        {
          for(int m=0; m<m_tot; m++)
          {
            if(nspin==1)
              ofs_sig << std::left << std::setw(5) << iomega+1 
              << std::setw(4) << m+1 << std::setw(4) << is+1
              << std::setw(20) << std::fixed << std::setprecision(12) 
              << sigma_a[iomega][0][m*m_tot+m].real()*Hartree_to_eV << " "
              << std::setw(20) << std::fixed << std::setprecision(12) 
              << sigma_a[iomega][0][m*m_tot+m].imag()*Hartree_to_eV << '\n';
            else
              ofs_sig << std::left << std::setw(5) << iomega+1
              << std::setw(4) << m+1 << std::setw(4) << is+1
              << std::setw(20) << std::fixed << std::setprecision(12) 
              << sigma_a[iomega][is][m*m_tot+m].real()/Hartree_to_eV << " "
              << std::setw(20) << std::fixed << std::setprecision(12) 
              << sigma_a[iomega][is][m*m_tot+m].imag()/Hartree_to_eV << '\n';       
          }//m
        }//is
      }//iomega
      ofs_sig.close();

    }//ineq

    return;
  }
  
}
