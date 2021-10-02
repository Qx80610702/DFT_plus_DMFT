#include "alps_cthyb_segment.h"

#include "../debug/debug.h"
#include "../timer.h"
#include "math_zone.h"
#include "../constants.h"

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
  void ALPS_CTHYB_SEGMENT::output(const int istep, 
        const double mu, DMFT::input_info& in, 
        DFT_output::atoms_info& atom, DFT_output::KS_bands& band,
        std::vector<std::vector<std::vector<std::complex<double>>>>& Eimp,
        std::vector<std::vector<std::vector<std::vector<std::complex<double>>>>>& Gf_in,
        std::vector<std::vector<std::vector<std::vector<std::complex<double>>>>>& Weiss)
  {
    //======================================================
    //        energy unit: eV
    //======================================================
    debug::codestamp("ALPS_CTHYB_SEGMENT::output");

    const int ineq_num = atom.inequ_atoms();
    const std::vector<int>& norb_sub = atom.iatom_norb();
    const bool soc = band.soc();
    const int nspin = band.nspins();
    const int ntau = *(int*)in.parameter("n_tau");
    const int nomega = *(int*)in.parameter("n_omega");

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
      const int m_tot=norb_sub[iatom];

      std::vector<std::vector<std::vector<std::complex<double>>>>& 
            hyb_taua = this->hyb_tau[ineq];
      
      std::vector<std::vector<std::complex<double>>>&
            hoppinga = Eimp[ineq];
      
      std::vector<std::vector<std::vector<std::complex<double>>>>&
            Gf_ina = Gf_in[ineq];

      std::vector<std::vector<std::vector<std::complex<double>>>>& 
            Weissa = Weiss[ineq];

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
        //           write  hyb.param
        //===========================================
        std::string param_file = current_dir+"/hyb.param";
        std::ofstream ofs_param(param_file.c_str(), std::ios::out);

        long long sweeps = *(long long*)in.parameter("mc_step")/20;

        ofs_param << "FLAVORS=" << 2*m_tot << '\n';
        ofs_param << "N_TAU=" << *(int*)in.parameter("n_tau") << '\n';        
        ofs_param << "BETA=" << std::fixed << std::setprecision(9) 
        << *(double*)in.parameter("beta")/Hartree_to_eV << '\n';  //Hartree to eV
        ofs_param << "THERMALIZATION=100000\n";
        ofs_param << "N_MEAS=20\n";
        ofs_param << "SWEEPS=" << sweeps << '\n';
        ofs_param << "DELTA=delta.dat\n";
        ofs_param << "MU_VECTOR=mu_vector.dat\n";
        ofs_param << "U_MATRIX=Uij.dat\n";
        ofs_param << "N_LEGENDRE=50\n";
        ofs_param << "N_HISTOGRAM_ORDERS=200\n";
        ofs_param << "NMATSUBARA=" << *(int*)in.parameter("n_omega") << '\n';
        
        ofs_param.close();

      //===========================================
      //           write  delta.txt
      //===========================================
      std::string delta_file = current_dir+"/delta.dat";
      std::ofstream ofs_delta(delta_file.c_str(), std::ios::out);

      for(int itau=0; itau<ntau+1; itau++)
      {
        ofs_delta << std::setw(5) << itau;

        for(int m=0; m<m_tot; m++)
        {
          if(nspin==2)
          {
            ofs_delta << std::setw(22) << std::fixed << std::setprecision(15) << 
            std::pow(Hartree_to_eV,2)*hyb_taua[0][itau][m*m_tot+m].real()
            << std::setw(22) << std::fixed << std::setprecision(15)
            << std::pow(Hartree_to_eV,2)*hyb_taua[1][itau][m*m_tot+m].real();
          }
          else
          {
            ofs_delta << std::setw(22) << std::fixed << std::setprecision(15) << 
            std::pow(Hartree_to_eV,2)*hyb_taua[0][itau][m*m_tot+m].real()
            << std::setw(22) << std::fixed << std::setprecision(15)
            << std::pow(Hartree_to_eV,2)*hyb_taua[0][itau][m*m_tot+m].real();
          }
        }//m_index

        ofs_delta << '\n';
      }//itau
      ofs_delta.close();

      //===========================================
      //           write  mu_vector.dat
      //===========================================
      std::string hopping_file = current_dir+"/mu_vector.dat";
      std::ofstream ofs_hopping(hopping_file.c_str(), std::ios::out);

      for(int m=0; m<m_tot; m++)
      {
        if(nspin==2)
          ofs_hopping << std::left << std::fixed << std::setprecision(15) 
          << Hartree_to_eV*(mu-hoppinga[0][m*m_tot+m].real()) << ' '
          << std::fixed << std::setprecision(15) 
          << Hartree_to_eV*(mu-hoppinga[1][m*m_tot+m].real()) << ' ';
        else
          ofs_hopping << std::left << std::fixed << std::setprecision(15) 
          << Hartree_to_eV*(mu-hoppinga[0][m*m_tot+m].real()) << ' '
          << std::fixed << std::setprecision(15) 
          << Hartree_to_eV*(mu-hoppinga[0][m*m_tot+m].real()) << ' ';
      }
      ofs_hopping << '\n';
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
          for(int m=0; m<m_tot; m++)
          {
            if(nspin==1)
              ofs_G0 << std::setw(5) << iomega << std::setw(3) << is << std::setw(3) << m
                     << std::setw(22) << std::fixed << std::setprecision(15) 
                     << Weissa[0][iomega][m*m_tot+m].real()/Hartree_to_eV << " "
                     << std::setw(22) << std::fixed << std::setprecision(15) 
                     << Weissa[0][iomega][m*m_tot+m].imag()/Hartree_to_eV << '\n' ;
            else
              ofs_G0 << std::setw(5) << iomega << std::setw(3) << is << std::setw(3) << m
                     << std::setw(22) << std::fixed << std::setprecision(15) 
                     << Weissa[is][iomega][m*m_tot+m].real()/Hartree_to_eV << " "
                     << std::setw(22) << std::fixed << std::setprecision(15) 
                     << Weissa[is][iomega][m*m_tot+m].imag()/Hartree_to_eV << '\n' ;   
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
              ofs_gf << std::setw(5) << iomega << std::setw(3) << is << std::setw(3) << m
                     << std::setw(22) << std::fixed << std::setprecision(15) 
                     << Gf_ina[0][iomega][m*m_tot+m].real()/Hartree_to_eV
                     << std::setw(22) << std::fixed << std::setprecision(15) 
                     << Gf_ina[0][iomega][m*m_tot+m].imag()/Hartree_to_eV << '\n';
            else
              ofs_gf << std::setw(5) << iomega << std::setw(3) << is << std::setw(3) << m
                     << std::setw(22) << std::fixed << std::setprecision(15) 
                     << Gf_ina[is][iomega][m*m_tot+m].real()/Hartree_to_eV
                     << std::setw(22) << std::fixed << std::setprecision(15) 
                     << Gf_ina[is][iomega][m*m_tot+m].imag()/Hartree_to_eV << '\n';
          }//m_index
        }//is
      }//iomega
      ofs_gf.close();

    }//ineq

    return;
  }

  void ALPS_CTHYB_SEGMENT::read_last_step(
          const int istep, 
          DFT_output::KS_bands& band,
          DMFT::input_info& in, 
          DFT_output::atoms_info& atom,
          std::vector<std::vector<std::
          vector<std::vector<
          std::complex<double>>>>>& hyb_save,
          std::vector<std::vector<std::
          vector<std::vector<
          std::complex<double>>>>>& Gf_qmc,
          std::vector<std::vector<std::
          vector<std::vector<
          std::complex<double>>>>>& Gf_save)
  {
    debug::codestamp("ALPS_CTHYB_SEGMENT::read_last_step");

    const int ineq_num = atom.inequ_atoms();
    const bool soc = band.soc();
    const int nspin = band.nspins();

    int iomega, m1, m2;
    double real, imag;

    //=================================================
    //    PART 1: allocation
    //=================================================
    this->hopping_matrix.resize(ineq_num);
    for(int ineq=0; ineq<ineq_num; ineq++)
    {
      const int iatom = atom.ineq_iatom(ineq);
      const int angular_L = atom.L(iatom);
      const int m_tot = 2*angular_L+1;

      this->hopping_matrix[ineq].resize(nspin);
      for(int is=0; is<nspin; is++)
      {
        this->hopping_matrix[ineq][is].resize(m_tot*m_tot);
      }//is
    }//ineq

    //=================================================
    //    PART 2: read
    //=================================================
    //directory impurity 
    std::string dir_impurity_solving = "impurity_solving";

    //directory step+num
    std::stringstream step_dir_ss;
    step_dir_ss << "/step" << istep-1;
    std::string step_dir= step_dir_ss.str();

    for(int ineq=0; ineq<ineq_num; ineq++)
    {
      const int iatom = atom.ineq_iatom(ineq);
      const int angular_L = atom.L(iatom);
      const int m_tot = 2*angular_L+1;

      std::stringstream site_dir_ss;
      site_dir_ss << "/impurity" << ineq;
      std::string site_dir= site_dir_ss.str();
      
      std::vector<std::vector<std::complex<double>>>&
            hoppinga = this->hopping_matrix[ineq];

      std::vector<std::vector<std::vector<std::complex<double>>>>&
            Gf_savea = Gf_save[ineq];

      std::vector<std::vector<std::vector<std::complex<double>>>>&
            Gf_qmca = Gf_qmc[ineq];

      std::vector<std::vector<std::vector<std::complex<double>>>>&
            hyb_savea =hyb_save[ineq];
                  
      for(int is=0; is<nspin; is++)
      {
        std::stringstream spin_dir_ss;
        spin_dir_ss << "/spin" << is;
        std::string spin_dir= spin_dir_ss.str();

        std::stringstream current_dir_ss;
        current_dir_ss << dir_impurity_solving << step_dir << site_dir << spin_dir;
        std::string current_dir = current_dir_ss.str();

        //===========================================
        //    Read hopping matrix
        //===========================================
        std::vector<std::complex<double>>&
            hoppingb = hoppinga[is];

        std::string hopping_file = current_dir+"/hopping.txt";
        std::ifstream ifs_Tij(hopping_file.c_str(), std::ios::in);

        if (!ifs_Tij)  
	      {
	      	std::cout << "Fail to oepn " << hopping_file.c_str() << std::endl;
          std::exit(EXIT_FAILURE);
        }

        ifs_Tij.seekg(0);   //set the position at the beginning of the file
        while(ifs_Tij.good())
        {
          ifs_Tij >> m1;
          if(ifs_Tij.eof()) break; //Check whether end of file is reached   
          ifs_Tij >> m2;
          ifs_Tij >> real;
          ifs_Tij >> imag;
          ifs_Tij.ignore(150,'\n');

          hoppingb[m1*m_tot+m2] = std::complex<double>(real,imag);

          hoppingb[m1*m_tot+m2] /= Hartree_to_eV;   //eV to Hartree

          if(ifs_Tij.eof()) break;//Check whether end of file is reached       
        }
        ifs_Tij.close();

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
        std::vector<std::vector<std::complex<double>>>&
            hyb_saveb = hyb_savea[is];

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

          hyb_saveb[iomega][m1*m_tot+m2] = std::complex<double>(real,imag);

          if(ifs_hybsave.eof()) break; //Check whether end of file is reached       
        }
        ifs_hybsave.close();

      }//is
    }//ineq

    return;
  }
  
}
