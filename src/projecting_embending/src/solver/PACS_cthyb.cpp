#include "PACS_cthyb.h"

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
  void PACS_CTHYB::output(
        const int char_step, const int DMFT_step,
        const double mu, DMFT::input_info& in, DFT_output::atoms_info& atom, 
        DFT_output::KS_bands& band, const std::vector<double>& freq,
        std::vector<std::vector<std::vector<std::complex<double>>>>& Eimp,
        std::vector<std::vector<std::vector<std::vector<std::complex<double>>>>>& Gf_in,
        std::vector<std::vector<std::vector<std::vector<std::complex<double>>>>>& Weiss,
        std::vector<std::vector<std::vector<std::vector<std::complex<double>>>>>& hyb_omega)
  {
    //======================================================
    //             energy unit: eV
    //======================================================
    debug::codestamp("PACS_CTHYB::output");

    const int ineq_num = atom.inequ_atoms();
    const std::vector<int>& norb_sub = atom.iatom_norb();
    const bool soc = band.soc();
    const int nspin = band.nspins();
    const int ntau = *(int*)in.parameter("n_tau");
    const int nomega = *(int*)in.parameter("n_omega");
    double zero=0.0;

        //Create directory impurity 
    std::string dir_impurity_solving = "dmft";
    std::stringstream make_dir1;
    make_dir1 << "test -d " << dir_impurity_solving << " || mkdir " << dir_impurity_solving;
    system(make_dir1.str().c_str());

    std::stringstream char_dir_ss;
    char_dir_ss << "/charge_step" << char_step;
    std::string char_step_dir= char_dir_ss.str();

    std::stringstream make_char_dir;
    make_char_dir << "test -d " << dir_impurity_solving << char_step_dir
            << " || mkdir " << dir_impurity_solving << char_step_dir;
    system(make_char_dir.str().c_str());

    std::stringstream step_dir_ss;
    step_dir_ss << "/dmft_step" << DMFT_step;
    std::string step_dir= step_dir_ss.str();

    std::stringstream make_dmft_dir;
    make_dmft_dir << "test -d " << dir_impurity_solving << char_step_dir << step_dir
            << " || mkdir " << dir_impurity_solving << char_step_dir << step_dir;
    system(make_dmft_dir.str().c_str());

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
      std::stringstream make_imp_dir;
      make_imp_dir << "test -d " << dir_impurity_solving << char_step_dir << step_dir << site_dir
              << " || mkdir " << dir_impurity_solving << char_step_dir << step_dir << site_dir;
      system(make_imp_dir.str().c_str());

      std::stringstream current_dir_ss;
      current_dir_ss << dir_impurity_solving << char_step_dir << step_dir << site_dir;
      std::string current_dir = current_dir_ss.str();

      //===========================================
      //           write  control.in
      //===========================================
      std::string input_file = current_dir+"/parameters.in";
      std::ofstream ofs_input(input_file.c_str(), std::ios::out);

      ofs_input << "&Model_Parameter\n";
      ofs_input << "      nOrbit = " << m_tot << '\n';
      ofs_input << "        Beta = " << std::fixed << std::setprecision(9) 
                << *(double*)in.parameter("beta")/Hartree_to_eV  << '\n';
      ofs_input << "          xU = " << atom.Uval(iatom)*Hartree_to_eV  << '\n';
      ofs_input << "          xJ = " << atom.Jval(iatom)*Hartree_to_eV  << '\n';
      ofs_input << "      nOmega = " << nomega << '\n';
      ofs_input << "        nTau = " << ntau << '\n';
      ofs_input << "/\n\n";

      ofs_input << "&MC_Parameter\n";
      ofs_input << "   nWarm = 100000\n";
      ofs_input << "   nMeasure = " << *(long long*)in.parameter("mc_step")/10 << '\n';
      ofs_input << "    nBin = 10\n";
      ofs_input << "/\n";

      ofs_input.close();

      //===========================================
      //           write  delta.dat
      //===========================================
      std::string delta_file = current_dir+"/delta.dat";
      std::ofstream ofs_delta(delta_file.c_str(), std::ios::out);

      for(int iomega=0; iomega<nomega; iomega++)
      {
        ofs_delta << std::setw(22) << std::fixed << std::setprecision(15)
          << Hartree_to_eV*freq[iomega];
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
      // //===========================================
      // std::string deltat_file = current_dir+"/delta_t.dat";
      // std::ofstream ofs_deltat(deltat_file.c_str(), std::ios::out);

      // for(int itau=0; itau<ntau+1; itau++)
      // {
      //   for(int is=0; is<2; is++)
      //   {
      //     for(int m=0; m<m_tot; m++)
      //     {
      //       if(nspin==1)
      //         ofs_deltat << std::left << std::setw(5) << itau
      //         << std::setw(4) << m+1 << std::setw(4) << is+1
      //         << std::setw(22) << std::fixed << std::setprecision(15) 
      //         << std::pow(Hartree_to_eV,2)*hybt[0][itau][m*m_tot+m].real() << '\n';
      //       else
      //         ofs_deltat << std::left << std::setw(5) << itau
      //         << std::setw(4) << m+1 << std::setw(4) << is+1
      //         << std::setw(22) << std::fixed << std::setprecision(15) 
      //         << std::pow(Hartree_to_eV,2)*hybt[is][itau][m*m_tot+m].real() << '\n';      
      //     }//m
      //   }//is
      // }//iomega
      // ofs_deltat.close();

      //===========================================
      //           write  mu_vector.dat
      //===========================================
      std::string hopping_file = current_dir+"/mu_vector.dat";
      std::ofstream ofs_hopping(hopping_file.c_str(), std::ios::out);

      int count1=0;
      for(int is=0; is<2; is++)
      {
        for(int m=0; m<m_tot; m++)
        {
          if(nspin==1)
            ofs_hopping << std::setw(22) << std::fixed << std::setprecision(15) 
              << Hartree_to_eV*(mu - hoppinga[0][m*m_tot+m].real());
          else
            ofs_hopping << std::setw(22) << std::fixed << std::setprecision(15) 
              << Hartree_to_eV*(mu - hoppinga[is][m*m_tot+m].real());
        }//is
      }//m
      ofs_hopping.close();

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
        ofs_gf << std::setw(5) << iomega;
        for(int is=0; is<2; is++)
        {
          for(int m=0; m<m_tot; m++)
          {
            if(nspin==1)
              ofs_gf << std::setw(22) << std::fixed << std::setprecision(15) 
                     << Gf_ina[0][iomega][m*m_tot+m].real()/Hartree_to_eV
                     << std::setw(22) << std::fixed << std::setprecision(15) 
                     << Gf_ina[0][iomega][m*m_tot+m].imag()/Hartree_to_eV;
            else
              ofs_gf << std::setw(22) << std::fixed << std::setprecision(15) 
                     << Gf_ina[is][iomega][m*m_tot+m].real()/Hartree_to_eV
                     << std::setw(22) << std::fixed << std::setprecision(15) 
                     << Gf_ina[is][iomega][m*m_tot+m].imag()/Hartree_to_eV ;      
          }//m
        }//is
        ofs_gf << std::endl;
      }//iomega
      ofs_gf.close();

    }//ineq

    return;
  }

  void PACS_CTHYB::read_last_step(
        const int char_step,
        const int DMFT_step,
        DFT_output::KS_bands& band,
        DMFT::input_info& in, 
        DFT_output::atoms_info& atom,
        std::vector<std::vector<std::
        vector<std::vector<
        std::complex<double>>>>>& Gw_qmc,
        std::vector<std::vector<std::
        vector<std::vector<
        std::complex<double>>>>>& Gw_save,
        std::vector<std::vector<std::
        vector<std::vector<
        std::complex<double>>>>>& Sw )
  {
    debug::codestamp("PACS_CTHYB::::read_last_step");

    const int ineq_num = atom.inequ_atoms();
    const std::vector<int>& norb_sub = atom.iatom_norb();
    const bool soc = band.soc();
    const int nspin = band.nspins();
    const int magnetism = *(int*)in.parameter("magnetism");

    double omega, real, imag;
    std::string str_tmp;

    std::string dir_dmft_solving = "dmft";

    std::stringstream char_dir_ss;
    char_dir_ss << "/charge_step" << char_step;
    std::string char_step_dir= char_dir_ss.str();

    std::stringstream dmft_dir_ss;
    dmft_dir_ss << "/dmft_step" << DMFT_step-1;
    std::string dmft_step_dir= dmft_dir_ss.str();

    for(int ineq=0; ineq<ineq_num; ineq++)
    {
      const int iatom = atom.ineq_iatom(ineq);
      const int m_tot=norb_sub[iatom];
      const int nomega = Gw_save[ineq][0].size();

      std::stringstream site_dir_ss;
      site_dir_ss << "/impurity" << ineq;
      std::string site_dir= site_dir_ss.str();

      std::vector<std::vector<std::vector<std::complex<double>>>>&
            Gw_savea = Gw_save[ineq];

      std::vector<std::vector<std::vector<std::complex<double>>>>&
            Gw_qmca = Gw_qmc[ineq];
      
      std::vector<std::vector<std::vector<std::complex<double>>>>&
            Swa = Sw[ineq];

      std::stringstream current_dir_ss;
      current_dir_ss << dir_dmft_solving << char_step_dir << dmft_step_dir << site_dir;
      std::string current_dir = current_dir_ss.str();

      //=================================================
      //    Read interacting Green function of last step
      //=================================================
      std::string Gf_file = current_dir+"/Gw.dat";
      std::ifstream ifs_gf(Gf_file.c_str(), std::ios::in);

      if (!ifs_gf)  
	    {
	    	std::cout << "Fail to oepn " << Gf_file.c_str() << std::endl;
        std::exit(EXIT_FAILURE);
      }

      std::vector<std::vector<double>> Gw_real(2);
      std::vector<std::vector<double>> Gw_im(2);
      for(int is=0; is<2; is++)
      {
        Gw_real[is].resize(m_tot);
        Gw_im[is].resize(m_tot);
      }

      ifs_gf.seekg(0);    //set the position at the beginning of the file
      ifs_gf >> str_tmp;
      ifs_gf.ignore(400,'\n');
      int count=0;
      while(ifs_gf.good())
      {     
        ifs_gf >> omega;
        if(ifs_gf.eof()) break; //Check whether end of file is reached
       
        for(int is=0; is<2; is++)
        {
          for(int m=0; m<m_tot; m++)
          {           
            ifs_gf >> Gw_real[is][m];
            ifs_gf >> Gw_im[is][m];
          }
        }
        ifs_gf.ignore(150,'\n');

        for(int is=0; is<nspin; is++)
        {
          for(int m=0; m<m_tot; m++)
          {
            if(magnetism==3 || magnetism==4)//none magnetic or paramamagnetic
            {
              if(nspin==1) 
              {
                Gw_real[0][m] = (Gw_real[0][m]+Gw_real[1][m])/2.0;
                Gw_im[0][m] = (Gw_im[0][m]+Gw_im[1][m])/2.0;
              }
              else
              {
                Gw_real[0][m] = (Gw_real[0][m]+Gw_real[1][m])/2.0;
                Gw_real[1][m] = Gw_real[0][m];
                Gw_im[0][m] = (Gw_im[0][m]+Gw_im[1][m])/2.0;
                Gw_im[1][m] = Gw_im[0][m];
              }
            }
          }       
          DFT_output::atoms_info::symmetry_operation_vector<double>(
                    atom.local_sym(), atom.L(ineq), 
                    m_tot, &Gw_real[is][0]);

          DFT_output::atoms_info::symmetry_operation_vector<double>(
                    atom.local_sym(), atom.L(ineq), 
                    m_tot, &Gw_im[is][0]);
                 
        }

        for(int is=0; is<nspin; is++)
          for(int m=0; m<m_tot; m++)
            Gw_qmca[is][count][m*m_tot+m] = Hartree_to_eV*
            std::complex<double>(Gw_real[is][m],Gw_im[is][m]);

        count++;
        if(ifs_gf.eof()) break;//Check whether end of file is reached       
      }
      ifs_gf.close();

      if(count<nomega)
      {
        std::cout << "The number of Matsubara points of Gw.dat is less than nomega\n";
        std::exit(EXIT_FAILURE);
      }

      //=============================================
      //    Read input Green function of last step
      //=============================================
      std::string Gf_save_file = current_dir+"/Gf.in";
      std::ifstream ifs_gfsave(Gf_save_file.c_str(), std::ios::in);

      if (!ifs_gfsave)  
	    {
	    	std::cout << "Fail to oepn " << Gf_save_file.c_str() << std::endl;
        std::exit(EXIT_FAILURE);
      }

      ifs_gfsave.seekg(0);    //set the position at the beginning of the file
      count=0;
      while(ifs_gfsave.good())
      {
        ifs_gfsave >> omega;
        if(ifs_gfsave.eof()) break; //Check whether end of file is reached 
        for(int is=0; is<2; is++)
        {
          for(int m=0; m<m_tot; m++)
          {
            ifs_gfsave >> Gw_real[is][m];
            ifs_gfsave >> Gw_im[is][m];
          }
        }
        ifs_gfsave.ignore(150,'\n');

        for(int is=0; is<nspin; is++)
          for(int m=0; m<m_tot; m++)
            Gw_savea[is][count][m*m_tot+m] = Hartree_to_eV*
            std::complex<double>(Gw_real[is][m],Gw_im[is][m]);

        count++;
        if(ifs_gfsave.eof()) break;//Check whether end of file is reached       
      }
      ifs_gfsave.close();

      //=================================================
      //    Read self-energy of last step
      //=================================================
      std::string Gw_file = current_dir+"/Sigma.dat";
      std::ifstream ifSw(Gw_file.c_str(), std::ios::in);

      if (!ifSw)  
	    {
	    	std::cout << "Fail to oepn " << Gw_file.c_str() << std::endl;
        std::exit(EXIT_FAILURE);
      }

      ifSw.seekg(0);    //set the position at the beginning of the file
      ifSw >> str_tmp;
      ifSw.ignore(400,'\n');
      count=0;
      while(ifSw.good())
      {
        ifSw >> omega;
        if(ifSw.eof()) break; //Check whether end of file is reached 
        for(int is=0; is<2; is++)
        {
          for(int m=0; m<m_tot; m++)
          {
            ifSw >> Gw_real[is][m];
            ifSw >> Gw_im[is][m];
          }
        }
        ifSw.ignore(150,'\n');

        for(int is=0; is<nspin; is++)
        {
          for(int m=0; m<m_tot; m++)
          {
            if(magnetism==3 || magnetism==4) //none magnetic or paramamagnetic
            {
              if(nspin==1) 
              {
                Gw_real[0][m] = (Gw_real[0][m]+Gw_real[1][m])/2.0;
                Gw_im[0][m] = (Gw_im[0][m]+Gw_im[1][m])/2.0;
              }
              else
              {
                Gw_real[0][m] = (Gw_real[0][m]+Gw_real[1][m])/2.0;
                Gw_real[1][m] = Gw_real[0][m];
                Gw_im[0][m] = (Gw_im[0][m]+Gw_im[1][m])/2.0;
                Gw_im[1][m] = Gw_im[0][m];
              }
            }
          }
          DFT_output::atoms_info::symmetry_operation_vector<double>(
              atom.local_sym(), atom.L(ineq), m_tot, &Gw_real[is][0] );

          DFT_output::atoms_info::symmetry_operation_vector<double>(
              atom.local_sym(), atom.L(ineq), m_tot, &Gw_im[is][0] );
        }

        for(int is=0; is<nspin; is++)
          for(int m=0; m<m_tot; m++)
            Swa[is][count][m*m_tot+m] = 
            std::complex<double>(Gw_real[is][m],Gw_im[is][m])/Hartree_to_eV;

        count++;
        if(ifSw.eof()) break;  //Check whether end of file is reached       
      }
      ifSw.close();

      if(count<nomega)
      {
        std::cout << "The number of Matsubara points of Sigma.dat is less than nomega" << std::endl;
        std::exit(EXIT_FAILURE);
      }

    }//ineq

    return;
  }

}
