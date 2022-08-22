#include "PACS_cthyb.h"

#include "../debug.h"
#include "../timer.h"
#include "math_zone.h"
#include "../constants.h"
#include "math_zone.h"
#include "../global_variables.h"
#include "../timer.h"
#include "../mpi_environment.h"
#include "../utilities.h"

#include <mpi.h>
#include <memory>
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>
#include <cstdlib>   //Use exit function
#include <unistd.h>

extern "C" void pacs_run_();

namespace DMFT
{
  void PACS_CTHYB::output(
        const int char_step, const int DMFT_step,
        const double mu, DMFT::input_info& in, DFT_output::atoms_info& atom, 
        DFT_output::KS_bands& band, const std::vector<double>& freq,
        std::vector<std::vector<std::vector<std::complex<double>>>>& Eimp,
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
    if(access(dir_impurity_solving.c_str(),0) != 0){
      if(mk_dir(dir_impurity_solving.c_str()) !=0){
        std::cout << "Fails to creat the directory " << dir_impurity_solving << std::endl;
        std::exit(EXIT_FAILURE);
      };
    }

    std::stringstream char_dir_ss;
    char_dir_ss << "/charge_step" << char_step;
    std::string char_step_dir= char_dir_ss.str();

    std::stringstream make_char_dir;
    make_char_dir << dir_impurity_solving << char_step_dir;
    if(access(make_char_dir.str().c_str(),0) != 0){
      if(mk_dir(make_char_dir.str().c_str()) !=0){
        std::cout << "Fails to creat the directory " << make_char_dir.str() << std::endl;
        std::exit(EXIT_FAILURE);
      };
    }

    std::stringstream step_dir_ss;
    step_dir_ss << "/dmft_step" << DMFT_step;
    std::string step_dir= step_dir_ss.str();

    std::stringstream make_dmft_dir;
    make_dmft_dir << dir_impurity_solving << char_step_dir << step_dir;
    if(access(make_dmft_dir.str().c_str(),0) != 0){
      if(mk_dir(make_dmft_dir.str().c_str()) !=0){
        std::cout << "Fails to creat the directory " << make_dmft_dir.str() << std::endl;
        std::exit(EXIT_FAILURE);
      };
    }
    
    for(int ineq=0; ineq<ineq_num; ineq++)
    {
      const int iatom = atom.ineq_iatom(ineq);
      const int m_tot = norb_sub[iatom];
      
      const std::vector<std::vector<std::complex<double>>>&
            hoppinga = Eimp[ineq];
      
      // const std::vector<std::vector<std::vector<std::complex<double>>>>&
      //       Sigma_ina = Sigma_in[ineq];

      const std::vector<std::vector<std::vector<std::complex<double>>>>& 
            hyb = hyb_omega[ineq];
      
      // const std::vector<std::vector<std::vector<std::complex<double>>>>& 
      //       G0 = Weiss[ineq];

      const std::vector<std::vector<std::vector<std::complex<double>>>>& 
            hybt = this->hyb_tau[ineq];

      std::stringstream site_dir_ss;
      site_dir_ss << "/impurity" << ineq;
      std::string site_dir= site_dir_ss.str();

      std::stringstream make_imp_dir;
      make_imp_dir << dir_impurity_solving << char_step_dir << step_dir << site_dir;
      if(access(make_imp_dir.str().c_str(),0) != 0){
        if(mk_dir(make_imp_dir.str().c_str()) !=0){
          std::cout << "Fails to creat the directory " << make_imp_dir.str() << std::endl;
          std::exit(EXIT_FAILURE);
        };
      }

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
                << *(double*)in.parameter("beta")/GLC::Hartree_to_eV  << '\n';
      ofs_input << "          xU = " << atom.Uval(iatom)*GLC::Hartree_to_eV  << '\n';
      ofs_input << "          xJ = " << atom.Jval(iatom)*GLC::Hartree_to_eV  << '\n';
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
          << GLC::Hartree_to_eV*freq[iomega];
        for(int is=0; is<2; is++)
        {
          for(int m=0; m<m_tot; m++)
          {
            if(nspin==1)
              ofs_delta << std::setw(22) << std::fixed << std::setprecision(15) 
              << GLC::Hartree_to_eV*hyb[0][iomega][m*m_tot+m].real() 
              << std::setw(22) << std::fixed << std::setprecision(15) 
              << GLC::Hartree_to_eV*hyb[0][iomega][m*m_tot+m].imag();
            else
              ofs_delta << std::setw(22) << std::fixed << std::setprecision(15) 
              << GLC::Hartree_to_eV*hyb[is][iomega][m*m_tot+m].real() 
              << std::setw(22) << std::fixed << std::setprecision(15) 
              << GLC::Hartree_to_eV*hyb[is][iomega][m*m_tot+m].imag();       
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
      //         << std::pow(GLC::Hartree_to_eV,2)*hybt[0][itau][m*m_tot+m].real() << '\n';
      //       else
      //         ofs_deltat << std::left << std::setw(5) << itau
      //         << std::setw(4) << m+1 << std::setw(4) << is+1
      //         << std::setw(22) << std::fixed << std::setprecision(15) 
      //         << std::pow(GLC::Hartree_to_eV,2)*hybt[is][itau][m*m_tot+m].real() << '\n';      
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
              << GLC::Hartree_to_eV*(mu - hoppinga[0][m*m_tot+m].real());
          else
            ofs_hopping << std::setw(22) << std::fixed << std::setprecision(15) 
              << GLC::Hartree_to_eV*(mu - hoppinga[is][m*m_tot+m].real());
        }//is
      }//m
      ofs_hopping.close();

      //==================================================
      //   write  Gf.in; the input self-energy
      //   of current step (in Matrsubara frequency), 
      //   which will be read by last step to judge whether
      //   the self-consistency is achieved
      //==================================================
      // std::string Sig_file = current_dir+"/Gf.in";
      // std::ofstream ofs_Sig(Sig_file.c_str(), std::ios::out);

      // for(int iomega=0; iomega<nomega; iomega++)
      // {
      //   ofs_Sig << std::setw(5) << iomega;
      //   for(int is=0; is<2; is++)
      //   {
      //     for(int m=0; m<m_tot; m++)
      //     {
      //       if(nspin==1)
      //         ofs_Sig << std::setw(22) << std::fixed << std::setprecision(15) 
      //                << Sigma_ina[0][iomega][m*m_tot+m].real()*GLC::Hartree_to_eV
      //                << std::setw(22) << std::fixed << std::setprecision(15) 
      //                << Sigma_ina[0][iomega][m*m_tot+m].imag()*GLC::Hartree_to_eV;
      //       else
      //         ofs_Sig << std::setw(22) << std::fixed << std::setprecision(15) 
      //                << Sigma_ina[is][iomega][m*m_tot+m].real()*GLC::Hartree_to_eV
      //                << std::setw(22) << std::fixed << std::setprecision(15) 
      //                << Sigma_ina[is][iomega][m*m_tot+m].imag()*GLC::Hartree_to_eV ;      
      //     }//m
      //   }//is
      //   ofs_Sig << std::endl;
      // }//iomega
      // ofs_Sig.close();

    }//ineq

    return;
  }

  void PACS_CTHYB::read_self_energy(
        const int char_step,
        const int DMFT_step,
        DFT_output::KS_bands& band,
        DMFT::input_info& in, 
        DFT_output::atoms_info& atom,
        std::vector<std::vector<std::
        vector<std::vector<
        std::complex<double>>>>>& Sw)
  {
    debug::codestamp("PACS_CTHYB::read_last_step");

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
    dmft_dir_ss << "/dmft_step" << DMFT_step;
    std::string dmft_step_dir= dmft_dir_ss.str();

    for(int ineq=0; ineq<ineq_num; ineq++)
    {
      const int iatom = atom.ineq_iatom(ineq);
      const int m_tot=norb_sub[iatom];
      const int nomega = Sw[ineq][0].size();

      int count = 0;

      std::stringstream site_dir_ss;
      site_dir_ss << "/impurity" << ineq;
      std::string site_dir= site_dir_ss.str();

      // std::vector<std::vector<std::vector<std::complex<double>>>>&
      //       Sw_savea = Sw_save[ineq];
      
      std::vector<std::vector<std::vector<std::complex<double>>>>&
            Swa = Sw[ineq];

      std::stringstream current_dir_ss;
      current_dir_ss << dir_dmft_solving << char_step_dir << dmft_step_dir << site_dir;
      std::string current_dir = current_dir_ss.str();

      //=================================================
      //    Read interacting Green function of last step
      //=================================================
      // std::string Gf_file = current_dir+"/Gw.dat";
      // std::ifstream ifs_gf(Gf_file.c_str(), std::ios::in);

      // if (!ifs_gf)  
	    // {
	    // 	std::cerr << "Fail to oepn " << Gf_file.c_str() << std::endl;
      //   std::exit(EXIT_FAILURE);
      // }

      // std::vector<std::vector<double>> Gw_real(2);
      // std::vector<std::vector<double>> Gw_im(2);
      // for(int is=0; is<2; is++)
      // {
      //   Gw_real[is].resize(m_tot);
      //   Gw_im[is].resize(m_tot);
      // }

      // ifs_gf.seekg(0);    //set the position at the beginning of the file
      // ifs_gf >> str_tmp;
      // ifs_gf.ignore(400,'\n');
      // int count=0;
      // while(ifs_gf.good())
      // {     
      //   ifs_gf >> omega;
      //   if(ifs_gf.eof()) break; //Check whether end of file is reached
       
      //   for(int is=0; is<2; is++)
      //   {
      //     for(int m=0; m<m_tot; m++)
      //     {           
      //       ifs_gf >> Gw_real[is][m];
      //       ifs_gf >> Gw_im[is][m];
      //     }
      //   }
      //   ifs_gf.ignore(150,'\n');

      //   for(int is=0; is<nspin; is++)
      //   {
      //     for(int m=0; m<m_tot; m++)
      //     {
      //       if(magnetism==3 || magnetism==4)//none magnetic or paramamagnetic
      //       {
      //         if(nspin==1) 
      //         {
      //           Gw_real[0][m] = (Gw_real[0][m]+Gw_real[1][m])/2.0;
      //           Gw_im[0][m] = (Gw_im[0][m]+Gw_im[1][m])/2.0;
      //         }
      //         else
      //         {
      //           Gw_real[0][m] = (Gw_real[0][m]+Gw_real[1][m])/2.0;
      //           Gw_real[1][m] = Gw_real[0][m];
      //           Gw_im[0][m] = (Gw_im[0][m]+Gw_im[1][m])/2.0;
      //           Gw_im[1][m] = Gw_im[0][m];
      //         }
      //       }
      //     }       
      //     DFT_output::atoms_info::symmetry_operation_vector<double>(
      //               atom.local_sym(), atom.L(ineq), 
      //               m_tot, &Gw_real[is][0]);

      //     DFT_output::atoms_info::symmetry_operation_vector<double>(
      //               atom.local_sym(), atom.L(ineq), 
      //               m_tot, &Gw_im[is][0]);
                 
      //   }

      //   for(int is=0; is<nspin; is++)
      //     for(int m=0; m<m_tot; m++)
      //       Gw_qmca[is][count][m*m_tot+m] = GLC::Hartree_to_eV*
      //       std::complex<double>(Gw_real[is][m],Gw_im[is][m]);

      //   count++;
      //   if(ifs_gf.eof()) break;//Check whether end of file is reached       
      // }
      // ifs_gf.close();

      // if(count<nomega)
      // {
      //   std::cerr << "The number of Matsubara points of Gw.dat is less than nomega\n";
      //   std::exit(EXIT_FAILURE);
      // }

      //=============================================
      //    Read input self-energy of last step
      //=============================================
      // std::vector<std::vector<double>> Gw_real(2);
      // std::vector<std::vector<double>> Gw_im(2);
      // for(int is=0; is<2; is++)
      // {
      //   Gw_real[is].resize(m_tot);
      //   Gw_im[is].resize(m_tot);
      // }

      // std::string Gf_save_file = current_dir+"/Sigma.in";
      // std::ifstream ifs_sigsave(Gf_save_file.c_str(), std::ios::in);

      // if (!ifs_sigsave)  
	    // {
	    // 	std::cerr << "Fail to oepn " << Gf_save_file.c_str() << std::endl;
      //   std::exit(EXIT_FAILURE);
      // }

      // ifs_sigsave.seekg(0);    //set the position at the beginning of the file
      // count=0;
      // while(ifs_sigsave.good())
      // {
      //   ifs_sigsave >> omega;
      //   if(ifs_sigsave.eof()) break; //Check whether end of file is reached 
      //   for(int is=0; is<2; is++)
      //   {
      //     for(int m=0; m<m_tot; m++)
      //     {
      //       ifs_sigsave >> Gw_real[is][m];
      //       ifs_sigsave >> Gw_im[is][m];
      //     }
      //   }
      //   ifs_sigsave.ignore(150,'\n');

      //   for(int is=0; is<nspin; is++)
      //     for(int m=0; m<m_tot; m++)
      //       Sw_savea[is][count][m*m_tot+m] = 
      //         std::complex<double>(Gw_real[is][m],Gw_im[is][m])/GLC::Hartree_to_eV;

      //   count++;
      //   if(ifs_sigsave.eof()) break;//Check whether end of file is reached       
      // }
      // ifs_sigsave.close();

      //=================================================
      //    Read self-energy of last step
      //=================================================
      std::vector<std::vector<double>> Gw_real(2);
      std::vector<std::vector<double>> Gw_im(2);
      for(int is=0; is<2; is++)
      {
        Gw_real[is].resize(m_tot);
        Gw_im[is].resize(m_tot);
      }

      std::string Gw_file = current_dir+"/Sigma.dat";
      std::ifstream ifSw(Gw_file.c_str(), std::ios::in);

      if (!ifSw)  
	    {
	    	std::cerr << "Fail to oepn " << Gw_file.c_str() << std::endl;
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
            std::complex<double>(Gw_real[is][m],Gw_im[is][m])/GLC::Hartree_to_eV;

        count++;
        if(ifSw.eof()) break;  //Check whether end of file is reached       
      }
      ifSw.close();

      if(count<nomega)
      {
        std::cerr << "The number of Matsubara points of Sigma.dat is less than nomega" << std::endl;
        std::exit(EXIT_FAILURE);
      }

    }//ineq

    return;
  }

  void PACS_CTHYB::impurities_solving(
          const int char_step,
          const int DMFT_step,
          DFT_output::atoms_info& atom )
  {
    debug::codestamp("PACS_CTHYB::impurities_solving");

    const int ineq_num = atom.inequ_atoms();

    std::string str_tmp;

    std::string dir_dmft_solving = "dmft";

    std::stringstream char_dir_ss;
    char_dir_ss << "/charge_step" << char_step;
    std::string char_step_dir= char_dir_ss.str();

    std::stringstream dmft_dir_ss;
    dmft_dir_ss << "/dmft_step" << DMFT_step;
    std::string dmft_step_dir= dmft_dir_ss.str();

    for(int ineq=0; ineq<ineq_num; ineq++)
    {
      const int iatom = atom.ineq_iatom(ineq);

      std::stringstream site_dir_ss;
      site_dir_ss << "/impurity" << ineq;
      std::string site_dir= site_dir_ss.str();
      
      std::stringstream current_dir_ss;
      current_dir_ss << dir_dmft_solving << char_step_dir << dmft_step_dir << site_dir;
      std::string current_dir = current_dir_ss.str();

      int ierr = chdir(current_dir.c_str());
      if(ierr != 0){
        std::cout << "Process " << mpi_rank() << " fails to enter to the directory " << current_dir << std::endl;
      }

      MPI_Barrier(MPI_COMM_WORLD);  //Blocks until all processes reach here

      std::string date;
      timer::get_date_time(date);
      GLV::ofs_running << "  impurity" << ineq << "     " << date;
      GLV::ofs_running.flush();
      
      double time, seconds;
      int hours, minutes;
      timer::timestamp(time);

      pacs_run_();

      timer::get_date_time(date);
      GLV::ofs_running << "       " << date;

      timer::get_time(time, seconds, minutes, hours);

      GLV::ofs_running << "             " << hours << "h " 
                           << minutes << "m "
                           << (int)seconds << "s" << std::endl;

      ierr = chdir("../../../../");
      if(ierr != 0){
        std::cout << "Process " << mpi_rank() << " fails to return to root directory "<< std::endl;
      }
      
      MPI_Barrier(MPI_COMM_WORLD);  //Blocks until all processes reach here
    }//ineq

    return;
  }

}
