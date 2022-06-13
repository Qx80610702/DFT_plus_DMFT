#include "alps_cthyb.h"

#include "../debug.h"
#include "../timer.h"
#include "math_zone.h"
#include "../constants.h"
#include "../global_variables.h"
#include "../utilities.h"

#include <omp.h>
#include <memory>
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>
#include <cstdlib>   //Use exit function
#include <unistd.h>

namespace DMFT
{
  void ALPS_CTHYB::output(
        const int char_step, const int DMFT_step, 
        const double mu, DMFT::input_info& in, 
        DFT_output::atoms_info& atom, DFT_output::KS_bands& band,
        std::vector<std::vector<std::vector<std::complex<double>>>>& Eimp,
        std::vector<std::vector<std::vector<std::vector<std::complex<double>>>>>& Sigma_in,
        std::vector<std::vector<std::vector<std::vector<std::complex<double>>>>>& Weiss,
        std::vector<std::vector<std::vector<std::vector<std::complex<double>>>>>& hyb_omega)
  {
    //======================================================
    //        energy unit: eV
    //======================================================
    debug::codestamp("ALPS_CTHYB::output");

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
      const int m_tot=norb_sub[iatom];

      const std::vector<std::vector<std::vector<std::complex<double>>>>& 
            hyb_taua = this->hyb_tau[ineq];
      
      const std::vector<std::vector<std::complex<double>>>&
            hoppinga = Eimp[ineq];
      
      const std::vector<std::vector<std::vector<std::complex<double>>>>&
            Sigma_ina = Sigma_in[ineq];

      const std::vector<std::vector<std::vector<std::complex<double>>>>& 
            Weissa = Weiss[ineq];
      
      const std::vector<std::vector<std::vector<std::complex<double>>>>& 
            hyb_omegaa = hyb_omega[ineq];

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

      int flavor;
      if(soc) flavor=2*m_tot;
      else flavor=2*m_tot;

      int timelimit;
      if(flavor<=6) timelimit=3600; //d system, non-soc
      else if(flavor>6 && flavor<=10) timelimit=7200; // f system, non-soc
      // else if(flavor>10 && flavor<11) timelimit=28800; //d system, soc
      // else if(flavor>=11 && flavor<15) timelimit=28800; //f system, soc

      std::stringstream current_dir_ss;
      current_dir_ss << dir_impurity_solving << char_step_dir << step_dir << site_dir;
      std::string current_dir = current_dir_ss.str();

      //===========================================
      //           write  input.ini
      //===========================================
      std::string input_file = current_dir+"/input.ini";
      std::ofstream ofs_input(input_file.c_str(), std::ios::out);

      ofs_input << "seed=20\n";
      ofs_input << "algorithm=complex-matrix\n";
      ofs_input << "timelimit=" << timelimit << '\n';
      ofs_input << "model.sites=" << m_tot << '\n';
      ofs_input << "model.spins=2\n";
      ofs_input << "model.coulomb_tensor_input_file=\"Uijkl.txt\"\n";
      ofs_input << "model.hopping_matrix_input_file=\"hopping.txt\"\n";
      ofs_input << "model.delta_input_file=\"delta.txt\"\n";
      ofs_input << "model.beta=" << std::fixed << std::setprecision(9) 
      << *(double*)in.parameter("beta")/GLC::Hartree_to_eV << '\n';  //Hartree to eV
      ofs_input << "model.n_tau_hyb=" << *(int*)in.parameter("n_tau") << '\n';
      ofs_input << "measurement.G1.n_tau=" << *(int*)in.parameter("n_tau") << '\n';
      ofs_input << "measurement.G1.n_matsubara=" << *(int*)in.parameter("n_omega") << '\n';
      ofs_input.close();

      //===========================================
      //           write  delta.txt
      //===========================================
      std::string delta_file = current_dir+"/delta.txt";
      std::ofstream ofs_delta(delta_file.c_str(), std::ios::out);

      for(int itau=0; itau<ntau+1; itau++)
      {
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
                  ofs_delta << std::left << std::setw(5) << itau 
                  << std::setw(4) << count1 << std::setw(4) << count2
                  << std::setw(22) << std::fixed << std::setprecision(15) 
                  << zero << " "
                  << std::setw(22) << std::fixed << std::setprecision(15) 
                  << zero << '\n';
                else
                  if(nspin==1)
                    ofs_delta << std::left << std::setw(5) << itau 
                    << std::setw(4) << count1 << std::setw(4) << count2
                    << std::setw(22) << std::fixed << std::setprecision(15) 
                    << std::pow(GLC::Hartree_to_eV,2)*hyb_taua[0][itau][m1*m_tot+m2].real() << " "
                    << std::setw(20) << std::fixed << std::setprecision(15) 
                    << std::pow(GLC::Hartree_to_eV,2)*hyb_taua[0][itau][m1*m_tot+m2].imag() << '\n';
                  else
                    ofs_delta << std::left << std::setw(5) << itau 
                    << std::setw(4) << count1 << std::setw(4) << count2
                    << std::setw(22) << std::fixed << std::setprecision(15) 
                    << std::pow(GLC::Hartree_to_eV,2)*hyb_taua[is1][itau][m1*m_tot+m2].real() << " "
                    << std::setw(22) << std::fixed << std::setprecision(15) 
                    << std::pow(GLC::Hartree_to_eV,2)*hyb_taua[is1][itau][m1*m_tot+m2].imag() << '\n';

                count2++;
              }//is2
            }//m2
            count1++;
          }//is1
        }//m1
      }//itau
      ofs_delta.close();

      //===========================================
      //           write  hopping.txt
      //===========================================
      std::string hopping_file = current_dir+"/hopping.txt";
      std::ofstream ofs_hopping(hopping_file.c_str(), std::ios::out);

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
                ofs_hopping << std::setw(4) << count1 << std::setw(4) << count2
                << std::setw(22) << std::fixed << std::setprecision(15) 
                << zero << " "
                << std::setw(22) << std::fixed << std::setprecision(15) 
                << zero << '\n';
              else
                if(nspin==1)
                  if(count1==count2)
                    ofs_hopping << std::setw(4) << count1 << std::setw(4) << count2
                    << std::setw(22) << std::fixed << std::setprecision(15) 
                    << GLC::Hartree_to_eV*(hoppinga[0][m1*m_tot+m2].real() - mu) << " "
                    << std::setw(22) << std::fixed << std::setprecision(15) 
                    << GLC::Hartree_to_eV*hoppinga[0][m1*m_tot+m2].imag() << '\n';
                  else
                    ofs_hopping << std::setw(4) << count1 << std::setw(4) << count2
                    << std::setw(22) << std::fixed << std::setprecision(15) 
                    << zero << " "
                    << std::setw(22) << std::fixed << std::setprecision(15) 
                    << zero << '\n';
                else
                  if(count1==count2)
                    ofs_hopping << std::setw(4) << count1 << std::setw(4) << count2
                    << std::setw(22) << std::fixed << std::setprecision(15) 
                    << GLC::Hartree_to_eV*(hoppinga[is1][m1*m_tot+m2].real() - mu) << " "
                    << std::setw(22) << std::fixed << std::setprecision(15) 
                    << GLC::Hartree_to_eV*hoppinga[is1][m1*m_tot+m2].imag() << '\n';
                  else
                    ofs_hopping << std::setw(4) << count1 << std::setw(4) << count2
                    << std::setw(22) << std::fixed << std::setprecision(15) 
                    << zero << " "
                    << std::setw(22) << std::fixed << std::setprecision(15) 
                    << zero << '\n';

              count2++;
            }//is2
          }//m2
          count1++;
        }//is1
      }//m1
      ofs_hopping.close();

      //==================================================
      //   write  Gf.in; the input self-energy 
      //   of current step (in Matrsubara frequency), 
      //   which will be read by last step to judge whether
      //   the self-consistency is achieved
      //==================================================
      std::string Sig_file = current_dir+"/Sigma.in";
      std::ofstream ofs_Sig(Sig_file.c_str(), std::ios::out);
  
      for(int iomega=0; iomega<nomega; iomega++)
      {
        for(int is=0; is<2; is++)
        {
          for(int m_index=0; m_index<m_tot*m_tot; m_index++)
          {
            int m1 = m_index/m_tot;
            int m2 = m_index%m_tot;

            if(nspin==1)
              ofs_Sig << std::left << std::setw(5) << iomega << std::setw(2) << is 
                  << std::setw(3) << m1 << std::setw(3) << m2
                  << std::setw(22) << std::fixed << std::setprecision(15) 
                  << Sigma_ina[0][iomega][m_index].real()*GLC::Hartree_to_eV << " "
                  << std::setw(22) << std::fixed << std::setprecision(15) 
                  << Sigma_ina[0][iomega][m_index].imag()*GLC::Hartree_to_eV << '\n';
            else
              ofs_Sig << std::left << std::setw(5) << iomega << std::setw(2) << is 
                  << std::setw(3) << m1 << std::setw(3) << m2
                  << std::setw(22) << std::fixed << std::setprecision(15) 
                  << Sigma_ina[is][iomega][m_index].real()*GLC::Hartree_to_eV << " "
                  << std::setw(22) << std::fixed << std::setprecision(15) 
                  << Sigma_ina[is][iomega][m_index].imag()*GLC::Hartree_to_eV << '\n';
          }//m_index
        }//iomega
      }
      ofs_Sig.close();

      //=====================================================
      //        write  Weiss function
      //=====================================================
      std::string Weiss_file = current_dir+"/Weiss.dat";
      std::ofstream ofs_Weiss(Weiss_file.c_str(), std::ios::out);

      for(int iomega=0; iomega<nomega; iomega++)
      {
        for(int is=0; is<2; is++)
        {
          for(int m_index=0; m_index<m_tot*m_tot; m_index++)
          {
            int m1 = m_index/m_tot;
            int m2 = m_index%m_tot;

            if(nspin==1)
              ofs_Weiss << std::left << std::setw(5) << iomega << std::setw(2) << is 
                  << std::setw(3) << m1 << std::setw(3) << m2
                  << std::setw(22) << std::fixed << std::setprecision(15) 
                  << Weissa[0][iomega][m_index].real()/GLC::Hartree_to_eV << " "
                  << std::setw(20) << std::fixed << std::setprecision(15) 
                  << Weissa[0][iomega][m_index].imag()/GLC::Hartree_to_eV << " " << '\n';
            else
              ofs_Weiss << std::left << std::setw(5) << iomega << std::setw(2) << is 
                  << std::setw(3) << m1 << std::setw(3) << m2
                  << std::setw(22) << std::fixed << std::setprecision(15) 
                  << Weissa[is][iomega][m_index].real()/GLC::Hartree_to_eV << " "
                  << std::setw(22) << std::fixed << std::setprecision(15) 
                  << Weissa[is][iomega][m_index].imag()/GLC::Hartree_to_eV << " " << '\n';

          }//m_index
        }//iomega
      }//is
      ofs_Weiss.close();

      //==========================================================
      //   write  hybridization function of Matsubara frequency
      //==========================================================
      std::string hyb_file = current_dir+"/delta_omega.dat";
      std::ofstream ofs_hyb(hyb_file.c_str(), std::ios::out);
    
      for(int iomega=0; iomega<nomega; iomega++)
      {
        for(int is=0; is<2; is++)
        {
          for(int m_index=0; m_index<m_tot*m_tot; m_index++)
          {
            int m1 = m_index/m_tot;
            int m2 = m_index%m_tot;

            if(nspin==1)
              ofs_hyb << std::left << std::setw(5) << iomega << std::setw(2) << is 
                  << std::setw(3) << m1 << std::setw(3) << m2
                  << std::setw(22) << std::fixed << std::setprecision(15) 
                  << hyb_omegaa[0][iomega][m_index].real()*GLC::Hartree_to_eV << " "
                  << std::setw(22) << std::fixed << std::setprecision(15) 
                  << hyb_omegaa[0][iomega][m_index].imag()*GLC::Hartree_to_eV << '\n';
            else
              ofs_hyb << std::left << std::setw(5) << iomega << std::setw(2) << is 
                  << std::setw(3) << m1 << std::setw(3) << m2
                  << std::setw(22) << std::fixed << std::setprecision(15) 
                  << hyb_omegaa[is][iomega][m_index].real()*GLC::Hartree_to_eV << " "
                  << std::setw(22) << std::fixed << std::setprecision(15) 
                  << hyb_omegaa[is][iomega][m_index].imag()*GLC::Hartree_to_eV << '\n';
          }//m_index
        }//iomega
      }//is
      ofs_hyb.close();

    }//ineq

    return;
  }

  void ALPS_CTHYB::read_last_step(
          const int char_step,
          const int DMFT_step, 
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

    int iomega, flavor1, flavor2;
    int ispin1, ispin2, m1, m2;
    double real, imag;
      
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

      std::stringstream site_dir_ss;
      site_dir_ss << "/impurity" << ineq;
      std::string site_dir= site_dir_ss.str();

      std::vector<std::vector<std::vector<std::complex<double>>>>&
            Gf_savea = Gf_save[ineq];

      std::vector<std::vector<std::vector<std::complex<double>>>>&
            Gf_qmca = Gf_qmc[ineq];

      std::vector<std::vector<std::vector<std::complex<double>>>>&
            Weissa =Weiss[ineq];

      std::stringstream current_dir_ss;
      current_dir_ss << dir_dmft_solving << char_step_dir << dmft_step_dir << site_dir;
      std::string current_dir = current_dir_ss.str();

      //=================================================
      //    Read interacting Green function of last step
      //=================================================
      std::string Gf_file = current_dir+"/G_omega.dat";
      std::ifstream ifs_gf(Gf_file.c_str(), std::ios::in);

      if (!ifs_gf)  
	    {
	    	GLV::ofs_error << "Fail to oepn " << Gf_file.c_str() << std::endl;
        std::exit(EXIT_FAILURE);
      }

      ifs_gf.seekg(0);    //set the position at the beginning of the file
      while(ifs_gf.good())
      {
        ifs_gf >> iomega;
        if(ifs_gf.eof()) break; //Check whether end of file is reached 
        ifs_gf >> flavor1;
        ifs_gf >> flavor2;
        ifs_gf >> real;
        ifs_gf >> imag;
        ifs_gf.ignore(150,'\n');

        ispin1 = flavor1%2;
        ispin2 = flavor2%2;
        m1 = flavor1/2;
        m2 = flavor2/2;

        if(ispin1==ispin2)
          if(nspin==2)
            Gf_qmca[ispin1][iomega][m1*m_tot+m2] = GLC::Hartree_to_eV*std::complex<double>(real,imag);
          else if(nspin==1 && ispin1==0)
            Gf_qmca[ispin1][iomega][m1*m_tot+m2] = GLC::Hartree_to_eV*std::complex<double>(real,imag);

        if(ifs_gf.eof()) break;//Check whether end of file is reached       
      }
      ifs_gf.close();

      //=============================================
      //    Read input Green function of last step
      //=============================================
      std::string Gf_save_file = current_dir+"/Gf.in";
      std::ifstream ifs_gfsave(Gf_save_file.c_str(), std::ios::in);

      if (!ifs_gfsave)  
	    {
	    	GLV::ofs_error << "Fail to oepn " << Gf_save_file.c_str() << std::endl;
        std::exit(EXIT_FAILURE);
      }

      ifs_gfsave.seekg(0);    //set the position at the beginning of the file
      while(ifs_gfsave.good())
      {
        ifs_gfsave >> iomega;
        if(ifs_gfsave.eof()) break; //Check whether end of file is reached 
        ifs_gfsave >> ispin1;
        ifs_gfsave >> m1;
        ifs_gfsave >> m2;
        ifs_gfsave >> real;
        ifs_gfsave >> imag;

        ifs_gfsave.ignore(150,'\n');

        if(nspin==2)
          Gf_savea[ispin1][iomega][m1*m_tot+m2] = GLC::Hartree_to_eV*std::complex<double>(real,imag);
        else if(nspin==1 && ispin1==0)
          Gf_savea[ispin1][iomega][m1*m_tot+m2] = GLC::Hartree_to_eV*std::complex<double>(real,imag);

        if(ifs_gfsave.eof()) break; //Check whether end of file is reached
      }
      ifs_gfsave.close();

      //===================================================
      //    Read Weiss Green's function of last step
      //===================================================
      std::string Weiss_file = current_dir+"/Weiss.dat";
      std::ifstream ifs_Weiss(Weiss_file.c_str(), std::ios::in);

      if (!ifs_Weiss)  
	    {
	    	GLV::ofs_error << "Fail to oepn " << Weiss_file.c_str() << std::endl;
        std::exit(EXIT_FAILURE);
      }

      ifs_Weiss.seekg(0);    //set the position at the beginning of the file
      while(ifs_Weiss.good())
      {
        ifs_Weiss >> iomega;
        if(ifs_Weiss.eof()) break; //Check whether end of file is reached
        ifs_Weiss >> ispin1;
        ifs_Weiss >> m1;
        ifs_Weiss >> m2;
        ifs_Weiss >> real;
        ifs_Weiss >> imag;
        ifs_Weiss.ignore(150,'\n');

        if(m1==m2)
        {
          if(nspin==2)
            Weissa[ispin1][iomega][m1*m_tot+m2] = GLC::Hartree_to_eV*std::complex<double>(real,imag);
          else if(nspin==1 && ispin1==0)
            Weissa[ispin1][iomega][m1*m_tot+m2] = GLC::Hartree_to_eV*std::complex<double>(real,imag);
        }

        if(ifs_Weiss.eof()) break; //Check whether end of file is reached       
      }
      ifs_Weiss.close();

    }//ineq

    return;
  }
  
  void ALPS_CTHYB::out_sigma_last_step(
        const int istep, DFT_output::KS_bands& band,
        DMFT::input_info& in, DFT_output::atoms_info& atom,
        const  std::vector<std::vector<std::vector<
        std::vector<std::complex<double>>>>>& sigma)
  {
    debug::codestamp("ALPS_CTHYB::out_sigma_last_step");

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
          for(int m_index=0; m_index<m_tot*m_tot; m_index++)
          {
            int m1 = m_index/m_tot;
            int m2 = m_index%m_tot;

            if(nspin==1)
              ofs_sig << std::left << std::setw(5) << iomega << std::setw(2) << is 
                  << std::setw(3) << m1 << std::setw(3) << m2
                  << std::setw(20) << std::fixed << std::setprecision(12) 
                  << sigma_a[0][iomega][m_index].real()*GLC::Hartree_to_eV << " "
                  << std::setw(20) << std::fixed << std::setprecision(12) 
                  << sigma_a[0][iomega][m_index].imag()*GLC::Hartree_to_eV << '\n';
            else
              ofs_sig << std::left << std::setw(5) << iomega << std::setw(2) << is 
                  << std::setw(3) << m1 << std::setw(3) << m2
                  << std::setw(20) << std::fixed << std::setprecision(12) 
                  << sigma_a[is][iomega][m_index].real()*GLC::Hartree_to_eV << " "
                  << std::setw(20) << std::fixed << std::setprecision(12) 
                  << sigma_a[is][iomega][m_index].imag()*GLC::Hartree_to_eV << '\n';
          }//m_index
        }//iomega
      }//is
      ofs_sig.close();

    }//ineq

    return;
  }
}
