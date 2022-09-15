#include "iQIST_narcissus.h"

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
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>
#include <cstdlib>   //Use exit function
#include <unistd.h>

extern "C" void iqist_narcissus_run_(); 

namespace DMFT
{
  void IQIST_NARCISSUS::output(
        const int char_step, const int DMFT_step, 
        const double mu, DMFT::input_info& in, DFT_output::atoms_info& atom, 
        DFT_output::KS_bands& band, const std::vector<double>& freq,
        std::vector<std::vector<std::vector<std::complex<double>>>>& Eimp,
        std::vector<std::vector<std::vector<std::vector<std::complex<double>>>>>& hyb_omega,
        DMFT::coulomb_tensor& Umat )
  {
    //======================================================
    //             energy unit: eV
    //======================================================
    debug::codestamp("IQIST_NARCISSUS::output");

    const int ineq_num = atom.inequ_atoms();
    const int symm = atom.local_sym();
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
      const int corr_L = atom.L(iatom);
      
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
      //        write solver.ctqmc.in
      //===========================================
      this->write_solver_ctqmc_in(current_dir+"/solver.ctqmc.in", m_tot, in);

      //===========================================
      //        write solver.eimp.in
      //===========================================
      std::vector<std::vector<std::complex<double>>> muvec = hoppinga;
      for(int is=0; is<nspin; is++)
        for(int m=0; m<m_tot; m++)
          muvec[is][m*m_tot+m] = muvec[is][m*m_tot+m]-mu;

      this->write_solver_eimp_in(
        current_dir+"/solver.eimp.in", 
        muvec, m_tot, symm, corr_L, nspin);

      //===========================================
      //        write solver.umat.in
      //===========================================
      this->write_solver_umat_in(current_dir+"/solver.umat.in", Umat.Coulomb_matrix()[ineq]);

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
      //                 << Sigma_ina[0][iomega][m*m_tot+m].real()*GLC::Hartree_to_eV
      //                 << std::setw(22) << std::fixed << std::setprecision(15) 
      //                 << Sigma_ina[0][iomega][m*m_tot+m].imag()*GLC::Hartree_to_eV;
      //       else
      //         ofs_Sig << std::setw(22) << std::fixed << std::setprecision(15) 
      //                 << Sigma_ina[is][iomega][m*m_tot+m].real()*GLC::Hartree_to_eV
      //                 << std::setw(22) << std::fixed << std::setprecision(15) 
      //                 << Sigma_ina[is][iomega][m*m_tot+m].imag()*GLC::Hartree_to_eV ;      
      //     }//m
      //   }//is
      //   ofs_Sig << std::endl;
      // }//iomega
      // ofs_Sig.close();

    }//ineq

    return;
  }

  void IQIST_NARCISSUS::write_solver_ctqmc_in(
        const std::string file, 
        const int nband,
        DMFT::input_info& in)
  {
    debug::codestamp("IQIST_NARCISSUS::write_solver_ctqmc_in");

    std::ofstream ofs(file.c_str(), std::ios::out);

    ofs << "###   setup general control flags  ###" << std::endl;
    ofs << "isscf = 1     #one-shot non-self-consistent scheme" << std::endl;
    ofs << "isscr = 1     #normal Hubbard model/Anderson impurity model" << std::endl;
    ofs << "isbnd = 2     #the bands are symmetrized according to symmetry matrix" << std::endl;
    ofs << "isspn = 1     #let spin up and spin down states evolve independently" << std::endl;
    ofs << "iswor = 2     #with worm algorithm, slow but more reliable" << std::endl;
    ofs << "isort = 3     #using singular value decomposition representation" << std::endl;
    ofs << "isobs = 1     #various physical observables(do nothing)" << std::endl;
    ofs << "issus = 1     #do not calculate charge/spin susceptibility" << std::endl;
    ofs << "isvrt = 1     #do not calculate two-particle green's functions" << std::endl;

    ofs << "\n###   setup common variables for quantum impurity model  ###" << std::endl;
    ofs << "niter = 1     #one-shot non-self-consistent scheme" << std::endl;
    ofs << "nband = " << nband << "     #number of correlated bands" << std::endl;
    ofs << "nspin = 2     #number of spin projections" << std::endl;
    ofs << "norbs = " << 2*nband << "     #number of correlated orbitals" << std::endl;
    ofs << "ncfgs = " << (int)std::pow(2,2*nband) << "     #number of atomic eigenstates" << std::endl;
    ofs << "mune = 0.0     #chemical potential" << std::endl;
    ofs << "beta = " << std::fixed << std::setprecision(9) 
        << *(double*)in.parameter("beta")/GLC::Hartree_to_eV 
        << "     #inverse temperature" << std::endl;
    
    ofs << "\n###   setup common variables for quantum impurity solver  ###" << std::endl;
    ofs << "lemax = 32        #maximum expansion order for legendre polynomial" << std::endl;
    ofs << "legrd = 20001     #number of mesh points for legendre polynomial" << std::endl;
    ofs << "svmax = 80        #maximum expansion order for svd polynomial" << std::endl;
    ofs << "svgrd = 2001      #number of mesh points for svd polynomial [-1,1]" << std::endl;
    ofs << "mkink = 1024      #maximum perturbation expansion order" << std::endl;
    ofs << "mfreq = " << *(int*)in.parameter("n_omega")   //100eV
        << "        #maximum number of matsubara frequency points" << std::endl;
    ofs << "nfreq = " << *(int*)in.parameter("n_omega")/4         //~25eV
        << "        #number of sampled matsubara frequency points" << std::endl;
    ofs << "ntime = " << *(int*)in.parameter("n_tau") 
        << "        #number of time slices" << std::endl;
    ofs << "nflip = 20000     #flip period for spin up and spin down states" << std::endl;
    ofs << "ntherm = 200000   #flip period for spin up and spin down states" << std::endl;
    ofs << "nsweep = " << *(long long*)in.parameter("mc_step") 
        << "        #number of Monte Carlo sweeping steps" << std::endl;
    ofs << "nwrite = 2000000  #output period" << std::endl;
    ofs << "nclean = 100000   #clean update period" << std::endl;
    ofs << "nmonte = 10       #how often to sample the observables" << std::endl;
    ofs << "ncarlo = 10       #how often to sample the observables" << std::endl;

    ofs.close();

    return;
  }

  void IQIST_NARCISSUS::write_solver_eimp_in(
        const std::string file, 
        std::vector<std::vector<
        std::complex<double>>>& muvec,
        const int nband, const int symm,
        const int corr_L, const int nspin )
  {
    debug::codestamp("IQIST_NARCISSUS::write_solver_eimp_in");

    std::ofstream ofs(file.c_str(), std::ios::out);

    int count=1;
    for(int is=0; is<2; is++)
    {
      for(int m=0; m<nband; m++)
      {
        ofs << std::setw(2) << count;
        if(nspin==1)
        {
          ofs << std::setw(22) << std::fixed << std::setprecision(15)
              << (muvec[0][m*nband+m].real())*GLC::Hartree_to_eV;
          
          if(symm==1 && corr_L==2) //cubic symmetry, d orbital
          {
            if(m==2 || m==4) ofs << std::setw(3) << 3 << std::endl;
            else ofs << std::setw(3) << 1 << std::endl;
          }
          else if(symm==2 || symm==3)//t2g or eg only
          {
            ofs << std::setw(3) << 1 << std::endl;
          }
          else //only spin symmetry
          {
            ofs << std::setw(3) << m+1 << std::endl;
          }
        }
        else if(nspin==2)
        {
          ofs << std::setw(22) << std::fixed << std::setprecision(15)
              << muvec[is][m*nband+m]*GLC::Hartree_to_eV;
          
          if(symm==1 && corr_L==2) //cubic symmetry, d orbital
          {
            if(m==2 || m==4) ofs << std::setw(3) << 3 + 5*is << std::endl;
            else ofs << std::setw(3) << 1 + 5*is << std::endl;
          }
          else if(symm==2)//t2g only
          {
            ofs << std::setw(3) << 1 + 3*is << std::endl;
          }
          else if(symm==3)//eg only
          {
            ofs << std::setw(3) << 1 + 2*is << std::endl;
          }
          else //no symmetry
          {
            ofs << std::setw(3) << m+1 << std::endl;
          }
        }
        
        count++;
      }
    }

    ofs.close();

    return;
  }

  void IQIST_NARCISSUS::write_solver_umat_in(
         const std::string file, 
         std::vector<std::vector<std::vector<
         std::vector<double>>>>& Umat )
  {
    debug::codestamp("IQIST_NARCISSUS::write_solver_umat_in");

    std::ofstream ofs(file.c_str(), std::ios::out);
    
    for(int iorb1=0; iorb1<Umat.size(); iorb1++)
      for(int iorb2=0; iorb2<Umat.size(); iorb2++)
        ofs << std::setw(2) << iorb1+1 << std::setw(4) << iorb2+1 
            << std::setw(12) << std::fixed << std::setprecision(6)
            << Umat[iorb1][iorb2][iorb2][iorb1]*GLC::Hartree_to_eV << std::endl;

    ofs.close();
    return;
  }

  void IQIST_NARCISSUS::read_self_energy(
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

    int count = 0;

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

  void IQIST_NARCISSUS::impurities_solving(
          const int char_step,
          const int DMFT_step,
          DFT_output::atoms_info& atom )
  {
    debug::codestamp("IQIST_NARCISSUS::impurities_solving");

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

      MPI_Barrier(MPI_COMM_WORLD);   //Blocks until all processes reach here

      std::string date;
      timer::get_date_time(date);
      GLV::ofs_running << "  impurity" << ineq << "     " << date;
      GLV::ofs_running.flush();

      double time, seconds;
      int hours, minutes;
      timer::timestamp(time);

      iqist_narcissus_run_();

      timer::get_date_time(date);
      GLV::ofs_running << "       " << date;

      timer::get_time(time, seconds, minutes, hours);

      GLV::ofs_running << "             " << hours << "h " 
                           << minutes << "m "
                           << (int)seconds << "s" << std::endl;

      ierr = chdir("../../../../");
      if(ierr != 0){
        std::cout << "Process " << mpi_rank() << " fails to return to root directory " << std::endl;
      }

      MPI_Barrier(MPI_COMM_WORLD);  //Blocks until all processes reach here
    }//ineq

    return;
  }

}
