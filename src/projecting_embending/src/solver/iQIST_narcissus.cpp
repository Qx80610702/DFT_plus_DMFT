#include "iQIST_narcissus.h"

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
  void IQIST_NARCISSUS::output(const int istep, 
        const double mu, DMFT::input_info& in, DFT_output::atoms_info& atom, 
        DFT_output::KS_bands& band, const std::vector<double>& freq,
        std::vector<std::vector<std::vector<std::complex<double>>>>& Eimp,
        std::vector<std::vector<std::vector<std::vector<std::complex<double>>>>>& Gf_in,
        std::vector<std::vector<std::vector<std::vector<std::complex<double>>>>>& Weiss,
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
      const int corr_L = atom.L(iatom);
      
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
        << *(double*)in.parameter("beta")/Hartree_to_eV 
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
        std::vector<std::vector<std::complex<double>>>& muvec,
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
              << muvec[0][m*nband+m].real()*Hartree_to_eV;
          
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
              << muvec[is][m*nband+m]*Hartree_to_eV;
          
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
            << Umat[iorb1][iorb2][iorb2][iorb1]*Hartree_to_eV << std::endl;

    ofs.close();
    return;
  }

  void IQIST_NARCISSUS::read_last_step(
        const int istep, 
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

    //directory impurity 
    std::string dir_impurity_solving = "impurity_solving";

    //directory step+num
    std::stringstream step_dir_ss;
    step_dir_ss << "/step" << istep-1;
    std::string step_dir= step_dir_ss.str();

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
      current_dir_ss << dir_impurity_solving << step_dir << site_dir;
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
