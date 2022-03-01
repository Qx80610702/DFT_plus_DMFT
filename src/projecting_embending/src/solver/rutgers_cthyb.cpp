#include "rutgers_cthyb.h"

#include "../debug.h"
#include "../constants.h"

#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>
#include <cstdlib>

namespace DMFT
{
  void Rutgers_CTHYB::output(
        const int char_step, const int DMFT_step, 
        const double mu, DMFT::input_info& in, 
        DFT_output::atoms_info& atom, DFT_output::KS_bands& band, const std::vector<double>& freq,
        std::vector<std::vector<std::vector<std::complex<double>>>>& Eimp,
        std::vector<std::vector<std::vector<std::vector<std::complex<double>>>>>& Gf_in,
        std::vector<std::vector<std::vector<std::vector<std::complex<double>>>>>& Weiss,
        std::vector<std::vector<std::vector<std::vector<std::complex<double>>>>>& hyb_omega,
        DMFT::coulomb_tensor& Umat)
  {
    //======================================================
    //        energy unit: eV
    //======================================================
    debug::codestamp("Rutgers_CTHYB::output");

    const int ineq_num = atom.inequ_atoms();
    const std::vector<int>& norb_sub = atom.iatom_norb();
    const bool soc = band.soc();
    const int nspin = band.nspins();
    const int ntau = *(int*)in.parameter("n_tau");
    const int nomega = *(int*)in.parameter("n_omega");

    //Create directory impurity 
    std::string dir_impurity_solving = "dmft_solving";
    std::stringstream make_dir1;
    make_dir1 << "test -d " << dir_impurity_solving << " || mkdir " << dir_impurity_solving;
    system(make_dir1.str().c_str());

    std::stringstream char_dir_ss;
    char_dir_ss << "/charge_step" << char_step;
    std::string char_step_dir= char_dir_ss.str();

    std::stringstream step_dir_ss;
    step_dir_ss << "/dmft_step" << DMFT_step;
    std::string step_dir= step_dir_ss.str();
    std::stringstream make_dir2;
    make_dir2 << "test -d " << dir_impurity_solving << char_step_dir << step_dir
            << " || mkdir " << dir_impurity_solving << char_step_dir << step_dir;
    system(make_dir2.str().c_str());

    for(int ineq=0; ineq<ineq_num; ineq++)
    {
      const int iatom = atom.ineq_iatom(ineq);
      const int m_tot=norb_sub[iatom];

      const std::vector<std::vector<std::vector<std::complex<double>>>>& 
            hyb = hyb_omega[ineq];
      
      std::vector<std::vector<std::complex<double>>>
            Eimpa = Eimp[ineq];

      const std::vector<std::vector<std::vector<std::vector<double>>>>&
            Utensor = Umat.Coulomb_matrix()[ineq];
      
      const std::vector<std::vector<std::vector<std::complex<double>>>>&
            Gf_ina = Gf_in[ineq];


      std::stringstream site_dir_ss;
      site_dir_ss << "/impurity" << ineq;
      std::string site_dir= site_dir_ss.str();
      std::stringstream make_dir3;
      make_dir3 << "test -d " << dir_impurity_solving << char_step_dir << step_dir << site_dir
              << " || mkdir " << dir_impurity_solving << char_step_dir << step_dir << site_dir;
      system(make_dir3.str().c_str());

      std::stringstream current_dir_ss;
      current_dir_ss << dir_impurity_solving << char_step_dir << step_dir << site_dir;
      std::string current_dir = current_dir_ss.str();

      //===========================================
      //           write PARAMS
      //===========================================
      double average_Eimp = 0.0;
      int count=0;
      for(int is=0; is<nspin; is++)
      {
        for(int m=0; m<m_tot; m++)
        {
          average_Eimp += Eimpa[is][m*m_tot+m].real();
          count++;
        }
      }
      average_Eimp /= count;

      for(int is=0; is<nspin; is++)
        for(int m=0; m<m_tot; m++)
          Eimpa[is][m*m_tot+m] -= average_Eimp;

      this->write_params(in, current_dir+"/PARAMS", mu-average_Eimp);

      //===========================================
      //           write  delta.dat
      //===========================================
      this->write_delta_omega(in, current_dir+"/delta.dat",
                          nspin, nomega, m_tot, freq, hyb);

      //===========================================
      //           write  ctqmc.cix
      //===========================================
      this->write_cix(nspin, m_tot, current_dir+"/ctqmc.cix", Eimpa, Utensor);

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
          }//m_index
        }//is
        ofs_gf << '\n';
      }//iomega
      ofs_gf.close();

    }//ineq

    return;
  }

  void Rutgers_CTHYB::write_params(DMFT::input_info& in, 
                  const std::string file, const double mu)
  {
    debug::codestamp("Rutgers_CTHYB::write_params");

    std::ofstream ofs(file.c_str(), std::ios::out);

    ofs << "nom  " << *(int*)in.parameter("n_omega") 
        << "    #Number of Matsubara frequency points sampled\n";
    ofs << "Ntau " << *(int*)in.parameter("n_tau") 
        << "    #Number of imaginary time τ points to spline input hybridization function\n";
    ofs << "beta " << std::fixed << std::setprecision(9) 
        << *(double*)in.parameter("beta")/Hartree_to_eV << "    #Inverse temperature\n";  //Hartree to eV
    ofs << "U 0.0    #Coulomb repulsion (F0), This information is wrtten in cix file\n";
    ofs << "mu " << std::setprecision(15) << mu*Hartree_to_eV 
        << "    #Chemical potential, This information is wrtten in cix file\n";
    ofs << "warmup 100000    #Warmup number of QMC steps\n";
    ofs << "M " << *(long long*)in.parameter("mc_step") << "    #Total number of Monte Carlo steps\n";
    ofs << "tsample 20    #How often to record measurements\n";
    ofs << "Nmax 200    #How often to record measurements\n";
    ofs << "GlobalFlip 500000    #how often to perform global flip\n";
    ofs << "svd_lmax 25    #We will use SVD basis to expand G, with this cutoff\n";
    ofs << "svd_Ntau " << *(int*)in.parameter("n_tau") + 1 << '\n';
    ofs << "aom 1    #number of frequency points to determin high frequency tail\n";
    ofs << "mode SH    #We will use self-energy sampling, and Hubbard I tail\n";
    ofs << "Delta delta.dat    #Input bath function hybridization\n";
    ofs << "cix ctqmc.cix    #Input file with atomic state\n";
    ofs << "Gf Gw.out    #Filename of the output Green’s function G(iw)\n";
    ofs << "Sig Sigma.dat    #Filename of the output self-energy \\sigma(iw)";

    ofs.close();
    return;
  }

  void Rutgers_CTHYB::write_delta_omega(DMFT::input_info& in, 
          const std::string file, const int nspin, const int nomega,
          const int m_tot, const std::vector<double>& freq,
          const std::vector<std::vector<
          std::vector<std::complex<double>>>>& hyb_omega)
  {
    debug::codestamp("Rutgers_CTHYB::write_params");

    std::ofstream ofs(file.c_str(), std::ios::out);

    for(int iomega=0; iomega<nomega; iomega++)
    {
      ofs << std::setw(22) << std::fixed << std::setprecision(15)
          << Hartree_to_eV*freq[iomega];
      for(int is=0; is<2; is++)
      {
        for(int m=0; m<m_tot; m++)
        {
          if(nspin==1)
            ofs << std::setw(22) << std::fixed << std::setprecision(15)
            << Hartree_to_eV*hyb_omega[0][iomega][m*m_tot+m].real() 
            << std::setw(22) << std::fixed << std::setprecision(15) 
            << Hartree_to_eV*hyb_omega[0][iomega][m*m_tot+m].imag();
          else
            ofs << std::setw(22) << std::fixed << std::setprecision(15)
            << Hartree_to_eV*hyb_omega[is][iomega][m*m_tot+m].real()
            << std::setw(22) << std::fixed << std::setprecision(15) 
            << Hartree_to_eV*hyb_omega[is][iomega][m*m_tot+m].imag();     
        }//m
      }//is
      ofs << '\n';
    }//iomega

    ofs.close();
    return;
  }

  void Rutgers_CTHYB::write_cix(const int nspin, const int m_tot, const std::string file, 
          const std::vector<std::vector<std::complex<double>>>& Esplit,
          const std::vector<std::vector<std::vector<std::vector<double>>>>& Utensor)
  {
    debug::codestamp("Rutgers_CTHYB::write_cix");

    int orb_count;
    std::ofstream ofs(file.c_str(), std::ios::out);

    ofs << "# CIX file for ctqmc!\n";
    ofs << "# cluster_size, number of states, number of baths, maximum_matrix_size\n";
    ofs << "1 " << (int)std::pow(2,2*m_tot) << ' ' << 2*m_tot << " 1\n";

    ofs << "# baths, dimension, symmetry, global flip\n";
    orb_count=0;
    for(int is=0; is<2; is++)
      for(int m=0; m<m_tot; m++)
      {
        ofs << orb_count << "  1  " << orb_count << "  " ;
        if(is==0) ofs << orb_count << '\n';
        else if(is==1) ofs << orb_count-m_tot << '\n';
        orb_count++;
      }

    ofs << "# cluster energies for unique baths, eps[k]\n";
    for(int is=0; is<2; is++)
      for(int m=0; m<m_tot; m++)
        if(nspin==2)
          ofs << std::setw(22) << std::fixed << std::setprecision(15) 
              << (Esplit[is][m*m_tot+m].real())*Hartree_to_eV;
        else if(nspin==1)
          ofs << std::setw(22) << std::fixed << std::setprecision(15) 
              << (Esplit[0][m*m_tot+m].real())*Hartree_to_eV;
    ofs << '\n';

    ofs << "#   N   K   Sz size\n";
    for(int dim=0; dim<(int)std::pow(2,2*m_tot); dim++)
    {
      ofs << std::setw(7) << dim+1 
          << std::setw(4) << site_occ(dim)
          << "  0  " << std::setw(6) << std::fixed << std::setprecision(3) << site_Sz(dim, 2*m_tot)
          << "  1  ";
      
      for(int ibath=0; ibath<2*m_tot; ibath++)
        ofs << std::setw(7) << F_dagger_state(dim, ibath, 2*m_tot);

      ofs << std::setw(25) << std::fixed << std::setprecision(15) 
          << state_energy(dim, m_tot, nspin, Utensor, Esplit) << "  0\n"; 
    }

    ofs << "# Matrix elements\n";
    for(int dim=0; dim<(int)std::pow(2,2*m_tot); dim++)
    {
      for(int ibath=0; ibath<2*m_tot; ibath++)
      {
        int end_state = F_dagger_state(dim, ibath, 2*m_tot);
        if(end_state==0)
        {
          ofs << std::setw(7) << dim+1
              << "   0   0  0\n";
        }
        else
        {
          ofs << std::setw(7) << dim+1
              << std::setw(7) << end_state
              << "   1  1"
              << std::setw(11) << std::fixed << std::setprecision(6)
              << F_dagger_sign(dim, ibath) << '\n';
        }
      }
    }

    ofs << "HB2   # Hubbard-I is used to determine high-frequency\n";

    ofs << "# UCoulomb : (m1,s1) (m2,s2) (m3,s2) (m4,s1)  Uc[m1,m2,m3,m4]\n";
    for(int ibath1=0; ibath1<2*m_tot; ibath1++)
    {
      for(int ibath2=0; ibath2<2*m_tot; ibath2++)
      {
        for(int ibath3=0; ibath3<2*m_tot; ibath3++)
        {
          for(int ibath4=0; ibath4<2*m_tot; ibath4++)
          {
            if(std::fabs(Utensor[ibath1][ibath2][ibath3][ibath4])>1.0e-3)
            {
              ofs << std::setw(3) << ibath1
                  << std::setw(3) << ibath2
                  << std::setw(3) << ibath3
                  << std::setw(3) << ibath4
                  << std::setw(12) << std::fixed << std::setprecision(6) 
                  << Utensor[ibath1][ibath2][ibath3][ibath4]*Hartree_to_eV << '\n';
            }
          }
        }
      }
    }

    ofs << "# number of operators needed\n";
    ofs << "0\n";

    ofs.close();
    return;
  }

  int Rutgers_CTHYB::site_occ(int n)
  {
    int count = 0;
    while (n)
    {
      count += n & 1;
      n >>= 1;
    }
    return count;
  }

  double Rutgers_CTHYB::site_Sz(int n, const int nbath)
  {
    double Sz = 0.0;

    for(int i=0; i<nbath; i++)
    {
      if(i<nbath/2) Sz += 0.5*(n & 1);
      else Sz -= 0.5*(n & 1);
      n >>= 1;
    }

    return Sz;
  }

  int Rutgers_CTHYB::F_dagger_state(int n, const int ibath, const int nbath)
  {
    int count=0;

    for(int i=0; i<nbath; i++)
    {
      if(i==ibath)
      {
        if(n&1) return 0;
        else count += (int)std::pow(2, i);
      }
      count += (int)std::pow(2, i)*(n&1);
      n >>= 1;
    }

    return count+1;
  }

  double Rutgers_CTHYB::F_dagger_sign(int n, const int ibath)
  {
    int count = 0;
    for(int j=0; j<ibath; j++)
    {
        count += ((n >> j) & 1);
    }
    return  std::pow(-1, count);
  }

  double Rutgers_CTHYB::state_energy(
    const int n, const int m_tot, const int nspin, 
    const std::vector<std::vector<std::vector<std::vector<double>>>>& Utensor, 
    const std::vector<std::vector<std::complex<double>>>& Esplit)
  {
    double Etot=0.0;

    int ibath=0;
    for(int is=0; is<2; is++)
    {
      for(int m=0; m<m_tot; m++)
      {
        if(nspin==2)
          Etot += ((n>>ibath)&1)*(Esplit[is][m*m_tot+m].real());
        else
          Etot += ((n>>ibath)&1)*(Esplit[0][m*m_tot+m].real());
        ibath++;
      }
    }

    for(int ibath1=0; ibath1<2*m_tot; ibath1++)
    {
      int is1=ibath1/m_tot;
      int m1 = ibath1<m_tot ? ibath1 : (ibath1-m_tot);
      int occ1 = ((n>>ibath1)&1);
      for(int ibath2=0; ibath2<ibath1; ibath2++)
      {
        int is2=ibath2/m_tot;
        int m2 = ibath2<m_tot ? ibath2 : (ibath2-m_tot);
        int occ2 = ((n>>ibath2)&1);

        // if(is1==is2)
        // {
        //   Etot += (U-3.0*J)*occ1*occ2;
        // }
        // else
        // {
        //   if(m1==m2) Etot += U*occ1*occ2;
        //   else Etot += (U-2.0*J)*occ1*occ2;
        // }
        Etot += Utensor[ibath1][ibath2][ibath2][ibath1]*occ1*occ2;
      }
    }

    return Etot*Hartree_to_eV;
  }

  void Rutgers_CTHYB::read_last_step(
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
    debug::codestamp("Rutgers_CTHYB::read_last_step");

    const int ineq_num = atom.inequ_atoms();
    const std::vector<int>& norb_sub = atom.iatom_norb();
    const bool soc = band.soc();
    const int nspin = band.nspins();
    const int magnetism = *(int*)in.parameter("magnetism");

    double omega, real, imag;
    std::string str_tmp;

    std::string dir_dmft_solving = "dmft_solving";

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
      std::string Gf_file = current_dir+"/Gw.out";
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
