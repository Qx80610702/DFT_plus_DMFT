#pragma once

#include <string>
#include <vector>

namespace DMFT
{
  class input_info
  {
    public:
    input_info(){;}
    ~input_info(){;}

    bool read();

    bool read_parameter(char* keyword, std::vector<std::string>& val, int count=1);

    static std::string strtolower(char *in_str);
    static std::string strtolower(std::string in_str);

    void test();    //test whether input parameter is correct
                    
    void out();     //test whether reading worked corrected            
    
    const void* parameter(char* keyword);    //   interfaces

    const std::vector<double>& projection_window(){return this->proj_window;}
    const std::vector<double>& spectra_window(){return this->DOS_window;}

    private:
    // Input parameter
    int flag_DFT_solver;            //1:aims, 2:ABACUS
    double temperature;             //unit:Kelvin; default:300
    int flag_magnetism;             //AFM:1, FM:2, paramagenetism:3, none:4
    int n_tau;                      //number of segments on the interval [0,beta] when using CTQMC impuroty solver;default:500
    int n_omega;                    //number of Matsubara frequency points
    //ALPS-CTHYB:1;  LPS-CTHYB-SEGMENT:2; (not supported now)
    //PACS:3;  Rutgers-CTHYB:4; iQIST:5; default:4
    int flag_impurity_solver;    
    int flag_double_counting;       //nominal:1;   default:1
    int charge_step_max;            //number of the maximum charge self-consistent interation loop; default:1
    int DMFT_step_max;              //number of the maximum DMFT interation step under each charge loop; default:1
    int DFT_step_max;               //number of the maximum DFT interation step under each charge loop; default:1
    long long MC_step;              //numuber of MC sampling steps; default:1000
    std::vector<double> proj_window;  //Energy window for hybridization function; default proj_window[0]=-5.0, proj_window[1]=5.0;
    std::vector<double> DOS_window;   //Energy window for DOS; default proj_window[0]=proj_window[0], proj_window[1]=proj_window[1];
    bool hyb_func;                  //Whether hybrid functional is used or not; default:false
    double hyf_xc_alpha;            //The mixing factor in hybrid funtional; default:0.25;
    double charge_mix_param;         //The charge density mixing parameter; default:0.05
    double delta_sigma;             //The convergency criteria of self-energy (unit eV); default:0.1
    double delta_rho;               //The convergency criteria of the charge density; default: 1.0e-3
    int calculation_type;           //DFT+DMFT scf:0; spectra:1; default:0
    int mixing_step;                //Number of past iteration include in the density mixing; default:8

    std::string dft_solver_exe;     //the executable of DFT solver

    //restart
    bool flag_restart;              
    int start_charge_step;          //The starting step of charge loop; default:1
    int start_DMFT_step;            //The stsrting step of DMFT loop; default:1
    int last_charge_step;           //The last one step before the starting step of charge loop; default:1
    int last_DMFT_step;             //The last one tsep before the stsrting step of DMFT loop; default:1

    //parameters determined by input temperature
    double beta;                 //1/(k_b*T); unit:1/Hartree

  };
}
