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

    static void strtolower(char *in_str, char *out_str);

    void test();    //test whether input parameter is correct
                    
    void out();     //test whether reading worked corrected            
    
    const void* parameter(char* keyword);    //   interfaces

    const std::vector<double>& window(){return this->E_window;}

    private:
    // Input parameter
    int flag_DFT_solver;         //1:aims, 2:ABACUS
    double temperature;          //unit:Kelvin; default:300
    int flag_magnetism=-1;       //AFM:1, FM:2, paramagenetism:3, none:4
    int n_tau;                   //number of segments on the interval [0,beta] when using CTQMC impuroty solver;default:500
    int n_omega;                 //number of Matsubara frequency points
    //ALPS-CTHYB:1;  LPS-CTHYB-SEGMENT:2; (not supported now)
    //PACS:3;  Rutgers-CTHYB:4; iQIST:5; default:3
    int flag_impurity_solver;    
    int flag_double_counting;    //nominal:1;   default:1
    int charge_step_max;             //number of the maximum charge self-consistent interation step
    int DMFT_step_max;               //number of the maximum DMFT interation step under each charge step
    long long MC_step;           //numuber of MC sampling steps; default:1000
    std::vector<double> E_window;  //Energy window for hybridization function; default E_window[0]=-2.0, E_window[1]=2.0;

    //parameters determined by input parameters
    double beta;                 //1/(k_b*T); unit:1/Hartree

  };
}
