#pragma once

#include <string>
#include <vector>

namespace DFT_output
{
  class KS_bands
  {
    public:
    KS_bands(){;}
    ~KS_bands(){;}

    bool read();
    void out();       //test whether reading worked corrected

    //=================
    //  interfaces
    //=================
    double& Fermi_level(){return Efermi;}
    int nband(){return nbands;}
    int nspins(){return nspin;}
    int nk(){return kpoints;}
    bool soc(){return flag_SOC;}
    double n_electrons(){return tot_elec_num;}
    const std::vector<std::vector<std::vector<double>>>& enk(){return eigen_values;}
    const std::vector<double>& kweight(){return k_weight;}
    const std::vector<std::vector<std::vector<double>>>& dft_occ(){return DFT_occ_numbers;}

    private:

    double Efermi;                  //Fermi level  given by DFT
    int nbands;                     //Total number of bands of each spin
    int kpoints;
    bool flag_SOC;                  //true for SOC, false for non_SOC
    double tot_elec_num;               //total number of electrons of the system

    //=============================
    //   non_SOC case
    //=============================
    int nspin;                   //Total number of spin
    std::vector<std::vector<std::vector<double>>> eigen_values;      //eigen_values[i_spin][ik][i_bands]

    std::vector<std::vector<std::vector<double>>> DFT_occ_numbers;   //DFT_occ_numbers[i_spin][ik][i_bands]

    std::vector<double> k_weight;       //k_weight[ik]

  };
  
}
