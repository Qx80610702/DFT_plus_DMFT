#include "input.h"
#include "../constants.h"
#include "../debug.h"
#include "../mpi_environment.h"
#include "../global_variables.h"

#include <fstream>
#include <iostream>
#include <cstring>
#include <iomanip>   //use setw and setprecison function
#include <cstdlib> 

namespace DMFT
{
  bool input_info::read()
  {
    debug::codestamp("input_info::read");

    GLV::ofs_running << "Reading DMFT.in ......" << std::endl;
    GLV::ofs_running << "=========== All parameters used in the running time ===========" << std::endl;

    //calculation
    try {
      std::vector<std::string> str_val;
      this->read_parameter("calculation", str_val);

      if(std::strcmp("scf", str_val[0].c_str())==0)
        this->calculation_type = 0;
      else if(std::strcmp("spectra", str_val[0].c_str())==0)
        this->calculation_type = 1;
      else{
        GLV::ofs_error << "Unsupported paramater for the key calculationi" << std::endl;
        std::exit(EXIT_FAILURE);
      }  
    }
    catch (const std::string messg) {
      GLV::ofs_error << messg << std::endl;
      std::exit(EXIT_FAILURE);
    }
    catch(const bool not_given){
      // GLV::ofs_running << "Warning: calculation is not given and set default value scf" << std::endl;
      this->calculation_type = 0;
    }

    GLV::ofs_running << "calculation  ";
    if(this->calculation_type==0) GLV::ofs_running << "scf" << std::endl;
    else if(this->calculation_type==1) GLV::ofs_running << "spectra" << std::endl;

    //dft_solver
    try {
      std::vector<std::string> str_val;
      this->read_parameter("dft_solver", str_val);

      if(std::strcmp("aims",str_val[0].c_str())==0){
        this->flag_DFT_solver = 1;
        #ifndef __FHIaims
          GLV::ofs_error << "FHI-aims has not been installed!!!  ";
          GLV::ofs_error << "Suggestion:Install FHI-aims and then re-compile the codes." << std::endl;
          std::exit(EXIT_FAILURE);
        #endif
      }
      else if(std::strcmp("abacus", str_val[0].c_str())==0){
        this->flag_DFT_solver = 2;
        #ifndef __ABACUS
          GLV::ofs_error << "ABACUS has not been installed!!!   ";
          GLV::ofs_error << "Suggestion:Install ABACUS and then re-compile the codes." << std::endl;
          std::exit(EXIT_FAILURE);
        #endif
      }
      else
      {
        GLV::ofs_error << "unsupported value of DFT_solver" << std::endl;
        std::exit(EXIT_FAILURE);
      }
    }
    catch (const std::string messg) {
      GLV::ofs_error << messg << std::endl;
      std::exit(EXIT_FAILURE);
    }
    catch(const bool not_given){
      // GLV::ofs_running << "Warning: dft_solver is not given and set default value aims" << std::endl;
      this->flag_DFT_solver = 1;   //FHI-aims
      #ifndef __FHIaims
        GLV::ofs_error << "FHI-aims has not been installed!!!  ";
        GLV::ofs_error << "Suggestion:Install FHI-aims and then re-compile the codes." << std::endl;
        std::exit(EXIT_FAILURE);
      #endif
    }

    GLV::ofs_running << "dft_solver  ";
    if(this->flag_DFT_solver==1) GLV::ofs_running << "aims" << std::endl;
    else if(this->flag_DFT_solver==2) GLV::ofs_running << "abacus" << std::endl;

    //temperature
    try {
      std::vector<std::string> str_val;
      this->read_parameter("temperature", str_val);

      this->temperature = atof(str_val[0].c_str());
    }
    catch (const std::string messg) {
      GLV::ofs_error << messg << std::endl;
      std::exit(EXIT_FAILURE);
    }
    catch(const bool not_given){
      // GLV::ofs_running << "Warning: temperature is not given and set default value 300 Kelven" << std::endl;
      this->temperature = 300.0;
    }
    GLV::ofs_running << "temperature  " << this->temperature << std::endl;

    this->beta = GLC::Hartree_to_eV/(this->temperature*GLC::K_BOLTZMAN_EV);
    // this->n_tau = 25/(this->temperature*K_BOLTZMAN_EV) > 500 ? 25/(this->temperature*K_BOLTZMAN_EV) : 500;
    this->n_tau = 25/(this->temperature*GLC::K_BOLTZMAN_EV);
    // this->n_omega = this->n_tau > 500 ? 500 : this->n_tau;
    //the energy maximum is 100eV, which may be a safe value to guarantee the Asymptotic behavior of hybridization function
    this->n_omega = 50/(GLC::PI*this->temperature*GLC::K_BOLTZMAN_EV); 
    // this->n_omega = 500; // ? 500 : this->n_tau;
    if(this->n_tau%2==1) this->n_tau++;   //n_tau+1 must be odd for Simpson integral
    if(this->n_omega%2==1) this->n_omega++;

    //magnetism
    try {
      std::vector<std::string> str_val;
      this->read_parameter("magnetism", str_val);
    
      if(std::strcmp("afm",str_val[0].c_str())==0) this->flag_magnetism = 1;
      else if(std::strcmp("fm",str_val[0].c_str())==0) this->flag_magnetism = 2;
      else if(std::strcmp("para",str_val[0].c_str())==0) this->flag_magnetism = 3;
      else if(std::strcmp("none",str_val[0].c_str())==0) this->flag_magnetism = 4;
      else
      {
        GLV::ofs_error << "unsupported value of magnetism" << std::endl;
        std::exit(EXIT_FAILURE);
      }
    }
    catch (const std::string messg) {
      GLV::ofs_error << messg << std::endl;
      std::exit(EXIT_FAILURE);
    }
    catch(const bool not_given){
      this->flag_magnetism = 3;
      // GLV::ofs_error << "Error: parameter magnetism must be given!!!" << std::endl;
      // std::exit(EXIT_FAILURE);
    }
    GLV::ofs_running << "magnetism  ";
    if(this->flag_magnetism==1) GLV::ofs_running << "afm" << std::endl;
    else if(this->flag_magnetism==2) GLV::ofs_running << "fm" << std::endl;
    else if(this->flag_magnetism==3) GLV::ofs_running << "para" << std::endl;
    else if(this->flag_magnetism==4) GLV::ofs_running << "none" << std::endl;

    //impurity_solver
    try {
      std::vector<std::string> str_val;
      this->read_parameter("impurity_solver", str_val);

      if(std::strcmp("alps-cthyb",str_val[0].c_str())==0) 
      {
        this->flag_impurity_solver = 1;
        GLV::ofs_error << "ALPS-CTHYB is unsupported now" << std::endl;
        std::exit(EXIT_FAILURE);
      }
      else if(std::strcmp("alps-cthyb-segment",str_val[0].c_str())==0)
      {
        this->flag_impurity_solver = 2;
        GLV::ofs_error << "ALPS-CTHYB-SEGMENT is unsupported now" << std::endl;
        std::exit(EXIT_FAILURE);
      }
      else if(std::strcmp("pacs", str_val[0].c_str())==0) this->flag_impurity_solver = 3;
      else if(std::strcmp("rutgers-cthyb",str_val[0].c_str())==0) this->flag_impurity_solver = 4;
      else if(std::strcmp("iqist",str_val[0].c_str())==0)
      {
        this->flag_impurity_solver = 5;
      }
      else
      {
        GLV::ofs_error << "unsupported value of impurity_solver " << str_val[0] << std::endl;
        std::exit(EXIT_FAILURE);
      }
    }
    catch (const std::string messg) {
      GLV::ofs_error << messg << std::endl;
      std::exit(EXIT_FAILURE);
    }
    catch(const bool not_given){
      GLV::ofs_running << "Warning: impurity_solver is not given and set default value Rutgers" << std::endl;
      this->flag_impurity_solver = 4;
    }
    GLV::ofs_running << "impurity_solver  ";
    if(this->flag_impurity_solver==1) GLV::ofs_running << "alps-cthyb" << std::endl;
    else if(this->flag_impurity_solver==2) GLV::ofs_running << "alps-cthyb-segment" << std::endl;
    else if(this->flag_impurity_solver==3) GLV::ofs_running << "pacs" << std::endl;
    else if(this->flag_impurity_solver==4) GLV::ofs_running << "rutgers-cthyb" << std::endl;
    else if(this->flag_impurity_solver==4) GLV::ofs_running << "iqist" << std::endl;

    //double_counting
    try {
      std::vector<std::string> str_val;
      this->read_parameter("double_counting", str_val);

      if(std::strcmp("nominal",str_val[0].c_str())==0) this->flag_double_counting = 1;
      else
      {
        GLV::ofs_error << "unsupported value of double_counting" << std::endl;
        std::exit(EXIT_FAILURE);
      }

    }
    catch (const std::string messg) {
      GLV::ofs_error << messg << std::endl;
      std::exit(EXIT_FAILURE);
    }
    catch(const bool not_given){
      // GLV::ofs_running << "Warning: double_counting is not given and set default value nominal" << std::endl;
      this->flag_double_counting = 1;
    }

    //max_charge_step
    try {
      std::vector<std::string> str_val;
      this->read_parameter("max_charge_step", str_val);

      this->charge_step_max = atoi(str_val[0].c_str());
    }
    catch (const std::string messg) {
      GLV::ofs_error << messg << std::endl;
      std::exit(EXIT_FAILURE);
    }
    catch(const bool not_given){
      // GLV::ofs_running << "Warning: max_charge_step is not given and set default value 1, i.e., one-shoot DFT+DMFT." << std::endl;
      this->charge_step_max = 1;
    }
    GLV::ofs_running << "max_charge_step  " << this->charge_step_max << std::endl;

    //delta_sigma
    try {
      std::vector<std::string> str_val;
      this->read_parameter("delta_sigma", str_val);

      this->delta_sigma = atof(str_val[0].c_str());
    }
    catch (const std::string messg) {
      GLV::ofs_error << messg << std::endl;
      std::exit(EXIT_FAILURE);
    }
    catch(const bool not_given){
      // GLV::ofs_running << "Warning: max_charge_step is not given and set default value 1, i.e., non self-consitent DFT+DMFT." << std::endl;
      this->delta_sigma = 0.1;
    }
    GLV::ofs_running << "delta_sigma  " 
            << std::setprecision(3) << std::setiosflags(std::ios::scientific)
            << this->delta_sigma << std::endl;

    //charge_mix_beta
    try {
      std::vector<std::string> str_val;
      this->read_parameter("charge_mix_beta", str_val);

      this->charge_mix_beta = atof(str_val[0].c_str());
    }
    catch (const std::string messg) {
      GLV::ofs_error << messg << std::endl;
      std::exit(EXIT_FAILURE);
    }
    catch(const bool not_given){
      // GLV::ofs_running << "Warning: charge mixing parameter is not given and set default value 0.05" << std::endl;
      this->charge_mix_beta = 0.4;
    }

    //max_dmft_step
    try {
      std::vector<std::string> str_val;
      this->read_parameter("max_dmft_step", str_val);

      this->DMFT_step_max = atoi(str_val[0].c_str());
    }
    catch (const std::string messg) {
      GLV::ofs_error << messg << std::endl;
      std::exit(EXIT_FAILURE);
    }
    catch(const bool not_given){
      // GLV::ofs_running << "Warning: max_dmft_step is not given and set default value 1" << std::endl;
      this->DMFT_step_max = 1;
    }
    GLV::ofs_running << "max_dmft_step  " << this->DMFT_step_max << std::endl;

    //max_dft_step
    try {
      std::vector<std::string> str_val;
      this->read_parameter("max_dft_step", str_val);

      this->DFT_step_max = atoi(str_val[0].c_str());
    }
    catch (const std::string messg) {
      GLV::ofs_error << messg << std::endl;
      std::exit(EXIT_FAILURE);
    }
    catch(const bool not_given){
      // GLV::ofs_running << "Warning: max_dmft_step is not given and set default value 1" << std::endl;
      this->DFT_step_max = 1;
    }
    GLV::ofs_running << "max_dft_step  " << this->DFT_step_max << std::endl;

    //mc_step
    try {
      std::vector<std::string> str_val;
      this->read_parameter("mc_step", str_val);

      this->MC_step = atoi(str_val[0].c_str());
    }
    catch (const std::string messg) {
      GLV::ofs_error << messg << std::endl;
      std::exit(EXIT_FAILURE);
    }
    catch(const bool not_given){
      // GLV::ofs_running << "Warning: mc_step is not given and set default value 5000000" << std::endl;
      this->MC_step = 5000000;
    }
    GLV::ofs_running << "mc_step  " << this->MC_step << std::endl;

    //energy_window
    try {
      std::vector<std::string> str_val;
      this->read_parameter("energy_window", str_val, 2);
      this->E_window.resize(2);
      this->E_window[0] = atof(str_val[0].c_str());
      this->E_window[1] = atof(str_val[1].c_str());
    }
    catch (const std::string messg) {
      GLV::ofs_error << messg << std::endl;
      std::exit(EXIT_FAILURE);
    }
    catch(const bool not_given){
      // GLV::ofs_running << "Warning: energy_window is not given and set default value -5.0~5.0" << std::endl;
      this->E_window[0] = -5.0;
      this->E_window[1] = 5.0;
    }
    GLV::ofs_running << "energy_window  " 
            << std::setprecision(3) << this->E_window[0] 
            << "  "
            << std::setprecision(3) << this->E_window[1] << std::endl;

    //restart
    try {
      std::vector<std::string> str_val;
      this->read_parameter("restart", str_val, 4);
      
      this->start_charge_step = atoi(str_val[0].c_str());
      if(this->start_charge_step<1){
        GLV::ofs_error << "The starting charge step must be greater than 0." << std::endl;
        std::exit(EXIT_FAILURE);
      }

      this->start_DMFT_step = atoi(str_val[1].c_str());
      if(this->start_DMFT_step<1){
        GLV::ofs_error << "The starting DMFT step must be greater than 0." << std::endl;
        std::exit(EXIT_FAILURE);
      }

      this->last_charge_step = atoi(str_val[2].c_str());
      if(this->last_charge_step<1){
        GLV::ofs_error << "The last charge step must be greater than 0." << std::endl;
        std::exit(EXIT_FAILURE);
      }

      this->last_DMFT_step = atoi(str_val[3].c_str());
      if(this->last_DMFT_step<1){
        GLV::ofs_error << "The last DMFT step must be greater than 0." << std::endl;
        std::exit(EXIT_FAILURE);
      }

      this->flag_restart = true;
    }
    catch (const std::string messg) {
      GLV::ofs_error << messg << std::endl;
      std::exit(EXIT_FAILURE);
    }
    catch(const bool not_given){
      this->start_charge_step = 1;
      this->start_DMFT_step = 1;
      this->last_charge_step = 1;
      this->last_DMFT_step = 1;
      this->flag_restart = false;
    }

    GLV::ofs_running << "restart  " 
            << this->start_charge_step << " "
            << this->start_DMFT_step << " "
            << this->last_charge_step << " "
            << this->last_DMFT_step << std::endl;

    //dft_xc
    try {
      std::vector<std::string> str_val;
      this->read_parameter("dft_xc", str_val);

      if(std::strcmp("hybrid", str_val[0].c_str())==0)
        this->hyb_func = true;
      else
        this->hyb_func = false;
      
    }
    catch (const std::string messg) {
      GLV::ofs_error << messg << std::endl;
      std::exit(EXIT_FAILURE);
    }
    catch(const bool not_given){
      this->hyb_func = false;
    }

    //hyf_xc_alpha
    try{
      std::vector<std::string> str_val;
      this->read_parameter("hyf_xc_alpha", str_val);

      this->hyf_xc_alpha = atof(str_val[0].c_str());
    }
    catch (const std::string messg) {
      GLV::ofs_error << messg << std::endl;
      std::exit(EXIT_FAILURE);
    }
    catch (const bool not_given){
      this->hyf_xc_alpha = 0.25;
    }

    int local_symmetry_tmp;
    try {
      std::vector<std::string> str_val;
      this->read_parameter("local_symmetry", str_val);
      local_symmetry_tmp = atoi(str_val[0].c_str());
    }
    catch (const std::string messg) {
      GLV::ofs_error << messg << std::endl;
      std::exit(EXIT_FAILURE);
    }
    catch(const bool not_given){
      // GLV::ofs_running << "Warning: local_symmetry is not given and set default value 0" << std::endl;
      local_symmetry_tmp = 0;
    }
    GLV::ofs_running << "local_symmetry  " << local_symmetry_tmp << std::endl;

    GLV::ofs_running << "===============================================================\n" << std::endl;

    return true;
  }
  
  bool input_info::read_parameter(char* keyword, std::vector<std::string>& val, int count)
  {
    debug::codestamp("input_info::read_parameter");

    bool not_given=true;
    std::string file;
    char word[100], word_low[100];

    std::ifstream ifs("DMFT.in", std::ios::in);

    if(!ifs)
    {
      GLV::ofs_error << "Fail to oepn file DMFT.in" << std::endl;
      std::exit(EXIT_FAILURE);
    }

    ifs.seekg(0);   //set the position at the beginning of the file

    while(ifs.good())
    {
      ifs.getline(word, 100);   
      if(ifs.eof()) break;
      
      DMFT::input_info::strtolower(word, word_low);

      std::string line(word_low);

      if(!line.empty())
      {
        line.erase(0, line.find_first_not_of(" "));
        line.erase(line.find_last_not_of(" ") + 1);
      }

      if(!line.empty() && line[0]!='#')
      {
        
        if(line.find_first_of('#') != std::string::npos) line.erase(line.find_first_of('#'));

        std::vector<std::string> key_val;
        while(!line.empty())
        {
          key_val.push_back(line.substr(0,line.find_first_of(' ')));
          line.erase(0,line.find_first_of(' '));
          line.erase(0,line.find_first_not_of(' '));
        }

        if(std::strcmp(key_val[0].c_str(), keyword)==0)
        {
          if(key_val.size() != count+1){
            std::string messg(word);
            throw messg;
          }

          for(int i=0; i<count; i++)
            val.push_back(key_val[i+1]);
          
          not_given=false;
          break;
        }
        else if(
          std::strcmp("calculation", key_val[0].c_str())==0 || 
          std::strcmp("dft_solver", key_val[0].c_str())==0 || 
          std::strcmp("temperature", key_val[0].c_str())==0 ||
          std::strcmp("magnetism", key_val[0].c_str())==0 ||
          std::strcmp("impurity_solver", key_val[0].c_str())==0 ||
          std::strcmp("double_counting", key_val[0].c_str())==0 ||
          std::strcmp("max_charge_step", key_val[0].c_str())==0 ||
          std::strcmp("max_dmft_step", key_val[0].c_str())==0 ||
          std::strcmp("max_dft_step", key_val[0].c_str())==0 ||
          std::strcmp("mc_step", key_val[0].c_str())==0 || 
          std::strcmp("energy_window", key_val[0].c_str())==0 ||
          std::strcmp("local_symmetry", key_val[0].c_str())==0 ||
          std::strcmp("restart", key_val[0].c_str())==0 ||
          std::strcmp("dft_xc", key_val[0].c_str())==0 ||
          std::strcmp("charge_mix_beta", key_val[0].c_str())==0 ||
          std::strcmp("delta_sigma", key_val[0].c_str())==0 ||
          std::strcmp("hybrid_xc_coeff", key_val[0].c_str())==0 )
        {;}
        else
        {
          std::string messg = "Unsupported parameter " + key_val[0];
          throw messg;
        }
      }

      if(ifs.eof()) break;
    }
    ifs.close();

    if(not_given) throw not_given;

    return true;
  }

  void input_info::strtolower(char *in_str, char *out_str)
  {
    char c;
    int len = std::strlen(in_str);
    for (int i = 0; i < len; i++)
    {
      c = in_str[i];
      out_str[i] = std::tolower(c);
    }
    out_str[len] = '\0';
  }

  const void* input_info::parameter(char* keyword)
  {
    char word[40];
    this->strtolower(keyword, word);

    if(std::strcmp("dft_solver", word)==0) return &flag_DFT_solver;
    else if(std::strcmp("calculation", word)==0) return &calculation_type;
    else if(std::strcmp("temperature", word)==0) return &temperature;
    else if(std::strcmp("beta", word)==0) return &beta;
    else if(std::strcmp("magnetism", word)==0) return &flag_magnetism;
    else if(std::strcmp("n_tau", word)==0) return &n_tau;
    else if(std::strcmp("n_omega", word)==0) return &n_omega;
    else if(std::strcmp("impurity_solver", word)==0) return &flag_impurity_solver;
    else if(std::strcmp("double_counting", word)==0) return &flag_double_counting;
    else if(std::strcmp("max_charge_step", word)==0) return &charge_step_max;
    else if(std::strcmp("max_dmft_step", word)==0) return &DMFT_step_max;
    else if(std::strcmp("max_dft_step", word)==0) return &DFT_step_max;
    else if(std::strcmp("mc_step", word)==0) return &MC_step;
    else if(std::strcmp("hyb_func", word)==0) return &hyb_func;
    else if(std::strcmp("hyf_xc_alpha", word)==0) return &hyf_xc_alpha;
    else if(std::strcmp("charge_mix_beta", word)==0) return &charge_mix_beta;
    else if(std::strcmp("delta_sigma", word)==0) return &delta_sigma;
    else if(std::strcmp("start_charge_step", word)==0) return &start_charge_step;
    else if(std::strcmp("start_dmft_step", word)==0) return &start_DMFT_step;
    else if(std::strcmp("last_charge_step", word)==0) return &last_charge_step;
    else if(std::strcmp("last_dmft_step", word)==0) return &last_DMFT_step;
    else if(std::strcmp("restart", word)==0) return &flag_restart;
    else
    {
      GLV::ofs_error << "No parameter " << keyword << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }

  void input_info::test()
  {
    if(this->flag_magnetism==-1)
    {
      GLV::ofs_error << "Please set the parameter of magnetism" << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }

  void input_info::out()
  {
    debug::codestamp("input_info::out");

    std::string file;
    file="DMFT.in";
    std::ofstream ofs(file.c_str(),std::ios::out);
    if(!ofs)
    {
      GLV::ofs_error << "Fail to oepn " << file.c_str() << std::endl;
      std::exit(EXIT_FAILURE);
    }

    ofs << "DFT_solver" << std::setw(5) << this->flag_DFT_solver << '\n';
    ofs << "temperature" << std::setw(10) << std::setprecision(2) << this->temperature << '\n';
    ofs << "magnetism" << std::setw(5) << this->flag_magnetism << '\n';
    ofs << "impurity_solver" << std::setw(5) << this->flag_impurity_solver << '\n';
    ofs << "n_tau" << std::setw(5) << this->n_tau << '\n';
    ofs << "n_omega" << std::setw(5) << this->n_omega << '\n';
    ofs << "double_counting" << std::setw(5) << this->flag_double_counting << std::endl;

    ofs.close();
  }

}
