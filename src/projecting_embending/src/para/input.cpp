#include "input.h"
#include "../constants.h"
#include "../debug/debug.h"
#include "../mpi_environment.h"

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

    if(mpi_rank()==0) std::cout << "Reading DMFT.in ......" << std::endl;

    //dft_solver
    try {
      std::vector<std::string> str_val;
      this->read_parameter("dft_solver", str_val);

      if(std::strcmp("aims",str_val[0].c_str())==0) this->flag_DFT_solver = 1;
      else if(std::strcmp("abacus", str_val[0].c_str())==0) this->flag_DFT_solver = 2;
      else
      {
        std::cout << "unsupported value of DFT_solver" << std::endl;
        std::exit(EXIT_FAILURE);
      }
    }
    catch (const std::string messg) {
      std::cout << messg << std::endl;
      std::exit(EXIT_FAILURE);
    }
    catch(const bool not_given){
      if(mpi_rank()==0) std::cout << "Warning: dft_solver is not given and set default value aims" << std::endl;
      this->flag_DFT_solver = 1;   //aims
    }

    //temperature
    try {
      std::vector<std::string> str_val;
      this->read_parameter("temperature", str_val);

      this->temperature = atof(str_val[0].c_str());
    }
    catch (const std::string messg) {
      std::cout << messg << std::endl;
      std::exit(EXIT_FAILURE);
    }
    catch(const bool not_given){
      if(mpi_rank()==0) std::cout << "Warning: temperature is not given and set default value 300 Kelven" << std::endl;
      this->temperature = 300.0;
    }
    this->beta = Hartree_to_eV/(this->temperature*K_BOLTZMAN_EV);
    // this->n_tau = 25/(this->temperature*K_BOLTZMAN_EV) > 500 ? 25/(this->temperature*K_BOLTZMAN_EV) : 500;
    this->n_tau = 25/(this->temperature*K_BOLTZMAN_EV);
    this->n_omega = this->n_tau > 500 ? 500 : this->n_tau;
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
        std::cout << "unsupported value of magnetism" << std::endl;
        std::exit(EXIT_FAILURE);
      }
    }
    catch (const std::string messg) {
      std::cout << messg << std::endl;
      std::exit(EXIT_FAILURE);
    }
    catch(const bool not_given){
      std::cout << "Error: parameter magnetism must be given!!!" << std::endl;
      std::exit(EXIT_FAILURE);
    }

    //impurity_solver
    try {
      std::vector<std::string> str_val;
      this->read_parameter("impurity_solver", str_val);

      if(std::strcmp("alps_cthyb",str_val[0].c_str())==0) this->flag_impurity_solver = 1;
      else if(std::strcmp("alps_cthyb_segment",str_val[0].c_str())==0) this->flag_impurity_solver = 2;
      else if(std::strcmp("lg_cthyb",str_val[0].c_str())==0) this->flag_impurity_solver = 3;
      else if(std::strcmp("rutgers_cthyb",str_val[0].c_str())==0) this->flag_impurity_solver = 4;
      else
      {
        std::cout << "unsupported value of impurity_solver" << std::endl;
        std::exit(EXIT_FAILURE);
      }
    }
    catch (const std::string messg) {
      std::cout << messg << std::endl;
      std::exit(EXIT_FAILURE);
    }
    catch(const bool not_given){
      if(mpi_rank()==0) std::cout << "Warning: impurity_solver is not given and set default value rutgers_cthyb" << std::endl;
        this->flag_impurity_solver = 3;
    }

    //double_counting
    try {
      std::vector<std::string> str_val;
      this->read_parameter("double_counting", str_val);

      if(std::strcmp("nominal",str_val[0].c_str())==0) this->flag_double_counting = 1;
      else
      {
        std::cout << "unsupported value of double_counting" << std::endl;
        std::exit(EXIT_FAILURE);
      }

    }
    catch (const std::string messg) {
      std::cout << messg << std::endl;
      std::exit(EXIT_FAILURE);
    }
    catch(const bool not_given){
      if(mpi_rank()==0) std::cout << "Warning: double_counting is not given and set default value nominal" << std::endl;
      this->flag_double_counting = 1;
    }

    //max_dmft_step
    try {
      std::vector<std::string> str_val;
      this->read_parameter("max_dmft_step", str_val);

      this->DMFT_step = atoi(str_val[0].c_str());
    }
    catch (const std::string messg) {
      std::cout << messg << std::endl;
      std::exit(EXIT_FAILURE);
    }
    catch(const bool not_given){
      if(mpi_rank()==0) std::cout << "Warning: max_dmft_step is not given and set default value 10" << std::endl;
      this->DMFT_step = 10;
    }

    //mc_step
    try {
      std::vector<std::string> str_val;
      this->read_parameter("mc_step", str_val);

      this->MC_step = atoi(str_val[0].c_str());
    }
    catch (const std::string messg) {
      std::cout << messg << std::endl;
      std::exit(EXIT_FAILURE);
    }
    catch(const bool not_given){
      if(mpi_rank()==0) std::cout << "Warning: mc_step is not given and set default value 5000000" << std::endl;
      this->MC_step = 5000000;
    }

    //energy_window
    try {
      std::vector<std::string> str_val;
      this->read_parameter("energy_window", str_val, 2);
      this->E_window.resize(2);
      this->E_window[0] = atof(str_val[0].c_str());
      this->E_window[1] = atof(str_val[1].c_str());
    }
    catch (const std::string messg) {
      std::cout << messg << std::endl;
      std::exit(EXIT_FAILURE);
    }
    catch(const bool not_given){
      if(mpi_rank()==0) std::cout << "Warning: energy_window is not given and set default value -5.0~5.0" << std::endl;
      this->E_window[0] = -5.0;
      this->E_window[1] = 5.0;
    }

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
      std::cout << "Fail to oepn file DMFT.in" << std::endl;
      std::exit(EXIT_FAILURE);
    }

    ifs.seekg(0);   //set the position at the beginning of the file

    while(ifs.good())
    {
      ifs.getline(word, 100);   
      if(ifs.eof()) break;
      
      this->strtolower(word, word_low);

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
          if(key_val.size() != count+1) throw "Error in number of values of parameter  " + key_val[0];

          for(int i=0; i<count; i++)
            val.push_back(key_val[i+1]);
          
          not_given=false;
          break;
        }
        else if(
          std::strcmp("dft_solver", key_val[0].c_str())==0 || 
          std::strcmp("temperature", key_val[0].c_str())==0 ||
          std::strcmp("magnetism", key_val[0].c_str())==0 ||
          std::strcmp("impurity_solver", key_val[0].c_str())==0 ||
          std::strcmp("double_counting", key_val[0].c_str())==0 ||
          std::strcmp("max_dmft_step", key_val[0].c_str())==0 ||
          std::strcmp("mc_step", key_val[0].c_str())==0 || 
          std::strcmp("energy_window", key_val[0].c_str())==0 ||
          std::strcmp("local_symmetry", key_val[0].c_str())==0 )
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
    else if(std::strcmp("temperature", word)==0) return &temperature;
    else if(std::strcmp("beta", word)==0) return &beta;
    else if(std::strcmp("magnetism", word)==0) return &flag_magnetism;
    else if(std::strcmp("n_tau", word)==0) return &n_tau;
    else if(std::strcmp("n_omega", word)==0) return &n_omega;
    else if(std::strcmp("impurity_solver", word)==0) return &flag_impurity_solver;
    else if(std::strcmp("double_counting", word)==0) return &flag_double_counting;
    else if(std::strcmp("max_dmft_step", word)==0) return &DMFT_step;
    else if(std::strcmp("mc_step", word)==0) return &MC_step;
    else
    {
      std::cout << "No parameter " << keyword << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }

  void input_info::test()
  {
    if(this->flag_magnetism==-1)
    {
      std::cout << "Please set the parameter of magnetism" << std::endl;
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
      std::cout << "Fail to oepn " << file.c_str() << std::endl;
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
