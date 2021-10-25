#include "mpi_environment.h"
#include "./solver/solver.h"
#include "timer.h"

#include <mpi.h>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <string>

bool parsing_commond_line(std::vector<std::string>& all_args, argument_lists& args);

int main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);

  //Standard output redirection
  std::streambuf* coutBuf = std::cout.rdbuf();
  std::ofstream stdout("DMFT_running.log", std::ios_base::app);
  std::streambuf* fileBuf = stdout.rdbuf();
  std::cout.rdbuf(fileBuf);

  // Parsing commond line
  argument_lists args_val ={1,false,false};
  if(argc>1)
  {
    std::vector<std::string> all_args(argv+1,argv+argc);
    bool success=parsing_commond_line(all_args, args_val);
    if(!success)
    {
      std::cout << "Error in parsing commond line" << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }

  DFT_DMFT_solver psolver(args_val);
  psolver.solve();

  std::cout.rdbuf(coutBuf);

  MPI_Finalize();

  return 0;
}

bool parsing_commond_line(std::vector<std::string>& all_args, argument_lists& args)
{
  for(std::vector<std::string>::iterator iter=all_args.begin(); iter!=all_args.end(); iter += 2)
  {
    std::string param = *iter;
    std::string value = *(iter+1);
    
    size_t pos_begin=0;
    if (param.substr(0,2)=="--") pos_begin=2;
    else if(param.substr(0,1)=="-") pos_begin=1;
    else
    {
      std::cout << "Error argument : " << param << '\n';
      std::exit(EXIT_FAILURE);
    }

    if(std::strcmp("dmft.step", param.substr(pos_begin).c_str())==0)
      args.global_step = atoi(value.c_str());
    else if(std::strcmp("eva.sigma_only", param.substr(pos_begin).c_str())==0)
      args.sigma_only = (bool)atoi(value.c_str());
    else if(std::strcmp("eva.spectrum", param.substr(pos_begin).c_str())==0)
      args.cal_spectrum = (bool)atoi(value.c_str());
    else
    {
      std::cout << "Illegal arguments : " << param << '\n';
      std::exit(EXIT_FAILURE);
    }
  }
  return true;
}