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
  for(const std::string& arg : all_args)
  {
    size_t pos_begin=0;
    if (arg.substr(0,2)=="--") pos_begin=2;
    else if(arg.substr(0,1)=="-") pos_begin=1;
    else
    {
      std::cout << "Error argument : " << arg << '\n';
      std::exit(EXIT_FAILURE);
    }

    size_t pos_seg = arg.find('=');
    if(pos_seg == std::string::npos)
    {
      std::cout << "Error argument : " << arg << '\n';
      std::exit(EXIT_FAILURE);
    }

    if(std::strcmp("dmft.step", arg.substr(pos_begin, pos_seg-pos_begin).c_str())==0)
      args.global_step = atoi(arg.substr(pos_seg+1).c_str());
    else if(std::strcmp("eva.sigma_only", arg.substr(pos_begin, pos_seg-pos_begin).c_str())==0)
      args.sigma_only = (bool)atoi(arg.substr(pos_seg+1).c_str());
    else if(std::strcmp("eva.spectrum", arg.substr(pos_begin, pos_seg-pos_begin).c_str())==0)
      args.cal_spectrum = (bool)atoi(arg.substr(pos_seg+1).c_str());
    else
    {
      std::cout << "Illegal arguments : " << arg << '\n';
      std::exit(EXIT_FAILURE);
    }
  }
  return true;
}