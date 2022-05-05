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

extern "C" void aims_(int* mpi_comm, int* unit, bool* mpi_switch);

int main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);

  //Standard output redirection
  // std::streambuf* coutBuf = std::cout.rdbuf();
  // std::ofstream stdout("DMFT_running.log", std::ios_base::app);
  // std::streambuf* fileBuf = stdout.rdbuf();
  // std::cout.rdbuf(fileBuf);

  // Parsing commond line
  // argument_lists args_val = {1, 1, 1, 1, false, false, false};  //default values
  // if(argc>1)
  // {
  //   std::vector<std::string> all_args(argv+1,argv+argc);
  //   bool success=parsing_commond_line(all_args, args_val);
  //   if(!success)
  //   {
  //     std::cout << "Error in parsing commond line" << std::endl;
  //     std::exit(EXIT_FAILURE);
  //   }
  // }

  // DFT_DMFT_solver psolver(args_val);
  DFT_plus_DMFT::solver DFT_DMFT_solver;
  DFT_DMFT_solver.solve();

  // std::cout.rdbuf(coutBuf);
  
  // int unit=6;
  // bool use_mpi=true;
  // int mpi_comm_global=MPI_COMM_WORLD;
  // aims_(&mpi_comm_global,&unit,&use_mpi);

  MPI_Finalize();

  return 0;
}

bool parsing_commond_line(std::vector<std::string>& all_args, argument_lists& args)
{
  for(std::vector<std::string>::iterator iter=all_args.begin(); iter!=all_args.end(); )
  {
    std::string param = *iter;
    
    size_t pos_begin=0;
    if (param.substr(0,2)=="--") pos_begin=2;
    else if(param.substr(0,1)=="-") pos_begin=1;
    else
    {
      std::cout << "Error argument : " << param << '\n';
      std::exit(EXIT_FAILURE);
    }

    if(std::strcmp("current_step", param.substr(pos_begin).c_str())==0){
      std::string value1 = *(iter+1);
      std::string value2 = *(iter+2);

      args.current_charge_step = atoi(value1.c_str());
      args.current_DMFT_step = atoi(value2.c_str());

      iter += 3;
    }
    else if(std::strcmp("last_step", param.substr(pos_begin).c_str())==0){
      std::string value1 = *(iter+1);
      std::string value2 = *(iter+2);

      args.last_charge_step = atoi(value1.c_str());
      args.last_DMFT_step = atoi(value2.c_str());

      iter += 3;
    }
    else if(std::strcmp("eva.sigma_only", param.substr(pos_begin).c_str())==0){
      std::string value = *(iter+1);
      args.sigma_only = (bool)atoi(value.c_str());
      iter += 2;
    }
    else if(std::strcmp("eva.spectrum", param.substr(pos_begin).c_str())==0){
      std::string value = *(iter+1);
      args.cal_spectrum = (bool)atoi(value.c_str());
      iter += 2;
    }
    else if(std::strcmp("eva.density", param.substr(pos_begin).c_str())==0){
      std::string value = *(iter+1);
      args.update_density = (bool)atoi(value.c_str());
      iter += 2;
    }
    else
    {
      std::cout << "Illegal arguments : " << param << '\n';
      std::exit(EXIT_FAILURE);
    }
  }

  return true;
}

