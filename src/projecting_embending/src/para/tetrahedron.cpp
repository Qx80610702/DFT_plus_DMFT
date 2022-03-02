#include "tetrahedron.h"
#include "../debug.h"
#include "../mpi_environment.h"

#include <fstream>
#include <iostream>
#include <cstdlib>   //Use exit function
#include <iomanip>   //use setw and setprecison function

namespace DFT_output
{
  bool tretrahedron::read()
  {
    debug::codestamp("tretrahedron::read");
    // if(mpi_rank()==0) std::cout << "Reading tetrahedron.dat ......" << std::endl;

    int i_x, i_y, i_z;

    std::ifstream ifs("dft/outputs/tetrahedron.dat", std::ios::in);

    if (!ifs)  
	  {
	  	std::cout << "Fail to oepn dft/outputs/tetrahedron.dat" << std::endl;
      std::exit(EXIT_FAILURE);
    }

    this->lattice_vector.resize(3);
    for(int i_vec=0; i_vec<3; i_vec++)
    {
      this->lattice_vector[i_vec].resize(3);
    }

    this->k_points_nosymmetry.resize(3);

    ifs.seekg(0);      //set the position at the beginning of the file
    for(int i_vec=0; i_vec<3; i_vec++)
    {
      ifs >> this->lattice_vector[i_vec][0];
      ifs >> this->lattice_vector[i_vec][1];
      ifs >> this->lattice_vector[i_vec][2];
      ifs.ignore(150,'\n');
    }

    ifs >> this->k_points_nosymmetry[0];
    ifs >> this->k_points_nosymmetry[1];
    ifs >> this->k_points_nosymmetry[2];
    ifs.ignore(150,'\n');

    this->ik_irred_map.resize(this->k_points_nosymmetry[0]);
    this->k_nosym_weight.resize(this->k_points_nosymmetry[0]);
    for(i_x=0; i_x<this->k_points_nosymmetry[0]; i_x++)
    {
      this->ik_irred_map[i_x].resize(this->k_points_nosymmetry[1]);
      this->k_nosym_weight[i_x].resize(this->k_points_nosymmetry[1]);
      for(i_y=0; i_y<this->k_points_nosymmetry[1]; i_y++)
      {
        this->ik_irred_map[i_x][i_y].resize(this->k_points_nosymmetry[2]);
        this->k_nosym_weight[i_x][i_y].resize(this->k_points_nosymmetry[2]);
      }
    }

    while(ifs.good())
    {
      ifs >> i_x;
      if(ifs.eof()) break;
      ifs >> i_y;
      ifs >> i_z;
      ifs >> this->ik_irred_map[i_x][i_y][i_z];
      ifs >> this->k_nosym_weight[i_x][i_y][i_z];

      ifs.ignore(150, '\n');

      if(ifs.eof()) break;  //Check whether end of file is reached 
    }
    ifs.close();

    return true;
  }

  void tretrahedron::out()
  {
    debug::codestamp("tretrahedron::out");

    std::string file;

    //===============================
    //  tetrahedron.dat
    //===============================
    file="tetrahedron.dat";
    std::ofstream ofs(file.c_str(),std::ios::out);
    if(!ofs)
    {
      std::cout << "Fail to oepn " << file.c_str() << std::endl;
      std::exit(EXIT_FAILURE);
    }

    for(int i_vec=0; i_vec<3; i_vec++)
    {
      ofs << std::setw(12) << std::fixed << std::setprecision(6) << this->lattice_vector[i_vec][0]
          << std::setw(12) << std::fixed << std::setprecision(6) << this->lattice_vector[i_vec][1]
          << std::setw(12) << std::fixed << std::setprecision(6) << this->lattice_vector[i_vec][2] << '\n';
    }

    ofs << std::setw(4) << this->k_points_nosymmetry[0]
        << std::setw(4) << this->k_points_nosymmetry[1]
        << std::setw(4) << this->k_points_nosymmetry[2] << '\n';

    for(int i_x=0; i_x<this->k_points_nosymmetry[0]; i_x++)
    {
      for(int i_y=0; i_y<this->k_points_nosymmetry[1]; i_y++)
      {
        for(int i_z=0; i_z<this->k_points_nosymmetry[2]; i_z++)
        {
          ofs << std::setw(4) << i_x << std::setw(4) << i_y << std::setw(4) << i_z 
          << std::setw(6)  << this->ik_irred_map[i_x][i_y][i_z] << std::endl; 
        }
      }
    }

    ofs.close();
  }
  
}
