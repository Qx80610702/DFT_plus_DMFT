#pragma once

#include <string>
#include <vector>

namespace DFT_output
{
  class tretrahedron
  {
    public:
    tretrahedron(){;}
    ~tretrahedron(){;}

    // void initialize();
    bool read();

    void out();       //test whether reading worked corrected

    private:
    //lattice_vector[i_vec][i_coord]
    std::vector<std::vector<double>> lattice_vector; 

    //Number of k-points in three direction before symmetrilization
    std::vector<int> k_points_nosymmetry;         
    
    //map between the k-points and their corresponding irreduced k-points
    // ik_irred_map[i_x_k][i_y_k][i_z_k]
    std::vector<std::vector<std::vector<int>>> ik_irred_map;

    std::vector<std::vector<std::vector<double>>> k_nosym_weight;

  };
  
}
