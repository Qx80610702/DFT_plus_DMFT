#pragma once
#include "input.h" 

#include <string>
#include <vector>

namespace DFT_output
{
  class atoms_info
  {
    public:
    atoms_info(){;}
    ~atoms_info(){;}

    void read(DMFT::input_info& in);
    bool read_structure_symmetry();
    bool read_local_symmetry(DMFT::input_info& in);
    bool read_correlated_atoms_info();

    void out();       //test whether reading worked corrected

    //=================
    //  interfaces
    //=================
    int total_atoms(){return n_DMFT_atoms;}
    int inequ_atoms(){return n_inequivalent_atoms;}
    int norb(){return n_DMFT_orb;}
    int equ_atom(int iatom);
    int L(const int iatom);
    int iorb2ibasis(const int iorb);
    double Uval(const int iatom);
    double Jval(const int iatom);
    double occ_num(const int iatom, const int ispin);
    std::vector<std::vector<std::vector<double>>>& occ_num_m(){return occ_number_m;}
    std::vector<std::vector<double>>& occ_num_ref(){return occ_number;}
    int ineq_iatom(const int ineq);
    int magnetic(const int iatom);
    const std::vector<int>& iatom_norb(){return sub_norb;}
    int local_sym(){return local_symmetry;}
    const std::vector<std::vector<int>>& Im2iorb(){return DMFT_orb_index;}

    public:
    //For matrix form quantity in the form as A[m*m]
    template <typename T>
    static inline void symmetry_operation_matrix(
    const int symm, const int corr_L, const int m_tot, T* matrix);

    //For vector form quantity in the form as A[m]
    template <typename T>
    static inline void symmetry_operation_vector(
    const int symm, const int corr_L, const  int m_tot, T* colum);

    private:
    int n_DMFT_atoms;            //Total number of correlated atoms requring DMFT correction
    int n_inequivalent_atoms;    //Total number of inquivalent correlated atoms requring DMFT correction
    std::vector<int> equivalent_atom;       //Maping between the correlated atoms and their equivalent atoms; equivalent_atom[n_DMFT_atoms]
    std::vector<int> ineq_atom_to_iatom;    //maping between the i-th inequivalent atom and its iatom index

    std::vector<int> angular_momment;       //angular_momment[iatom]
    std::vector<std::vector<int>> magnetic_number;      //magnetic_number[iatom][m]
    std::vector<std::vector<int>> DMFT_orb_index;       //Mapping between correlated atoms I, orbit m and orb index;DMFT_orb_index[iatom][m]
    int n_DMFT_orb;
    std::vector<int> iorb_ibasis;           //Mapping between iorb and ibasis

    std::vector<std::vector<double>> occ_number;        //occ_number[iatom][spin up(0) or down(1)]
    std::vector<std::vector<std::vector<double>>> occ_number_m; //occ_number_m[iatom][is][m]
    std::vector<double> Hubbard_U;          //Parameter Hubbard U; unit:Hartree
    std::vector<double> Hund_J;             //Parameter Hund J; unit:Hartree
    std::vector<int> mag_moment;            //the direction of magnetic moment

    // local orbital symmetry; 
    // 0:no symmetry and no rotation
    // 1:cubic symmetry
    // 2:t2g only
    // 3:eg only
    // ...
    int local_symmetry;                     

    //=======sub-shell=========
    std::vector<int> sub_norb;             //number of orbitals in the sub-shell; sub_norb[iatom]

  };

  template <typename T>
  void atoms_info::symmetry_operation_vector(
            const int symm, const int corr_L, 
            const  int m_tot, T* colum)
  {
    if(symm==0){;}
    else if(symm==1 && corr_L==2)//cubic symmetry, d orbital
    {
      //Real harmonics for l=2
      //m=-2(dxy); m=-1(dyz); m=0(dz^2); m=1(dxz); m=2(dx^2-y^2)
      T t2g = (colum[0]+colum[1]+colum[3])/3.0;
      T eg = (colum[2]+colum[4])/2.0;

      for(int m=0; m<m_tot; m++)
        if(m==2 || m==4) colum[m]=eg;
        else colum[m]=t2g;
    }
    else if(symm==2)//t2g only
    {
      T t2g = (colum[0]+colum[1]+colum[2])/3.0;
      for(int m=0; m<m_tot; m++) colum[m]=t2g;
    }
    else if(symm==3)//eg only
    {
      T eg = (colum[0]+colum[1])/2.0;
      for(int m=0; m<m_tot; m++) colum[m]=eg;
    }
    return;
  }

  template <typename T>
  void atoms_info::symmetry_operation_matrix(
            const int symm, const int corr_L, 
            const  int m_tot, T* matrix)
  {
    if(symm==0){;}
    else if(symm==1 && corr_L==2)//cubic symmetry, d orbital
    {
      T t2g = (matrix[0]+matrix[1*m_tot+1]+matrix[3*m_tot+3])/3.0;
      T eg = (matrix[2*m_tot+2]+matrix[4*m_tot+4])/2.0;

      for(int m=0; m<m_tot; m++)
        if(m==2 || m==4) matrix[m*m_tot+m]=eg;
        else matrix[m*m_tot+m]=t2g;
    }
    else if(symm==2)//t2g only
    {
      T t2g = (matrix[0]+matrix[1*m_tot+1]+matrix[2*m_tot+2])/3.0;
      for(int m=0; m<m_tot; m++) matrix[m*m_tot+m]=t2g;
    }
    else if(symm==3)//eg only
    {
      T eg = (matrix[0*m_tot+0]+matrix[1*m_tot+1])/2.0;
      for(int m=0; m<m_tot; m++) matrix[m*m_tot+m]=eg;
    }
    return;
  }

}
