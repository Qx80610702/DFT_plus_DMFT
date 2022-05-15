#include "correlated_atoms.h"
#include "../debug.h"
#include "../constants.h"
#include "input.h"
#include "../mpi_environment.h"
#include "../global_variables.h"

#include <memory>
#include <fstream>
#include <iostream>
#include <cstring>
#include <cstdlib>   
#include <iomanip>   //use setw and setprecison function

namespace DFT_output
{ 
  int atoms_info::equ_atom(int iatom)
  {
    if(iatom<0 || iatom>=this->n_DMFT_atoms)
    {
      GLV::ofs_error << "Array equivalent_atom is out of range" << std::endl;
      std::exit(EXIT_FAILURE);
    }

    return this->equivalent_atom[iatom];
  }

  int atoms_info::L(const int iatom)
  {
    if(iatom<0 || iatom>=this->n_DMFT_atoms)
    {
      GLV::ofs_error << "Array angular_momment is out of range" << std::endl;
      std::exit(EXIT_FAILURE);
    }

    return this->angular_momment[iatom];
  }

  int atoms_info::iorb2ibasis(const int iorb)
  {
    if(iorb<0 || iorb>=this->n_DMFT_orb)
    {
      GLV::ofs_error << "Array iorb_ibasis is out of range" <<std::endl;
      std::exit(EXIT_FAILURE);
    }
    
    return this->iorb_ibasis[iorb];
  }

  double atoms_info::Uval(const int iatom)
  {
    if(iatom<0 || iatom>=this->n_DMFT_atoms)
    {
      GLV::ofs_error << "Array Hubbard_U is out of range" <<std::endl;
      std::exit(EXIT_FAILURE);
    }
    return this->Hubbard_U[iatom];
  }

  double atoms_info::Jval(const int iatom)
  {
    if(iatom<0 || iatom>=this->n_DMFT_atoms)
    {
      GLV::ofs_error << "Array Hund_J is out of range" <<std::endl;
      std::exit(EXIT_FAILURE);
    }
    return this->Hund_J[iatom];
  }

  double atoms_info::occ_num(const int iatom, const int ispin)
  {
    if(iatom<0 || iatom>=this->n_DMFT_atoms
      || ispin<0 || ispin>=2)
    {
      GLV::ofs_error << "Array Hund_J is out of range" <<std::endl;
      std::exit(EXIT_FAILURE);
    }
    return this->occ_number[iatom][ispin];
  }

  int atoms_info::ineq_iatom(const int ineq)
  {
    if(ineq<0 || ineq>=this->n_inequivalent_atoms)
    {
      GLV::ofs_error << "Array ineq_atom_to_iatom is out of range" <<std::endl;
      std::exit(EXIT_FAILURE);
    }
    return this->ineq_atom_to_iatom[ineq];
  }

  int atoms_info::magnetic(const int iatom)
  {
    if(iatom<0 || iatom>=this->n_DMFT_atoms)
    {
      GLV::ofs_error << "Array mag_moment is out of range" <<std::endl;
      std::exit(EXIT_FAILURE);
    }
    return this->mag_moment[iatom];
  }

  void atoms_info::read(DMFT::input_info& in)
  {
    GLV::ofs_running << "Reading correlated atoms information ......" << std::endl;

    this->read_structure_symmetry();
    this->read_local_symmetry(in);
    this->read_correlated_atoms_info();
    return;
  }

  bool atoms_info::read_structure_symmetry()
  {
    debug::codestamp("atoms_info::read_structure_symmetry");

    int i_atom_tmp1, i_atom_tmp2; 
    
    std::ifstream if_symmetry("dft/outputs_to_DMFT/symmetry.dat", std::ios::in);

    if (!if_symmetry) 
	  {
	  	GLV::ofs_error << "Fail to oepn file dft/outputs_to_DMFT/symmetry.dat" << std::endl;
      std::exit(EXIT_FAILURE);
    }

    if_symmetry.seekg(0); //set the position at the beginning of the file

    if_symmetry >> this->n_DMFT_atoms;
    if_symmetry.ignore(150, '\n');

    if_symmetry >> this->n_inequivalent_atoms;
    if_symmetry.ignore(150, '\n');

    this->equivalent_atom.resize(this->n_DMFT_atoms);
    while(if_symmetry.good())
    {
      if_symmetry >> i_atom_tmp1;
      if(if_symmetry.eof()) break;   //Check whether end of file is reached
      if_symmetry >> i_atom_tmp2;
      this->equivalent_atom[i_atom_tmp1] = i_atom_tmp2;
      if_symmetry.ignore(150, '\n');

      if(if_symmetry.eof()) break;   //Check whether end of file is reached 
    }

    if_symmetry.close();

    this->ineq_atom_to_iatom.resize(this->n_inequivalent_atoms);
    int count=-1;
    for(int iatom=0; iatom<this->n_DMFT_atoms; iatom++)
    {
      if(iatom==this->equivalent_atom[iatom])
      {
        count++;
        this->ineq_atom_to_iatom[count] = iatom;
      }
    }

    return true;	  	
  }

  bool atoms_info::read_correlated_atoms_info()
  {
    debug::codestamp("atoms_info::read_correlated_atoms_info");

    int iatom, magnetic_num, ibasis;
    char word[40];
    std::unique_ptr<std::unique_ptr<int[]>[]> tmp;
    
    this->angular_momment.resize(this->n_DMFT_atoms);
    this->magnetic_number.resize(this->n_DMFT_atoms);
    this->DMFT_orb_index.resize(this->n_DMFT_atoms);
    this->Hubbard_U.resize(this->n_DMFT_atoms);
    this->Hund_J.resize(this->n_DMFT_atoms);
    this->mag_moment.resize(this->n_DMFT_atoms);
    this->sub_norb.resize(this->n_DMFT_atoms);

    this->occ_number.resize(this->n_DMFT_atoms);
    this->occ_number_m.resize(this->n_DMFT_atoms);
    for(int iatom=0; iatom<this->n_DMFT_atoms; iatom++)
    {
      this->occ_number[iatom].resize(2);
      this->occ_number_m[iatom].resize(2);
    }

    tmp = std::make_unique<std::unique_ptr<int[]>[]>(this->n_DMFT_atoms);

    // if(mpi_rank()==0) GLV::ofs_error << "Reading correlated_atoms.info ......" << std::endl;

    std::ifstream ifs("dft/outputs_to_DMFT/correlated_atoms.info", std::ios::in);
    if (!ifs) 
	  {
	  	GLV::ofs_error << "Fail to oepn file dft/outputs_to_DMFT/correlated_atoms.info" << std::endl;
      std::exit(EXIT_FAILURE);
    }

    ifs.seekg(0);     //set the position at the beginning of the file

    this->n_DMFT_orb =0;
    int iorb=0;
    while(ifs.good())
    {
      ifs >> word;
      if(ifs.eof()) break;
      ifs >> iatom;
      ifs.ignore(150,'\n');
      ifs >> word;
      ifs >> this->angular_momment[iatom];
    
      ifs.ignore(150,'\n');

      ifs >> this->Hubbard_U[iatom];    
      ifs >> this->Hund_J[iatom];
      this->Hubbard_U[iatom] /= GLC::Hartree_to_eV;   //eV to Hartree
      this->Hund_J[iatom] /= GLC::Hartree_to_eV;      //eV to Hartree

      ifs >> this->mag_moment[iatom];
      
      const int ineq = this->equivalent_atom[iatom];
      if(this->local_symmetry==0
        || this->local_symmetry==1)
      {
        this->sub_norb[iatom] = 2*this->angular_momment[iatom]+1;
      }
      else if(this->local_symmetry==2 
        && this->angular_momment[iatom]==2) //t2g only
      {
        this->sub_norb[iatom] = 3;
      }
      else if(this->local_symmetry==3 
        && this->angular_momment[iatom]==2) //eg only
      {
        this->sub_norb[iatom] = 2;
      }
      else
      {
        GLV::ofs_error << "Unsupported local symmetry!!!" << std::endl;
        std::exit(EXIT_FAILURE); 
      }
      
      this->n_DMFT_orb += this->sub_norb[iatom];
      ifs.ignore(150,'\n');

      this->occ_number_m[iatom][0].resize(this->sub_norb[iatom]); //spin up
      this->occ_number_m[iatom][1].resize(this->sub_norb[iatom]); //spin down

      this->DMFT_orb_index[iatom].resize(this->sub_norb[iatom]);
      this->magnetic_number[iatom].resize(this->sub_norb[iatom]);
      tmp[iatom] = std::make_unique<int[]>(this->sub_norb[iatom]);
      int m_count=0;
      for(int m=0; m<2*this->angular_momment[iatom]+1; m++)
      {
        ifs >> magnetic_num;
        ifs >> ibasis;
        
        if(this->local_symmetry==0
          || this->local_symmetry==1 )
        {
          this->magnetic_number[iatom][m_count] = magnetic_num;
          this->DMFT_orb_index[iatom][m_count] = iorb;
          tmp[iatom][m_count] = ibasis;
          m_count++;
          iorb++;
        }
        else if(this->local_symmetry==2 
          && this->angular_momment[iatom]==2) //t2g only
        {
          if(magnetic_num==-2 || magnetic_num==-1 || magnetic_num==1)
          {
            this->magnetic_number[iatom][m_count] = magnetic_num;
            this->DMFT_orb_index[iatom][m_count] = iorb;
            tmp[iatom][m_count] = ibasis;
            m_count++;
            iorb++;
          }
        }
        else if(this->local_symmetry==3 
          && this->angular_momment[iatom]==2) //eg only
        {
          if(magnetic_num==0 || magnetic_num==2)
          {
            this->magnetic_number[iatom][m_count] = magnetic_num;
            this->DMFT_orb_index[iatom][m_count] = iorb;
            tmp[iatom][m_count] = ibasis;
            m_count++;
            iorb++;
          }
        }
      }
      ifs.ignore(150,'\n');

      if(ifs.eof()) break;
    }
    ifs.close();

    this->iorb_ibasis.resize(this->n_DMFT_orb);
    for(int iatom=0; iatom<this->n_DMFT_atoms; iatom++)
    {
      for(int m=0; m<this->sub_norb[iatom]; m++)
      {
        this->iorb_ibasis[this->DMFT_orb_index[iatom][m]] = tmp[iatom][m];
      }
    }

    return true;	  	
  }

  bool atoms_info::read_local_symmetry(DMFT::input_info& in)
  {
    debug::codestamp("input_info::read_local_symmetry");

    try {
      std::vector<std::string> str_val;
      in.read_parameter("local_symmetry", str_val);
      this->local_symmetry = atoi(str_val[0].c_str());
    }
    catch (const std::string messg) {
      GLV::ofs_error << messg << std::endl;
      std::exit(EXIT_FAILURE);
    }
    catch(const bool not_given){
      // GLV::ofs_running << "Warning: local_symmetry is not given and set default value 0" << std::endl;
      this->local_symmetry = 0;
    }

    return true;
  }

  void atoms_info::out()
  {
    debug::codestamp("atoms_info::out");

    std::string file;

    //================================
    //    symmetry.dat
    //================================
    file="symmetry.dat";
    std::ofstream ofs(file.c_str(), std::ios::out);

    if(!ofs)
    {
      GLV::ofs_error << "Fail to oepn " << file.c_str() << std::endl;
      std::exit(EXIT_FAILURE);
    }

    ofs << this->n_DMFT_atoms << '\n';
    ofs << this->n_inequivalent_atoms << '\n';
    for(int iatom=0; iatom<this->n_DMFT_atoms; iatom++)
    {
      ofs << std::setw(4) << iatom+1 << std::setw(4) << this->equivalent_atom[iatom] << std::endl;
    }
    ofs.close();

    //===============================
    // correlated_atoms.info
    //==============================
    file="correlated_atoms.info";
    ofs.open(file.c_str(),std::ios::out);
    if(!ofs)
    {
      GLV::ofs_error << "Fail to oepn " << file.c_str() << std::endl;
      std::exit(EXIT_FAILURE);
    }
    
    for(int iatom=0; iatom<this->n_DMFT_atoms; iatom++)
    {
      ofs << "atom  " << iatom << std::endl;
      ofs << "angular_moment  " << this->angular_momment[iatom] << std::endl;
      ofs << std::setw(6) << std::fixed << std::setprecision(2) << this->Hubbard_U[iatom]*GLC::Hartree_to_eV
      << std::setw(6) << std::fixed << std::setprecision(2) << this->Hund_J[iatom]*GLC::Hartree_to_eV
      << std::setw(4) << this->mag_moment[iatom]
      << std::setw(10) << std::fixed << std::setprecision(5) << this->occ_number[iatom][0]
      << std::setw(10) << std::fixed << std::setprecision(5) << this->occ_number[iatom][1] << std::endl;

      for(int m=0; m<2*this->angular_momment[iatom]+1; m++)
      {
        ofs << std::setw(2) << m-this->angular_momment[iatom] << "  " 
        << std::setw(2) << this->DMFT_orb_index[iatom][m] 
        << std::setw(6) << this->iorb_ibasis[this->DMFT_orb_index[iatom][m]] << std::endl;
      }
    }
    ofs.close();
  }

}
