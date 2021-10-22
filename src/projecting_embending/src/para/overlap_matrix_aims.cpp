#include "overlap_matrix_aims.h"
#include "../debug.h"

#include <fstream>
#include <iostream>
#include <cstdlib>   //Use exit function
#include <iomanip>   //use setw and setprecison function
#include <sstream>
#include <cmath>

namespace DFT_output
{
  namespace aims
  { 
    int overlap_matrix::ncells(){return n_cells_in_hamiltonian;}

    int overlap_matrix::kpoints(){return n_k_points;}

    std::complex<double> overlap_matrix::kphase(int icell, int ik)
    {
      if(icell<0 || icell>=this->n_cells_in_hamiltonian
        || ik<0 || ik>=this->n_k_points)
      {
        std::cout << "Array kphase is out of range" << std::endl;
        std::exit(EXIT_FAILURE);
      }

      return this->k_phase[icell][ik];
    }

    int overlap_matrix::index_ham(int rc, int icell, int i_basis)
    {
      if(rc<0 || rc>=2
        || icell<0 || icell>=this->n_cells_in_hamiltonian
        || i_basis<0 || i_basis>=this->n_basis)
      {
        std::cout << "Array index_hamiltonian is out of range" << std::endl;
        std::exit(EXIT_FAILURE);
      }

      return this->index_hamiltonian[rc][icell][i_basis];
    }

    int overlap_matrix::index_col_ham(int index)
    {
      if(index<0 || index>=size_colum_index_hamiltonian)
      {
        std::cout << "Array colum_index_hamiltonian is out of range" << std::endl;
        std::exit(EXIT_FAILURE);
      }

      return this->colum_index_hamiltonian[index];
    }

    double overlap_matrix::ovlp_mat_real(int index)
    {
      if(index<0 || index>=this->size_overlap_matrix)
      {
        std::cout << "Array ovlp_matrix_real is out of range" << std::endl;
        std::exit(EXIT_FAILURE);
      }

      return this->ovlp_matrix_real[index];
    }

    bool overlap_matrix::read(const std::string dir)
    { 
      debug::codestamp("overlap_matrix::read");

      int icell, ik, i_basis;
      int ival, ival1;
      double real, imag;
      std::string file;

      //=======================================
      //         k_phase
      //=======================================
      file=dir+"/k_phase.dat";
      std::ifstream if_kphase(file.c_str(), std::ios::in);

      if (!if_kphase) 
	    {
	    	std::cout << "Fail to oepn " << file.c_str() << std::endl;
        std::exit(EXIT_FAILURE);
      }

      if_kphase.seekg(0); //set the position at the beginning of the file

      if_kphase >> this->n_cells_in_hamiltonian;
      this->k_phase.resize(this->n_cells_in_hamiltonian);
      if_kphase >> this->n_k_points;
      for(icell=0; icell<this->n_cells_in_hamiltonian; icell++)
      {
        this->k_phase[icell].resize(this->n_k_points);
      }
      if_kphase.ignore(150,'\n');

      while(if_kphase.good())
      {
        if_kphase >> icell;
        if(if_kphase.eof()) break;
        if_kphase >> ik;
        if_kphase >> real;
        if_kphase >> imag;
        this->k_phase[icell][ik] = std::complex<double>(real,imag);
        if_kphase.ignore(150,'\n');

        if(if_kphase.eof()) break;   //Check whether end of file is reached 
      }
      if_kphase.close();

      //============================================
      //  index_hamiltonian
      //============================================
      file=dir+"/index_hamiltonian.dat";
      std::ifstream ifs(file.c_str(), std::ios::in);

      ifs.seekg(0); //set the position at the beginning of the file
      ifs >> ival;
      ifs >> ival;
      ifs >> this->n_basis;

      this->index_hamiltonian.resize(2);
      for(int i=0; i<2; i++)
      {
        this->index_hamiltonian[i].resize(this->n_cells_in_hamiltonian);
        for(icell=0; icell<this->n_cells_in_hamiltonian; icell++)
        {
          this->index_hamiltonian[i][icell].resize(this->n_basis);
        }
      }
      ifs.ignore(150,'\n');

      while(ifs.good())
      {
        ifs >> ival;
        if(ifs.eof()) break;
        ifs >> icell;
        ifs >> i_basis;
        ifs >> ival1;
        this->index_hamiltonian[ival][icell][i_basis] = ival1;
        ifs.ignore(150,'\n');

        if(ifs.eof()) break;
      }
      ifs.close();

      //====================================================
      //    colum_index_hamiltonian
      //====================================================
      file=dir+"/colum_index_hamiltonian.dat";
      ifs.open(file.c_str(), std::ios::in);

      ifs >> this->size_colum_index_hamiltonian;
      this->colum_index_hamiltonian.resize(this->size_colum_index_hamiltonian);
      ifs.ignore(150,'\n');

      while(ifs.good())
      {
        ifs >> ival;
        if(ifs.eof()) break;
        ifs >> ival1;
        this->colum_index_hamiltonian[ival] = ival1;
        ifs.ignore(150,'\n');

        if(ifs.eof()) break;
      }
      ifs.close();

      //==================================================
      //   overlap matrix
      //==================================================
      file=dir+"/overlap_matrix.dat";
      ifs.open(file.c_str(), std::ios::in);

      ifs >> this->size_overlap_matrix;
      this->ovlp_matrix_real.resize(this->size_overlap_matrix);
      ifs.ignore(150,'\n');

      while(ifs.good())
      {
        ifs >> ival;
        if(ifs.eof()) break;
        ifs >> this->ovlp_matrix_real[ival];
        ifs.ignore(150,'\n');

        if(ifs.eof()) break;
      }
      ifs.close();

      return true;	  	
    }
  
    void overlap_matrix::evaluate_ovlp_k(const int ik, atoms_info& at_info)
    {
      debug::codestamp("overlap_matrix::evaluate_ovlp_k");

      const int norb=at_info.norb();
      const std::complex<double> zero(0.0,0.0);

      if(this->ovlp_matrix_k.empty()) this->ovlp_matrix_k.resize(this->n_basis*norb);
      if(this->ovlp_matrix_work.empty()) this->ovlp_matrix_work.resize(this->n_basis*this->n_basis);
      if(this->ovlp_localorb_k.empty()) this->ovlp_localorb_k.resize(norb*norb);

      for(std::complex<double>& iter : this->ovlp_matrix_work) iter = zero;

      for(int ibasis_row=0; ibasis_row<this->n_basis; ibasis_row++)
      {
        for(int icell=0; icell<this->n_cells_in_hamiltonian-1; icell++)
        {
          if( this->index_ham(0,icell,ibasis_row) > -1 )
          {
            int index_real = this->index_ham(0,icell,ibasis_row)-1;

            for(int i_size=index_real+1; i_size<=this->index_ham(1,icell,ibasis_row); i_size++)
            {
              index_real += 1;
              int ibasis_col= this->index_col_ham(index_real);

              this->ovlp_matrix_work[ibasis_col*this->n_basis+ibasis_row] += 
                  this->ovlp_matrix_real[index_real]*this->k_phase[icell][ik];
            }//i_size
          }
        }//icell
      }//ibasis_row

      for(int ibasis_row=0; ibasis_row<this->n_basis; ibasis_row++)
      {
        for(int ibasis_col=ibasis_row+1; ibasis_col<this->n_basis; ibasis_col++)
        {
          this->ovlp_matrix_work[ibasis_col*this->n_basis+ibasis_row] = 
             std::conj(this->ovlp_matrix_work[ibasis_row*this->n_basis+ibasis_col]);
        }
      }

      for(int iorb=0; iorb<norb; iorb++)
      {
        const int ibasis_col = at_info.iorb2ibasis(iorb);
        
        for(int ibasis_row=0; ibasis_row<this->n_basis; ibasis_row++)
        {
          this->ovlp_matrix_k[ibasis_row*norb+iorb] = 
            this->ovlp_matrix_work[ibasis_row*this->n_basis+ibasis_col];
        }//icell
      }//iorb

      for(int iorb1=0; iorb1<norb; iorb1++)
      {   
        for(int iorb2=0; iorb2<norb; iorb2++)
        {
          this->ovlp_localorb_k[iorb1*norb+iorb2] = 
            this->ovlp_matrix_work[ at_info.iorb2ibasis(iorb1)*
            this->n_basis + at_info.iorb2ibasis(iorb2) ];
        }//iorb2
      }//iorb1

// std::stringstream ss;
// ss << "overlap_matrix/ovlp_ik" << ik << ".dat";
// std::string file=ss.str();
// std::ofstream ofs(file.c_str(),std::ios::out);
// for(int ibasis_row=0; ibasis_row<this->n_basis; ibasis_row++)
// {
//   for(int ibasis_col=0; ibasis_col<this->n_basis; ibasis_col++)
//   {
//     ofs << std::setw(5) << ibasis_row << std::setw(5) <<  ibasis_col 
//     << std::setw(15) << std::fixed << std::setprecision(9) << this->ovlp_matrix_work[ibasis_col*this->n_basis+ibasis_row].real()
//     << std::setw(15) << std::fixed << std::setprecision(9) << this->ovlp_matrix_work[ibasis_col*this->n_basis+ibasis_row].imag() << '\n';
//   }
// }
// ofs.close();

      return;
    }
    
    void overlap_matrix::out()
    {
      debug::codestamp("overlap_matrix::out");

      std::string dir="./overlap_matrix";
      std::string file;

      std::stringstream ss;
      ss << "test -d " << dir << " || mkdir " << dir;
      system(ss.str().c_str());

      //=======================================
      //         k_phase
      //=======================================
      file="./overlap_matrix/k_phase.dat";
      std::ofstream ofs(file.c_str(), std::ios::out);
      if(!ofs)
      {
        std::cout << "Fail to oepn " << file.c_str() << std::endl;
        std::exit(EXIT_FAILURE);
      }
      
      ofs << std::setw(5) << this->n_cells_in_hamiltonian << std::setw(6) << this->n_k_points << '\n';
      for(int icell=0; icell<this->n_cells_in_hamiltonian; icell++)
      {
        for(int ik=0; ik<this->n_k_points; ik++)
        {
          ofs << std::setw(5) << icell << std::setw(6) << ik
          << std::setw(15) << std::fixed << std::setprecision(9) << this->k_phase[icell][ik].real()
          << std::setw(15) << std::fixed << std::setprecision(9) << this->k_phase[icell][ik].imag() << '\n';
        }
      }
      ofs.close();

      //============================================
      //  index_hamiltonian
      //============================================
      file="./overlap_matrix/index_hamiltonian.dat";
      ofs.open(file.c_str(),std::ios::out);
      if(!ofs)
      {
        std::cout << "Fail to oepn " << file.c_str() << std::endl;
        std::exit(EXIT_FAILURE);
      }
      
      ofs << std::setw(2) << 2 << std::setw(5) << this->n_cells_in_hamiltonian 
      << std::setw(6) << this->n_basis << '\n';

      for(int i=0; i<2; i++)
      {
        for(int icell=0; icell<this->n_cells_in_hamiltonian; icell++)
        {
          for(int ibasis=0; ibasis<this->n_basis; ibasis++)
          {
            ofs << std::setw(2) << i << std::setw(5) << icell 
            << std::setw(6) << ibasis
            << std::setw(12) << this->index_hamiltonian[i][icell][ibasis] << '\n';
          }
        }
      }
      ofs.close();

      //============================================
      //  index_hamiltonian
      //============================================
      file="./overlap_matrix/colum_index_hamiltonian.dat";
      ofs.open(file.c_str(),std::ios::out);
      if(!ofs)
      {
        std::cout << "Fail to oepn " << file.c_str() << std::endl;
        std::exit(EXIT_FAILURE);
      }

      ofs << this->size_colum_index_hamiltonian << '\n';

      for(long i=0; i<this->size_colum_index_hamiltonian; i++)
      {
        ofs << std::setw(12) << i << std::setw(6) << this->colum_index_hamiltonian[i] << '\n';
      }
      ofs.close();

      //============================================
      //    overlap matrix
      //============================================
      file="./overlap_matrix/overlap_matrix.dat";
      ofs.open(file.c_str(),std::ios::out);
      if(!ofs)
      {
        std::cout << "Fail to oepn " << file.c_str() << std::endl;
        std::exit(EXIT_FAILURE);
      }

      ofs << std::setw(12) << this->size_overlap_matrix << '\n';

      for(long i=0; i<this->size_overlap_matrix; i++)
      {
        ofs << std::setw(12) << i
        << std::setw(18) << std::fixed << std::setprecision(12) << this->ovlp_matrix_real[i] << std::endl;
      }
      ofs.close();

    }

  }
}
