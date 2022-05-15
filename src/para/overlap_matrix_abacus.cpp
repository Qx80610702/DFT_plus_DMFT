#include "overlap_matrix_abacus.h"
#include "../constants.h"
#include "../debug.h"
#include "../global_variables.h"

#include <fstream>
#include <cstdlib>
#include <cmath>
#include <iostream>

//tests
#include <sstream>
#include <iomanip>

namespace DFT_output
{
  namespace abacus
  {
    overlap_matrix::overlap_matrix(const std::string dir)
    {
      this->read(dir);
    }

    void overlap_matrix::read(const std::string dir)
    {
      debug::codestamp("DFT_output::abacus::overlap_matrix::read");

      int nkpoints, ik_point;
      char words[150];
      int size;

      //=======================================
      //         k_vector
      //=======================================
      std::string kfile = dir+"/k_vector.dat";
      std::ifstream ifk(kfile.c_str(), std::ios::in);

      if (!ifk) 
	    {
	    	GLV::ofs_error << "Fail to oepn " << kfile.c_str() << std::endl;
        std::exit(EXIT_FAILURE);
      }

      ifk.seekg(0);  //set the position at the beginning of the file

      ifk >> nkpoints;
      ifk.ignore(150,'\n');

      this->kvector.resize(nkpoints);
      for(auto& iter : this->kvector) iter.resize(3);

      while(ifk.good())
      {
        ifk >> ik_point;
        if(ifk.eof()) break;
        ifk >> this->kvector[ik_point][0];
        ifk >> this->kvector[ik_point][1];
        ifk >> this->kvector[ik_point][2];
        ifk.ignore(150,'\n');

        if(ifk.eof()) break;    //Check whether the end of file is reached 
      }
      ifk.close();

      //===============================================
      //        overlap matrix in R space: SR.csr
      //===============================================
      std::string Srf = dir + "/SR.csr";
      std::ifstream ifSr(Srf.c_str(), std::ios::in);

      if (!ifSr) 
	    {
	    	GLV::ofs_error << "Fail to oepn " << Srf.c_str() << std::endl;
        std::exit(EXIT_FAILURE);
      }

      ifSr.getline(words, 150);
      std::string line1(words);
      this->get_last_number_in_line(line1, this->nbasis);

      ifSr.getline(words, 150);
      std::string line2(words);
      this->get_last_number_in_line(line2, this->ncells);

      this->R_vector.resize(this->ncells);
      for(auto& iter : this->R_vector) iter.resize(3);

      this->ovlp_csr.resize(this->ncells);
      this->col_indices.resize(this->ncells);

      this->row_indptr.resize(this->ncells);
      for(auto& iter : this->row_indptr) iter.resize(this->nbasis+1);

      for(int icell=0; icell<this->ncells; icell++)
      {
        ifSr >> this->R_vector[icell][0];
        ifSr >> this->R_vector[icell][1];
        ifSr >> this->R_vector[icell][2];
        ifSr >> size;
        ifSr.ignore(150,'\n');

        this->ovlp_csr[icell].resize(size);
        for(int i=0; i<size; i++)
          ifSr >> this->ovlp_csr[icell][i];
        ifSr.ignore(150,'\n');

        this->col_indices[icell].resize(size);
        for(int i=0; i<size; i++)
          ifSr >> this->col_indices[icell][i];
        ifSr.ignore(150,'\n');

        for(int i=0; i<this->nbasis+1; i++)
          ifSr >> this->row_indptr[icell][i];
        ifSr.ignore(150,'\n');
      }

      return;
    }

    void overlap_matrix::get_last_number_in_line(std::string& line, int& val)
    {
      line.erase(0, line.find_first_not_of(" "));
      line.erase(line.find_last_not_of(" ") + 1);

      line.erase(0, line.find_last_of(' '));
      
      val = atoi(line.c_str());
    }

    void overlap_matrix::evaluate_ovlp_k(const int ik, atoms_info& atom)
    {
      debug::codestamp("DFT_output::abacus::overlap_matrix::evaluate_ovlp_k");

      const int norb = atom.norb();
      const std::complex<double> zero(0.0,0.0);

std::vector<std::complex<double>> Sk(this->nbasis*this->nbasis, zero);

      if(this->ovlp_matrix_k.empty()) this->ovlp_matrix_k.resize(this->nbasis*norb);
      if(this->ovlp_localorb_k.empty()) this->ovlp_localorb_k.resize(norb*norb);

      for(std::complex<double>& iter : this->ovlp_matrix_k) iter = zero;

      for(int icell=0; icell<this->ncells; icell++)
      {
        double arg = 2.0*GLC::PI*( this->kvector[ik][0]*this->R_vector[icell][0] +
                              this->kvector[ik][1]*this->R_vector[icell][1] +
                              this->kvector[ik][2]*this->R_vector[icell][2] );
        std::complex<double> kphase = std::complex<double>(std::cos(arg), -std::sin(arg));

        for(int iorb=0; iorb<norb; iorb++)
        {
          int ibasis_row = atom.iorb2ibasis(iorb);
          for(int index=this->row_indptr[icell][ibasis_row]; index<this->row_indptr[icell][ibasis_row+1]; index++)
            this->ovlp_matrix_k[ this->col_indices[icell][index]*norb + iorb] 
              += kphase*this->ovlp_csr[icell][index];
        }

        // for(int ibasis_row=0; ibasis_row<this->nbasis; ibasis_row++)
        //   for(int index=this->row_indptr[icell][ibasis_row]; index<this->row_indptr[icell][ibasis_row+1]; index++)
        //     this->ovlp_matrix_work[ ibasis_row*this->nbasis + this->col_indices[icell][index]] 
        //           += kphase*this->ovlp_csr[icell][index];

//tests
// std::stringstream ss;
// ss << "overlap_matrix/Sk" << ik << ".dat";
// std::ofstream ofs(ss.str().c_str(), std::ios::out);
// for(int ir=0; ir<this->nbasis; ir++)
//   for(int ic=0; ic<norb; ic++)
//     ofs << std::setw(20) << std::fixed << std::setprecision(12) 
//         << this->ovlp_matrix_work[ir*norb+ic].real()
//         << std::setw(20) << std::fixed << std::setprecision(12) 
//         << this->ovlp_matrix_work[ir*norb+ic].imag() << std::endl;
// ofs.close();

      }

      for(int iorb1=0; iorb1<norb; iorb1++)
        for(int iorb2=0; iorb2<norb; iorb2++)
          this->ovlp_localorb_k[iorb1*norb+iorb2] = 
            this->ovlp_matrix_k[ atom.iorb2ibasis(iorb1)*norb + iorb2 ];

      return;
    }

  }
}
