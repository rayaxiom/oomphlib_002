//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC//           Version 0.90. August 3, 2009.
//LIC// 
//LIC// Copyright (C) 2006-2009 Matthias Heil and Andrew Hazel
//LIC// 
//LIC// This library is free software; you can redistribute it and/or
//LIC// modify it under the terms of the GNU Lesser General Public
//LIC// License as published by the Free Software Foundation; either
//LIC// version 2.1 of the License, or (at your option) any later version.
//LIC// 
//LIC// This library is distributed in the hope that it will be useful,
//LIC// but WITHOUT ANY WARRANTY; without even the implied warranty of
//LIC// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//LIC// Lesser General Public License for more details.
//LIC// 
//LIC// You should have received a copy of the GNU Lesser General Public
//LIC// License along with this library; if not, write to the Free Software
//LIC// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
//LIC// 02110-1301  USA.
//LIC// 
//LIC// The authors may be contacted at oomph-lib@maths.man.ac.uk.
//LIC// 
//LIC//====================================================================
#include "block_preconditioner.h"

namespace oomph
{

 /// \short Static boolean to allow block_matrix_test(...) to be run.
 /// Defaults to false.
 template<typename MATRIX> 
 bool BlockPreconditioner<MATRIX>::Run_block_matrix_test=false;
 
//=============================================================================
/// \short Gets block (i,j) from the original matrix and returns it in
/// block_matrix_pt (Specialisation for CCDoubleMatrix)
//=============================================================================
 template<> 
 void BlockPreconditioner<CCDoubleMatrix>:: 
 get_block(const unsigned& i, const unsigned& j, 
	    CCDoubleMatrix*& block_pt) const
 {

  // pointers for the jacobian matrix is compressed column sparse format 
  int* j_column_start;
  int* j_row_index;
  double* j_value;
  
    // Cast the matrix pointer
    CCDoubleMatrix* cc_matrix_pt = dynamic_cast<CCDoubleMatrix*>(matrix_pt());

  // sets pointers to jacobian matrix
    j_column_start = cc_matrix_pt->column_start();
    j_row_index = cc_matrix_pt->row_index();
    j_value = cc_matrix_pt->value();

  // get the block dimensions
  unsigned block_nrow = this->block_dimension(i);
  unsigned block_ncol = this->block_dimension(j);

  // allocate temporary storage for the component vectors of block (i,j)
  // temp_ptr is used to point to an element in each column - required as
  // cannot assume that order of block's rows in jacobian and the block
  // matrix will be the same
  Vector<int> temp_row_index;
  Vector<int> temp_column_start(block_ncol+1);
  Vector<int> temp_ptr(block_ncol+1);
  Vector<double> temp_value;
  int block_nnz = 0;
  
  // get number of rows in source matrix
  unsigned master_nrow = this->master_nrow();
  
  // determine how many non zeros there are in the block (i,j)
  // also determines how many non zeros are stored in each row or column - 
  // stored in temp_ptr temporarily
  for (unsigned k = 0; k < master_nrow; k++)
   {
    if (block_number(k) == static_cast<int>(j))
     {
      for (int l = j_column_start[k]; 
           l < j_column_start[k+1]; l++)
       {
        if (block_number(j_row_index[l]) == 
            static_cast<int>(i))
         {
          block_nnz++;
          temp_ptr[index_in_block(k)+1]++;
         }
       }
     }
   }

  // if the block matrix is not empty
  if (block_nnz > 0)
   {

    // uses number of elements in each column of block to determine values
    // for the block column start (temp_column_start)
    temp_column_start[0] = 0;
    for (unsigned k = 1; k <= block_ncol; k++)
     {
      
      temp_column_start[k] = temp_column_start[k-1]+temp_ptr[k];
      temp_ptr[k] = temp_column_start[k];
     }
    
    // resizes the block row index and value to store all non zeros
    temp_row_index.resize(block_nnz);
    temp_value.resize(block_nnz);
    
    // copies the relevant elements of the jacobian to the correct entries 
    // of the block matrix
    for (unsigned k = 0; k < master_nrow; k++)
     {
      if (block_number(k) == static_cast<int>(j))
       {
        for (int l = j_column_start[k]; 
             l < j_column_start[k+1]; l++)
         {
          if (block_number(j_row_index[l]) == 
              static_cast<int>(i))
           {
            int kk = temp_ptr[index_in_block(k)]++;
            temp_value[kk] = j_value[l];
            temp_row_index[kk] = index_in_block(j_row_index[l]); 
           }
         }
       }
     }								
    
    // creates a new compressed column sparse matrix for the pointer for
    // the current block 
    block_pt = new CCDoubleMatrix(temp_value,temp_row_index,
                                  temp_column_start,block_nrow,block_ncol);

#ifdef PARANOID
    if (Run_block_matrix_test)
     {
      // checks to see if block matrix has been set up correctly 
	    block_matrix_test(i,j,block_pt);
     }
#endif
   }
 
  // else the matrix is empty
  else
   {
    block_pt = 0;
   }
 }


//=============================================================================
/// \short Gets block (i,j) from the original matrix and returns it in
/// block_matrix_pt (Specialisation for CRDoubleMatrix)
//=============================================================================
 template<> 
 void BlockPreconditioner<CRDoubleMatrix>:: 
 get_block(const unsigned& block_i, const unsigned& block_j, 
	    CRDoubleMatrix*& block_pt) const
 {

#ifdef PARANOID
  // the number of blocks
  unsigned n_blocks = this->nblock_types();

  // paranoid check that block i is in this block preconditioner
  if (block_i >= n_blocks || block_j >= n_blocks)
   {
    std::ostringstream error_message;
    error_message << "Requested block (" << block_i << "," << block_j   
                  << "), however this preconditioner has nblock_types() "
                  << "= " << nblock_types() << std::endl;
    throw OomphLibError(error_message.str(),
                        "BlockPreconditioner::get_block(...)",
                        OOMPH_EXCEPTION_LOCATION);
   }
#endif

    // Cast the pointer
    CRDoubleMatrix* cr_matrix_pt = dynamic_cast<CRDoubleMatrix*>(matrix_pt());

  // if + only one processor
  //    + more than one processor but matrix_pt is not distributed
  // then use the serial get_block method
    if (cr_matrix_pt->distribution_pt()->communicator_pt()->nproc() == 1 ||
	!cr_matrix_pt->distribution_pt()->distributed())
   {
    // pointers for the jacobian matrix is compressed row sparse format 
    int* j_row_start;
    int* j_column_index;
    double* j_value;
    
    // sets pointers to jacobian matrix
	j_row_start = cr_matrix_pt->row_start();
	j_column_index = cr_matrix_pt->column_index();
	j_value = cr_matrix_pt->value();
    
    // get the block dimensions
    unsigned block_nrow = this->block_dimension(block_i);
    unsigned block_ncol = this->block_dimension(block_j);
    
    // allocate temporary storage for the component vectors of block (i,j)
    // temp_ptr is used to point to an element in each column - required as
    // cannot assume that order of block's rows in jacobian and the block
    // matrix will be the same
    int* temp_row_start = new int[block_nrow+1];
    for (unsigned i = 0; i <= block_nrow; i++)
     {
      temp_row_start[i] = 0;
     }
    Vector<int> temp_ptr(block_nrow+1);
    int block_nnz = 0;
    
    // get number of rows in source matrix
    unsigned master_nrow = this->master_nrow();
    
    // determine how many non zeros there are in the block (i,j)
    // also determines how many non zeros are stored in each row or column - 
    // stored in temp_ptr temporarily
    for (unsigned k = 0; k < master_nrow; k++)
     {
      if (block_number(k) == static_cast<int>(block_i))
       {
        for (int l = j_row_start[k]; 
             l < j_row_start[k+1]; l++)
         {
          if (block_number(j_column_index[l]) == 
              static_cast<int>(block_j))
           {
            block_nnz++;
            temp_ptr[index_in_block(k)+1]++;
           }
         }
       }
     }
    
    // if the matrix is not empty
    int* temp_column_index = new int[block_nnz];
    double* temp_value = new double[block_nnz];
    if (block_nnz > 0)
     {
      
      // uses number of elements in each column of block to determine values
      // for the block column start (temp_row_start)
      temp_row_start[0] = 0;
      for (unsigned k = 1; k <= block_nrow; k++)
       {
        temp_row_start[k] = temp_row_start[k-1]+temp_ptr[k];
        temp_ptr[k] = temp_row_start[k];
       }
      
      // copies the relevant elements of the jacobian to the correct entries 
      // of the block matrix
      for (unsigned k = 0; k < master_nrow; k++)
       {
        if (block_number(k) == static_cast<int>(block_i))
         {
          for (int l = j_row_start[k]; 
               l < j_row_start[k+1]; l++)
           {
            if (block_number(j_column_index[l]) == 
                static_cast<int>(block_j))
             {
              int kk = temp_ptr[index_in_block(k)]++;
              temp_value[kk] = j_value[l];
              temp_column_index[kk] = 
               index_in_block(j_column_index[l]); 
             }
           }
         }
       }
     }
      
      
    // creates a new compressed column sparse matrix for the pointer for
    // the current block 
    if (block_pt == 0)
     {
      block_pt = new CRDoubleMatrix(Block_distribution_pt[block_i]);
     }
    block_pt->build_without_copy(block_ncol,block_nnz,
                                 temp_value,temp_column_index,
                                 temp_row_start);
 
#ifdef PARANOID
    // checks to see if block matrix has been set up correctly 
    //   block_matrix_test(matrix_pt,block_i,block_j,block_pt);
    if (Run_block_matrix_test)
     {
      // checks to see if block matrix has been set up correctly 
	    block_matrix_test(block_i,block_j,block_pt);
     }
#endif 
   }


  // otherwise we are dealing with a distributed matrix
  else
   {
#ifdef OOMPH_HAS_MPI
    // number of processors
    unsigned nproc = this->distribution_pt()->communicator_pt()->nproc();

    // my rank
    unsigned my_rank = this->distribution_pt()->communicator_pt()->my_rank();

    // sets pointers to jacobian matrix
	int* j_row_start = cr_matrix_pt->row_start();
	int* j_column_index = cr_matrix_pt->column_index();
	double* j_value = cr_matrix_pt->value();

    // number of non zeros in each row to be sent
    Vector<int*> nnz_send(nproc,0);

    // number of non zeros in each row to be received
    Vector<int*> nnz_recv(nproc,0);

    // storage for data to be sent
    Vector<int*> column_index_for_proc(nproc,0);
    Vector<double*> values_for_proc(nproc,0);

    // number of non zeros to be sent to each processor
    Vector<unsigned> total_nnz_send(nproc,0);

    // number of rows of the block matrix on this processor
    unsigned nrow_local = Block_distribution_pt[block_i]->nrow_local();

    // resize the nnz storage and compute nnz_send
    // and send and recv the nnz
    Vector<MPI_Request> send_req;
    Vector<MPI_Request> recv1_req;
    for (unsigned p = 0; p < nproc; p++)
     {
      int nrow_send = Nrows_to_send_for_get_block(block_i,p);
      int nrow_recv = Nrows_to_recv_for_get_block(block_i,p);

      // assemble nnz recv
      nnz_recv[p] = new int[nrow_recv];

      // assemble the storage to send
      if (nrow_send > 0 && p != my_rank)
       {
        nnz_send[p] = new int[nrow_send];
       }

      // compute the number of nnzs in each row and the total number
      // of nnzs
      for (int i = 0; i < nrow_send; i++)
       {
        unsigned row = Rows_to_send_for_get_block(block_i,p)[i];
        int c = 0;
        for (int r = j_row_start[row]; r < j_row_start[row+1]; r++)
         {
          if (block_number(j_column_index[r]) == int(block_j))
           {
            c++;
           }
         }
        if (p != my_rank)
         {
          nnz_send[p][i] = c;
         }
        else
         {
          nnz_recv[p][i] = c;
         }
        total_nnz_send[p] += c;
       }

      // send
      if (p != my_rank)
       {
        if (nrow_send)
         {
          MPI_Request req;
          MPI_Isend(nnz_send[p],nrow_send,MPI_INT,p,0,
                    this->distribution_pt()->communicator_pt()->mpi_comm(),
                    &req);
          send_req.push_back(req);
         }

        // recv
        if (nrow_recv)
         {
          MPI_Request req;
          MPI_Irecv(nnz_recv[p],nrow_recv,MPI_INT,p,0,
                    this->distribution_pt()->communicator_pt()->mpi_comm(),
                    &req);
          recv1_req.push_back(req);
         }
       }
     }

    // next assemble the values and row_start data to be sent for each
    // processor
    for (unsigned p = 0; p < nproc; p++)
     {
      int nrow_send = Nrows_to_send_for_get_block(block_i,p);

      // assemble the storage for the values and column indices to be sent
      if (p != my_rank)
       {
        if (total_nnz_send[p] > 0)
         {
          values_for_proc[p] = new double[total_nnz_send[p]];
          column_index_for_proc[p] = new int[total_nnz_send[p]];
          
          // copy the values and column indices to the storage
          unsigned ptr = 0;
          for (int i = 0; i < nrow_send; i++)
           {
            unsigned row = Rows_to_send_for_get_block(block_i,p)[i];
            for (int r = j_row_start[row]; r < j_row_start[row+1]; r++)
             {
              if (block_number(j_column_index[r]) == int(block_j))
               {
                values_for_proc[p][ptr] = j_value[r];
                column_index_for_proc[p][ptr] = 
                 index_in_block(j_column_index[r]);
                ptr++;
               }
             }
           }
       
          // create the datatypes
          MPI_Datatype types[2];
          MPI_Type_contiguous(total_nnz_send[p],MPI_DOUBLE,&types[0]);
          MPI_Type_commit(&types[0]);
          MPI_Type_contiguous(total_nnz_send[p],MPI_INT,&types[1]);
          MPI_Type_commit(&types[1]);
          
          // get the start address of the vectors
          MPI_Aint displacement[2];
          MPI_Address(values_for_proc[p],&displacement[0]);
          MPI_Address(column_index_for_proc[p],&displacement[1]);
          
          // compute the displacements
          displacement[1] -= displacement[0];
          displacement[0] -= displacement[0];

          // compute the block lengths
          int length[2];
          length[0] = length[1] = 1;

          // build the struct data type
          MPI_Datatype final_type;
          MPI_Type_struct(2,length,displacement,types,&final_type);
          MPI_Type_commit(&final_type);
          MPI_Type_free(&types[0]);
          MPI_Type_free(&types[1]);

          // and send
          MPI_Request req;
          MPI_Isend(values_for_proc[p],1,final_type,p,1,
                    this->distribution_pt()->communicator_pt()->mpi_comm(),
                    &req);
          send_req.push_back(req);
          MPI_Type_free(&final_type);
         }
       }
     }

    // wait for the recv to complete (the row_start recv which actually
    // contains the number of nnzs in each row)
    int c_recv = recv1_req.size();
    if (c_recv != 0)
     {
      MPI_Waitall(c_recv,&recv1_req[0],MPI_STATUS_IGNORE);
     }

    // compute the total number of nnzs to be received
    Vector<int> total_nnz_recv_from_proc(nproc);
    int local_block_nnz = 0;
    for (unsigned p = 0; p < nproc; p++)
     {
      // compute the total nnzs
      for (unsigned i = 0; i < Nrows_to_recv_for_get_block(block_i,p); i++)
       {
        total_nnz_recv_from_proc[p] += nnz_recv[p][i];

       }
      local_block_nnz += total_nnz_recv_from_proc[p];
     }

    // compute the offset for each block of nnzs (a matrix row) in the 
    // values_recv and column_index_recv vectors

    // fisrt determine how many blocks of rows are to be recv
    Vector<int> n_recv_block(nproc,0);
    for (unsigned p = 0; p < nproc; p++)
     {
      if (Nrows_to_recv_for_get_block(block_i,p) > 0)
       {
        n_recv_block[p] = 1;
       }
      for (unsigned i = 1; i < Nrows_to_recv_for_get_block(block_i,p); i++)
       {
        if (Rows_to_recv_for_get_block(block_i,p)[i] !=
            Rows_to_recv_for_get_block(block_i,p)[i-1] + 1)
         {
          n_recv_block[p]++;
         }
       }
     }

    // next assemble row start recv
    int* row_start_recv = new int[nrow_local+1];
    for (unsigned i = 0; i <= nrow_local; i++)
     {
      row_start_recv[i] = 0;
     }
    for (unsigned p = 0; p < nproc; p++)
     {
      for (unsigned i = 0; i < Nrows_to_recv_for_get_block(block_i,p); i++)   
       {
        row_start_recv[Rows_to_recv_for_get_block(block_i,p)[i]] 
         = nnz_recv[p][i];
       }
     }
    int g = row_start_recv[0];
    row_start_recv[0] = 0;
    for (unsigned i = 1; i < nrow_local; i++)
     {
      int temp_g = g;
      g = row_start_recv[i];
      row_start_recv[i] = row_start_recv[i-1] + temp_g;
     }    
    row_start_recv[nrow_local] = row_start_recv[nrow_local-1] + g;

    // next assemble the offset and the number of nzs in each recv block
    Vector<int*> offset_recv_block(nproc,0);
    Vector<int*> nnz_recv_block(nproc,0);
    for (unsigned p = 0; p < nproc; p++)
     {
      if (Nrows_to_recv_for_get_block(block_i,p) > 0)
       {
        offset_recv_block[p] = new int[n_recv_block[p]];
        offset_recv_block[p][0] = 0;
        nnz_recv_block[p] = new int[n_recv_block[p]];
        for (int i = 0; i < n_recv_block[p]; i++)
         {
          nnz_recv_block[p][i] = 0;
         }
        unsigned ptr = 0;
        nnz_recv_block[p][ptr] += nnz_recv[p][0];
        offset_recv_block[p][0] 
         = row_start_recv[Rows_to_recv_for_get_block(block_i,p)[0]];
        for (unsigned i = 1; i < Nrows_to_recv_for_get_block(block_i,p); i++)
         {
          if (Rows_to_recv_for_get_block(block_i,p)[i] !=
              Rows_to_recv_for_get_block(block_i,p)[i-1] + 1)
           {
            ptr++;
            offset_recv_block[p][ptr] 
             = row_start_recv[Rows_to_recv_for_get_block(block_i,p)[i]];
           }
          nnz_recv_block[p][ptr] += nnz_recv[p][i];
         }
       }
      delete[] nnz_recv[p];
     }

    // post the receives
    int* column_index_recv = new int[local_block_nnz];
    double* values_recv = new double[local_block_nnz];
    Vector<MPI_Request> recv2_req;
    for (unsigned p = 0; p < nproc; p++)
     {
      if (p != my_rank)
       {
        if (total_nnz_recv_from_proc[p] != 0)
         {
          // create the datatypes
          MPI_Datatype types[2];
          MPI_Type_indexed(n_recv_block[p],nnz_recv_block[p],         
                           offset_recv_block[p],MPI_DOUBLE,&types[0]);
          MPI_Type_commit(&types[0]);
          MPI_Type_indexed(n_recv_block[p],nnz_recv_block[p],         
                           offset_recv_block[p],MPI_INT,&types[1]);
          MPI_Type_commit(&types[1]);
          
          // compute the displacements
          MPI_Aint displacements[2];
          MPI_Address(values_recv,&displacements[0]);
          MPI_Address(column_index_recv,&displacements[1]);
          displacements[1] -= displacements[0];
          displacements[0] -= displacements[0];
          
          // compute the block lengths
          int length[2];
          length[0] = length[1] = 1;
          
          // create the final datatype
          MPI_Datatype final_type;
          MPI_Type_struct(2,length,displacements,types,&final_type);
          MPI_Type_commit(&final_type);
	  MPI_Type_free(&types[0]);
	  MPI_Type_free(&types[1]);
          
          // and the recv
          MPI_Request req;
          MPI_Irecv(values_recv,1,final_type,p,1,
                    this->distribution_pt()->communicator_pt()->mpi_comm(),
                    &req);
          recv2_req.push_back(req);
	  MPI_Type_free(&final_type);
         }
       }
      else
       {
        // next send the values and column indices to self
        unsigned block_ptr = 0;
        unsigned counter = 0;
        int nrow_send = Nrows_to_send_for_get_block(block_i,my_rank);
        if (nrow_send > 0)
         {
          unsigned offset = offset_recv_block[my_rank][0];
          for (int i = 0; i < nrow_send; i++)
           {
            if (i > 0)
             {
              if (Rows_to_recv_for_get_block(block_i,p)[i] !=
                  Rows_to_recv_for_get_block(block_i,p)[i-1] + 1)
               {
                counter = 0;
                block_ptr++;
                offset = offset_recv_block[my_rank][block_ptr];
               }
             }
            unsigned row = Rows_to_send_for_get_block(block_i,my_rank)[i];
            for (int r = j_row_start[row]; r < j_row_start[row+1]; r++)
             {
              if (block_number(j_column_index[r]) == int(block_j))
               {
                values_recv[offset+counter] = j_value[r];
                column_index_recv[offset + counter] = 
                 index_in_block(j_column_index[r]);
                counter++;
               }
             }
           }
         }
       }
     }
       
    // wait for the recv to complete (for the column_index and the values_
    c_recv = recv2_req.size();
    if (c_recv != 0)
     {
      MPI_Waitall(c_recv,&recv2_req[0],MPI_STATUS_IGNORE);   
     }

    // create the matrix
    // creates a new compressed column sparse matrix for the pointer for
    // the current block 
    if (block_pt == 0)
     {
      block_pt = new CRDoubleMatrix(Block_distribution_pt[block_i]);
     }
    block_pt->build_without_copy(this->block_dimension(block_j),
                                 local_block_nnz,
                                 values_recv,column_index_recv,
                                 row_start_recv);
    
    // wait for the send to complete (nnz / row_start)
    int c_send = send_req.size();
    if (c_send)
     {
      MPI_Waitall(c_send,&send_req[0],MPI_STATUS_IGNORE);
     }

    // delete temp storage used for assembling data for communication
    for (unsigned p = 0; p < nproc; p++)
     {
      delete[] nnz_send[p];
      delete[] column_index_for_proc[p];
      delete[] values_for_proc[p];
      delete[] offset_recv_block[p];
      delete[] nnz_recv_block[p];
     }
#else
    // throw error
    std::ostringstream error_message;
    error_message << "The matrix is distributed and on more than one "
                  << "processor. MPI is required.";
    throw OomphLibError(error_message.str(),
                        "BlockPreconditioner<MATRIX>::get_block(...)",
                        OOMPH_EXCEPTION_LOCATION);
#endif 
   }
 }
 

//=============================================================================
/// \short test function to check that every element in the block matrix
/// (block_i,block_j) matches the corresponding element in the original matrix
//=============================================================================
  template<typename MATRIX> void BlockPreconditioner<MATRIX>::
  block_matrix_test(const unsigned& block_i, const unsigned& block_j,
		    const MATRIX* block_matrix_pt) const
 {

  // boolean flag to indicate whether test is passed
  bool check = true;
  
  // number of rows in matrix
    unsigned n_row = matrix_pt()->nrow();
  
  // number of columns in matrix
    unsigned n_col = matrix_pt()->ncol();
  
  // loop over rows of original matrix
  for (unsigned i = 0; i < n_row; i++)
   {
    
    // if this coefficient is associated with a block in this block 
    // preconditioner
    if (static_cast<int>(block_i) == this->block_number(i))
     {
      
      // loop over columns of original matrix
      for (unsigned j = 0; j < n_col; j++)
       {
        
        // if the coeeficient is associated with a block in this block
        // preconditioner
        if (static_cast<int>(block_j) == this->block_number(j))
         {
          
          // check whether elements in original matrix and matrix of block 
          // pointers match
		    if ( matrix_pt()->operator()(i,j) !=
               block_matrix_pt
               ->operator()(index_in_block(i),index_in_block(j)) )
           {
            check = false;
           }
         }
       }
     }
   }
  
  // throw error
  if (!check)
   {
    std::ostringstream error_message;
    error_message << "The require elements have not been successfully copied"
                  << " from the original matrix to the block matrices";
    throw OomphLibError(error_message.str(),
                        "BlockPreconditioner<MATRIX>::block_matrix_test()",
                        OOMPH_EXCEPTION_LOCATION);
   }
 }


//=============================================================================
/// Assemble the block preconditioner as a single matrix. 
/// In some cases the block preconditioner cannot be applied in 
/// individual components; this function takes the matrix
/// of block pointers and returns a single matrix containing all the 
/// blocks of the matrix of blocks in a single matrix that can
/// be solved directly. Specialised version for CCDoubleMatrix.
//=============================================================================
 template<> 
  void BlockPreconditioner<CCDoubleMatrix>::
  build_preconditioner_matrix(DenseMatrix<CCDoubleMatrix*>& block_matrix_pt,
			      CCDoubleMatrix*& preconditioner_matrix_pt) const
 {


  // number of non zeros in the final preconditioner matrix
  int p_nnz = 0;
  for (unsigned j = 0; j < Nblock_types; j++)
   {
    for (unsigned i = 0; i < Nblock_types; i++)
     {
      if (block_matrix_pt(i,j) != 0)
       {
        p_nnz += block_matrix_pt(i,j)->nnz();	
       }
     }
   }
    
  // vector indicating the offset required for a block
  Vector<unsigned> block_start(Nblock_types,0);
  for (unsigned a = 1; a < Nblock_types; a++)
   {
    block_start[a] = block_start[a-1]+block_dimension(a-1);
   }
    
  // determine p matrix n_row
  unsigned p_nrow = Nrow;
    
  // temporary storage for complete block preconditioner matrix
  Vector<int> p_column_start(p_nrow+1);
  Vector<int> p_row_index(p_nnz);
  Vector<double> p_value(p_nnz);
  int p_ptr = 0;
    
    
  p_column_start[0] = 0;
    
  // loops of the block columns 
  for (unsigned j = 0; j < Nblock_types; j++)
   {
      
    // determines the block column offset
    int j_block_start = block_start[j];
      
    // loop over the columns of the current block
    for (unsigned k = 0; k < block_dimension(j); k++)
     {
        
      //
      p_column_start[j_block_start + k+1] = p_column_start[j_block_start+k];
        
      // loop over the block rows
      for (unsigned i = 0; i < Nblock_types; i++)
       {
          
        // if block(i,j) pointer no null then
        if (block_matrix_pt(i,j) != 0)
         {
            
          // sets the block row offset
          unsigned i_block_start = block_start[i];			
            
          // creates pointers for the elements in the current block
          int* temp_column_start= block_matrix_pt(i,j)->column_start();
          int* temp_row_index = block_matrix_pt(i,j)->row_index();
          double* temp_value = block_matrix_pt(i,j)->value();
            
          // sets the next column start
          p_column_start[j_block_start + k + 1] += temp_column_start[k+1] -
           temp_column_start[k];
            
          // adds of the current row of the current block to the preconditioner
          // matrix
          for (int l = temp_column_start[k]; l < temp_column_start[k+1]; l++)
           {
            p_row_index[p_ptr] = temp_row_index[l] + i_block_start;
            p_value[p_ptr] = temp_value[l];
            p_ptr++;
           }
         }	
       }			
     }
   }
    
  // creates a new compressed column sparse matrix for the pointer for
  // the current block 
  if (preconditioner_matrix_pt != 0)
   {
    delete preconditioner_matrix_pt;
   }
  preconditioner_matrix_pt = new CCDoubleMatrix(p_value, 
						p_row_index, 
						p_column_start,
						p_nrow,p_nrow);
 }



//=============================================================================
/// Assemble the block preconditioner as a single matrix. 
/// In some cases the block preconditioner cannot be applied in 
/// individual components; this function takes the matrix
/// of block pointers and returns a single matrix containing all the 
/// blocks of the matrix of blocks in a single matrix that can
/// be solved directly. Specialised version for CRDoubleMatrix.
//=============================================================================
 template<>
  void BlockPreconditioner<CRDoubleMatrix>::
  build_preconditioner_matrix(DenseMatrix<CRDoubleMatrix*>& block_matrix_pt,
			      CRDoubleMatrix*& preconditioner_matrix_pt) const
 {
  // the number of blocks
  unsigned n_blocks = this->nblock_types();

#ifdef PARANOID

  // paranoid check that block i is in this block preconditioner
  if (block_matrix_pt.nrow() != n_blocks || block_matrix_pt.ncol() != n_blocks)
   {
    std::ostringstream error_message;
    error_message << "Block_matrix_pt must be of dimension nblock x nblock "
		  << "(" << n_blocks << " x " << n_blocks << "), but it is "
		  << block_matrix_pt.nrow() << " x " << block_matrix_pt.ncol()
		  << std::endl;
    throw OomphLibError(error_message.str(),
                        "BlockPreconditioner::get_block(...)",
                        OOMPH_EXCEPTION_LOCATION);
   }
#endif

  // determine the number of processors
  unsigned nproc = Problem_pt->communicator_pt()->nproc();

   // determine the number of nnzs
   unsigned p_nnz = 0;
  for (unsigned i = 0; i < Nblock_types; i++)
   {
    for (unsigned j = 0; j < Nblock_types; j++)
     {
      if (block_matrix_pt(i,j) !=0)
       {
        p_nnz += block_matrix_pt(i,j)->nnz();
       }
     }
   }

  // construct block offset
  DenseMatrix<unsigned> col_offset(nproc,n_blocks,0);
  unsigned off = 0;
  for (unsigned p = 0; p < nproc; p++)
   {
    for (unsigned b = 0; b < n_blocks; b++)
     {
      col_offset(p,b) = off;
      off += Block_distribution_pt[b]->nrow_local(p);
     }
   }

  //  nrow local
  unsigned p_nrow_local = 
   this->preconditioner_matrix_distribution_pt()->nrow_local();

  // storage for the new matrix
  int* p_row_start = new int[p_nrow_local+1]; 
  int* p_column_index = new int[p_nnz];
  double* p_value = new double[p_nnz];

  // initialise the zero-th entry
  p_row_start[0] = 0;

  // loop over the block rows
  unsigned p_pt = 0;
  unsigned r_pt = 0;
  for (unsigned i = 0; i < Nblock_types; i++)
   {
    
    // loop over the rows of the current block row
    unsigned block_nrow = Block_distribution_pt[i]->nrow_local();
    for (unsigned k = 0; k < block_nrow; k++)
     {

      // initialise p_row_start
      p_row_start[r_pt+1] = p_row_start[r_pt];

      // loop over the block columns
      for (unsigned j = 0; j < Nblock_types; j++)
       {
        
        // if block(i,j) pointer not null then
        if (block_matrix_pt(i,j) != 0)
         {
          
          // creates pointers for the elements in the current block
          int* b_row_start = block_matrix_pt(i,j)->row_start();
          int* b_column_index =block_matrix_pt(i,j)->column_index();
          double* b_value = block_matrix_pt(i,j)->value();
          
          // 
          for (int l = b_row_start[k]; l < b_row_start[k+1]; l++)
           {

            // if b_column_index[l] was a row index, what processor
            // would it be on
            unsigned p = 0;
            int b_first_row = this->block_distribution_pt(j)->first_row(0);
            int b_nrow_local = this->block_distribution_pt(j)->nrow_local(0);
            while (b_column_index[l] < b_first_row || 
                   b_column_index[l] >= b_nrow_local+b_first_row)
             {
              p++;
              b_first_row = this->block_distribution_pt(j)->first_row(p);
              b_nrow_local = this->block_distribution_pt(j)->nrow_local(p);
             }

            // determine the local equation number in the block j/processor p
            // "column block"
            unsigned eqn = b_column_index[l]-b_first_row;

            // add to the preconditioner matrix
            p_value[p_pt] = b_value[l];
            p_column_index[p_pt] = col_offset(p,j)+eqn;
            p_row_start[r_pt+1]++;
            p_pt++;
           }
         }
       }
      
      // increment the row pt
      r_pt++;
     }
   }
 
  // creates a new compressed column sparse matrix for the pointer for
  // the current block 
  if (preconditioner_matrix_pt == 0)
   {
    preconditioner_matrix_pt = 
     new CRDoubleMatrix(this->preconditioner_matrix_distribution_pt());
   }
  unsigned p_nrow = this->preconditioner_matrix_distribution_pt()->nrow();
  preconditioner_matrix_pt->build_without_copy(p_nrow,p_nnz,
                                               p_value,
                                               p_column_index, 
                                               p_row_start);   
 }


}

