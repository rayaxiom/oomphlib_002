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
#include "constrained_navier_stokes_preconditioners.h"


namespace oomph
{


///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////

/*
//======start_of_namespace============================================
/// Namespace for exact solution for pressure advection diffusion 
/// problem
//====================================================================
 namespace PressureAdvectionDiffusionValidation
 {
  

  /// Flag for solution
  unsigned Flag=0;
  
  /// Peclet number -- overwrite with actual Reynolds number
  double Peclet=0.0;
  
  /// Wind
  void wind_function(const Vector<double>& x, Vector<double>& wind)
  {
   if (Flag==0)
    {
     wind[0]=sin(6.0*x[1]);
     wind[1]=cos(6.0*x[0]);
    }
   else
    {
     wind[0]=1.0;
     wind[1]=0.0;
    }
  }
  
  /// Exact solution as a Vector
  void get_exact_u(const Vector<double>& x, Vector<double>& u)
  {
   u.resize(3);
   wind_function(x,u);
   if (Flag==0)
    {
     u[2]=x[0]*x[0]*pow(1.0-x[0],2.0)*x[1]*x[1]*pow(1.0-x[1],2.0);
    }
   else
    {
     u[2]=0.1E1-Peclet*x[0]*(0.1E1-0.5*x[0]);
    }
  }
  
  /// Exact solution as a scalar
  void get_exact_u(const Vector<double>& x, double& u)
  {
   if (Flag==0)
    {
     u=x[0]*x[0]*pow(1.0-x[0],2.0)*x[1]*x[1]*pow(1.0-x[1],2.0);
    }
   else
    {
     u=0.1E1-Peclet*x[0]*(0.1E1-0.5*x[0]);
    }
  }
  
  /// Source function required to make the solution above an exact solution 
  double source_function(const Vector<double>& x_vect)
  {

   double x[2];
   x[0]=x_vect[0];
   x[1]=x_vect[1];

   
   double source=0.0;

   if (Flag==0)
    {
     source=
      Peclet*(sin(0.6E1*x[1])*(2.0*x[0]*pow(1.0-x[0],2.0)*x[1]*x[1]*pow(
                                1.0-x[1],2.0)-2.0*x[0]*x[0]*(1.0-x[0])*x[1]*x[1]*pow(1.0-x[1],2.0))+cos(0.6E1*x
                                                                                                        [0])*(2.0*x[0]*x[0]*pow(1.0-x[0],2.0)*x[1]*pow(1.0-x[1],2.0)-2.0*x[0]*x[0]*pow(
                                                                                                               1.0-x[0],2.0)*x[1]*x[1]*(1.0-x[1])))-2.0*pow(1.0-x[0],2.0)*x[1]*x[1]*pow(1.0-x
                                                                                                                                                                                        [1],2.0)+8.0*x[0]*(1.0-x[0])*x[1]*x[1]*pow(1.0-x[1],2.0)-2.0*x[0]*x[0]*x[1]*x
      [1]*pow(1.0-x[1],2.0)-2.0*x[0]*x[0]*pow(1.0-x[0],2.0)*pow(1.0-x[1],2.0)+8.0*x
      [0]*x[0]*pow(1.0-x[0],2.0)*x[1]*(1.0-x[1])-2.0*x[0]*x[0]*pow(1.0-x[0],2.0)*x[1]
      *x[1];
    }
   else
    {
     source=Peclet*(-0.1E1*Peclet*(0.1E1-0.5*x[0])+0.5*Peclet*x[0])-0.1E1*Peclet
;

    }
   
   return source;
  }
  
  
 } // end of namespace
*/ 
 
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
 

//===========================================================================
/// Setup the least-squares commutator Navier Stokes preconditioner. This
/// extracts blocks corresponding to the velocity and pressure unknowns,
/// creates the matrices actually needed in the application of the
/// preconditioner and deletes what can be deleted... Note that
/// this preconditioner needs a CRDoubleMatrix.
//============================================================================
 void ConstrainedNavierStokesSchurComplementPreconditioner::
 setup(Problem* problem_pt, DoubleMatrixBase* matrix_pt)
 {

  // For debugging...
  bool doc_block_matrices=false;

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // NOTE: In the interest of minimising memory usage, several containers
  //       are recycled, therefore their content/meaning changes
  //       throughout this function. The code is carefully annotated
  //       but you'll have to read it line by line!
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  // make sure any old data is deleted
  clean_up_memory();

#ifdef PARANOID
  // paranoid check that the navier stokes mesh pt has been set
  if (Navier_stokes_mesh_pt == 0)
   {
    std::ostringstream error_message;
    error_message << "The navier stokes elements mesh pointer must be set.\n"
                  << "Use method set_navier_stokes_mesh(...)";
    throw OomphLibError(error_message.str(),
                     	"ConstrainedNavierStokesSchurComplementPreconditioner::setup()",
                        OOMPH_EXCEPTION_LOCATION);
   }
//  if (Prec_blocks.size() == 0)
//   {
//    std::ostringstream error_message;
//    error_message << "set_prec_blocks(...) function has not been called.";
//    throw OomphLibError(error_message.str(),
//                     	"ConstrainedNavierStokesSchurComplementPreconditioner::setup()",
//                        OOMPH_EXCEPTION_LOCATION);
//   }
//  if (Master_doftype_order.size() == 0)
//   {
//    std::ostringstream error_message;
//    error_message << "set_master_doftype_ordering(...) function has not been called.";
//    throw OomphLibError(error_message.str(),
//                     	"ConstrainedNavierStokesSchurComplementPreconditioner::setup()",
//                        OOMPH_EXCEPTION_LOCATION);
//   }
//
#endif

  // set the mesh
  this->set_mesh(0,problem_pt,Navier_stokes_mesh_pt);

  // Get blocks
  // ----------

  // In comes the current Jacobian. Recast it to a CR double matrix;
  // shout if that can't be done.
  CRDoubleMatrix* cr_matrix_pt = dynamic_cast<CRDoubleMatrix*>(matrix_pt);


#ifdef PARANOID
  if (cr_matrix_pt==0)
   {
    std::ostringstream error_message;
    error_message 
     << "ConstrainedNavierStokesSchurComplementPreconditioner only works with "
     << "CRDoubleMatrix matrices" << std::endl;
    throw OomphLibError(error_message.str(),
                     	"ConstrainedNavierStokesSchurComplementPreconditioner::setup()",
                        OOMPH_EXCEPTION_LOCATION);
   }
#endif

  if (doc_block_matrices)
   {
    std::stringstream junk;
    junk << "j_matrix" << problem_pt->communicator_pt()->my_rank()
         << ".dat";
    oomph_info << "About to output " << junk.str() << std::endl;
    cr_matrix_pt->sparse_indexed_output_with_offset(junk.str());
    oomph_info << "Done output of " << junk.str() << std::endl;
   }


  // Set up block look up schemes (done automatically in the
  // BlockPreconditioner base class, based on the information 
  // provided in the block-preconditionable elements in the problem)

  // this preconditioner has two types of block:
  // type 0: velocity - corresponding to DOFs 0 to n-2
  // type 1: pressure - corresponding to DOF n-1
  double t_block_start = TimingHelpers::timer();
  unsigned ndof_types = 0;
  if (this->is_subsidiary_block_preconditioner())
   {
    ndof_types = this->ndof_types();
   }
  else
   {
    ndof_types = this->ndof_types_in_mesh(0);
   }

  Dim = Navier_stokes_mesh_pt->finite_element_pt(0)->dim();

  this->block_setup(problem_pt,matrix_pt);

  for (unsigned i = 0; i < 5; i++) 
  {
    CRDoubleMatrix* temp_mat_pt = 0;
    this->get_block(0,i,cr_matrix_pt,temp_mat_pt);
    unsigned temp_mat_ncol = temp_mat_pt->ncol();
    std::cout << "ncol:::: " << temp_mat_ncol << std::endl; 
  }

  std::cout << "new testsss: \n" << std::endl;
  unsigned pb_nrow = Prec_blocks.nrow();
  unsigned pb_ncol = Prec_blocks.ncol();
  for (unsigned i = 0; i < pb_ncol; i++) 
  {
    unsigned subblock_ncol = Prec_blocks(0,i)->ncol();
    std::cout << "ncol from pb:" << subblock_ncol << std::endl; 
    
  }
  
pause("BUBBLE POP");

  double t_block_finish = TimingHelpers::timer();
  double block_setup_time = t_block_finish - t_block_start;
  if(Doc_time)
   {
    oomph_info << "Time for block_setup(...) [sec]: "
               << block_setup_time << "\n";
   }

  // determine whether the F preconditioner is a block preconditioner (and
  // therefore a subsidiary preconditioner)
  BlockPreconditioner<CRDoubleMatrix>* F_block_preconditioner_pt = 
   dynamic_cast<BlockPreconditioner<CRDoubleMatrix>* >(F_preconditioner_pt);
  F_preconditioner_is_block_preconditioner = true;
  if (F_block_preconditioner_pt == 0)
   {
    F_preconditioner_is_block_preconditioner = false;
   }
  
  // Get B (the divergence block)
  // RAYRAY
  //double t_get_B_start = TimingHelpers::timer();
  
  // RAYHERE
  //CRDoubleMatrix* b_pt = Prec_blocks[1];
  CRDoubleMatrix* b_pt = 0;
  this->get_block(1,0,cr_matrix_pt,b_pt);
  //double t_get_B_finish = TimingHelpers::timer();
  //if(Doc_time)
  // {
  //double get_B_time = t_get_B_finish - t_get_B_start;
  //  oomph_info << "Time to get B [sec]: "
  //             << get_B_time << "\n";
  // }


  if (doc_block_matrices)
   {
    std::stringstream junk;
    junk << "b_matrix" << problem_pt->communicator_pt()->my_rank()
         << ".dat";
    b_pt->sparse_indexed_output_with_offset(junk.str());
    oomph_info << "Done output of " << junk.str() << std::endl;
   }
  
  /////////////////////////////////////////////////////////////////////////////
  // get the inverse velocity and pressure mass matrices
  CRDoubleMatrix* inv_p_mass_pt = 0;
  CRDoubleMatrix* inv_v_mass_pt = 0;

  // We assume that there is only one pressure dof type, the others are 
  // velocity dof type.
  N_velocity_doftypes = ndof_types - 1;
  

  DenseMatrix<CRDoubleMatrix*> inv_vmm_pts(1,N_velocity_doftypes,0);

  // We only need the velocity mass matrix
  bool do_both=false;
     
  double ivmm_assembly_start_t = TimingHelpers::timer();
  for (unsigned vdof_i = 0; vdof_i < N_velocity_doftypes; vdof_i++) 
  {
    //unsigned required_block_number = Master_doftype_order[vdof_i];
    //std::cout << "vodf_i: " << vdof_i << " re. block: "<< required_block_number<< std::endl; 
    
 //    assemble_inv_press_and_veloc_mass_matrix_diagonal
 //     (inv_p_mass_pt, inv_vmm_pts(0,vdof_i), do_both, 
 //      required_block_number,problem_pt);
     assemble_inv_press_and_veloc_mass_matrix_diagonal
      (inv_p_mass_pt, inv_vmm_pts(0,vdof_i), do_both, 
       vdof_i,problem_pt);

  }
   // pause("doumooo"); 
    
    
  double ivmm_assembly_finish_t = TimingHelpers::timer();
  if (Doc_time)
   {
    double 
     ivmm_assembly_time = ivmm_assembly_finish_t - ivmm_assembly_start_t;
    oomph_info << "Time to assemble inverse diagonal velocity and pressure"
               << "mass matrices) [sec]: "
               << ivmm_assembly_time << "\n";
   }
  
  // Now cat the matrices in inv_vmm_pts along the diagonal by doing:
  // Get the total number of rows.
  // loop through all of the rows and get the values. 
  //   Set the row_start and column.
  double t_cat_ivmm_start = TimingHelpers::timer();
  {
  bool distributed = this->master_distribution_pt()->distributed();
  unsigned long total_nrow_global = 0;
  for(unsigned doftype_i = 0; doftype_i < N_velocity_doftypes; doftype_i++)
  {
    total_nrow_global += inv_vmm_pts(0,doftype_i)->nrow();
    //cout << "nrow = " << inv_v_mass_ptrs(0,doftype_i)->nrow() << endl;
  }
  
  //cout << " total_global_nrow = " << total_nrow_global << endl;
  
  LinearAlgebraDistribution* new_distribution_pt 
    = new LinearAlgebraDistribution(problem_pt->communicator_pt(),
                                    total_nrow_global,distributed);

  // no. values = total_nrow_global since:
  // 1) This is not disributed
  // 2) We only extract the diagonals.
  Vector<double> new_values(total_nrow_global);
  Vector<int> new_column_indices(total_nrow_global);
  Vector<int> new_row_start(total_nrow_global + 1);

  // Loop through the inv_v_mm_pts and populate the
  // new_values, new_column_indices and new_row_start
  unsigned long row_i = 0;
  for(unsigned doftype_i = 0; doftype_i < N_velocity_doftypes; doftype_i++)
  {
    CRDoubleMatrix* current_block_pt = inv_vmm_pts(0,doftype_i);
    double* current_block_values = current_block_pt->value();
    unsigned long block_nrow_local = current_block_pt->nrow_local();
    for(unsigned long block_nrow_local_i = 0;
        block_nrow_local_i < block_nrow_local;
        block_nrow_local_i++)
    {
      new_values[row_i] = current_block_values[block_nrow_local_i];
      new_column_indices[row_i] = row_i;
      new_row_start[row_i] = row_i;
      
      row_i++;
    }
  }

  new_row_start[total_nrow_global] = total_nrow_global;
  inv_v_mass_pt = new CRDoubleMatrix(new_distribution_pt);
  // This is a square matrix
  inv_v_mass_pt->build(total_nrow_global,
                       new_values,
                       new_column_indices,
                       new_row_start);
  }
  
  double t_cat_ivmm_finish = TimingHelpers::timer();
  if(Doc_time)
   {
     double t_cat_ivmm_time = t_cat_ivmm_finish - t_cat_ivmm_start;
     cout << "t_cat_ivmm_time: " << t_cat_ivmm_time << "\n";
   } 
  
  // Delete the pointers to the sub matrices.
  for(unsigned mass_i = 0; mass_i < N_velocity_doftypes; mass_i++)
  {
    delete inv_vmm_pts(0,mass_i);
    inv_vmm_pts(0,mass_i) = 0;
  }


  if (doc_block_matrices)
   {
    std::stringstream junk;
    junk << "inv_v_mass_matrix" 
         << problem_pt->communicator_pt()->my_rank()
         << ".dat";
    inv_v_mass_pt->sparse_indexed_output_with_offset(junk.str());
    oomph_info << "Done output of " << junk.str() << std::endl;
   }
  /////////////////////////////////////////////////////////////////////////////
  
  // Get gradient matrix Bt
  // RAYHERE
  //CRDoubleMatrix* bt_pt = Prec_blocks[2];
  CRDoubleMatrix* bt_pt = 0;
  double t_get_Bt_start = TimingHelpers::timer();
  this->get_block(0,1,cr_matrix_pt,bt_pt);
  double t_get_Bt_finish = TimingHelpers::timer();
  if(Doc_time)
   {
    double t_get_Bt_time = t_get_Bt_finish - t_get_Bt_start;
    oomph_info << "Time to get Bt [sec]: "
               << t_get_Bt_time << std::endl;
   }

  if (doc_block_matrices)   
   {
    std::stringstream junk;
    junk << "bt_matrix" << problem_pt->communicator_pt()->my_rank()
         << ".dat";
    bt_pt->sparse_indexed_output_with_offset(junk.str());
    oomph_info << "Done output of " << junk.str() << std::endl;
   }
  
  // Build pressure Poisson matrix 
  CRDoubleMatrix* p_matrix_pt = new CRDoubleMatrix;

  // Multiply inverse velocity mass matrix by gradient matrix B^T
  double t_QBt_matrix_start = TimingHelpers::timer();
  CRDoubleMatrix* qbt_pt = new CRDoubleMatrix;
  inv_v_mass_pt->multiply(*bt_pt, *qbt_pt);
  //delete bt_pt;
  
  //Prec_blocks[2]->sparse_indexed_output("Bt_after_delete");
    
  // Store product in bt_pt 
  bt_pt = qbt_pt;
  double t_QBt_matrix_finish = TimingHelpers::timer();
  if(Doc_time)
   {
    double t_QBt_time = t_QBt_matrix_finish - t_QBt_matrix_start;
    oomph_info << "Time to generate QBt [sec]: "
               << t_QBt_time << std::endl;
   }
  delete inv_v_mass_pt;
  
  
  // Multiply B from left by divergence matrix B and store result in 
  // pressure Poisson matrix.
  double t_p_matrix_start = TimingHelpers::timer();
  b_pt->multiply(*bt_pt, *p_matrix_pt);
  double t_p_matrix_finish = TimingHelpers::timer();
  if(Doc_time)
   {
    double t_p_time = t_p_matrix_finish - t_p_matrix_start;
    oomph_info << "Time to generate P [sec]: "
               << t_p_time << std::endl;
   }
  // Kill divergence matrix because we don't need it any more
  //delete b_pt;
  

  double t_QBt_MV_start = TimingHelpers::timer();
  QBt_mat_vec_pt = new MatrixVectorProduct;
  QBt_mat_vec_pt->setup(bt_pt);
  double t_QBt_MV_finish = TimingHelpers::timer();
  if(Doc_time)
   {
    double t_p_time = t_QBt_MV_finish - t_QBt_MV_start;
    oomph_info << "Time to build QBt matrix vector operator [sec]: "
               << t_p_time << std::endl;
   }
  // Kill gradient matrix B^T (it's been overwritten anyway and
  // needs to be recomputed afresh below)
  //delete bt_pt;
  

  // Do we need the Fp stuff?
  /* 
  if (!Use_LSC)
   {
    // Get pressure advection diffusion matrix Fp and store in 
    // a "big" matrix (same size as the problem's Jacobian)
    double t_get_Fp_start = TimingHelpers::timer();
    CRDoubleMatrix full_fp_matrix;
    get_pressure_advection_diffusion_matrix(full_fp_matrix);
    
    // Now extract the pressure pressure block
    CRDoubleMatrix* fp_matrix_pt=0;
    this->get_block(1,1,&full_fp_matrix,fp_matrix_pt);
    double t_get_Fp_finish = TimingHelpers::timer();
    if(Doc_time)
     {
      double t_get_Fp_time = t_get_Fp_finish - t_get_Fp_start;
      oomph_info << "Time to get Fp [sec]: "
                 << t_get_Fp_time << std::endl;
     }
    
    // Build vector product of pressure advection diffusion matrix with
    // inverse pressure mass matrix
    CRDoubleMatrix* fp_qp_inv_pt = new CRDoubleMatrix; 
    fp_matrix_pt->multiply(*inv_p_mass_pt, *fp_qp_inv_pt);
    
    // Build the matvec operator for E = F_p Q_p^{-1}
    double t_Fp_Qp_inv_MV_start = TimingHelpers::timer();
    E_mat_vec_pt = new MatrixVectorProduct;
    E_mat_vec_pt->setup(fp_qp_inv_pt);
    double t_Fp_Qp_inv_MV_finish = TimingHelpers::timer();
    if(Doc_time)
     {
      double t_p_time = t_Fp_Qp_inv_MV_finish - t_Fp_Qp_inv_MV_start;
      oomph_info << "Time to build Fp Qp^{-1} matrix vector operator [sec]: "
                 << t_p_time << std::endl;
     }
    // Kill pressure advection diffusion and inverse pressure mass matrices
    delete inv_p_mass_pt;
    inv_p_mass_pt=0;
    delete fp_qp_inv_pt;
    fp_qp_inv_pt=0;
   }
   */


  // Get momentum block F
  CRDoubleMatrix* f_pt = 0; //Prec_blocks[0];
  double t_get_F_start = TimingHelpers::timer();
  this->get_block(0,0,cr_matrix_pt,f_pt);
  double t_get_F_finish = TimingHelpers::timer();
  if(Doc_time)
   {
    double t_get_F_time = t_get_F_finish - t_get_F_start;
    oomph_info << "Time to get F [sec]: "
               << t_get_F_time << std::endl;
   }
  
  // form the matrix vector product helper
  double t_F_MV_start = TimingHelpers::timer();
  F_mat_vec_pt = new MatrixVectorProduct;
  F_mat_vec_pt->setup(f_pt);
  double t_F_MV_finish = TimingHelpers::timer();
  if(Doc_time)
   {
    double t_F_MV_time = t_F_MV_finish - t_F_MV_start;
    oomph_info << "Time to build F Matrix Vector Operator [sec]: "
               << t_F_MV_time << std::endl;
   }
  
  // if F is a block preconditioner then we can delete the F matrix
  if (F_preconditioner_is_block_preconditioner)
   {
    delete f_pt;
   }
  
  // Rebuild Bt (remember that we temporarily overwrote
  // it by its product with the inverse velocity mass matrix)
  t_get_Bt_start = TimingHelpers::timer();
  //bt_pt = Prec_blocks[2];
  bt_pt->sparse_indexed_output("Bt_from_constns");
  //pause("Before get bt block"); 
  
  this->get_block(0,1,cr_matrix_pt,bt_pt);
  t_get_Bt_finish = TimingHelpers::timer();
  
  if(Doc_time)
   {
    double t_get_Bt_time = t_get_Bt_finish - t_get_Bt_start;
    oomph_info << "Time to get Bt [sec]: "
               << t_get_Bt_time << std::endl;
   }
 

  // form the matrix vector operator for Bt
  double t_Bt_MV_start = TimingHelpers::timer();
  Bt_mat_vec_pt = new MatrixVectorProduct;

  std::cout << "Bt matrix built: " << bt_pt->built()<< std::endl; 
  
// pause("after get bt block"); 
  Bt_mat_vec_pt->setup(bt_pt);
//pause("after bt mat vec setup"); 
  double t_Bt_MV_finish = TimingHelpers::timer();
  if(Doc_time)
   {
    double t_Bt_MV_time = t_Bt_MV_finish - t_Bt_MV_start;
    oomph_info << "Time to build Bt Matrix Vector Operator [sec]: "
               << t_Bt_MV_time << std::endl;
   }
  //delete bt_pt;

  // if the P preconditioner has not been setup
  if (P_preconditioner_pt == 0)
   {
  //  pause("Setting up SuperLu for P"); 
    P_preconditioner_pt = new SuperLUPreconditioner;
    Using_default_p_preconditioner = true;
   }
 

  // Setup the preconditioner for the Pressure matrix
  double t_p_prec_start = TimingHelpers::timer();

  if (doc_block_matrices)
   {
    std::stringstream junk;
    junk << "p_matrix" << problem_pt->communicator_pt()->my_rank()
         << ".dat";
    p_matrix_pt->sparse_indexed_output_with_offset(junk.str());
    oomph_info << "Done output of " << junk.str() << std::endl;
   }
 //  pause("P_prec setup"); 
   
  P_preconditioner_pt->setup(problem_pt, p_matrix_pt);
  delete p_matrix_pt;
  p_matrix_pt=0;
  double t_p_prec_finish = TimingHelpers::timer();
  if(Doc_time)
   {
    double t_p_prec_time = t_p_prec_finish - t_p_prec_start;
    oomph_info << "P sub-preconditioner setup time [sec]: "
               << t_p_prec_time << "\n";
   }
  
  // Set up solver for solution of system with momentum matrix
  // ----------------------------------------------------------

  // if the F preconditioner has not been setup
  if (F_preconditioner_pt == 0)
   {
//    pause("Creating SuperLU for F_prec"); 
    
    F_preconditioner_pt = new SuperLUPreconditioner;
    Using_default_f_preconditioner = true;
   }

  // if F is a block preconditioner
  double t_f_prec_start = TimingHelpers::timer();
  if (F_preconditioner_is_block_preconditioner)
   {
    unsigned ndof_types = this->ndof_types();
    ndof_types--;
    Vector<unsigned> dof_map(ndof_types);
    for (unsigned i = 0; i < ndof_types; i++)
     {
      dof_map[i] = i;
     }
    F_block_preconditioner_pt->
     turn_into_subsidiary_block_preconditioner(this,dof_map);
    // NOTE: I HAVE TO CHANGE THIS SO IT USES MY AUGMENTED MATRIX
 //   pause("F_prec setup"); 
    
    F_block_preconditioner_pt->setup(problem_pt,matrix_pt);
   }
  // otherwise F is not a block preconditioner
  else
   {
    F_preconditioner_pt->setup(problem_pt,f_pt);
    //delete f_pt;
   }
  double t_f_prec_finish = TimingHelpers::timer();
  if(Doc_time)
   {
    double t_f_prec_time = t_f_prec_finish - t_f_prec_start;
    oomph_info << "F sub-preconditioner setup time [sec]: "
               << t_f_prec_time << "\n";
   }
  
  // Remember that the preconditioner has been setup so
  // the stored information can be wiped when we
  // come here next...
  Preconditioner_has_been_setup = true;
 }



//=======================================================================
 /// Apply preconditioner to r.
//=======================================================================
 void ConstrainedNavierStokesSchurComplementPreconditioner:: 
 preconditioner_solve(const DoubleVector &r, DoubleVector &z)
 {

#ifdef PARANOID
  if (Preconditioner_has_been_setup==false)
   {
    std::ostringstream error_message;
    error_message << "setup must be called before using preconditioner_solve";
    throw OomphLibError(
     error_message.str(),
     "ConstrainedNavierStokesSchurComplementPreconditioner::preconditioner_solve()",
     OOMPH_EXCEPTION_LOCATION);
   }
  if (z.built())
   {
    if (z.nrow() != r.nrow())
     {
      std::ostringstream error_message;
      error_message << "The vectors z and r must have the same number of "
                    << "of global rows";
      throw OomphLibError(
       error_message.str(),
       "ConstrainedNavierStokesSchurComplementPreconditioner::preconditioner_solve()",
       OOMPH_EXCEPTION_LOCATION);      
     }
   }
#endif

  // if z is not setup then give it the same distribution
  if (!z.distribution_pt()->built())
   {
    z.build(r.distribution_pt(),0.0);
   }

  // Step 1 - apply approximate Schur inverse to pressure unknowns (block 1)
  // -----------------------------------------------------------------------

  // Working vectors
  DoubleVector temp_vec;
  DoubleVector another_temp_vec;
  DoubleVector yet_another_temp_vec;

  // Copy pressure values from residual vector to temp_vec:
  // Loop over all entries in the global vector (this one
  // includes velocity and pressure dofs in some random fashion)
  //this->get_block_vector(1,r,temp_vec);
  //this->get_block_vector(Master_doftype_order[N_velocity_doftypes],r,temp_vec);
  this->get_block_vector(N_velocity_doftypes,r,temp_vec);
  
  // NOTE: The vector temp_vec now contains the vector r_p.


  // LSC version
  if (Use_LSC)
   {
    
    // Solve first pressure Poisson system
#ifdef PARANOID
    // check a solver has been set
    if (P_preconditioner_pt==0)
     {
      std::ostringstream error_message;
      error_message << "P_preconditioner_pt has not been set.";
      throw OomphLibError(
       error_message.str(),
       "ConstrainedNavierStokesSchurComplementPreconditioner::preconditioner_solve()",
       OOMPH_EXCEPTION_LOCATION);
     }
#endif
    
    // use some Preconditioner's preconditioner_solve function
    P_preconditioner_pt->preconditioner_solve(temp_vec, another_temp_vec);
    
    // NOTE: The vector another_temp_vec now contains the vector P^{-1} r_p
    
    // Multiply another_temp_vec by matrix E and stick the result into temp_vec
    temp_vec.clear();  
    QBt_mat_vec_pt->multiply(another_temp_vec, temp_vec);
    //pause("Le multiply"); 
    another_temp_vec.clear();
    F_mat_vec_pt->multiply(temp_vec,another_temp_vec);
    //pause("Le duex multiply"); 
    temp_vec.clear();
    QBt_mat_vec_pt->multiply_transpose(another_temp_vec, temp_vec);
    //pause("Le Tois multiply"); 
    
    
    // NOTE: The vector temp_vec now contains E P^{-1} r_p
    
    // Solve second pressure Poisson system using preconditioner_solve
    another_temp_vec.clear();
    P_preconditioner_pt->preconditioner_solve(temp_vec, another_temp_vec);
    
    // NOTE: The vector another_temp_vec now contains z_p = P^{-1} E P^{-1} r_p
    //       as required (apart from the sign which we'll fix in the
    //       next step.

   }
  // Fp version
/*else
{
      
    // Multiply temp_vec by matrix E and stick the result into 
    // yet_another_temp_vec
    E_mat_vec_pt->multiply(temp_vec,yet_another_temp_vec);
        
    // NOTE: The vector yet_another_temp_vec now contains Fp Qp^{-1} r_p

    // Solve pressure Poisson system
#ifdef PARANOID
    // check a solver has been set
    if (P_preconditioner_pt==0)
     {
      std::ostringstream error_message;
      error_message << "P_preconditioner_pt has not been set.";
      throw OomphLibError(
       error_message.str(),
       "ConstrainedNavierStokesSchurComplementPreconditioner::preconditioner_solve()",
       OOMPH_EXCEPTION_LOCATION);
     }
#endif
    
    // Solve second pressure Poisson system using preconditioner_solve
    another_temp_vec.clear();
    P_preconditioner_pt->preconditioner_solve(yet_another_temp_vec, 
                                              another_temp_vec);
    
    // NOTE: The vector another_temp_vec now contains 
    //       z_p = P^{-1} Fp Qp^{-1} r_p
    //       as required (apart from the sign which we'll fix in the
    //       next step.

   }*/

  // Now copy another_temp_vec (i.e. z_p) back into the global vector z.
  // Loop over all entries in the global results vector z:
  temp_vec.build(another_temp_vec.distribution_pt(),0.0);
  temp_vec -= another_temp_vec;
//  this->return_block_vector(Master_doftype_order[N_velocity_doftypes],temp_vec,z);
  this->return_block_vector(N_velocity_doftypes,temp_vec,z);
  

    
  // Step 2 - apply preconditioner to velocity unknowns (block 0)
  // ------------------------------------------------------------
  
  // Recall that another_temp_vec (computed above) contains the
  // negative of the solution of the Schur complement systen, -z_p.
  // Multiply by G (stored in Block_matrix_pt(0,1) and store
  // result in temp_vec (vector resizes itself).
  temp_vec.clear();
  //pause("Before 4th damn multiply"); 
  
  Bt_mat_vec_pt->multiply(another_temp_vec, temp_vec);
//pause("The 4th damn multiply"); 


  // NOTE: temp_vec now contains -G z_p

  // The vector another_temp_vec is no longer needed -- re-use it to store
  // velocity quantities:
  another_temp_vec.clear();

  // Concat the velocity RHS:
  unsigned long v_nrow = 0;
  for (unsigned i = 0; i < N_velocity_doftypes; i++) 
  {
//    this->get_block_vector(Master_doftype_order[i],r,another_temp_vec);
    this->get_block_vector(i,r,another_temp_vec);
    v_nrow += another_temp_vec.nrow();
    another_temp_vec.clear();
  }
  
  bool distributed = this->master_distribution_pt()->distributed();
  LinearAlgebraDistribution* new_distribution_pt
    = new LinearAlgebraDistribution(problem_pt()->communicator_pt(),
                                    v_nrow,distributed);

  yet_another_temp_vec.build(new_distribution_pt,0.0);

  unsigned long v_nrow_i = 0;
  for (unsigned i = 0; i < N_velocity_doftypes; i++) 
  {
//    this->get_block_vector(Master_doftype_order[i],r,another_temp_vec);
    this->get_block_vector(i,r,another_temp_vec);
    unsigned long current_block_nrow = another_temp_vec.nrow();

    for (unsigned current_block_i = 0; current_block_i < current_block_nrow; 
         current_block_i++) 
    {
      yet_another_temp_vec[v_nrow_i] = another_temp_vec[current_block_i];
      v_nrow_i++;
    }
    another_temp_vec.clear();
  }
  
  

  // Loop over all enries in the global vector and find the
  // entries associated with the velocities:
  //get_block_vector(0,r,another_temp_vec);
  //another_temp_vec += temp_vec;
 
  //yet_another_temp_vec contained the merged blocks in the order we want.
  yet_another_temp_vec += temp_vec;
  temp_vec.clear();


  // NOTE:  The vector another_temp_vec now contains r_u - G z_p

  // Solve momentum system
#ifdef PARANOID
  // check a solver has been set
  if (F_preconditioner_pt==0)
   {
    std::ostringstream error_message;
    error_message << "F_preconditioner_pt has not been set."; 
    throw OomphLibError(
     error_message.str(),
     "ConstrainedNavierStokesSchurComplementPreconditioner::preconditioner_solve()",
     OOMPH_EXCEPTION_LOCATION);
   }
#endif

  // use some Preconditioner's preconditioner solve
  // and return
  if (F_preconditioner_is_block_preconditioner)
   {
    pause("should not get here.... from NS prec solve...");
    return_block_vector(0,yet_another_temp_vec,z);
    F_preconditioner_pt->preconditioner_solve(z,z);
   }
  else
   {
    F_preconditioner_pt->preconditioner_solve(yet_another_temp_vec, temp_vec);

    unsigned long global_row_i = 0;
    for(unsigned i = 0; i< N_velocity_doftypes; i++)
    {
      another_temp_vec.clear();
     
//      this->get_block_vector(Master_doftype_order[i],r,another_temp_vec);
      this->get_block_vector(i,r,another_temp_vec);
     
      unsigned long current_block_nrow = another_temp_vec.nrow();
      
      for(unsigned long current_block_nrow_i = 0; 
          current_block_nrow_i < current_block_nrow;
          current_block_nrow_i++)
      {
        another_temp_vec[current_block_nrow_i] = temp_vec[global_row_i];
        global_row_i++;
      }
     
//      this->return_block_vector(Master_doftype_order[i],another_temp_vec,z);
      this->return_block_vector(i,another_temp_vec,z);
    } //  loop though the block vectors

    //this->return_block_vector(Master_doftype_order[i],another_temp_vec,z);
   }
 }




//========================================================================
/// Helper function to assemble the diagonal of the velocity
/// mass matrix from the elemental contributions defined in
/// NavierStokesEquations<DIM>::get_mass_matrix_diagonal(...).
/// If do_both=true, both are computed, otherwise only the velocity
/// mass matrix (the LSC version of the preconditioner only needs
/// that one)
//========================================================================
// void ConstrainedNavierStokesSchurComplementPreconditioner:: 
// assemble_inv_press_and_veloc_mass_matrix_diagonal(
//  CRDoubleMatrix*& inv_p_mass_pt,
//  CRDoubleMatrix*& inv_v_mass_pt,
//  const bool& do_both)
// {
//
//  // determine the velocity rows required by this processor
//  unsigned v_first_row = this->block_distribution_pt(0)->first_row();
//  unsigned v_nrow_local = this->block_distribution_pt(0)->nrow_local();
//  unsigned v_nrow = this->block_distribution_pt(0)->nrow();
//  
//  // create storage for the diagonals
//  double* v_values = new double[v_nrow_local];
//  for (unsigned i = 0; i < v_nrow_local; i++)
//   {
//    v_values[i] = 0.0;
//   }
//
//  // Equivalent information for pressure mass matrix (only needed for 
//  // Fp version)
//  unsigned p_first_row=0;
//  unsigned p_nrow_local=0;
//  unsigned p_nrow=0;
//  double* p_values = 0;
//  if (!Use_LSC)
//   {
//    // determine the pressure rows required by this processor
//    p_first_row = this->block_distribution_pt(1)->first_row();
//    p_nrow_local = this->block_distribution_pt(1)->nrow_local();
//    p_nrow = this->block_distribution_pt(1)->nrow();
//  
//    // create storage for the diagonals
//    p_values = new double[p_nrow_local];
//    for (unsigned i = 0; i < p_nrow_local; i++)
//     {
//      p_values[i] = 0.0;
//     }
//   }
//
//  // store the problem pt
//  const Problem* problem_pt = this->problem_pt();
//
//  // if the problem is distributed
//  bool distributed = false;
//#ifdef OOMPH_HAS_MPI
//  if (problem_pt->distributed() || 
//      this->master_distribution_pt()->distributed())
//   {
//    distributed = true;
//   }
//#endif
//
//  // next we get the diagonal velocity mass matrix data
//  if (distributed)
//   {
//
//#ifdef OOMPH_HAS_MPI
//
//    // the number of processors
//    unsigned nproc = problem_pt->communicator_pt()->nproc();
//
//    // and my rank
//    unsigned my_rank = problem_pt->communicator_pt()->my_rank();
//
//    // determine the rows for which we have lookup rows
//
//    // if the problem is NOT distributed then we only classify global equations
//    // on this processor to avoid duplication (as every processor holds 
//    // every element)
//    unsigned first_lookup_row = 0; 
//    unsigned last_lookup_row = 0;
//    first_lookup_row = this->master_distribution_pt()->first_row();
//    last_lookup_row = first_lookup_row + 
//     this->master_distribution_pt()->nrow_local() - 1;
//
//    // find number of local elements
//    unsigned n_el = Navier_stokes_mesh_pt->nelement();
//    
//    // get the master distribution pt
//    const LinearAlgebraDistribution* master_distribution_pt = 
//     this->master_distribution_pt();
//
//    // Do the two blocks (0: veloc; 1: press)
//    unsigned max_block=0;
//    if (!Use_LSC) max_block=1;
//    for (unsigned block_index=0;block_index<=max_block;block_index++)
//     {
//      
//      // Local working variables: Default to velocity
//      unsigned v_or_p_first_row=v_first_row;
//      double* v_or_p_values=v_values;
//      // Switch to pressure
//      if (block_index==1)
//       {
//        v_or_p_first_row=p_first_row;
//        v_or_p_values=p_values;
//       }
//    
//
//      // the diagonal mass matrix contributions that have been
//      // classified and should be sent to another processor
//      Vector<double>* classified_contributions_send 
//       = new Vector<double>[nproc];
//      
//      // the corresponding block indices
//      Vector<unsigned>* classified_indices_send
//       = new Vector<unsigned>[nproc];
//      
//      // the matrix contributions that cannot be classified by this processor
//      // and therefore must be sent to another for classification
//      Vector<double>* unclassified_contributions_send
//       = new Vector<double>[nproc];
//      
//      // the corresponding global indices that require classification
//      Vector<unsigned>* unclassified_indices_send
//       = new Vector<unsigned>[nproc];
//      
//      // get the velocity or pressure distribution pt
//      const LinearAlgebraDistribution* velocity_or_press_dist_pt 
//       = this->block_distribution_pt(block_index);
//      
//      // get the contribution for each element
//      for (unsigned e = 0; e < n_el; e++)
//       {
//        
//        // Get element
//        GeneralisedElement* el_pt=Navier_stokes_mesh_pt->element_pt(e);
//
//        // check that the element is not halo
//        if (!el_pt->is_halo())
//         {
//          
//          // find number of degrees of freedom in the element
//          // (this is slightly too big because it includes the
//          // pressure dofs but this doesn't matter)
//          unsigned el_dof = el_pt->ndof();
//          
//          // Allocate local storage for the element's contribution to the
//          // mass matrix diagonal
//          Vector<double> el_vmm_diagonal(el_dof,0.0);
//          Vector<double> el_pmm_diagonal(el_dof,0.0);
//          
//          unsigned which_one=2;
//          if (block_index==1) which_one=1;
//
//          NavierStokesElementWithDiagonalMassMatrices* cast_el_pt=0;
//          cast_el_pt=dynamic_cast<NavierStokesElementWithDiagonalMassMatrices*>
//           (el_pt);
//          if (cast_el_pt!=0)
//           {
//            cast_el_pt->get_pressure_and_velocity_mass_matrix_diagonal( 
//             el_pmm_diagonal,el_vmm_diagonal,which_one);
//           }
//
//          // get the contribution for each dof
//          for (unsigned i = 0; i < el_dof; i++)
//           {
//            
//            //Get the equation number
//            unsigned eqn_number = el_pt->eqn_number(i);
//            
//            // if I have lookup information on this processor
//            if ((eqn_number >= first_lookup_row) && 
//                (eqn_number <= last_lookup_row)    )
//             {
//              
//              // Only use the dofs that we're dealing with here
//              if ( this->block_number(eqn_number)==int(block_index) )
//               {
//                
//                // get the index in the block
//                unsigned index = this->index_in_block(eqn_number);
//                
//                // determine which processor requires the block index
//                for (unsigned p = 0; p < nproc; p++)
//                 {
//                  if ( (index >= velocity_or_press_dist_pt->first_row(p)) &&
//                       (index < (velocity_or_press_dist_pt->first_row(p)
//                                 +velocity_or_press_dist_pt->nrow_local(p)) ) )
//                   {
//                    
//                    // if it is required by this processor then add the 
//                    // contribution
//                    if (p == my_rank)
//                     {
//                      if (block_index==0)
//                       {
//                        v_or_p_values[index-v_or_p_first_row] 
//                         += el_vmm_diagonal[i];
//                       }
//                      else if (block_index==1)
//                       {
//                        v_or_p_values[index-v_or_p_first_row] 
//                         += el_pmm_diagonal[i];
//                       }
//                     }
//                    // otherwise store it for communication
//                    else
//                     {
//                      if (block_index==0)
//                       {
//                        classified_contributions_send[p]
//                         .push_back(el_vmm_diagonal[i]);
//                        classified_indices_send[p].push_back(index);
//                       }
//                      else if (block_index==1)
//                       {
//                        classified_contributions_send[p]
//                         .push_back(el_pmm_diagonal[i]);
//                        classified_indices_send[p].push_back(index);
//                       }
//                     }
//                   }
//                 }
//               }
//             }
//            // if we do not have the lookup information on this processor
//            // then we send the mass matrix contribution to a processor
//            // which we know to have the lookup information
//            // the assumption: the processor for which the master block
//            // preconditioner distribution will definitely hold the lookup
//            // data for eqn_number (although others may)
//            else if (problem_pt->distributed())
//             {
//              
//              // determine which processor requires the block index
//              unsigned p = 0;
//              while (!(eqn_number >=master_distribution_pt->first_row(p) &&
//                       (eqn_number<(master_distribution_pt->first_row(p)
//                                    +master_distribution_pt->nrow_local(p)))))
//               {
//                p++;
//               }
//              
//              // store the data
//              if (block_index==0)
//               {
//                unclassified_contributions_send[p]
//                 .push_back(el_vmm_diagonal[i]);
//                unclassified_indices_send[p].push_back(eqn_number);
//               }
//              else if (block_index==1)
//               {
//                unclassified_contributions_send[p]
//                 .push_back(el_pmm_diagonal[i]);
//                unclassified_indices_send[p].push_back(eqn_number);
//               }
//              
//             }
//           }
//         }
//       }
//      
//      //next the unclassified contributions are communicated to 
//      //processors that can classify them
//      
//      //first determine how many unclassified rows are to be sent to
//      //each processor
//      unsigned* n_unclassified_send = new unsigned[nproc];
//      for (unsigned p = 0; p < nproc; p++)
//       {
//        if (p == my_rank)
//         {
//          n_unclassified_send[p] = 0;
//         }
//        else
//         {
//          n_unclassified_send[p] 
//           = unclassified_contributions_send[p].size();
//         }
//       }
//      
//      //then all-to-all com number of unclassified to be sent / recv
//      unsigned* n_unclassified_recv = new unsigned[nproc];
//      MPI_Alltoall(n_unclassified_send,1,MPI_UNSIGNED,
//                   n_unclassified_recv,1,MPI_UNSIGNED,
//                   problem_pt->communicator_pt()->mpi_comm());
//      
//      //the base displacement for the sends
//      MPI_Aint base_displacement;
//      MPI_Address(v_or_p_values,&base_displacement); 
//      
//      //allocate storage for the data to be received
//      //and post the sends and recvs
//      Vector<double*> unclassified_contributions_recv(nproc);
//      Vector<unsigned*> unclassified_indices_recv(nproc);
//      Vector<MPI_Request> unclassified_recv_requests;
//      Vector<MPI_Request> unclassified_send_requests;
//      Vector<unsigned> unclassified_recv_proc;
//      for (unsigned p = 0; p < nproc; p++)
//       {
//        if (p != my_rank)
//         {
//          //recv
//          if (n_unclassified_recv[p] > 0)
//           {
//            unclassified_contributions_recv[p] 
//             = new double[n_unclassified_recv[p]];
//            unclassified_indices_recv[p] = new 
//             unsigned[n_unclassified_recv[p]];
//              
//            //data for the struct data type
//            MPI_Datatype recv_types[2];
//            MPI_Aint recv_displacements[2];
//            int recv_sz[2];
//              
//            //contributions
//            MPI_Type_contiguous(n_unclassified_recv[p],MPI_DOUBLE,
//                                &recv_types[0]);
//            MPI_Type_commit(&recv_types[0]);
//            MPI_Address(unclassified_contributions_recv[p],
//                        &recv_displacements[0]);
//            recv_displacements[0] -= base_displacement;
//            recv_sz[0] = 1;
//              
//            //indices
//            MPI_Type_contiguous(n_unclassified_recv[p],MPI_UNSIGNED,
//                                &recv_types[1]);
//            MPI_Type_commit(&recv_types[1]);
//            MPI_Address(unclassified_indices_recv[p],
//                        &recv_displacements[1]);
//            recv_displacements[1] -= base_displacement;
//            recv_sz[1] = 1;
//              
//            //build the final recv type
//            MPI_Datatype final_recv_type;
//            MPI_Type_struct(2,recv_sz,recv_displacements,recv_types,
//                            &final_recv_type);
//            MPI_Type_commit(&final_recv_type);
//              
//            //and recv
//            MPI_Request req;
//            MPI_Irecv(v_or_p_values,1,final_recv_type,p,0,
//                      problem_pt->communicator_pt()->mpi_comm(),&req);
//            unclassified_recv_requests.push_back(req);
//            unclassified_recv_proc.push_back(p); 
//            MPI_Type_free(&recv_types[0]);
//            MPI_Type_free(&recv_types[1]);
//            MPI_Type_free(&final_recv_type);
//           }
//            
//          //send
//          if (n_unclassified_send[p] > 0)
//           {
//            //data for the struct data type
//            MPI_Datatype send_types[2];
//            MPI_Aint send_displacements[2];
//            int send_sz[2];
//              
//            //contributions
//            MPI_Type_contiguous(n_unclassified_send[p],MPI_DOUBLE,
//                                &send_types[0]);
//            MPI_Type_commit(&send_types[0]);
//            MPI_Address(&unclassified_contributions_send[p][0],
//                        &send_displacements[0]);
//            send_displacements[0] -= base_displacement;
//            send_sz[0] = 1;
//              
//            //indices
//            MPI_Type_contiguous(n_unclassified_send[p],MPI_UNSIGNED,
//                                &send_types[1]);
//            MPI_Type_commit(&send_types[1]);
//            MPI_Address(&unclassified_indices_send[p][0],
//                        &send_displacements[1]);
//            send_displacements[1] -= base_displacement;
//            send_sz[1] = 1;
//              
//            //build the final send type
//            MPI_Datatype final_send_type;
//            MPI_Type_struct(2,send_sz,send_displacements,send_types,
//                            &final_send_type);
//            MPI_Type_commit(&final_send_type);
//              
//            //and send
//            MPI_Request req;
//            MPI_Isend(v_or_p_values,1,final_send_type,p,0,
//                      problem_pt->communicator_pt()->mpi_comm(),&req);
//            unclassified_send_requests.push_back(req);
//            MPI_Type_free(&send_types[0]);
//            MPI_Type_free(&send_types[1]);
//            MPI_Type_free(&final_send_type);
//           }
//         }
//       }
//        
//      //next classify the data as it is received
//      unsigned n_unclassified_recv_req = unclassified_recv_requests.size();
//      while (n_unclassified_recv_req > 0)
//       {
//        //get the processor number and remove the completed request
//        //for the vector of requests
//        int req_num;
//        MPI_Waitany(n_unclassified_recv_req,&unclassified_recv_requests[0],
//                    &req_num,MPI_STATUS_IGNORE);
//        unsigned p = unclassified_recv_proc[req_num];
//        unclassified_recv_requests.erase(unclassified_recv_requests.begin()
//                                         +req_num);    
//        unclassified_recv_proc.erase(unclassified_recv_proc.begin()+req_num);
//        n_unclassified_recv_req--;
//          
//        //next classify the dofs 
//        //and store them for sending to other processors if required
//        unsigned n_recv = n_unclassified_recv[p];
//        for (unsigned i = 0; i < n_recv; i++)
//         {
//          unsigned eqn_number = unclassified_indices_recv[p][i];
//          //Only deal with our block unknowns
//          if ( this->block_number(eqn_number)==int(block_index) )
//           {
//              
//            //get the index in the block
//            unsigned index = this->index_in_block(eqn_number);
//              
//            //determine which processor requires the block index
//            for (unsigned pp = 0; pp < nproc; pp++)
//             {
//                
//                
//              if ( (index >= velocity_or_press_dist_pt->first_row(pp)) && 
//                   (index < (velocity_or_press_dist_pt->first_row(pp)          
//                             +velocity_or_press_dist_pt->nrow_local(pp)) ) )
//               {
//                
//                //if it is required by this processor then add the 
//                //contribution
//                if (pp == my_rank)
//                 {
//                  v_or_p_values[index-v_or_p_first_row] 
//                   += unclassified_contributions_recv[p][i];
//                 }
//                //otherwise store it for communication
//                else
//                 {
//                  double v = unclassified_contributions_recv[p][i];
//                  classified_contributions_send[pp].push_back(v);
//                  classified_indices_send[pp].push_back(index);
//                 }
//               }
//             }
//           }
//         }
//          
//        //clean up
//        delete[] unclassified_contributions_recv[p];
//        delete[] unclassified_indices_recv[p];
//       }
//      delete[] n_unclassified_recv;
//        
//      //now all indices have been classified
//        
//      //next the classified contributions are communicated to 
//      //processors that require them
//        
//      //first determine how many classified rows are to be sent to
//      //each processor
//      unsigned* n_classified_send = new unsigned[nproc];
//      for (unsigned p = 0; p < nproc; p++)
//       {
//        if (p == my_rank)
//         {
//          n_classified_send[p] = 0;
//         }
//        else
//         {
//          n_classified_send[p] 
//           = classified_contributions_send[p].size();
//         }
//       }
//        
//      //then all-to-all number of classified to be sent / recv
//      unsigned* n_classified_recv = new unsigned[nproc];
//      MPI_Alltoall(n_classified_send,1,MPI_UNSIGNED,
//                   n_classified_recv,1,MPI_UNSIGNED,
//                   problem_pt->communicator_pt()->mpi_comm());
//        
//      //allocate storage for the data to be received
//      //and post the sends and recvs
//      Vector<double*> classified_contributions_recv(nproc);
//      Vector<unsigned*> classified_indices_recv(nproc);
//      Vector<MPI_Request> classified_recv_requests;
//      Vector<MPI_Request> classified_send_requests;
//      Vector<unsigned> classified_recv_proc;
//      for (unsigned p = 0; p < nproc; p++)
//       {
//        if (p != my_rank)
//         {
//          //recv
//          if (n_classified_recv[p] > 0)
//           {
//            classified_contributions_recv[p] 
//             = new double[n_classified_recv[p]];
//            classified_indices_recv[p] = new unsigned[n_classified_recv[p]];
//              
//            //data for the struct data type
//            MPI_Datatype recv_types[2];
//            MPI_Aint recv_displacements[2];
//            int recv_sz[2];
//              
//            //contributions
//            MPI_Type_contiguous(n_classified_recv[p],MPI_DOUBLE,
//                                &recv_types[0]);
//            MPI_Type_commit(&recv_types[0]);
//            MPI_Address(classified_contributions_recv[p],
//                        &recv_displacements[0]);
//            recv_displacements[0] -= base_displacement;
//            recv_sz[0] = 1;
//              
//            //indices
//            MPI_Type_contiguous(n_classified_recv[p],MPI_UNSIGNED,
//                                &recv_types[1]);
//            MPI_Type_commit(&recv_types[1]);
//            MPI_Address(classified_indices_recv[p],
//                        &recv_displacements[1]);
//            recv_displacements[1] -= base_displacement;
//            recv_sz[1] = 1;
//              
//            //build the final recv type
//            MPI_Datatype final_recv_type;
//            MPI_Type_struct(2,recv_sz,recv_displacements,recv_types,
//                            &final_recv_type);
//            MPI_Type_commit(&final_recv_type);
//              
//            //and recv
//            MPI_Request req;
//            MPI_Irecv(v_or_p_values,1,final_recv_type,p,0,
//                      problem_pt->communicator_pt()->mpi_comm(),&req);
//            classified_recv_requests.push_back(req);
//            classified_recv_proc.push_back(p);
//            MPI_Type_free(&recv_types[0]);
//            MPI_Type_free(&recv_types[1]);
//            MPI_Type_free(&final_recv_type);
//           }
//            
//          //send
//          if (n_classified_send[p] > 0)
//           {
//            //data for the struct data type
//            MPI_Datatype send_types[2];
//            MPI_Aint send_displacements[2];
//            int send_sz[2];
//              
//            //contributions
//            MPI_Type_contiguous(n_classified_send[p],MPI_DOUBLE,
//                                &send_types[0]);
//            MPI_Type_commit(&send_types[0]);
//            MPI_Address(&classified_contributions_send[p][0],
//                        &send_displacements[0]);
//            send_displacements[0] -= base_displacement;
//            send_sz[0] = 1;
//              
//            //indices
//            MPI_Type_contiguous(n_classified_send[p],MPI_UNSIGNED,
//                                &send_types[1]);
//            MPI_Type_commit(&send_types[1]);
//            MPI_Address(&classified_indices_send[p][0],
//                        &send_displacements[1]);
//            send_displacements[1] -= base_displacement;
//            send_sz[1] = 1;
//              
//            //build the final send type
//            MPI_Datatype final_send_type;
//            MPI_Type_struct(2,send_sz,send_displacements,send_types,
//                            &final_send_type);
//            MPI_Type_commit(&final_send_type);
//              
//            //and send
//            MPI_Request req;
//            MPI_Isend(v_or_p_values,1,final_send_type,p,0,
//                      problem_pt->communicator_pt()->mpi_comm(),&req);
//            classified_send_requests.push_back(req);
//            MPI_Type_free(&send_types[0]);
//            MPI_Type_free(&send_types[1]);
//            MPI_Type_free(&final_send_type);
//           }
//         }
//       }
//        
//      //next classify the data as it is received
//      unsigned n_classified_recv_req = classified_recv_requests.size();
//      while (n_classified_recv_req > 0)
//       {
//        //get the processor number and remove the completed request
//        //for the vector of requests
//        int req_num;
//        MPI_Waitany(n_classified_recv_req,&classified_recv_requests[0],
//                    &req_num,MPI_STATUS_IGNORE);
//        unsigned p = classified_recv_proc[req_num];
//        classified_recv_requests.erase(classified_recv_requests.begin()
//                                       +req_num);    
//        classified_recv_proc.erase(classified_recv_proc.begin()+req_num);
//        n_classified_recv_req--;
//          
//        //next classify the dofs 
//        //and store them for sending to other processors if required
//        unsigned n_recv = n_classified_recv[p];
//        for (unsigned i = 0; i < n_recv; i++)
//         {
//          v_or_p_values[classified_indices_recv[p][i]-v_or_p_first_row] 
//           += classified_contributions_recv[p][i];
//         }
//          
//        //clean up
//        delete[] classified_contributions_recv[p];
//        delete[] classified_indices_recv[p];
//       }
//        
//      //wait for the unclassified sends to complete
//      unsigned n_unclassified_send_req = unclassified_send_requests.size();
//      if (n_unclassified_send_req > 0)
//       {
//        MPI_Waitall(n_unclassified_send_req,&unclassified_send_requests[0],
//                    MPI_STATUS_IGNORE);
//       }
//      delete[] unclassified_contributions_send;
//      delete[] unclassified_indices_send;
//      delete[] n_unclassified_send;
//        
//      //wait for the classified sends to complete
//      unsigned n_classified_send_req = classified_send_requests.size();
//      if (n_classified_send_req > 0) 
//       {
//        MPI_Waitall(n_classified_send_req,&classified_send_requests[0],
//                    MPI_STATUS_IGNORE);
//       }
//      delete[] classified_indices_send;
//      delete[] classified_contributions_send;
//      delete[] n_classified_recv;
//      delete[] n_classified_send;
//
//      // Copy the values back where they belong
//      if (block_index==0)
//       {
//        v_values=v_or_p_values;
//       }
//      else if (block_index==1)
//       {
//        p_values=v_or_p_values;
//       }
//        
//     }
//    
//#endif
//    
//   }
//  // or if the problem is not distributed
//  else
//   {
//
//    // find number of elements
//    unsigned n_el = Navier_stokes_mesh_pt->nelement();
//    
//    // Fp needs pressure and velocity mass matrices
//    unsigned which_one=0;
//    if (Use_LSC) which_one=2;
//
//    // get the contribution for each element
//    for (unsigned e = 0; e < n_el; e++)
//     {
//
//      // Get element
//      GeneralisedElement* el_pt=Navier_stokes_mesh_pt->element_pt(e);
//      
//      // find number of degrees of freedom in the element
//      // (this is slightly too big because it includes the
//      // pressure dofs but this doesn't matter)
//      unsigned el_dof = el_pt->ndof();
//      
//      // allocate local storage for the element's contribution to the
//      // pressure and velocity mass matrix diagonal
//      Vector<double> el_vmm_diagonal(el_dof,0.0);
//      Vector<double> el_pmm_diagonal(el_dof,0.0);
//      
//      NavierStokesElementWithDiagonalMassMatrices* cast_el_pt=0;
//      cast_el_pt=dynamic_cast<NavierStokesElementWithDiagonalMassMatrices*>(
//       el_pt);
//      if (cast_el_pt!=0)
//       {
//        cast_el_pt->get_pressure_and_velocity_mass_matrix_diagonal( 
//         el_pmm_diagonal,el_vmm_diagonal,which_one);
//       }
//      else
//       {
//#ifdef PARANOID
//        std::ostringstream error_message;
//        error_message 
//         << "Navier-Stokes mesh contains element that is not of type \n"
//         << "NavierStokesElementWithDiagonalMassMatrices. \n"
//         << "The element is in fact of type " 
//         << typeid(*el_pt).name() 
//         << "\nWe'll assume that it does not make a used contribution \n" 
//         << "to the inverse diagonal mass matrix used in the preconditioner\n" 
//         << "and (to avoid divisions by zero) fill in dummy unit entries.\n"
//         << "[This case currently arises with flux control problems\n"
//         << "where for simplicity the NetFluxControlElement has been added \n"
//         << "to the Navier Stokes mesh -- this should probably be changed at\n"
//         << "some point -- if you get this warning in any other context\n"
//         << "you should check your code very carefully]\n";
//        OomphLibWarning(
//         error_message.str(),
//         "ConstrainedNavierStokesSchurComplementPreconditioner::assemble_inv_press_and_veloc_mass_matrix_diagonal()",
//         OOMPH_EXCEPTION_LOCATION);       
//#endif 
//
//        // Fill in dummy entries to stop division by zero below
//        for (unsigned j=0;j<el_dof;j++)
//         {
//          el_vmm_diagonal[j]=1.0;
//          el_pmm_diagonal[j]=1.0;
//         }
//       }
//      
//      // Get the contribution for each dof
//      for (unsigned i = 0; i < el_dof; i++)
//       {
//        //Get the equation number
//        unsigned eqn_number = el_pt->eqn_number(i);
//        
//        // Get the velocity dofs
//        if (this->block_number(eqn_number)==0)
//         {
//          // get the index in the block
//          unsigned index = this->index_in_block(eqn_number);
//          
//          // if it is required on this processor
//          if ((index >= v_first_row) &&
//              (index < (v_first_row + v_nrow_local) ) )
//           {
//            v_values[index-v_first_row] += el_vmm_diagonal[i];
//           }
//         }
//        // Get the pressure dofs
//        else if (this->block_number(eqn_number)==1)
//         {
//          if (!Use_LSC)
//           {
//            // get the index in the block
//            unsigned index = this->index_in_block(eqn_number);
//            
//            // if it is required on this processor
//            if ((index >= p_first_row)&&
//                (index < (p_first_row + p_nrow_local)) )
//             {
//              p_values[index-p_first_row] += el_pmm_diagonal[i];
//             }
//           }
//         }
//       } 
//     }
//   }
//  
//  // Create column index and row start for velocity mass matrix
//  int* v_column_index = new int[v_nrow_local]; 
//  int* v_row_start = new int[v_nrow_local+1];
//  for (unsigned i = 0; i < v_nrow_local; i++)
//   {
//#ifdef PARANOID
//    if (v_values[i]==0.0)
//     {
//      std::ostringstream error_message;
//      error_message << "Zero entry in diagonal of velocity mass matrix\n"
//                    << "Index: " << i << std::endl;
//      throw OomphLibError(
//       error_message.str(),
//       "ConstrainedNavierStokesSchurComplementPreconditioner::assemble_inv_press_and_veloc_mass_matrix_diagonal()",
//       OOMPH_EXCEPTION_LOCATION);
//     }
//#endif
//    v_values[i] = 1.0/v_values[i];   
//    v_column_index[i] = v_first_row + i;
//    v_row_start[i] = i;
//   }
//  v_row_start[v_nrow_local] = v_nrow_local;
//  
//  // Build the velocity mass matrix
//  inv_v_mass_pt = new CRDoubleMatrix(this->block_distribution_pt(0));
//  inv_v_mass_pt->build_without_copy(v_nrow,v_nrow_local,
//                                    v_values,v_column_index,
//                                    v_row_start);
//  
//  // Create pressure mass matrix
//  if (!Use_LSC)
//   {
//    // Create column index and row start for pressure mass matrix
//    int* p_column_index = new int[p_nrow_local];
//    int* p_row_start = new int[p_nrow_local+1];
//    for (unsigned i = 0; i < p_nrow_local; i++)
//     {
//      
//#ifdef PARANOID
//      if (p_values[i]==0.0)
//       {
//        std::ostringstream error_message;
//        error_message << "Zero entry in diagonal of pressure mass matrix\n"
//                      << "Index: " << i << std::endl;
//        throw OomphLibError(
//         error_message.str(),
//         "ConstrainedNavierStokesSchurComplementPreconditioner::assemble_inv_press_and_veloc_mass_matrix_diagonal()",
//         OOMPH_EXCEPTION_LOCATION);
//       }
//#endif
//      p_values[i] = 1.0/p_values[i];
//      
//      p_column_index[i] = p_first_row + i;
//      p_row_start[i] = i;
//     }
//    p_row_start[p_nrow_local] = p_nrow_local;
//    
//    // Build the pressure mass matrix
//    inv_p_mass_pt = new CRDoubleMatrix(this->block_distribution_pt(1));
//    inv_p_mass_pt->build_without_copy(p_nrow,p_nrow_local,
//                                      p_values,p_column_index,
//                                      p_row_start);
//
//   }
//  
// }

 void ConstrainedNavierStokesSchurComplementPreconditioner::
     assemble_inv_press_and_veloc_mass_matrix_diagonal(
  CRDoubleMatrix*& inv_p_mass_pt,
  CRDoubleMatrix*& inv_v_mass_pt,
  const bool& do_both,
  const unsigned& procnumber,Problem* problem_pt)
 {
  int pronumber = (int)procnumber;
  // determine the velocity rows required by this processor
  // RRR_type
  unsigned v_first_row = this->block_distribution_pt(pronumber)->first_row();
  unsigned v_nrow_local = this->block_distribution_pt(pronumber)->nrow_local();
  unsigned v_nrow = this->block_distribution_pt(pronumber)->nrow();
 
 //cout << "v_first_row = " << v_first_row << endl;
 //cout << "v_nrow_local = " << v_nrow_local << endl;
 //cout << "v_nrow = " << v_nrow << endl;
 
  // create storage for the diagonals
  double* v_values = new double[v_nrow_local];
  for (unsigned i = 0; i < v_nrow_local; i++)
   {
    v_values[i] = 0.0;
   }
 
  // Equivalent information for pressure mass matrix (only needed for 
  // Fp version)
  unsigned p_first_row=0;
  unsigned p_nrow_local=0;
  unsigned p_nrow=0;
  double* p_values = 0;

  if (!Use_LSC)
   {
    // determine the pressure rows required by this processor
    p_first_row = this->block_distribution_pt(1)->first_row();
    p_nrow_local = this->block_distribution_pt(1)->nrow_local();
    p_nrow = this->block_distribution_pt(1)->nrow();
  
    // create storage for the diagonals
    p_values = new double[p_nrow_local];
    for (unsigned i = 0; i < p_nrow_local; i++)
     {
      p_values[i] = 0.0;
     }
   } // if (!Use_LSC)
   
  // store the problem pt, this seems to be not used.
  //const Problem* problem_pt = this->problem_pt();

  // if the problem is distributed
  bool distributed = false;
  
  if (distributed)
   {
   // To do
   }
  else
   {
    // find number of elements
    unsigned n_el = Navier_stokes_mesh_pt->nelement();
    
    //unsigned n_el = problem_pt->mesh_pt(0)->nelement();
    // Fp needs pressure and velocity mass matrices
    unsigned which_one=0;
    if (Use_LSC) 
      which_one=2;
    
    // get the contribution for each element
    for (unsigned e = 0; e < n_el; e++)
     {
      // Get element
      //GeneralisedElement* el_pt=problem_pt->mesh_pt(0)->element_pt(e);
      GeneralisedElement* el_pt=Navier_stokes_mesh_pt->element_pt(e);

      // find number of degrees of freedom in the element
      // (this is slightly too big because it includes the
      // pressure dofs but this doesn't matter)
      unsigned el_dof = el_pt->ndof();
      
      // allocate local storage for the element's contribution to the
      // pressure and velocity mass matrix diagonal
      Vector<double> el_vmm_diagonal(el_dof);
      Vector<double> el_pmm_diagonal(el_dof);
      
      dynamic_cast<TemplateFreeNavierStokesEquationsBase*>(el_pt)->
       get_pressure_and_velocity_mass_matrix_diagonal( 
        el_pmm_diagonal,el_vmm_diagonal,which_one);
        
      // Get the contribution for each dof
      for (unsigned i = 0; i < el_dof; i++)
       {
        //Get the equation number
        unsigned eqn_number = el_pt->eqn_number(i);
        
        // Get the velocity dofs
        if (this->block_number(eqn_number)==pronumber) // RAY_TOCHANGE
         {
         //cout << "GOT HERE!!!!" << endl;
          // get the index in the block
          unsigned index = this->index_in_block(eqn_number);
          
          // if it is required on this processor
          if ((index >= v_first_row) &&
              (index < (v_first_row + v_nrow_local) ) )
           {
           //cout << "ZOMG Got in if" << endl;
           //cout << "index-v_first_row = " << index-v_first_row << endl;
            v_values[index-v_first_row] += el_vmm_diagonal[i];
           }
         }
        // Get the pressure dofs
        // NOTE: This is not used for the LSC case.
        else if (this->block_number(eqn_number)==1) // RAY_TOCHANGE
         {
          if (!Use_LSC)
           {
            // get the index in the block
            unsigned index = this->index_in_block(eqn_number);
            
            // if it is required on this processor
            if ((index >= p_first_row)&&
                (index < (p_first_row + p_nrow_local)) )
             {
              p_values[index-p_first_row] += el_pmm_diagonal[i];
             }
           } // if (!Use_LSC)
         }
       }
     } // for (unsigned e = 0; e < n_el; e++)
   }// if (distributed), else
   
  // Create column index and row start for velocity mass matrix
  int* v_column_index = new int[v_nrow_local]; 
  int* v_row_start = new int[v_nrow_local+1];
  for (unsigned i = 0; i < v_nrow_local; i++)
   {
#ifdef PARANOID
    if (v_values[i]==0.0)
     {
      std::ostringstream error_message;
      error_message << "Zero entry in diagonal of velocity mass matrix\n"
                    << "Index: " << i << std::endl;
      throw OomphLibError(
       error_message.str(),
       "ConstrainedNavierStokesSchurComplementPreconditioner::assemble_inv_press_and_veloc_mass_matrix_diagonal()",
       OOMPH_EXCEPTION_LOCATION);
     }
#endif
    v_values[i] = 1.0/v_values[i];   
    v_column_index[i] = v_first_row + i;
    v_row_start[i] = i;
   }
  v_row_start[v_nrow_local] = v_nrow_local;
  
  // Build the velocity mass matrix
  inv_v_mass_pt = new CRDoubleMatrix(this->block_distribution_pt(pronumber));
  inv_v_mass_pt->build_without_copy(v_nrow,v_nrow_local,
                                    v_values,v_column_index,
                                    v_row_start);
   
  // Create pressure mass matrix
  if (!Use_LSC)
   {
    // Create column index and row start for pressure mass matrix
    int* p_column_index = new int[p_nrow_local];
    int* p_row_start = new int[p_nrow_local+1];
    for (unsigned i = 0; i < p_nrow_local; i++)
     {
      
#ifdef PARANOID
      if (p_values[i]==0.0)
       {
        std::ostringstream error_message;
        error_message << "Zero entry in diagonal of pressure mass matrix\n"
                      << "Index: " << i << std::endl;
        throw OomphLibError(
         error_message.str(),
         "ConstrainedNavierStokesSchurComplementPreconditioner::assemble_inv_press_and_veloc_mass_matrix_diagonal()",
         OOMPH_EXCEPTION_LOCATION);
       }
#endif
      p_values[i] = 1.0/p_values[i];
      
      p_column_index[i] = p_first_row + i;
      p_row_start[i] = i;
     }
    p_row_start[p_nrow_local] = p_nrow_local;
    
    // Build the pressure mass matrix
    inv_p_mass_pt = new CRDoubleMatrix(this->block_distribution_pt(1)); // This also needs to change.
    inv_p_mass_pt->build_without_copy(p_nrow,p_nrow_local,
                                      p_values,p_column_index,
                                      p_row_start);

   } // if (!Use_LSC)
   
 }

//=========================================================================
/// Helper function to delete preconditioner data.
//=========================================================================
 void ConstrainedNavierStokesSchurComplementPreconditioner::clean_up_memory()
 {
  if (Preconditioner_has_been_setup)
   {
    // delete matvecs
    delete Bt_mat_vec_pt; 
    Bt_mat_vec_pt = 0;

    delete F_mat_vec_pt;
    F_mat_vec_pt = 0;

    delete QBt_mat_vec_pt;
    QBt_mat_vec_pt = 0;

    delete E_mat_vec_pt;
    E_mat_vec_pt = 0;
    
    // delete stuff from velocity solve
    if (Using_default_f_preconditioner)
     {
      delete F_preconditioner_pt;
      F_preconditioner_pt = 0;
     }
    
    // delete stuff from Schur complement approx
    if (Using_default_p_preconditioner)
     {
      delete P_preconditioner_pt;
      P_preconditioner_pt = 0;
     }
   }
 }


}// end oomph namespace
