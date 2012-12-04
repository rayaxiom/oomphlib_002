//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC//           Version 0.85. June 9, 2008.
//LIC// 
//LIC// Copyright (C) 2006-2008 Matthias Heil and Andrew Hazel
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

#include <sstream>

// Oomph-lib includes
#include "generic.h"
#include "navier_stokes.h"
#include "raymon.h"

// The 3D mesh
#include "meshes/simple_cubic_mesh.h"
#include "meshes/simple_cubic_tet_mesh.h"

using namespace std;
using namespace oomph;


double raycos(const double &rad)
{
  double temp = std::cos(rad);
  if(std::abs(temp) < 1e-10)
  {
    return 0;
  }
  else
  {  
    return temp;
  }
}


double raysin(const double &rad)
{
  double temp = std::sin(rad);
  if(std::abs(temp) < 1e-10)
  {
    return 0;
  }
  else
  {
    return temp;
  }
}




//===start_of_namespace=================================================
/// Namepspace for global parameters
//======================================================================
namespace Global_Parameters
{
 /// Current Reynolds number
 double Re = 0.0;
 
 double Ang_x = 0.0;
 double Ang_y = 0.0;
 double Ang_z = 0.0;
 
 /// Lagrange Multiplier ID 
 const unsigned Lagrange_multiplier_po = 42;
 const unsigned Lagrange_multiplier_ib = 43;
 /// Constant for pi
 const double Pi = 4.0*atan(1.0);


 /// sigma for the problem
 double Sigma = 0.0;
 
 string Prec_str = ""; 
 string Dim_str = "";
 string Prob_str = "";
 string Vis_str = "";
 string Ang_x_str = "";
 string Ang_y_str = "";
 string Ang_z_str = "";
 string Rey_str = "";
 string Noel_str = "";
 string Current_settings = "";
 
 /// Storage for number of iterations during Newton steps 
 Vector<unsigned> Iterations;
 
 /// Storage for linear solver times during Newton steps 
 Vector<double> Linear_solver_time;
 
 // counter for the current newtown steps
 unsigned before_newton_step_counter = 0;
 // flag to determine if this is the first iteration of the newton solve
 bool first_solve = true;
 // counter for iterations within each newton step
 unsigned its = 0;


 bool Dump_matrices = false;
 bool Doc_solution = false;
 bool Use_axnorm = false;

 bool Use_lsc = false;
} // Global_Parameters


namespace oomph {

//=============================================================================
/// \short rayupdate
//=============================================================================
 class ParallelOutflowPreconditioner 
  : public BlockPreconditioner<CRDoubleMatrix>
 {

   public:
  
  /// \short rayupdate
  typedef Preconditioner* (*SubsidiaryPreconditionerFctPt)();
  
  /// \short The augmented fluid system can be solved in
  /// one of four ways. \n
  /// 0 - Exact preconditioner \n
  /// 1 - Exact LSC preconditioner \n
  /// The type of preconditioning used will determine how
  /// the sub-blocks are merged together.
  enum Lagrange_preconditioner_type {Exact_block_preconditioner,
                                     Exact_lsc_block_preconditioner};

  /// \short Default (and only) constructor.
  ParallelOutflowPreconditioner()
   {
    // How many meshes do we have?
    //Lagrange_multiplier_mesh_pt = 0;
    //Fluid_mesh_pt = 0;

  // Pointer to the 'preconditioner' for the pressure matrix
  P_preconditioner_pt = 0;
  
  // Pointer to the 'preconditoner' for the F matrix
  F_preconditioner_pt = 0;
  NS_preconditioner_pt=0;
  
  // Pointer to the 'preconditoner' for the W matrix
  // This has been changed to a vector of vectors.
  //W_preconditioner_pt = 0;
  
  // flag to indicate whether the default F preconditioner is used
  Using_default_f_preconditioner = true;
  
  // flag to indicate whether the default p preconditioner is used
  Using_default_p_preconditioner = true;
  
  // flag to indicate whether the default w preconditioner is used
  Using_default_w_preconditioner = true;
  
  // flag to indicate LSC preconditioner 
  L_prec_type = Exact_block_preconditioner;
  Use_default_norm_of_f_scaling = true;
  Scaling_sigma = 0.0;

  N_dof_types = 0;
  N_lagrange_dof_types = 0;
  N_fluid_dof_types = 0;
  
  N_velocity_dof_types = 0;

  Fluid_block_size = 0;
  Pressure_block_size = 0;
  Velocity_block_size = 0;
   }
  
  /// destructor
  virtual ~ParallelOutflowPreconditioner()
   {
    this->clean_up_memory();
   }
  
  /// Broken copy constructor
  ParallelOutflowPreconditioner
   (const ParallelOutflowPreconditioner&)
   { 
    BrokenCopy::broken_copy("ParallelOutflowPreconditioner");
   } 
  
  /// Broken assignment operator
  void operator=
   (const ParallelOutflowPreconditioner&) 
   {
    BrokenCopy::broken_assign(" ParallelOutflowPreconditioner");
   }
  
  /// Dump out some useful stuff
  void dump_it(Problem* problem_pt)
   {
    cout << "ndof tyes: " << ndof_types() << std::endl;
   }

  /// Setup method for the ParallelOutflowPreconditioner.
  void setup(Problem* problem_pt, DoubleMatrixBase* matrix_pt);
  
  /// \short Get the infinity norm of a matrix.
  ///  The matrix may be composed of several matrices.
  void get_inf_norm(DenseMatrix<CRDoubleMatrix* > &matrix_pt, 
                    double &max_row_value);
  
  /// Merge submatrices into a single matrix.
  void merge(DenseMatrix<CRDoubleMatrix* > &matrix_pt, 
             CRDoubleMatrix *&block_pt);
  
  /// Add a scalar to each of the diagonal entry of a matrix.
  void add_scaling_to_diag(double &Scaling, CRDoubleMatrix *&block_pt);
  
  /// Extract the diagonal entries of a matrix. 
  void get_diag(CRDoubleMatrix *&block_pt, Vector<double>& diag);
  
  /// Element-wise addition of two matrices.
  void add_matrices(CRDoubleMatrix *&block_pt1, CRDoubleMatrix *&block_pt2);
                                                 
  /// \short Apply the preconditioner. Method implemented in two
  /// other methods (elastic and lagrange multiplier subsidiary
  /// preocnditioner) for the PseudoElasticFSIPreconditioner
  /// r is the residual (rhs), z will contain the solution.
  void preconditioner_solve(const DoubleVector& r, DoubleVector& z)
  {
   std::stringstream NewtonStepStream;
   unsigned currentNS = Global_Parameters::before_newton_step_counter - 1;
   NewtonStepStream << "NS" << currentNS; // Because the newton step counter
                                          // has already been incremented before
                                          // this current Newton solve
   string newton_step_counter = NewtonStepStream.str();
   
   std::stringstream ItsStream;
   ItsStream << "i" << Global_Parameters::its;
   string its_counter = ItsStream.str();
   //cout << "RAYRAY: NS: " << currentNS << "its: " << Global_Parameters::its << endl;
   
   string currentsetting = Global_Parameters::Current_settings
                           + newton_step_counter
                           + its_counter;
                           
                           
// NEW STUFF
   DoubleVector temp_vec;
   DoubleVector another_temp_vec;
   DoubleVector yet_another_temp_vec;
   bool Doc_time = false;
   // Note: We have:
   // 0  1   2  3    4  5
   // u  u_c v  v_c  p  L
   // n_dof_types.
   
////////////////////////////////////////////////////////////////////////////////
   // First we solve all the r_l blocks:
   
   for(unsigned l_i = 0; l_i < N_lagrange_dof_types; l_i++)
   {
     unsigned l_ii = N_fluid_dof_types + l_i;
     this->get_block_vector(Dof_type_list[l_ii],r,temp_vec);
     W_preconditioner_pt[l_i]->preconditioner_solve(temp_vec, 
                                                    another_temp_vec);
     this->return_block_vector(Dof_type_list[l_ii],another_temp_vec,z);
     temp_vec.clear();
     another_temp_vec.clear();
   }                           
// At this point, all vectors are cleared.

if(L_prec_type == Exact_block_preconditioner)
{
  // Merge the fluid block vectors
  double t_get_merge_fluid_rhs_start = TimingHelpers::timer();

  bool distributed = this->master_distribution_pt()->distributed();
  LinearAlgebraDistribution* new_distribution_pt 
    = new LinearAlgebraDistribution(problem_pt()->communicator_pt(),
                                    Fluid_block_size,distributed);
  
  temp_vec.build(new_distribution_pt,0.0);
  // v_nrow_i, for indexing the new merged vector.
  unsigned long v_nrow_i = 0;
  // Loop through the fluid rhs block vectors
  for(unsigned i = 0; i < N_fluid_dof_types; i++)
  {
    this->get_block_vector(Dof_type_list[i],r,another_temp_vec);
    unsigned long current_block_nrow = another_temp_vec.nrow();

    // Loop through the entries in this block vector
    for(unsigned long current_block_i = 0; 
        current_block_i < current_block_nrow;
        current_block_i++)
    {
      temp_vec[v_nrow_i] = another_temp_vec[current_block_i];
      v_nrow_i++;
    } // for
    another_temp_vec.clear();
  } // for
   double t_get_merge_fluid_rhs_finish = TimingHelpers::timer();
   if(Doc_time && Global_Parameters::first_solve)
    {
      double t_get_merge_fluid_rhs_time = t_get_merge_fluid_rhs_finish 
                                    - t_get_merge_fluid_rhs_start;
      std::cout << "t_get_merge_fluid_rhs_time: " 
                << t_get_merge_fluid_rhs_time << endl;
    }

    // temp_vec contains the fluid rhs.
    NS_preconditioner_pt->preconditioner_solve(temp_vec,another_temp_vec);
    temp_vec.clear();

    // We now have to put the block vectors in another_temp_vec back!
    double t_put_back_rhs_start = TimingHelpers::timer();
    unsigned long merged_vec_row_i = 0;
    // loop through the fluid block vectors
    // to extract entries from the merged vector.
    for(unsigned i = 0; i < N_fluid_dof_types; i++)
    {
      this->get_block_vector(Dof_type_list[i],r,temp_vec);
      unsigned long current_block_nrow = temp_vec.nrow();

      // loop through the entries of this block vector
      for(unsigned long current_block_i = 0;
          current_block_i < current_block_nrow;
          current_block_i++)
      {
        temp_vec[current_block_i] = another_temp_vec[merged_vec_row_i];
        merged_vec_row_i++;
      }
      this->return_block_vector(Dof_type_list[i],temp_vec,z);
      temp_vec.clear();
    }//for
    double t_put_back_rhs_finish = TimingHelpers::timer();
    if(Doc_time && Global_Parameters::first_solve)
    {
      double t_put_back_rhs_time = t_put_back_rhs_finish
                                   - t_put_back_rhs_start;
      std::cout << "t_put_back_rhs_time: "
                << t_put_back_rhs_time << std::endl;
    }


}//if(L_prec_type == Exact_block_preconditioner)
else if(L_prec_type == Exact_lsc_block_preconditioner)
{
  // Solve the pressure block, then velocity.
  ////////////////////////////////////////////////////////////////////////////////
   // Now solve for the pressure
   
   // Get r_p
   this->get_block_vector(Dof_type_list[N_velocity_dof_types],r,temp_vec);
   // passed
   //temp_vec.output("temp_vec");
   
   P_preconditioner_pt->preconditioner_solve(temp_vec,another_temp_vec);
   temp_vec.clear();
   //passed
   //another_temp_vec.output("another_temp_vec");
   
   QBt_mat_vec_pt->multiply(another_temp_vec,temp_vec);
   another_temp_vec.clear();
   //passed
   //temp_vec.output("temp_vec");
   
   F_mat_vec_pt->multiply(temp_vec, another_temp_vec);
   // passed
   //another_temp_vec.output("another_temp_vec");
   
   temp_vec.clear();
   
   QBt_mat_vec_pt->multiply_transpose(another_temp_vec,temp_vec);
   // passed
   //temp_vec.output("temp_vec");
   
   another_temp_vec.clear();
   P_preconditioner_pt->preconditioner_solve(temp_vec,another_temp_vec);
   
   //passed
   //another_temp_vec.output("another_temp_vec");
   
   temp_vec.build(another_temp_vec.distribution_pt(),0.0);
   temp_vec -= another_temp_vec;
   // passed
   //temp_vec.output("temp_vec");
   
   
   this->return_block_vector(Dof_type_list[N_velocity_dof_types],temp_vec,z);
   //pause("testo!");
///////////////////////////////////////////////////////////////////////////////
  
   temp_vec.clear();
   Bt_mat_vec_pt->multiply(another_temp_vec,temp_vec);
   // temp_vec now contained -BtZp
   
   another_temp_vec.clear();
   double t_get_merge_rhs_start = TimingHelpers::timer();
   // merge the rhs vector:
   unsigned long v_nrow = 0;
   for(unsigned i = 0; i < N_velocity_dof_types; i++)
   {
     this->get_block_vector(Dof_type_list[i],r,another_temp_vec);
     v_nrow += another_temp_vec.nrow();
     another_temp_vec.clear();
   }
   
   bool distributed = this->master_distribution_pt()->distributed();
   LinearAlgebraDistribution* new_distribution_pt 
     = new LinearAlgebraDistribution(problem_pt()->communicator_pt(),
                                     v_nrow,distributed);
   
   yet_another_temp_vec.build(new_distribution_pt,0.0);
   //yet_another_temp_vec.output("yet_another_temp_vec");
   
   
   unsigned long v_nrow_i = 0;
   for(unsigned i = 0; i < N_velocity_dof_types; i++)
   {
     this->get_block_vector(Dof_type_list[i],r,another_temp_vec);
     unsigned long current_block_nrow = another_temp_vec.nrow();
     
     for(unsigned long current_block_i = 0; current_block_i<current_block_nrow;
         current_block_i++)
     {
       yet_another_temp_vec[v_nrow_i] = another_temp_vec[current_block_i];
       v_nrow_i++;
     }
     //std::stringstream blockstringsteam;
     //blockstringsteam << "vec_block_" << Dof_type_list[i];
     //another_temp_vec.output(blockstringsteam.str());
     another_temp_vec.clear();
   }
   double t_get_merge_rhs_finish = TimingHelpers::timer();
   
   if(Doc_time && Global_Parameters::first_solve)
     {
       double t_get_merge_rhs_time = t_get_merge_rhs_finish 
                                     - t_get_merge_rhs_start;
       std::cout << "t_get_merge_rhs_time: " << t_get_merge_rhs_time << endl;
     }
   // passed
   //yet_another_temp_vec.output("yet_another_temp_vec");
   
   //yet_another_temp_vec contained the merged blocks in the order we want.
   
   yet_another_temp_vec += temp_vec;
   
   // yet_another_temp_vec now contains r_u - G z_p
   temp_vec.clear();
   
   F_preconditioner_pt->preconditioner_solve(yet_another_temp_vec,temp_vec);
   
   // Now we have to put temp_vec back.
   double t_put_back_rhs_start = TimingHelpers::timer();
   unsigned long global_row_i = 0;
   for(unsigned i = 0; i< N_velocity_dof_types; i++)
   {
     another_temp_vec.clear();
     
     this->get_block_vector(Dof_type_list[i],r,another_temp_vec);
     
     unsigned long current_block_nrow = another_temp_vec.nrow();
     
     for(unsigned long current_block_nrow_i = 0; 
         current_block_nrow_i < current_block_nrow;
         current_block_nrow_i++)
     {
       another_temp_vec[current_block_nrow_i] = temp_vec[global_row_i];
       global_row_i++;
     }
     
     this->return_block_vector(Dof_type_list[i],another_temp_vec,z);
   } //  loop though the block vectors
   double t_put_back_rhs_finish = TimingHelpers::timer();
   
   if(Doc_time && Global_Parameters::first_solve)
     {
       double t_put_back_rhs_time = t_put_back_rhs_finish 
                                    - t_put_back_rhs_start;
       std::cout << "t_put_back_rhs_time: " 
                 << t_put_back_rhs_time << std::endl;
     }
 

}// Exact_lsc_block_preconditioner



////////////////////////////////////////////////////////////////////////////////
  
  DoubleVector x; // Will contain the re-ordered rhs
  this->get_block_ordered_preconditioner_vector(r,x);
   
    
//   DoubleVector y; // Will contain the solution
   
//   NS_preconditioner_pt->preconditioner_solve(x,y);
   // RRR_DUMP
   //*
   // Set up the string to record the current newton step
   if(Global_Parameters::Dump_matrices)
   {
     if(Global_Parameters::first_solve)
     {                       
       std::stringstream rhsxstream;
       rhsxstream << "rhsx_" << currentsetting;
       x.output(rhsxstream.str());
   
       //std::stringstream precstream;
       //precstream << "prec_" << currentsetting;
       //s_prec_pt->sparse_indexed_output(precstream.str());
   
       // Setting up the submatrix dimension file.
       //ofstream precdimensionfile;
       //string precdimensionstring = "precdim_" + currentsetting;
        
       //precdimensionfile.open(precdimensionstring.c_str());
      //int nrows = s_prec_pt->nrow ();
      //int ncols = s_prec_pt->ncol ();
      //precdimensionfile << nrows << " " << ncols << endl;
      //precdimensionfile.close();
   
     }
   } // if(Dump_matrices)
   // */
   //cout << "before J_prec_pt->solve(x,y)" << endl;
   //J_prec_pt->solve(x,y); // x is rhs, y is soln.
   //cout << "after  J_prec_pt->solve(x,y)" << endl; 

   // putting y into z but blocked.
  // this->return_block_ordered_preconditioner_vector(y,z);
   
   // RRR_DUMP
   /*
   if(Global_Parameters::first_solve)
   {                       
   std::stringstream solnystream;
   solnystream << "solny_" << currentsetting;
   y.output(solnystream.str());
   
   std::stringstream solnzstream;
   solnzstream << "solnz_" << currentsetting;
   z.output(solnzstream.str());
   }
   // */
   
   //we do not want to output the rhs and solution for subsequent solves
   Global_Parameters::first_solve = false;
   //cout << "preconditioner_solve(const DoubleVector& r, DoubleVector& z)" 
   //     << endl;
   ////Preconditioner_pt->preconditioner_solve(x,y);
   //cout << Preconditioner_pt->iterations() << endl;
   //this->elastic_preconditioner_solve(r,z);
   // this->lagrange_multiplier_preconditioner_solve(r,z);
   
   Global_Parameters::its++;
   
  } // end of preconditioner_solve
  
  /// \short Access function to mesh containing the block-preconditionable
  /// fluid elements
  //void set_fluid_mesh(Mesh* mesh_pt) 
  //{
  // Fluid_mesh_pt = mesh_pt;
  //}
  
  /// \short Access function to mesh containing the block-preconditionable
  /// lagrange multiplier elements 
  //void set_lagrange_multiplier_mesh(Mesh* mesh_pt) 
  //{
  // Lagrange_multiplier_mesh_pt = mesh_pt;
  //}

  void set_meshes(Vector<Mesh*> meshes_pt)
  {
    Meshes_pt = meshes_pt;
    
    unsigned nmeshes = Meshes_pt.size();

    this->set_nmesh(nmeshes);

  }

  /// \short Access function to the Scaling sigma of the preconditioner
  double& scaling_sigma()
  {
    Use_default_norm_of_f_scaling = false;
    return Scaling_sigma;
  }
/// \short Function to get the scaling Sigma of the preconditioner  
  double scaling_sigma() const
  {
    return Scaling_sigma;
  }

  void use_default_norm_of_f_scaling()
  {
    Use_default_norm_of_f_scaling = true;
  }
   /// \short Helper function to assemble the diagonal of the pressure
   /// and velocity mass matrices from the elemental contributions defined in
   /// NavierStokesEquations<DIM>.
   /// If do_both=true, both are computed, otherwise only the velocity
   /// mass matrix (the LSC version of the preconditioner only needs
   /// that one)
   void assemble_inv_press_and_veloc_mass_matrix_diagonal(
    CRDoubleMatrix*& inv_p_mass_pt, 
    CRDoubleMatrix*& inv_v_mass_pt, 
    const bool& do_both,
    const unsigned& procnumber);
  
///raydo comment
  void use_lsc()
  {
    L_prec_type = Exact_lsc_block_preconditioner;
  }
///raydo comment
  void use_exact()
  {
    L_prec_type = Exact_block_preconditioner;
  }
  /// \short Clears the memory.
  void clean_up_memory();
  
   private:
  
  /// \short the dimension of the problem
  unsigned Dim;
  
  /// \short the Scaling_sigma variable of this preconditioner
  double Scaling_sigma;
  
  Vector<Mesh*> Meshes_pt;

  /// Pointer to the Navier Stokes preconditioner
 // NavierStokesLSCPreconditioner * Navier_stokes_preconditioner_pt;
  
  /// Pointer to the W block preconditioner
 // Preconditioner* Preconditioner_pt;
  
  
  ////////// NEW STUFF
  bool Use_default_norm_of_f_scaling;
  
  Lagrange_preconditioner_type L_prec_type;
 
  // P_prec and F_prec are for LSC.
  // Pointer to the 'preconditioner' for the pressure matrix
  Preconditioner* P_preconditioner_pt;
  
  // Pointer to the 'preconditoner' for the F matrix
  Preconditioner* F_preconditioner_pt;
  
  // Pointer to the 'preconditioner' for the Navier-Stokes block
  Preconditioner* NS_preconditioner_pt;

  // Pointer to the 'preconditoner' for the W matrix
  Vector<Preconditioner*> W_preconditioner_pt;
  
  // flag to indicate whether the default F preconditioner is used
  bool Using_default_f_preconditioner;
  
  // flag to indicate whether the default p preconditioner is used
  bool Using_default_p_preconditioner;
  
  // flag to indicate whether the default ns preconditioner is used
  
  // flag to indicate whether the default w preconditioner is used
  bool Using_default_w_preconditioner;
  
  //

  bool Preconditioner_has_been_setup;
  
  bool F_preconditioner_is_block_preconditioner;
  
  MatrixVectorProduct* QBt_mat_vec_pt;
  
  MatrixVectorProduct* Bt_mat_vec_pt;
  
  MatrixVectorProduct* F_mat_vec_pt;
  
  // the re-arraned dof types.
  Vector<unsigned> Dof_type_list;
 
  // These are assigned in the setup but used in
  // preconditioner_solve() to re-arrange blocks.
  Vector<unsigned long> Dof_type_block_size;
  unsigned long Fluid_block_size;
  unsigned long Pressure_block_size;
  unsigned long Velocity_block_size;
  
  //
  unsigned N_dof_types;
  unsigned N_lagrange_dof_types;
  unsigned N_fluid_dof_types;
  unsigned N_velocity_dof_types;
  
 }; // end of ParallelOutflowPreconditioner class



//========================================================================
/// Helper function to assemble the diagonal of the velocity
/// mass matrix from the elemental contributions defined in
/// NavierStokesEquations<DIM>::get_velocity_mass_matrix_diagonal(...).
/// If do_both=true, both are computed, otherwise only the velocity
/// mass matrix (the LSC version of the preconditioner only needs
/// that one)
//========================================================================
 void ParallelOutflowPreconditioner:: 
  assemble_inv_press_and_veloc_mass_matrix_diagonal(
  CRDoubleMatrix*& inv_p_mass_pt,
  CRDoubleMatrix*& inv_v_mass_pt,
  const bool& do_both,
  const unsigned& procnumber)
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

  if (L_prec_type != Exact_lsc_block_preconditioner)
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
    unsigned n_el = Meshes_pt[0]->nelement();
    
    // Fp needs pressure and velocity mass matrices
    unsigned which_one=0;
    if (L_prec_type == Exact_lsc_block_preconditioner) 
      which_one=2;
    
    // get the contribution for each element
    for (unsigned e = 0; e < n_el; e++)
     {
      // Get element
      GeneralisedElement* el_pt=Meshes_pt[0]->element_pt(e);
      
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
          if (L_prec_type != Exact_lsc_block_preconditioner)
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
       "NavierStokesSchurComplementPreconditioner::assemble_inv_press_and_veloc_mass_matrix_diagonal()",
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
  if (L_prec_type != Exact_lsc_block_preconditioner)
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
         "NavierStokesSchurComplementPreconditioner::assemble_inv_press_and_veloc_mass_matrix_diagonal()",
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
// void ParallelOutflowPreconditioner::assemble_inv_press_and_veloc_mass_matrix_diagonal



/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////

 //========================================================================
 /// add the scaled identity matrix to the specified block
 //========================================================================
  void ParallelOutflowPreconditioner::add_matrices(CRDoubleMatrix *&block_pt1, 
                                                   CRDoubleMatrix *&block_pt2)
  {
    double* values1 = block_pt1->value();
    int* column_index1 = block_pt1->column_index();
    int* row_start1 = block_pt1->row_start();
    unsigned nrow_local1 = block_pt1->nrow_local();
    //unsigned first_row1 = block_pt1->first_row();
    unsigned nrow_global1 = block_pt1->nrow();
    
    double* values2 = block_pt2->value();
    int* column_index2 = block_pt2->column_index();
    int* row_start2 = block_pt2->row_start();
    //unsigned nrow_local2 = block_pt2->nrow_local();
    //unsigned first_row2 = block_pt2->first_row();
    
    // nrow_local should be the same.
    Vector< set<int> > column_set(nrow_local1);
    Vector<int> new_row_start(nrow_local1+1, 0);
    unsigned total_entries = 0;
    
    // Check which columns appear twice.
    pair<set<int>::iterator,bool> ret;
    
    for(unsigned i = 0; i < nrow_local1; i++)
    {
      // Union of both column indices
      for (int j = row_start1[i]; j < row_start1[i+1]; j++)
      {
        column_set[i].insert(column_index1[j]);
      }
      
      for (int j = row_start2[i]; j < row_start2[i+1]; j++)
      {
        ret = column_set[i].insert(column_index2[j]);
      }
      
      total_entries += column_set[i].size();
      new_row_start[i+1] = total_entries;
    }
    
    //Vector<double> ray_values33(9);
    //Vector<int> ray_column_indices33(9);
    //Vector<int> ray_row_start33(6);
  
  //cout << "Total entries: " << total_entries << endl;
  Vector<int> new_column_index(total_entries, 0);
  Vector<double> new_values(total_entries, 0);
  
  unsigned nci = 0;
  for(unsigned i = 0; i < nrow_local1; i++)
  {
    for(set<int>::const_iterator p = column_set[i].begin();
        p != column_set[i].end(); p++)
    {
      new_column_index[nci] =  *p;
      nci++;
      //cout << *p << " ";
    }
    //cout << endl;
  }
  
  
  
  // Now we add all the entries.
  
    for (unsigned i = 0; i < nrow_local1; i++)
      {
        // Adding the first matrix to the new matrix
        // Loop through the entries on ith row of matrix 1
        for (int j = row_start1[i]; j < row_start1[i+1]; j++)
        {
          bool found = false;
          // loop through the entries on the ith row of result matrix
          for (int jj = new_row_start[i]; jj  < new_row_start[i+1] && !found; jj++)
          {
            if(column_index1[j] == new_column_index[jj])
            {
              new_values[jj] += values1[j];
              found =true;
            }
          }
        }
        
        // Loop through the entries on ith row of matrix 2
        for (int j = row_start2[i]; j < row_start2[i+1]; j++)
        {
          bool found = false;
          // loop through the entries on the ith row of result matrix
          for (int jj = new_row_start[i]; jj  < new_row_start[i+1] && !found; jj++)
          {
            if(column_index2[j] == new_column_index[jj])
            {
              new_values[jj] += values2[j];
              found =true;
            }
          }
        }
      }
  
  // FIX THIS! This currently only works with square matrices
  block_pt2->build(nrow_global1,new_values,new_column_index,new_row_start);
  
  /// ENDS HERE
  /*
  cout << "new_column_index = " << endl;
  for(int i = 0; i < total_entries; i++)
  {
    cout << new_column_index[i] << " ";
  }
  cout << endl;
  
  cout << "new_row_start = ";
  for(int i = 0; i < nrow_local1 + 1; i++)
  {
    cout << new_row_start[i] << " ";
  }
  cout << endl;
  */
  } // add_matrices
 
 
  void ParallelOutflowPreconditioner::get_inf_norm(DenseMatrix<CRDoubleMatrix* > 
                                                     &matrix_pt, 
                                                   double &max_row_value)
  {
    unsigned long matrix_nrow = matrix_pt.nrow();
    unsigned long matrix_ncol = matrix_pt.ncol();
    
    max_row_value = 0.0;
    //unsigned test_i = 0;
    // Loop through the block rows.
    for(unsigned row_block = 0; row_block < matrix_nrow; row_block++)
    {
      // Get the number of rows in this row_block from the first block.
      unsigned long block_nrow_local = matrix_pt(row_block,0)->nrow_local();
      
      // Loop through the number of local rows
      for(unsigned i = 0; i < block_nrow_local; i++)
      {
        double current_row_total = 0.0;
        // Loop through the column blocks on this row:
        for(unsigned column_block = 0; column_block < matrix_ncol; column_block++)
        {
          CRDoubleMatrix* current_block_pt = matrix_pt(row_block,column_block);
          double* current_block_values = current_block_pt->value();
          //int* current_block_column_indicies = current_block_pt->column_index();
          int* current_block_row_start = current_block_pt->row_start();
          
          for(int j = current_block_row_start[i]; 
              j < current_block_row_start[i+1]; j++)
          {
            current_row_total += fabs(current_block_values[j]);
          } // loop through the 
        } // loop through the columns
        
        max_row_value = max(max_row_value, current_row_total);
      } // loop through the local rows
    } // loop through the block rows
  } // void ParallelOutflowPreconditioner::get_inf_norm
  
  void ParallelOutflowPreconditioner::merge(DenseMatrix<CRDoubleMatrix* > 
                                              &matrix_pt,
                                            CRDoubleMatrix *&block_pt)
  {
    bool distributed = this->master_distribution_pt()->distributed();
    unsigned long matrix_nrow = matrix_pt.nrow();
    unsigned long matrix_ncol = matrix_pt.ncol();
    Vector<unsigned> block_cols(matrix_ncol);
    
    // get the block cols for offset
    for(unsigned col_block_i = 0; col_block_i < matrix_ncol; col_block_i++)
    {
      CRDoubleMatrix* current_block_pt = matrix_pt(0,col_block_i);
      block_cols[col_block_i] = current_block_pt->ncol();
    }
    
    // Get the ncol global
    unsigned total_ncol_global = 0;
    for(unsigned col_i = 0; col_i < matrix_ncol; col_i++)
    {
      CRDoubleMatrix* current_block_pt = matrix_pt(0,col_i);
      total_ncol_global += current_block_pt->ncol();
    }
    
    // Get the nrow global.
    unsigned total_nrow_global = 0;
    for(unsigned row_i = 0; row_i < matrix_nrow; row_i++)
    {
      CRDoubleMatrix* current_block_pt = matrix_pt(row_i,0);
      total_nrow_global += current_block_pt->nrow();
    } // for
    
    LinearAlgebraDistribution* new_distribution_pt 
      = new LinearAlgebraDistribution(problem_pt()->communicator_pt(),
                                      total_nrow_global,distributed);
    
    // Get the nnz in all matrices.
    unsigned long total_nnz = 0;
    // Loop through the block rows.
    for(unsigned row_block = 0; row_block < matrix_nrow; row_block++)
    {
      // Loop through the blocks on this row:
      for(unsigned column_block = 0; column_block < matrix_ncol; column_block++)
      {
        CRDoubleMatrix* current_block_pt = matrix_pt(row_block,column_block);
        total_nnz += current_block_pt->nnz();
      }
    }
    
    Vector<double> new_values(total_nnz);
    Vector<int> new_column_indices(total_nnz);
    Vector<int> new_row_start(total_nrow_global + 1);
    
    // Loop through the block rows.
    unsigned long new_val_i = 0;
    unsigned long new_row_start_i = 0;
    unsigned long n_empty_rows = 0;
    for(unsigned row_block = 0; row_block < matrix_nrow; row_block++)
    {
      // Get the number of rows in this row_block from the first block.
      unsigned long block_nrow_local = matrix_pt(row_block,0)->nrow_local();
      
      // Loop through the number of local rows
      for(unsigned i = 0; i < block_nrow_local; i++)
      {
        bool first_ele_in_row = true;
        // Loop through the column blocks on this row:
        for(unsigned column_block = 0; 
            column_block < matrix_ncol; column_block++)
        {
          CRDoubleMatrix* current_block_pt 
            = matrix_pt(row_block,column_block);
          double* current_block_values 
            = current_block_pt->value();
          int* current_block_column_indicies 
            = current_block_pt->column_index();
          int* current_block_row_start 
            = current_block_pt->row_start();
          
          // calculate the off-set
          unsigned long offset = 0;
          for(unsigned pre_col_block = 0; 
              pre_col_block < column_block; pre_col_block++)
          {
            offset += block_cols[pre_col_block];
          }
          // Because indices starts at zero in C++
          
          
          // We only go in here if there is a non-empty row.
          for(int j = current_block_row_start[i]; 
              j < current_block_row_start[i+1]; j++)
          {
            //cout << "new_val_i " << new_val_i << endl;
            new_values[new_val_i] = current_block_values[j];
            new_column_indices[new_val_i] = current_block_column_indicies[j] 
                                            + offset;
            
            // Filling in the row_start
            if(first_ele_in_row)
            {
              if(n_empty_rows != 0)
              {
                // We have to fill in the row start for all the zero rows.
                for(unsigned long empty_row_i = 0; 
                    empty_row_i < n_empty_rows; empty_row_i++)
                {
                  //cout << "new_row_start_i (empty row)" << new_row_start_i << endl;
                  // We fill it with the current row start
                  new_row_start[new_row_start_i] = new_val_i;
                  new_row_start_i++;
                } // for
                // reset the number of empty rows.
                n_empty_rows = 0;
              } // if
              //cout << "new_row_start_i (non-empty)" << new_row_start_i << endl;
              new_row_start[new_row_start_i] = new_val_i;
              new_row_start_i++;
              first_ele_in_row = false;
            } // if
            
            new_val_i++;
          } // for - looping the values

        } // for - looping through block columns
          // At the end of looping through all of the column blocks,
          // If true => no first element is reached => no row start => empty row.
          if(first_ele_in_row)
          {
            n_empty_rows++;
          } // if
      } // for - looping throug local rows
    } // for - looping through blockrows
    
    // If there are empty rows, then we fill them!
    if(n_empty_rows != 0)
    {
      // We have to fill in the row start for all the zero rows.
      for(unsigned long empty_row_i = 0; 
          empty_row_i < n_empty_rows; empty_row_i++)
      {
        // We fill it with the current row start
        new_row_start[new_row_start_i] = total_nnz;
        new_row_start_i++;
      } // for
      
      // reset the number of empty rows.
      n_empty_rows = 0;
    } // if
    
    new_row_start[total_nrow_global] = total_nnz;
    
    
    if(block_pt == 0)
    {
      block_pt = new CRDoubleMatrix(new_distribution_pt );
    }
    
    block_pt->build(total_ncol_global, new_values, 
                    new_column_indices, new_row_start);
    /*
    unsigned block_nrow = total_nrow_global;
    unsigned block_ncol = total_nrow_global;
    unsigned block_nnz = 0;
    double* temp_value = new double[block_nnz];
    int* temp_column_index = new int[block_nnz];
    
    int* temp_row_start = new int[block_nrow+1];
    for (unsigned i = 0; i <= block_nrow; i++)
    {
      temp_row_start[i] = 0;
    }
    
    block_pt->build_without_copy(block_ncol,block_nnz,
                                 temp_value,temp_column_index,
                                 temp_row_start);
    */
  } // void ParallelOutflowPreconditioner::merge
  
 //========================================================================
 /// add the scaled identity matrix to the specified block
 //========================================================================
  void ParallelOutflowPreconditioner::get_diag(CRDoubleMatrix *&block_pt, 
                                                       Vector<double>& diag)
  {
    // Note that diag_sqrd is nrow_local long.
    double* values = block_pt->value();
    int* column_index = block_pt->column_index();
    int* row_start = block_pt->row_start();
    unsigned nrow_local = block_pt->nrow_local();
    unsigned first_row = block_pt->first_row();
    
      for (unsigned i = 0; i < nrow_local; i++)
      {
        bool found = false;
        for (int j = row_start[i]; j < row_start[i+1] && !found; j++)
        {
          if (column_index[j] == (int)(i + first_row))
          {
            diag[i] = values[j];
            found = true;
          }
        }
        if(!found)
        {
          diag[i] = 0.0;
        }
      }
  } // end_of_get_diag(CRDoubleMatrix *&block_pt,Vector<double>& diag)
  
 //========================================================================
 /// add the scaled identity matrix to the specified block
 //========================================================================
  void ParallelOutflowPreconditioner::add_scaling_to_diag(double &Scaling, 
                                                 CRDoubleMatrix *&block_pt)
  {
    double* values = block_pt->value();
    int* column_index = block_pt->column_index();
    int* row_start = block_pt->row_start();
    unsigned nrow_local = block_pt->nrow_local();
    unsigned first_row = block_pt->first_row();
  
    // We check if there are any zero entries on the diagonal of
    // the block.
    long unsigned nzero_diag = 0;
    Vector<unsigned> diag_pres(nrow_local, 0);
    for (unsigned i = 0; i < nrow_local; i++)
    {
      bool found = false;
      for (int j = row_start[i];
           j < row_start[i+1] && !found; j++)
      {
        if (column_index[j] == (int)(i + first_row))
        {
          diag_pres[i] = 1;
          found = true;
        }
      }
      if(!found)
      {
        nzero_diag++;
      }
    }
    
    if(nzero_diag != 0)
    {
      long unsigned nnz = block_pt->nnz();
      Vector<double> new_values(nnz + nzero_diag);
      Vector<int> new_column_index(nnz + nzero_diag);
      Vector<int> new_row_start(nrow_local+1);
      long unsigned offset = 0;
      unsigned nrow_global = block_pt->nrow();
      new_row_start[0] = row_start[0];
      
      // Loop through each row
      for(unsigned i = 0; i < nrow_local; i++)
      {
        // Have we set the artificial zero yet?
        bool set_zero = false;
        
        if(diag_pres[i])
        {
          // The diagonal is present. We just copy over the 
          // column index and value.
          for(int j = row_start[i]; j < row_start[i+1]; j++)
          {
            if (column_index[j] == (int)(i + first_row))
            {
              new_column_index[j + offset] = column_index[j];
              new_values[j + offset] = values[j] + Scaling; // added scaling.
            }
            else
            {
              new_column_index[j + offset] = column_index[j];
              new_values[j + offset] = values[j];
            }
          }
        }
        else
        {
          // Diagonal is zero. There are 3 cases:
          // (1) Zero row - there are no non-zero entries on this row.
          // (2) Right entry - there are non-zero entries to the right 
          //                   of the diagonal
          // (3) Left entry only - there are entries to the left 
          //                       of the diagonal ONLY.
          
          // Case (1):
          if(row_start[i] == row_start[i+1])
          {
            // Zero row, set the artificial zero.
            new_column_index[row_start[i] + offset] = i + first_row;
            new_values[row_start[i] + offset] = Scaling; // added scaling
            offset++;
            set_zero = true;
          }
          else
          {
            // Case (2) and (3): There are some values on this row.
            // We loop through all of these values.
            for(int j = row_start[i]; j < row_start[i+1]; j++)
            {
              if(column_index[j] < (int)(i + first_row))
              {
                // First we copy all the values to the left of the diagonal.
                new_column_index[j + offset] = column_index[j];
                new_values[j + offset] = values[j];
              }
              else
              {
                // Now we are at either the diagonal or 
                // to the right of the diagonal.
                
                if(!set_zero)
                {
                  // We have reached the diagonal.
                  new_column_index[j + offset] = i + first_row;
                  new_values[j + offset] = Scaling; // added scaling
                  offset++;
                  set_zero = true;
                }
                
                // Values to the right of the diagonal.
                new_column_index[j + offset] = column_index[j];
                new_values[j + offset] = values[j];
              }
            }
            // Case (3): If there are values to the left of the diagonal only,
            // we never reach the diagonal. We set the diagonal in this case.
            if(!set_zero)
            {
              new_column_index[row_start[i+1] + offset] = i + first_row;
              new_values[row_start[i+1] + offset] = Scaling; // added scaling
              offset++;
              set_zero = true;
            }
            
          }
        }
        new_row_start[i+1] = row_start[i+1]+offset;
      }
      // This assumes that this is a square matrix
      block_pt->build(nrow_global,new_values,new_column_index,new_row_start);
    }
    else
    {
      for (unsigned i = 0; i < nrow_local; i++)
      {
        bool found = false;
        for (int j = row_start[i]; j < row_start[i+1] && !found; j++)
        {
          if (column_index[j] == (int)(i + first_row))
          {
            values[j] += Scaling;
            found = true;
          }
        }
      }
    }
  }// end_of_add_scaling_on_diag
  
 //========================================================================
 /// Setup method for the ParallelOutflowPreconditioner.
 //========================================================================
 void ParallelOutflowPreconditioner::setup(Problem* problem_pt, 
                                           DoubleMatrixBase* matrix_pt)
 {
  
  // For debugging
  //bool doc_block_matrices = false;
  bool Doc_time = false;

  // clean
  this->clean_up_memory();
  
  unsigned nmeshes = Meshes_pt.size();
//*
#ifdef PARANOID
  // paranoid check that meshes have been set
  for(unsigned mesh_i = 0; mesh_i < nmeshes; mesh_i++)
  {
    if (Meshes_pt[mesh_i]==0)
    {
      std::ostringstream error_message;
      error_message << "The Meshes_pt[" << mesh_i << "] must be set.";
      throw OomphLibError(error_message.str(),
                          "ParallelOutflowPreconditioner",
                          OOMPH_EXCEPTION_LOCATION);
    }
  }
#endif
  
  // set the mesh
  for(unsigned mesh_i = 0; mesh_i < nmeshes; mesh_i++)
  {
    this->set_mesh(mesh_i,problem_pt,Meshes_pt[mesh_i]);
  }
 

  if (this->is_master_block_preconditioner())
   {
    // get the number of fluid dof types from the first element
    N_fluid_dof_types = this->ndof_types_in_mesh(0);
    std::cout << "N_fluid_dof_types: " << N_fluid_dof_types << std::endl;
    // reset the N_lagrange_dof_types to 0
    N_lagrange_dof_types = 0;

    // the rest of the meshes are Lagrange_multiplier blocks.
    for(unsigned mesh_i = 1; mesh_i < nmeshes; mesh_i++)
    {
      N_lagrange_dof_types += this->ndof_types_in_mesh(mesh_i);
      std::cout << "mesh_i = " << mesh_i 
                << "N_lagrange+dof_types = " << N_lagrange_dof_types << std::endl;
    }
    
    // get the total number of dof types
    N_dof_types = N_fluid_dof_types + N_lagrange_dof_types;
   }
  else
   {
     // RAYEDIT - I'm not sure what to do here.
     // I'll cross the bridge when I come to it.
    N_dof_types = this->ndof_types();
    N_fluid_dof_types = 5; // rayupdate update this! (int)(((double)2*n_dof_types)/3);
   }
   
 cout << "nmeshes: " << nmeshes << endl;
  
 // We get the spatial dimension straight from the sauce!
  FiniteElement* el_pt =  
    dynamic_cast<FiniteElement*> (Meshes_pt[0]->element_pt(0));
  Dim = el_pt->dim();
  cout << "Dim: " << Dim << endl;
  // determine the number of velocity dof types
  N_velocity_dof_types = N_fluid_dof_types - 1;
  
  // Call block setup for this preconditioner
  this->block_setup(problem_pt,matrix_pt);

  // Recast Jacobian matrix to CRDoubleMatrix
  CRDoubleMatrix* cr_matrix_pt 
    = dynamic_cast<CRDoubleMatrix*>(matrix_pt);
#ifdef PARANOID
  if (cr_matrix_pt==0)
   {
    std::ostringstream error_message;
    error_message << "PrallelOutflowPreconditioner only works with"
                  << " CRDoubleMatrix matrices" << std::endl;
    throw OomphLibError(error_message.str(),
                        "ParallelOutflowPreconditioner",
                        OOMPH_EXCEPTION_LOCATION);
   }
#endif
  
  // Set up the string to record the current newton step
  std::stringstream NewtonStepStream;
  NewtonStepStream << "NS" 
                   << Global_Parameters::before_newton_step_counter;
  
  Global_Parameters::before_newton_step_counter++;
  Global_Parameters::its=0;
  
  string newton_step_counter = NewtonStepStream.str();
  string currentsetting = Global_Parameters::Current_settings
                          + newton_step_counter;

//cout << "currentsetting from setup: " << currentsetting << endl;
  
  // We begin a new newton step
  Global_Parameters::first_solve = true;
  
  // Setting up the submatrix dimension file.
  if(Global_Parameters::Dump_matrices)
  {
  ofstream submatrixdimensions;
  string subdim = "subdim_" + currentsetting;
  char *charsubdim = new char[subdim.length()+1];
  strcpy(charsubdim, subdim.c_str());
  
  //RRR_DUMP
  submatrixdimensions.open(charsubdim);
  
  //std::stringstream jacobianfile;
  //jacobianfile << "jacobian_" << currentsetting;
  
  //RRR_DUMP
  //cr_matrix_pt->sparse_indexed_output(jacobianfile.str());
  

  
  //RRR_DUMP
  //*
  for(unsigned Mi=0; Mi<N_dof_types; Mi++)
  {
    for(unsigned Mj=0; Mj<N_dof_types; Mj++)
     {  
      CRDoubleMatrix* sub_matrix_pt = 0;
      this->get_block(Mi,Mj,cr_matrix_pt,sub_matrix_pt);
      std::stringstream blockname;
      blockname << "subjacblock_"<< currentsetting<< "_" << Mi << Mj;
      int nrows = sub_matrix_pt->nrow ();
      int ncols = sub_matrix_pt->ncol ();
      //cout << blockname.str() <<": "<< nrows << "x" << ncols << endl;
      submatrixdimensions << nrows << " " << ncols << endl;
      sub_matrix_pt->sparse_indexed_output(blockname.str());
     }//for
  }//for
   submatrixdimensions.close();
   delete []charsubdim;
 // */
  } // if(Global_Parameters::Dump_matrices)

//pause("dfdsfdsfdsfdsfes");
//////////////////////////////////////////////////////////////////////////////
//////////////////////////// NEW STUFF ///////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
  

  
  ///////////////////////////////////////////////////////////////////////////
  
  // Re-order the dof_types.
  // By default the dof_types are in the order (in 2D):
  // U V P U_c V_c L
  // We re-arranging the dof_types to:
  // U U_c V V_c P L
  // by supplying the vector Dof_type_list = 0 3 1 4 2 5
  //
  // in Dof_type_list,
  // index i = 0; i < 2*Dim is the fluid dof_types.
  // 2*Dim is the pressure block.
  // i = 2*Dim+1; i < N_dof_types is the lagrange multiplier block
  
  Dof_type_list.resize(N_dof_types);
  
  // fluid dof_types
  for(unsigned i = 0; i < Dim; i++)
  {
    Dof_type_list[i*2] = i;
    Dof_type_list[i*2 + 1] = i + Dim + 1;
  }
  
  // pressure
  Dof_type_list[Dim*2] = Dim;
  
  // lagrange multipliers
  for(unsigned i = 2*Dim+1; i < N_dof_types; i++)
  {
    Dof_type_list[i] = i;
  }
  
  /* //RRR_DUMP
  cout << "Dof_type_list[i]:"<< endl;
  for(unsigned i = 0; i < N_dof_types; i++)
  {
    cout << Dof_type_list[i] << endl;
  }
  // */
 
  // Store the size of get type of block.
  // This will be used in the preconditioner_solve() function.
  // (per Newton iteration, per Newton Step) when re-arranging
  // the block vectors.
  Dof_type_block_size.assign(N_dof_types,0);
  for(unsigned i = 0; i < N_dof_types; i++)
  {
    Dof_type_block_size[i] = block_distribution_pt(Dof_type_list[i])->nrow_local();
  }
  
  Fluid_block_size = 0;
  Pressure_block_size = 0;
  Velocity_block_size = 0;
  for(unsigned i = 0; i < N_velocity_dof_types; i++)
  {
    Velocity_block_size += Dof_type_block_size[i];
  }
  
  Pressure_block_size = Dof_type_block_size[N_velocity_dof_types];

  Fluid_block_size = Velocity_block_size + Pressure_block_size;

  ///////////////////////////////////////////////////////////////////////////
/*
  cout << "Velocity_block_size" << Velocity_block_size << endl;
  cout << "Pressure_block_size" << Pressure_block_size << endl;
  cout << "Fluid_block_size" << Fluid_block_size << endl;
  pause("done done");
*/
  ////////////////////////////////////////////////////////////////////////////////
  // Need to create the norms, used for Sigma, if required

  if(Use_default_norm_of_f_scaling)
  {  
  // Create the scaling, norm of the momentum block.
  // Creating A, Ax and Ay to see the infinity norm, we re-arrange the matrix:
  // A = 0-4, Ax = 0-3, Ay = 1-4
  //    0    3     1    4     2    5
  //  0 Axx  Axxo  Axy  Axyo  Bxt  0
  //  3 Axox Axoxo Axoy Axoyo Bxot Mx
  //  1 Ayx  Ayxo  Ayy  Ayyo  Byt  0
  //  4 Ayox Ayoxo Ayoy Ayoyo Byot My
  //  2 Bx   Bxo   By   Byo   0    0
  //  5 0    Mx    0    My    0    0
  
  double ax_norm = 0.0;
  /*
  double a_norm = 0.0;
  DenseMatrix<CRDoubleMatrix* > a_pts(N_velocity_dof_types,
                                      N_velocity_dof_types,0);
  for(unsigned row_i = 0; row_i < N_velocity_dof_types; row_i++)
  {
    for(unsigned col_i = 0; col_i < N_velocity_dof_types; col_i++)
    {
      this->get_block(Dof_type_list[row_i], Dof_type_list[col_i],
                      cr_matrix_pt,a_pts(row_i,col_i));
    }
  }
  get_inf_norm(a_pts,a_norm);
  // */
  
  double t_norm_ax_start = TimingHelpers::timer();
  // Get the momentum block corresponding to the x or u block.
  DenseMatrix<CRDoubleMatrix* > ax_pts(2,2,0);
  for(unsigned row_i = 0; row_i < 2; row_i++)
  {
    for(unsigned col_i = 0; col_i < 2; col_i++)
    {
      this->get_block(Dof_type_list[row_i], Dof_type_list[col_i],
                      cr_matrix_pt,ax_pts(row_i,col_i));
    }
  }


  // Get the norm
  get_inf_norm(ax_pts,ax_norm);
  double t_norm_ax_finish = TimingHelpers::timer();
  
  if(Doc_time)
  {
    double t_norm_ax_time = t_norm_ax_finish - t_norm_ax_start;
    cout << "t_norm_ax_time: " << t_norm_ax_time << std::endl;
  }
  
  // Set Scaling_sigma
  Scaling_sigma =  -ax_norm;
  streamsize cout_precision = cout.precision();
  cout << "RAYNORM: " << currentsetting << " " << setprecision(15)
       << ax_norm << setprecision(cout_precision) << endl;
  } // if(Use_default_f_scaling)
  
#ifdef PARANOID
  if(Scaling_sigma == 0.0)
  {
    std::ostringstream warning_stream;
    warning_stream << "WARNING: " << std::endl
                   << "The scaling (Scaling_sigma) is " << Scaling_sigma << std::endl
                   << "Division by Scaling_sigma = 0 will implode the world."
                   << std::endl;
    OomphLibWarning(warning_stream.str(),
                    "ParallelOutflowPreconditioner::setup()",
                    OOMPH_EXCEPTION_LOCATION);
  }
  if(Scaling_sigma > 0.0)
  {
    std::ostringstream warning_stream;
    warning_stream << "WARNING: " << std::endl
                   << "The scaling (Scaling_sigma) is " << Scaling_sigma << std::endl
                   << "Performance may be degraded."
                   << std::endl;
    OomphLibWarning(warning_stream.str(),
                    "ParallelOutflowPreconditioner::setup()",
                    OOMPH_EXCEPTION_LOCATION);
  }
#endif  

  streamsize cout_precision = cout.precision();
  cout << "RAYSIGMA: " << setprecision(15) << Scaling_sigma 
                       << setprecision(cout_precision) 
                       << endl;

////////////////////////////////////////////////////////////////////////////////

  // We will now extract all of the mass matrices
  // which corresponds to the L block in the Jacobian
  // Create the d(W) block and then the augmentation block.

  // Extract all the mass matrices by
  // looping through the lagrange multipliers
  DenseMatrix<CRDoubleMatrix* > mm_pointers(N_lagrange_dof_types,Dim,0);
  double t_get_massmatrices_start = TimingHelpers::timer();
  // loop through the number of Lagrange multiplier rows
  for(unsigned l_i = 0; l_i < N_lagrange_dof_types; l_i++)
  {
    // The corresponding block row in the Jacobian
    unsigned row_i = N_fluid_dof_types + l_i;
    
    // loop through the constrained columns (same as number of spatial dim).
    for(unsigned d_i = 0; d_i < Dim; d_i++)
    {
      // The corresponding block column in the Jacobian
      unsigned col_i = 2*d_i+1;
      this->get_block(Dof_type_list[row_i], Dof_type_list[col_i],
                      cr_matrix_pt,mm_pointers(l_i,d_i));
    }
  }

  double t_get_massmatrices_finish = TimingHelpers::timer();
  if(Doc_time)
    {
      double t_get_massmatrices_time = t_get_massmatrices_finish 
                                       - t_get_massmatrices_start;
      cout << "t_get_massmatrices_time: " << t_get_massmatrices_time << std::endl;
    }
  /* //RRR_DUMP
  //CRDoubleMatrix* l_pt = 0;
  //merge(mm_pointers,l_pt);
  //l_pt->sparse_indexed_output("merged_l_block");
  // */
  
////////////////////////////////////////////////////////////////////////////////
  // We new create all of the d(W) sub-blocks d(W_i)
  // Then merge them into final d(W), or we can just solve each block!
  double t_create_dW_start = TimingHelpers::timer();
  // The number of W_i is the same as the number of Lagrange mulipliers.
  DenseMatrix<CRDoubleMatrix* > w_pointers(1,N_lagrange_dof_types,0);
  DenseMatrix<CRDoubleMatrix* > invw_pointers(1,N_lagrange_dof_types,0);
  // loop through the number of lagrange multipliers
  for(unsigned l_i = 0; l_i < N_lagrange_dof_types; l_i++)
  { 
    // Get the number of local rows for this lagrange block.
    // We shall use the block in the first column.
    unsigned long l_i_nrow_local = mm_pointers(l_i,0)->nrow_local();
    // A Vector of Vectors containing the diagonals of 
    // each mass matrix of this Lagrange mult. constraint.
    Vector<Vector<double> > m_diag(Dim,Vector<double>(l_i_nrow_local,0.0));
    
    // Extract all the diagonals of the mass matrices for this
    // block lagrange row (l_i)
    for(unsigned d_i = 0; d_i < Dim; d_i++)
    {
      get_diag(mm_pointers(l_i,d_i),m_diag[d_i]);
    }
    
    // A vector to contain the results of mass matrices squared.
    Vector<double> w_i_diag_values(l_i_nrow_local,0);
    Vector<double> invw_i_diag_values(l_i_nrow_local,0);
    
    Vector<int> w_i_column_indices(l_i_nrow_local);
    Vector<int> w_i_row_start(l_i_nrow_local+1);
  
    // Component-wise, square and add all the diagonals.
    for(unsigned d_i = 0; d_i < Dim; d_i++)
    {
      for(unsigned long row_i = 0; row_i < l_i_nrow_local; row_i++)
      {
        w_i_diag_values[row_i] += m_diag[d_i][row_i]*m_diag[d_i][row_i];
      }
    }
    
    // Divide by Scaling_sigma and create the inverse of w.
    for(unsigned long row_i = 0; row_i < l_i_nrow_local; row_i++)
    {
      w_i_diag_values[row_i] /= Scaling_sigma;

      // w_i is a diagonal matrix, so take the inverse to
      // invert the matrix.
      invw_i_diag_values[row_i] = 1/w_i_diag_values[row_i];
      w_i_column_indices[row_i] = row_i;
      w_i_row_start[row_i] = row_i;
    }
    w_i_row_start[l_i_nrow_local] = l_i_nrow_local;
    
    // need to put the thing in here:
    bool distributed = this->master_distribution_pt()->distributed();
    // Get
    
    unsigned long l_i_nrow_global = mm_pointers(l_i,0)->nrow();
    
    LinearAlgebraDistribution* new_distribution_pt 
      = new LinearAlgebraDistribution(problem_pt->communicator_pt(),
                                      l_i_nrow_global,distributed);
    
    w_pointers(0,l_i) = new CRDoubleMatrix(new_distribution_pt);
    invw_pointers(0,l_i) = new CRDoubleMatrix(new_distribution_pt);
    
    // Theses are square matrices
    w_pointers(0,l_i)->build(l_i_nrow_global,
                             w_i_diag_values,
                             w_i_column_indices,
                             w_i_row_start);
    
    invw_pointers(0,l_i)->build(l_i_nrow_global,
                                invw_i_diag_values,
                                w_i_column_indices,
                                w_i_row_start);
  
  } // Loop through all of the lagrange multipliers.
    // for(unsigned l_i = 0; l_i < N_lagrange_dof_types; l_i++)
  double t_create_dW_finish = TimingHelpers::timer();
  if(Doc_time)
    {
      double t_create_dW_time = t_create_dW_finish - t_create_dW_start;
      cout << "t_create_dW_time: " << t_create_dW_time << std::endl;
    }
////////////////////////////////////////////////////////////////////////////////
// CREATING THE AUGMENTED F BLOCK

  // What we have done already:
  // Dof_type_list(N_dof_types): constains the permutations
  // mm_pointers(N_lagrange_dof_types,Dim): the mm from L block
  // w_pointers(1,N_lagrange_dof_types): the w_i's
  // invw_pointers(1,N_lagrange_dof_types): inverse of above

 // Note that we shall use the re-arranged order of the block dof types.
 // i.e. we have U Uc V Vc W Wc P L1, etc... 
    
  // Setup the fluid subsidiary precoditioner
  if(L_prec_type == Exact_block_preconditioner)
  {
  
  // Get the fluid (Navier-Stokes) block matrix.
  double t_get_F_start = TimingHelpers::timer();
  DenseMatrix<CRDoubleMatrix* > f_aug_ptrs(N_fluid_dof_types,
                                           N_fluid_dof_types,0);
  // Extract the fluid blocks (including pressure)
  for(unsigned row_i = 0; row_i < N_fluid_dof_types; row_i++)
  {
    for(unsigned col_i = 0; col_i < N_fluid_dof_types; col_i++)
    {
      this->get_block(Dof_type_list[row_i], Dof_type_list[col_i],
                      cr_matrix_pt,f_aug_ptrs(row_i,col_i));
    }//for 
  }//for
  
  double t_get_F_finish = TimingHelpers::timer();
  if(Doc_time)
    {
      double t_get_F_time = t_get_F_finish - t_get_F_start;
      cout << "t_get_F_time: " << t_get_F_time << std::endl;
    }//if
  
  // Create and add the augmentation.
  // The constrained blocks are located at: 2*i + 1, 
  // using the new Dof_type_list, where i is 
  // the index for spatial dimension.
  // Outter most loop: Go through the lagrange multipliers
  double t_create_add_aug_start = TimingHelpers::timer();
  // There are l_i blocks to add.
  for(unsigned l_i = 0; l_i < N_lagrange_dof_types; l_i++)
  {
    // loop through the constrained rows
    for(unsigned row_i = 0; row_i < Dim; row_i++)
    {
      // Dim = 3 -> row_i  = 0,1,2
      // constrained_row_i = 1,3,5 (re-arranged dofs)
      unsigned constrained_row_i = 2*row_i + 1; 

      // loop through the constrained columns
      for(unsigned col_i = 0; col_i < Dim; col_i++)
      {
        // Dim = 3 -> col_i  = 0,1,2
        // constrained_col_i = 1,3,5
        unsigned constrained_col_i = 2*col_i + 1; 
        
        // A temp pointer to store the intermediate results.
        CRDoubleMatrix* aug_pt = 0;

        // being lazy, need to fix this... I am extracting the block
        // unnecessarily. This is because the multiply method requires
        // a matrix, not just a pointer.
        this->get_block(Dof_type_list[constrained_row_i],
                        Dof_type_list[constrained_col_i],
                        cr_matrix_pt,aug_pt);    
        
        mm_pointers(l_i,row_i)->multiply((*invw_pointers(0,l_i)),(*aug_pt));
        aug_pt->multiply((*mm_pointers(l_i,col_i)),(*aug_pt));
        

        add_matrices(aug_pt,f_aug_ptrs(constrained_row_i,
                                       constrained_col_i));
        
      } // loop through the constrained columns
    } // loop through the constrained rows
  } // loop through the lagrange multipliers
  
  double t_create_add_aug_finish = TimingHelpers::timer();
  if(Doc_time)
    {
      double t_create_add_aug_time = t_create_add_aug_finish 
                                     - t_create_add_aug_start;
      std::cout << "t_create_add_aug_time: " << t_create_add_aug_time << std::endl;
    }//if




  double t_merge_F_aug_start = TimingHelpers::timer();
  CRDoubleMatrix* f_aug_pt = 0;

  merge(f_aug_ptrs,f_aug_pt);
  double t_merge_F_aug_finish = TimingHelpers::timer();
  if(Doc_time)
    {
      double t_merge_F_aug_time = t_merge_F_aug_finish-t_merge_F_aug_start;
      cout << "t_merge_F_aug_time: " << t_merge_F_aug_time << std::endl;
    }

  // delete the sub F pointers
  for(unsigned row_i = 0; row_i < N_velocity_dof_types; row_i++)
  {
    for(unsigned col_i = 0; col_i < N_velocity_dof_types; col_i++)
    {
      delete f_aug_ptrs(row_i, col_i);
      f_aug_ptrs(row_i, col_i) = 0;
    }
  }
  NS_preconditioner_pt = new SuperLUPreconditioner;
  NS_preconditioner_pt->setup(problem_pt,f_aug_pt);

  delete f_aug_pt;
  f_aug_pt = 0;
  } // if(L_prec_type == Exact_block_preconditioner)
  else if(L_prec_type == Exact_lsc_block_preconditioner)
  {
    // We solve the pressure block first, then velocity block.
  // Get B (the divergence block)
  // This is done in three steps. 
  // 1) Extract all the sub-blocks
  // 2) Merge the sub-blocks
  // 3) Delete the pointers to the sub-blocks.
  //
  double t_get_B_start = TimingHelpers::timer();

  // This has to be a matrix, not a vector, because the merge function
  // requires a DenseMatrix. This is so we can generalise it...
  DenseMatrix<CRDoubleMatrix* > b_pointers(1,N_velocity_dof_types,0);
  
  // Encapsulation of the variable row_i
  {
    unsigned row_i = N_velocity_dof_types;
    for(unsigned col_i = 0; col_i < N_velocity_dof_types; col_i++)
    {
      this->get_block(Dof_type_list[row_i], Dof_type_list[col_i],
                      cr_matrix_pt,b_pointers(0,col_i));
    }//for(unsigned col_i = 0; col_i < N_velocity_dof_types; col_i++)
  }
  
  
  double t_get_B_finish = TimingHelpers::timer();
  
  if(Doc_time)
    {
      double t_get_B_time = t_get_B_finish - t_get_B_start;
      cout << "t_get_B_time: " << t_get_B_time << "\n";
    }

  // merge the sub-blocks.
  double t_merge_B_start = TimingHelpers::timer();
  CRDoubleMatrix* b_pt = 0;
  merge(b_pointers,b_pt);
  double t_merge_B_finish = TimingHelpers::timer();
  if(Doc_time)
    {
      double t_merge_B_time = t_merge_B_finish - t_merge_B_start;
      cout << "t_merge_B_time: " << t_merge_B_time << "\n";
    }
  
  // delete the sub-blocks
  for(unsigned col_i = 0; col_i < N_velocity_dof_types; col_i++)
  {
    delete b_pointers(0,col_i);
    b_pointers(0,col_i) = 0;
  }


  
////////////////////////////////////////////////////////////////////////////////
  
  // Compute and merge the velocity mass matrix.
  // The velocity dof types are located at from 0 to 2*Dim-1 in Dof_type_list.
  
  // Store the mass matrices in here:
  DenseMatrix<CRDoubleMatrix* > inv_v_mass_ptrs(1,N_velocity_dof_types,0);
  
  // get the inverse velocity and pressure mass matrices
  CRDoubleMatrix* inv_v_mass_pt = 0;
  CRDoubleMatrix* inv_p_mass_pt = 0; // This is not required for LSC.
  
  // Extract all of the inv_v_mass.
  bool do_both=false;
  double t_get_ivmm_start = TimingHelpers::timer();
  for(unsigned dof_type_i = 0; 
      dof_type_i < N_velocity_dof_types; dof_type_i++)
  {
    unsigned required_block = Dof_type_list[dof_type_i];

    assemble_inv_press_and_veloc_mass_matrix_diagonal
      (inv_p_mass_pt, inv_v_mass_ptrs(0, dof_type_i), 
       do_both, required_block);

    if(Global_Parameters::Dump_matrices)
    {
    // RRR_DUMP
    std::stringstream blockname;
    blockname << "submassblock_"<< currentsetting<< "_" 
              << required_block << required_block;
    inv_v_mass_ptrs(0,dof_type_i)->sparse_indexed_output(blockname.str());
    }
  }

  double t_get_ivmm_finish = TimingHelpers::timer();
  
  if(Doc_time)
    {
      double t_get_ivmm_time = t_get_ivmm_finish - t_get_ivmm_start;
      std::cout << "t_get_ivmm_time: " << t_get_ivmm_time << "\n";
    }
  
  /*
  inv_v_mass_ptrs(0,0)->sparse_indexed_output("submassblock_Noel3_00_t");
  inv_v_mass_ptrs(0,1)->sparse_indexed_output("submassblock_Noel3_33_t");
  inv_v_mass_ptrs(0,2)->sparse_indexed_output("submassblock_Noel3_11_t");
  inv_v_mass_ptrs(0,3)->sparse_indexed_output("submassblock_Noel3_44_t");
  */
  
  // Now merge the matrices in inv_v_mass_ptrs(1,N_velocity_dof_types) 
  // along the diagonal.
  
  // Get the total number of rows.
  // loop through all of the rows and get the values. 
  //   Set the row_start and column.
  double t_merge_ivmm_start = TimingHelpers::timer();
  {
    bool distributed = this->master_distribution_pt()->distributed();
    
    unsigned long total_nrow_global = 0;
    for(unsigned dof_type_i = 0; 
        dof_type_i < N_velocity_dof_types; dof_type_i++)
    {
      total_nrow_global += inv_v_mass_ptrs(0,dof_type_i)->nrow();
      //cout << "nrow = " << inv_v_mass_ptrs(0,dof_type_i)->nrow() << endl;
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
    
    // Loop through the inv_v_mass_ptrs and populate the
    // new_values, new_column_indices and new_row_start
    unsigned long row_i = 0;
    for(unsigned dof_type_i = 0; 
        dof_type_i < N_velocity_dof_types; dof_type_i++)
    {
      CRDoubleMatrix* current_block_pt = inv_v_mass_ptrs(0,dof_type_i);
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
    
    //inv_v_mass_pt->sparse_indexed_output("inv_v_mass_pt");
  }

  double t_merge_ivmm_finish = TimingHelpers::timer();
  
  if(Doc_time)
    {
      double t_merge_ivmm_time = t_merge_ivmm_finish - t_merge_ivmm_start;
      cout << "t_merge_ivmm_time: " << t_merge_ivmm_time << "\n";
    }
  
  // Delete the pointers to the sub matrices.
  for(unsigned mass_i = 0; mass_i < N_velocity_dof_types; mass_i++)
  {
    delete inv_v_mass_ptrs(0,mass_i);
    inv_v_mass_ptrs(0,mass_i) = 0;
  }

////////////////////////////////////////////////////////////////////////////////
  // Get Bt (the discrete something block)
  DenseMatrix<CRDoubleMatrix* > bt_pointers(N_velocity_dof_types,1,0);
  
  double t_get_Bt_start = TimingHelpers::timer();
  {
    unsigned col_i = N_velocity_dof_types;
    for(unsigned row_i = 0; row_i < N_velocity_dof_types; row_i++)
    {
      this->get_block(Dof_type_list[row_i], Dof_type_list[col_i],
                      cr_matrix_pt,bt_pointers(row_i,0));
    }//for(unsigned row_i = 0; row_i < N_velocity_dof_types; row_i++)
  }
  
  double t_get_Bt_finish = TimingHelpers::timer();
  
  if(Doc_time)
    {
      double t_get_Bt_time = t_get_Bt_finish - t_get_Bt_start;
      cout << "t_get_Bt_time: " << t_get_Bt_time << "\n";
    }
  
  CRDoubleMatrix* bt_pt = 0;
  
  double t_merge_Bt_start = TimingHelpers::timer();
  merge(bt_pointers,bt_pt);
  double t_merge_Bt_finish = TimingHelpers::timer();
  
  
  if(Doc_time)
    {
      double t_merge_Bt_time = t_merge_Bt_finish - t_merge_Bt_start;
      cout << "t_merge_Bt_time: " << t_merge_Bt_time << "\n";
    }

  // Delete sub-blocks
  for(unsigned row_i = 0; row_i < N_velocity_dof_types; row_i++)
  {
    delete bt_pointers(row_i,0);
    bt_pointers(row_i,0) = 0;
  }

////////////////////////////////////////////////////////////////////////////////
  // Build the pressure poisson matrix
  CRDoubleMatrix* p_matrix_pt = new CRDoubleMatrix;
  
  // Multiply the inverse velocity mass matrix by the gradient matrix B^T
  double t_QBt_matrix_start = TimingHelpers::timer();
  CRDoubleMatrix* qbt_pt = new CRDoubleMatrix;
  inv_v_mass_pt->multiply(*bt_pt, *qbt_pt);
  //cout << "inv_v_mass_pt = " 
  //     << inv_v_mass_pt->nrow() << ", " 
  //     << inv_v_mass_pt->ncol() << endl;
  //
  // NOTE: qbt_pt = Q^{-1} * Bt
  // RRR_DUMP - passed
  //qbt_pt->sparse_indexed_output("qbt_pt");

  delete bt_pt;
  bt_pt = 0;

  // Store the product in bt_pt
  bt_pt = qbt_pt;
  double t_QBt_matrix_finish = TimingHelpers::timer();
  if(Doc_time)
   {
    double t_mult_QBt_time = t_QBt_matrix_finish - t_QBt_matrix_start;
    cout << "t_mult_QBt_time: "
               << t_mult_QBt_time << std::endl;
   }
  
  // delete the inverse of the velocity mass matrix.
  // when working out B* Q^-1, we can transpose qbt_pt
  delete inv_v_mass_pt;
  inv_v_mass_pt = 0;

  // Multiply Q^-1 * B^T from the left by the divergence matrix B and store
  // the result in the pressure poisson matrix.
  double t_p_matrix_start = TimingHelpers::timer();
  b_pt->multiply(*bt_pt,*p_matrix_pt);
  double t_p_matrix_finish = TimingHelpers::timer();
  // RRR_DUMP - passed
  // p_matrix_pt->sparse_indexed_output("p_matrix_pt");
  // NOTE: p_matrix_pt = B * Q^-1 * B^T
    
  if(Doc_time)
   {
    double t_mult_p_time = t_p_matrix_finish - t_p_matrix_start;
    cout << "t_mult_p_time: "
               << t_mult_p_time << std::endl;
   }
    
  // Kill divergence matrix because we don't need it any more
  delete b_pt;

  // Build the matvec operator for QBt
  // NOTE: bt_pt = Q^-1 * B^T
  double t_QBt_MV_start = TimingHelpers::timer();
  QBt_mat_vec_pt = new MatrixVectorProduct;
  QBt_mat_vec_pt->setup(bt_pt);
  double t_QBt_MV_finish = TimingHelpers::timer();
  if(Doc_time)
   {
    double t_MV_QBt_time = t_QBt_MV_finish - t_QBt_MV_start;
    oomph_info << "t_MV_QBt_time: "
               << t_MV_QBt_time << std::endl;
   }
  
  delete bt_pt;
  bt_pt = 0;

////////////////////////////////////////////////////////////////////////////////
// CREATING THE AUGMENTED F BLOCK

  // We now have:
  // mass matrices in: mm_pointers(N_lagrange_dof_types,Dim);
  // dw in: w_pointers(1,N_lagrange_dof_types)
  // inv(dw) in: invw_pointers(1,N_lagrange_dof_types)
  
  // We now form the augmented block:
  
  // Get the momentum matrix.
  // This does not include the pressure block.
  double t_get_F_start = TimingHelpers::timer();
  DenseMatrix<CRDoubleMatrix* > f_aug_ptrs(N_velocity_dof_types,
                                           N_velocity_dof_types,0);
  
  for(unsigned row_i = 0; row_i < N_velocity_dof_types; row_i++)
  {
    for(unsigned col_i = 0; col_i < N_velocity_dof_types; col_i++)
    {
      this->get_block(Dof_type_list[row_i], Dof_type_list[col_i],
                      cr_matrix_pt,f_aug_ptrs(row_i,col_i));
    }
  }
  double t_get_F_finish = TimingHelpers::timer();
  if(Doc_time)
    {
      double t_get_F_time = t_get_F_finish - t_get_F_start;
      cout << "t_get_F_time: " << t_get_F_time << std::endl;
    }
  
  // Add the augmentation.
  // The constrained blocks are located at: 2*i + 1
  // loop through the lagrange multipliers
  double t_create_add_aug_start = TimingHelpers::timer();
  for(unsigned l_i = 0; l_i < N_lagrange_dof_types; l_i++)
  {
    // loop through the constrained rows
    for(unsigned row_i = 0; row_i < Dim; row_i++)
    {
      //
      //unsigned row_i = 1;
      //
      unsigned constrained_row_i = 2*row_i + 1;
      // loop through the constrained columns
      for(unsigned col_i = 0; col_i < Dim; col_i++)
      {
        //
        //unsigned col_i = 1;
        //
        unsigned constrained_col_i = 2*col_i + 1;
        CRDoubleMatrix* aug_pt = 0;
        // being lazy, need to fix this... I am extracting the block
        // unnecessarily.
        
        //cout << "Dof_type_list[constrained_row_i] = " << Dof_type_list[constrained_row_i]
        //     << "Dof_type_list[constrained_col_i] = " << Dof_type_list[constrained_col_i]
        //     << " " << endl;
        this->get_block(Dof_type_list[constrained_row_i],
                        Dof_type_list[constrained_col_i],
                        cr_matrix_pt,
                        aug_pt);
        
        mm_pointers(l_i,row_i)->multiply((*invw_pointers(0,l_i)),(*aug_pt));
        aug_pt->multiply((*mm_pointers(l_i,col_i)),(*aug_pt));
        
        //aug_pt->sparse_indexed_output("aug_pt");

        add_matrices(aug_pt,f_aug_ptrs(constrained_row_i, constrained_col_i));
        
        //f_aug_ptrs(constrained_row_i, constrained_col_i)
        //  ->sparse_indexed_output("f_aug_pt");
        
        //pause("newtesto");
        //cout << "l_i = " << l_i << ", i = " << constrained_row_i 
        //<< ". j = " << constrained_col_i << endl;
      } // loop through the constrained columns
    } // loop through the constrained rows
  } // loop through the lagrange multipliers
  double t_create_add_aug_finish = TimingHelpers::timer();
  if(Doc_time)
    {
      double t_create_add_aug_time = t_create_add_aug_finish 
                                     - t_create_add_aug_start;
      std::cout << "t_create_add_aug_time: " << t_create_add_aug_time << std::endl;
    }
  // f_aug_ptrs(N_velocity_dof_types,N_velocity_dof_types) has now been augmented.
  /*
  f_aug_ptrs(1,1)->sparse_indexed_output("f_aug_33");
  f_aug_ptrs(1,3)->sparse_indexed_output("f_aug_34");
  f_aug_ptrs(3,1)->sparse_indexed_output("f_aug_43");
  f_aug_ptrs(3,3)->sparse_indexed_output("f_aug_44");
  */

  // delete the mm_pointers
  // loop through the Lagrange multiplier blocks (constrained rows)
  for(unsigned n_i = 0; n_i < N_lagrange_dof_types; n_i++)
  {
    // loop through the columns (no. of dimensions)
    for(unsigned d_i = 0; d_i < Dim; d_i++)
    {
      delete mm_pointers(n_i,d_i);
      mm_pointers(n_i,d_i) = 0;   
    }
  }


  double t_merge_F_aug_start = TimingHelpers::timer();
  CRDoubleMatrix* f_aug_pt = 0;
  merge(f_aug_ptrs,f_aug_pt);
  double t_merge_F_aug_finish = TimingHelpers::timer();
  if(Doc_time)
    {
      double t_merge_F_aug_time = t_merge_F_aug_finish-t_merge_F_aug_start;
      cout << "t_merge_F_aug_time: " << t_merge_F_aug_time << std::endl;
    }
  
  // delete the sub F pointers
  for(unsigned row_i = 0; row_i < N_velocity_dof_types; row_i++)
  {
    for(unsigned col_i = 0; col_i < N_velocity_dof_types; col_i++)
    {
      delete f_aug_ptrs(row_i, col_i);
      f_aug_ptrs(row_i, col_i) = 0;
    }
  }

  // RRR_DUMP
  //f_aug_pt-> sparse_indexed_output("f_aug_pt");
  //pause("test f_aug_pt");
  
  double t_F_MV_start = TimingHelpers::timer();
  F_mat_vec_pt = new MatrixVectorProduct;
  F_mat_vec_pt->setup(f_aug_pt);
  double t_F_MV_finish = TimingHelpers::timer();
  if(Doc_time)
   {
    double t_F_MV_time = t_F_MV_finish - t_F_MV_start;
    oomph_info << "t_F_MV_time: "
               << t_F_MV_time << std::endl;
   }

  // Rebuild Bt
  // recall that we have deleted the Bt sub matrices.
  // So we need to extract them, merge, then delete them.
  //*
  t_merge_Bt_start = TimingHelpers::timer();
  bt_pt = 0;
  
  // re-extract the sub Bt blocks.
  {
    unsigned col_i = N_velocity_dof_types;
    for(unsigned row_i = 0; row_i < N_velocity_dof_types; row_i++)
    {
      this->get_block(Dof_type_list[row_i], Dof_type_list[col_i],
                      cr_matrix_pt,bt_pointers(row_i,0));
    }
  }

  // merge
  merge(bt_pointers,bt_pt);
  t_merge_Bt_finish = TimingHelpers::timer();
  if(Doc_time)
    {
      double t_merge_Bt_time = t_merge_Bt_finish - t_merge_Bt_start;
      cout << "t_merge_Bt_time: " << t_merge_Bt_time << std::endl;
    }
  
  // delete the sub Bt blocks again
  for(unsigned row_i = 0; row_i < N_velocity_dof_types; row_i++)
  {
    delete bt_pointers(row_i,0);
    bt_pointers(row_i,0) = 0;
  } 

  // */
  // form the matrix vector operator for Bt
  //*
  double t_Bt_MV_start = TimingHelpers::timer();
  Bt_mat_vec_pt = new MatrixVectorProduct;
  Bt_mat_vec_pt->setup(bt_pt);
  double t_Bt_MV_finish = TimingHelpers::timer();
  if(Doc_time)
   {
    double t_Bt_MV_time = t_Bt_MV_finish - t_Bt_MV_start;
    oomph_info << "t_Bt_MV_time: "
               << t_Bt_MV_time << std::endl;
   }
  delete bt_pt;
  //*/

  double t_p_prec_start = TimingHelpers::timer();
  P_preconditioner_pt = new SuperLUPreconditioner;
  //cout << "NROW AND NCOL = " << p_matrix_pt->nrow() << " " 
  //     << p_matrix_pt->ncol() << endl;
  //p_matrix_pt->sparse_indexed_output("p_matrix_pt2");
  P_preconditioner_pt->setup(problem_pt, p_matrix_pt);
  double t_p_prec_finish = TimingHelpers::timer();
  if(Doc_time)
   {
    double t_p_prec_time = t_p_prec_finish - t_p_prec_start;
    cout << "t_p_prec_time: "
               << t_p_prec_time << "\n";
   }


  // Solver for the momentum matrix.
  double t_f_prec_start = TimingHelpers::timer();
  F_preconditioner_pt = new SuperLUPreconditioner;
  F_preconditioner_pt->setup(problem_pt,f_aug_pt);
  delete f_aug_pt;
  double t_f_prec_finish = TimingHelpers::timer();
  if(Doc_time)
   {
    double t_f_prec_time = t_f_prec_finish - t_f_prec_start;
    cout << "t_f_prec_time: "
               << t_f_prec_time << "\n";
   }



  } // if else (L_prec_type == Exact_lsc_block_preconditioner)
 
  // Solver for the W block.
  double t_w_prec_start = TimingHelpers::timer();
  W_preconditioner_pt.resize(N_lagrange_dof_types);
  for(unsigned l_i = 0; l_i < N_lagrange_dof_types; l_i++)
  {
    W_preconditioner_pt[l_i] = new SuperLUPreconditioner;
    W_preconditioner_pt[l_i]->setup(problem_pt,w_pointers(0,l_i));
  }
  double t_w_prec_finish = TimingHelpers::timer();
  if(Doc_time)
   {
    double t_w_prec_time = t_w_prec_finish - t_w_prec_start;
    cout << "t_w_prec_time: "
               << t_w_prec_time << "\n";
   }


 
 
 } // end of ParallelOutflowPreconditioner::setup
 

 //========================================================================
 /// \short Clears the memory.
 //========================================================================
 void ParallelOutflowPreconditioner::clean_up_memory()
 {
  // clean the block preconditioner base class memory
  this->clear_block_preconditioner_base();  
 } // end of ParallelOutflowPreconditioner::clean_up_memory



/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////

//=========================================================================
/// Wrapper class 
//========================================================================
template <class ELEMENT>
class BulkFluidSubjectToLagrangeMultiplierElement : public virtual ELEMENT
{
 
public:
 
 /// Default constructor
 BulkFluidSubjectToLagrangeMultiplierElement() : ELEMENT() {}
 
 /// \short Returns the number of DOF types associated with this element:
 /// Twice the number of its spatial dimension plus one (for the pressure)
 unsigned ndof_types()
  {
   return 2*ELEMENT::dim()+1;
  }
 
 /// \short Create a list of pairs for all unknowns in this element,
 /// so that the first entry in each pair contains the global equation
 /// number of the unknown, while the second one contains the number
 /// of the "DOF" that this unknown is associated with.
 /// (Function can obviously only be called if the equation numbering
 /// scheme has been set up.)\n
 /// E.g. in a 2D problem there are 5 types of DOF:\n
 /// 0 - x velocity (without lagr mult )\n
 /// 1 - y velocity (without lagr mult )\n
 /// 2 - pressure
 /// 3 - x velocity (with lagr mult )\n
 /// 4 - y velocity (with lagr mult )\n
 void get_dof_numbers_for_unknowns(
  std::list<std::pair<unsigned long,unsigned> >& block_lookup_list)
  {
   // temporary pair (used to store block lookup prior to being added to list
   std::pair<unsigned,unsigned> block_lookup;
   
   // number of nodes
   const unsigned n_node = this->nnode();
   
   //Get the dimension of the node
   const unsigned nodal_dim = ELEMENT::nodal_dimension();
   
   //Integer storage for local unknown
   int local_unknown=0;
   
   //Loop over the nodes
   for(unsigned n=0;n<n_node;n++)
    {
     // Check if node has been resized, i.e. if Lagrange multipliers
     // have been attached; in that case the veloc dofs get
     // a different dof type.
     unsigned offset = 0;
     if (this->node_pt(n)->nvalue() != this->required_nvalue(n))
      {
       offset = ELEMENT::dim()+1;
      } // if
     
     //Loop over dimension for velocity components
     for(unsigned i=0;i<nodal_dim;i++)
      {
       //If the variable is free
       local_unknown = ELEMENT::nodal_local_eqn(n,i);
       // ignore pinned values
       if (local_unknown >= 0)
        {
         // store block lookup in temporary pair: First entry in pair
         // is global equation number; second entry is block type
         block_lookup.first = this->eqn_number(local_unknown);
         block_lookup.second = offset+i;
         
         // add to list
         block_lookup_list.push_front(block_lookup);  
         if (this->node_pt(n)->nvalue() != this->required_nvalue(n))
         {
//         cout << "g_eqn: " << block_lookup.first << ", type: " 
//               << block_lookup.second << std::endl;
         }
        } // if local_unknown >= 0
      } // for i - loop over dimensions
       
       // Pressure (Taylor Hood only!)
       if (this->required_nvalue(n)==(ELEMENT::dim()+1))
        {         
         //If the variable is free
         local_unknown = ELEMENT::nodal_local_eqn(n,ELEMENT::dim());
         
         // ignore pinned values
         if (local_unknown >= 0)
          {
           // store block lookup in temporary pair: First entry in pair
           // is global equation number; second entry is block type
           block_lookup.first = this->eqn_number(local_unknown);
           block_lookup.second = ELEMENT::dim();
           
           // add to list
           block_lookup_list.push_front(block_lookup);    
         if (this->node_pt(n)->nvalue() != this->required_nvalue(n))
         {
            cout << "Eqn: " << block_lookup.first << " classed as " 
                 << block_lookup.second << std::endl;
         }
          } // if local_unknown >= 0
        } // if this->required_nvalue(n)==(ELEMENT::dim()+1) checking pressure

    } // for loop over node n
  } // end of function get_dof_numbers_for unknowns

}; // end of BulkFluidSubjectToLagrangeMultiplierElement
   
//===========start_face_geometry==============================================
/// FaceGeometry of wrapped element is the same as the underlying element
//============================================================================
template<class ELEMENT>
class FaceGeometry<BulkFluidSubjectToLagrangeMultiplierElement<ELEMENT> > :
 public virtual FaceGeometry<ELEMENT>
{
};





//========================================================================
/// \short A Sloping Mesh  class.
///
/// derived from SimpleCubicMesh:
/// the same mesh rotated with angles ang_x, ang_y, and ang_z
/// about the x-axis, y-axis and z-axis respectively.
//========================================================================
 template<class ELEMENT> 
 class SlopingCubicMesh : public SimpleCubicMesh<ELEMENT>
 {
 public:

  /// Constructor.
  SlopingCubicMesh(const unsigned& nx, const unsigned& ny, const unsigned& nz,
                   const double& lx,  const double& ly, const double& lz,
                   const double& ang_x, 
                   const double& ang_y, 
                   const double& ang_z) :
   SimpleCubicMesh<ELEMENT>(nx,ny,nz,lx,ly,lz)
   {
    // Find out how many nodes there are
    unsigned n_node=this->nnode();

    // Loop over all nodes
    for (unsigned n=0;n<n_node;n++)
     {
      // Pointer to node:
      Node* nod_pt=this->node_pt(n);

      // Get the x/y coordinates
      double x=nod_pt->x(0);
      double y=nod_pt->x(1);
      double z=nod_pt->x(2);

      // Set new nodal coordinates by premultiplying by R_xyz.
      // R_xyz = R_x * R_y * R_z where R_x is the standard
      // "counter-clockwise" in three dimensions.
      
      nod_pt->x(0)=x*cos(ang_y)*cos(ang_z)
       -y*cos(ang_y)*sin(ang_z)
       +z*sin(ang_y);
      nod_pt->x(1)=x*(cos(ang_x)*sin(ang_z) + cos(ang_z)*sin(ang_x)*sin(ang_y))
       +y*(cos(ang_x)*cos(ang_z) - sin(ang_x)*sin(ang_y)*sin(ang_z))
       -z*(cos(ang_y)*sin(ang_x));
      nod_pt->x(2)=x*(sin(ang_x)*sin(ang_z) - cos(ang_x)*cos(ang_z)*sin(ang_y))
       +y*(cos(ang_z)*sin(ang_x) + cos(ang_x)*sin(ang_y)*sin(ang_z))
       +z*(cos(ang_x)*cos(ang_y));
     } // for, loop over nodes
   } // SimpleCubicMesh
 };
}



//==start_of_problem_class============================================
/// Driven cavity problem in rectangular domain
//====================================================================
template<class ELEMENT>
class ParallelOutflowBoundaryProblem : public Problem
{

public:

 /// Constructor
 ParallelOutflowBoundaryProblem(const unsigned& n_element);

 /// Update before solve is empty
 void actions_before_newton_solve() 
 {
   // Initialise counters for each newton solve.
   Global_Parameters::Iterations.clear();
   Global_Parameters::Linear_solver_time.clear();
 } // actions_before_newton_solve

 void actions_after_newton_step()
 {
   Global_Parameters::Iterations.push_back(
       dynamic_cast<IterativeLinearSolver*>
       (this->linear_solver_pt())->iterations());

   Global_Parameters::Linear_solver_time.push_back(
       linear_solver_pt()->linear_solver_solution_time());
 }
 

 ///Fix pressure in element e at pressure dof pdof and set to pvalue
 void fix_pressure(const unsigned &e, const unsigned &pdof, 
                   const double &pvalue)
  {
   //Cast to full element type and fix the pressure at that element
   dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e))->
    fix_pressure(pdof,pvalue);
  }

 /// Doc the solution
 void doc_solution(DocInfo& doc_info);

 /// \short Create lagrange elements on boundary b of the Mesh pointed
 /// to by bulk_mesh_pt and add them to the Mesh object pointed to by 
 /// surface_mesh_pt
 void create_parall_outflow_lagrange_elements(const unsigned &b, 
                                              Mesh* const &bulk_mesh_pt,
                                              Mesh* const &surface_mesh_pt);
 /// \short Create lagrange elements on boundary b of the Mesh pointed
 /// to by bulk_mesh_pt and add them to the Mesh object pointed to by 
 /// surface_mesh_pt
 void create_impenetrable_lagrange_elements(const unsigned &b, 
                                              Mesh* const &bulk_mesh_pt,
                                              Mesh* const &surface_mesh_pt);
 /// Pointer to the "bulk" mesh
 Mesh*& bulk_mesh_pt() {return Bulk_mesh_pt;}

private:

 /// ID of imposed flow boundary
 unsigned Imposed_flow_boundary;

 /// ID of neumann boundary
 unsigned Neumann_boundary;

 /// Pointer to the "bulk" mesh
 Mesh* Bulk_mesh_pt;
 
 /// Pointer to the "surface" mesh for Lagrange multiplier elements
 Mesh* Surface_mesh_PO_pt;
 Mesh* Surface_mesh_IB_pt;

 /// Solver
 Preconditioner* Prec_pt;

 /// Solver
 IterativeLinearSolver* Solver_pt;

};



//==start_of_constructor==================================================
/// Constructor for DrivenCavity problem 
//========================================================================
template<class ELEMENT> 
ParallelOutflowBoundaryProblem<ELEMENT>::ParallelOutflowBoundaryProblem(
 const unsigned& n_el)
{ 

 // Setup mesh
 
 // # of elements in x-direction
 unsigned n_x=n_el;
 
 // # of elements in y-direction
 unsigned n_y=n_el;

 // # of elements in z-direction
 unsigned n_z=n_el;
 
 // Domain length in x-direction
 double l_x=1.0;
 
 // Domain length in y-direction
 double l_y=1.0;
 
 // Domain length in y-direction
 double l_z=1.0;
 
 // Build and assign mesh
 Bulk_mesh_pt = 
  new SlopingCubicMesh<ELEMENT >(n_x,n_y,n_z,l_x,l_y,l_z,
                                 Global_Parameters::Ang_x,
                                 Global_Parameters::Ang_y,
                                 Global_Parameters::Ang_z);

// Create "surface mesh" that will contain only the Lagrange multiplier 
 // elements.
 Surface_mesh_PO_pt = new Mesh;
 Surface_mesh_IB_pt = new Mesh;

 //Imposed_flow_boundary=5;
 //Neumann_boundary=0;
 
 // 0 0 0
 Imposed_flow_boundary=4;
 unsigned PO_b=3;
 unsigned IB_b=1;

 // remember to change this function as well.
// create_impenetrable_lagrange_elements(IB_b,
//                                         Bulk_mesh_pt,Surface_mesh_IB_pt);
 create_parall_outflow_lagrange_elements(PO_b,
                                         Bulk_mesh_pt,Surface_mesh_PO_pt);

 // Add the two sub meshes to the problem
 add_sub_mesh(Bulk_mesh_pt);
 add_sub_mesh(Surface_mesh_PO_pt);
 //add_sub_mesh(Surface_mesh_IB_pt);

 // Combine all submeshes into a single Mesh
 build_global_mesh();
 
 // Set the boundary conditions for this problem: All nodes are
 // free by default -- just pin the ones that have Dirichlet conditions
 // here. 
 unsigned num_bound = Bulk_mesh_pt->nboundary();
 for(unsigned ibound=0;ibound<num_bound;ibound++)
  {
   //if((ibound != PO_b)&&(ibound!=IB_b))
   if(ibound != PO_b)
    {
     unsigned num_nod= Bulk_mesh_pt->nboundary_node(ibound);
     for (unsigned inod=0;inod<num_nod;inod++)
      {
       // Loop over values (u, v and w velocities)
       for (unsigned i=0;i<3;i++)
        {
         Bulk_mesh_pt->boundary_node_pt(ibound,inod)->pin(i);
        }
      }
    }
  }
 
 // Impose inflow and negative inflow (outflow).
  unsigned num_nod= Bulk_mesh_pt->nboundary_node(Imposed_flow_boundary);
  for (unsigned inod=0;inod<num_nod;inod++)
   {
    Node* nod_pt=Bulk_mesh_pt->boundary_node_pt(Imposed_flow_boundary,inod);
    double x=nod_pt->x(0);
    double y=nod_pt->x(1);
    double z=nod_pt->x(2);
      
    // Reverse the tilting
    // We apply R_zyx
    double tilt_back_x = 0;
    double tilt_back_y = 0;
    double tilt_back_z = 0;

    tilt_back_x = x*(cos(-Global_Parameters::Ang_y)
                     *cos(-Global_Parameters::Ang_z))
                  +y*(cos(-Global_Parameters::Ang_z)
                      *sin(-Global_Parameters::Ang_x)
                    *sin(-Global_Parameters::Ang_y) 
                    - cos(-Global_Parameters::Ang_x)
                    *sin(-Global_Parameters::Ang_z))
                  +z*(sin(-Global_Parameters::Ang_x)
                      *sin(-Global_Parameters::Ang_z) 
                      + cos(-Global_Parameters::Ang_x)
                      *cos(-Global_Parameters::Ang_z)
                      *sin(-Global_Parameters::Ang_y));

    tilt_back_y = x*(cos(-Global_Parameters::Ang_y)
                            *sin(-Global_Parameters::Ang_z))
     + y*(cos(-Global_Parameters::Ang_x)
          *cos(-Global_Parameters::Ang_z) 
          + sin(-Global_Parameters::Ang_x)
          *sin(-Global_Parameters::Ang_y)
          *sin(-Global_Parameters::Ang_z))
     + z*(cos(-Global_Parameters::Ang_x)
          *sin(-Global_Parameters::Ang_y)
          *sin(-Global_Parameters::Ang_z) 
          - cos(-Global_Parameters::Ang_z)
          *sin(-Global_Parameters::Ang_x));
          
    tilt_back_z = -x*sin(-Global_Parameters::Ang_y)
     +y*cos(-Global_Parameters::Ang_y)
     *sin(-Global_Parameters::Ang_x)
     +z*cos(-Global_Parameters::Ang_x)
     *cos(-Global_Parameters::Ang_y);

    // The imposed velocity
    double ref_u_x = 0.0;
    double ref_u_y = 0.0;
    double ref_u_z = 0.0;
    
    //*
    // 0 0 0
    //if(tilt_back_y > 0.5)
    {
      ref_u_x=(tilt_back_y)*(1.0-tilt_back_y)
               *(tilt_back_z)*(1.0-tilt_back_z);
    }
   // else
   // {
   //   ref_u_x=-(tilt_back_y-0.0)*(0.5-tilt_back_y)
   //            *(tilt_back_z)*(1.0-tilt_back_z);
   // }

    // Now apply Rxyz to u, using rotation matrices.
    // We have velocity in the x direction only.
    // Thus the vector to rotate is [u,0,0] since the imposed
    // velocity in the y and z direction is 0.
    
    double imposed_u_x = 0.0;
    double imposed_u_y = 0.0;
    double imposed_u_z = 0.0;

     imposed_u_x = ref_u_x*cos(Global_Parameters::Ang_y)
                          *cos(Global_Parameters::Ang_z)
                   -ref_u_y*cos(Global_Parameters::Ang_y)
                           *sin(Global_Parameters::Ang_z)
                   +ref_u_z*sin(Global_Parameters::Ang_y);
     imposed_u_y = ref_u_x*(cos(Global_Parameters::Ang_x)
                          *sin(Global_Parameters::Ang_z) 
                          +cos(Global_Parameters::Ang_z)
                          *sin(Global_Parameters::Ang_x)
                          *sin(Global_Parameters::Ang_y))
                   +ref_u_y*(cos(Global_Parameters::Ang_x)
                           *cos(Global_Parameters::Ang_z) 
                           -sin(Global_Parameters::Ang_x)
                           *sin(Global_Parameters::Ang_y)
                           *sin(Global_Parameters::Ang_z))
                   -ref_u_z*(cos(Global_Parameters::Ang_y)
                           *sin(Global_Parameters::Ang_x));
     imposed_u_z = ref_u_x*(sin(Global_Parameters::Ang_x)
                          *sin(Global_Parameters::Ang_z) 
                          -cos(Global_Parameters::Ang_x)
                          *cos(Global_Parameters::Ang_z)
                          *sin(Global_Parameters::Ang_y))
                   +ref_u_y*(cos(Global_Parameters::Ang_z)
                           *sin(Global_Parameters::Ang_x) 
                           +cos(Global_Parameters::Ang_x)
                           *sin(Global_Parameters::Ang_y)
                           *sin(Global_Parameters::Ang_z))
                   +ref_u_z*(cos(Global_Parameters::Ang_x)
                           *cos(Global_Parameters::Ang_y));
  
    nod_pt->set_value(0,imposed_u_x);
    nod_pt->set_value(1,imposed_u_y);
    nod_pt->set_value(2,imposed_u_z);
   } // for
   
 // Complete the build of all elements so they are fully functional

 //Find number of elements in mesh
 unsigned n_element = Bulk_mesh_pt->nelement();

 // Loop over the elements to set up element-specific 
 // things that cannot be handled by constructor
 for(unsigned e=0;e<n_element;e++)
  {
   // Upcast from GeneralisedElement to the present element
   ELEMENT* el_pt = 
    dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e));
 // cout << "testdimmmm " << el_pt->dim() << endl;
 //pause("Dalle!"); 
   //Set the Reynolds number
   el_pt->re_pt() = &Global_Parameters::Re;
  }
 
 // Setup equation numbering scheme
 cout <<"Number of equations: " << assign_eqn_numbers() << std::endl; 
///* 
 // Build the preconditioner
 ParallelOutflowPreconditioner* prec_pt=new ParallelOutflowPreconditioner;
 Prec_pt = prec_pt;
 
Vector<Mesh*> meshes_pt;
meshes_pt.resize(2);
meshes_pt[0] = Bulk_mesh_pt;
meshes_pt[1] = Surface_mesh_PO_pt;
//meshes_pt[2] = Surface_mesh_IB_pt;
prec_pt->set_meshes(meshes_pt);

 //prec_pt->set_fluid_mesh(Bulk_mesh_pt);
 //prec_pt->set_lagrange_multiplier_mesh(Surface_mesh_pt);
 
 if(!Global_Parameters::Use_axnorm)
 {
   prec_pt->scaling_sigma() = Global_Parameters::Sigma;
 }
 
 if(Global_Parameters::Use_lsc)
 {
   prec_pt->use_lsc(); 
 }
 else
 {
   prec_pt->use_exact();
 }
 // Build solve and preconditioner
 Solver_pt = new GMRES<CRDoubleMatrix>;
 
 // We use RHS preconditioning. Note that by default, 
 // left hand preconditioning is used.
 dynamic_cast<GMRES<CRDoubleMatrix>*>(Solver_pt)->set_preconditioner_RHS();
 
 // Set solver and preconditioner
 Solver_pt->preconditioner_pt() = Prec_pt;
 linear_solver_pt() = Solver_pt;
//*/
}



//==start_of_doc_solution=================================================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void ParallelOutflowBoundaryProblem<ELEMENT>::doc_solution(DocInfo& doc_info)
{ 
 ofstream some_file;
 char filename[100];
 if(Global_Parameters::Doc_solution)
 {
   // Number of plot points
 unsigned npts;
 npts=5; 

 // Output solution 
 sprintf(filename,"%s/%s.dat",doc_info.directory().c_str(),
         Global_Parameters::Current_settings.c_str());
 some_file.open(filename);
 Bulk_mesh_pt->output(some_file,npts);
 some_file.close();
 }
 /*
 some_file.open("pressure_dofs.dat");
 unsigned nnod=Bulk_mesh_pt->nnode();
 for (unsigned j=0;j<nnod;j++)
  {
   if (Bulk_mesh_pt->node_pt(j)->nvalue()==3)
    {
     some_file << Bulk_mesh_pt->node_pt(j)->x(0) << " " 
               << Bulk_mesh_pt->node_pt(j)->x(1) << " " 
               << Bulk_mesh_pt->node_pt(j)->eqn_number(2) << " " 
               << std::endl;
    }
  }
 some_file.close();
 */
}


//============start_of_create_parall_outflow_lagrange_elements===========
/// Create ImposeParallelOutflowElement  on the b-th boundary of the
/// Mesh object pointed to by bulk_mesh_pt and add the elements to the 
/// Mesh object pointeed to by surface_mesh_pt.
//=======================================================================
template<class ELEMENT>
void ParallelOutflowBoundaryProblem<ELEMENT>::
create_impenetrable_lagrange_elements(const unsigned &b,
                                        Mesh* const &bulk_mesh_pt,
                                        Mesh* const &surface_mesh_pt)
{
 // How many bulk elements are adjacent to boundary b?
 unsigned n_element = bulk_mesh_pt->nboundary_element(b);

 // Loop over the bulk elements adjacent to boundary b?
 for(unsigned e=0;e<n_element;e++)
  {
   // Get pointer to the bulk element that is adjacent to boundary b
   ELEMENT* bulk_elem_pt = dynamic_cast<ELEMENT*>(
    bulk_mesh_pt->boundary_element_pt(b,e));
   
   //What is the index of the face of element e along boundary b
   int face_index = bulk_mesh_pt->face_index_at_boundary(b,e);
   
   // Build the corresponding my_lagrange_element
   ImposeImpenetrabilityElement<ELEMENT>* flux_element_pt = new 
    ImposeImpenetrabilityElement<ELEMENT>
     (bulk_elem_pt,face_index,Global_Parameters::Lagrange_multiplier_ib);
   
   //Add the prescribed-flux element to the surface mesh
   surface_mesh_pt->add_element_pt(flux_element_pt);

   // Loop over the nodes 
   unsigned nnod=flux_element_pt->nnode();
   for (unsigned j=0;j<nnod;j++)
    {
     Node* nod_pt = flux_element_pt->node_pt(j);
     
     // Is the node also on boundary 0, 1, 3 or 5?
     // 0 0 0 
     if ((nod_pt->is_on_boundary(0))||(nod_pt->is_on_boundary(2))
         || (nod_pt->is_on_boundary(4))||(nod_pt->is_on_boundary(5)))
     {
       // How many nodal values were used by the "bulk" element
       // that originally created this node?
       // Cast to a boundary node
       BoundaryNodeBase *bnod_pt =
        dynamic_cast<BoundaryNodeBase*>(nod_pt);
       
       // Get the index of the first Lagrange multiplier
       unsigned first_lmi=bnod_pt->
        index_of_first_value_assigned_by_face_element(
         Global_Parameters::Lagrange_multiplier_ib);
       
       // There is only one lagrange multiplier for the single 
       // outer unit normal vector.
       for (unsigned j=0;j<1;j++)
        {
         nod_pt->pin(first_lmi+j);
        }
      }
    }
  }
}



//============start_of_create_parall_outflow_lagrange_elements===========
/// Create ImposeParallelOutflowElement  on the b-th boundary of the
/// Mesh object pointed to by bulk_mesh_pt and add the elements to the 
/// Mesh object pointeed to by surface_mesh_pt.
//=======================================================================
template<class ELEMENT>
void ParallelOutflowBoundaryProblem<ELEMENT>::
create_parall_outflow_lagrange_elements(const unsigned &b,
                                        Mesh* const &bulk_mesh_pt,
                                        Mesh* const &surface_mesh_pt)
{
 // How many bulk elements are adjacent to boundary b?
 unsigned n_element = bulk_mesh_pt->nboundary_element(b);

 // Loop over the bulk elements adjacent to boundary b?
 for(unsigned e=0;e<n_element;e++)
  {
   // Get pointer to the bulk element that is adjacent to boundary b
   ELEMENT* bulk_elem_pt = dynamic_cast<ELEMENT*>(
    bulk_mesh_pt->boundary_element_pt(b,e));
   
   //What is the index of the face of element e along boundary b
   int face_index = bulk_mesh_pt->face_index_at_boundary(b,e);
   
   // Build the corresponding my_lagrange_element
   ImposeParallelOutflowElement<ELEMENT>* flux_element_pt = new 
    ImposeParallelOutflowElement<ELEMENT>
     (bulk_elem_pt,face_index,Global_Parameters::Lagrange_multiplier_po);

   flux_element_pt->boundary_id() = 1;
   //cout << "testing b_id = " << flux_element_pt->boundary_id() << endl;
   //Add the prescribed-flux element to the surface mesh
   surface_mesh_pt->add_element_pt(flux_element_pt);

   // Loop over the nodes 
   unsigned nnod=flux_element_pt->nnode();
   for (unsigned j=0;j<nnod;j++)
    {
     Node* nod_pt = flux_element_pt->node_pt(j);
     
     // Is the node also on boundary 0, 1, 3 or 5?
     // 0 0 0 
     if ((nod_pt->is_on_boundary(0)) || (nod_pt->is_on_boundary(2))
         || (nod_pt->is_on_boundary(4))||(nod_pt->is_on_boundary(5)))
     {
       // How many nodal values were used by the "bulk" element
       // that originally created this node?
       // Cast to a boundary node
       BoundaryNodeBase *bnod_pt =
        dynamic_cast<BoundaryNodeBase*>(nod_pt);
       
       // Get the index of the first Lagrange multiplier
       unsigned first_lmi=bnod_pt->
        index_of_first_value_assigned_by_face_element(
         Global_Parameters::Lagrange_multiplier_po);
       
       // There are two lagrange multipliers (two tangent vectors).
       for (unsigned j=0;j<2;j++)
        {
         nod_pt->pin(first_lmi+j);
        }
      }
    }
  }
}

int str2int(const string &str)
{
  stringstream ss(str);
  int n;
  ss >> n;
  return n;
}

unsigned str2unsigned(const string &str)
{
  stringstream ss(str);
  unsigned n;
  ss >> n;
  return n;
}

double str2double(const string &str)
{
  stringstream ss(str);
  double n;
  ss >> n;
  return n;
}

void set_dim_str(const string& Dim_str)
{
  Global_Parameters::Dim_str = Dim_str;
}

void set_prob_str(const string& Prob_str)
{
  Global_Parameters::Prob_str = Prob_str;
}

void set_vis_str(const string& do_stress_div)
{
  if (do_stress_div.compare("Str") == 0)
    {
     NavierStokesEquations<3>::Gamma[0]=1.0;
     NavierStokesEquations<3>::Gamma[1]=1.0;
     NavierStokesEquations<3>::Gamma[2]=1.0;
     Global_Parameters::Vis_str = "Str";
    }
   else
    {
     NavierStokesEquations<3>::Gamma[0]=0.0;
     NavierStokesEquations<3>::Gamma[1]=0.0;
     NavierStokesEquations<3>::Gamma[2]=0.0;
     Global_Parameters::Vis_str = "Sim";
    } // else - setting viscuous term.
}
  
void set_ang_str(const double& Ang_degree, const string& prefix_str,
                 string& Ang_str, double& Ang_radian)
{
  // Setting the string.
  stringstream Ang_str_stream;
  Ang_str_stream << prefix_str << Ang_degree;
  Ang_str = Ang_str_stream.str();

  // Setting the angle in radians.
  Ang_radian = Ang_degree * (Global_Parameters::Pi / 180.0);
}

void set_rey_str(const double& Re)
{
  Global_Parameters::Re = Re;
  stringstream Rey_str_stream;
  Rey_str_stream << "Rey" << Re;
  Global_Parameters::Rey_str = Rey_str_stream.str();
}

void set_noel_str(const unsigned& Noel)
{
  stringstream Noel_stream;
  Noel_stream << "Noel" << Noel;
  Global_Parameters::Noel_str = Noel_stream.str();
}

void set_current_settings()
{
 stringstream Current_settings_stream;

 Current_settings_stream << Global_Parameters::Prec_str
                         << Global_Parameters::Dim_str
                         << Global_Parameters::Prob_str
                         << Global_Parameters::Vis_str
                         << Global_Parameters::Ang_x_str
                         << Global_Parameters::Ang_y_str
                         << Global_Parameters::Ang_z_str
                         << Global_Parameters::Rey_str
                         << Global_Parameters::Noel_str;

 Global_Parameters::Current_settings = Current_settings_stream.str();
}

//===start_of_main======================================================
/// Driver code 
//======================================================================
int main(int argc, char* argv[]) 
{

#ifdef OOMPH_HAS_MPI
 MPI_Helpers::init(argc,argv);
#endif

 // Set up doc info
 DocInfo doc_info;
 doc_info.number()=0;
 doc_info.set_directory("RESLT");
 
 //Doc number of gmres iterations
 ofstream out_file;


 Global_Parameters::Dump_matrices = true;
 Global_Parameters::Doc_solution = true;

 string prec_str=argv[8];
 if(prec_str.compare("Exact")==0)
 {
   Global_Parameters::Use_lsc=false;
   Global_Parameters::Prec_str="Exact";
 }
 else
 {
   Global_Parameters::Use_lsc=true;
   Global_Parameters::Prec_str="Lsc";
 }
 string sigma_str = argv[7];
 if(sigma_str.compare("axnorm") == 0)
 {
   Global_Parameters::Use_axnorm = true;
 }
 else
 {
   Global_Parameters::Use_axnorm = false;
   Global_Parameters::Sigma = str2double(argv[7]);
 }

 // These naming parameters are hard coded.
 // Change these per file. 
 set_dim_str("3DIb");
 set_prob_str("PressureOut");

 set_vis_str(argv[1]);

 set_ang_str(str2double(argv[2]),"Angx",
             Global_Parameters::Ang_x_str,
             Global_Parameters::Ang_x);

 set_ang_str(str2double(argv[3]),"y",
             Global_Parameters::Ang_y_str,
             Global_Parameters::Ang_y);

 set_ang_str(str2double(argv[4]),"z",
             Global_Parameters::Ang_z_str,
             Global_Parameters::Ang_z);

 set_rey_str(str2double(argv[5]));
 
 unsigned n_el_1d = str2unsigned(argv[6]);
 set_noel_str(n_el_1d);

 set_current_settings();
 
 std::cout << Global_Parameters::Current_settings << std::endl;

 // The rest of the parameters are as follows: 
 // 1 - str, dimension of the problem, "3D" or "2D"
 // 2 - str, "PO" or "IB"
 // 3 - 
 // Set the tilting angle about the x/y/z axis
 //Global_Parameters::Ang_x = Global_Parameters::Pi/6;
 //Global_Parameters::Ang_y = Global_Parameters::Pi/6;
 //Global_Parameters::Ang_z = Global_Parameters::Pi/6;

 // Set the number of elements in one dimension
 //unsigned n_el_1d = 4;

         time_t rawtime;
         time ( &rawtime );
         
         cout << "RAYDOING: " 
              << Global_Parameters::Current_settings
              << " on " << ctime (&rawtime) << endl;

 ParallelOutflowBoundaryProblem<
   BulkFluidSubjectToLagrangeMultiplierElement<
     QTaylorHoodElement<3> > > problem(n_el_1d);

 // Solve the problem 
 problem.newton_solve();
           
 // Doc solution
 problem.doc_solution(doc_info);
 doc_info.number()++;

         // Output number of iterations
         unsigned iter = Global_Parameters::Iterations.size();
         double total_its = 0;
         cout << "RAYITS: ";
         for (unsigned j = 0; j < iter; j++)
         {
           total_its += Global_Parameters::Iterations[j];
           cout << Global_Parameters::Iterations[j] << " ";
         }
         double average_its = total_its/iter;
          
          // Print to one decimal place if the average its is not an exact
          // integer. Otherwise we print normally.
         ((int(average_its*10))%10)?
         cout << "\t"<< fixed << setprecision(1)
              << average_its << "(" << iter << ")" << endl:
         cout << "\t"<< average_its << "(" << iter << ")" << endl;


         // Output linear solver time 
         double total_time = 0;
         cout << "RAYTIME: " << setprecision(15);
         for (unsigned j = 0; j < iter; j++)
         {
           total_time += Global_Parameters::Linear_solver_time[j];
           cout << Global_Parameters::Linear_solver_time[j] << " ";
         }
         double average_time = total_time/iter;
          
          // Print to one decimal place if the average its is not an exact
          // integer. Otherwise we print normally.
         cout << "\t" << average_time << endl;
         
         time ( &rawtime );
         cout << "RAYDONE: " 
              << Global_Parameters::Current_settings
              << " on " << ctime (&rawtime) << endl;

#ifdef OOMPH_HAS_MPI
 MPI_Helpers::finalize();
#endif
}
