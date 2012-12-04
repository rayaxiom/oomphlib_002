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
#ifndef OOMPH_CONSTRAINED_NAVIER_STOKES_PRECONDITIONERS_HEADER
#define OOMPH_CONSTRAINED_NAVIER_STOKES_PRECONDITIONERS_HEADER


// Config header generated by autoconfig
#ifdef HAVE_CONFIG_H
  #include <oomph-lib-config.h>
#endif


// oomphlib headers
#include "../generic/matrices.h"
#include "../generic/assembly_handler.h"
#include "../generic/problem.h"
#include "../generic/block_preconditioner.h"
#include "../generic/preconditioner.h"
#include "../generic/SuperLU_preconditioner.h"
#include "../generic/matrix_vector_product.h"
#include "navier_stokes_elements.h"
#include "refineable_navier_stokes_elements.h"
#include "navier_stokes_preconditioners.h"


namespace oomph
{

//===========================================================================
/// \short The least-squares commutator (LSC; formerly BFBT) Navier Stokes 
/// preconditioner. It uses blocks corresponding to the velocity 
/// and pressure unknowns, i.e. there are a total of 2x2 blocks, 
/// and all velocity components are treated as a single block of unknowns.
/// \n\n
/// Here are the details: An "ideal" Navier-Stokes preconditioner
/// would solve the system
/// \f[
/// \left( 
/// \begin{array}{cc}
/// {\bf F} & {\bf G} \\ {\bf D} & {\bf 0} 
/// \end{array} 
/// \right)
/// \left( 
/// \begin{array}{c}
/// {\bf z}_u \\ {\bf z}_p
/// \end{array} 
/// \right) =
/// \left( 
/// \begin{array}{c}
/// {\bf r}_u \\ {\bf r}_p
/// \end{array} 
/// \right)
/// \f]
/// where \f$ {\bf F}\f$ is the momentum block,  \f$ {\bf G} \f$ the
/// discrete gradient operator, and \f$ {\bf D}\f$ the discrete
/// divergence operator. (For unstabilised elements, we have 
/// \f$ {\bf D} = {\bf G}^T \f$ and in much of the literature
/// the divergence matrix is denoted by \f$ {\bf B} \f$ .)
/// The use of this preconditioner would ensure the convergence
/// of any iterative linear solver in a single iteration but its
/// application is, of course, exactly as expensive as a direct solve.
/// The LSC/BFBT preconditioner replaces the exact Jacobian by 
/// a block-triangular approximation
/// \f[
/// \left( 
/// \begin{array}{cc}
/// {\bf F} & {\bf G} \\ {\bf 0} & -{\bf M}_s 
/// \end{array} 
/// \right) 
/// \left( 
/// \begin{array}{c}
/// {\bf z}_u \\ {\bf z}_p
/// \end{array} 
/// \right) =
/// \left( 
/// \begin{array}{c}
/// {\bf r}_u \\ {\bf r}_p
/// \end{array} 
/// \right),
/// \f]
/// where \f${\bf M}_s\f$ is an approximation to the pressure 
/// Schur-complement \f$ {\bf S} = {\bf D} {\bf F}^{-1}{\bf G}. \f$
/// This system can be solved in two steps:
/// -# Solve the second row for \f$ {\bf z}_p\f$ via
///    \f[ 
///    {\bf z}_p = - {\bf M}_s^{-1} {\bf r}_p
///    \f]
/// -# Given \f$ {\bf z}_p \f$ , solve the first row for \f$ {\bf z}_u\f$ via
///    \f[ 
///    {\bf z}_u = {\bf F}^{-1} \big( {\bf r}_u - {\bf G} {\bf z}_p \big)
///    \f]
/// .
/// In the LSC/BFBT preconditioner, the action of the inverse pressure
/// Schur complement 
/// \f[
/// {\bf z}_p = - {\bf M}_s^{-1} {\bf r}_p
/// \f]
/// is approximated by
/// \f[
/// {\bf z}_p = - 
/// \big({\bf D} \widehat{\bf Q}^{-1}{\bf G} \big)^{-1}
/// \big({\bf D} \widehat{\bf Q}^{-1}{\bf F} \widehat{\bf Q}^{-1}{\bf G}\big) 
/// \big({\bf D} \widehat{\bf Q}^{-1}{\bf G} \big)^{-1}
/// {\bf r}_p,
/// \f]
/// where  \f$ \widehat{\bf Q} \f$ is the diagonal of the velocity
/// mass matrix. The evaluation of this expression involves
/// two linear solves involving the matrix
/// \f[
/// {\bf P} = \big({\bf D} \widehat{\bf Q}^{-1}{\bf G} \big)
/// \f]
/// which has the character of a matrix arising from the discretisation 
/// of a Poisson problem on the pressure space. We also have
/// to evaluate matrix-vector products with the matrix 
/// \f[
/// {\bf E}={\bf D}\widehat{\bf Q}^{-1}{\bf F}\widehat{\bf Q}^{-1}{\bf G}
/// \f]
/// Details of the theory can be found in "Finite Elements and 
/// Fast Iterative Solvers with Applications in Incompressible Fluid 
/// Dynamics" by Howard C. Elman, David J. Silvester, and Andrew J. Wathen,
/// published by Oxford University Press, 2006.
/// \n\n
/// In our implementation of the preconditioner, the linear systems
/// can either be solved "exactly", using SuperLU (in its incarnation
/// as an exact preconditioner; this is the default) or by any 
/// other Preconditioner (inexact solver) specified via the access functions
/// \code
/// ConstrainedNavierStokesSchurComplementPreconditioner::set_f_preconditioner(...)
/// \endcode
/// or 
/// \code
/// ConstrainedNavierStokesSchurComplementPreconditioner::set_p_preconditioner(...)
/// \endcode
//===========================================================================
 class ConstrainedNavierStokesSchurComplementPreconditioner :
 public BlockPreconditioner<CRDoubleMatrix>
 {
  
   public :
  
    /// Constructor - sets defaults for control flags
    ConstrainedNavierStokesSchurComplementPreconditioner() : 
     BlockPreconditioner<CRDoubleMatrix>()
    {
     // Use Robin BC elements for Fp preconditioner -- yes by default
     Use_robin_for_fp=true;

     // Flag to indicate that the preconditioner has been setup
     // previously -- if setup() is called again, data can
     // be wiped.
     Preconditioner_has_been_setup = false;

     // By default we use the LSC version
     Use_LSC=true;

     // By default we use SuperLU for both p and f blocks
     Using_default_p_preconditioner=true;
     Using_default_f_preconditioner=true;

     // Pin pressure dof in press adv diff problem for Fp precond
     Pin_first_pressure_dof_in_press_adv_diff=true;

     // resize the mesh pt
     // note: meaningless if subsidiary preconditioner
     this->set_nmesh(1);
     Navier_stokes_mesh_pt = 0;

     // Set default preconditioners (inexact solvers) -- they are 
     // members of this class!
     P_preconditioner_pt = 0;
     F_preconditioner_pt = 0;

     // set Doc_time to false
     Doc_time = false;

     // null the off diagonal Block matrix pt
     Bt_mat_vec_pt = 0;

     // null the F matrix vector product helper
     F_mat_vec_pt = 0;

     // null the QBt matrix vector product pt
     QBt_mat_vec_pt = 0;

     // null the E matrix vector product helper in Fp
     E_mat_vec_pt = 0;

     // RAYRAY initialise all my variables.
     Dim = 0;

     N_velocity_doftypes = 0;
    }

   /// Destructor
   ~ConstrainedNavierStokesSchurComplementPreconditioner()
    {
     clean_up_memory();
    }

   /// Broken copy constructor
   ConstrainedNavierStokesSchurComplementPreconditioner(
    const ConstrainedNavierStokesSchurComplementPreconditioner&)
    {
     BrokenCopy::broken_copy("ConstrainedNavierStokesSchurComplementPreconditioner");
    }
   
   /// Broken assignment operator
   void operator=(const ConstrainedNavierStokesSchurComplementPreconditioner&)
    {
     BrokenCopy::broken_assign("ConstrainedNavierStokesSchurComplementPreconditioner");
    }
   
//   // RAYRAY Set the blocks for the Navier Stokes preconditioner
//   void set_prec_blocks(Vector<CRDoubleMatrix*> &required_prec_blocks)
//   { 
//#ifdef PARANOID
//     if(required_prec_blocks.size() != 3)                                                   
//     {                                                                           
//       std::ostringstream error_message;                                         
//       error_message << "There must be three blocks for the\n"
//                     << "LSC preconditioner, for F, B and Bt" << std::endl;          
//       throw OomphLibError(error_message.str(),                                  
//                          "ConstrainedNavierStokesSchurComplementPreconditioner",                      
//                          OOMPH_EXCEPTION_LOCATION);                            
//     }
//     for (unsigned block_i = 0; block_i < 3; block_i++) 
//     {
//       if (required_prec_blocks[block_i] == 0) 
//       {
//       std::ostringstream error_message;                                         
//       error_message << "Block " << block_i << " is not set." << std::endl;
//       throw OomphLibError(error_message.str(),                                  
//                          "ConstrainedNavierStokesSchurComplementPreconditioner",                      
//                          OOMPH_EXCEPTION_LOCATION); 
//       }
//     }
//#endif
//     Prec_blocks = required_prec_blocks;
//
//     // set bool Prec_blocks_has_been_set = true - will I ever have to set it to
//     // false?
//   }
//
//   // RAYRAY 
//   void set_master_doftype_ordering(Vector<unsigned> &block_ordering)
//   {
//     unsigned ndof_types = this->ndof_types();
//     unsigned nblocks = block_ordering.size();
//     // check if this has the same ndof as the ndoftype of this preconditioner.
//#ifdef PARANOID
//     if(ndof_types > nblocks)
//     {
//       std::ostringstream error_message;                                         
//       error_message << "The number of blocks in the master preconditioner\n"
//                     << "must be equal to, or more than, the number of \n"
//                     << "blocks in the subsidiary preconditioner." 
//                     << std::endl;
//       throw OomphLibError(error_message.str(),
//                          "ConstrainedNavierStokesSchurComplementPreconditioner",
//                          OOMPH_EXCEPTION_LOCATION); 
//     }
//#endif
//     // Set the Master_doftype_order vector.
//     Master_doftype_order = block_ordering;
//     
//     // Resize it, to cut off the dof types we do not require.
//     Master_doftype_order.resize(ndof_types);
//   }
   
   /// Setup the preconditioner
   void setup(Problem* problem_pt, DoubleMatrixBase* matrix_pt);

   /// Apply preconditioner to Vector r
   void preconditioner_solve(const DoubleVector&r, DoubleVector &z);
   
   /// specify the mesh containing the mesh containing the 
   /// block-preconditionable Navier-Stokes elements. The dimension of the
   /// problem must also be specified.
   void set_navier_stokes_mesh(Mesh* mesh_pt)
    { 
     Navier_stokes_mesh_pt = mesh_pt;
    }

   /// Function to set a new pressure matrix preconditioner (inexact solver)
   void set_p_preconditioner(Preconditioner* new_p_preconditioner_pt)
   {
    // If the default preconditioner has been used
    // clean it up now...
    if (Using_default_p_preconditioner)
     {
      delete P_preconditioner_pt;
     }
    P_preconditioner_pt = new_p_preconditioner_pt;
    Using_default_p_preconditioner = false;
   }

   /// \short Function to (re-)set pressure matrix preconditioner  (inexact 
   /// solver) to SuperLU
   void set_p_superlu_preconditioner()
   {
    if (!Using_default_p_preconditioner)
     {
      delete P_preconditioner_pt;
      P_preconditioner_pt = new SuperLUPreconditioner;
      Using_default_p_preconditioner = true;
     }
   }

   /// Function to set a new momentum matrix preconditioner (inexact solver)
   void set_f_preconditioner(Preconditioner* new_f_preconditioner_pt)
   {
    // If the default preconditioner has been used
    // clean it up now...
    if (Using_default_f_preconditioner)
     {
      delete F_preconditioner_pt;
     }
    F_preconditioner_pt = new_f_preconditioner_pt;
    Using_default_f_preconditioner = false;
   }

   /// Use LSC version of the preconditioner
   void use_lsc(){Use_LSC=true;}
   
   /// Use Fp version of the preconditioner
   void use_fp(){Use_LSC=false;}

   ///\short Function to (re-)set momentum matrix preconditioner (inexact 
   /// solver) to SuperLU
   void set_f_superlu_preconditioner()
   {
    if (!Using_default_f_preconditioner)
     {
      delete F_preconditioner_pt;
      F_preconditioner_pt = new SuperLUPreconditioner;
      Using_default_f_preconditioner = true;
     }
   }

   ///Enable documentation of time
   void enable_doc_time() {Doc_time = true;}

   ///Disable documentation of time
   void disable_doc_time() {Doc_time = false;}


   /// \short Helper function to delete preconditioner data.
   void clean_up_memory();

   /// \short Use  Robin BC elements for the Fp preconditioner
   void enable_robin_for_fp() {Use_robin_for_fp = true;}

   /// \short Don't use Robin BC elements for the Fp preconditioenr
   void disable_robin_for_fp() {Use_robin_for_fp = false;}

   /// \short Set boolean indicating that we want to pin first pressure 
   /// dof in Navier Stokes mesh when
   /// assembling the pressure advection diffusion system for
   /// Fp preconditoner -- needed at zero Reynolds number
   /// for non-enclosed flows but seems harmless in any case
   void pin_first_pressure_dof_in_press_adv_diff()
   {Pin_first_pressure_dof_in_press_adv_diff=true;}

    /// \short Set boolean indicating that we do not want to pin first pressure 
   /// dof in Navier Stokes mesh when
   /// assembling the pressure advection diffusion system for
   /// Fp preconditoner -- needed at zero Reynolds number
   /// for non-enclosed flows but seems harmless in any case
   void unpin_first_pressure_dof_in_press_adv_diff()
   {Pin_first_pressure_dof_in_press_adv_diff=false;}

   /// \short Validate auxiliary pressure advection diffusion problem
   /// in 2D
   template<class ELEMENT>
   void validate(DocInfo& doc_info, Problem* orig_problem_pt)
   {
    FpPressureAdvectionDiffusionProblem<ELEMENT> 
     aux_problem(Navier_stokes_mesh_pt,orig_problem_pt);
    aux_problem.validate(doc_info); 
   }
   
   /// \short Pin all non-pressure dofs
   void pin_all_non_pressure_dofs()
   {
    // Backup pin status and then pin all non-pressure degrees of freedom
    unsigned nelem=Navier_stokes_mesh_pt->nelement();
    for (unsigned e=0;e<nelem;e++)
     {
      // Get pointer to the bulk element that is adjacent to boundary b
      TemplateFreeNavierStokesEquationsBase* el_pt = 
       dynamic_cast<TemplateFreeNavierStokesEquationsBase*>(
        Navier_stokes_mesh_pt->element_pt(e));
      
      // Pin 
      el_pt->pin_all_non_pressure_dofs(Eqn_number_backup);
     }
   }
   
   /// Get the pressure advection diffusion matrix
   void get_pressure_advection_diffusion_matrix(CRDoubleMatrix& fp_matrix)
   {
    // Pin all non-pressure dofs 
    pin_all_non_pressure_dofs(); 
    
    // Get the spatial dimension
    unsigned ndim=Navier_stokes_mesh_pt->finite_element_pt(0)->dim();


    // Backup assembly handler and set new one
    AssemblyHandler* backed_up_assembly_handler_pt=
     problem_pt()->assembly_handler_pt();
    problem_pt()->assembly_handler_pt()=
     new FpPreconditionerAssemblyHandler(ndim);

    // Backup submeshes (if any)
    unsigned n_sub_mesh=problem_pt()->nsub_mesh();
    Vector<Mesh*> backed_up_sub_mesh_pt(n_sub_mesh);
    for (unsigned i=0;i<n_sub_mesh;i++)
     {
      backed_up_sub_mesh_pt[i]=problem_pt()->mesh_pt(i);
     }
    // Flush sub-meshes but don't call rebuild_global_mesh()
    // so we can simply re-instate the sub-meshes below.
    problem_pt()->flush_sub_meshes();
    
    // Backup the problem's own mesh pointer 
    Mesh* backed_up_mesh_pt=problem_pt()->mesh_pt();
    problem_pt()->mesh_pt()=Navier_stokes_mesh_pt;

#ifdef OOMPH_HAS_MPI

    // Backup the start and end elements for the distributed
    // assembly process
    Vector<unsigned> backed_up_first_el_for_assembly;
    Vector<unsigned> backed_up_last_el_for_assembly;
    problem_pt()->get_first_and_last_element_for_assembly(
     backed_up_first_el_for_assembly,
     backed_up_last_el_for_assembly);

    // Now re-assign
    problem_pt()->set_default_first_and_last_element_for_assembly();

#endif

    // Identify pinned pressure dof
    int pinned_pressure_eqn=-2;
    if (Pin_first_pressure_dof_in_press_adv_diff)
     {
      // Loop over all Navier Stokes elements to find first non-pinned 
      // pressure dof
      unsigned nel=Navier_stokes_mesh_pt->nelement();
      for (unsigned e=0;e<nel;e++)
       {
        TemplateFreeNavierStokesEquationsBase* bulk_elem_pt = 
         dynamic_cast<TemplateFreeNavierStokesEquationsBase*>(
          Navier_stokes_mesh_pt->element_pt(e));
        int local_eqn=bulk_elem_pt->p_local_eqn(0);
        if (local_eqn>=0)
         {
          pinned_pressure_eqn=bulk_elem_pt->eqn_number(local_eqn);
          break;
         }
       }



#ifdef OOMPH_HAS_MPI
      if (problem_pt()->distributed())
       {
        int pinned_pressure_eqn_local=pinned_pressure_eqn;
        int pinned_pressure_eqn_global=pinned_pressure_eqn;
        MPI_Allreduce(&pinned_pressure_eqn_local,
                      &pinned_pressure_eqn_global,1,MPI_INT,MPI_MIN,
                      this->problem_pt()->communicator_pt()->mpi_comm());
        pinned_pressure_eqn=pinned_pressure_eqn_global;
       }    
      
#endif
      
     }
    
    
    // Loop over all Navier Stokes elements
    unsigned nel=Navier_stokes_mesh_pt->nelement();
    for (unsigned e=0;e<nel;e++)
     {
      TemplateFreeNavierStokesEquationsBase* bulk_elem_pt = 
       dynamic_cast<TemplateFreeNavierStokesEquationsBase*>(
        Navier_stokes_mesh_pt->element_pt(e));
      
      // Set pinned pressure equation
      bulk_elem_pt->pinned_fp_pressure_eqn()=pinned_pressure_eqn;      
     }
    
 
    // Attach Robin BC elements
    unsigned nbound=Navier_stokes_mesh_pt->nboundary(); 
    if (Use_robin_for_fp)
     {
      // Loop over all boundaries of Navier Stokes mesh
      for (unsigned b=0;b<nbound;b++)
       {
        // How many bulk elements are adjacent to boundary b?
        unsigned n_element = Navier_stokes_mesh_pt->nboundary_element(b);
        
        // Loop over the bulk elements adjacent to boundary b?
        for(unsigned e=0;e<n_element;e++)
         {
          TemplateFreeNavierStokesEquationsBase* bulk_elem_pt =
           dynamic_cast<TemplateFreeNavierStokesEquationsBase*>(
            Navier_stokes_mesh_pt->boundary_element_pt(b,e));
          
          //What is the index of the face of the bulk element e on bondary b
          int face_index=Navier_stokes_mesh_pt->face_index_at_boundary(b,e);
          
          // Build face element
          bulk_elem_pt->build_fp_press_adv_diff_robin_bc_element(face_index);
          
         } //end of loop over bulk elements adjacent to boundary b
       }
     }
        
    // Get "Jacobian" of the modified system
    DoubleVector dummy_residuals;
    problem_pt()->get_jacobian(dummy_residuals,fp_matrix);

    // Kill Robin BC elements
    if (Use_robin_for_fp)
     {
      // Loop over all boundaries of Navier Stokes mesh
      for (unsigned b=0;b<nbound;b++)
       {
        // How many bulk elements are adjacent to boundary b?
        unsigned n_element = Navier_stokes_mesh_pt->nboundary_element(b);
        
        // Loop over the bulk elements adjacent to boundary b?
        for(unsigned e=0;e<n_element;e++)
         {
          
          TemplateFreeNavierStokesEquationsBase* bulk_elem_pt = 
           dynamic_cast<TemplateFreeNavierStokesEquationsBase*>(
            Navier_stokes_mesh_pt->boundary_element_pt(b,e));
          
          // Kill associated Robin elements
          bulk_elem_pt->delete_pressure_advection_diffusion_robin_elements();

         } //end of loop over bulk elements adjacent to boundary b
       }
     }

    // Reset pin status
    reset_pin_status();

#ifdef OOMPH_HAS_MPI

    // Reset start and end elements for the distributed
    // assembly process
    problem_pt()->set_first_and_last_element_for_assembly(
     backed_up_first_el_for_assembly,
     backed_up_last_el_for_assembly);
    
#endif

    // Cleanup and reset assembly handler
    delete problem_pt()->assembly_handler_pt();
    problem_pt()->assembly_handler_pt()=backed_up_assembly_handler_pt;

    // Re-instate submeshes. (No need to call rebuild_global_mesh()
    // as it was never unbuilt).
    for (unsigned i=0;i<n_sub_mesh;i++)
     {
      problem_pt()->add_sub_mesh(backed_up_sub_mesh_pt[i]);
     }

    
    // Reset the problem's mesh pointer 
    problem_pt()->mesh_pt()= backed_up_mesh_pt;    
   }

   
   /// \short Reset pin status of all values
   void reset_pin_status()
   {

    // Reset pin status for nodes
    unsigned nnod=Navier_stokes_mesh_pt->nnode();
    for (unsigned j=0;j<nnod;j++)
     {
      Node* nod_pt=Navier_stokes_mesh_pt->node_pt(j);
      unsigned nval=nod_pt->nvalue();
      for (unsigned i=0;i<nval;i++)
       {
        nod_pt->eqn_number(i)=Eqn_number_backup[nod_pt][i];
       }

      // If it's a solid node deal with its positional data too
      SolidNode* solid_nod_pt=dynamic_cast<SolidNode*>(nod_pt);
      if (solid_nod_pt!=0)
       {
        Data* solid_posn_data_pt=solid_nod_pt->variable_position_pt();
        unsigned nval=solid_posn_data_pt->nvalue();
        for (unsigned i=0;i<nval;i++)
         {
          solid_posn_data_pt->eqn_number(i)=
           Eqn_number_backup[solid_posn_data_pt][i];
         }
       }
      
     }
    
    // ... and internal data
    unsigned nelem=Navier_stokes_mesh_pt->nelement();
    for (unsigned e=0;e<nelem;e++)
     {
      // Pointer to element
      GeneralisedElement* el_pt=Navier_stokes_mesh_pt->element_pt(e);
      
      // Pin/unpin internal data
      unsigned nint=el_pt->ninternal_data();
      for (unsigned j=0;j<nint;j++)
       {
        Data* data_pt=el_pt->internal_data_pt(j);
        unsigned nvalue=data_pt->nvalue();
        for (unsigned i=0;i<nvalue;i++)
         {
          data_pt->eqn_number(i)=Eqn_number_backup[data_pt][i];
         }
       }
     }
    
    // Free up storage
    Eqn_number_backup.clear();
   }
  
    private:

   // oomph-lib objects
   // -----------------

   // RAYRAY: Used to determine the number of velocity directions
   // when building the Master_doftype_order vector. Any left over will
   // be pressure dof types.
   unsigned Dim;
   
   // RAYRAY: used to know when the velocity dof types end and pressure begins.
   // since we may have more than one pressure dof type.
   unsigned N_velocity_doftypes;
   
   
   // Pointers to preconditioner (=inexact solver) objects
   // -----------------------------------------------------
   /// Pointer to the 'preconditioner' for the pressure matrix
   Preconditioner* P_preconditioner_pt;

   /// Pointer to the 'preconditioner' for the F matrix
   Preconditioner* F_preconditioner_pt;

   /// flag indicating whether the default F preconditioner is used
   bool Using_default_f_preconditioner;

   /// flag indicating whether the default P preconditioner is used
   bool Using_default_p_preconditioner;

   /// \short Control flag is true if the preconditioner has been setup
   /// (used so we can wipe the data when the preconditioner is
   /// called again)
   bool Preconditioner_has_been_setup;

   /// \short Helper function to assemble the diagonal of the pressure
   /// and velocity mass matrices from the elemental contributions defined in
   /// NavierStokesEquations<DIM>.
   /// If do_both=true, both are computed, otherwise only the velocity
   /// mass matrix (the LSC version of the preconditioner only needs
   /// that one)
   //void assemble_inv_press_and_veloc_mass_matrix_diagonal(
   // CRDoubleMatrix*& inv_p_mass_pt, 
   // CRDoubleMatrix*& inv_v_mass_pt, 
   // const bool& do_both);
   void assemble_inv_press_and_veloc_mass_matrix_diagonal(
     CRDoubleMatrix*& inv_p_mass_pt, 
     CRDoubleMatrix*& inv_v_mass_pt, 
     const bool& do_both,
     const unsigned& procnumber,
     Problem* problem_pt);

   
   /// \short Boolean indicating whether the momentum system preconditioner 
   /// is a block preconditioner
   bool F_preconditioner_is_block_preconditioner;

   /// Set Doc_time to true for outputting results of timings
   bool Doc_time;

   /// MatrixVectorProduct operator for Qv^{-1} Bt 
   MatrixVectorProduct* QBt_mat_vec_pt;

   /// MatrixVectorProduct operator for Bt
   MatrixVectorProduct* Bt_mat_vec_pt;

   /// MatrixVectorProduct operator for F
   MatrixVectorProduct* F_mat_vec_pt;

   /// MatrixVectorProduct operator for E = Fp Qp^{-1} (only for Fp variant)
   MatrixVectorProduct* E_mat_vec_pt;

   /// \short the pointer to the mesh of block preconditionable Navier
   /// Stokes elements.
   Mesh* Navier_stokes_mesh_pt;

   /// Boolean to indicate use of LSC (true) or Fp (false) variant
   bool Use_LSC;
   
   /// Use Robin BC elements for Fp preconditoner?
   bool Use_robin_for_fp;

   /// \short Map to store original eqn numbers of various Data values when 
   /// assembling pressure advection diffusion matrix
   std::map<Data*,std::vector<int> > Eqn_number_backup;

   /// \short Boolean indicating if we want to pin first pressure 
   /// dof in Navier Stokes mesh when
   /// assembling the pressure advection diffusion system for
   /// Fp preconditoner -- needed at zero Reynolds number
   /// for non-enclosed flows but seems harmless in any case
   bool Pin_first_pressure_dof_in_press_adv_diff;
   
   };




///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////


/*
//============================================================================
/// \short The exact Navier Stokes preconditioner. This extracts 2x2 blocks
/// (corresponding to the velocity and pressure unknowns) and uses these to
/// build a single preconditioner matrix for testing purposes.
/// Iterative solvers should converge in a single step if this is used.
/// If it doesn't something is wrong in the setup of the block matrices.
//=============================================================================
 template<typename MATRIX>
  class NavierStokesExactPreconditioner : public BlockPreconditioner<MATRIX>
  {

   public :
     
    /// Constructor - do nothing
    NavierStokesExactPreconditioner() : BlockPreconditioner<MATRIX>(){}
   
   
   /// Destructor - do nothing
   ~NavierStokesExactPreconditioner(){}

   
   /// Broken copy constructor
   NavierStokesExactPreconditioner(const NavierStokesExactPreconditioner&)
    {
     BrokenCopy::broken_copy("NavierStokesExactPreconditioner");
    }
   

   /// Broken assignment operator   
   void operator=(const NavierStokesExactPreconditioner&)
    {
     BrokenCopy::broken_assign("NavierStokesExactPreconditioner");
    }
   
   
   /// Setup the preconditioner
   void setup(Problem* problem_pt, DoubleMatrixBase* matrix_pt);
   
   /// Apply preconditioner to r
   void preconditioner_solve(const Vector<double>&r,
                             Vector<double> &z);

   protected :
    
    /// Preconditioner matrix
    MATRIX P_matrix;

  };

 */

}
#endif
