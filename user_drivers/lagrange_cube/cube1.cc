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
#include "ray.h"

// The 3D mesh
#include "meshes/simple_cubic_mesh.h"
#include "meshes/simple_cubic_tet_mesh.h"

using namespace std;
using namespace oomph;


struct CubeLagrangeVariables
{

  /// Lagrange Multiplier ID 
  unsigned Lagrange_multiplier_po;
  unsigned Lagrange_multiplier_ib;


  unsigned W_solver; // CL
  unsigned NS_solver; // CL
  unsigned F_solver; // CL
  unsigned P_solver; // CL
  unsigned Vis; // CL
  double Ang_x; // CL
  double Ang_y; // CL
  double Ang_z; // CL
  mutable double Rey; // CL
  unsigned Noel; // CL
  double Scaling_sigma; // CL

  std::string Prob_str; // To set, no CL
  std::string W_str; // To set from CL
  std::string NS_str; // To set from CL
  std::string F_str; // To set from CL
  std::string P_str; // To set from CL
  std::string Vis_str; // To set from CL
  std::string Ang_x_str; // To set from CL
  std::string Ang_y_str; // To set from CL
  std::string Ang_z_str; // To set from CL
  std::string Rey_str; // To set from CL
  std::string Noel_str; // To set from CL
  std::string Sigma_str; // Set from CL
  std::string W_approx_str;

  bool Use_axnorm; // To set from CL
  bool Use_diagonal_w_block;

  // Setting the defaults:
  CubeLagrangeVariables()
  {
    /// Lagrange Multiplier ID 
    Lagrange_multiplier_po = 42;
    Lagrange_multiplier_ib = 43;

    W_solver = 0; // 0 - Exact, no other solver coded in driver.
    NS_solver = 0; // 0 - Exact, 1 - LSC
    F_solver = 0; // 0 - Exact, 1 - multigrid.
    P_solver = 0; // 0 - Exact, 1 - multigrid.
    Vis = 0; // 0 - Sim, 1 - Str.
    Ang_x = 30.0;
    Ang_y = 30.0;
    Ang_z = 30.0;
    Rey = 100.0;
    Noel = 4;
    Scaling_sigma = 0; // This will be changed.

    Use_diagonal_w_block = true; // This will be set from the commandline.
    
    Use_axnorm = true; // Will be reset if Scaling sigma is provided.

    Prob_str = "SETPROBSTR"; // To set, no CL
    W_str = "We"; // To set from CL
    NS_str = "Nl"; // To set from CL
    F_str = "Fe"; // To set from CL
    P_str = "Pe"; // To set from CL
    Vis_str = "Sim"; // To set from CL
    Ang_x_str = "x30"; // To set from CL
    Ang_y_str = "y30"; // To set from CL
    Ang_z_str = "z30"; // To set from CL
    Rey_str = "r100"; // To set from CL
    Noel_str = "n4"; // To set from CL
    Sigma_str = ""; // Set from CL
    W_approx_str = "dw";
  }
};





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
 
 unsigned NS_solver = 1;
 unsigned F_solver = 0;
 unsigned P_solver = 0; 
 unsigned Use_diagonal_w_block = true;

} // Global_Parameters


namespace oomph {
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
 ParallelOutflowBoundaryProblem(CubeLagrangeVariables&, DocInfo&);

 /// Update before solve is empty
 void actions_before_newton_solve() 
 {
   // Initialise counters for each newton solve.
   Doc_info_pt->setup_new_time_step();
 } // actions_before_newton_solve

 void actions_after_newton_step()
 {
   
   unsigned iters = 0;
   double solver_time = 0.0;

#ifdef PARANOID
   IterativeLinearSolver* iterative_solver_pt
     = dynamic_cast<IterativeLinearSolver*>
       (this->linear_solver_pt());
   if(iterative_solver_pt == 0)
   {
     std::ostringstream error_message;
     error_message << "Cannot cast the solver pointer." << std::endl;

     throw OomphLibError(error_message.str(),
                         "TiltedCavityProblem",
                         OOMPH_EXCEPTION_LOCATION);
   }
   else
   {
     iters = iterative_solver_pt->iterations();
   }
#else
   iters = static_cast<IterativeLinearSolver*>
             (this->linear_solver_pt())->iterations();
#endif

   solver_time = linear_solver_pt()->linear_solver_solution_time();

   Doc_info_pt->add_iteration_and_time(iters,solver_time);

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

 DocInfo* Doc_info_pt;

 CubeLagrangeVariables* Myvar_pt;

};



//==start_of_constructor==================================================
/// Constructor for DrivenCavity problem 
//========================================================================
template<class ELEMENT> 
ParallelOutflowBoundaryProblem<ELEMENT>::ParallelOutflowBoundaryProblem(
 CubeLagrangeVariables &myvar, DocInfo &doc_info)
{ 

  Doc_info_pt = &doc_info;
  Myvar_pt = &myvar;


 // Setup mesh
 
 // # of elements in x-direction
 unsigned n_x=myvar.Noel;
 
 // # of elements in y-direction
 unsigned n_y=myvar.Noel;

 // # of elements in z-direction
 unsigned n_z=myvar.Noel;
 
 // Domain length in x-direction
 double l_x=1.0;
 
 // Domain length in y-direction
 double l_y=1.0;
 
 // Domain length in y-direction
 double l_z=1.0;
 
 // Build and assign mesh
 Bulk_mesh_pt = 
  new SlopingCubicMesh<ELEMENT >(n_x,n_y,n_z,l_x,l_y,l_z,
                                 myvar.Ang_x,
                                 myvar.Ang_y,
                                 myvar.Ang_z);

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
 create_impenetrable_lagrange_elements(IB_b,
                                         Bulk_mesh_pt,Surface_mesh_IB_pt);
 create_parall_outflow_lagrange_elements(PO_b,
                                         Bulk_mesh_pt,Surface_mesh_PO_pt);

 // Add the two sub meshes to the problem. The order we add the meshes are
 // important for preconditioning. The bulk mesh must go first, then the 
 // surface meshes.
 add_sub_mesh(Bulk_mesh_pt);
 add_sub_mesh(Surface_mesh_PO_pt);
 add_sub_mesh(Surface_mesh_IB_pt);

 // Combine all submeshes into a single Mesh
 build_global_mesh();
 
 // Set the boundary conditions for this problem: All nodes are
 // free by default -- just pin the ones that have Dirichlet conditions
 // here. 
 unsigned num_bound = Bulk_mesh_pt->nboundary();
 for(unsigned ibound=0;ibound<num_bound;ibound++)
  {
   if((ibound != PO_b)&&(ibound!=IB_b))
   //if(ibound != PO_b)
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

    tilt_back_x = x*(cos(-myvar.Ang_y)
                     *cos(-myvar.Ang_z))
                  +y*(cos(-myvar.Ang_z)
                      *sin(-myvar.Ang_x)
                    *sin(-myvar.Ang_y) 
                    - cos(-myvar.Ang_x)
                    *sin(-myvar.Ang_z))
                  +z*(sin(-myvar.Ang_x)
                      *sin(-myvar.Ang_z) 
                      + cos(-myvar.Ang_x)
                      *cos(-myvar.Ang_z)
                      *sin(-myvar.Ang_y));

    tilt_back_y = x*(cos(-myvar.Ang_y)
                            *sin(-myvar.Ang_z))
     + y*(cos(-myvar.Ang_x)
          *cos(-myvar.Ang_z) 
          + sin(-myvar.Ang_x)
          *sin(-myvar.Ang_y)
          *sin(-myvar.Ang_z))
     + z*(cos(-myvar.Ang_x)
          *sin(-myvar.Ang_y)
          *sin(-myvar.Ang_z) 
          - cos(-myvar.Ang_z)
          *sin(-myvar.Ang_x));
          
    tilt_back_z = -x*sin(-myvar.Ang_y)
     +y*cos(-myvar.Ang_y)
     *sin(-myvar.Ang_x)
     +z*cos(-myvar.Ang_x)
     *cos(-myvar.Ang_y);

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

     imposed_u_x = ref_u_x*cos(myvar.Ang_y)
                          *cos(myvar.Ang_z)
                   -ref_u_y*cos(myvar.Ang_y)
                           *sin(myvar.Ang_z)
                   +ref_u_z*sin(myvar.Ang_y);
     imposed_u_y = ref_u_x*(cos(myvar.Ang_x)
                          *sin(myvar.Ang_z) 
                          +cos(myvar.Ang_z)
                          *sin(myvar.Ang_x)
                          *sin(myvar.Ang_y))
                   +ref_u_y*(cos(myvar.Ang_x)
                           *cos(myvar.Ang_z) 
                           -sin(myvar.Ang_x)
                           *sin(myvar.Ang_y)
                           *sin(myvar.Ang_z))
                   -ref_u_z*(cos(myvar.Ang_y)
                           *sin(myvar.Ang_x));
     imposed_u_z = ref_u_x*(sin(myvar.Ang_x)
                          *sin(myvar.Ang_z) 
                          -cos(myvar.Ang_x)
                          *cos(myvar.Ang_z)
                          *sin(myvar.Ang_y))
                   +ref_u_y*(cos(myvar.Ang_z)
                           *sin(myvar.Ang_x) 
                           +cos(myvar.Ang_x)
                           *sin(myvar.Ang_y)
                           *sin(myvar.Ang_z))
                   +ref_u_z*(cos(myvar.Ang_x)
                           *cos(myvar.Ang_y));
  
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
   el_pt->re_pt() = &myvar.Rey;
  }
 
 // Setup equation numbering scheme
 cout <<"Number of equations: " << assign_eqn_numbers() << std::endl; 
///* 
 // Build the preconditioner
 LagrangeEnforcedflowPreconditioner* prec_pt=new LagrangeEnforcedflowPreconditioner;
 Prec_pt = prec_pt;
 
 Vector<Mesh*> meshes_pt;
 meshes_pt.resize(3);

 meshes_pt[0] = Bulk_mesh_pt;
 meshes_pt[1] = Surface_mesh_PO_pt;
 meshes_pt[2] = Surface_mesh_IB_pt;
 prec_pt->set_meshes(meshes_pt);


 if(!Global_Parameters::Use_axnorm)
 {
   prec_pt->scaling_sigma() = Global_Parameters::Sigma;
 }


 // W solver. Use SuperLU
 /* 
 if(myvar.W_solver == 0)
 {
 }
 else
 {
   std::cout << "Other W solvers not complemented yet. Using default SuperLU"
             << std::endl;
 }
*/

 // The preconditioner for the fluid block:
 ConstrainedNavierStokesSchurComplementPreconditioner* ns_preconditioner_pt =
 new ConstrainedNavierStokesSchurComplementPreconditioner;
 
 if(myvar.NS_solver == 0) // Exact solve.
 {
 }
 else if(myvar.NS_solver == 1) // LSC
 {
   prec_pt->set_navier_stokes_lsc_preconditioner(ns_preconditioner_pt);
   ns_preconditioner_pt->set_navier_stokes_mesh(Bulk_mesh_pt);

   // F block solve
   // myvar.F_solver == 0 is default, so do nothing.
   if(myvar.F_solver == 1)
   {
#ifdef OOMPH_HAS_HYPRE
     // LSC takes type "Preconditioner".
     Preconditioner* f_preconditioner_pt = new HyprePreconditioner;

     // Cast it so we can set the HYPRE properties.
     HyprePreconditioner* hypre_preconditioner_pt =
       static_cast<HyprePreconditioner*>(f_preconditioner_pt);

     // Set the HYPRE properties. See hypre_solver.cc for settings.
     Hypre_default_settings::
     set_defaults_for_2D_poisson_problem(hypre_preconditioner_pt);

     // Set the preconditioner in the LSC preconditioner.
     ns_preconditioner_pt->set_f_preconditioner(f_preconditioner_pt);
#endif
   }

   // P block solve
   //myvar.P_solver == 0 is default, so do nothing.
   if(myvar.P_solver == 1)
   {
#ifdef OOMPH_HAS_HYPRE
     Preconditioner* p_preconditioner_pt = new HyprePreconditioner;

     HyprePreconditioner* hypre_preconditioner_pt =
       static_cast<HyprePreconditioner*>(p_preconditioner_pt);

     Hypre_default_settings::
     set_defaults_for_2D_poisson_problem(hypre_preconditioner_pt);

     ns_preconditioner_pt->set_p_preconditioner(p_preconditioner_pt);
#endif
   }
 }
 else
 {
   pause("There is no solver for NS.");
 }


 // Set the doc info for book keeping purposes.
 //prec_pt->set_doc_info(&doc_info);

 if(myvar.Use_diagonal_w_block)
 {
   prec_pt->use_diagonal_w_block();
 }
 else
 {
   prec_pt->use_block_diagonal_w_block();
 }

 //if(doc_info.is_doc_prec_data_enabled())
 //{
 //  prec_pt->enable_doc_prec();
 //}



 // Build solve and preconditioner
 Solver_pt = new GMRES<CRDoubleMatrix>;
 
 // We use RHS preconditioning. Note that by default, 
 // left hand preconditioning is used.
 //dynamic_cast<GMRES<CRDoubleMatrix>*>(Solver_pt)->preconditioner_LHS()=false;
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
 //if(doc_info.is_doc_prec_data_enabled())
 {
   // Number of plot points
 unsigned npts;
 npts=5; 

 // Output solution 
 sprintf(filename,"%s/%s.dat",doc_info.directory().c_str(),
         doc_info.label().c_str());
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

   //flux_element_pt->boundary_id() = 1;
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

CubeLagrangeVariables myvar = CubeLagrangeVariables();

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
 set_dim_str("3DPo");
 set_prob_str("Cube");

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
   QTaylorHoodElement<3> > problem(myvar,doc_info);

 // Solve the problem 
 problem.newton_solve();
           
 // Doc solution
 problem.doc_solution(doc_info);
// doc_info.number()++;


 // Get the iteration count and time from the doc_info
 Vector<Vector<pair<unsigned, double> > > iters_times
   = doc_info.iterations_and_times();

 // Below outputs the iteration counts and time.
 // Output the number of iterations
 // Since this is a steady state problem, there is only
 // one "time step".
 //*

 // Loop over the time step:
 unsigned ntimestep = iters_times.size();
 for(unsigned intimestep = 0; intimestep < ntimestep; intimestep++)
 {
   // New timestep:
   std::cout << "RAYITS:\t" << intimestep << "\t";
   // Loop through the Newtom Steps
   unsigned nnewtonstep = iters_times[intimestep].size();
   unsigned sum_of_newtonstep_iters = 0;
   for(unsigned innewtonstep = 0; innewtonstep < nnewtonstep;
       innewtonstep++)
   {
      sum_of_newtonstep_iters += iters_times[intimestep][innewtonstep].first;
      std::cout << iters_times[intimestep][innewtonstep].first << " ";
   }
   double average_its = ((double)sum_of_newtonstep_iters)
                        / ((double)nnewtonstep);

   // Print to one decimal place if the average is not an exact
   // integer. Ohterwise we print normally.
   ((unsigned(average_its*10))%10)?
   std::cout << "\t"<< fixed << setprecision(1)
             << average_its << "(" << nnewtonstep << ")" << endl:
   std::cout << "\t"<< average_its << "(" << nnewtonstep << ")" << std::endl;

 }

 // Now doing the times
 for(unsigned intimestep = 0; intimestep < ntimestep; intimestep++)
 {
   // New timestep:
   std::cout << "RAYTIME:\t" << intimestep << "\t";
   // Loop through the Newtom Steps
   unsigned nnewtonstep = iters_times[intimestep].size();
   double sum_of_newtonstep_times = 0;
   for(unsigned innewtonstep = 0; innewtonstep < nnewtonstep;
       innewtonstep++)
   {
      sum_of_newtonstep_times += iters_times[intimestep][innewtonstep].second;
      std::cout << iters_times[intimestep][innewtonstep].second << " ";
   }
   double average_time = ((double)sum_of_newtonstep_times)
                        / ((double)nnewtonstep);

   // Print to one decimal place if the average is not an exact
   // integer. Ohterwise we print normally.
   std::cout << "\t"<< average_time << "(" << nnewtonstep << ")" << endl;
 }




//         // Output number of iterations
//         unsigned iter = Global_Parameters::Iterations.size();
//         double total_its = 0;
//         cout << "RAYITS: ";
//         for (unsigned j = 0; j < iter; j++)
//         {
//           total_its += Global_Parameters::Iterations[j];
//           cout << Global_Parameters::Iterations[j] << " ";
//         }
//         double average_its = total_its/iter;
//          
//          // Print to one decimal place if the average its is not an exact
//          // integer. Otherwise we print normally.
//         ((int(average_its*10))%10)?
//         cout << "\t"<< fixed << setprecision(1)
//              << average_its << "(" << iter << ")" << endl:
//         cout << "\t"<< average_its << "(" << iter << ")" << endl;
//
//
//         // Output linear solver time 
//         double total_time = 0;
//         cout << "RAYTIME: " << setprecision(15);
//         for (unsigned j = 0; j < iter; j++)
//         {
//           total_time += Global_Parameters::Linear_solver_time[j];
//           cout << Global_Parameters::Linear_solver_time[j] << " ";
//         }
//         double average_time = total_time/iter;
//          
//          // Print to one decimal place if the average its is not an exact
//          // integer. Otherwise we print normally.
//         cout << "\t" << average_time << endl;
         
         time ( &rawtime );
         cout << "RAYDONE: " 
              << Global_Parameters::Current_settings
              << " on " << ctime (&rawtime) << endl;

#ifdef OOMPH_HAS_MPI
 MPI_Helpers::finalize();
#endif
}
