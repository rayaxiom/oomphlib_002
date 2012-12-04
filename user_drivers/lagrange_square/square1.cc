// Working tangential flow.
// ./square1 --w_solver 0 --ns_solver 0 --visc Sim --ang 30 --rey 100 --noel 8 --diagw --doc_soln
//LIC//====================================================================
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
//#include "raymon.h"
#include "ray.h"

// The 2D mesh
#include "meshes/simple_rectangular_quadmesh.h"
#include "meshes/rectangular_quadmesh.h"

//using namespace std;
using namespace oomph;

struct SquareLagrangeVariables{
  // 0 = Exact, 1 = Lsc Exact, 2 = Lsc multigrid.
  unsigned W_solver; // CL
  unsigned NS_solver; // CL
  unsigned F_solver; // CL
  unsigned P_solver; // CL
  unsigned Vis; // CL
  double Ang; // CL
  double Rey; // CL
  unsigned Noel; // CL
  double Scaling_sigma; // CL

  std::string Prob_str; // To set, no CL
  std::string W_str; // To set from CL
  std::string NS_str; // To set from CL
  std::string F_str; // To set from CL
  std::string P_str; // To set from CL
  std::string Vis_str; // CL
  std::string Ang_str; // To set from CL
  std::string Rey_str; // To set from CL
  std::string Noel_str; // To set from CL
  std::string Sigma_str; // Set from CL
  std::string W_approx_str;

  bool Use_axnorm; // To set from CL
  bool Use_diagonal_w_block;

  // Setting the defaults:
  SquareLagrangeVariables() :
    W_solver(0), // 0 = Exact, no other W solver coded.
    NS_solver(1), // 0 = Exact, 1 = LSC
    F_solver(0), // 0 = Exact, 1 = AMG
    P_solver(0), // 0 = Exact, 1 = AMG
    Vis(0), // 0 - Simple, 1 - Stress divergence.
    Ang(30.0), // Angle, in degrees
    Rey(100.0), // Reynolds number
    Noel(4), // Number of elements.
    Scaling_sigma(0), // Scaling Sigma
    Prob_str("SqTf"), // This is set inside the code, not from commandline.
    W_str("We"), // e - Exact, no other solver.
    NS_str("Nl"), // e - Exact, l - LSC
    F_str("Fe"), // e - Exact, a - AMG
    P_str("Pe"), // e - Exact, a - AMG
    Vis_str("Sim"), // Sim - Simple, Str = Stress Divergence
    Ang_str("A30"), // z - angle of rotation, about the z axis.
    Rey_str("R100"), // Reynolds number.
    Noel_str("N4"), // Number of elements.
    Sigma_str(""), // Sigma string.
    W_approx_str(""), // diagonal approximation
    Use_axnorm(true), // Use norm of velocity in the x direction for Sigma.
    Use_diagonal_w_block(true) // Use the diagonal approximation for W.
  {}

};

namespace oomph
{

//========================================================================
/// \short A Sloping Mesh  class.
///
/// derived from RectangularQuadMesh:
/// the same mesh rotated with an angle phi
//========================================================================
 template<class ELEMENT>
 class SlopingQuadMesh : public RectangularQuadMesh<ELEMENT>
 {
 public:

  /// Constructor.
  SlopingQuadMesh(const unsigned& nx, const unsigned& ny,
                  const double& lx,  const double& ly, const double& phi ) :
   RectangularQuadMesh<ELEMENT>(nx,ny,lx,ly)
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

      // Set new nodal coordinates
      nod_pt->x(0)=x*cos(phi)-y*sin(phi);
      nod_pt->x(1)=x*sin(phi)+y*cos(phi);
     }
   }
 };




} // end of namespace


//===start_of_problem_class=============================================
//======================================================================

template<class ELEMENT>
class TiltedCavityProblem : public Problem
{
public:

 /// \short Constructor: Pass number of elements in x and y directions and
 /// lengths
 TiltedCavityProblem(SquareLagrangeVariables&,
                     DocInfo&);

 /// Update before solve is empty
 void actions_before_newton_solve()
 {
   // Initialise counters for each newton solve.
   Doc_info_pt->setup_new_time_step();
 }

 /// \short Update after solve is empty
 void actions_after_newton_solve()
 {
 }

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
 /// Doc the solution
 void doc_solution(DocInfo&);

 /// \short Create lagrange elements on boundary b of the Mesh pointed
 /// to by bulk_mesh_pt and add them to the Mesh object pointed to by
 /// surface_mesh_pt
 void create_parall_outflow_lagrange_elements(const unsigned &b,
                                              Mesh* const &bulk_mesh_pt,
                                              Mesh* const &surface_mesh_pt);

 void create_impenetrable_lagrange_elements(const unsigned &b,
                                            Mesh* const &bulk_mesh_pt,
                                            Mesh* const &surface_mesh_pt);

private:

 /// Pointer to the "bulk" mesh
 SlopingQuadMesh<ELEMENT>* Bulk_mesh_pt;

 /// Pointer to the "surface" mesh
 Mesh* Surface_mesh_T_pt;
 Mesh* Surface_mesh_P_pt;

 // Preconditioner
 Preconditioner* Prec_pt;
 // Solver
 IterativeLinearSolver* Solver_pt;

 DocInfo* Doc_info_pt;
};



//===start_of_constructor=============================================
/// Problem constructor
//====================================================================
template<class ELEMENT> // rrrback - changed here.
TiltedCavityProblem<ELEMENT>::TiltedCavityProblem
(SquareLagrangeVariables &myvar, DocInfo &doc_info)
{

 Doc_info_pt = &doc_info;

 // Assign the boundaries:
 unsigned if_b=3;
 unsigned tf_b=1;
 //unsigned po_b=1;

 /// Setup the mesh
 // # of elements in x-direction
 unsigned nx=myvar.Noel;

 // # of elements in y-direction
 unsigned ny=myvar.Noel;

 // Domain length in x-direction
 double lx=1.0;

 // Domain length in y-direction
 double ly=1.0;

 Bulk_mesh_pt =
  new SlopingQuadMesh<ELEMENT>(nx,ny,lx,ly,myvar.Ang);

 // Create a "surface mesh" that will contain only
 // ImposeParallelOutflowElements in boundary 1
 // The constructor just creates the mesh without
 // giving it any elements, nodes, etc.
 //Surface_mesh_P_pt = new Mesh;
 Surface_mesh_T_pt = new Mesh;

 // Create ImposeParallelOutflowElement from all elements that are
 // adjacent to the Neumann boundary.
 //create_parall_outflow_lagrange_elements(po_b,
 //                                        Bulk_mesh_pt,Surface_mesh_P_pt);
 create_impenetrable_lagrange_elements(tf_b,
                                         Bulk_mesh_pt,Surface_mesh_T_pt);

 // Add the two sub meshes to the problem
 add_sub_mesh(Bulk_mesh_pt);
 //add_sub_mesh(Surface_mesh_P_pt);
 add_sub_mesh(Surface_mesh_T_pt);
 // Combine all submeshes into a single Mesh
 build_global_mesh();

 unsigned num_bound=Bulk_mesh_pt->nboundary();

 // Set the boundary conditions for this problem: All nodes are
 // free by default -- just pin the ones that have Dirichlet conditions
 // here.
 for(unsigned ibound=0;ibound<num_bound;ibound++)
  {
  //if((ibound != po_b)&&(ibound != tf_b))
  if(ibound != tf_b)
  {
   unsigned num_nod=Bulk_mesh_pt->nboundary_node(ibound);

   for (unsigned inod=0;inod<num_nod;inod++)
    {
     // Get node
     Node* nod_pt=Bulk_mesh_pt->boundary_node_pt(ibound,inod);

     nod_pt->pin(0);
     nod_pt->pin(1);

    }
  }
  }

 // Inflow on upper half of the Imposed_flow_boundary.
 // Outflow on lower half of the Imposed_flow_boundary.
 unsigned num_nod= mesh_pt()->nboundary_node(if_b);
 for(unsigned inod=0;inod<num_nod;inod++)
 {
   Node* nod_pt=mesh_pt()->boundary_node_pt(if_b,inod);
   double x=nod_pt->x(0);
   double y=nod_pt->x(1);

   // Tilt it back
   double ytiltedback = x*sin(-myvar.Ang)
                        +y*cos(-myvar.Ang);
   double u=0.0;
   //u=(ytiltedback-0.0)*(1-ytiltedback);


   //*
   if(ytiltedback > 0.5)
    {
     // Impose inflow velocity
     u=(ytiltedback-0.5)*(1-ytiltedback);
    }
   else
    {
     // Impose outflow velocity
     u=(ytiltedback-0.5)*ytiltedback;
    }
    // */

   // Now apply the rotation to u, using rotation matrices.
   // with x = u and y = 0, i.e. R*[u;0] since we have the
   // velocity in the x direction only. There is no velocity
   // in the y direction.
   double ux=u*cos(myvar.Ang);
   double uy=u*sin(myvar.Ang);

   nod_pt->set_value(0,ux);
   nod_pt->set_value(1,uy);
 }

 //Complete the problem setup to make the elements fully functional

 //Loop over the elements
 unsigned n_el = Bulk_mesh_pt->nelement();
 for(unsigned e=0;e<n_el;e++)
  {
   //Cast to a fluid element
   ELEMENT *el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e));

   //Set the Reynolds number, etc
   el_pt->re_pt() = &myvar.Rey;

  } // for(unsigned e=0;e<n_el;e++)

 //Assgn equation numbers
 cout << "\n equation numbers : "<< assign_eqn_numbers() << std::endl;

 ////// Build the preconditioner
 LagrangeEnforcedflowPreconditioner* prec_pt=new LagrangeEnforcedflowPreconditioner;
 Prec_pt = prec_pt;

 Vector<Mesh*> meshes_pt;
 meshes_pt.resize(2);
 meshes_pt[0] = Bulk_mesh_pt;
 //meshes_pt[1] = Surface_mesh_P_pt;
 meshes_pt[1] = Surface_mesh_T_pt;
 prec_pt->set_meshes(meshes_pt);
 

 if(!myvar.Use_axnorm)
 {
   prec_pt->scaling_sigma() = myvar.Scaling_sigma;
 }

 //////////////////////////////////////////////////////////////////////////////
 // Setting up the solver an preconditioners.

 // W solver. Use SuperLU
 if(myvar.W_solver == 0)
 {
 }
 else
 {
   std::cout << "Other W solvers not complemented yet. Using default SuperLU"
             << std::endl;
 }

 ConstrainedNavierStokesSchurComplementPreconditioner* ns_preconditioner_pt =
 new ConstrainedNavierStokesSchurComplementPreconditioner;

 // The preconditioner for the fluid block:

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
 prec_pt->set_doc_info(&doc_info);

 if(myvar.Use_diagonal_w_block)
 {
   prec_pt->use_diagonal_w_block();
 }
 else
 {
   prec_pt->use_block_diagonal_w_block();
 }

 if(doc_info.is_doc_prec_data_enabled())
 {
   prec_pt->enable_doc_prec();
 }

 // Build solve and preconditioner
 Solver_pt = new GMRES<CRDoubleMatrix>;

 // We use RHS preconditioning. Note that by default,
 // left hand preconditioning is used.
 static_cast<GMRES<CRDoubleMatrix>*>(Solver_pt)->set_preconditioner_RHS();

 // Set solver and preconditioner
 Solver_pt->preconditioner_pt() = Prec_pt;
 linear_solver_pt() = Solver_pt;

}



//==start_of_doc_solution=================================================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void TiltedCavityProblem<ELEMENT>::doc_solution(DocInfo& doc_info)
{
  if(doc_info.is_doc_solution_enabled())
  {
    ofstream some_file;
    char filename[100];

    // Number of plot points
    unsigned npts=5;

    // Output solution
    sprintf(filename,"%s/%s.dat",doc_info.directory().c_str(),
            doc_info.label().c_str());
    some_file.open(filename);
    Bulk_mesh_pt->output(some_file,npts);
    some_file.close();
  } // if
}



//============start_of_create_parall_outflow_lagrange_elements===========
/// Create ImposeParallelOutflowElement on the b-th boundary of the
/// Mesh object pointed to by bulk_mesh_pt and add the elements to the
/// Mesh object pointeed to by surface_mesh_pt.
//=======================================================================
template<class ELEMENT>
void TiltedCavityProblem<ELEMENT>::
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

   // Build the corresponding impose_impenetrability_element
   ImposeParallelOutflowElement<ELEMENT>* flux_element_pt = new
    ImposeParallelOutflowElement<ELEMENT>(bulk_elem_pt,
                                          face_index);


   // Add the prescribed-flux element to the surface mesh
   surface_mesh_pt->add_element_pt(flux_element_pt);

   // Loop over the nodes
   unsigned nnod=flux_element_pt->nnode();
   for (unsigned j=0;j<nnod;j++)
    {
     Node* nod_pt = flux_element_pt->node_pt(j);

     // Is the node also on boundary 0 or 2?
     if ((nod_pt->is_on_boundary(0))||(nod_pt->is_on_boundary(2)))
      {
       // How many nodal values were used by the "bulk" element
       // that originally created this node?
       unsigned n_bulk_value=flux_element_pt->nbulk_value(j);

       // The remaining ones are Lagrange multipliers and we pin them.
       unsigned nval=nod_pt->nvalue();
       for (unsigned j=n_bulk_value;j<nval;j++)
        {
         nod_pt->pin(j);
        }
      }
    }
  }
}

//============start_of_create_parall_outflow_lagrange_elements===========
/// Create ImposeParallelOutflowElement on the b-th boundary of the
/// Mesh object pointed to by bulk_mesh_pt and add the elements to the
/// Mesh object pointeed to by surface_mesh_pt.
//=======================================================================
template<class ELEMENT>
void TiltedCavityProblem<ELEMENT>::
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

   // Build the corresponding impose_impenetrability_element
   ImposeImpenetrabilityElement<ELEMENT>* flux_element_pt = new
    ImposeImpenetrabilityElement<ELEMENT>(bulk_elem_pt,
                                          face_index);

   // Add the prescribed-flux element to the surface mesh
   surface_mesh_pt->add_element_pt(flux_element_pt);

   // Loop over the nodes
   unsigned nnod=flux_element_pt->nnode();
   for (unsigned j=0;j<nnod;j++)
    {
     Node* nod_pt = flux_element_pt->node_pt(j);

     // Is the node also on boundary 0 or 2?
     if ((nod_pt->is_on_boundary(0))||(nod_pt->is_on_boundary(2)))
      {
       // How many nodal values were used by the "bulk" element
       // that originally created this node?
       unsigned n_bulk_value=flux_element_pt->nbulk_value(j);

       // The remaining ones are Lagrange multipliers and we pin them.
       unsigned nval=nod_pt->nvalue();
       for (unsigned j=n_bulk_value;j<nval;j++)
        {
         nod_pt->pin(j);
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

//===start_of_main======================================================
/// Driver code
//======================================================================
int main(int argc, char* argv[])
{
 #ifdef OOMPH_HAS_MPI
 // Initialise MPI
 MPI_Helpers::init(argc,argv);
 #endif

 /*
 SquareLagrangeVariables problemvars;
 problemvars.Dim = 2;
 std::cout << "probvar: " << problemvars.Dim << endl;

 // Store commandline arguments
 CommandLineArgs::setup(argc,argv);
 CommandLineArgs::specify_command_line_flag("--dim", &problemvars.Dim);
cout << "new dim: " << problemvars.Dim << endl;
CommandLineArgs::parse_and_assign();
CommandLineArgs::doc_specified_flags();
cout << "new dim: " << problemvars.Dim << endl;

pause("closer");

 */
//cout << "MathematicalConstants::Pi " << MathematicalConstants::Pi << endl;
//pause("damn");

 SquareLagrangeVariables myvar = SquareLagrangeVariables();

 // Set up doc info
 DocInfo doc_info;
 //doc_info.number()=0;
 doc_info.set_directory("RESLT");

 //Doc number of gmres iterations
 //ofstream out_file;

 // Store commandline arguments
 CommandLineArgs::setup(argc,argv);


 // Flag to output the solution.
 CommandLineArgs::specify_command_line_flag("--doc_soln");
 // Flag to output the preconditioner, used for debugging.
 CommandLineArgs::specify_command_line_flag("--doc_prec");

 CommandLineArgs::specify_command_line_flag("--w_solver", &myvar.W_solver);
 CommandLineArgs::specify_command_line_flag("--ns_solver", &myvar.NS_solver);
 CommandLineArgs::specify_command_line_flag("--p_solver", &myvar.P_solver);
 CommandLineArgs::specify_command_line_flag("--f_solver", &myvar.F_solver);
 CommandLineArgs::specify_command_line_flag("--visc", &myvar.Vis);
 CommandLineArgs::specify_command_line_flag("--ang", &myvar.Ang);
 CommandLineArgs::specify_command_line_flag("--rey", &myvar.Rey);
 CommandLineArgs::specify_command_line_flag("--noel", &myvar.Noel);
 CommandLineArgs::specify_command_line_flag("--sigma",
                                            &myvar.Scaling_sigma);
 CommandLineArgs::specify_command_line_flag("--bdw");


/*
cout << "Ang: " << myvar.Ang << " Rey: " << myvar.Rey << endl;
cout << "noel: " << myvar.Noel << " Sigma: " << myvar.Scaling_sigma << endl;
cout << "Visc: " << myvar.Vis_str << " Prec: " << myvar.Prec << endl;
*/


 CommandLineArgs::parse_and_assign();
 CommandLineArgs::doc_specified_flags();
/*
 cout << "doc_soln: "  << doc_soln << ", doc_prec: "<< doc_prec << endl;

 cout << "Ang: " << myvar.Ang << " Rey: " << myvar.Rey << endl;
cout << "noel: " << myvar.Noel << " Sigma: " << myvar.Scaling_sigma << endl;
cout << "Visc: " << myvar.Vis_str << " Prec: " << myvar.Prec << endl;
*/

 // Set a string to identify the problem. This is unique to each problem,
 // so we hard code this. 2DStrPo = 2 dimension, straight parallel outflow.
 // straight describes the velocity flow field. Po = Parallel outflow
 // describes the boundary type.
 myvar.Prob_str = "SqTf";

 // Set the string to identify the preconditioning,
 // This is used purely for book keeping purposes.
 if(CommandLineArgs::command_line_flag_has_been_set("--w_solver"))
 {
  switch(myvar.W_solver)
  {
    case 0:
      myvar.W_str = "We";
      break;
    case 1:
      myvar.W_str = "Wa";
      break;
    default:
      std::cout << "Do not recognise W: " << myvar.W_solver << "\n"
                << "Exact preconditioning = 0\n"
                << "AMG = 1\n"
                << "Using default: Exact (W_solver = 0)"<< endl;
  }  // switch
 } // if
 // Set the string to identify the preconditioning,
 // This is used purely for book keeping purposes.
 if(CommandLineArgs::command_line_flag_has_been_set("--ns_solver"))
 {
  switch(myvar.NS_solver)
  {
    case 0:
      myvar.NS_str = "Ne";
      break;
    case 1:
      myvar.NS_str = "Nl";
      break;
    default:
      std::cout << "Do not recognise NS: " << myvar.NS_solver << "\n"
                << "Exact solve = 0\n"
                << "LSC = 1\n"
                << "Using default: LSC for NS block (NS_solver = 1)"<<std::endl;
  }  // switch
 } // if
 // Set the string to identify the preconditioning,
 // This is used purely for book keeping purposes.
 if(CommandLineArgs::command_line_flag_has_been_set("--p_solver"))
 {
  if(myvar.NS_solver == 0)
  {
    pause("NS solve is exact. There cannot be a P solver.");
  }
  switch(myvar.P_solver)
  {
    case 0:
      myvar.P_str = "Pe";
      break;
    case 1:
      myvar.P_str = "Pa";
      break;
    default:
      std::cout << "Do not recognise P: " << myvar.P_solver << "\n"
                << "Exact preconditioning = 0\n"
                << "AMG = 1\n"
                << "Using default: Exact P solve (P_solver = 0)"<< endl;
  }  // switch
 } // if
 // Set the string to identify the preconditioning,
 // This is used purely for book keeping purposes.
 if(CommandLineArgs::command_line_flag_has_been_set("--f_solver"))
 {
  if(myvar.NS_solver == 0)
  {
    pause("NS solve is exact. There cannot be an F solver.");
  }
  switch(myvar.F_solver)
  {
    case 0:
      myvar.F_str = "Fe";
      break;
    case 1:
      myvar.F_str = "Fa";
      break;
    default:
      std::cout << "Do not recognise F: " << myvar.F_solver << "\n"
                << "Exact preconditioning = 0\n"
                << "AMG = 1\n"
                << "Using default: Exact F solve (F_solver = 0)"<< endl;
  }  // switch
 } // if


 // Set the viscuous term.
 if(CommandLineArgs::command_line_flag_has_been_set("--visc"))
 {
   if (myvar.Vis == 0)
   {
     myvar.Vis_str = "Sim";
     NavierStokesEquations<2>::Gamma[0]=0.0;
     NavierStokesEquations<2>::Gamma[1]=0.0;

   }
   else if (myvar.Vis == 1)
   {
     myvar.Vis_str = "Str";
     NavierStokesEquations<2>::Gamma[0]=1.0;
     NavierStokesEquations<2>::Gamma[1]=1.0;
   } // else - setting viscuous term.
   else
   {
     std::cout << "There is no such Viscous term, using 0 = simple." 
               << std::endl; 
   }
 }

 // Set Ang_str, used for book keeping.
 if(CommandLineArgs::command_line_flag_has_been_set("--ang"))
 {
   std::ostringstream strs;
   strs << "A" << myvar.Ang;
   myvar.Ang_str = strs.str();

   // Now we need to convert Ang into radians.
   myvar.Ang = myvar.Ang * (MathematicalConstants::Pi / 180.0);

 }

 // Set Rey_str, used for book keeping.
 if(CommandLineArgs::command_line_flag_has_been_set("--rey"))
 {
   std::ostringstream strs;
   strs << "R" << myvar.Rey;
   myvar.Rey_str = strs.str();
 }

 // Set Noel_str, used for book keeping.
 if(CommandLineArgs::command_line_flag_has_been_set("--noel"))
 {
   std::ostringstream strs;
   strs << "N" << myvar.Noel;
   myvar.Noel_str = strs.str();
 }

 // Set Use_axnorm, if sigma has not been set, norm os momentum block is used.
 if(CommandLineArgs::command_line_flag_has_been_set("--sigma"))
 {
   myvar.Use_axnorm = false;

   std::ostringstream strs;
   strs << "S" << myvar.Scaling_sigma;
   myvar.Sigma_str = strs.str();
 }

 // use the diagonal or block diagonal approximation for W block.
 if(CommandLineArgs::command_line_flag_has_been_set("--bdw"))
 {
   myvar.Use_diagonal_w_block = false;
   myvar.W_approx_str = "bdw";
 }
 else
 {
   myvar.Use_diagonal_w_block = true;
   myvar.W_approx_str = "";
 }


 // Setup the label. Used for doc solution and preconditioner.
 if(myvar.NS_solver == 0)
 {
   doc_info.label() = myvar.Prob_str + myvar.W_str + myvar.NS_str + myvar.Vis_str
                      + myvar.Ang_str + myvar.Rey_str + myvar.Noel_str
                      + myvar.W_approx_str + myvar.Sigma_str;
 }
 else if(myvar.NS_solver == 1)
 {
   doc_info.label() = myvar.Prob_str
                      + myvar.W_str + myvar.NS_str + myvar.F_str + myvar.P_str
                      + myvar.Vis_str + myvar.Ang_str + myvar.Rey_str
                      + myvar.Noel_str + myvar.W_approx_str + myvar.Sigma_str;
 }
 else
 {
   pause("There is no such NS preconditioner");
 }


 // Document the solution?.
 if(CommandLineArgs::command_line_flag_has_been_set("--doc_soln"))
 {
   doc_info.enable_doc_solution();
 }
 else
 {
   doc_info.disable_doc_solution();
 }

 // Document the preconditioner?.
 if(CommandLineArgs::command_line_flag_has_been_set("--doc_prec"))
 {
   doc_info.enable_doc_prec_data();
 }
 else
 {

   doc_info.disable_doc_prec_data();
 }
 

 time_t rawtime;
 time(&rawtime);

 std::cout << "RAYDOING: "
           << doc_info.label()
           << " on " << ctime(&rawtime) << std::endl;
//*
 // Solve with Taylor-Hood element, set up problem
 TiltedCavityProblem< QTaylorHoodElement<2> > problem(myvar,doc_info);

 // Solve the problem
 problem.newton_solve();

 //Output solution
 problem.doc_solution(doc_info);


 // We now output the iteration and time.
 Vector<Vector<pair<unsigned, double> > > iters_times
   = doc_info.iterations_and_times();
/*
 for(unsigned intimestep = 0;
     intimestep < test_thing.size(); intimestep++)
 {
   for(unsigned innewtonstep = 0;
       innewtonstep < test_thing[intimestep].size();
       innewtonstep++)
   {
     cout << test_thing[intimestep][innewtonstep].first
          << " " << test_thing[intimestep][innewtonstep].second
          << endl;
   }
 }
*/

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
/*
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
 std::cout << "\t" << average_time << std::endl;

 time(&rawtime);
 std::cout << "RAYDONE: "
           << Global_Parameters::Current_settings
           << " on " << ctime (&rawtime) << endl;
*/

#ifdef OOMPH_HAS_MPI
// finalize MPI
MPI_Helpers::finalize();
#endif
 return(EXIT_SUCCESS);
} // end_of_main
