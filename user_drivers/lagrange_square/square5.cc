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

// The 2D mesh
#include "meshes/simple_rectangular_quadmesh.h"
#include "meshes/rectangular_quadmesh.h"

using namespace std;
using namespace oomph;

//===start_of_namespace=================================================
/// Namepspace for global parameters
//======================================================================       
namespace Global_Parameters
{
 /// Current Reynolds number
 double Re = 0.0;
 
 /// Current tilting angle
 double Ang = 0.0;
 
 /// Constant for pi
 const double Pi = 4.0*atan(1.0);
 
} // end of namespace

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
                  const double& lx,  const double& ly, double& phi ) :
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
}



//===start_of_problem_class=============================================
//======================================================================

template<class ELEMENT>
class TiltedCavityProblem : public Problem
{
public:

 /// \short Constructor: Pass number of elements in x and y directions and 
 /// lengths
 TiltedCavityProblem(const unsigned&);

 /// Update before solve is empty
 void actions_before_newton_solve() {}

 /// \short Update after solve is empty
 void actions_after_newton_solve() {}
 
 /// Doc the solution
 void doc_solution(DocInfo&);

 /// \short Create lagrange elements on boundary b of the Mesh pointed
 /// to by bulk_mesh_pt and add them to the Mesh object pointed to by 
 /// surface_mesh_pt
 void create_parall_outflow_lagrange_elements(const unsigned &b, 
                                              Mesh* const &bulk_mesh_pt,
                                              Mesh* const &surface_mesh_pt);

private:

  /// ID of imposed flow boundary
 unsigned Imposed_flow_boundary;

 /// ID of neumann boundary
 unsigned Neumann_boundary;

 /// Pointer to the "bulk" mesh
 SlopingQuadMesh<ELEMENT>* Bulk_mesh_pt;
 
 /// Pointer to the "surface" mesh
 Mesh* Surface_mesh_pt;

};



//===start_of_constructor=============================================
/// Problem constructor
//====================================================================
template<class ELEMENT> // rrrback - changed here.
TiltedCavityProblem<ELEMENT>::TiltedCavityProblem
(const unsigned &n_el_1d)
{

 // Assign the boundaries:
 Imposed_flow_boundary=3;
 Neumann_boundary=1;
 
 /// Setup the mesh
 // # of elements in x-direction
 unsigned nx=n_el_1d;

 // # of elements in y-direction
 unsigned ny=n_el_1d;

 // Domain length in x-direction
 double lx=1.0;

 // Domain length in y-direction
 double ly=1.0;
 Bulk_mesh_pt = 
  new SlopingQuadMesh<ELEMENT>(nx,ny,lx,ly,Global_Parameters::Ang);
 
 // Create a "surface mesh" that will contain only 
 // ImposeParallelOutflowElements in boundary 1
 // The constructor just creates the mesh without
 // giving it any elements, nodes, etc.
 Surface_mesh_pt = new Mesh;

 // Create ImposeParallelOutflowElement from all elements that are 
 // adjacent to the Neumann boundary.
 create_parall_outflow_lagrange_elements(Neumann_boundary,
                                         Bulk_mesh_pt,Surface_mesh_pt);
 
 // Add the two sub meshes to the problem
 add_sub_mesh(Bulk_mesh_pt);
 add_sub_mesh(Surface_mesh_pt);

 // Combine all submeshes into a single Mesh
 build_global_mesh();

 unsigned num_bound=Bulk_mesh_pt->nboundary();

 // Set the boundary conditions for this problem: All nodes are
 // free by default -- just pin the ones that have Dirichlet conditions
 // here.
 for(unsigned ibound=0;ibound<num_bound;ibound++)
  {
  if(ibound != Neumann_boundary)
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
 unsigned num_nod= mesh_pt()->nboundary_node(Imposed_flow_boundary);
 for (unsigned inod=0;inod<num_nod;inod++)
  {
   Node* nod_pt=mesh_pt()->boundary_node_pt(Imposed_flow_boundary,inod);
   double x=nod_pt->x(0);
   double y=nod_pt->x(1);
    
   // Tilt it back
   double ytiltedback = x*sin(-Global_Parameters::Ang)
                        +y*cos(-Global_Parameters::Ang);
   double u=0.0;
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
    
   // Now apply the rotation to u, using rotation matrices.
   // with x = u and y = 0, i.e. R*[u;0] since we have the
   // velocity in the x direction only. There is no velocity
   // in the y direction.
   double ux=u*cos(Global_Parameters::Ang);
   double uy=u*sin(Global_Parameters::Ang);
    
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
   el_pt->re_pt() = &Global_Parameters::Re;

  } // for(unsigned e=0;e<n_el;e++)

 //Assgn equation numbers
 cout << "\n equation numbers : "<< assign_eqn_numbers() << std::endl;
 
}



//==start_of_doc_solution=================================================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void TiltedCavityProblem<ELEMENT>::doc_solution(DocInfo& doc_info)
{ 
 ofstream some_file;
 char filename[100];

 // Number of plot points
 unsigned npts=5; 

 // Output solution 
 sprintf(filename,"%s/soln%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 Bulk_mesh_pt->output(some_file,npts);
 some_file.close();
 
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



//===start_of_main======================================================
/// Driver code 
//======================================================================
int main(int argc, char* argv[]) 
{
 // Set up doc info
 DocInfo doc_info;
 doc_info.number()=0;
 doc_info.set_directory("RESLT");
 
 //Doc number of gmres iterations
 ofstream out_file;

 // Setting the tilt angle in radians
 Global_Parameters::Ang = Global_Parameters::Pi/6;

 // Setting the reynolds number
 Global_Parameters::Re=0.0;

 // Setting the no_el number
 unsigned n_el_1d = 4;
         
 // Solve with Taylor-Hood element, set up problem
 TiltedCavityProblem< QTaylorHoodElement<2> > problem(n_el_1d);
         
 // Solve the problem
 problem.newton_solve();
         
 //Output solution
 problem.doc_solution(doc_info);
 doc_info.number()++;

} // end_of_main
