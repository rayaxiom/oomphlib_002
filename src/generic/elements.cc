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
//Non-inline member functions for generic elements

#include<float.h>

//oomph-lib includes
#include "elements.h"
#include "timesteppers.h"
#include "integral.h"
#include "shape.h"
#include "oomph_definitions.h"
#include "element_with_external_element.h"

namespace oomph
{

/// Static boolean to suppress warnings about repeated internal
/// data. Defaults to false
bool GeneralisedElement::Suppress_warning_about_repeated_internal_data=false;


/// Static boolean to suppress warnings about repeated external
/// data. Defaults to true
bool GeneralisedElement::Suppress_warning_about_repeated_external_data=true;


///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
//  Functions for generalised elements
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////

//=======================================================================
/// Add a (pointer to an) internal data object to the element and
/// return the index required to obtain it from the access
/// function \c internal_data_pt()
//=======================================================================
 unsigned GeneralisedElement::add_internal_data(Data* const &data_pt,
                                                const bool &fd)
 {
  //Local cache of numbers of internal and external data
  const unsigned n_internal_data = Ninternal_data;
  const unsigned n_external_data = Nexternal_data;

  //Find out whether the data is already stored in the array
  
  //Loop over the number of internals
  //The internal data are stored at the beginning of the array
  for(unsigned i=0;i<n_internal_data;i++)
   {
    //If the passed pointer is stored in the array
    if(internal_data_pt(i) == data_pt)
     {
#ifdef PARANOID
      if (!Suppress_warning_about_repeated_internal_data)
       {
        oomph_info << std::endl << std::endl;
        oomph_info 
         << "-----------------------------------------------------------------"
         << std::endl
         << "Info: Data is already included in this element's internal storage." 
         << std::endl
         << "It's stored as entry " << i <<  " and I'm not adding it again." 
         << std::endl<< std::endl
         << "Note: You can suppress this message by recompiling without"
         << "\n      PARANOID or setting the boolean \n"
         << "\n      GeneralisedElement::Suppress_warning_about_repeated_internal_data"
         << "\n\n      to true."
         << std::endl
         << "-----------------------------------------------------------------"
         << std::endl << std::endl;
       }
#endif
      //Return the index to the data object
      return i;
     }
   }
  
  //Allocate new storage for the pointers to data
  Data** new_data_pt = new Data*[n_internal_data + n_external_data + 1];

  //Copy the old iternal values across to the beginning of the array
  for(unsigned i=0;i<n_internal_data;i++) {new_data_pt[i] = Data_pt[i];}
  
  //Now add the new value to the end of the internal data
  new_data_pt[n_internal_data] = data_pt;
  
  //Copy the external values across
  for(unsigned i=0;i<n_external_data;i++)
   {new_data_pt[n_internal_data + 1 + i] = Data_pt[n_internal_data + i];}
  
  //Delete the storage associated with the previous values
  delete[] Data_pt;
  
  //Set the pointer to the new storage
  Data_pt = new_data_pt;
  
  //Resize the array of boolean flags
  Data_fd.resize(n_internal_data + n_external_data + 1);
  //Shuffle the flags for the external data to the end of the array
  for(unsigned i=n_external_data;i>0;i--)
   {
    Data_fd[n_internal_data + i] = Data_fd[n_internal_data + i-1];
   }
  //Now add the new flag to the end of the internal data
  Data_fd[n_internal_data] = fd;
  
  //Increase the number of internals
  ++Ninternal_data;
  
  //Return the final index to the new internal data
  return n_internal_data;
 }
 
//=======================================================================
/// Add the contents of the queue global_eqn_numbers to
/// the local storage for the equation numbers, which represents the 
/// local-to-global translation scheme. It is essential that the entries
/// are added in order, i.e. from the front.
//=======================================================================
 void GeneralisedElement::add_global_eqn_numbers(
  std::deque<unsigned long> const &global_eqn_numbers)
 {
  //Find the number of dofs
  const unsigned n_dof = Ndof;
  //Find the number of additional dofs
  const unsigned n_additional_dof = global_eqn_numbers.size();
  //If there are none, return immediately
  if(n_additional_dof==0) {return;}
  
  //Find the new total number of equation numbers
  const unsigned new_n_dof = n_dof + n_additional_dof;
  //Create storage for all equations, initialised to NULL
  unsigned long *new_eqn_number = new unsigned long[new_n_dof];

  //Copy over the exisiting values to the start new storage
  for(unsigned i=0;i<n_dof;i++) {new_eqn_number[i] = Eqn_number[i];}
  
  //Set an index to the next position in the new storage
  unsigned index = n_dof;
  //Loop over the queue and add it's entries to our new storage
  for(std::deque<unsigned long>::const_iterator it=global_eqn_numbers.begin();
      it!=global_eqn_numbers.end();++it)
   {
    //Add the value to the storage
    new_eqn_number[index] = *it;
    //Increase the array index
    ++index;
   }
  
  //Now delete the old storage
  delete[] Eqn_number;
  //Set the pointer to address the new storage
  Eqn_number = new_eqn_number; 
  //Finally update the number of degrees of freedom
  Ndof = new_n_dof;
 }

//========================================================================
/// Empty dense matrix used as a dummy argument to combined
/// residual and jacobian functions in the case when only the residuals
/// are being assembled
//========================================================================
DenseMatrix<double> GeneralisedElement::Dummy_matrix;
 
//=========================================================================
/// Default value used as the increment for finite difference calculations
/// of the jacobian matrices
//=========================================================================
 double GeneralisedElement::Default_fd_jacobian_step=1.0e-8;

//=======================================================================
/// Return the global time, accessed via the time pointer
//======================================================================
 double GeneralisedElement::time() const
 {
  //If no Time_pt, return 0.0
  if(Time_pt==0) {return 0.0;}
  else {return Time_pt->time();}
 }

//==========================================================================
/// Destructor for generalised elements: Wipe internal data. Pointers 
/// to external data get NULLed  but are not deleted because they 
/// are (generally) shared by lots of elements.
//==========================================================================
 GeneralisedElement::~GeneralisedElement()
 {
  //Delete each of the objects stored as internal data
  for(unsigned i=Ninternal_data;i>0;i--)
   {
    //The objects are stored at the beginning of the Data_pt array
    delete Data_pt[i-1];
    Data_pt[i-1] = 0;
   }

  //Now delete the storage for internal and external data
  delete[] Data_pt;
  
  //Now if we have allocated storage for the local equation for
  //the internal and external data, delete it.
  if(Data_local_eqn)
   {
    delete[] Data_local_eqn[0];
    delete[] Data_local_eqn;
   }

  //Delete the storage for the global equation numbers
  delete[] Eqn_number;
 }


//=======================================================================
/// Add a (pointer to an) external data object to the element and
/// return the index required to obtain it from the access
/// function \c external_data_pt()
//=======================================================================
 unsigned GeneralisedElement::add_external_data(Data* const &data_pt,
                                                const bool &fd)
 {
  //Find the numbers of internal and external data
  const unsigned n_internal_data = Ninternal_data;
  const unsigned n_external_data = Nexternal_data;
  //Find out whether the data is already stored in the array
  
  //Loop over the number of externals
  //The external data are stored at the end of the array Data_pt
  for(unsigned i=0;i<n_external_data;i++)
   {
    //If the passed pointer is stored in the array
    if(external_data_pt(i) == data_pt)
     {
#ifdef PARANOID
      if (!Suppress_warning_about_repeated_external_data)
       {
        oomph_info << std::endl << std::endl;
        oomph_info 
         << "-----------------------------------------------------------------"
         << std::endl
         << "Info: Data is already included in this element's external storage." 
         << std::endl
         << "It's stored as entry " << i <<  " and I'm not adding it again" 
         << std::endl << std::endl
         << "Note: You can suppress this message by recompiling without"
         << "\n      PARANOID or setting the boolean \n"
         << "\n      GeneralisedElement::Suppress_warning_about_repeated_external_data"
         << "\n\n      to true."
         << std::endl
         << "-----------------------------------------------------------------"
         << std::endl << std::endl;
       }
#endif
      //Return the index to the data object
      return i;
     }
   }
  
  //Allocate new storage for the pointers to data
  Data** new_data_pt = new Data*[n_internal_data + n_external_data + 1];
  
  //Copy the old iternal and external values across to the new array
  for(unsigned i=0;i<(n_internal_data + n_external_data);i++) 
   {new_data_pt[i] = Data_pt[i];}
  
  //Add the new data pointer to the end of the array
  new_data_pt[n_internal_data + n_external_data] = data_pt;
  
  //Delete the storage associated with the previous values
  delete[] Data_pt;
  
  //Set the pointer to the new storage
  Data_pt = new_data_pt;
  
  //Resize the array of boolean flags
  Data_fd.resize(n_internal_data + n_external_data + 1);
  //Now add the new flag to the end of the external data
  Data_fd[n_internal_data + n_external_data] = fd;

  //Increase the number of externals
  ++Nexternal_data;
 
  //Return the final index to the new external data
  return n_external_data;
 }

//========================================================================
/// Flush all the external data, i.e. clear the pointers to external
/// data from the internal storage.
//========================================================================
 void GeneralisedElement::flush_external_data()
 {
  //Get the numbers of internal and external data
  const unsigned n_external_data = Nexternal_data;
  //If there is external data
  if(n_external_data>0)
   {
    //Storage for new data, initialised to NULL
    Data** new_data_pt=0;

    //Find the number of internal data
    const unsigned n_internal_data = Ninternal_data;
    //If there is internal data resize Data_pt and copy over
    //the pointers
    if(n_internal_data  > 0)
     {
      //The new data pointer should only be the size of the internal data
      new_data_pt = new Data*[n_internal_data];
      //Copy over the internal data only
      for(unsigned i=0;i<n_internal_data;i++) {new_data_pt[i] = Data_pt[i];}
     }
    
    //Delete the old storage
    delete[] Data_pt;
    //Set the new storage, this will be NULL if there is no internal data
    Data_pt = new_data_pt;
    //Set the number of externals to zero
    Nexternal_data = 0;

    //Resize the array of boolean flags to the number of internals
    Data_fd.resize(n_internal_data);
   }
 }



//=========================================================================
/// Remove the object addressed by data_pt from the external data array
/// Note that this could mess up the numbering of other external data
//========================================================================
 void GeneralisedElement::flush_external_data(Data* const &data_pt)
 {
  //Get the current numbers of external data
  const unsigned n_external_data = Nexternal_data;
  //Index of the data to be removed
  //We initialise this to be n_external_data, and it will be smaller than
  //n_external_data if the data pointer is found in the array
  unsigned index = n_external_data;
  //Loop over the external data and find the argument
  for(unsigned i=0;i<n_external_data;i++)
   {
    //If we have found the data pointer, set the index and break
    if(external_data_pt(i) == data_pt)
     {
      index=i;
      break;
     }
   }

  //If we have found an index less than Nexternal_data, then we have found
  //the data in the array
  if(index < n_external_data)
   {
    //Initialise a new array to NULL
    Data** new_data_pt = 0;
    
    //Find the number of internals
    const unsigned n_internal_data = Ninternal_data;

    //Find the new number of total data items (one fewer)
    const unsigned n_total_data = n_internal_data + n_external_data - 1;
 
   //Create a new array containing the data items, if we have any
    if(n_total_data > 0) {new_data_pt = new Data*[n_total_data];}

    //Copy over all the internal data
    for(unsigned i=0;i<n_internal_data;i++) {new_data_pt[i] = Data_pt[i];}

    //Now copy over the un-flushed data
    unsigned counter=0;
    for(unsigned i=0;i<n_external_data;i++)
     {
      //If we are not at the deleted index 
      if(i!=index)
       {
        //Copy the undeleted entry into the next entry in our new 
        //array of pointers to Data
        new_data_pt[n_internal_data + counter] = Data_pt[i];
        //Increase the counter
        ++counter;
       }
     }

    //Delete the storage associated with the previous values
    delete[] Data_pt;

    //Set pointers to the new storage, will be NULL if no data left
    Data_pt = new_data_pt;

    //Remove the entry from the array of boolean flags
    Data_fd.erase(Data_fd.begin()+n_internal_data+index);

    //Decrease the number of externals
    --Nexternal_data;
    
    //Issue a warning if there will be external data remaining
    if(Nexternal_data > 1)
     {
      std::ostringstream warning_stream;
      warning_stream << "Data removed from element's external data   "
                     << std::endl
                     << "You may have to update the indices for remaning data "
                     << std::endl
                     << "This can be achieved by using add_external_data()    "
                     << std::endl;
      OomphLibWarning(warning_stream.str(),
                      "GeneralisedElement::flush_external_data()",
                      OOMPH_EXCEPTION_LOCATION);
     }
   }
 }



//==========================================================================
/// This function loops over the internal data of the element and assigns
/// GLOBAL equation numbers to the data objects.
///
/// Pass: 
/// - the current total number of dofs, global_number, which gets
///   incremented. 
/// - Dof_pt, the Vector of pointers to the global dofs 
///   (to which any internal dofs are added).
//==========================================================================
 void GeneralisedElement::
 assign_internal_eqn_numbers(unsigned long &global_number,
                             Vector<double *> &Dof_pt)
 {
  //Loop over the internal data and assign the equation numbers
  //The internal data are stored at the beginning of the Data_pt array
  for(unsigned i=0;i<Ninternal_data;i++)
   {internal_data_pt(i)->assign_eqn_numbers(global_number,Dof_pt);}
 }

 //==========================================================================
 /// This function loops over the internal data of the element and add 
 /// pointers to their unconstrained values to a map indexed by the global
 /// equation number.
 //==========================================================================
 void GeneralisedElement:: add_internal_value_pt_to_map(
  std::map<unsigned,double*> &map_of_value_pt)
 {
  //Loop over the internal data and add their data to the map
  //The internal data are stored at the beginning of the Data_pt array
  for(unsigned i=0;i<Ninternal_data;i++)
   {internal_data_pt(i)->add_value_pt_to_map(map_of_value_pt);}
 }




#ifdef OOMPH_HAS_MPI
 //=========================================================================
 /// Add all internal data and time history values to the vector in 
 /// the internal storage order
 //=========================================================================
 void GeneralisedElement::add_internal_data_values_to_vector(
  Vector<double> &vector_of_values)
 {
  for(unsigned i=0;i<Ninternal_data;i++)
   {internal_data_pt(i)->add_values_to_vector(vector_of_values);}
 }

 //========================================================================
 /// Read all internal data and time history values from the vector
 /// starting from index. On return the index will be
 /// set to the value at the end of the data that has been read in
 //========================================================================
 void GeneralisedElement::read_internal_data_values_from_vector(
  const Vector<double> &vector_of_values, 
  unsigned &index)
 {
  for(unsigned i=0;i<Ninternal_data;i++)
   {internal_data_pt(i)->read_values_from_vector(vector_of_values,index);}
 }

 //=========================================================================
 /// Add all equation numbers associated with internal data 
 /// to the vector in the internal storage order
 //=========================================================================
 void GeneralisedElement::add_internal_eqn_numbers_to_vector(
  Vector<long> &vector_of_eqn_numbers)
 {
  for(unsigned i=0;i<Ninternal_data;i++)
   {internal_data_pt(i)->add_eqn_numbers_to_vector(vector_of_eqn_numbers);}
 }

 //=========================================================================
 ///  Read all equation numbers associated with internal data
 /// from the vector
 /// starting from index. On return the index will be
 /// set to the value at the end of the data that has been read in
 //=========================================================================
 void GeneralisedElement::read_internal_eqn_numbers_from_vector(
  const Vector<long> & vector_of_eqn_numbers, unsigned &index)
 {
  for(unsigned i=0;i<Ninternal_data;i++)
   {internal_data_pt(i)->read_eqn_numbers_from_vector(vector_of_eqn_numbers,
                                                      index);}
 }
 
#endif


 //====================================================================
 /// Setup the arrays of local equation numbers for the element. 
 /// In general, this function should not need to be overloaded. Instead
 /// the function assign_all_generic_local_eqn_numbers() will be overloaded
 /// by different types of Element.
 /// That said, the function is virtual so that it 
 /// may be overloaded by the user if they *really* know what they are doing.
 //==========================================================================
 void GeneralisedElement::assign_local_eqn_numbers() 
 {
  clear_global_eqn_numbers();
  assign_all_generic_local_eqn_numbers();
  assign_additional_local_eqn_numbers();

  //Check that no global equation numbers are repeated
#ifdef PARANOID

  std::ostringstream error_stream;

  //Loop over the array of equation numbers and add to set to assess
  //uniqueness
  std::map<unsigned,bool> is_repeated;
  std::set<unsigned long> set_of_global_eqn_numbers;
  unsigned old_ndof=0;
  for(unsigned n=0;n<Ndof;++n) 
   {
    set_of_global_eqn_numbers.insert(Eqn_number[n]);
    if (set_of_global_eqn_numbers.size()==old_ndof)
     {
      error_stream << "Repeated global eqn: " << Eqn_number[n] << std::endl;
      is_repeated[Eqn_number[n]]=true;
     }
    old_ndof=set_of_global_eqn_numbers.size();
   }
  
  //If the sizes do not match we have a repeat, throw an error
  if(set_of_global_eqn_numbers.size() != Ndof)
   {
#ifdef OOMPH_HAS_MPI
    error_stream << "Element is ";
    if (!is_halo()) error_stream << "not ";
    error_stream << "a halo element\n\n";
#endif
    error_stream << "\nLocal/lobal equation numbers: " << std::endl;
    for(unsigned n=0;n<Ndof;++n)
     {
      error_stream << "  " << n << "        " << Eqn_number[n] << std::endl;
     }
    
    // It's helpful for debugging purposes to output more about this element
    error_stream << std::endl << " element address is " << this << std::endl;
    
    // Check if the repeated dofs are among the internal Data values
    unsigned nint=this->ninternal_data();
    error_stream << "Number of internal data " << nint << std::endl;
    for (unsigned i=0;i<nint;i++)
     {
      Data* data_pt=this->internal_data_pt(i);
      unsigned nval=data_pt->nvalue();
      for (unsigned j=0;j<nval;j++)
       {
        int eqn_no=data_pt->eqn_number(j);
        error_stream << "Internal dof: " << eqn_no << std::endl;
        if (is_repeated[unsigned(eqn_no)])
         {
          error_stream << "Repeated internal dof: " << eqn_no << std::endl;
         }
       }
     }
    
    // Check if the repeated dofs are among the external Data values
    unsigned next=this->nexternal_data();
    error_stream << "Number of external data " << next << std::endl;
    for (unsigned i=0;i<next;i++)
     {
      Data* data_pt=this->external_data_pt(i);
      unsigned nval=data_pt->nvalue();
      for (unsigned j=0;j<nval;j++)
       {
        int eqn_no=data_pt->eqn_number(j);
        error_stream << "External dof: " << eqn_no << std::endl;
        if (is_repeated[unsigned(eqn_no)])
         {
          error_stream << "Repeated external dof: " << eqn_no;            
          Node* nod_pt=dynamic_cast<Node*>(data_pt);
          if (nod_pt!=0)
           {
            error_stream << " (is a node at: ";
            unsigned ndim=nod_pt->ndim();
            for (unsigned ii=0;ii<ndim;ii++)
             {
              error_stream << nod_pt->x(i) << " ";
             }
           }
          error_stream << ")\n";
         }
       }
     }
    
    
    // If it's an element with external element check the associated
    // Data
    ElementWithExternalElement* e_el_pt=
     dynamic_cast<ElementWithExternalElement*>(this);
    if (e_el_pt!=0)
     {
      // Check if the repeated dofs are among the external Data values
      {
       unsigned next=e_el_pt->nexternal_interaction_field_data();
       error_stream << "Number of external element data " << next << std::endl;
       Vector<Data*> data_pt(e_el_pt->external_interaction_field_data_pt());
       for (unsigned i=0;i<next;i++)
        {
         unsigned nval=data_pt[i]->nvalue();
         for (unsigned j=0;j<nval;j++)
          {
           int eqn_no=data_pt[i]->eqn_number(j);
           error_stream << "External element dof: " << eqn_no << std::endl;
           if (is_repeated[unsigned(eqn_no)])
            {
             error_stream << "Repeated external element dof: " 
                          << eqn_no;            
             Node* nod_pt=dynamic_cast<Node*>(data_pt[i]);
             if (nod_pt!=0)
              {
               error_stream << " (is a node at: ";
               unsigned ndim=nod_pt->ndim();
               for (unsigned ii=0;ii<ndim;ii++)
                {
                 error_stream << nod_pt->x(ii) << " ";
                }
              }
             error_stream << ")\n";
            }
          }
        }
      }
       
      
      // Check if the repeated dofs are among the external geom Data values
      {
       unsigned next=e_el_pt->nexternal_interaction_geometric_data();
       error_stream << "Number of external element geom data " 
                    << next << std::endl;
       Vector<Data*> data_pt(e_el_pt->
                             external_interaction_geometric_data_pt());
       for (unsigned i=0;i<next;i++)
        {
         unsigned nval=data_pt[i]->nvalue();
         for (unsigned j=0;j<nval;j++)
          {
           int eqn_no=data_pt[i]->eqn_number(j);
           error_stream << "External element geometric dof: " 
                        << eqn_no << std::endl;
           if (is_repeated[unsigned(eqn_no)])
            {
             error_stream << "Repeated external element geometric dof: " 
                          << eqn_no << " " 
                          << data_pt[i]->value(j);            
             Node* nod_pt=dynamic_cast<Node*>(data_pt[i]);
             if (nod_pt!=0)
              {
               error_stream << " (is a node at: ";
               unsigned ndim=nod_pt->ndim();
               for (unsigned ii=0;ii<ndim;ii++)
                {
                 error_stream << nod_pt->x(i) << " ";
                }
               error_stream << ")";
              }
             error_stream << "\n";
            }
          }
        }
      }

     }

    // If it's a FiniteElement then output its nodes
    FiniteElement* f_el_pt=dynamic_cast<FiniteElement*>(this);
    if (f_el_pt!=0)
     {
      unsigned n_node=f_el_pt->nnode();
      for (unsigned n=0;n<n_node;n++)
       {
        Node* nod_pt=f_el_pt->node_pt(n);
        unsigned nval=nod_pt->nvalue();
        for (unsigned j=0;j<nval;j++)
         {
          int eqn_no=nod_pt->eqn_number(j);
          error_stream << "Node " << n 
                       << ": Nodal dof: " 
                       << eqn_no;
          if (eqn_no>=0)
           {
            if (is_repeated[unsigned(eqn_no)])
             {
              error_stream << "Node " << n 
                           << ": Repeated nodal dof: " 
                           << eqn_no;
              if (j>=f_el_pt->required_nvalue(n))
               {
                error_stream << " (resized)";
               }
              error_stream << std::endl;
             }
           }
         }       
        SolidNode* solid_nod_pt=dynamic_cast<SolidNode*>(nod_pt);
        if (solid_nod_pt!=0)
         {
          Data* data_pt=solid_nod_pt->variable_position_pt();
          unsigned nval=data_pt->nvalue();
          for (unsigned j=0;j<nval;j++)
           {
            int eqn_no=data_pt->eqn_number(j);
            error_stream <<  "Node " << n << ": Positional dof: " 
                         << eqn_no << std::endl;
            if (is_repeated[unsigned(eqn_no)])
             {
              error_stream << "Repeated positional dof: " 
                           << eqn_no << " " << data_pt->value(j) << std::endl;
             }
           }
         }        
       }
      
      // Output nodal coordinates of element
      n_node=f_el_pt->nnode();
      for (unsigned n=0;n<n_node;n++)
       {
        Node* nod_pt=f_el_pt->node_pt(n);
        unsigned n_dim=nod_pt->ndim();
        error_stream << "Node " << n << " at ( ";
        for (unsigned i=0;i<n_dim;i++)
         {
          error_stream << nod_pt->x(i) << " ";
         }
        error_stream << ")" << std::endl;
       }

     }


    throw OomphLibError(error_stream.str(),
                        "GeneralisedElement::assign_local_eqn_numbers()",
                        OOMPH_EXCEPTION_LOCATION);
   }
#endif
 }
 

//==========================================================================
/// This function loops over the internal and external data of the element, 
/// adds the GLOBAL equation numbers to the local-to-global look-up scheme and 
/// fills in the look-up schemes for the local equation 
/// numbers
//==========================================================================
 void GeneralisedElement::assign_internal_and_external_local_eqn_numbers()
 {
  //Find the number of internal and external data
  const unsigned n_internal_data = Ninternal_data;
  const unsigned n_external_data = Nexternal_data;
  //Find the total number of data
  const unsigned n_total_data = n_internal_data + n_external_data;

  //If there is data
  if(n_total_data  > 0)
   {
    //Find the number of local equations assigned so far
    unsigned local_eqn_number = ndof();

    //We need to find the total number of values stored in all the
    //internal and external data
    //Initialise to the number stored in the first data object
    unsigned n_total_values = Data_pt[0]->nvalue();
    //Loop over the other data and add the number of values stored
    for(unsigned i=1;i<n_total_data;++i)
     {n_total_values += Data_pt[i]->nvalue();}

    //If allocated delete the old storage
    if(Data_local_eqn) 
     {
      delete[] Data_local_eqn[0];
      delete[] Data_local_eqn;
     }
        
    //If there are no values then we are done, null out the storage and
    //return
    if(n_total_values==0) {Data_local_eqn=0; return;}

    //Allocate the storage for the local equation numbers
    //The idea is that we store internal equation numbers followed by
    //external equation numbers

    //Firstly allocate pointers to the rows for each data item
    Data_local_eqn = new int*[n_total_data];
    //Now allocate storage for all the equation numbers
    Data_local_eqn[0] = new int[n_total_values];
    //Set all values to be unclassified
    for(unsigned i=0;i<n_total_values;++i)
     {Data_local_eqn[0][i] = Data::Is_unclassified;}

    //Loop over the remaining rows and set their pointers
    for(unsigned i=1;i<n_total_data;++i)
     {
      //Initially set the pointer to the i-th row to the pointer
      //to the i-1th row
      Data_local_eqn[i] = Data_local_eqn[i-1];
      //Now increase the row pointer by the number of values 
      //stored at the i-1th data object
      Data_local_eqn[i] += Data_pt[i-1]->nvalue();
     }
    
    //A local queue to store the global equation numbers
    std::deque<unsigned long> global_eqn_number_queue;

    //Now loop over the internal data and assign local equation numbers
    for(unsigned i=0;i<n_internal_data;i++)
     {
      //Find the number of values stored at the internal data
      unsigned n_value = internal_data_pt(i)->nvalue();
     
      //Loop over the number of values
      for(unsigned j=0;j<n_value;j++)
       {
        //Get the GLOBAL equation number
        long eqn_number = internal_data_pt(i)->eqn_number(j);
        //If the GLOBAL equation number is positive (a free variable)
        if(eqn_number >= 0)
         {
          //Add the GLOBAL equation number to the queue
          global_eqn_number_queue.push_back(eqn_number);
          //Add the local equation number to the storage scheme
          Data_local_eqn[i][j] = local_eqn_number;
          //Increase the local number
          local_eqn_number++;
         }
        else
         {
          //Set the local scheme to be pinned
          Data_local_eqn[i][j] = Data::Is_pinned;
         } 
       }
     } //End of loop over internal data

    
    //Now loop over the external data and assign local equation numbers
    for(unsigned i=0;i<n_external_data;i++)
     {
      //Find the number of values stored at the external data
      unsigned n_value = external_data_pt(i)->nvalue();
     
      //Loop over the number of values
      for(unsigned j=0;j<n_value;j++)
       {
        //Get the GLOBAL equation number
        long eqn_number = external_data_pt(i)->eqn_number(j);
        //If the GLOBAL equation number is positive (a free variable)
        if(eqn_number >= 0)
         {
          //Add the GLOBAL equation number to the queue
          global_eqn_number_queue.push_back(eqn_number);
          //Add the local equation number to the local scheme
          Data_local_eqn[n_internal_data + i][j] = local_eqn_number;
          //Increase the local number
          local_eqn_number++;
         }
        else
         {
          //Set the local scheme to be pinned
          Data_local_eqn[n_internal_data + i][j] = Data::Is_pinned;
         }
       }
     }

    //Now add our global equations numbers to the internal element storage
    add_global_eqn_numbers(global_eqn_number_queue);
   }
 }


//============================================================================
/// This function calculates the entries of Jacobian matrix, used in 
/// the Newton method, associated with the internal degrees of freedom.
/// It does this using finite differences, 
/// rather than an analytical formulation, so can be done in total generality.
/// If the boolean argument is true, the finite
/// differencing will be performed for all internal data, irrespective of
/// the information in Data_fd. The default value (false) 
/// uses the information in Data_fd to selectively difference only
/// certain data.
//==========================================================================
 void GeneralisedElement::
 fill_in_jacobian_from_internal_by_fd(Vector<double> &residuals,
                                      DenseMatrix<double> &jacobian,
                                      const bool &fd_all_data)
 {
  //Locally cache the number of internal data
  const unsigned n_internal_data = Ninternal_data;

  //If there aren't any internal data, then return straight away
  if(n_internal_data == 0) {return;}

  //Call the update function to ensure that the element is in
  //a consistent state before finite differencing starts
  update_before_internal_fd();

  //Find the number of dofs in the element
  const unsigned n_dof = ndof();

  //Create newres vector
  Vector<double> newres(n_dof);

  //Integer storage for local unknown
  int local_unknown=0;
  
  //Use the default finite difference step
  const double fd_step = Default_fd_jacobian_step;

  //Loop over the internal data
  for(unsigned i=0;i<n_internal_data;i++)
   {
    //If we are doing all finite differences or 
    //the variable is included in the finite difference loop, do it
    if(fd_all_data || internal_data_fd(i))
     {
      //Get the number of value at the internal data
      const unsigned n_value = internal_data_pt(i)->nvalue();
      
      //Loop over the number of values
      for(unsigned j=0;j<n_value;j++)
       {
        //Get the local equation number
        local_unknown = internal_local_eqn(i,j);
        //If it's not pinned
        if(local_unknown >= 0)
         {
          //Get a pointer to the value of the internal data
          double* const value_pt = internal_data_pt(i)->value_pt(j);
          
          //Save the old value of the Internal data
          const double old_var = *value_pt;
          
          //Increment the value of the Internal data
          *value_pt += fd_step;
          
          //Now update any slaved variables
          update_in_internal_fd(i);
          
          //Calculate the new residuals
          get_residuals(newres);
          
          //Do finite differences
          for(unsigned m=0;m<n_dof;m++)
           {
            double sum = (newres[m] - residuals[m])/fd_step;
            //Stick the entry into the Jacobian matrix
            jacobian(m,local_unknown) = sum;
           }
          
          //Reset the Internal data
          *value_pt = old_var;
          
          //Reset any slaved variables
          reset_in_internal_fd(i);
         }
       }
     } //End of finite-differencing for datum (if block)
   }

  //End of finite difference loop
  //Final reset of any slaved data
  reset_after_internal_fd();
 }

//============================================================================
/// This function calculates the entries of Jacobian matrix, used in 
/// the Newton method, associated with the external degrees of freedom.
/// It does this using finite differences, 
/// rather than an analytical formulation, so can be done in total generality.
/// If the boolean argument is true, the finite
/// differencing will be performed for all internal data, irrespective of
/// the information in Data_fd. The default value (false) 
/// uses the information in Data_fd to selectively difference only
/// certain data.
//==========================================================================
 void GeneralisedElement::
 fill_in_jacobian_from_external_by_fd(Vector<double> &residuals,
                                      DenseMatrix<double> &jacobian,
                                      const bool &fd_all_data)
 {
  //Locally cache the number of external data
  const unsigned n_external_data = Nexternal_data;
  //If there aren't any external data, then return straight awayy
  if(n_external_data == 0) {return;}

  //Call the update function to ensure that the element is in
  //a consistent state before finite differencing starts
  update_before_external_fd();

  //Find the number of dofs in the element
  const unsigned n_dof = ndof();

  //Create newres vector
  Vector<double> newres(n_dof);

  //Integer storage for local unknown
  int local_unknown=0;
  
  //Use the default finite difference step
  const double fd_step = Default_fd_jacobian_step;

  //Loop over the external data
  for(unsigned i=0;i<n_external_data;i++)
   {
    //If we are doing all finite differences or 
    //the variable is included in the finite difference loop, do it
    if(fd_all_data || external_data_fd(i))
     {
      //Get the number of value at the external data
      const unsigned n_value = external_data_pt(i)->nvalue();
   
      //Loop over the number of values
      for(unsigned j=0;j<n_value;j++)
       {
        //Get the local equation number
        local_unknown = external_local_eqn(i,j);
        //If it's not pinned
        if(local_unknown >= 0)
         {
          //Get a pointer to the External data value
          double* const value_pt = external_data_pt(i)->value_pt(j);
          
          //Save the old value of the External data
          const double old_var = *value_pt;
          
          //Increment the value of the External data
          *value_pt += fd_step;
          
          //Now update any slaved variables
          update_in_external_fd(i);
          
          //Calculate the new residuals
          get_residuals(newres);
          
          //Do finite differences
          for(unsigned m=0;m<n_dof;m++)
           {
            double sum = (newres[m] - residuals[m])/fd_step;
            //Stick the entry into the Jacobian matrix
            jacobian(m,local_unknown) = sum;
           }
          
          //Reset the External data
          *value_pt = old_var;
          
          //Reset any slaved variables
          reset_in_external_fd(i);
         }
       }
     } //End of finite differencing for datum (if block)
   }

  //End of finite difference loop
  //Final reset of any slaved data
  reset_after_external_fd();
 }

 //=====================================================================
 /// \short Add the elemental contribution to the mass matrix
 /// and the residuals vector. Note that
 /// this function will NOT initialise the residuals vector or the mass
 /// matrix. It must be called after the residuals vector and 
 /// jacobian matrix have been initialised to zero. The default
 /// is deliberately broken.
 //=====================================================================
 void GeneralisedElement:: 
 fill_in_contribution_to_mass_matrix(Vector<double> &residuals,
                                     DenseMatrix<double> &mass_matrix) 
 {
  std::string error_message =
   "Empty fill_in_contribution_to_mass_matrix() has been called.\n";
  error_message +=
   "This function is called from the default implementations of\n";
  error_message += "get_mass_matrix();\n";
  error_message += 
    "and must calculate the residuals vector and mass matrix ";
  error_message += "without initialising any of their entries.\n\n";
  
  error_message += 
   "If you do not wish to use these defaults, you must overload\n";
  error_message += 
   "get_mass_matrix(), which must initialise the entries\n";
  error_message +=
   "of the residuals vector and mass matrix to zero.\n";
  
  throw 
   OomphLibError(error_message,
                 "GeneralisedElement::fill_in_contribution_to_mass_matrix()",
                 OOMPH_EXCEPTION_LOCATION);
 }

 //=====================================================================
 /// \short Add the elemental contribution to the jacobian matrix,
 /// mass matrix and the residuals vector. Note that
 /// this function should NOT initialise any entries.
 /// It must be called after the residuals vector and 
 /// matrices have been initialised to zero. The default is deliberately
 /// broken.
 //=====================================================================
 void GeneralisedElement::fill_in_contribution_to_jacobian_and_mass_matrix(
  Vector<double> &residuals,
  DenseMatrix<double> &jacobian, DenseMatrix<double> &mass_matrix)
 {
  std::string error_message =
   "Empty fill_in_contribution_to_jacobian_and_mass_matrix() has been ";
  error_message += "called.\n";
  error_message +=
   "This function is called from the default implementations of\n";
  error_message += "get_jacobian_and_mass_matrix();\n";
  error_message += 
    "and must calculate the residuals vector and mass and jacobian matrices ";
  error_message += "without initialising any of their entries.\n\n";
  
  error_message += 
   "If you do not wish to use these defaults, you must overload\n";
  error_message += 
   "get_jacobian_and_mass_matrix(), which must initialise the entries\n";
  error_message +=
   "of the residuals vector, jacobian and mass matrix to zero.\n";
  
  throw 
   OomphLibError(
    error_message,
    "GeneralisedElement::fill_in_contribution_to_jacobian_and_mass_matrix()",
    OOMPH_EXCEPTION_LOCATION);
 }


 //=====================================================================
 /// Add the elemental contribution to the derivatives of 
 /// the residuals with respect to a parameter. This function should
 /// NOT initialise any entries and must be called after the entries
 /// have been initialised to zero
 /// The default implementation is deliberately broken
 //========================================================================
 void GeneralisedElement::fill_in_contribution_to_dresiduals_dparameter(
  double* const &parameter_pt, Vector<double> &dres_dparam)
 {
  std::string error_message =
   "Empty fill_in_contribution_to_dresiduals_dparameter() has been ";
  error_message += "called.\n";
  error_message +=
   "This function is called from the default implementations of\n";
  error_message += "get_dresiduals_dparameter();\n";
  error_message += 
   "and must calculate the derivatives of the residuals vector with respect\n";
  error_message += "to the passed parameter ";
  error_message += "without initialising any values.\n\n";

  error_message += 
   "If you do not wish to use these defaults, you must overload\n";
  error_message += 
   "get_dresiduals_dparameter(), which must initialise the entries\n";
  error_message +=
   "of the returned vector to zero.\n";
  
  error_message += 
   "This function is intended for use instead of the default (global) \n";
  error_message += 
   "finite-difference routine when analytic expressions are to be used\n";
  error_message += "in continuation and bifurcation tracking problems.\n\n";
  error_message += "This function is only called when the function\n";
  error_message += 
   "Problem::set_analytic_dparameter() has been called in the driver code\n";
  
  throw 
   OomphLibError(
    error_message,
    "GeneralisedElement::fill_in_contribution_to_dresiduals_dparameter()",
    OOMPH_EXCEPTION_LOCATION);
 }

 //======================================================================
 /// Add the elemental contribution to the derivatives of 
 /// the elemental Jacobian matrix 
 /// and residuals with respect to a parameter. This function should
 /// NOT initialise any entries and must be called after the entries
 /// have been initialised to zero
 /// The default implementation is to use finite differences to calculate 
 /// the derivatives.
 //========================================================================
 void GeneralisedElement::fill_in_contribution_to_djacobian_dparameter(
  double* const &parameter_pt,
  Vector<double> &dres_dparam, 
  DenseMatrix<double> &djac_dparam)
 {
  std::string error_message =
   "Empty fill_in_contribution_to_djacobian_dparameter() has been ";
  error_message += "called.\n";
  error_message +=
   "This function is called from the default implementations of\n";
  error_message += "get_djacobian_dparameter();\n";
  error_message += 
   "and must calculate the derivatives of residuals vector and jacobian ";
   error_message += "matrix\n";
  error_message += "with respect to the passed parameter ";
  error_message += "without initialising any values.\n\n";

  error_message += 
   "If you do not wish to use these defaults, you must overload\n";
  error_message += 
   "get_djacobian_dparameter(), which must initialise the entries\n";
  error_message +=
   "of the returned vector and matrix to zero.\n\n";
  
  error_message += 
   "This function is intended for use instead of the default (global) \n";
  error_message += 
   "finite-difference routine when analytic expressions are to be used\n";
  error_message += "in continuation and bifurcation tracking problems.\n\n";
  error_message += "This function is only called when the function\n";
  error_message += 
   "Problem::set_analytic_dparameter() has been called in the driver code\n";

  
  throw 
   OomphLibError(
    error_message,
    "GeneralisedElement::fill_in_contribution_to_djacobian_dparameter()",
    OOMPH_EXCEPTION_LOCATION);
 }


 //=====================================================================
 /// \short Add the elemental contribution to the derivative of the 
 /// jacobian matrix,
 /// mass matrix and the residuals vector with respect to a parameter. 
 /// Note that
 /// this function should NOT initialise any entries.
 /// It must be called after the residuals vector and 
 /// matrices have been initialised to zero. The default is deliberately
 /// broken.
 //=====================================================================
 void GeneralisedElement::
 fill_in_contribution_to_djacobian_and_dmass_matrix_dparameter(
  double* const &parameter_pt,
  Vector<double> &dres_dparam,
  DenseMatrix<double> &djac_dparam, 
  DenseMatrix<double> &dmass_matrix_dparam)
 {
  std::string error_message =
   "Empty fill_in_contribution_to_djacobian_and_dmass_matrix_dparameter() has";
  error_message += " been called.\n";
  error_message +=
   "This function is called from the default implementations of\n";
  error_message += "get_djacobian_and_dmass_matrix_dparameter();\n";
  error_message += 
    "and must calculate the residuals vector and mass and jacobian matrices ";
  error_message += "without initialising any of their entries.\n\n";
  
  error_message += 
   "If you do not wish to use these defaults, you must overload\n";
  error_message += 
   "get_djacobian_and_dmass_matrix_dparameter(), which must initialise the\n";
    error_message += 
   "entries of the returned vector and  matrices to zero.\n";


  error_message += 
   "This function is intended for use instead of the default (global) \n";
  error_message += 
   "finite-difference routine when analytic expressions are to be used\n";
  error_message += "in continuation and bifurcation tracking problems.\n\n";
  error_message += "This function is only called when the function\n";
  error_message += 
   "Problem::set_analytic_dparameter() has been called in the driver code\n";

  

  throw 
   OomphLibError(
    error_message,
    "GeneralisedElement::fill_in_contribution_to_djacobian_and_dmass_matrix_dparameter()",
    OOMPH_EXCEPTION_LOCATION);
 }



 //========================================================================
 /// Fill in contribution to the product of the Hessian 
 /// (derivative of Jacobian with
 /// respect to all variables) an eigenvector, Y, and 
 /// other specified vectors, C
 /// (d(J_{ij})/d u_{k}) Y_{j} C_{k}
 /// At the moment the dof pointer is passed in to enable
 /// easy calculation of finite difference default
//==========================================================================
 void GeneralisedElement::fill_in_contribution_to_hessian_vector_products(
  Vector<double> const &Y,
  DenseMatrix<double> const &C,
  DenseMatrix<double> &product)
 {
  std::string error_message =
   "Empty fill_in_contribution_to_hessian_vector_products() has been ";
  error_message += "called.\n";
  error_message +=
   "This function is called from the default implementations of\n";
  error_message += "get_hessian_vector_products(); ";
  error_message += " and must calculate the products \n";
  error_message += "of the hessian matrix with the (global) ";
  error_message += "vectors Y and C\n";
  error_message += "without initialising any values.\n\n";

  error_message += 
   "If you do not wish to use these defaults, you must overload\n";
  error_message += 
   "get_hessian_vector_products(), which must initialise the entries\n";
  error_message +=
   "of the returned vector to zero.\n\n";
  
  error_message += 
   "This function is intended for use instead of the default (global) \n";
  error_message += 
   "finite-difference routine when analytic expressions are to be used.\n";
  error_message += "This function is only called when the function\n";
  error_message += 
   "Problem::set_analytic_hessian_products() has been called in the driver code\n";
  
  throw 
   OomphLibError(
    error_message,
    "GeneralisedElement::fill_in_contribution_to_hessian_vector_product()",
    OOMPH_EXCEPTION_LOCATION);
 }


//==========================================================================
/// Self-test: Have all internal values been classified as 
/// pinned/unpinned? Return 0 if OK. 
//==========================================================================
 unsigned GeneralisedElement::self_test()
 {
  // Initialise
  bool passed=true;
 
  //Loop over internal Data
  for(unsigned i=0;i<Ninternal_data;i++)
   {
    if(internal_data_pt(i)->self_test()!=0)
     {
      passed=false;
      oomph_info 
       << "\n ERROR: Failed GeneralisedElement::self_test()!" << std::endl;
      oomph_info 
       << "for internal data object number: " << i << std::endl;
     }
   }

  //Loop over external Data
  for(unsigned i=0;i<Nexternal_data;i++)
   {
    if(external_data_pt(i)->self_test()!=0)
     {
      passed=false;
      oomph_info 
       << "\n ERROR: Failed GeneralisedElement::self_test()!" << std::endl;
      oomph_info 
       << "for external data object number: " << i << std::endl;
     }
   }

  // Return verdict
  if (passed) {return 0;}
  else {return 1;}
 }


//======================================================================
/// Helper namespace for tolerances, number of iterations, etc
/// used in the locate_zeta function in FiniteElement
//======================================================================
namespace Locate_zeta_helpers
{
 /// Convergence tolerance for the newton solver
 double Newton_tolerance = 1.0e-7;	
 
 /// Maximum number of newton iterations
 unsigned Max_newton_iterations = 10;
 
 /// Rounding tolerance for whether coordinate is in element or not
 double Rounding_tolerance = 1.0e-12;

 /// Number of points along one dimension of each element used
 /// to create coordinates within the element in order to see
 /// which has the smallest initial residual (and is therefore used
 /// as the initial guess in the Newton method when locating coordinate)
 unsigned N_local_points = 2;
}


///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
//  Functions for finite elements
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////



//========================================================================
/// \short Internal function used to check for singular or negative values
/// of the determinant of the Jacobian of the mapping between local and 
/// global or lagrangian coordinates. Negative jacobians are allowed if the 
/// Accept_negative_jacobian flag is set to true. 
//========================================================================
 void FiniteElement::check_jacobian(const double &jacobian) const
 {
  //First check for a zero jacobian
  if(std::fabs(jacobian) < 1.0e-16)
   {
    if (FiniteElement::Suppress_output_while_checking_for_inverted_elements)
     {
      throw OomphLibQuietException();
     }
    else
     {
      throw OomphLibError(
       "Determinant of Jacobian matrix is zero --- singular mapping!",
       "FiniteElement::check_jacobian()",
       OOMPH_EXCEPTION_LOCATION);
     }
   }
  //Now check for negative jacobians, if we're not allowing them (default)
  if((Accept_negative_jacobian==false) && (jacobian < 0.0))
   {

    if (FiniteElement::Suppress_output_while_checking_for_inverted_elements)
     {
      throw OomphLibQuietException();
     }
    else
     {
      std::ostringstream error_stream;
      error_stream 
       << "Negative Jacobian in transform from "
       << "local to global coordinates"
       << std::endl 
       << "       You have an inverted coordinate system!" 
       << std::endl << std::endl;
      error_stream 
       << "Here are the nodal coordinates of the inverted element\n" 
       << "in the form \n\n       x,y[,z], hang_status\n\n" 
       << "where hang_status = 1 or 2 for non-hanging or hanging\n" 
       << "nodes respectively (useful for automatic sizing of\n" 
       << "tecplot markers to identify the hanging nodes). \n\n" ;
      unsigned n_dim_nod=node_pt(0)->ndim();
      unsigned n_nod=nnode();
      unsigned hang_count=0;
      for (unsigned j=0;j<n_nod;j++)
       {
        for (unsigned i=0;i<n_dim_nod;i++)
         {
          error_stream << node_pt(j)->x(i) << " ";
         }
        if (node_pt(j)->is_hanging())
         {
          error_stream << " 2";
          hang_count++;
         }
        else
         {
          error_stream << " 1";
         }
        error_stream << std::endl;
       }
      error_stream << std::endl << std::endl;
      if ((Macro_elem_pt!=0)&&(0!=hang_count))
       {
        error_stream 
         << "NOTE: The inverted element is associated with a MacroElement\n"
         << "       AND the element has hanging nodes! \n"
         << "       If an element is thin and highly curved, the \n"
         << "       constraints imposed by\n \n"
         << "       (1) inter-element continuity (imposed by the hanging\n"
         << "           node constraints) and \n\n"
         << "       (2) the need to respect curvilinear domain boundaries\n"
         << "           during mesh refinement (imposed by the element's\n"
         << "           macro element mapping)\n\n"
         << "       may be irreconcilable!  \n \n"
         << "       You may have to re-design your base mesh to avoid \n"
         << "       the creation of thin, highly curved elements during\n"
         << "       the refinement process.\n"
         << std::endl;
       }
      
      error_stream 
       << std::endl << std::endl
       << "If you believe that inverted elements do not cause any\n" 
       << "problems in your specific application you can \n " 
       << "suppress this test by: " << std::endl
       << "  i) setting the (static) flag " 
       << "FiniteElement::Accept_negative_jacobian to be true" << std::endl;
      error_stream << " ii) switching OFF the PARANOID flag" 
                   << std::endl << std::endl;
      
      throw OomphLibError(error_stream.str(),
                          "FiniteElement::check_jacobian()",
                          OOMPH_EXCEPTION_LOCATION);
     }
   }
 }

//=========================================================================
/// Internal function that is used to assemble the jacobian of the mapping
/// from local coordinates (s) to the eulerian coordinates (x), given the
/// derivatives of the shape functions. 
//=========================================================================
void FiniteElement::
assemble_local_to_eulerian_jacobian(const DShape &dpsids,
                                    DenseMatrix<double> &jacobian) const
{
 //Locally cache the elemental dimension
 const unsigned el_dim = dim();
 //The number of shape functions must be equal to the number
 //of nodes (by definition)
 const unsigned n_shape = nnode();
 //The number of shape function types must be equal to the number
 //of nodal position types (by definition)
 const unsigned n_shape_type = nnodal_position_type();

#ifdef PARANOID
 //Check for dimensional compatibility
 if(Elemental_dimension != Nodal_dimension)
  {
   std::ostringstream error_message;
   error_message << "Dimension mismatch" << std::endl;
   error_message << "The elemental dimension: " << Elemental_dimension
                 << " must equal the nodal dimension: " 
                 << Nodal_dimension 
                 << " for the jacobian of the mapping to be well-defined"
                 << std::endl;
   throw OomphLibError(error_message.str(),
                       "FiniteElement::assemble_local_to_eulerian_jacobian()",
                       OOMPH_EXCEPTION_LOCATION);
  }
#endif
   

 //Loop over the rows of the jacobian
 for(unsigned i=0;i<el_dim;i++)
  {
   //Loop over the columns of the jacobian
   for(unsigned j=0;j<el_dim;j++)
    {
     //Zero the entry
     jacobian(i,j) = 0.0;
     //Loop over the shape functions
     for(unsigned l=0;l<n_shape;l++)
      {
       for(unsigned k=0;k<n_shape_type;k++)
        {
         //Jacobian is dx_j/ds_i, which is represented by the sum
         //over the dpsi/ds_i of the nodal points X j
         //Call the Non-hanging version of positions
         //This is overloaded in refineable elements
         jacobian(i,j) += raw_nodal_position_gen(l,k,j)*dpsids(l,k,i);
        }
      }
    }
  }
}

//=========================================================================
/// Internal function that is used to assemble the jacobian of second
/// derivatives of the the mapping from local coordinates (s) to the 
/// eulerian coordinates (x), given the second derivatives of the 
/// shape functions. 
//=========================================================================
 void FiniteElement::
 assemble_local_to_eulerian_jacobian2(const DShape &d2psids,
                                      DenseMatrix<double> &jacobian2) const
 {
  //Find the the dimension of the element
  const unsigned el_dim = dim();
  //Find the number of shape functions and shape functions types
  //Must be equal to the number of nodes and their position types
  //by the definition of the shape function.
  const unsigned n_shape = nnode();
  const unsigned n_shape_type = nnodal_position_type();
  //Find the number of second derivatives
  const unsigned n_row = N2deriv[el_dim];
 
  //Assemble the "jacobian" (d^2 x_j/ds_i^2) for second derivatives of 
  //shape functions
  //Loop over the rows (number of second derivatives)
  for(unsigned i=0;i<n_row;i++)
   {
    //Loop over the columns (element dimension
    for(unsigned j=0;j<el_dim;j++)
     {
      //Zero the entry
      jacobian2(i,j) = 0.0;
      //Loop over the shape functions
      for(unsigned l=0;l<n_shape;l++)
       {
        //Loop over the shape function types
        for(unsigned k=0;k<n_shape_type;k++)
         {
          //Add the terms to the jacobian entry
          //Call the Non-hanging version of positions
          //This is overloaded in refineable elements
          jacobian2(i,j) += raw_nodal_position_gen(l,k,j)*d2psids(l,k,i);
         }
       }
     }
   }
 }

 //=====================================================================
 /// Assemble the covariant Eulerian base vectors and return them in the
 /// matrix interpolated_G. The derivatives of the shape functions with
 /// respect to the local coordinate should already have been calculated
 /// before calling this function
 //=====================================================================
 void FiniteElement::assemble_eulerian_base_vectors(
  const DShape &dpsids, DenseMatrix<double> &interpolated_G) const
 {
  //Find the number of nodes and position types
  const unsigned n_node = nnode();
  const unsigned n_position_type = nnodal_position_type();
  //Find the dimension of the node and element
  const unsigned n_dim_node = nodal_dimension();
  const unsigned n_dim_element = dim();

  //Loop over the dimensions and compute the entries of the 
  //base vector matrix
  for(unsigned i=0;i<n_dim_element;i++)
   {
    for(unsigned j=0;j<n_dim_node;j++)
     {
      //Initialise the j-th component of the i-th base vector to zero
      interpolated_G(i,j) = 0.0;
      for(unsigned l=0;l<n_node;l++)
       {
        for(unsigned k=0;k<n_position_type;k++)
         {
          interpolated_G(i,j) += raw_nodal_position_gen(l,k,j)*dpsids(l,k,i);
         }
       }
     }
   }
 }


//============================================================================
/// Zero-d specialisation of function to calculate inverse of jacobian mapping
//============================================================================
 template<>
 double FiniteElement::
 invert_jacobian<0>(const DenseMatrix<double> &jacobian,
                    DenseMatrix<double> &inverse_jacobian) const
 {
  //Issue a warning
  oomph_info << "\nWarning: You are trying to invert the jacobian for "
             << "a 'point' element" << std::endl
             << "This makes no sense and is almost certainly an error"
             << std::endl << std::endl;

  //Dummy return
  return(1.0);
 }


//===========================================================================
/// One-d specialisation of function to calculate inverse of jacobian mapping
//===========================================================================
 template<>
 double FiniteElement::
 invert_jacobian<1>(const DenseMatrix<double> &jacobian,
                    DenseMatrix<double> &inverse_jacobian) const
 {
  //Calculate the determinant of the matrix
  const double det = jacobian(0,0);

//Report if Matrix is singular or negative
#ifdef PARANOID
  check_jacobian(det);
#endif 

  //Calculate the inverse --- trivial in 1D
  inverse_jacobian(0,0) = 1.0/jacobian(0,0);
 
  //Return the determinant
  return(det);
 }

//===========================================================================
/// Two-d specialisation of function to calculate inverse of jacobian mapping
//===========================================================================
 template<>
 double FiniteElement::
 invert_jacobian<2>(const DenseMatrix<double> &jacobian,
                    DenseMatrix<double> &inverse_jacobian) const
 {
  //Calculate the determinant of the matrix
  const double det = jacobian(0,0)*jacobian(1,1) - jacobian(0,1)*jacobian(1,0);

//Report if Matrix is singular or negative
#ifdef PARANOID
  check_jacobian(det);
#endif 

  //Calculate the inverset of the 2x2 matrix
  inverse_jacobian(0,0) =  jacobian(1,1)/det;
  inverse_jacobian(0,1) = -jacobian(0,1)/det;
  inverse_jacobian(1,0) = -jacobian(1,0)/det;
  inverse_jacobian(1,1) =  jacobian(0,0)/det;

  //Return the jacobian
  return(det);
 }

//=============================================================================
/// Three-d specialisation of function to calculate inverse of jacobian mapping
//=============================================================================
 template<>
 double FiniteElement::
 invert_jacobian<3>(const DenseMatrix<double> &jacobian,
                    DenseMatrix<double> &inverse_jacobian) const
 {
  //Calculate the determinant of the matrix
  const double det = jacobian(0,0)*jacobian(1,1)*jacobian(2,2) 
   + jacobian(0,1)*jacobian(1,2)*jacobian(2,0) 
   + jacobian(0,2)*jacobian(1,0)*jacobian(2,1) 
   - jacobian(0,0)*jacobian(1,2)*jacobian(2,1) 
   - jacobian(0,1)*jacobian(1,0)*jacobian(2,2) 
   - jacobian(0,2)*jacobian(1,1)*jacobian(2,0); 
 
  //Report if Matrix is singular or negative
#ifdef PARANOID
  check_jacobian(det);
#endif

  //Calculate the inverse of the 3x3 matrix
  inverse_jacobian(0,0) =  (jacobian(1,1)*jacobian(2,2) 
                            - jacobian(1,2)*jacobian(2,1))/det; 
  inverse_jacobian(0,1) = -(jacobian(0,1)*jacobian(2,2) 
                            - jacobian(0,2)*jacobian(2,1))/det; 
  inverse_jacobian(0,2) =  (jacobian(0,1)*jacobian(1,2) 
                            - jacobian(0,2)*jacobian(1,1))/det; 
  inverse_jacobian(1,0) = -(jacobian(1,0)*jacobian(2,2) 
                            - jacobian(1,2)*jacobian(2,0))/det; 
  inverse_jacobian(1,1) =  (jacobian(0,0)*jacobian(2,2) 
                            - jacobian(0,2)*jacobian(2,0))/det; 
  inverse_jacobian(1,2) = -(jacobian(0,0)*jacobian(1,2) 
                            - jacobian(0,2)*jacobian(1,0))/det; 
  inverse_jacobian(2,0) =  (jacobian(1,0)*jacobian(2,1) 
                            - jacobian(1,1)*jacobian(2,0))/det; 
  inverse_jacobian(2,1) = -(jacobian(0,0)*jacobian(2,1) 
                            - jacobian(0,1)*jacobian(2,0))/det; 
  inverse_jacobian(2,2) =  (jacobian(0,0)*jacobian(1,1) 
                            - jacobian(0,1)*jacobian(1,0))/det; 

  //Return the determinant
  return(det);
 }

//========================================================================
/// Template-free interface for inversion of the jacobian of a mapping.
/// This is slightly inefficient, given that it uses a switch statement.
/// It can always be overloaded in specific geometric elements, 
/// for efficiency reasons.
//========================================================================
 double FiniteElement::
 invert_jacobian_mapping(const DenseMatrix<double> &jacobian,
                         DenseMatrix<double> &inverse_jacobian) const
 {
  //Find the spatial dimension of the element
  const unsigned el_dim = dim();
  //Call the appropriate templated function, depending on the 
  //element dimension
  switch(el_dim)
   {
   case 0:
    return invert_jacobian<0>(jacobian,inverse_jacobian);
    break;
   case 1:
    return invert_jacobian<1>(jacobian,inverse_jacobian);
    break;
   case 2:
    return invert_jacobian<2>(jacobian,inverse_jacobian);
    break;
   case 3:
    return invert_jacobian<3>(jacobian,inverse_jacobian);
    break;
    //Catch-all default case: issue warning and die
   default:
    std::ostringstream error_stream;
    error_stream 
     << "Dimension of the element must be 0,1,2 or 3, not " << el_dim 
     << std::endl;
  
    throw OomphLibError(error_stream.str(),
                        "FiniteElement::invert_jacobian_mapping(..)",
                        OOMPH_EXCEPTION_LOCATION);
   }
  //Dummy return for Intel compiler
  return 1.0;
 }

//============================================================================
/// Zero-d specialisation of function to calculate the derivative of the
/// jacobian of a mapping with respect to the nodal coordinates X_ij.
//============================================================================
 template<>
 void FiniteElement::dJ_eulerian_dnodal_coordinates_templated_helper<0>(
  const DenseMatrix<double> &jacobian,const DShape &dpsids,
  DenseMatrix<double> &djacobian_dX) const
 {
  // Issue a warning
  oomph_info << "\nWarning: You are trying to calculate derivatives of "
             << "a jacobian w.r.t. nodal coordinates for a 'point' "
             << "element." << std::endl
             << "This makes no sense and is almost certainly an error."
             << std::endl << std::endl;
 }

//===========================================================================
/// One-d specialisation of function to calculate the derivative of the
/// jacobian of a mapping with respect to the nodal coordinates X_ij.
//===========================================================================
 template<>
 void FiniteElement::dJ_eulerian_dnodal_coordinates_templated_helper<1>(
  const DenseMatrix<double> &jacobian,const DShape &dpsids,
  DenseMatrix<double> &djacobian_dX) const
 {
  // Determine the number of nodes in the element
  const unsigned n_node = nnode();
  
  // Loop over nodes
  for(unsigned j=0;j<n_node;j++)
   {
    djacobian_dX(0,j) = dpsids(j,0);
   }
 }

//===========================================================================
/// Two-d specialisation of function to calculate the derivative of the
/// jacobian of a mapping with respect to the nodal coordinates X_ij.
//===========================================================================
 template<>
 void FiniteElement::dJ_eulerian_dnodal_coordinates_templated_helper<2>(
  const DenseMatrix<double> &jacobian,const DShape &dpsids,
  DenseMatrix<double> &djacobian_dX) const
 {
  // Determine the number of nodes in the element
  const unsigned n_node = nnode();

  // Loop over nodes
  for(unsigned j=0;j<n_node;j++)
   {
    // i=0
    djacobian_dX(0,j) = dpsids(j,0)*jacobian(1,1) - dpsids(j,1)*jacobian(0,1);

    // i=1
    djacobian_dX(1,j) = dpsids(j,1)*jacobian(0,0) - dpsids(j,0)*jacobian(1,0);
   }
 }

//=============================================================================
/// Three-d specialisation of function to calculate the derivative of the
/// jacobian of a mapping with respect to the nodal coordinates X_ij.
//=============================================================================
 template<>
 void FiniteElement::dJ_eulerian_dnodal_coordinates_templated_helper<3>(
  const DenseMatrix<double> &jacobian,const DShape &dpsids,
  DenseMatrix<double> &djacobian_dX) const
 {
  // Determine the number of nodes in the element
  const unsigned n_node = nnode();

  // Loop over nodes
  for(unsigned j=0;j<n_node;j++)
   {
    // i=0
    djacobian_dX(0,j)
     = dpsids(j,0)*(jacobian(1,1)*jacobian(2,2)
                    - jacobian(1,2)*jacobian(2,1))
     + dpsids(j,1)*(jacobian(0,2)*jacobian(2,1)
                    - jacobian(0,1)*jacobian(2,2))
     + dpsids(j,2)*(jacobian(0,1)*jacobian(1,2)
                    - jacobian(0,2)*jacobian(1,1));
    
    // i=1
    djacobian_dX(1,j)
     = dpsids(j,0)*(jacobian(1,2)*jacobian(2,0)
                    - jacobian(1,0)*jacobian(2,2))
     + dpsids(j,1)*(jacobian(0,0)*jacobian(2,2)
                    - jacobian(0,2)*jacobian(2,0))
     + dpsids(j,2)*(jacobian(0,2)*jacobian(1,0)
                    - jacobian(0,0)*jacobian(1,2));

    // i=2
    djacobian_dX(2,j)
     = dpsids(j,0)*(jacobian(1,0)*jacobian(2,1)
                    - jacobian(1,1)*jacobian(2,0))
     + dpsids(j,1)*(jacobian(0,1)*jacobian(2,0)
                    - jacobian(0,0)*jacobian(2,1))
     + dpsids(j,2)*(jacobian(0,0)*jacobian(1,1)
                    - jacobian(0,1)*jacobian(1,0));
   }
 }

//============================================================================
/// Zero-d specialisation of function to calculate the derivative w.r.t. the
/// nodal coordinates \f$ X_{pq} \f$ of the derivative of the shape functions
/// w.r.t. the global eulerian coordinates \f$ x_i \f$.
//============================================================================
 template<>
 void FiniteElement::d_dshape_eulerian_dnodal_coordinates_templated_helper<0>(
  const double &det_jacobian,
  const DenseMatrix<double> &jacobian,
  const DenseMatrix<double> &djacobian_dX,
  const DenseMatrix<double> &inverse_jacobian,
  const DShape &dpsids,
  RankFourTensor<double> &d_dpsidx_dX) const
 {
  // Issue a warning
  oomph_info << "\nWarning: You are trying to calculate derivatives of "
             << "eulerian derivatives of shape functions w.r.t. nodal "
             << "coordinates for a 'point' element." << std::endl
             << "This makes no sense and is almost certainly an error."
             << std::endl << std::endl;
 }

//===========================================================================
/// One-d specialisation of function to calculate the derivative w.r.t. the
/// nodal coordinates \f$ X_{pq} \f$ of the derivative of the shape functions
/// w.r.t. the global eulerian coordinates \f$ x_i \f$.
//===========================================================================
 template<>
 void FiniteElement::d_dshape_eulerian_dnodal_coordinates_templated_helper<1>(
  const double &det_jacobian,
  const DenseMatrix<double> &jacobian,
  const DenseMatrix<double> &djacobian_dX,
  const DenseMatrix<double> &inverse_jacobian,
  const DShape &dpsids,
  RankFourTensor<double> &d_dpsidx_dX) const
 {
  // Find inverse of determinant of jacobian of mapping
  const double inv_det_jac = 1.0/det_jacobian;

  // Determine the number of nodes in the element
  const unsigned n_node = nnode();

  // Loop over the shape functions
  for(unsigned q=0;q<n_node;q++)
   {
    // Loop over the shape functions
    for(unsigned j=0;j<n_node;j++)
     {
      d_dpsidx_dX(0,q,j,0)
       = - djacobian_dX(0,q)*dpsids(j,0)*inv_det_jac*inv_det_jac;
     }
   }
 }

//===========================================================================
/// Two-d specialisation of function to calculate the derivative w.r.t. the
/// nodal coordinates \f$ X_{pq} \f$ of the derivative of the shape functions
/// w.r.t. the global eulerian coordinates \f$ x_i \f$.
//===========================================================================
 template<>
 void FiniteElement::d_dshape_eulerian_dnodal_coordinates_templated_helper<2>(
  const double &det_jacobian,
  const DenseMatrix<double> &jacobian,
  const DenseMatrix<double> &djacobian_dX,
  const DenseMatrix<double> &inverse_jacobian,
  const DShape &dpsids,
  RankFourTensor<double> &d_dpsidx_dX) const
 {
  // Find inverse of determinant of jacobian of mapping
  const double inv_det_jac = 1.0/det_jacobian;

  // Determine the number of nodes in the element
  const unsigned n_node = nnode();

  // Loop over the spatial dimension (this must be 2)
  for(unsigned p=0;p<2;p++)
   {
    // Loop over the shape functions
    for(unsigned q=0;q<n_node;q++)
     {
      // Loop over the shape functions
      for(unsigned j=0;j<n_node;j++)
       {
        // i=0
        d_dpsidx_dX(p,q,j,0) = - djacobian_dX(p,q)*
         (inverse_jacobian(0,0)*dpsids(j,0)
          + inverse_jacobian(0,1)*dpsids(j,1));
        
        if(p==1)
         {
          d_dpsidx_dX(p,q,j,0)
           += dpsids(j,0)*dpsids(q,1) - dpsids(j,1)*dpsids(q,0);
         }
        d_dpsidx_dX(p,q,j,0) *= inv_det_jac;

        // i=1
        d_dpsidx_dX(p,q,j,1) = - djacobian_dX(p,q)*
         (inverse_jacobian(1,1)*dpsids(j,1)
          + inverse_jacobian(1,0)*dpsids(j,0));
        
        if(p==0)
         {
          d_dpsidx_dX(p,q,j,1)
           += dpsids(j,1)*dpsids(q,0) - dpsids(j,0)*dpsids(q,1);
         }
        d_dpsidx_dX(p,q,j,1) *= inv_det_jac;
       }
     }
   }
 }

//=============================================================================
/// Three-d specialisation of function to calculate the derivative w.r.t. the
/// nodal coordinates \f$ X_{pq} \f$ of the derivative of the shape functions
/// w.r.t. the global eulerian coordinates \f$ x_i \f$.
//=============================================================================
 template<>
 void FiniteElement::d_dshape_eulerian_dnodal_coordinates_templated_helper<3>(
  const double &det_jacobian,
  const DenseMatrix<double> &jacobian,
  const DenseMatrix<double> &djacobian_dX,
  const DenseMatrix<double> &inverse_jacobian,
  const DShape &dpsids,
  RankFourTensor<double> &d_dpsidx_dX) const
 {
  // Find inverse of determinant of jacobian of mapping
  const double inv_det_jac = 1.0/det_jacobian;
  
  // Determine the number of nodes in the element
  const unsigned n_node = nnode();
  
  // Loop over the spatial dimension (this must be 3)
  for(unsigned p=0;p<3;p++)
   {
    // Loop over the shape functions
    for(unsigned q=0;q<n_node;q++)
     {
      // Loop over the shape functions
      for(unsigned j=0;j<n_node;j++)
       {
        // Terms not multiplied by delta function
        for(unsigned i=0;i<3;i++)
         {
          d_dpsidx_dX(p,q,j,i)
           = - djacobian_dX(p,q)*(inverse_jacobian(i,0)*dpsids(j,0)
                                  + inverse_jacobian(i,1)*dpsids(j,1)
                                  + inverse_jacobian(i,2)*dpsids(j,2));
         }

        // Delta function terms
        switch(p)
         {
         case 0:
          d_dpsidx_dX(p,q,j,1)+=((dpsids(q,2)*jacobian(1,2)
                                  - dpsids(q,1)*jacobian(2,2))*dpsids(j,0)
                                 + (dpsids(q,0)*jacobian(2,2)
                                    - dpsids(q,2)*jacobian(0,2))*dpsids(j,1)
                                 + (dpsids(q,1)*jacobian(0,2)
                                    - dpsids(q,0)*jacobian(1,2))*dpsids(j,2));

          d_dpsidx_dX(p,q,j,2)+=((dpsids(q,1)*jacobian(2,1)
                                  - dpsids(q,2)*jacobian(1,1))*dpsids(j,0)
                                 + (dpsids(q,2)*jacobian(0,1)
                                    - dpsids(q,0)*jacobian(2,1))*dpsids(j,1)
                                 + (dpsids(q,0)*jacobian(1,1)
                                    - dpsids(q,1)*jacobian(0,1))*dpsids(j,2));
          break;

         case 1:
        
          d_dpsidx_dX(p,q,j,0)+=((dpsids(q,1)*jacobian(2,2)
                                  - dpsids(q,2)*jacobian(1,2))*dpsids(j,0)
                                 + (dpsids(q,2)*jacobian(0,2)
                                    - dpsids(q,0)*jacobian(2,2))*dpsids(j,1)
                                 + (dpsids(q,0)*jacobian(1,2)
                                    - dpsids(q,1)*jacobian(0,2))*dpsids(j,2));

          d_dpsidx_dX(p,q,j,2)+=((dpsids(q,2)*jacobian(1,0)
                                  - dpsids(q,1)*jacobian(2,0))*dpsids(j,0)
                                 + (dpsids(q,0)*jacobian(2,0)
                                    - dpsids(q,2)*jacobian(0,0))*dpsids(j,1)
                                 + (dpsids(q,1)*jacobian(0,0)
                                    - dpsids(q,0)*jacobian(1,0))*dpsids(j,2));
          break;
         
         case 2:

          d_dpsidx_dX(p,q,j,0)+=((dpsids(q,2)*jacobian(1,1)
                                  - dpsids(q,1)*jacobian(2,1))*dpsids(j,0)
                                 + (dpsids(q,0)*jacobian(2,1)
                                    - dpsids(q,2)*jacobian(0,1))*dpsids(j,1)
                                 + (dpsids(q,1)*jacobian(0,1)
                                    - dpsids(q,0)*jacobian(1,1))*dpsids(j,2));

          d_dpsidx_dX(p,q,j,1)+=((dpsids(q,1)*jacobian(2,0)
                                  - dpsids(q,2)*jacobian(1,0))*dpsids(j,0)
                                 + (dpsids(q,2)*jacobian(0,0)
                                    - dpsids(q,0)*jacobian(2,0))*dpsids(j,1)
                                 + (dpsids(q,0)*jacobian(1,0)
                                    - dpsids(q,1)*jacobian(0,0))*dpsids(j,2));
          break;
         }
        
        // Divide through by the determinant of the Jacobian mapping
        for(unsigned i=0;i<3;i++)
         {
          d_dpsidx_dX(p,q,j,i) *= inv_det_jac;
         }
       }
     }
   }
 }

//=======================================================================
/// Default value for the number of values at a node
//=======================================================================
 const unsigned FiniteElement::Default_Initial_Nvalue = 0;

//======================================================================
/// \short Default value that is used for the tolerance required when 
/// locating nodes via local coordinates
 const double FiniteElement::Node_location_tolerance = 1.0e-14;

//======================================================================
/// Set the default value of the Accept_negative_jacobian flag to be 
/// false
//======================================================================
 bool FiniteElement::Accept_negative_jacobian = false;


//======================================================================
///  Set default for static boolean to suppress output while checking 
/// for inverted elements
//======================================================================
bool FiniteElement::Suppress_output_while_checking_for_inverted_elements=false;

//========================================================================
/// Static array that holds the number of rows in the second derivative
/// matrix as a function of spatial dimension. In one-dimension, there is
/// only one possible second derivative. In two-dimensions, there are three,
/// the two second derivatives and the mixed derivatives. In three 
/// dimensions there are six.
//=========================================================================
 const unsigned FiniteElement::N2deriv[4]={0,1,3,6};

//==========================================================================
/// Calculate the mapping from local to eulerian coordinates
/// assuming that the coordinates are aligned in the direction of the local 
/// coordinates, i.e. there are no cross terms and the jacobian is diagonal.
/// The local derivatives are passed as dpsids and the jacobian and 
/// inverse jacobian are returned.
//==========================================================================
 double FiniteElement::
 local_to_eulerian_mapping_diagonal(const DShape &dpsids,
                                    DenseMatrix<double> &jacobian,
                                    DenseMatrix<double> &inverse_jacobian)
  const
 {
  //Find the dimension of the element
  const unsigned el_dim = dim();
  //Find the number of shape functions and shape functions types
  //Equal to the number of nodes and their position types by definition
  const unsigned n_shape = nnode();
  const unsigned n_shape_type = nnodal_position_type();

#ifdef PARANOID
  //Check for dimension compatibility
 if(Elemental_dimension != Nodal_dimension)
  {
   std::ostringstream error_message;
   error_message << "Dimension mismatch" << std::endl;
   error_message << "The elemental dimension: " << Elemental_dimension
                 << " must equal the nodal dimension: " 
                 << Nodal_dimension 
                 << " for the jacobian of the mapping to be well-defined"
                 << std::endl;
   throw OomphLibError(error_message.str(),
                       "FiniteElement::local_to_eulerian_jacobian_diagonal()",
                       OOMPH_EXCEPTION_LOCATION);
  }
#endif
 
  //In this case we assume that there are no cross terms, that is
  //global coordinate 0 is always in the direction of local coordinate 0
 
  //Loop over the coordinates
  for(unsigned i=0;i<el_dim;i++)
   {
    //Zero the jacobian and inverse jacobian entries
    for(unsigned j=0;j<el_dim;j++) 
     {jacobian(i,j) = 0.0; inverse_jacobian(i,j) = 0.0;}
   
    //Loop over the shape functions
    for(unsigned l=0;l<n_shape;l++)
     {
      //Loop over the types of dof
      for(unsigned k=0;k<n_shape_type;k++)
       {
        //Derivatives are always dx_{i}/ds_{i}
        jacobian(i,i) += raw_nodal_position_gen(l,k,i)*dpsids(l,k,i);
       }
     }
   }
 
  //Now calculate the determinant of the matrix
  double det = 1.0;
  for(unsigned i=0;i<el_dim;i++) {det *= jacobian(i,i);}
 
//Report if Matrix is singular, or negative
#ifdef PARANOID
  check_jacobian(det);
#endif
 
  //Calculate the inverse mapping (trivial in this case)
  for(unsigned i=0;i<el_dim;i++) 
   {inverse_jacobian(i,i) = 1.0/jacobian(i,i);}
 
  //Return the value of the Jacobian
  return(det);
 }

//========================================================================
/// Template-free interface calculating the derivative of the jacobian
/// of a mapping with respect to the nodal coordinates X_ij. This is
/// slightly inefficient, given that it uses a switch statement. It can
/// always be overloaded in specific geometric elements, for efficiency
/// reasons.
//========================================================================
void FiniteElement::dJ_eulerian_dnodal_coordinates(
 const DenseMatrix<double> &jacobian,const DShape &dpsids,
 DenseMatrix<double> &djacobian_dX) const
{
 // Determine the spatial dimension of the element
 const unsigned el_dim = dim();
 
#ifdef PARANOID
 // Determine the number of nodes in the element
 const unsigned n_node = nnode();
 
 // Check that djacobian_dX has the correct number of rows (= el_dim)
 if(djacobian_dX.nrow()!=el_dim)
  {
   std::ostringstream error_message;
   error_message << "djacobian_dX must have the same number of rows as the"
                 << "\nspatial dimension of the element.";
   throw OomphLibError(error_message.str(),
                       "FiniteElement::dJ_eulerian_dnodal_coordinates()",
                       OOMPH_EXCEPTION_LOCATION);
  }
 // Check that djacobian_dX has the correct number of columns (= n_node)
 if(djacobian_dX.ncol()!=n_node)
  {
   std::ostringstream error_message;
   error_message << "djacobian_dX must have the same number of columns as the"
                 << "\nnumber of nodes in the element.";
   throw OomphLibError(error_message.str(),
                       "FiniteElement::dJ_eulerian_dnodal_coordinates()",
                       OOMPH_EXCEPTION_LOCATION);
  }
#endif

 // Call the appropriate templated function, depending on the 
 // element dimension
 switch(el_dim)
  {
  case 0:
   dJ_eulerian_dnodal_coordinates_templated_helper<0>(jacobian,dpsids,
                                                      djacobian_dX);
   break;
  case 1:
   dJ_eulerian_dnodal_coordinates_templated_helper<1>(jacobian,dpsids,
                                                      djacobian_dX);
   break;
  case 2:
   dJ_eulerian_dnodal_coordinates_templated_helper<2>(jacobian,dpsids,
                                                      djacobian_dX);
   break;
  case 3:
   dJ_eulerian_dnodal_coordinates_templated_helper<3>(jacobian,dpsids,
                                                      djacobian_dX);
   break;
   // Catch-all default case: issue warning and die
  default:
   std::ostringstream error_stream;
   error_stream 
    << "Dimension of the element must be 0,1,2 or 3, not " << el_dim 
    << std::endl;
   
   throw OomphLibError(error_stream.str(),
                       "FiniteElement::dJ_eulerian_dnodal_coordinates(..)",
                       OOMPH_EXCEPTION_LOCATION);
  }
}

//========================================================================
/// Template-free interface calculating the derivative w.r.t. the nodal
/// coordinates \f$ X_{pq} \f$ of the derivative of the shape functions
/// \f$ \psi_j \f$ w.r.t. the global eulerian coordinates \f$ x_i \f$.
/// I.e. this function calculates
/// \f[
/// \frac{\partial}{\partial X_{pq}}
/// \left( \frac{\partial \psi_j}{\partial x_i} \right).
/// \f]
/// To do this it requires the determinant of the jacobian mapping, its
/// derivative w.r.t. the nodal coordinates \f$ X_{pq} \f$, the inverse
/// jacobian and the derivatives of the shape functions w.r.t. the local
/// coordinates. The result is returned as a tensor of rank four.
/// \n\n Numbering: \n
/// d_dpsidx_dX(p,q,j,i) = \f$ \frac{\partial}{\partial X_{pq}}
/// \left( \frac{\partial \psi_j}{\partial x_i} \right) \f$ \n
/// This function is slightly inefficient, given that it uses a switch
/// statement. It can always be overloaded in specific geometric elements,
/// for efficiency reasons.
//========================================================================
void FiniteElement::d_dshape_eulerian_dnodal_coordinates(
 const double &det_jacobian,
 const DenseMatrix<double> &jacobian,
 const DenseMatrix<double> &djacobian_dX,
 const DenseMatrix<double> &inverse_jacobian,
 const DShape &dpsids,
 RankFourTensor<double> &d_dpsidx_dX) const
{
 // Determine the spatial dimension of the element
 const unsigned el_dim = dim();

#ifdef PARANOID
 // Determine the number of nodes in the element
 const unsigned n_node = nnode();
 
 // Check that d_dpsidx_dX is of the correct size
 if(d_dpsidx_dX.nindex1()!=el_dim || d_dpsidx_dX.nindex2()!=n_node
    || d_dpsidx_dX.nindex3()!=n_node || d_dpsidx_dX.nindex4()!=el_dim)
  {
   std::ostringstream error_message;
   error_message << "d_dpsidx_dX must be of the following dimensions:"
                 << "\nd_dpsidx_dX(el_dim,n_node,n_node,el_dim)";
   throw OomphLibError(error_message.str(),
                       "FiniteElement::d_dshape_eulerian_dnodal_coordinates(...)",
                       OOMPH_EXCEPTION_LOCATION);
  }
#endif

 // Call the appropriate templated function, depending on the 
 // element dimension
 switch(el_dim)
  {
  case 0:
   d_dshape_eulerian_dnodal_coordinates_templated_helper<0>(
    det_jacobian,jacobian,djacobian_dX,inverse_jacobian,dpsids,d_dpsidx_dX);
   break;
  case 1:
   d_dshape_eulerian_dnodal_coordinates_templated_helper<1>(
    det_jacobian,jacobian,djacobian_dX,inverse_jacobian,dpsids,d_dpsidx_dX);
   break;
  case 2:
   d_dshape_eulerian_dnodal_coordinates_templated_helper<2>(
    det_jacobian,jacobian,djacobian_dX,inverse_jacobian,dpsids,d_dpsidx_dX);
   break;
  case 3:
   d_dshape_eulerian_dnodal_coordinates_templated_helper<3>(
    det_jacobian,jacobian,djacobian_dX,inverse_jacobian,dpsids,d_dpsidx_dX);
   break;
   // Catch-all default case: issue warning and die
  default:
   std::ostringstream error_stream;
   error_stream 
    << "Dimension of the element must be 0,1,2 or 3, not " << el_dim 
    << std::endl;
   
   throw OomphLibError(error_stream.str(),
                       "FiniteElement::d_dshape_eulerian_dnodal_coordinates(...)",
                       OOMPH_EXCEPTION_LOCATION);
  }
}

//=====================================================================
/// Convert derivatives w.r.t local coordinates to derivatives w.r.t the
/// coordinates used to assemble the inverse jacobian mapping passed as
/// inverse_jacobian. The derivatives passed in dbasis will be
/// modified in this function from dbasisds to dbasisdX.
//======================================================================
 void FiniteElement::transform_derivatives(const DenseMatrix<double> 
                                           &inverse_jacobian,
                                           DShape &dbasis) const
 {
  //Find the number of basis functions and basis function types
  const unsigned n_basis = dbasis.nindex1();
  const unsigned n_basis_type = dbasis.nindex2();
  //Find the dimension of the array (Must be the elemental dimension)
  const unsigned n_dim = dim();
 
  //Allocate temporary (stack) storage of the dimension of the element
  double new_derivatives[n_dim];

  //Loop over the number of basis functions
  for(unsigned l=0;l<n_basis;l++)
   {
    //Loop over type of basis functions
    for(unsigned k=0;k<n_basis_type;k++)
     {
      //Loop over the coordinates
      for(unsigned j=0;j<n_dim;j++)
       {
        //Zero the new transformed derivatives
        new_derivatives[j] = 0.0;
        //Do premultiplication by inverse jacobian
        for(unsigned i=0;i<n_dim;i++)
         {
          new_derivatives[j] += inverse_jacobian(j,i)*dbasis(l,k,i);
         }
       }
      //We now copy the new derivatives into the shape functions
      for(unsigned j=0;j<n_dim;j++) {dbasis(l,k,j) = new_derivatives[j];}
     }
   }
 }

//=====================================================================
/// Convert derivatives w.r.t local coordinates to derivatives w.r.t the
/// coordinates used to assemble the inverse jacobian mapping passed as
/// inverse_jacobian, assuming that the mapping is diagonal. This merely
/// saves a few loops, but is probably worth it.
//======================================================================
 void FiniteElement::transform_derivatives_diagonal(const DenseMatrix<double> 
                                                    &inverse_jacobian,
                                                    DShape &dbasis) const
 {
  //Find the number of basis functions and basis function types
  const unsigned n_basis = dbasis.nindex1();
  const unsigned n_basis_type = dbasis.nindex2();
  //Find the dimension of the array (must be the elemental dimension)
  const unsigned n_dim = dim();
 
  //Loop over the number of basis functions
  for(unsigned l=0;l<n_basis;l++)
   {
    //Loop over type of basis functions
    for(unsigned k=0;k<n_basis_type;k++)
     {
      //Loop over the coordinates
      for(unsigned j=0;j<n_dim;j++)
       {
        dbasis(l,k,j) *= inverse_jacobian(j,j);
       }
     }
   }
 }

//=======================================================================
/// Convert derivatives and second derivatives w.r.t local coordinates to 
/// derivatives w.r.t. the coordinates used to assemble the jacobian,
/// inverse_jacobian and jacobian 2 passed. This must be specialised
/// for each dimension, otherwise it gets very ugly
/// Specialisation to one dimension.
//=======================================================================
 template<>
 void FiniteElement::
 transform_second_derivatives_template<1>(const DenseMatrix<double>
                                          &jacobian,
                                          const DenseMatrix<double> 
                                          &inverse_jacobian,
                                          const DenseMatrix<double> 
                                          &jacobian2,
                                          DShape &dbasis,
                                          DShape &d2basis) const
 {
  //Find the number of basis functions and basis function types
  const unsigned n_basis = dbasis.nindex1();
  const unsigned n_basis_type = dbasis.nindex2();
  
  //The second derivatives are easy, because there can be no mixed terms
  //Loop over number of basis functions
  for(unsigned l=0;l<n_basis;l++) 
   {
    //Loop over number of basis function types
    for(unsigned k=0;k<n_basis_type;k++)
     {
      d2basis(l,k,0) = d2basis(l,k,0)/(jacobian(0,0)*jacobian(0,0))
       //Second term comes from taking d/dx of (dpsi/ds / dx/ds)
       - dbasis(l,k,0)*jacobian2(0,0)/
       (jacobian(0,0)*jacobian(0,0)*jacobian(0,0));
     }
   }

  //Assemble the first derivatives (do this last so that we don't
  //overwrite the dphids before we use it in the above)
  transform_derivatives(inverse_jacobian,dbasis);

 }

//=======================================================================
/// Convert derivatives and second derivatives w.r.t local coordinates to 
/// derivatives w.r.t. the coordinates used to assemble the jacobian,
/// inverse_jacobian and jacobian 2 passed. This must be specialised
/// for each dimension, otherwise it gets very ugly.
/// Specialisation to two spatial dimensions
//=======================================================================
 template<>
 void FiniteElement::
 transform_second_derivatives_template<2>(const DenseMatrix<double>
                                          &jacobian,
                                          const DenseMatrix<double> 
                                          &inverse_jacobian,
                                          const DenseMatrix<double> 
                                          &jacobian2,
                                          DShape &dbasis,
                                          DShape &d2basis) const
 {
  //Find the number of basis functions and basis function types
  const unsigned n_basis = dbasis.nindex1();
  const unsigned n_basis_type = dbasis.nindex2();

  //Calculate the determinant
  const double det = jacobian(0,0)*jacobian(1,1) - jacobian(0,1)*jacobian(1,0);

  //Second derivatives ... the approach taken here is to construct
  //dphidX/ds which can then be used to calculate the second derivatives
  //using the relationship d/dX = inverse_jacobian*d/ds

  double ddetds[2];

  ddetds[0] = jacobian2(0,0)*jacobian(1,1) + jacobian(0,0)*jacobian2(2,1) 
   - jacobian2(0,1)*jacobian(1,0) - jacobian(0,1)*jacobian2(2,0);
  ddetds[1] = jacobian2(2,0)*jacobian(1,1) + jacobian(0,0)*jacobian2(1,1)
   - jacobian2(2,1)*jacobian(1,0) - jacobian(0,1)*jacobian2(1,0);

  //Calculate the derivatives of the inverse jacobian wrt the local coordinates
  double dinverse_jacobiands[2][2][2];

  dinverse_jacobiands[0][0][0] = jacobian2(2,1)/det - 
   jacobian(1,1)*ddetds[0]/(det*det);
  dinverse_jacobiands[0][1][0] = -jacobian2(0,1)/det +
   jacobian(0,1)*ddetds[0]/(det*det);
  dinverse_jacobiands[1][0][0] = -jacobian2(2,0)/det +
   jacobian(1,0)*ddetds[0]/(det*det);
  dinverse_jacobiands[1][1][0] = jacobian2(0,0)/det -
   jacobian(0,0)*ddetds[0]/(det*det);

  dinverse_jacobiands[0][0][1] = jacobian2(1,1)/det - 
   jacobian(1,1)*ddetds[1]/(det*det);
  dinverse_jacobiands[0][1][1] = -jacobian2(2,1)/det +
   jacobian(0,1)*ddetds[1]/(det*det);
  dinverse_jacobiands[1][0][1] = -jacobian2(1,0)/det +
   jacobian(1,0)*ddetds[1]/(det*det);
  dinverse_jacobiands[1][1][1] = jacobian2(2,0)/det -
   jacobian(0,0)*ddetds[1]/(det*det);

  //Set up derivatives of dpsidx wrt local coordinates
  DShape dphidXds0(n_basis,n_basis_type,2), dphidXds1(n_basis,n_basis_type,2);
  
  for(unsigned l=0;l<n_basis;l++)
   {
    for(unsigned k=0;k<n_basis_type;k++)
     {
      for(unsigned j=0;j<2;j++)
       {
        //Note that we can't have an inner loop because of the
        //convention I've chosen for the mixed derivatives!
        dphidXds0(l,k,j) = dinverse_jacobiands[j][0][0]*dbasis(l,k,0)
         + dinverse_jacobiands[j][1][0]*dbasis(l,k,1)
         + inverse_jacobian(j,0)*d2basis(l,k,0)
         + inverse_jacobian(j,1)*d2basis(l,k,2);
       
        dphidXds1(l,k,j) = dinverse_jacobiands[j][0][1]*dbasis(l,k,0)
         + dinverse_jacobiands[j][1][1]*dbasis(l,k,1)
         + inverse_jacobian(j,0)*d2basis(l,k,2)
         + inverse_jacobian(j,1)*d2basis(l,k,1);
       }
     }
   }
  
  //Now calculate the DShape d2phidX
  for(unsigned l=0;l<n_basis;l++)
   {
    for(unsigned k=0;k<n_basis_type;k++)
     {
      //Zero dpsidx 
      for(unsigned j=0;j<3;j++) {d2basis(l,k,j) = 0.0;}
      
      //Do premultiplication by inverse jacobian
      for(unsigned i=0;i<2;i++)
       {
        d2basis(l,k,i) = inverse_jacobian(i,0)*dphidXds0(l,k,i)+
         inverse_jacobian(i,1)*dphidXds1(l,k,i);
       }
      //Calculate mixed derivative term
      d2basis(l,k,2) += inverse_jacobian(0,0)*dphidXds0(l,k,1)
       + inverse_jacobian(0,1)*dphidXds1(l,k,1);
     }
   }
 
  //Assemble the first derivatives second, so that we don't
  //overwrite dphids
  transform_derivatives(inverse_jacobian,dbasis);
 }


//=======================================================================
/// Convert derivatives and second derivatives w.r.t local coordinates to 
/// derivatives w.r.t. the coordinates used to assemble the jacobian,
/// inverse_jacobian and jacobian 2 passed. This must be specialised
/// for each dimension, otherwise it gets very ugly
/// Specialisation to one dimension.
//=======================================================================
 template<>
 void FiniteElement::
 transform_second_derivatives_diagonal<1>(const DenseMatrix<double>
                                          &jacobian,
                                          const DenseMatrix<double> 
                                          &inverse_jacobian,
                                          const DenseMatrix<double> 
                                          &jacobian2,
                                          DShape &dbasis,
                                          DShape &d2basis) const
 {
  FiniteElement::
   transform_second_derivatives_template<1>(jacobian,
                                            inverse_jacobian,
                                            jacobian2,dbasis,d2basis);
 }


//=========================================================================
/// Convert second derivatives w.r.t. local coordinates to 
/// second derivatives w.r.t. the coordinates passed in the tensor
/// coordinate. Specialised to two spatial dimension
//=========================================================================
 template<>
 void FiniteElement::
 transform_second_derivatives_diagonal<2>(const DenseMatrix<double>
                                          &jacobian,
                                          const DenseMatrix<double> 
                                          &inverse_jacobian,
                                          const DenseMatrix<double> 
                                          &jacobian2,
                                          DShape &dbasis,
                                          DShape &d2basis) const
 {
  //Find the number of basis functions and basis function types
  const unsigned n_basis = dbasis.nindex1();
  const unsigned n_basis_type = dbasis.nindex2();
 
  //Again we assume that there are no cross terms and that  coordinate
  //i depends only upon local coordinate i

  //Now calculate the DShape d2phidx
  for(unsigned l=0;l<n_basis;l++)
   {
    for(unsigned k=0;k<n_basis_type;k++)
     {
      //Second derivatives
      d2basis(l,k,0) = d2basis(l,k,0)/(jacobian(0,0)*jacobian(0,0))
       - dbasis(l,k,0)*jacobian2(0,0)/
       (jacobian(0,0)*jacobian(0,0)*jacobian(0,0));

      d2basis(l,k,1) = d2basis(l,k,1)/(jacobian(1,1)*jacobian(1,1))
       - dbasis(l,k,1)*jacobian2(1,1)/
       (jacobian(1,1)*jacobian(1,1)*jacobian(1,1));
     
      d2basis(l,k,2) = d2basis(l,k,2)/(jacobian(0,0)*jacobian(1,1));
     }
   }


  //Assemble the first derivatives
  transform_derivatives_diagonal(inverse_jacobian,dbasis);
 }


//=============================================================================
/// \short Convert derivatives and second derivatives w.r.t. local coordiantes
/// to derivatives and second derivatives w.r.t. the coordinates used to
/// assemble the jacobian, inverse jacobian and jacobian2 passed to the
/// function. This is a template-free general interface, that should be
/// overloaded for efficiency
//============================================================================
 void FiniteElement::transform_second_derivatives(const DenseMatrix<double>
                                                  &jacobian,
                                                  const DenseMatrix<double> 
                                                  &inverse_jacobian,
                                                  const DenseMatrix<double> 
                                                  &jacobian2,
                                                  DShape &dbasis,
                                                  DShape &d2basis) const
 {
  //Find the dimension of the element
  const unsigned el_dim = dim();
  //Choose the appropriate function based on the dimension of the element
  switch(el_dim)
   {
   case 1:
    transform_second_derivatives_template<1>(jacobian,
                                             inverse_jacobian,jacobian2,
                                             dbasis,d2basis);
    break;
   case 2:
    transform_second_derivatives_template<2>(jacobian,
                                             inverse_jacobian,jacobian2,
                                             dbasis,d2basis);
    break;
         
   case 3:
    throw OomphLibError("Not implemented yet ... maybe one day",
                        "FiniteElement::transform_second_derivative()",
                        OOMPH_EXCEPTION_LOCATION);

    //transform_second_derivatives_template<3>(dphids,d2phids,jacobian,
    //                                         inverse_jacobian,jacobian2,
    //                                         dphidX,d2phidX);
    break;
    //Catch-all default case: issue warning and die
   default:
    std::ostringstream error_stream;
    error_stream
     << "Dimension of the element must be 1,2 or 3, not " << el_dim
     << std::endl;

    throw OomphLibError(error_stream.str(),
                        "FiniteElement::transform_second_derivatives(..)",
                        OOMPH_EXCEPTION_LOCATION);
   }
 }


//======================================================================
/// \short The destructor cleans up the memory allocated
/// for storage of pointers to nodes. Internal and external data get
/// wiped by the GeneralisedElement destructor; nodes get
/// killed in mesh destructor. 
//=======================================================================
 FiniteElement::~FiniteElement() 
 {
  //Delete the storage allocated for the pointers to the loca nodes
  delete[] Node_pt;
  
  //Delete the storage allocated for the nodal numbering schemes
  if(Nodal_local_eqn)
   {
    delete[] Nodal_local_eqn[0];
    delete[] Nodal_local_eqn;
   }
 }


 //==============================================================
 /// Get the local fraction of the node j in the element; 
 /// vector sets its own size
 //==============================================================
 void FiniteElement::local_fraction_of_node(const unsigned &j, 
                                            Vector<double> &s_fraction)
 {
  //Default implementation is rather dumb
  //Get the local coordinate and scale by local coordinate range
  local_coordinate_of_node(j,s_fraction);
  unsigned n_coordinates = s_fraction.size();
  for(unsigned i=0;i<n_coordinates;i++)
   {
    s_fraction[i] = (s_fraction[i] - s_min())/(s_max() - s_min());
   }
 }

//=======================================================================
/// Set the spatial integration scheme and also calculate the values of the
/// shape functions and their derivatives w.r.t. the local coordinates,
/// placing the values into storage so that they may be re-used, 
/// without recalculation
//=======================================================================
 void FiniteElement::set_integration_scheme(Integral* const &integral_pt) 
 {
  //Assign the integration scheme
  Integral_pt = integral_pt;
 }

//=========================================================================
/// \short Return the shape function stored at the ipt-th integration
/// point.
//=========================================================================
 void FiniteElement::shape_at_knot(const unsigned &ipt, Shape &psi) const
 {
  //Find the dimension of the element
  const unsigned el_dim = dim();
  //Storage for the local coordinates of the integration point
  Vector<double> s(el_dim); 
  //Set the local coordinate
  for(unsigned i=0;i<el_dim;i++) {s[i] = integral_pt()->knot(ipt,i);}
  //Get the shape function
  shape(s,psi);
 }

//=========================================================================
/// \short Return the shape function and its derivatives w.r.t. the local
/// coordinates at the ipt-th integration point.
//=========================================================================
 void FiniteElement::dshape_local_at_knot(const unsigned &ipt, Shape &psi,
                                          DShape &dpsids) const
 {
  //Find the dimension of the element
  const unsigned el_dim = dim();
  //Storage for the local coordinates of the integration point
  Vector<double> s(el_dim); 
  //Set the local coordinate
  for(unsigned i=0;i<el_dim;i++) {s[i] = integral_pt()->knot(ipt,i);}
  //Get the shape function and derivatives
  dshape_local(s,psi,dpsids);
 }

//=========================================================================
/// Calculate the shape function and its first and second derivatives 
/// w.r.t. local coordinates at the ipt-th integration point.
/// \n\n Numbering:
/// \n \b 1D: \n
/// d2psids(i,0) = \f$ d^2 \psi_j / d s^2 \f$
/// \n \b 2D: \n
/// d2psids(i,0) = \f$ \partial^2 \psi_j / \partial s_0^2 \f$ \n
/// d2psids(i,1) = \f$ \partial^2 \psi_j / \partial s_1^2 \f$ \n
/// d2psids(i,2) = \f$ \partial^2 \psi_j / \partial s_0 \partial s_1 \f$ \n
/// \n \b 3D: \n
/// d2psids(i,0) = \f$ \partial^2 \psi_j / \partial s_0^2 \f$ \n
/// d2psids(i,1) = \f$ \partial^2 \psi_j / \partial s_1^2 \f$ \n
/// d2psids(i,2) = \f$ \partial^2 \psi_j / \partial s_2^2 \f$ \n
/// d2psids(i,3) = \f$ \partial^2 \psi_j / \partial s_0 \partial s_1 \f$ \n
/// d2psids(i,4) = \f$ \partial^2 \psi_j / \partial s_0 \partial s_2 \f$ \n
/// d2psids(i,5) = \f$ \partial^2 \psi_j / \partial s_1 \partial s_2 \f$ \n
//=========================================================================
 void FiniteElement::d2shape_local_at_knot(const unsigned &ipt, Shape &psi, 
                                           DShape &dpsids, DShape &d2psids) 
  const
 {
  //Find the dimension of the element
  const unsigned el_dim = dim();
  //Storage for the local coordinates of the integration point
  Vector<double> s(el_dim); 
  //Set the local coordinate
  for(unsigned i=0;i<el_dim;i++) {s[i] = integral_pt()->knot(ipt,i);}
  //Get the shape function and first and second derivatives
  d2shape_local(s,psi,dpsids,d2psids);
 }

//=========================================================================
/// \short Compute the geometric shape functions and also
/// first derivatives w.r.t. global coordinates at local coordinate s;
/// Returns Jacobian of mapping from global to local coordinates.
/// Most general form of the function, but may be over-loaded, if desired
//=========================================================================
 double FiniteElement::dshape_eulerian(const Vector<double> &s, 
                                       Shape &psi,
                                       DShape &dpsi) const
 {
  //Find the element dimension
  const unsigned el_dim = dim();
 
  //Get the values of the shape functions and their local derivatives
  //Temporarily stored in dpsi
  dshape_local(s,psi,dpsi);
  
  //Allocate memory for the inverse jacobian
  DenseMatrix<double> inverse_jacobian(el_dim);
  //Now calculate the inverse jacobian
  const double det = local_to_eulerian_mapping(dpsi,inverse_jacobian);
  
  //Now set the values of the derivatives to be dpsidx
  transform_derivatives(inverse_jacobian,dpsi);
  //Return the determinant of the jacobian
  return det;
 }
 
//========================================================================
/// \short Compute the geometric shape functions and also first
/// derivatives w.r.t. global coordinates at integration point ipt.
/// Most general form of function, but may be over-loaded if desired
//========================================================================
 double FiniteElement::dshape_eulerian_at_knot(const unsigned &ipt,
                                               Shape &psi, 
                                               DShape &dpsi) const
 {
  //Find the element dimension
  const unsigned el_dim = dim();
 
  //Get the values of the shape function and local derivatives
  //Temporarily store it in dpsi
  dshape_local_at_knot(ipt,psi,dpsi);
 
  //Allocate memory for the inverse jacobian
  DenseMatrix<double> inverse_jacobian(el_dim);
  //Now calculate the inverse jacobian
  const double det = local_to_eulerian_mapping(dpsi,inverse_jacobian);

  //Now set the values of the derivatives to dpsidx
  transform_derivatives(inverse_jacobian,dpsi);
  //Return the determinant of the jacobian
  return det;
 }


//========================================================================
/// \short Compute the geometric shape functions (psi) and first
/// derivatives w.r.t. global coordinates (dpsidx) at integration point
/// ipt. Return the determinant of the jacobian of the mapping (detJ).
/// Additionally calculate the derivatives of both "detJ" and "dpsidx"
/// w.r.t. the nodal coordinates.
/// Most general form of function, but may be over-loaded if desired.
//========================================================================
 double FiniteElement::dshape_eulerian_at_knot(
  const unsigned &ipt,
  Shape &psi, 
  DShape &dpsi,
  DenseMatrix<double> &djacobian_dX,
  RankFourTensor<double> &d_dpsidx_dX) const
 {
  // Find the element dimension
  const unsigned el_dim = dim();
 
  // Get the values of the shape function and local derivatives
  // Temporarily store in dpsi
  dshape_local_at_knot(ipt,psi,dpsi);
 
  // Allocate memory for the jacobian and the inverse of the jacobian
  DenseMatrix<double> jacobian(el_dim), inverse_jacobian(el_dim);

  // Now calculate the inverse jacobian
  const double det = local_to_eulerian_mapping(dpsi,jacobian,inverse_jacobian);

  // Calculate the derivative of the jacobian w.r.t. nodal coordinates
  // Note: must call this before "transform_derivatives(...)" since this
  // function requires dpsids rather than dpsidx
  dJ_eulerian_dnodal_coordinates(jacobian,dpsi,djacobian_dX);

  // Calculate the derivative of dpsidx w.r.t. nodal coordinates
  // Note: this function also requires dpsids rather than dpsidx
  d_dshape_eulerian_dnodal_coordinates(det,jacobian,djacobian_dX,
                                       inverse_jacobian,dpsi,d_dpsidx_dX);
                                                            
  // Now set the values of the derivatives to dpsidx
  transform_derivatives(inverse_jacobian,dpsi);

  // Return the determinant of the jacobian
  return det;
 }



//===========================================================================
/// \short Compute the geometric shape functions and also first
/// and second derivatives w.r.t. global coordinates at local coordinate s;
/// Also returns Jacobian of mapping from global to local coordinates.
/// \n\n Numbering:
/// \n \b 1D: \n
/// d2psidx(i,0) = \f$ d^2 \psi_j / d x^2 \f$
/// \n \b 2D: \n
/// d2psidx(i,0) = \f$ \partial^2 \psi_j / \partial x_0^2 \f$ \n
/// d2psidx(i,1) = \f$ \partial^2 \psi_j / \partial x_1^2 \f$ \n
/// d2psidx(i,2) = \f$ \partial^2 \psi_j / \partial x_0 \partial x_1 \f$ \n
/// \n \b 3D: \n
/// d2psidx(i,0) = \f$ \partial^2 \psi_j / \partial x_0^2 \f$ \n
/// d2psidx(i,1) = \f$ \partial^2 \psi_j / \partial x_1^2 \f$ \n
/// d2psidx(i,2) = \f$ \partial^2 \psi_j / \partial x_2^2 \f$ \n
/// d2psidx(i,3) = \f$ \partial^2 \psi_j / \partial x_0 \partial x_1 \f$ \n
/// d2psidx(i,4) = \f$ \partial^2 \psi_j / \partial x_0 \partial x_2 \f$ \n
/// d2psidx(i,5) = \f$ \partial^2 \psi_j / \partial x_1 \partial x_2 \f$ \n
//===========================================================================
 double FiniteElement::d2shape_eulerian(const Vector<double> &s, 
                                        Shape &psi, 
                                        DShape &dpsi, 
                                        DShape &d2psi) const
 {
  //Find the values of the indices of the shape functions
  //Locally cached.
  //Find the element dimension
  const unsigned el_dim = dim();
  //Find the number of second derivatives required
  const unsigned n_deriv = N2deriv[el_dim];
 
  //Get the values of the shape function and local derivatives
  d2shape_local(s,psi,dpsi,d2psi);
 
  //Allocate memory for the jacobian and inverse jacobian 
  DenseMatrix<double> jacobian(el_dim), inverse_jacobian(el_dim);
  //Calculate the jacobian and inverse jacobian
  const double det = 
   local_to_eulerian_mapping(dpsi,jacobian,inverse_jacobian);

  //Allocate memory for the jacobian of second derivatives
  DenseMatrix<double> jacobian2(n_deriv,el_dim);
  //Assemble the jacobian of second derivatives
  assemble_local_to_eulerian_jacobian2(d2psi,jacobian2);

  //Now set the value of the derivatives
  transform_second_derivatives(jacobian,inverse_jacobian,
                               jacobian2,dpsi,d2psi);
  //Return the determinant of the mapping
  return det;
 }

//===========================================================================
/// \short Compute the geometric shape functions and also first
/// and second derivatives w.r.t. global coordinates at ipt-th integration
/// point
/// Returns Jacobian of mapping from global to local coordinates.
/// This is the most general version, may be overloaded, if desired.
/// \n\n Numbering:
/// \n \b 1D: \n
/// d2psidx(i,0) = \f$ d^2 \psi_j / d x^2 \f$
/// \n \b 2D: \n
/// d2psidx(i,0) = \f$ \partial^2 \psi_j / \partial x_0^2 \f$ \n
/// d2psidx(i,1) = \f$ \partial^2 \psi_j / \partial x_1^2 \f$ \n
/// d2psidx(i,2) = \f$ \partial^2 \psi_j / \partial x_0 \partial x_1 \f$ \n
/// \n \b 3D: \n
/// d2psidx(i,0) = \f$ \partial^2 \psi_j / \partial x_0^2 \f$ \n
/// d2psidx(i,1) = \f$ \partial^2 \psi_j / \partial x_1^2 \f$ \n
/// d2psidx(i,2) = \f$ \partial^2 \psi_j / \partial x_2^2 \f$ \n
/// d2psidx(i,3) = \f$ \partial^2 \psi_j / \partial x_0 \partial x_1 \f$ \n
/// d2psidx(i,4) = \f$ \partial^2 \psi_j / \partial x_0 \partial x_2 \f$ \n
/// d2psidx(i,5) = \f$ \partial^2 \psi_j / \partial x_1 \partial x_2 \f$ \n
//==========================================================================
 double FiniteElement::d2shape_eulerian_at_knot(const unsigned &ipt, 
                                                Shape &psi, 
                                                DShape &dpsi, 
                                                DShape &d2psi) const
 {
  //Find the values of the indices of the shape functions
  //Locally cached
  //Find the element dimension
  const unsigned el_dim = dim();
  //Find the number of second derivatives required
  const unsigned n_deriv = N2deriv[el_dim];

  //Get the values of the shape function and local derivatives
  d2shape_local_at_knot(ipt,psi,dpsi,d2psi);
 
  //Allocate memory for the jacobian and inverse jacobian 
  DenseMatrix<double> jacobian(el_dim), inverse_jacobian(el_dim);
  //Calculate the jacobian and inverse jacobian
  const double det = 
   local_to_eulerian_mapping(dpsi,jacobian,inverse_jacobian);

  //Allocate memory for the jacobian of second derivatives
  DenseMatrix<double> jacobian2(n_deriv,el_dim);
  //Assemble the jacobian of second derivatives
  assemble_local_to_eulerian_jacobian2(d2psi,jacobian2);

  //Now set the value of the derivatives
  transform_second_derivatives(jacobian,inverse_jacobian,
                               jacobian2,dpsi,d2psi);
  //Return the determinant of the mapping
  return det;
 }



//==========================================================================
/// This function loops over the nodal data of the element, adds the
/// GLOBAL equation numbers to the local-to-global look-up scheme and 
/// fills in the Nodal_local_eqn look-up scheme for the local equation 
/// numbers
//==========================================================================
 void FiniteElement::assign_nodal_local_eqn_numbers()
 {
  //Find the number of nodes
  const unsigned n_node = nnode();
  //If there are nodes
  if(n_node > 0)
   {
    //Find the number of local equations assigned so far
    unsigned local_eqn_number = ndof();
   
    //We need to find the total number of values stored at the node
    //Initialise to the number of values stored at the first node
    unsigned n_total_values = node_pt(0)->nvalue();
    //Loop over the other nodes and add the values stored
    for(unsigned n=1;n<n_node;n++)
     {n_total_values += node_pt(n)->nvalue();}

    //If allocated delete the old storage
    if(Nodal_local_eqn)
     {
      delete[] Nodal_local_eqn[0];
      delete[] Nodal_local_eqn;
     }
    
    //If there are no values, we are done, null out the storage and
    //return
    if(n_total_values==0) {Nodal_local_eqn=0; return;}

    //Resize the storage for the nodal local equation numbers
    //Firstly allocate pointers to rows for each node
    Nodal_local_eqn = new int*[n_node];
    //Now allocate storage for the equation numbers
    Nodal_local_eqn[0] = new int[n_total_values];
    //initially all local equations are unclassified
    for(unsigned i=0;i<n_total_values;i++)
     {Nodal_local_eqn[0][i] = Data::Is_unclassified;}
    
    //Loop over the remaining rows and set their pointers
    for(unsigned n=1;n<n_node;++n)
     {
      //Initially set the pointer to the i-th row to the pointer
      //to the i-1th row
      Nodal_local_eqn[n] = Nodal_local_eqn[n-1];
      //Now increase the row pointer by the number of values 
      //stored at the i-1th node
      Nodal_local_eqn[n] += Node_pt[n-1]->nvalue();
     }
    

    //A local queue to store the global equation numbers
    std::deque<unsigned long> global_eqn_number_queue;

    //Now loop over the nodes again and assign local equation numbers
    for(unsigned n=0;n<n_node;n++)
     {
      //Find the number of values stored at the node
      unsigned n_value = node_pt(n)->nvalue();
     
      //Loop over the number of values
      for(unsigned j=0;j<n_value;j++)
       {
        //Get the GLOBAL equation number
        long eqn_number = node_pt(n)->eqn_number(j);
        //If the GLOBAL equation number is positive (a free variable)
        if(eqn_number >= 0)
         {
          //Add the GLOBAL equation number to the queue
          global_eqn_number_queue.push_back(eqn_number);
          //Add the local equation number to the local scheme
          Nodal_local_eqn[n][j] = local_eqn_number;
          //Increase the local number
          local_eqn_number++;
         }
        else
         {
          //Set the local scheme to be pinned
          Nodal_local_eqn[n][j] = Data::Is_pinned;
         }
       }
     }
    
    //Now add our global equations numbers to the internal element storage
    add_global_eqn_numbers(global_eqn_number_queue);
   }
 }


//============================================================================
/// This function calculates the entries of Jacobian matrix, used in 
/// the Newton method, associated with the nodal degrees of freedom.
/// It does this using finite differences, 
/// rather than an analytical formulation, so can be done in total generality.
//==========================================================================
 void FiniteElement::
 fill_in_jacobian_from_nodal_by_fd(Vector<double> &residuals,
                                   DenseMatrix<double> &jacobian)
 {
  //Find the number of nodes
  const unsigned n_node = nnode();
  //If there aren't any nodes, then return straight awayy
  if(n_node == 0) {return;}

  //Call the update function to ensure that the element is in
  //a consistent state before finite differencing starts
  update_before_nodal_fd();

  //Find the number of dofs in the element
  const unsigned n_dof = ndof();
  //Create newres vector
  Vector<double> newres(n_dof);

  //Integer storage for local unknown
  int local_unknown=0;
  
  //Use the default finite difference step
  const double fd_step = Default_fd_jacobian_step;

  //Loop over the nodes
  for(unsigned n=0;n<n_node;n++)
   {
    //Get the number of values stored at the node
    const unsigned n_value = node_pt(n)->nvalue();
   
    //Loop over the number of values
    for(unsigned i=0;i<n_value;i++)
     {
      //Get the local equation number
      local_unknown = nodal_local_eqn(n,i);
      //If it's not pinned
      if(local_unknown >= 0)
       {
        //Store a pointer to the nodal data value
        double* const value_pt = node_pt(n)->value_pt(i);
        
        //Save the old value of the Nodal data
        const double old_var = *value_pt;
       
        //Increment the value of the Nodal data
        *value_pt += fd_step;

        //Now update any slaved variables
        update_in_nodal_fd(i);
        
        //Calculate the new residuals
        get_residuals(newres);
       
        //Do finite differences
        for(unsigned m=0;m<n_dof;m++)
         {
          double sum = (newres[m] - residuals[m])/fd_step;
          //Stick the entry into the Jacobian matrix
          jacobian(m,local_unknown) = sum;
         }
       
        //Reset the Nodal data
        *value_pt = old_var;

        //Reset any slaved variables
        reset_in_nodal_fd(i);
       }
     }
   }
  
  //End of finite difference loop
  //Final reset of any slaved data
  reset_after_nodal_fd();
 }




//=======================================================================
/// Compute derivatives of elemental residual vector with respect
/// to nodal coordinates. Default implementation by FD can be overwritten
/// for specific elements. 
/// dresidual_dnodal_coordinates(l,i,j) = d res(l) / dX_{ij}
////=======================================================================
void FiniteElement::get_dresidual_dnodal_coordinates(
 RankThreeTensor<double>& dresidual_dnodal_coordinates)
{
 // Number of nodes
 unsigned n_nod=nnode();
 
 // If the element has no nodes (why??!!) return straightaway
 if (n_nod==0) return;
 
 // Get dimension from first node
 unsigned dim_nod=node_pt(0)->ndim();
 
 // Number of dofs
 unsigned n_dof=ndof();
 
 // Get reference residual
 Vector<double> res(n_dof);
 Vector<double> res_pls(n_dof);
 get_residuals(res);
  
  // FD step 
  double eps_fd=GeneralisedElement::Default_fd_jacobian_step;

  // Do FD loop
  for (unsigned j=0;j<n_nod;j++)
   {
    // Get node
    Node* nod_pt=node_pt(j);

    // Loop over coordinate directions
    for (unsigned i=0;i<dim_nod;i++)
     {
      // Make backup
      double backup=nod_pt->x(i);

      // Do FD step. No node update required as we're
      // attacking the coordinate directly...
      nod_pt->x(i)+=eps_fd;

      // Perform auxiliary node update function
      nod_pt->perform_auxiliary_node_update_fct();
  
      // Get advanced residual
      get_residuals(res_pls);
      
      // Fill in FD entries [Loop order is "wrong" here as l is the
      // slow index but this is in a function that's costly anyway
      // and gives us the fastest loop outside where these tensor
      // is actually used.]
      for (unsigned l=0;l<n_dof;l++)
       {
        dresidual_dnodal_coordinates(l,i,j)=(res_pls[l]-res[l])/eps_fd;
       }

      // Reset coordinate. No node update required as we're
      // attacking the coordinate directly...
      nod_pt->x(i)=backup;

      // Perform auxiliary node update function
      nod_pt->perform_auxiliary_node_update_fct();
      
     }
   }
  
 }


//===============================================================
/// Return the number of the node located at *node_pt 
/// if this node is in the element, else return -1; 
//===============================================================
 int FiniteElement::get_node_number(Node* const &global_node_pt)
 {
   //Initialise the number to -1
   int number = -1;
   //Find the number of nodes
   unsigned n_node = nnode();
#ifdef PARANOID
   {
    //Error check that node does not appear in element more than once
    unsigned count=0;
    //Storage for the local node numbers of the element
    std::vector<int> local_node_number;
    //Loop over the nodes
    for(unsigned i=0;i<n_node;i++)
     {
      //If the node is present increase the counter
      //and store the local node number
      if(node_pt(i)==global_node_pt) 
       {
        ++count;
        local_node_number.push_back(i);
       }
     }
    
    //If the node appears more than once, complain
    if(count > 1)
     {
      std::ostringstream error_stream;
      error_stream << "Node " << global_node_pt << " appears " 
                   << count << " times in an element." << std::endl
                   << "In positions: ";
      for(std::vector<int>::iterator it=local_node_number.begin();
          it!=local_node_number.end();++it)
       {
        error_stream << *it << " ";
       }
      error_stream << std::endl
                   << "That seems very odd." << std::endl;
      
      throw OomphLibError(error_stream.str(),
                          "FiniteElement::get_node_number()",
                          OOMPH_EXCEPTION_LOCATION);
     }
   }
#endif
   
   //Loop over the nodes
   for(unsigned i=0;i<n_node;i++)
    {
     //If the passed node pointer is present in the element
     //set number to be its local node number
     if(node_pt(i)==global_node_pt)
      {      
       number=i;
       break;
      }
    }
   
   //Return the node number
   return number;
 }
 

 //==========================================================================
 /// \short If there is a node at the local coordinate, s, return the pointer 
 /// to  the node. If not return 0. Note that this is a default, brute
 /// force implementation, can almost certainly be made more efficient for
 /// specific elements.
 //==========================================================================
 Node* FiniteElement::get_node_at_local_coordinate(const Vector<double> &s)
  {
   //Locally cache the tolerance
   const double tol = Node_location_tolerance;
   Vector<double> s_node;
   //Locally cache the member data
   const unsigned el_dim = Elemental_dimension;
   const unsigned n_node = Nnode;
   //Loop over the nodes
   for(unsigned n=0;n<n_node;n++)
    {
     bool Match = true;
     //Find the local coordinate of the node
     local_coordinate_of_node(n,s_node);
     for(unsigned i=0;i<el_dim;i++)
      {
       //Calculate the difference between coordinates
       //and if it's bigger than our tolerance 
       //break out of the (inner)loop
       if(std::fabs(s[i] - s_node[i]) > tol)
        {
         Match = false;
         break;
        }
      }
     //If we haven't complained then we have a match
     if(Match) {return node_pt(n);}
    }
   //If we get here, we have no match
   return 0;
  }
   


//======================================================================
/// Return FE interpolated coordinate x[i] at local coordinate s
//======================================================================
 double FiniteElement::interpolated_x(const Vector<double> &s, 
                                      const unsigned &i) const
 {
  //Find the number of nodes
  const unsigned n_node = nnode();
  //Find the number of positional types
  const unsigned n_position_type = nnodal_position_type();
  //Assign storage for the local shape function
  Shape psi(n_node,n_position_type);
  //Find the values of shape function
  shape(s,psi);

  //Initialise value of x
  double interpolated_x = 0.0;
  //Loop over the local nodes
  for(unsigned l=0;l<n_node;l++) 
   {
    //Loop over the number of dofs
    for(unsigned k=0;k<n_position_type;k++)
     {
      interpolated_x += nodal_position_gen(l,k,i)*psi(l,k);
     }
   }
   
  return(interpolated_x);
 }

//=========================================================================
/// Return FE interpolated coordinate x[i] at local coordinate s
/// at previous timestep t (t=0: present; t>0: previous timestep)
//========================================================================
 double FiniteElement::interpolated_x(const unsigned &t,
                                      const Vector<double> &s, 
                                      const unsigned &i) const
 {
  //Find the number of nodes
  const unsigned n_node = nnode();
  //Find the number of positional types
  const unsigned n_position_type = nnodal_position_type();

  //Assign storage for the local shape function
  Shape psi(n_node,n_position_type);
  //Find the values of shape function
  shape(s,psi);

  //Initialise value of x
  double interpolated_x = 0.0;
  //Loop over the local nodes
  for(unsigned l=0;l<n_node;l++) 
   {
    //Loop over the number of dofs
    for(unsigned k=0;k<n_position_type;k++)
     {
      interpolated_x += nodal_position_gen(t,l,k,i)*psi(l,k);
     }
   }
   
  return(interpolated_x);
 }

//=======================================================================
/// Return FE interpolated position x[] at local coordinate s as Vector
//=======================================================================
 void FiniteElement::interpolated_x(const Vector<double> &s, Vector<double> &x)
  const
 {
   //Find the number of nodes
  const unsigned n_node = nnode();
  //Find the number of positional types
  const unsigned n_position_type = nnodal_position_type();
  //Find the dimension stored in the node
  const unsigned nodal_dim = nodal_dimension();

  //Assign storage for the local shape function
  Shape psi(n_node,n_position_type);
  //Find the values of shape function
  shape(s,psi);

  //Loop over the dimensions
  for(unsigned i=0;i<nodal_dim;i++)
   {
    //Initilialise value of x[i] to zero
    x[i] = 0.0;
    //Loop over the local nodes
    for(unsigned l=0;l<n_node;l++) 
     {
      //Loop over the number of dofs
      for(unsigned k=0;k<n_position_type;k++)
       {
        x[i] += nodal_position_gen(l,k,i)*psi(l,k);
       }
     }
   }
 }

//==========================================================================
/// Return FE interpolated position x[] at local coordinate s
/// at previous timestep t as Vector (t=0: present; t>0: previous timestep)
//==========================================================================
 void FiniteElement::interpolated_x(const unsigned &t, 
                                    const Vector<double> &s, 
                                    Vector<double>& x) const
 {
  //Find the number of nodes
  const unsigned n_node = nnode();
  //Find the number of positional types
  const unsigned n_position_type = nnodal_position_type();
  //Find the dimensions of the nodes
  const unsigned nodal_dim = nodal_dimension();

  //Assign storage for the local shape function
  Shape psi(n_node,n_position_type);
  //Find the values of shape function
  shape(s,psi);

  //Loop over the dimensions
  for(unsigned i=0;i<nodal_dim;i++)
   {
    //Initilialise value of x[i] to zero
    x[i] = 0.0;
    //Loop over the local nodes
    for(unsigned l=0;l<n_node;l++) 
     {
      //Loop over the number of dofs
      for(unsigned k=0;k<n_position_type;k++)
       {
        x[i] += nodal_position_gen(t,l,k,i)*psi(l,k);
       }
     }
   }
 }

//========================================================================
/// \short Calculate the Jacobian of the mapping between local and global
/// coordinates at the position s
//========================================================================
 double FiniteElement::J_eulerian(const Vector<double> &s) const
 {
  //Find the number of nodes and position types
  const unsigned n_node = nnode();
  const unsigned n_position_type = nnodal_position_type();
  //Find the dimension of the node and element
  const unsigned n_dim_node = nodal_dimension();
  const unsigned n_dim_element = dim();

  //Set up dummy memory for the shape functions
  Shape psi(n_node,n_position_type);
  DShape dpsids(n_node,n_position_type,n_dim_element);
  //Get the shape functions and local derivatives
  dshape_local(s,psi,dpsids);
  
  //Right calculate the base vectors
  DenseMatrix<double> interpolated_G(n_dim_element,n_dim_node);
  assemble_eulerian_base_vectors(dpsids,interpolated_G);
 
  //Calculate the metric tensor of the element
  DenseMatrix<double> G(n_dim_element,n_dim_element,0.0);
  for(unsigned i=0;i<n_dim_element;i++)
   {
    for(unsigned j=0;j<n_dim_element;j++)
     {
      for(unsigned k=0;k<n_dim_node;k++) 
       {
        G(i,j) += interpolated_G(i,k)*interpolated_G(j,k);
       }
     }
   }

  //Calculate the determinant of the metric tensor
  double det = 0.0;
  switch(n_dim_element)
   {
   case 0:
    throw OomphLibError("Cannot calculate J_eulerian() for point element\n",
                        "FiniteElement::J_eulerian()",
                        OOMPH_EXCEPTION_LOCATION);
    break;
   case 1:
    det = G(0,0);
    break;
   case 2:
    det = G(0,0)*G(1,1) - G(0,1)*G(1,0);
    break;
   case 3:
    det = G(0,0)*G(1,1)*G(2,2) + G(0,1)*G(1,2)*G(2,0) + G(0,2)*G(1,0)*G(2,1)
     - G(0,0)*G(1,2)*G(2,1) - G(0,1)*G(1,0)*G(2,2) - G(0,2)*G(1,1)*G(2,0);
    break;
   default:
    oomph_info << "More than 3 dimensions in J_eulerian()" << std::endl;
    break;
   }

#ifdef PARANOID
  check_jacobian(det);
#endif

  //Return the Jacobian (square-root of the determinant of the metric tensor)
  return sqrt(det);
 }

//========================================================================
/// \short Compute the Jacobian of the mapping between the local and global
/// coordinates at the ipt-th integration point
//========================================================================
 double FiniteElement::J_eulerian_at_knot(const unsigned &ipt)
  const
 {
  //Find the number of nodes
  const unsigned n_node = nnode();
  //Find the number of position types
  const unsigned n_position_type = nnodal_position_type();
  //Find the dimension of the node and element
  const unsigned n_dim_node = nodal_dimension();
  const unsigned n_dim_element = dim();

  //Set up dummy memory for the shape functions
  Shape psi(n_node,n_position_type);
  DShape dpsids(n_node,n_position_type,n_dim_element);
  //Get the shape functions and local derivatives at the knot
  //This call may use the stored versions, which is why this entire
  //function doesn't just call J_eulerian(s), after reading out s from
  //the knots.
  dshape_local_at_knot(ipt,psi,dpsids);

  //Right calculate the base vectors
  DenseMatrix<double> interpolated_G(n_dim_element,n_dim_node);
  assemble_eulerian_base_vectors(dpsids,interpolated_G);
 
  //Calculate the metric tensor of the element
  DenseMatrix<double> G(n_dim_element,n_dim_element,0.0);
  for(unsigned i=0;i<n_dim_element;i++)
   {
    for(unsigned j=0;j<n_dim_element;j++)
     {
      for(unsigned k=0;k<n_dim_node;k++) 
       {
        G(i,j) += interpolated_G(i,k)*interpolated_G(j,k);
       }
     }
   }

  //Calculate the determinant of the metric tensor
  double det = 0.0;
  switch(n_dim_element)
   {
   case 0:
    throw OomphLibError("Cannot calculate J_eulerian() for point element\n",
                        "FiniteElement::J_eulerian_at_knot()",
                        OOMPH_EXCEPTION_LOCATION);
    break;
   case 1:
    det = G(0,0);
    break;
   case 2:
    det = G(0,0)*G(1,1) - G(0,1)*G(1,0);
    break;
   case 3:
    det = G(0,0)*G(1,1)*G(2,2) + G(0,1)*G(1,2)*G(2,0) + G(0,2)*G(1,0)*G(2,1)
     - G(0,0)*G(1,2)*G(2,1) - G(0,1)*G(1,0)*G(2,2) - G(0,2)*G(1,1)*G(2,0);
    break;
   default:
    oomph_info << "More than 3 dimensions in J_eulerian()" << std::endl;
    break;
   }
 
  //Return the Jacobian (square-root of the determinant of the metric tensor)
  return sqrt(det);
 }

//====================================================================
/// Calculate the size of the element.
//====================================================================
 double FiniteElement::size() const
 {
  //Initialise the sum to zero
  double sum = 0.0;
 
  //Loop over the integration points
  const unsigned n_intpt = integral_pt()->nweight();
  for(unsigned ipt=0;ipt<n_intpt;ipt++)
   {
    //Get the integral weight
    double w = integral_pt()->weight(ipt);
    //Get the value of the Jacobian of the mapping to global coordinates
    double J = J_eulerian_at_knot(ipt);

    //Add the product to the sum
    sum += w*J;
   }
 
  //Return the answer
  return(sum);
 }

//==================================================================
/// Integrate Vector-valued time-dep function over element
//==================================================================
 void FiniteElement::
 integrate_fct(FiniteElement::UnsteadyExactSolutionFctPt integrand_fct_pt,
               Vector<double>& integral)
 {
  //Initialise all components of integral Vector and setup integrand vector
  const unsigned ncomponents=integral.size();
  Vector<double> integrand(ncomponents);
  for (unsigned i=0;i<ncomponents;i++) {integral[i]=0.0;}

  // Figure out the global (Eulerian) spatial dimension of the
  // element 
  const unsigned n_dim_eulerian = nodal_dimension();
  
  // Allocate Vector of global Eulerian coordinates
  Vector<double> x(n_dim_eulerian);

  // Set the value of n_intpt
  const unsigned n_intpt = integral_pt()->nweight();

  // Vector of local coordinates
  const unsigned n_dim = dim();
  Vector<double> s(n_dim);

  //Loop over the integration points
  for(unsigned ipt=0;ipt<n_intpt;ipt++)
   {
    //Assign the values of s 
    for(unsigned i=0;i<n_dim;i++) {s[i] = integral_pt()->knot(ipt,i);}

    //Assign the values of the global Eulerian coordinates
    for(unsigned i=0;i<n_dim_eulerian;i++) {x[i] = interpolated_x(s,i);}

    //Get the integral weight
    double w = integral_pt()->weight(ipt);
     
    //Get Jacobian of mapping
    double J = J_eulerian(s);
     
    // Evaluate the integrand at the knot points
    integrand_fct_pt(time(),x,integrand);

    //Add to components of integral Vector
    for (unsigned i=0;i<ncomponents;i++) {integral[i]+=integrand[i]*w*J;}
   }
 }

//==================================================================
/// Integrate Vector-valued function over element
//==================================================================
 void FiniteElement::
 integrate_fct(FiniteElement::SteadyExactSolutionFctPt integrand_fct_pt,
               Vector<double>& integral)
 {
  //Initialise all components of integral Vector
  const unsigned ncomponents=integral.size();
  Vector<double> integrand(ncomponents);
  for (unsigned i=0;i<ncomponents;i++) {integral[i]=0.0;}

  // Figure out the global (Eulerian) spatial dimension of the
  // element by checking the Eulerian dimension of the nodes
  const unsigned n_dim_eulerian = nodal_dimension();
  
  // Allocate Vector of global Eulerian coordinates
  Vector<double> x(n_dim_eulerian);

  // Set the value of n_intpt
  const unsigned n_intpt = integral_pt()->nweight();

  // Vector of local coordinates
  const unsigned n_dim = dim();
  Vector<double> s(n_dim);

  //Loop over the integration points
  for(unsigned ipt=0;ipt<n_intpt;ipt++)
   {
    //Assign the values of s 
    for(unsigned i=0;i<n_dim;i++) {s[i] = integral_pt()->knot(ipt,i);}

    //Assign the values of the global Eulerian coordinates
    for(unsigned i=0;i<n_dim_eulerian;i++) {x[i] = interpolated_x(s,i);}

    //Get the integral weight
    double w = integral_pt()->weight(ipt);
     
    //Get Jacobian of mapping
    double J = J_eulerian(s);
     
    // Evaluate the integrand at the knot points
    integrand_fct_pt(x,integrand);

    //Add to components of integral Vector
    for (unsigned i=0;i<ncomponents;i++) {integral[i]+=integrand[i]*w*J;}
   }
 }

//==========================================================================
/// Self-test: Have all internal values been classified as 
/// pinned/unpinned? Has pointer to spatial integration scheme
/// been set? Return 0 if OK. 
//==========================================================================
 unsigned FiniteElement::self_test()
 {
  // Initialise
  bool passed=true;

  if(GeneralisedElement::self_test()!=0) {passed=false;}

  // Check that pointer to spatial integration scheme has been set
  if(integral_pt()==0)
   {
    passed=false;

    OomphLibWarning(
     "Pointer to spatial integration scheme has not been set.",
     "FiniteElement::self_test()",
     OOMPH_EXCEPTION_LOCATION);
   }

  //If the dimension of the element is zero (point element), there
  //is not jacobian
  const unsigned dim_el = dim();

  if(dim_el> 0)
   {
    // Loop over integration points to check sign of Jacobian
    //-------------------------------------------------------
   
    //Set the value of n_intpt
    const unsigned n_intpt = integral_pt()->nweight();
 
    //Set the Vector to hold local coordinates
    Vector<double> s(dim_el);

    //Find the number of local nodes
    const unsigned n_node = nnode();
    const unsigned n_position_type = nnodal_position_type();
    //Set up memory for the shape and test functions
    Shape psi(n_node,n_position_type);
    DShape dpsidx(n_node,dim_el);
            
    // Jacobian
    double jacobian;
   

    // Two ways of testing for negative Jacobian for non-FaceElements
    unsigned ntest=1; 
   
    // For FaceElements checking the Jacobian via dpsidx doesn't
    // make sense
    FiniteElement* tmp_pt=const_cast<FiniteElement*>(this);
    FaceElement* face_el_pt=dynamic_cast<FaceElement*>(tmp_pt);
    //oomph_info << "ntest face_el_pt: " << ntest << " " << face_el_pt << std::endl;
    if (face_el_pt==0)
     {
      ntest=2;
      //oomph_info << "Changed to ntest=" << ntest << std::endl;
     }
   
    // For now overwrite -- the stuff above fails for Bretherton.
    // Not sure why.
    ntest=1;

    //Loop over the integration points
    for(unsigned ipt=0;ipt<n_intpt;ipt++)
     {
      //Assign values of s
      for(unsigned i=0;i<dim_el;i++)
       {
        s[i] = integral_pt()->knot(ipt,i);
       }
     

      // Do tests
      for (unsigned test=0;test<ntest;test++)
       {
   
        switch (test)
         {
         
         case 0:
         
          // Get the jacobian from the mapping between local and Eulerian
          // coordinates
          jacobian = J_eulerian(s);
         
          break;
         
         case 1:

          //Call the geometrical shape functions and derivatives  
          //This also computes the Jacobian by a slightly different 
          //method
          jacobian = dshape_eulerian_at_knot(ipt,psi,dpsidx);

          break;

         default:

          throw OomphLibError("Never get here",
                              "FiniteElement::self_test()",
                              OOMPH_EXCEPTION_LOCATION);
         }

       
        //Check for a singular jacobian
        if(std::fabs(jacobian) < 1.0e-16)
         {
          std::ostringstream warning_stream;
          warning_stream << "Determinant of Jacobian matrix is zero at ipt "
                         << ipt << std::endl;
          OomphLibWarning(warning_stream.str(),
                          "FiniteElement::self_test()",
                          OOMPH_EXCEPTION_LOCATION);
          passed = false;
          //Skip the next test
          continue;
         }
       
        //Check sign of Jacobian
        if ((Accept_negative_jacobian==false) && (jacobian < 0.0))
         {
          std::ostringstream warning_stream;
          warning_stream << "Jacobian negative at integration point ipt=" 
                         << ipt << std::endl;
          warning_stream 
           << "If you think that this is what you want you may: " << std::endl
           << "set the (static) flag " 
           << "FiniteElement::Accept_negative_jacobian to be true" << std::endl;
         
          OomphLibWarning(warning_stream.str(),
                          "FiniteElement::self_test()",
                          OOMPH_EXCEPTION_LOCATION);
          passed=false;
         }
       
       } // end of loop over two tests

     }
   } //End of non-zero dimension check

   
  // Return verdict
  if (passed) {return 0;}
  else {return 1;}
 }




//=======================================================================
/// Return the t-th time-derivative of the 
/// i-th FE-interpolated Eulerian coordinate at 
/// local coordinate s.
//=======================================================================
 double FiniteElement::interpolated_dxdt(const Vector<double> &s, 
                                         const unsigned &i,
                                         const unsigned &t_deriv)
 {
  //Find the number of nodes and positions (locally cached)
  const unsigned n_node = nnode();
  const unsigned n_position_type = nnodal_position_type();
  // Get shape functions: specify # of nodes, # of positional dofs
  Shape psi(n_node,n_position_type);
  shape(s,psi);
 
  // Initialise 
  double drdt=0.0;
 
  // Assemble time derivative
  //Loop over nodes
  for(unsigned l=0;l<n_node;l++) 
   {
    //Loop over types of dof
    for(unsigned k=0;k<n_position_type;k++)
     {
      drdt+=dnodal_position_gen_dt(t_deriv,l,k,i)*psi(l,k);  
     }
   }
  return drdt;
 }



//=======================================================================
/// Compute t-th time-derivative of the
/// FE-interpolated Eulerian coordinate vector at
/// local coordinate s. 
//=======================================================================
 void FiniteElement::interpolated_dxdt(const Vector<double> &s, 
                                       const unsigned &t_deriv,
                                       Vector<double>& dxdt)
 {
  //Find the number of nodes and positions (locally cached)
  const unsigned n_node = nnode();
  const unsigned n_position_type = nnodal_position_type();
  const unsigned nodal_dim = nodal_dimension();

  // Get shape functions: specify # of nodes, # of positional dofs
  Shape psi(n_node,n_position_type);
  shape(s,psi);
 
  // Loop over directions
  for (unsigned i=0;i<nodal_dim;i++)
   {
    // Initialise 
    dxdt[i]=0.0;
    
    // Assemble time derivative
    //Loop over nodes
    for(unsigned l=0;l<n_node;l++) 
     {
      //Loop over types of dof
      for(unsigned k=0;k<n_position_type;k++)
       {
        dxdt[i]+=dnodal_position_gen_dt(t_deriv,l,k,i)*psi(l,k);  
       }
     }
   }
 }

//============================================================================
/// Calculate the interpolated value of zeta, the intrinsic coordinate
/// of the element when viewed as a compound geometric object within a Mesh
/// as a function of the local coordinate of the element, s.
///  The default 
/// assumption is the zeta is interpolated using the shape functions of
/// the element with the values given by zeta_nodal().
/// A MacroElement representation of the intrinsic coordinate parametrised
/// by the local coordinate s is used if available. 
/// Choosing the MacroElement representation of zeta (Eulerian x by default)
/// allows a correspondence to be established between elements on different
/// Meshes covering the same curvilinear domain in cases where one element
/// is much coarser than the other.
//==========================================================================
 void FiniteElement::interpolated_zeta(const Vector<double> &s,  
                                       Vector<double> &zeta) const       
{
 //If there is a macro element use it
 if(Macro_elem_pt!=0) {this->get_x_from_macro_element(s,zeta);}
 //Otherwise interpolate zeta_nodal using the shape functions
 else
  {
   //Find the number of nodes
   const unsigned n_node = this->nnode();
   //Find the number of positional types
   const unsigned n_position_type = this->nnodal_position_type();
   //Storage for the shape functions
   Shape psi(n_node,n_position_type);
   //Get the values of the shape functions at the local coordinate s
   this->shape(s,psi);
   
   //Find the number of coordinates
   const unsigned ncoord = this->dim();
   //Initialise the value of zeta to zero
   for(unsigned i=0;i<ncoord;i++) {zeta[i] = 0.0;}

   //Add the contributions from each nodal dof to the interpolated value
   //of zeta.
   for(unsigned l=0;l<n_node;l++)
    {
     for(unsigned k=0;k<n_position_type;k++)
      {
       //Locally cache the value of the shape function
       const double psi_ = psi(l,k);
       for(unsigned i=0;i<ncoord;i++)
        {
         zeta[i] += this->zeta_nodal(l,k,i)*psi_;
        }
      }
    }
  }
}

//==========================================================================
/// For a given value of zeta, the "global" intrinsic coordinate of
/// a mesh of FiniteElements represented as a compound geometric object,
/// find the local coordinate in this element that corresponds to the 
/// requested value of zeta. 
/// This is achieved in generality by using Newton's method to find the value 
/// of the local coordinate, s, such that
/// interpolated_zeta(s) is equal to the requested value of zeta.
/// If zeta cannot be located in this element, geom_object_pt is set
/// to NULL. If zeta is located in this element, we return its "this"
/// pointer.
/// Setting the optional bool argument to true means that the coordinate
/// argument "s" is used as the initial guess. (Default is false).
//=========================================================================
void FiniteElement::locate_zeta(const Vector<double> &zeta,
                                GeomObject*& geom_object_pt, Vector<double> &s,
                                const bool& use_coordinate_as_initial_guess)
    {
     //Find the number of coordinates, the dimension of the element
     //This must be the same for the local and intrinsic global coordinate
     unsigned ncoord = this->dim();

     //Assign storage for the vector and matrix used in Newton's method
     Vector<double> dx(ncoord,0.0);
     DenseDoubleMatrix jacobian(ncoord,ncoord,0.0);

     // Make a list of (equally-spaced) local coordinates inside the element
     unsigned n_list=Locate_zeta_helpers::N_local_points; 
     double list_space=(1.0/(double(n_list)-1.0))*(s_max()-s_min());
     Vector<Vector<double> > s_list;

     // If the boolean argument use_coordinate_as_initial_guess was set
     // to true then we don't need to initialise s
     if (!use_coordinate_as_initial_guess)
      {
       // If there is no macro element
       if(Macro_elem_pt==0)
        {
         // Default to the centre of the element
         for(unsigned i=0;i<ncoord;i++)
          {
           s[i] = 0.5*(s_max()+s_min());
          }
        }
       else
        {
         // Create a list of coordinates within the element
         if (ncoord==1)
          {
           for (unsigned i=0;i<n_list;i++)
            {
             Vector<double> s_c(ncoord);
             s_c[0]=s_min()+(double(i)*list_space);
             s_list.push_back(s_c);
            }
          }
         else if (ncoord==2)
          {
           for (unsigned i=0;i<n_list;i++)
            {
             for (unsigned j=0;j<n_list;j++)
              {
               Vector<double> s_c(ncoord);
               s_c[0]=s_min()+(double(i)*list_space);
               s_c[1]=s_min()+(double(j)*list_space);
               s_list.push_back(s_c);
              }
            }
          }
         else if (ncoord==3)
          {
           for (unsigned i=0;i<n_list;i++)
            {
             for (unsigned j=0;j<n_list;j++)
              {
               for (unsigned k=0;k<n_list;k++)
                {
                 Vector<double> s_c(ncoord);
                 s_c[0]=s_min()+(double(i)*list_space);
                 s_c[1]=s_min()+(double(j)*list_space);
                 s_c[2]=s_min()+(double(k)*list_space);
                 s_list.push_back(s_c);
                }
              }
            }
          }
         else
          {
           // Shouldn't get in here...
           std::ostringstream error_stream;
           error_stream << "Element dimension is not equal to 1, 2 or 3 - ?\n";
           throw OomphLibError(error_stream.str(),
                               "FiniteElement::locate_zeta()",
                               OOMPH_EXCEPTION_LOCATION);
          }
        }
      }

     //Counter for the number of Newton steps
     unsigned count=0;

     //Control flag for the Newton loop
     bool keep_going=true;
   
     //Storage for the interpolated value of x
     Vector<double> inter_x(ncoord);
   
     // If no macro element, or we have specified the coordinate already
     if ((macro_elem_pt()==0) || (use_coordinate_as_initial_guess))
      {
       //Get the value of x at the initial guess
       this->interpolated_zeta(s,inter_x);
   
       //Set up the residuals
       for(unsigned i=0;i<ncoord;i++) {dx[i] = zeta[i] - inter_x[i];}
      }
     else
      {
       // Find the smallest residual from the list of coordinates made earlier
       double my_min_resid=DBL_MAX;
       Vector<double> s_local(ncoord);
       Vector<double> work_x(ncoord);
       Vector<double> work_dx(ncoord);

       unsigned n_list_coord=s_list.size();

       for (unsigned i_coord=0; i_coord<n_list_coord; i_coord++)
        {
         for (unsigned i=0;i<ncoord;i++)
          {
           s_local[i]=s_list[i_coord][i];
          }
         // get_x for this coordinate
         this->interpolated_zeta(s_local,work_x);

         // calculate residuals
         for(unsigned i=0;i<ncoord;i++)
          {
           work_dx[i] = zeta[i] - work_x[i];
          }

         double maxres =
          std::fabs(*std::max_element(work_dx.begin(),work_dx.end(),
                                     AbsCmp<double>()));

         // test against previous residuals
         if (maxres<my_min_resid)
          {
           my_min_resid=maxres;
           dx=work_dx;
           inter_x=work_x;
           s=s_local;
          }

        }

      }

     //Main Newton Loop
     do    // start of do while loop
      {
       //Increase loop counter
       count++;

       //Bail out if necessary (without an error for now...)
       if(count > Locate_zeta_helpers::Max_newton_iterations)
        {
         keep_going=false;
         continue;
        }
	     
       //If it's the first time round the loop, check the initial residuals
       if(count==1)
        {
         double maxres =
          std::fabs(*std::max_element(dx.begin(),dx.end(),AbsCmp<double>()));

         //If it's small enough exit
         if(maxres < Locate_zeta_helpers::Newton_tolerance)
          {
           keep_going=false;
           continue;
          }
        }

       //Is there a macro element? If so, assemble the Jacobian by FD-ing
       if (macro_elem_pt()!=0)
        {
         // Assemble jacobian on the fly by finite differencing
         Vector<double> work_s=s;
         Vector<double> r=inter_x; // i.e. the result of previous call to get_x

         // Finite difference step
         double fd_step=GeneralisedElement::Default_fd_jacobian_step;

         // Storage for calculated r from incremented s
         Vector<double> work_r(ncoord,0.0);

         // Loop over s coordinates
         for (unsigned i=0; i<ncoord; i++)
          {
           // Increment work_s by a small amount
           work_s[i]+=fd_step;

           // Calculate work_r from macro element
           this->interpolated_zeta(work_s,work_r);

           // Loop over r to fill Jacobian
           for (unsigned j=0; j<ncoord; j++)
            {
             jacobian(j,i)=-(work_r[j]-r[j])/fd_step;
            }
         
           // Reset work_s
           work_s[i]=s[i];

          }

        }
       else // no macro element, so compute Jacobian with shape functions etc.
        {
         //Compute the entries of the Jacobian matrix
         unsigned n_node = this->nnode();
         unsigned n_position_type = this->nnodal_position_type();
         Shape psi(n_node,n_position_type);
         DShape dpsids(n_node,n_position_type,ncoord);

         //Get the local shape functions and their derivatives
         dshape_local(s,psi,dpsids);
     
         //Calculate the values of dxds
         DenseMatrix<double> interpolated_dxds(ncoord,ncoord,0.0);

// MH: No longer needed
//          //This implementation will only work for n_position_type=1
//          //since the function nodal_position_gen does not yet exist
// #ifdef PARANOID
//          if (n_position_type!=1)
//           {
//            std::ostringstream error_stream;
//            error_stream << "This implementation does not exist yet;\n"
//                         << "it currently uses raw_nodal_position_gen\n"
//                         << "which does not take hangingness into account\n"
//                         << "It will work if n_position_type=1\n";
//            throw OomphLibError(error_stream.str(),
//                                "FiniteElement::locate_zeta()",
//                                OOMPH_EXCEPTION_LOCATION);
//           }
// #endif

         // Loop over the nodes
         for(unsigned l=0;l<n_node;l++)
          {
           // Loop over position type even though it should be 1 - the
           // functionality for n_position_type>1 will exist in the future
           for(unsigned k=0;k<n_position_type;k++)
            {
             // Add the contribution from the nodal coordinates to the matrix
             for(unsigned i=0;i<ncoord;i++)
              {
               for(unsigned j=0;j<ncoord;j++)
                {
                 interpolated_dxds(i,j) +=
                  this->zeta_nodal(l,k,i)*dpsids(l,k,j);
                }
              }
            }
          }

         //The entries of the Jacobian matrix are merely dresiduals/ds
         //i.e. - dx/ds
         for(unsigned i=0;i<ncoord;i++)
          {
           for(unsigned j=0;j<ncoord;j++)
            {
             jacobian(i,j) = - interpolated_dxds(i,j);
            }
          }

        }

       //Now solve the damn thing
       try
        {
         jacobian.solve(dx);
        }
       catch(OomphLibError &error)
        {
         oomph_info << "Error in linear solve for "
                    << "FiniteElement::locate_zeta"
                    << std::endl;
         oomph_info << "Should not affect the result!" << std::endl;
        }

       //Add the correction to the local coordinates
       for(unsigned i=0;i<ncoord;i++) {s[i] -= dx[i];}
     
       //Get the new residuals
       this->interpolated_zeta(s,inter_x);
       for(unsigned i=0;i<ncoord;i++) {dx[i] = zeta[i] - inter_x[i];}
     
       //Get the maximum residuals
       double maxres =
        std::fabs(*std::max_element(dx.begin(),dx.end(),AbsCmp<double>()));

       //If we have converged jump straight to the test at the end of the loop
       if(maxres < Locate_zeta_helpers::Newton_tolerance)
        {
         keep_going=false;
         continue;
        }
      }
     while(keep_going);

     //Test whether the local coordinate are valid or not
     bool valid=local_coord_is_valid(s,
                                     Locate_zeta_helpers::Rounding_tolerance);
     if (!valid)
      {
       geom_object_pt=0;
       return;
      }

     // It is also possible now that it may not have converged "correctly", 
     // i.e. count is greater than Max_newton_iterations
     if (count > Locate_zeta_helpers::Max_newton_iterations)
      {
       // Don't trust the current answer, return null
       geom_object_pt=0;
       return;
      }

     //Otherwise the required point is located in "this" element:
     geom_object_pt = this;

    }


//=======================================================================
/// Loop over all nodes in the element and update their positions
/// using each node's (algebraic) update function
//=======================================================================
 void FiniteElement::node_update()
 {
  const unsigned n_node = nnode();
  for(unsigned n=0;n<n_node;n++) {node_pt(n)->node_update();}
 }

//======================================================================
/// The purpose of this function is to identify all possible
/// Data that can affect the fields interpolated by the FiniteElement.
/// The information will typically be used in interaction problems in
/// which the FiniteElement provides a forcing term for an 
/// ElementWithExternalElement. The Data must be provided as 
/// \c paired_load data containing
///  - the pointer to a Data object
/// and
/// - the index of the value in that Data object
/// .
/// The generic implementation (should be overloaded in more specific
/// applications) is to include all nodal and internal Data stored in
/// the FiniteElement. Note that the geometric data, 
/// which includes the positions
/// of SolidNodes, is treated separately by the function 
/// \c identify_geometric_data()
//======================================================================
void FiniteElement::identify_field_data_for_interactions(
 std::set<std::pair<Data*,unsigned> > &paired_field_data)
{
 //Loop over all internal data
 const unsigned n_internal = this->ninternal_data();
 for(unsigned n=0;n<n_internal;n++)
  {
   //Cache the data pointer
   Data* const dat_pt = this->internal_data_pt(n);
   //Find the number of data values stored in the data object
   const unsigned n_value = dat_pt->nvalue();
   //Add the index of each data value and the pointer to the set
   //of pairs
   for(unsigned i=0;i<n_value;i++)
    {
     paired_field_data.insert(std::make_pair(dat_pt,i));
    }
  }

 //Loop over all the nodes
 const unsigned n_node = this->nnode();
 for(unsigned n=0;n<n_node;n++)
  {
   //Cache the node pointer
   Node* const nod_pt = this->node_pt(n);
   //Find the number of values stored at the node
   const unsigned n_value = nod_pt->nvalue();
   //Add the index of each data value and the pointer to the set
   //of pairs
   for(unsigned i=0;i<n_value;i++) 
    {
     paired_field_data.insert(std::make_pair(nod_pt,i));
    }
  }
}

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////


//========================================================================
/// \short Calculate the Jacobian of the mapping between local and global
/// coordinates at the position s. Overloaded from FiniteElement.
//========================================================================
 double FaceElement::J_eulerian(const Vector<double> &s) const
 {

  //Find out the sptial dimension of the element
  unsigned n_dim_el = this->dim();

  // Bail out if we're in a point element -- not sure what
  // J_eulerian actually is, but this is harmless 
  if (n_dim_el==0) return 1.0;
   
  //Find out how many nodes there are
  unsigned n_node = nnode();
  
  //Find out how many positional dofs there are
  unsigned n_position_type = this->nnodal_position_type();
  
  //Find out the dimension of the node
  unsigned n_dim = this->nodal_dimension();
    
  //Set up memory for the shape functions
  Shape psi(n_node,n_position_type);
  DShape dpsids(n_node,n_position_type,n_dim_el); 

  //Only need to call the local derivatives
  dshape_local(s,psi,dpsids);
  
  //Also calculate the surface Vectors (derivatives wrt local coordinates)
  DenseMatrix<double> interpolated_A(n_dim_el,n_dim,0.0);   
  
  //Calculate positions and derivatives
  for(unsigned l=0;l<n_node;l++) 
   {
    //Loop over positional dofs
    for(unsigned k=0;k<n_position_type;k++)
     {
      //Loop over coordinates
      for(unsigned i=0;i<n_dim;i++)
       {
        //Loop over LOCAL derivative directions, to calculate the tangent(s)
        for(unsigned j=0;j<n_dim_el;j++)
         {
          interpolated_A(j,i) += 
           nodal_position_gen(l,bulk_position_type(k),i)*dpsids(l,k,j);
         }
       }
     }
   }
  //Now find the local deformed metric tensor from the tangent Vectors
  DenseMatrix<double> A(n_dim_el,n_dim_el,0.0);
  for(unsigned i=0;i<n_dim_el;i++)
   {
    for(unsigned j=0;j<n_dim_el;j++)
     {
      //Take the dot product
      for(unsigned k=0;k<n_dim;k++)
       { 
        A(i,j) += interpolated_A(i,k)*interpolated_A(j,k);
       }
     }
   }
  
  //Find the determinant of the metric tensor
  double Adet =0.0;
  switch(n_dim_el) 
   {
   case 1:
    Adet = A(0,0);
    break;
   case 2:
    Adet = A(0,0)*A(1,1) - A(0,1)*A(1,0);
    break;
   default:
    throw 
     OomphLibError("Wrong dimension in FaceElement",
                   "FaceElement::J_eulerian()",
                   OOMPH_EXCEPTION_LOCATION);
   }
  
  // Return 
  return sqrt(Adet);
 }


//========================================================================
/// \short Compute the Jacobian of the mapping between the local and global
/// coordinates at the ipt-th integration point. Overloaded from 
/// FiniteElement.
//========================================================================
 double FaceElement::J_eulerian_at_knot(const unsigned &ipt)
  const
 {
  //Find the number of nodes
  const unsigned n_node = nnode();

  //Find the number of position types
  const unsigned n_position_type = nnodal_position_type();

  //Find the dimension of the node and element
  const unsigned n_dim = nodal_dimension();
  const unsigned n_dim_el = dim();

  //Set up dummy memory for the shape functions
  Shape psi(n_node,n_position_type);
  DShape dpsids(n_node,n_position_type,n_dim_el);

  //Get the shape functions and local derivatives at the knot
  //This call may use the stored versions, which is why this entire
  //function doesn't just call J_eulerian(s), after reading out s from
  //the knots.
  dshape_local_at_knot(ipt,psi,dpsids);

  //Also calculate the surface Vectors (derivatives wrt local coordinates)
  DenseMatrix<double> interpolated_A(n_dim_el,n_dim,0.0);   
  
  //Calculate positions and derivatives
  for(unsigned l=0;l<n_node;l++) 
   {
    //Loop over positional dofs
    for(unsigned k=0;k<n_position_type;k++)
     {
      //Loop over coordinates
      for(unsigned i=0;i<n_dim;i++)
       {
        //Loop over LOCAL derivative directions, to calculate the tangent(s)
        for(unsigned j=0;j<n_dim_el;j++)
         {
          interpolated_A(j,i) += 
           nodal_position_gen(l,bulk_position_type(k),i)*dpsids(l,k,j);
         }
       }
     }
   }

  //Now find the local deformed metric tensor from the tangent Vectors
  DenseMatrix<double> A(n_dim_el,n_dim_el,0.0);
  for(unsigned i=0;i<n_dim_el;i++)
   {
    for(unsigned j=0;j<n_dim_el;j++)
     {
      //Take the dot product
      for(unsigned k=0;k<n_dim;k++)
       { 
        A(i,j) += interpolated_A(i,k)*interpolated_A(j,k);
       }
     }
   }
  
  //Find the determinant of the metric tensor
  double Adet =0.0;
  switch(n_dim_el) 
   {
   case 1:
    Adet = A(0,0);
    break;
   case 2:
    Adet = A(0,0)*A(1,1) - A(0,1)*A(1,0);
    break;
   default:
    throw 
     OomphLibError("Wrong dimension in FaceElement",
                   "FaceElement::J_eulerian_at_knot()",
                   OOMPH_EXCEPTION_LOCATION);
   }
  
  // Return 
  return sqrt(Adet);
 }

//////////////// RAYRAY

//=======================================================================
/// Compute the tangent vector(s) at the specified local coordinate
//=======================================================================
void FaceElement::tangent(const Vector<double> &s,
                          Vector<Vector<double> > &tang_vec) const
{
  std::cout << "Hi from tangent (vector)" << std::endl;
  
  //Find the spatial dimension of the FaceElement
  const unsigned element_dim = dim();

  //Find the overall dimension of the problem 
  //(assume that it's the same for all nodes)
  const unsigned spatial_dim = nodal_dimension();

#ifdef PARANOID
  //Check the number of local coordinates passed
  if(s.size()!=element_dim)
   {
    std::ostringstream error_stream;
    error_stream
     << "Local coordinate s passed to tangent() has dimension " 
     << s.size() << std::endl
     << "but element dimension is " << element_dim << std::endl;
   
    throw OomphLibError(error_stream.str(),
                        "FaceElement::tangent()",
                        OOMPH_EXCEPTION_LOCATION);
   }

  //Check the dimension of the normal vector
  if(tang_vec[0].size()!=spatial_dim)
   {
    std::ostringstream error_stream;
    error_stream
     << "Unit normal passed to tangent() has dimension " 
     << tang_vec[0].size() << std::endl
     << "but spatial dimension is " << spatial_dim << std::endl;
   
    throw OomphLibError(error_stream.str(),
                        "FaceElement::tangent()",
                        OOMPH_EXCEPTION_LOCATION);
   }
#endif   

  //Now let's consider the different element dimensions
  switch(element_dim)
   {
    //Point element, derived from a 1D element, in this case
    //the tangent vector is merely the tangent to the bulk element
    //and there is only one free coordinate in the bulk element
    //Hence we will need to calculate the derivatives wrt the
    //local coordinates in the BULK element.
     case 0:
     {
       // Find the number of nodes in the Bulk element
       const unsigned n_node_bulk = Bulk_element_pt->nnode();
       //Find the number of position types in the bulk element
       const unsigned n_position_type_bulk = 
        Bulk_element_pt->nnodal_position_type();

       //Construct the local coordinate in the bulk element
       Vector<double> s_bulk(1);

       //Get the local coordinates in the bulk element
       get_local_coordinate_in_bulk(s,s_bulk);
       
      //Allocate storage for the shape functions and their derivatives wrt
      //local coordinates
      Shape psi(n_node_bulk,n_position_type_bulk);
      DShape dpsids(n_node_bulk,n_position_type_bulk,1);
      //Get the value of the shape functions at the given local coordinate
      Bulk_element_pt->dshape_local(s_bulk,psi,dpsids);

      
    //Calculate all derivatives of the spatial coordinates wrt 
    //local coordinates
    DenseMatrix<double> interpolated_dxds(1,spatial_dim);
    //Initialise to zero
    for(unsigned i=0;i<spatial_dim;i++) {interpolated_dxds(0,i) = 0.0;}
    
    //Loop over all parent nodes
    for(unsigned l=0;l<n_node_bulk;l++)
     {
      //Loop over all position types in the bulk
      for(unsigned k=0;k<n_position_type_bulk;k++)
       {
        //Loop over coordinate directions
        for(unsigned i=0;i<spatial_dim;i++)
         {
          //Compute the spatial derivative
          interpolated_dxds(0,i) += 
             Bulk_element_pt->nodal_position_gen(l,k,i)*dpsids(l,k,0);
         } // for
       } // for
     } // for

    //Now the unit normal is just the derivative of the position vector
    //with respect to the single coordinate
    for(unsigned i=0;i<spatial_dim;i++)
     {tang_vec[0][i] = interpolated_dxds(0,i);}
     } // case 0
     break;
   
    //Line element, derived from a 2D element, in this case
    //the normal is a mess of cross products
    //We need an interior direction, so we must find the local
    //derivatives in the BULK element
   case 1:
   {
    //Find the number of nodes in the bulk element
    const unsigned n_node_bulk = Bulk_element_pt->nnode();
    //Find the number of position types in the bulk element
    const unsigned n_position_type_bulk = 
     Bulk_element_pt->nnodal_position_type();
    
    //Construct the local coordinate in the bulk element
    Vector<double> s_bulk(2);
    //Get the local coordinates in the bulk element
    get_local_coordinate_in_bulk(s,s_bulk);
    
    //Allocate storage for the shape functions and their derivatives wrt
    //local coordinates
    Shape psi(n_node_bulk,n_position_type_bulk);
    DShape dpsids(n_node_bulk,n_position_type_bulk,2);
    //Get the value of the shape functions at the given local coordinate
    Bulk_element_pt->dshape_local(s_bulk,psi,dpsids);
 
    //Calculate all derivatives of the spatial coordinates 
    //wrt local coordinates
    DenseMatrix<double> interpolated_dxds(2,spatial_dim);
    //Initialise to zero
    for(unsigned j=0;j<2;j++)
     {for(unsigned i=0;i<spatial_dim;i++) {interpolated_dxds(j,i) = 0.0;}}
    
    //Loop over all parent nodes
    for(unsigned l=0;l<n_node_bulk;l++)
     {
      //Loop over all position types in the bulk
      for(unsigned k=0;k<n_position_type_bulk;k++)
       {
        //Loop over derivative direction
        for(unsigned j=0;j<2;j++)
         {
          //Loop over coordinate directions
          for(unsigned i=0;i<spatial_dim;i++)
           {
            //Compute the spatial derivative
            interpolated_dxds(j,i) += 
             Bulk_element_pt->nodal_position_gen(l,k,i)*dpsids(l,k,j);
           }
         }
       }
     }

    // RAYRAY I NEED HELP TO UNDERSTAND THIS... maybe

    //Initialise the tangent, interior tangent and normal vectors to zero
    //The idea is that even if the element is in a two-dimensional space,
    //the normal cannot be calculated without embedding the element in three
    //dimensions, in which case, the tangent and interior tangent will have
    //zero z-components.
    Vector<double> tangent(3,0.0), interior_tangent(3,0.0), normal(3,0.0);

    //We must get the relationship between the coordinate along the face
    //and the local coordinates in the bulk element
    //We must also find an interior direction
    DenseMatrix<double> dsbulk_dsface(2,1,0.0);
    unsigned interior_direction=0;
    get_ds_bulk_ds_face(s,dsbulk_dsface,interior_direction);
    //Load in the values for the tangents
    for(unsigned i=0;i<spatial_dim;i++)
     {
      //Tangent to the face is the derivative wrt to the face coordinate
      //which is calculated using dsbulk_dsface and the chain rule
      tangent[i] = interpolated_dxds(0,i)*dsbulk_dsface(0,0)
       + interpolated_dxds(1,i)*dsbulk_dsface(1,0);
      //Interior tangent to the face is the derivative in the interior 
      //direction
      interior_tangent[i] = interpolated_dxds(interior_direction,i);
     }

    //Now the (3D) normal to the element is the interior tangent 
    //crossed with the tangent
    normal[0] = 
     interior_tangent[1]*tangent[2] - interior_tangent[2]*tangent[1];
    normal[1] = 
     interior_tangent[2]*tangent[0] - interior_tangent[0]*tangent[2];
    normal[2] = 
     interior_tangent[0]*tangent[1] - interior_tangent[1]*tangent[0];
   
    //We find the line normal by crossing the element normal with the tangent
    Vector<double> full_normal(3);
    full_normal[0] = normal[1]*tangent[2] - normal[2]*tangent[1];
    full_normal[1] = normal[2]*tangent[0] - normal[0]*tangent[2];
    full_normal[2] = normal[0]*tangent[1] - normal[1]*tangent[0];

    //Copy the appropriate entries into the unit normal
    //Two or Three depending upon the spatial dimension of the system
    for(unsigned i=0;i<spatial_dim;i++) {tang_vec[0][i] = tangent[i];}



   } // case 1
   break;

   //Plane element, derived from 3D element, in this case the normal
   //is just the cross product of the two surface tangents
   //We assume, therefore, that we have three spatial coordinates
   //and two surface coordinates
   //Then we need only to get the derivatives wrt the local coordinates
   //in this face element
   case 2:
   {
#ifdef PARANOID
    //Check that we actually have three spatial dimensions
    if(spatial_dim != 3)
     {
      std::ostringstream error_stream;
      error_stream << "There are only " << spatial_dim
                   << "coordinates at the nodes of this 2D FaceElement,\n"
                   << "which must have come from a 3D Bulk element\n";
      throw OomphLibError(error_stream.str(),
                          "FaceElement::tangent()",
                          OOMPH_EXCEPTION_LOCATION);
     }
#endif


    //Find the number of nodes in the element
    const unsigned n_node  = this->nnode();
    //Find the number of position types
    const unsigned n_position_type = this->nnodal_position_type();
    
    //Allocate storage for the shape functions and their derivatives wrt
    //local coordinates
    Shape psi(n_node,n_position_type);
    DShape dpsids(n_node,n_position_type,2);
    //Get the value of the shape functions at the given local coordinate
    this->dshape_local(s,psi,dpsids);
 
    //Calculate all derivatives of the spatial coordinates 
    //wrt local coordinates
    DenseMatrix<double> interpolated_dxds(2,3);
    //Initialise to zero
    for(unsigned j=0;j<2;j++)
     {for(unsigned i=0;i<3;i++) {interpolated_dxds(j,i) = 0.0;}}
    
    //Loop over all nodes
    for(unsigned l=0;l<n_node;l++)
     {
      //Loop over all position types
      for(unsigned k=0;k<n_position_type;k++)
       {
        //Loop over derivative directions
        for(unsigned j=0;j<2;j++)
         {
          //Loop over coordinate directions
          for(unsigned i=0;i<3;i++)
           {
            //Compute the spatial derivative
            //Remember that we need to translate the position type
            //to its location in the bulk node
            interpolated_dxds(j,i) += 
             this->nodal_position_gen(l,bulk_position_type(k),i)*dpsids(l,k,j);
           }
         }
       }
     }//for

    tang_vec[0][0] = interpolated_dxds(0,0);
    tang_vec[0][1] = interpolated_dxds(0,1);
    tang_vec[0][2] = interpolated_dxds(0,2);

    tang_vec[1][0] = interpolated_dxds(1,0);
    tang_vec[1][1] = interpolated_dxds(1,1);
    tang_vec[1][2] = interpolated_dxds(1,2);

   }//case 2
    break;

   default:

    throw OomphLibError(
     "Cannot have a FaceElement with dimension higher than 2",
     "FaceElement::outer_unit_normal()",
     OOMPH_EXCEPTION_LOCATION);
    break;
   } // switch(element_dim)
}

//=======================================================================
/// Compute the tangent vector(s) at the ipt-th integration point
//=======================================================================
void FaceElement::tangent(const unsigned &ipt,
                          Vector<Vector<double> > &tang_vec) const
{
  std::cout << "Hi from tangent (ipt)" << std::endl;
  
  // Fine the dimension of the element
  const unsigned element_dim = dim();
  // Find the local coordinates of the ipt-th integration point
  Vector<double> s(element_dim);
  for(unsigned i=0;i<element_dim;i++) {s[i] = integral_pt()->knot(ipt,i);}
  // Call the tangent function
  tangent(s,tang_vec);
}



//=======================================================================
/// Compute the tangent and outer unit normal at the specified local coordinate
//=======================================================================
void FaceElement::tangent_and_outer_unit_normal(const Vector<double> &s,
                                                Vector<Vector<double> > &tang_vec,
                                                Vector<double> &unit_normal) const
{
 //Find the spatial dimension of the FaceElement
  const unsigned element_dim = dim();

  //Find the overall dimension of the problem 
  //(assume that it's the same for all nodes)
  const unsigned spatial_dim = nodal_dimension();
 
#ifdef PARANOID
  //Check the number of local coordinates passed
  if(s.size()!=element_dim)
   {
    std::ostringstream error_stream;
    error_stream
     << "Local coordinate s passed to outer_unit_normal() has dimension " 
     << s.size() << std::endl
     << "but element dimension is " << element_dim << std::endl;
   
    throw OomphLibError(error_stream.str(),
                        "FaceElement::outer_unit_normal()",
                        OOMPH_EXCEPTION_LOCATION);
   }
 
  //Check the dimension of the normal vector
  if(unit_normal.size()!=spatial_dim)
   {
    std::ostringstream error_stream;
    error_stream
     << "Unit normal passed to outer_unit_normal() has dimension " 
     << unit_normal.size() << std::endl
     << "but spatial dimension is " << spatial_dim << std::endl;
   
    throw OomphLibError(error_stream.str(),
                        "FaceElement::outer_unit_normal()",
                        OOMPH_EXCEPTION_LOCATION);
   }
#endif   


  //Now let's consider the different element dimensions
  switch(element_dim)
   {
    //Point element, derived from a 1D element, in this case
    //the normal is merely the tangent to the bulk element
    //and there is only one free coordinate in the bulk element
    //Hence we will need to calculate the derivatives wrt the
    //local coordinates in the BULK element.
   case 0:
   {
    //Find the number of nodes in the bulk element
    const unsigned n_node_bulk = Bulk_element_pt->nnode();
    //Find the number of position types in the bulk element
    const unsigned n_position_type_bulk = 
     Bulk_element_pt->nnodal_position_type();

    //Construct the local coordinate in the bulk element
    Vector<double> s_bulk(1);

    //Get the local coordinates in the bulk element
    get_local_coordinate_in_bulk(s,s_bulk);
    
    //Allocate storage for the shape functions and their derivatives wrt
    //local coordinates
    Shape psi(n_node_bulk,n_position_type_bulk);
    DShape dpsids(n_node_bulk,n_position_type_bulk,1);
    //Get the value of the shape functions at the given local coordinate
    Bulk_element_pt->dshape_local(s_bulk,psi,dpsids);
 
    //Calculate all derivatives of the spatial coordinates wrt 
    //local coordinates
    DenseMatrix<double> interpolated_dxds(1,spatial_dim);
    //Initialise to zero
    for(unsigned i=0;i<spatial_dim;i++) {interpolated_dxds(0,i) = 0.0;}
    
    //Loop over all parent nodes
    for(unsigned l=0;l<n_node_bulk;l++)
     {
      //Loop over all position types in the bulk
      for(unsigned k=0;k<n_position_type_bulk;k++)
       {
        //Loop over coordinate directions
        for(unsigned i=0;i<spatial_dim;i++)
         {
          //Compute the spatial derivative
          interpolated_dxds(0,i) += 
             Bulk_element_pt->nodal_position_gen(l,k,i)*dpsids(l,k,0);
         }
       }
     }
    
    // RAYGETN
    //Now the unit normal is just the derivative of the position vector
    //with respect to the single coordinate
    for(unsigned i=0;i<spatial_dim;i++)
     {unit_normal[i] = interpolated_dxds(0,i);}

    for(unsigned i=0;i<spatial_dim;i++)
     {tang_vec[0][i] = interpolated_dxds(0,i);}


   }
    break;

    //Line element, derived from a 2D element, in this case
    //the normal is a mess of cross products
    //We need an interior direction, so we must find the local
    //derivatives in the BULK element
   case 1:
   {
    //Find the number of nodes in the bulk element
    const unsigned n_node_bulk = Bulk_element_pt->nnode();
    //Find the number of position types in the bulk element
    const unsigned n_position_type_bulk = 
     Bulk_element_pt->nnodal_position_type();
    
    //Construct the local coordinate in the bulk element
    Vector<double> s_bulk(2);
    //Get the local coordinates in the bulk element
    get_local_coordinate_in_bulk(s,s_bulk);
    
    //Allocate storage for the shape functions and their derivatives wrt
    //local coordinates
    Shape psi(n_node_bulk,n_position_type_bulk);
    DShape dpsids(n_node_bulk,n_position_type_bulk,2);
    //Get the value of the shape functions at the given local coordinate
    Bulk_element_pt->dshape_local(s_bulk,psi,dpsids);
 
    //Calculate all derivatives of the spatial coordinates 
    //wrt local coordinates
    DenseMatrix<double> interpolated_dxds(2,spatial_dim);
    //Initialise to zero
    for(unsigned j=0;j<2;j++)
     {for(unsigned i=0;i<spatial_dim;i++) {interpolated_dxds(j,i) = 0.0;}}
    
    //Loop over all parent nodes
    for(unsigned l=0;l<n_node_bulk;l++)
     {
      //Loop over all position types in the bulk
      for(unsigned k=0;k<n_position_type_bulk;k++)
       {
        //Loop over derivative direction
        for(unsigned j=0;j<2;j++)
         {
          //Loop over coordinate directions
          for(unsigned i=0;i<spatial_dim;i++)
           {
            //Compute the spatial derivative
            interpolated_dxds(j,i) += 
             Bulk_element_pt->nodal_position_gen(l,k,i)*dpsids(l,k,j);
           }
         }
       }
     }
    
    //Initialise the tangent, interior tangent and normal vectors to zero
    //The idea is that even if the element is in a two-dimensional space,
    //the normal cannot be calculated without embedding the element in three
    //dimensions, in which case, the tangent and interior tangent will have
    //zero z-components.
    Vector<double> tangent(3,0.0), interior_tangent(3,0.0), normal(3,0.0);
    
    //We must get the relationship between the coordinate along the face
    //and the local coordinates in the bulk element
    //We must also find an interior direction
    DenseMatrix<double> dsbulk_dsface(2,1,0.0);
    unsigned interior_direction=0;
    get_ds_bulk_ds_face(s,dsbulk_dsface,interior_direction);
    //Load in the values for the tangents
    for(unsigned i=0;i<spatial_dim;i++)
     {
      //Tangent to the face is the derivative wrt to the face coordinate
      //which is calculated using dsbulk_dsface and the chain rule
      tangent[i] = interpolated_dxds(0,i)*dsbulk_dsface(0,0)
       + interpolated_dxds(1,i)*dsbulk_dsface(1,0);
      //Interior tangent to the face is the derivative in the interior 
      //direction
      interior_tangent[i] = interpolated_dxds(interior_direction,i);
     }

    //Now the (3D) normal to the element is the interior tangent 
    //crossed with the tangent
    normal[0] = 
     interior_tangent[1]*tangent[2] - interior_tangent[2]*tangent[1];
    normal[1] = 
     interior_tangent[2]*tangent[0] - interior_tangent[0]*tangent[2];
    normal[2] = 
     interior_tangent[0]*tangent[1] - interior_tangent[1]*tangent[0];
   
    //We find the line normal by crossing the element normal with the tangent
    Vector<double> full_normal(3);
    full_normal[0] = normal[1]*tangent[2] - normal[2]*tangent[1];
    full_normal[1] = normal[2]*tangent[0] - normal[0]*tangent[2];
    full_normal[2] = normal[0]*tangent[1] - normal[1]*tangent[0];

    //Copy the appropriate entries into the unit normal
    //Two or Three depending upon the spatial dimension of the system
    for(unsigned i=0;i<spatial_dim;i++) {unit_normal[i] = full_normal[i];}

    for(unsigned i=0;i<spatial_dim;i++) {tang_vec[0][i] = tangent[i];}
    
    // Normalise the tangent
    double tang_length = 0.0; 
    unsigned vec_i=0;
    for(unsigned dim_i=0;dim_i<spatial_dim;dim_i++)
      {tang_length += tang_vec[vec_i][dim_i]*tang_vec[vec_i][dim_i];}
    
    tang_length = sqrt(tang_length);

    for(unsigned dim_i=0;dim_i<spatial_dim;dim_i++)
      {tang_vec[vec_i][dim_i] /= tang_length;}
   }
   break;

   //Plane element, derived from 3D element, in this case the normal
   //is just the cross product of the two surface tangents
   //We assume, therefore, that we have three spatial coordinates
   //and two surface coordinates
   //Then we need only to get the derivatives wrt the local coordinates
   //in this face element
   case 2:
   {
#ifdef PARANOID
    //Check that we actually have three spatial dimensions
    if(spatial_dim != 3)
     {
      std::ostringstream error_stream;
      error_stream << "There are only " << spatial_dim
                   << "coordinates at the nodes of this 2D FaceElement,\n"
                   << "which must have come from a 3D Bulk element\n";
      throw OomphLibError(error_stream.str(),
                          "FaceElement::outer_unit_normal()",
                          OOMPH_EXCEPTION_LOCATION);
     }
#endif

    //Find the number of nodes in the element
    const unsigned n_node  = this->nnode();
    //Find the number of position types
    const unsigned n_position_type = this->nnodal_position_type();
    
    //Allocate storage for the shape functions and their derivatives wrt
    //local coordinates
    Shape psi(n_node,n_position_type);
    DShape dpsids(n_node,n_position_type,2);
    //Get the value of the shape functions at the given local coordinate
    this->dshape_local(s,psi,dpsids);
 
    //Calculate all derivatives of the spatial coordinates 
    //wrt local coordinates
    DenseMatrix<double> interpolated_dxds(2,3);
    //Initialise to zero
    for(unsigned j=0;j<2;j++)
     {for(unsigned i=0;i<3;i++) {interpolated_dxds(j,i) = 0.0;}}
    
    //Loop over all nodes
    for(unsigned l=0;l<n_node;l++)
     {
      //Loop over all position types
      for(unsigned k=0;k<n_position_type;k++)
       {
        //Loop over derivative directions
        for(unsigned j=0;j<2;j++)
         {
          //Loop over coordinate directions
          for(unsigned i=0;i<3;i++)
           {
            //Compute the spatial derivative
            //Remember that we need to translate the position type
            //to its location in the bulk node
            interpolated_dxds(j,i) += 
             this->nodal_position_gen(l,bulk_position_type(k),i)*dpsids(l,k,j);
           }
         }
       }
     }

    //We now take the cross product of the two normal vectors
    unit_normal[0] = 
     interpolated_dxds(0,1)*interpolated_dxds(1,2) -
     interpolated_dxds(0,2)*interpolated_dxds(1,1);
    unit_normal[1] = 
     interpolated_dxds(0,2)*interpolated_dxds(1,0) -
     interpolated_dxds(0,0)*interpolated_dxds(1,2);
    unit_normal[2] = 
     interpolated_dxds(0,0)*interpolated_dxds(1,1) -
     interpolated_dxds(0,1)*interpolated_dxds(1,0);
  
   /* 
    streamsize cout_precision = cout.precision();
    cout << setprecision(15) << interpolated_dxds(0,0) << " " 
                             << interpolated_dxds(0,1) << " "
                             << interpolated_dxds(0,2) << endl;

    cout << setprecision(15) << interpolated_dxds(1,0) << " "
                             << interpolated_dxds(1,1) << " "
                             << interpolated_dxds(1,2) << endl;

    cout << setprecision(cout_precision) << std::endl;
*/




    tang_vec[0][0] = interpolated_dxds(0,0);
    tang_vec[0][1] = interpolated_dxds(0,1);
    tang_vec[0][2] = interpolated_dxds(0,2);
    tang_vec[1][0] = interpolated_dxds(1,0);
    tang_vec[1][1] = interpolated_dxds(1,1);
    tang_vec[1][2] = interpolated_dxds(1,2);

//*
  // normalise 
  // Loop through the two vectors
  for(unsigned vec_i=0; vec_i<2; vec_i++)
  {
    // Get the length...
    double tang_length = 0.0;
    for(unsigned dim_i=0;dim_i<spatial_dim;dim_i++) 
     {tang_length += tang_vec[vec_i][dim_i]*tang_vec[vec_i][dim_i];}
    
    tang_length = sqrt(tang_length);

    for(unsigned dim_i=0;dim_i<spatial_dim;dim_i++) 
     {tang_vec[vec_i][dim_i] /= tang_length;}

  }
// */

/*
    streamsize cout_precision = cout.precision();
    cout << setprecision(15) << tang_vec[0][0] << " " 
                             << tang_vec[0][1] << " "
                             << tang_vec[0][2] << endl;

    cout << setprecision(15) << tang_vec[1][0] << " "
                             << tang_vec[1][1] << " "
                             << tang_vec[1][2] << endl;

    cout << setprecision(cout_precision) << std::endl;


pause("test done!");
// */
   }
    break;

   default:

    throw OomphLibError(
     "Cannot have a FaceElement with dimension higher than 2",
     "FaceElement::outer_unit_normal()",
     OOMPH_EXCEPTION_LOCATION);
    break;
   }
   
  //Finally normalise unit normal
  double length = 0.0;
  for(unsigned i=0;i<spatial_dim;i++) 
   {length += unit_normal[i]*unit_normal[i];}
  for(unsigned i=0;i<spatial_dim;i++) 
   {unit_normal[i] *= Normal_sign/sqrt(length);}

  // Lets normalise the tangent vectors as well... why not?



}

//=======================================================================
/// Compute the tangent and outer unit normal 
/// at the ipt-th integration point
//=======================================================================
void FaceElement::tangent_and_outer_unit_normal(const unsigned &ipt,
                                                Vector<Vector<double> > &tang_vec,
                                                Vector<double> &unit_normal) const
{
  // Fine the dimension of the element
  const unsigned element_dim = dim();
  // Find the local coordinates of the ipt-th integration point
  Vector<double> s(element_dim);
  for(unsigned i=0;i<element_dim;i++) {s[i] = integral_pt()->knot(ipt,i);}
  // Call the tangent function
  tangent_and_outer_unit_normal(s,tang_vec,unit_normal);

}





//=======================================================================
/// Compute the outer unit normal at the specified local coordinate
//=======================================================================
 void FaceElement::outer_unit_normal(const Vector<double> &s,
                                     Vector<double> &unit_normal) const
 {
  //Find the spatial dimension of the FaceElement
  const unsigned element_dim = dim();

  //Find the overall dimension of the problem 
  //(assume that it's the same for all nodes)
  const unsigned spatial_dim = nodal_dimension();
 
#ifdef PARANOID
  //Check the number of local coordinates passed
  if(s.size()!=element_dim)
   {
    std::ostringstream error_stream;
    error_stream
     << "Local coordinate s passed to outer_unit_normal() has dimension " 
     << s.size() << std::endl
     << "but element dimension is " << element_dim << std::endl;
   
    throw OomphLibError(error_stream.str(),
                        "FaceElement::outer_unit_normal()",
                        OOMPH_EXCEPTION_LOCATION);
   }
 
  //Check the dimension of the normal vector
  if(unit_normal.size()!=spatial_dim)
   {
    std::ostringstream error_stream;
    error_stream
     << "Unit normal passed to outer_unit_normal() has dimension " 
     << unit_normal.size() << std::endl
     << "but spatial dimension is " << spatial_dim << std::endl;
   
    throw OomphLibError(error_stream.str(),
                        "FaceElement::outer_unit_normal()",
                        OOMPH_EXCEPTION_LOCATION);
   }
#endif   

/*  //The spatial dimension of the bulk element will be element_dim+1
  const unsigned bulk_dim = element_dim + 1;

  //Find the number of nodes in the bulk element
  const unsigned n_node_bulk = Bulk_element_pt->nnode();
  //Find the number of position types in the bulk element
  const unsigned n_position_type_bulk = 
   Bulk_element_pt->nnodal_position_type();

  //Construct the local coordinate in the bulk element
  Vector<double> s_bulk(bulk_dim);
  //Set the value of the bulk coordinate that is fixed on the face
  //s_bulk[s_fixed_index()] = s_fixed_value();

  //Set the other bulk coordinates
  //for(unsigned i=0;i<element_dim;i++) {s_bulk[bulk_s_index(i)] = s[i];} 

  get_local_coordinate_in_bulk(s,s_bulk);

  //Allocate storage for the shape functions and their derivatives wrt
  //local coordinates
  Shape psi(n_node_bulk,n_position_type_bulk);
  DShape dpsids(n_node_bulk,n_position_type_bulk,bulk_dim);
  //Get the value of the shape functions at the given local coordinate
  Bulk_element_pt->dshape_local(s_bulk,psi,dpsids);
 
  //Calculate all derivatives of the spatial coordinates wrt local coordinates
  DenseMatrix<double> interpolated_dxds(bulk_dim,spatial_dim);
  //Initialise to zero
  for(unsigned j=0;j<bulk_dim;j++)
   {for(unsigned i=0;i<spatial_dim;i++) {interpolated_dxds(j,i) = 0.0;}}

  //Loop over all parent nodes
  for(unsigned l=0;l<n_node_bulk;l++)
   {
    //Loop over all position types in the bulk
    for(unsigned k=0;k<n_position_type_bulk;k++)
     {
      //Loop over derivative direction
      for(unsigned j=0;j<bulk_dim;j++)
       {
        //Loop over coordinate directions
        for(unsigned i=0;i<spatial_dim;i++)
         {
          //Compute the spatial derivative
          interpolated_dxds(j,i) += 
           Bulk_element_pt->nodal_position_gen(l,k,i)*dpsids(l,k,j);
         }
       }
     }
     }*/

  //Now let's consider the different element dimensions
  switch(element_dim)
   {
    //Point element, derived from a 1D element, in this case
    //the normal is merely the tangent to the bulk element
    //and there is only one free coordinate in the bulk element
    //Hence we will need to calculate the derivatives wrt the
    //local coordinates in the BULK element.
   case 0:
   {
    //Find the number of nodes in the bulk element
    const unsigned n_node_bulk = Bulk_element_pt->nnode();
    //Find the number of position types in the bulk element
    const unsigned n_position_type_bulk = 
     Bulk_element_pt->nnodal_position_type();

    //Construct the local coordinate in the bulk element
    Vector<double> s_bulk(1);

    //Get the local coordinates in the bulk element
    get_local_coordinate_in_bulk(s,s_bulk);
    
    //Allocate storage for the shape functions and their derivatives wrt
    //local coordinates
    Shape psi(n_node_bulk,n_position_type_bulk);
    DShape dpsids(n_node_bulk,n_position_type_bulk,1);
    //Get the value of the shape functions at the given local coordinate
    Bulk_element_pt->dshape_local(s_bulk,psi,dpsids);
 
    //Calculate all derivatives of the spatial coordinates wrt 
    //local coordinates
    DenseMatrix<double> interpolated_dxds(1,spatial_dim,0.0);
    
    //Loop over all parent nodes
    for(unsigned l=0;l<n_node_bulk;l++)
     {
      //Loop over all position types in the bulk
      for(unsigned k=0;k<n_position_type_bulk;k++)
       {
        //Loop over coordinate directions
        for(unsigned i=0;i<spatial_dim;i++)
         {
          //Compute the spatial derivative
          interpolated_dxds(0,i) += 
             Bulk_element_pt->nodal_position_gen(l,k,i)*dpsids(l,k,0);
         }
       }
     }

    //Now the unit normal is just the derivative of the position vector
    //with respect to the single coordinate
    for(unsigned i=0;i<spatial_dim;i++)
     {unit_normal[i] = interpolated_dxds(0,i);}
   }
    break;

    //Line element, derived from a 2D element, in this case
    //the normal is a mess of cross products
    //We need an interior direction, so we must find the local
    //derivatives in the BULK element
   case 1:
   {
    //Find the number of nodes in the bulk element
    const unsigned n_node_bulk = Bulk_element_pt->nnode();
    //Find the number of position types in the bulk element
    const unsigned n_position_type_bulk = 
     Bulk_element_pt->nnodal_position_type();
    
    //Construct the local coordinate in the bulk element
    Vector<double> s_bulk(2);
    //Get the local coordinates in the bulk element
    get_local_coordinate_in_bulk(s,s_bulk);
    
    //Allocate storage for the shape functions and their derivatives wrt
    //local coordinates
    Shape psi(n_node_bulk,n_position_type_bulk);
    DShape dpsids(n_node_bulk,n_position_type_bulk,2);
    //Get the value of the shape functions at the given local coordinate
    Bulk_element_pt->dshape_local(s_bulk,psi,dpsids);
 
    //Calculate all derivatives of the spatial coordinates 
    //wrt local coordinates
    DenseMatrix<double> interpolated_dxds(2,spatial_dim,0.0);
    
    //Loop over all parent nodes
    for(unsigned l=0;l<n_node_bulk;l++)
     {
      //Loop over all position types in the bulk
      for(unsigned k=0;k<n_position_type_bulk;k++)
       {
        //Loop over derivative direction
        for(unsigned j=0;j<2;j++)
         {
          //Loop over coordinate directions
          for(unsigned i=0;i<spatial_dim;i++)
           {
            //Compute the spatial derivative
            interpolated_dxds(j,i) += 
             Bulk_element_pt->nodal_position_gen(l,k,i)*dpsids(l,k,j);
           }
         }
       }
     }
    
    //Initialise the tangent, interior tangent and normal vectors to zero
    //The idea is that even if the element is in a two-dimensional space,
    //the normal cannot be calculated without embedding the element in three
    //dimensions, in which case, the tangent and interior tangent will have
    //zero z-components.
    Vector<double> tangent(3,0.0), interior_tangent(3,0.0), normal(3,0.0);
    
    //We must get the relationship between the coordinate along the face
    //and the local coordinates in the bulk element
    //We must also find an interior direction
    DenseMatrix<double> dsbulk_dsface(2,1,0.0);
    unsigned interior_direction=0;
    get_ds_bulk_ds_face(s,dsbulk_dsface,interior_direction);
    //Load in the values for the tangents
    for(unsigned i=0;i<spatial_dim;i++)
     {
      //Tangent to the face is the derivative wrt to the face coordinate
      //which is calculated using dsbulk_dsface and the chain rule
      tangent[i] = interpolated_dxds(0,i)*dsbulk_dsface(0,0)
       + interpolated_dxds(1,i)*dsbulk_dsface(1,0);
      //Interior tangent to the face is the derivative in the interior 
      //direction
      interior_tangent[i] = interpolated_dxds(interior_direction,i);
     }

    //Now the (3D) normal to the element is the interior tangent 
    //crossed with the tangent
    normal[0] = 
     interior_tangent[1]*tangent[2] - interior_tangent[2]*tangent[1];
    normal[1] = 
     interior_tangent[2]*tangent[0] - interior_tangent[0]*tangent[2];
    normal[2] = 
     interior_tangent[0]*tangent[1] - interior_tangent[1]*tangent[0];
   
    //We find the line normal by crossing the element normal with the tangent
    Vector<double> full_normal(3);
    full_normal[0] = normal[1]*tangent[2] - normal[2]*tangent[1];
    full_normal[1] = normal[2]*tangent[0] - normal[0]*tangent[2];
    full_normal[2] = normal[0]*tangent[1] - normal[1]*tangent[0];

    //Copy the appropriate entries into the unit normal
    //Two or Three depending upon the spatial dimension of the system
    for(unsigned i=0;i<spatial_dim;i++) {unit_normal[i] = full_normal[i];}
   }
   break;

   //Plane element, derived from 3D element, in this case the normal
   //is just the cross product of the two surface tangents
   //We assume, therefore, that we have three spatial coordinates
   //and two surface coordinates
   //Then we need only to get the derivatives wrt the local coordinates
   //in this face element
   case 2:
   {
#ifdef PARANOID
    //Check that we actually have three spatial dimensions
    if(spatial_dim != 3)
     {
      std::ostringstream error_stream;
      error_stream << "There are only " << spatial_dim
                   << "coordinates at the nodes of this 2D FaceElement,\n"
                   << "which must have come from a 3D Bulk element\n";
      throw OomphLibError(error_stream.str(),
                          "FaceElement::outer_unit_normal()",
                          OOMPH_EXCEPTION_LOCATION);
     }
#endif

    //Find the number of nodes in the element
    const unsigned n_node  = this->nnode();
    //Find the number of position types
    const unsigned n_position_type = this->nnodal_position_type();
    
    //Allocate storage for the shape functions and their derivatives wrt
    //local coordinates
    Shape psi(n_node,n_position_type);
    DShape dpsids(n_node,n_position_type,2);
    //Get the value of the shape functions at the given local coordinate
    this->dshape_local(s,psi,dpsids);
 
    //Calculate all derivatives of the spatial coordinates 
    //wrt local coordinates
    DenseMatrix<double> interpolated_dxds(2,3,0.0);
    
    //Loop over all nodes
    for(unsigned l=0;l<n_node;l++)
     {
      //Loop over all position types
      for(unsigned k=0;k<n_position_type;k++)
       {
        //Loop over derivative directions
        for(unsigned j=0;j<2;j++)
         {
          //Loop over coordinate directions
          for(unsigned i=0;i<3;i++)
           {
            //Compute the spatial derivative
            //Remember that we need to translate the position type
            //to its location in the bulk node
            interpolated_dxds(j,i) += 
             this->nodal_position_gen(l,bulk_position_type(k),i)*dpsids(l,k,j);
           }
         }
       }
     }

    //We now take the cross product of the two normal vectors
    unit_normal[0] = 
     interpolated_dxds(0,1)*interpolated_dxds(1,2) -
     interpolated_dxds(0,2)*interpolated_dxds(1,1);
    unit_normal[1] = 
     interpolated_dxds(0,2)*interpolated_dxds(1,0) -
     interpolated_dxds(0,0)*interpolated_dxds(1,2);
    unit_normal[2] = 
     interpolated_dxds(0,0)*interpolated_dxds(1,1) -
     interpolated_dxds(0,1)*interpolated_dxds(1,0);
   }
    break;

   default:

    throw OomphLibError(
     "Cannot have a FaceElement with dimension higher than 2",
     "FaceElement::outer_unit_normal()",
     OOMPH_EXCEPTION_LOCATION);
    break;
   }
   
  //Finally normalise unit normal
  double length = 0.0;
  for(unsigned i=0;i<spatial_dim;i++) 
   {length += unit_normal[i]*unit_normal[i];}
  for(unsigned i=0;i<spatial_dim;i++) 
   {unit_normal[i] *= Normal_sign/sqrt(length);}
 }

//=======================================================================
/// Compute the outer unit normal at the ipt-th integration point
//=======================================================================
 void FaceElement::outer_unit_normal(const unsigned &ipt,
                                     Vector<double> &unit_normal) const
 {
  //Find the dimension of the element
  const unsigned element_dim = dim();
  //Find the local coordiantes of the ipt-th integration point
  Vector<double> s(element_dim);
  for(unsigned i=0;i<element_dim;i++) {s[i] = integral_pt()->knot(ipt,i);}
  //Call the outer unit normal function
  outer_unit_normal(s,unit_normal);
 }



//=======================================================================
/// Return vector of local coordinates in bulk element, 
/// given the local coordinates in this FaceElement
//=======================================================================
Vector<double> FaceElement::local_coordinate_in_bulk(
  const Vector<double>& s) const
 {
  //Find the dimension of the bulk element
  unsigned dim_bulk = Bulk_element_pt->dim();

  // Vector of local coordinates in bulk element
  Vector<double> s_bulk(dim_bulk);

  //Use the function pointer if it is set
  if(Face_to_bulk_coordinate_fct_pt)
   {
    //Call the translation function
    (*Face_to_bulk_coordinate_fct_pt)(s,s_bulk);
   }
  else
   {
    throw OomphLibError("Face_to_bulk_coordinate mapping not set",
                        "FaceElement::local_coordinate_in_bulk()",
                        OOMPH_EXCEPTION_LOCATION);
   }
  
  // Return it
  return s_bulk;
 }



//=======================================================================
/// Calculate the  vector of local coordinates in bulk element, 
/// given the local coordinates in this FaceElement
//=======================================================================
void FaceElement::get_local_coordinate_in_bulk(
  const Vector<double>& s, Vector<double> &s_bulk) const
 {
  //Use the function pointer if it is set
  if(Face_to_bulk_coordinate_fct_pt)
   {
    //Call the translation function
    (*Face_to_bulk_coordinate_fct_pt)(s,s_bulk);
   }
  //Otherwise use the existing (not general) interface
  else
   {
        throw OomphLibError("Face_to_bulk_coordinate mapping not set",
                        "FaceElement::local_coordinate_in_bulk()",
                        OOMPH_EXCEPTION_LOCATION);
   }
 }


//=======================================================================
///  Calculate the derivatives of the local coordinates in the
/// bulk element with respect to the local coordinates in this FaceElement.
/// In addition return the index of a bulk local coordinate that varies away
/// from the face.
//=======================================================================
 void FaceElement::get_ds_bulk_ds_face(const Vector<double> &s,
                                       DenseMatrix<double> &dsbulk_dsface,
                                       unsigned &interior_direction) const
 {
  //Use the function pointer if it is set
  if(Bulk_coordinate_derivatives_fct_pt)
   {
    //Call the translation function
    (*Bulk_coordinate_derivatives_fct_pt)(s,dsbulk_dsface,interior_direction);
   }
  //Otherwise throw an error
  else
   {
    throw OomphLibError(
     "No function for derivatives of bulk coords wrt face coords set",
     "FaceElement::get_ds_bulk_ds_face()",
     OOMPH_EXCEPTION_LOCATION);
   }
 }



///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
//  Functions for elastic general elements
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////

//=========================================================================
/// Internal function that is used to assemble the jacobian of the mapping
/// from local coordinates (s) to the lagrangian coordinates (xi), given the
/// derivatives of the shape functions. 
//=========================================================================
 void SolidFiniteElement::
 assemble_local_to_lagrangian_jacobian(const DShape &dpsids,
                                       DenseMatrix<double> &jacobian) const
 {
  //Find the the dimension of the element
  const unsigned el_dim = dim();
  //Find the number of shape functions and shape functions types
  //We shall ASSUME (ENFORCE) that Lagrangian coordinates must 
  //be interpolated through the nodes
  const unsigned n_shape = nnode();
  const unsigned n_shape_type = nnodal_lagrangian_type();

#ifdef PARANOID
 //Check for dimensional compatibility
 if(el_dim != Lagrangian_dimension)
  {
   std::ostringstream error_message;
   error_message << "Dimension mismatch" << std::endl;
   error_message << "The elemental dimension: " << el_dim
                 << " must equal the nodal Lagrangian dimension: " 
                 << Lagrangian_dimension 
                 << " for the jacobian of the mapping to be well-defined"
                 << std::endl;
   throw OomphLibError(error_message.str(),
    "SolidFiniteElement::assemble_local_to_lagrangian_jacobian()",
                       OOMPH_EXCEPTION_LOCATION);
  }
#endif

  //Loop over the rows of the jacobian
  for(unsigned i=0;i<el_dim;i++)
   {
    //Loop over the columns of the jacobian
    for(unsigned j=0;j<el_dim;j++)
     {
      //Zero the entry
      jacobian(i,j) = 0.0;
      //Loop over the shape functions
      for(unsigned l=0;l<n_shape;l++)
       {
        for(unsigned k=0;k<n_shape_type;k++)
         {
          //Jacobian is dx_j/ds_i, which is represented by the sum
          //over the dpsi/ds_i of the nodal points X j
          //Call the Non-hanging version of positions
          //This is overloaded in refineable elements
          jacobian(i,j) += raw_lagrangian_position_gen(l,k,j)*dpsids(l,k,i);
         }
       }
     }
   }
 }

//=========================================================================
/// Internal function that is used to assemble the jacobian of second
/// derivatives of the the mapping from local coordinates (s) to the 
/// lagrangian coordinates (xi), given the second derivatives of the 
/// shape functions. 
//=========================================================================
 void SolidFiniteElement::
 assemble_local_to_lagrangian_jacobian2(const DShape &d2psids,
                                        DenseMatrix<double> &jacobian2) const
 {
  //Find the the dimension of the element
  const unsigned el_dim = dim();
  //Find the number of shape functions and shape functions types
  //We ENFORCE that Lagrangian coordinates must be interpolated through
  //the nodes
  const unsigned n_shape = nnode();
  const unsigned n_shape_type = nnodal_lagrangian_type();
  //Find the number of second derivatives
  const unsigned n_row = N2deriv[el_dim];
 
  //Assemble the "jacobian" (d^2 x_j/ds_i^2) for second derivatives of 
  //shape functions
  //Loop over the rows (number of second derivatives)
  for(unsigned i=0;i<n_row;i++)
   {
    //Loop over the columns (element dimension
    for(unsigned j=0;j<el_dim;j++)
     {
      //Zero the entry
      jacobian2(i,j) = 0.0;
      //Loop over the shape functions
      for(unsigned l=0;l<n_shape;l++)
       {
        //Loop over the shape function types
        for(unsigned k=0;k<n_shape_type;k++)
         {
          //Add the terms to the jacobian entry
          //Call the Non-hanging version of positions
          //This is overloaded in refineable elements
          jacobian2(i,j) += raw_lagrangian_position_gen(l,k,j)*d2psids(l,k,i);
         }
       }
     }
   }
 }

//============================================================================
/// \short Destructor for SolidFiniteElement:
//============================================================================
 SolidFiniteElement:: ~SolidFiniteElement() 
 {
  //Delete the storage allocated for the positional local equations
  delete[] Position_local_eqn;
 }


//==========================================================================
/// Calculate the mapping from local to lagrangian coordinates
/// assuming that the coordinates are aligned in the direction of the local 
/// coordinates, i.e. there are no cross terms and the jacobian is diagonal.
/// The local derivatives are passed as dpsids and the jacobian and 
/// inverse jacobian are returned.
//==========================================================================
 double SolidFiniteElement::
 local_to_lagrangian_mapping_diagonal(const DShape &dpsids,
                                      DenseMatrix<double> &jacobian,
                                      DenseMatrix<double> &inverse_jacobian)
  const
 {
  //Find the dimension of the element
  const unsigned el_dim = dim();
  //Find the number of shape functions and shape functions types
 //We shall ASSUME (ENFORCE) that Lagrangian coordinates must 
  //be interpolated through the nodes
  const unsigned n_shape = nnode();
  const unsigned n_shape_type = nnodal_lagrangian_type();
 
  //In this case we assume that there are no cross terms, that is
  //global coordinate 0 is always in the direction of local coordinate 0

  //Loop over the coordinates
  for(unsigned i=0;i<el_dim;i++)
   {
    //Zero the jacobian and inverse jacobian entries
    for(unsigned j=0;j<el_dim;j++) 
     {jacobian(i,j) = 0.0; inverse_jacobian(i,j) = 0.0;}
   
    //Loop over the shape functions
    for(unsigned l=0;l<n_shape;l++)
     {
      //Loop over the types of dof
      for(unsigned k=0;k<n_shape_type;k++)
       {
        //Derivatives are always dx_{i}/ds_{i}
        jacobian(i,i) += raw_lagrangian_position_gen(l,k,i)*dpsids(l,k,i);
       }
     }
   }
 
  //Now calculate the determinant of the matrix
  double det = 1.0;
  for(unsigned i=0;i<el_dim;i++) {det *= jacobian(i,i);}
 
//Report if Matrix is singular, or negative
#ifdef PARANOID
  check_jacobian(det);
#endif
 
  //Calculate the inverse mapping (trivial in this case)
  for(unsigned i=0;i<el_dim;i++) 
   {inverse_jacobian(i,i) = 1.0/jacobian(i,i);}
 
  //Return the value of the Jacobian
  return(det);
 }

//========================================================================
/// Calculate shape functions and derivatives w.r.t. Lagrangian 
/// coordinates at local coordinate s. Returns the Jacobian of the mapping
/// from Lagrangian to local coordinates.
/// General case, may be overloaded
//========================================================================
 double SolidFiniteElement::
 dshape_lagrangian(const Vector<double> &s, Shape &psi, DShape &dpsi) const
 {
  //Find the element dimension
  const unsigned el_dim = dim();
 
  //Get the values of the shape function and local derivatives
  //Temporarily stored in dpsi
  dshape_local(s,psi,dpsi);
 
  //Allocate memory for the inverse jacobian
  DenseMatrix<double> inverse_jacobian(el_dim);
  //Now calculate the inverse jacobian
  const double det = local_to_lagrangian_mapping(dpsi,inverse_jacobian);
 
  //Now set the values of the derivatives to be dpsidxi
  transform_derivatives(inverse_jacobian,dpsi);
  //Return the determinant of the jacobian
  return det;
 }

//========================================================================= 
/// Compute the geometric shape functions and also first
/// derivatives w.r.t. Lagrangian coordinates at integration point ipt.
/// Most general form of function, but may be over-loaded if desired
//========================================================================
 double SolidFiniteElement::dshape_lagrangian_at_knot(const unsigned &ipt,
                                                      Shape &psi, 
                                                      DShape &dpsi) const
 {
  //Find the element dimension
  const unsigned el_dim = dim();
 
  //Shape function for the local derivatives
  //Again we ASSUME (insist) that the lagrangian coordinates
  //are interpolated through the nodes
  //Get the values of the shape function and local derivatives
  dshape_local_at_knot(ipt,psi,dpsi);
 
  //Allocate memory for the inverse jacobian
  DenseMatrix<double> inverse_jacobian(el_dim);
  //Now calculate the inverse jacobian
  const double det = local_to_lagrangian_mapping(dpsi,inverse_jacobian);

  //Now set the values of the derivatives
  transform_derivatives(inverse_jacobian,dpsi);
  //Return the determinant of the jacobian
  return det;
 }

//========================================================================
/// \short Compute the geometric shape functions and also first
/// and second derivatives w.r.t. Lagrangian coordinates at 
/// local coordinate s;
/// Returns Jacobian of mapping from Lagrangian to local coordinates.
/// \n\n Numbering:
/// \n \b 1D: \n
/// d2pidxi(i,0) = \f$ d^2 \psi_j / d \xi^2 \f$
/// \n \b 2D: \n
/// d2psidxi(i,0) = \f$ \partial^2 \psi_j / \partial \xi_0^2 \f$ \n
/// d2psidxi(i,1) = \f$ \partial^2 \psi_j / \partial \xi_1^2 \f$ \n
/// d2psidxi(i,2) = \f$ \partial^2 \psi_j / \partial \xi_0 \partial \xi_1 \f$ \n
/// \n \b 3D: \n
/// d2psidxi(i,0) = \f$ \partial^2 \psi_j / \partial \xi_0^2 \f$ \n
/// d2psidxi(i,1) = \f$ \partial^2 \psi_j / \partial \xi_1^2 \f$ \n
/// d2psidxi(i,2) = \f$ \partial^2 \psi_j / \partial \xi_2^2 \f$ \n
/// d2psidxi(i,3) = \f$ \partial^2 \psi_j/\partial \xi_0 \partial \xi_1 \f$ \n
/// d2psidxi(i,4) = \f$ \partial^2 \psi_j/\partial \xi_0 \partial \xi_2 \f$ \n
/// d2psidxi(i,5) = \f$ \partial^2 \psi_j/\partial \xi_1 \partial \xi_2 \f$ \n
//========================================================================
 double SolidFiniteElement::
 d2shape_lagrangian(const Vector<double> &s, Shape &psi, 
                    DShape &dpsi, DShape &d2psi) const
 {
  //Find the element dimension
  const unsigned el_dim = dim();
  //Find the number of second derivatives required
  const unsigned n_deriv = N2deriv[el_dim];

  //Get the values of the shape function and local derivatives
  d2shape_local(s,psi,dpsi,d2psi);
 
  //Allocate memory for the jacobian and inverse jacobian 
  DenseMatrix<double> jacobian(el_dim), inverse_jacobian(el_dim);
  //Calculate the jacobian and inverse jacobian
  const double det = 
   local_to_lagrangian_mapping(dpsi,jacobian,inverse_jacobian);

  //Allocate memory for the jacobian of second derivatives
  DenseMatrix<double> jacobian2(n_deriv,el_dim);
  //Assemble the jacobian of second derivatives
  assemble_local_to_lagrangian_jacobian2(d2psi,jacobian2);

  //Now set the value of the derivatives
  transform_second_derivatives(jacobian,inverse_jacobian,
                               jacobian2,dpsi,d2psi);
  //Return the determinant of the mapping
  return det;
 }

//========================================================================== 
/// \short Compute the geometric shape functions and also first
/// and second derivatives w.r.t. Lagrangian coordinates at 
/// the ipt-th integration point
/// Returns Jacobian of mapping from Lagrangian to local coordinates.
/// \n\n Numbering:
/// \n \b 1D: \n
/// d2pidxi(i,0) = \f$ d^2 \psi_j / d \xi^2 \f$
/// \n \b 2D: \n
/// d2psidxi(i,0) = \f$ \partial^2 \psi_j / \partial \xi_0^2 \f$ \n
/// d2psidxi(i,1) = \f$ \partial^2 \psi_j / \partial \xi_1^2 \f$ \n
/// d2psidxi(i,2) = \f$ \partial^2 \psi_j/\partial \xi_0 \partial \xi_1 \f$ \n
/// \n \b 3D: \n
/// d2psidxi(i,0) = \f$ \partial^2 \psi_j / \partial \xi_0^2 \f$ \n
/// d2psidxi(i,1) = \f$ \partial^2 \psi_j / \partial \xi_1^2 \f$ \n
/// d2psidxi(i,2) = \f$ \partial^2 \psi_j / \partial \xi_2^2 \f$ \n
/// d2psidxi(i,3) = \f$ \partial^2 \psi_j/\partial \xi_0 \partial \xi_1 \f$ \n
/// d2psidxi(i,4) = \f$ \partial^2 \psi_j/\partial \xi_0 \partial \xi_2 \f$ \n
/// d2psidxi(i,5) = \f$ \partial^2 \psi_j/\partial \xi_1 \partial \xi_2 \f$ \n
//========================================================================
 double SolidFiniteElement::d2shape_lagrangian_at_knot(const unsigned &ipt, 
                                                       Shape &psi, 
                                                       DShape &dpsi, 
                                                       DShape &d2psi) const
 {
  //Find the values of the indices of the shape functions
  //Find the element dimension
  const unsigned el_dim = dim();
  //Find the number of second derivatives required
  const unsigned n_deriv = N2deriv[el_dim];
 
  //Get the values of the shape function and local derivatives
  d2shape_local_at_knot(ipt,psi,dpsi,d2psi);
 
  //Allocate memory for the jacobian and inverse jacobian 
  DenseMatrix<double> jacobian(el_dim), inverse_jacobian(el_dim);
  //Calculate the jacobian and inverse jacobian
  const double det = 
   local_to_lagrangian_mapping(dpsi,jacobian,inverse_jacobian);

  //Allocate memory for the jacobian of second derivatives
  DenseMatrix<double> jacobian2(n_deriv,el_dim);
  //Assemble the jacobian of second derivatives
  assemble_local_to_lagrangian_jacobian2(d2psi,jacobian2);

  //Now set the value of the derivatives
  transform_second_derivatives(jacobian,inverse_jacobian,
                               jacobian2,dpsi,d2psi);
  //Return the determinant of the mapping
  return det;
 }

//============================================================================
/// Assign local equation numbers for the solid equations in the element.
//  This can be done at a high level assuming, as I am, that the equations will
//  always be formulated in terms of nodal positions. 
//============================================================================
 void SolidFiniteElement::assign_solid_local_eqn_numbers()
 {
  //Find the number of nodes
  const unsigned n_node = nnode();
 
  //Check there are nodes!
  if(n_node > 0)
   {
    //Find the number of position types and dimensions of the nodes
    //Local caching
    const unsigned n_position_type = nnodal_position_type();
    const unsigned nodal_dim = nodal_dimension();
     
    //Delete the existing storage
    delete[] Position_local_eqn;
    //Resize the storage for the positional equation numbers
    Position_local_eqn = new int[n_node*n_position_type*nodal_dim];
    
    //A local queue to store the global equation numbers
    std::deque<unsigned long> global_eqn_number_queue;

    //Get the number of dofs so far, this must be outside both loops
    //so that both can use it
    unsigned local_eqn_number = ndof();
    
    //Loop over the nodes
    for(unsigned n=0;n<n_node;n++)
     {
      //Cast to a solid node
      SolidNode* cast_node_pt = static_cast<SolidNode*>(node_pt(n));

      //Loop over the number of position dofs
      for(unsigned j=0;j<n_position_type;j++)
       {
        //Loop over the dimension of each node
        for(unsigned k=0;k<nodal_dim;k++)
         {
          //Get equation number
          //Note eqn_number is long !
          long eqn_number = cast_node_pt->position_eqn_number(j,k);
          //If equation_number positive add to array
          if(eqn_number >= 0)
           {
            //Add to global array
            global_eqn_number_queue.push_back(eqn_number);
            //Add to look-up scheme
            Position_local_eqn[(n*n_position_type + j)*nodal_dim + k] = 
             local_eqn_number;
            //Increment the local equation number
            local_eqn_number++;
           }
          else
           {
            Position_local_eqn[(n*n_position_type + j)*nodal_dim + k] = 
             Data::Is_pinned;
           }
         }
       }
     } //End of loop over nodes

    //Now add our global equations numbers to the internal element storage
    add_global_eqn_numbers(global_eqn_number_queue);
   } //End of the case when there are nodes
 }


//============================================================================
/// This function calculates the entries of Jacobian matrix, used in 
/// the Newton method, associated with the elastic problem in which the
/// nodal position is a variable. It does this using finite differences, 
/// rather than an analytical formulation, so can be done in total generality.
//==========================================================================
 void SolidFiniteElement::
 fill_in_jacobian_from_solid_position_by_fd(Vector<double> &residuals,
                                            DenseMatrix<double> &jacobian)
 {

  // Flag to indicate if we use first or second order FD
  // bool use_first_order_fd=false;

  //Find the number of nodes
  const unsigned n_node = nnode();

  //If there aren't any nodes, then return straight away
  if(n_node == 0) {return;}

  //Call the update function to ensure that the element is in
  //a consistent state before finite differencing starts
  update_before_solid_position_fd();

  //Get the number of position dofs and dimensions at the node
  const unsigned n_position_type = nnodal_position_type();
  const unsigned nodal_dim = nodal_dimension();

  //Find the number of dofs in the element
  const unsigned n_dof = this->ndof();

  //Create newres vectors
  Vector<double> newres(n_dof);
  // Vector<double> newres_minus(n_dof);
  
  //Integer storage for local unknown
  int local_unknown=0;
 
  //Should probably give this a more global scope
  const double fd_step = Default_fd_jacobian_step;
 
  //Loop over the nodes
  for(unsigned n=0;n<n_node;n++)
   {
    //Loop over position dofs
    for(unsigned k=0;k<n_position_type;k++)
     {
      //Loop over dimension
      for(unsigned i=0;i<nodal_dim;i++)
       {
        //If the variable is free
        local_unknown = position_local_eqn(n,k,i);
        if(local_unknown >= 0)
         {
          //Store a pointer to the (generalised) Eulerian nodal position
          double* const value_pt = &(node_pt(n)->x_gen(k,i));

          //Save the old value of the (generalised) Eulerian nodal position
          const double old_var = *value_pt;
           
          //Increment the (generalised) Eulerian nodal position
          *value_pt += fd_step;
          
          // Perform any auxialiary node updates
          node_pt(n)->perform_auxiliary_node_update_fct();

          // Update any other slaved variables
          update_in_solid_position_fd(n);

          //Calculate the new residuals
          get_residuals(newres);
         
//          if (use_first_order_fd)
           {
            //Do forward finite differences
            for(unsigned m=0;m<n_dof;m++)
             {
              //Stick the entry into the Jacobian matrix
              jacobian(m,local_unknown) = (newres[m] - residuals[m])/fd_step;
             }
           }
//           else
//            {
//             //Take backwards step for the  (generalised) Eulerian nodal 
//             // position
//             node_pt(n)->x_gen(k,i) = old_var-fd_step;
           
//             //Calculate the new residuals at backward position
//             get_residuals(newres_minus);

//             //Do central finite differences
//             for(unsigned m=0;m<n_dof;m++)
//              {
//               //Stick the entry into the Jacobian matrix
//               jacobian(m,local_unknown) = 
//                (newres[m] - newres_minus[m])/(2.0*fd_step);
//              }
//            }

          //Reset the (generalised) Eulerian nodal position
          *value_pt = old_var;


          // Perform any auxialiary node updates
          node_pt(n)->perform_auxiliary_node_update_fct();
          
          // Reset any other slaved variables
          reset_in_solid_position_fd(n);
         }
       }
     }
   }

  //End of finite difference loop
  //Final reset of any slaved data
  reset_after_solid_position_fd();
 }

//============================================================================
/// Return i-th FE-interpolated Lagrangian coordinate at
/// local coordinate s.
//============================================================================
 double SolidFiniteElement::interpolated_xi(const Vector<double> &s, 
                                            const unsigned &i) const
 { 
  //Find the number of nodes
  const unsigned n_node = nnode();

  //Find the number of lagrangian types from the first node
  const unsigned n_lagrangian_type = nnodal_lagrangian_type();

  //Assign storage for the local shape function
  Shape psi(n_node,n_lagrangian_type);

  //Find the values of shape function
  shape(s,psi);

  //Initialise value of xi
  double interpolated_xi = 0.0;

  //Loop over the local nodes
  for(unsigned l=0;l<n_node;l++) 
   {
    //Loop over the number of dofs
    for(unsigned k=0;k<n_lagrangian_type;k++)
     {
      interpolated_xi += lagrangian_position_gen(l,k,i)*psi(l,k);
     }
   }
   
  return(interpolated_xi);
 }

//============================================================================
/// Compute FE-interpolated Lagrangian coordinate vector xi[] at
/// local coordinate s.
//============================================================================
 void SolidFiniteElement::interpolated_xi(const Vector<double> &s, 
                                          Vector<double> &xi) const
 { 
  //Find the number of nodes
  const unsigned n_node = nnode();

  //Find the number of lagrangian types from the first node
  const unsigned n_lagrangian_type = nnodal_lagrangian_type();

  //Assign storage for the local shape function
  Shape psi(n_node,n_lagrangian_type);

  //Find the values of shape function
  shape(s,psi);

  //Read out the number of lagrangian coordinates from the node
  const unsigned n_lagrangian = lagrangian_dimension();

  //Loop over the number of lagrangian coordinates
  for(unsigned i=0;i<n_lagrangian;i++)
   {
    //Initialise component to zero
    xi[i] = 0.0;

    //Loop over the local nodes
    for(unsigned l=0;l<n_node;l++) 
     {
      //Loop over the number of dofs
      for(unsigned k=0;k<n_lagrangian_type;k++)
       {
        xi[i] += lagrangian_position_gen(l,k,i)*psi(l,k);
       }
     }
   }
 }





//============================================================================
/// Compute derivatives of FE-interpolated Lagrangian coordinates xi
/// with respect to local coordinates: dxids[i][j]=dxi_i/ds_j
//============================================================================
void SolidFiniteElement::interpolated_dxids(const Vector<double> &s, 
                                            DenseMatrix<double> &dxids) const
{ 
 //Find the number of nodes
 const unsigned n_node = nnode();
 
 //Find the number of lagrangian types from the first node
 const unsigned n_lagrangian_type = nnodal_lagrangian_type();
 
 // Dimension of the element =  number of local coordinates
 unsigned el_dim=dim();
 
 //Assign storage for the local shape function
 Shape psi(n_node,n_lagrangian_type);
 DShape dpsi(n_node,n_lagrangian_type,el_dim);
 
 //Find the values of shape function and its derivs w.r.t. to local coords
 dshape_local(s,psi,dpsi);
 
 //Read out the number of lagrangian coordinates from the node
 const unsigned n_lagrangian = lagrangian_dimension();
 
 //Loop over the number of lagrangian and local coordinates 
 for(unsigned i=0;i<n_lagrangian;i++)
  {
   for(unsigned j=0;j<el_dim;j++)
    {
     //Initialise component to zero
     dxids(i,j) = 0.0;
     
     //Loop over the local nodes
     for(unsigned l=0;l<n_node;l++) 
      {
       //Loop over the number of dofs
       for(unsigned k=0;k<n_lagrangian_type;k++)
        {
         dxids(i,j) += lagrangian_position_gen(l,k,i)*dpsi(l,k,j);
        }
      }
    }
  }
}

//=======================================================================
/// Add jacobian and residuals for consistent assignment of 
/// initial "accelerations" in Newmark scheme. Jacobian is the mass matrix.
//=======================================================================
void SolidFiniteElement::
fill_in_jacobian_for_newmark_accel(DenseMatrix<double> &jacobian)
{
 
#ifdef PARANOID
 // Check that we're computing the real residuals, not the
 // ones corresponding to the assignement of a prescribed displacement
 if ((Solid_ic_pt!=0)||(!Solve_for_consistent_newmark_accel_flag))
  {
   std::ostringstream error_stream;
   error_stream << "Can only call fill_in_jacobian_for_newmark_accel()\n" 
                << "With Solve_for_consistent_newmark_accel_flag:" 
                << Solve_for_consistent_newmark_accel_flag << std::endl;
   error_stream << "Solid_ic_pt: " << Solid_ic_pt << std::endl;

   throw OomphLibError(error_stream.str(),
                       "SolidFinitElement::fill_in_jacobian_for_newmark_accel()",
                       OOMPH_EXCEPTION_LOCATION);
  }
#endif

 //Find the number of nodes
 const unsigned n_node = nnode();
 const unsigned n_position_type = nnodal_position_type();
 const unsigned nodal_dim = nodal_dimension();

 //Set the number of Lagrangian coordinates
 const unsigned n_lagrangian = dim();
  
 //Integer storage for local equation and unknown
 int local_eqn=0, local_unknown=0;

 // Set up memory for the shape functions:
 // # of nodes, # of positional dofs
 Shape psi(n_node,n_position_type);

 // # of nodes, # of positional dofs, # of lagrangian coords (for deriv)
 // Not needed but they come for free when compute the Jacobian.
 DShape dpsidxi(n_node,n_position_type,n_lagrangian);

 //Set # of integration points
 const unsigned n_intpt = integral_pt()->nweight();

 //Set the Vector to hold local coordinates
 Vector<double> s(n_lagrangian);

 // Get positional timestepper from first nodal point
 TimeStepper* timestepper_pt= node_pt(0)->position_time_stepper_pt();

#ifdef PARANOID
 // Of course all this only works if we're actually using a
 // Newmark timestepper!
 if (timestepper_pt->type()!="Newmark")
  {
   std::ostringstream error_stream;
   error_stream
    << "Assignment of Newmark accelerations obviously only works\n"
    << "for Newmark timesteppers. You've called me with " 
    << timestepper_pt->type() << std::endl;

   throw OomphLibError(error_stream.str(),
                       "SolidFinitElement::fill_in_jacobian_for_newmark_accel()",
                       OOMPH_EXCEPTION_LOCATION);
  }
#endif

 // "Acceleration" is the last history value:
 const unsigned ntstorage=timestepper_pt->ntstorage();
 const double accel_weight=timestepper_pt->weight(2,ntstorage-1);

 //Loop over the integration points
 for(unsigned ipt=0;ipt<n_intpt;ipt++)
  {
     
   //Assign values of s
   for(unsigned i=0;i<n_lagrangian;i++) {s[i] = integral_pt()->knot(ipt,i);}
   
   //Get the integral weight
   double w = integral_pt()->weight(ipt);

   // Jacobian of mapping between local and Lagrangian coords and shape
   // functions
   double J = dshape_lagrangian(s,psi,dpsidxi);

   // Get Lagrangian coordinate
   Vector<double> xi(n_lagrangian);
   interpolated_xi(s,xi);

   // Get multiplier for inertia terms
   double factor=multiplier(xi);

   //Premultiply the weights and the Jacobian
   double W = w*J;

   //Loop over the nodes
   for(unsigned l0=0;l0<n_node;l0++)
    {
     //Loop over the positional dofs: 'Type': 0: displacement; 1: deriv
     for(unsigned k0=0;k0<n_position_type;k0++)
      {
       // Loop over Eulerian directions
       for(unsigned i0=0;i0<nodal_dim;i0++)
        {
         local_eqn = position_local_eqn(l0,k0,i0);
         //If it's a degree of freedom, add contribution
         if(local_eqn >= 0)
          {
           //Loop over the nodes
           for(unsigned l1=0;l1<n_node;l1++)
            {
             //Loop over the positional dofs: 'Type': 0: displacement; 
             //1: deriv
             for(unsigned k1=0;k1<n_position_type;k1++)
              {
               // Loop over Eulerian directions
               unsigned i1=i0;
               {
                local_unknown = position_local_eqn(l1,k1,i1);
                //If it's a degree of freedom, add contribution
                if(local_unknown >= 0)
                 {
                  // Add contribution: Mass matrix, multiplied by
                  // weight for "acceleration" in Newmark scheme
                  // and general multiplier 
                  jacobian(local_eqn,local_unknown) +=
                   factor*accel_weight*psi(l0,k0)*psi(l1,k1)*W;
                 }
               }
              }
            }
          }
        }
      }
    }
   
  } //End of loop over the integration points

}





//=======================================================================
/// Helper function to fill in the residuals and (if flag==1) the Jacobian 
/// for the setup of an initial condition. The global equations are:
/// \f[
/// 0 = \int \left( \sum_{j=1}^N \sum_{k=1}^K X_{ijk} \psi_{jk}(\xi_n) 
/// - \frac{\partial^D R^{(IC)}_i(\xi_n)}{\partial t^D}
/// \right) \psi_{lm}(\xi_n) \ dv
/// \mbox{ \ \ \ \ for \ \ \ $l=1,...,N, \ \ m=1,...,K$}
/// \f]
/// where \f$ N \f$ is the number of nodes in the mesh and \f$ K \f$
/// the number of generalised nodal coordinates. The initial shape
/// of the solid body, \f$ {\bf R}^{(IC)},\f$ and its time-derivatives
/// are specified via the \c GeomObject that is stored in the 
/// \c SolidFiniteElement::SolidInitialCondition object. The latter also
/// stores the order of the time-derivative \f$ D \f$ to be assigned.
//=======================================================================
void SolidFiniteElement::fill_in_generic_jacobian_for_solid_ic(
 Vector<double> &residuals,
 DenseMatrix<double> &jacobian,
 const unsigned& flag)
{
 //Find the number of nodes and position types
 const unsigned n_node = nnode();
 const unsigned n_position_type = nnodal_position_type();

 //Set the dimension of the global coordinates
 const unsigned nodal_dim = nodal_dimension();
 
 //Find the number of lagragian types from the first node
 const unsigned n_lagrangian_type = nnodal_lagrangian_type();

 //Set the number of lagrangian coordinates
 const unsigned n_lagrangian = dim();
  
 //Integer storage for local equation number
 int local_eqn=0;
 int local_unknown=0;

 // # of nodes, # of positional dofs
 Shape psi(n_node,n_position_type);
 
 // # of nodes, # of positional dofs, # of lagrangian coords (for deriv)
 // not needed but they come for free when we compute the Jacobian
 //DShape dpsidxi(n_node,n_position_type,n_lagrangian);
 
 //Set # of integration points
 const unsigned n_intpt = integral_pt()->nweight();
 
 //Set the Vector to hold local coordinates
 Vector<double> s(n_lagrangian);
 
 //Loop over the integration points
 for(unsigned ipt=0;ipt<n_intpt;ipt++)
  {

   //Assign values of s
   for(unsigned i=0;i<n_lagrangian;i++) 
    {s[i] = integral_pt()->knot(ipt,i);}

   //Get the integral weight
   double w = integral_pt()->weight(ipt);

   // Shape fcts
   shape(s,psi);

   // Get Lagrangian coordinate
   Vector<double> xi(n_lagrangian,0.0);

   //Loop over the number of lagrangian coordinates
   for(unsigned i=0;i<n_lagrangian;i++)
    {
     //Loop over the local nodes
     for(unsigned l=0;l<n_node;l++) 
      {
       //Loop over the number of dofs
       for(unsigned k=0;k<n_lagrangian_type;k++)
        {
         xi[i] += lagrangian_position_gen(l,k,i)*psi(l,k);
        }
      }
    }

   // Get initial condition
   Vector<double> drdt_ic(nodal_dim);
   Solid_ic_pt->geom_object_pt()->
    dposition_dt(xi,Solid_ic_pt->ic_time_deriv(),drdt_ic);
   
   // Weak form of assignment of initial guess
   
   //Loop over the number of node
   for (unsigned l=0;l<n_node;l++)
    {
     //Loop over the type of degree of freedom
     for (unsigned k=0;k<n_position_type;k++)
      {
       //Loop over the coordinate directions
       for (unsigned i=0;i<nodal_dim;i++)
        {
         local_eqn = position_local_eqn(l,k,i);

         //If it's not a boundary condition
         if(local_eqn >= 0)
          {
           // Note we're ignoring the mapping between local and
           // global Lagrangian coords -- doesn't matter here;
           // we're just enforcing a slightly different
           // weak form.
           residuals[local_eqn] += (interpolated_x(s,i)-drdt_ic[i])*
            psi(l,k)*w;


           // Do Jacobian too?
           if (flag==1)
            {     

             //Loop over the number of node
             for (unsigned ll=0;ll<n_node;ll++)
              {                
               //Loop over the type of degree of freedom
               for (unsigned kk=0;kk<n_position_type;kk++)
                {    

                 // Only diagonal term
                 unsigned ii=i;
                                 
                 local_unknown = position_local_eqn(ll,kk,ii);
                 
                 //If it's not a boundary condition
                 if(local_unknown >= 0)
                  {
                   // Note we're ignoring the mapping between local and
                   // global Lagrangian coords -- doesn't matter here;
                   // we're just enforcing a slightly different
                   // weak form.
                   jacobian(local_eqn,local_unknown) +=
                    psi(ll,kk)*psi(l,k)*w;
                  }
                 else
                  {
                   oomph_info 
                    << "WARNING: You should really free all Data" 
                    << std::endl
                    << "         before setup of initial guess" 
                    << std::endl
                    << "ll, kk, ii " << ll << " " << kk << " " << ii 
                    << std::endl;
                  }
                }
              }
             
            }
          }
         else
          {
           oomph_info 
            << "WARNING: You should really free all Data" << std::endl
            << "         before setup of initial guess" << std::endl
            << "l, k, i " << l << " " << k << " " << i << std::endl;
          }
        }
      }
    }
  
  } //End of loop over the integration points
}


//===============================================================
/// Return the geometric shape function at the local coordinate s
//===============================================================
 void PointElement::shape(const Vector<double> &s, Shape &psi) const
  {
   // Single shape function always has value 1
   psi[0]=1.0;
  }
 
//=======================================================================
/// Assign the static Default_integration_scheme
//=======================================================================
PointIntegral PointElement::Default_integration_scheme;


}
