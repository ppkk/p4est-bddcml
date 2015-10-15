#ifndef P4_TO_P8
#include <p4est_bits.h>
#include <p4est_ghost.h>
#include <p4est_lnodes.h>
#include <p4est_vtk.h>
#else
#include <p8est_bits.h>
#include <p8est_ghost.h>
#include <p8est_lnodes.h>
#include <p8est_vtk.h>
#endif

#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <assert.h>

#include "arrays.h"
#include "bddcml_structs.h"
#include "my_p4est_interface.h"
#include "assemble.h"
#include "mesh.h"
#include "femspace.h"
#include "quadrature.h"
#include "element.h"

using namespace std;

void print_complete_matrix_rhs(BddcmlFemSpace *femsp, BddcmlDimensions *global_dims, SparseMatrix *matrix, RealArray *rhss, MPI_Comm mpicomm)
{
   real* compl_rhs = (real*) malloc(global_dims->n_dofs * sizeof(real));
   real* compl_mat = (real*) malloc(global_dims->n_dofs * global_dims->n_dofs * sizeof(real));

   if(!compl_rhs || !compl_mat)
   {
      printf("Unable to allocate dense matrix. Exiting print_complete_matrix_rhs....\n");
      return;
   }

   memset(compl_rhs, 0, global_dims->n_dofs * sizeof(real));
   memset(compl_mat, 0, global_dims->n_dofs * global_dims->n_dofs * sizeof(real));

   for(int loc_dof = 0; loc_dof < femsp->subdomain_dims->n_dofs; loc_dof++)
   {
      int glob_dof = femsp->dofs_global_map.val[loc_dof];
      compl_rhs[glob_dof] = rhss->val[loc_dof];
   }

   for(int idx = 0; idx < matrix->nnz; idx++)
   {
      int glob_dof_j = femsp->dofs_global_map.val[matrix->j[idx]];
      int glob_dof_i = femsp->dofs_global_map.val[matrix->i[idx]];
      compl_mat[global_dims->n_dofs * glob_dof_j + glob_dof_i] += matrix->val[idx];

      // adding also symmetric counterpart
      if(matrix->type != MatrixType::GENERAL)
         compl_mat[global_dims->n_dofs * glob_dof_i + glob_dof_j] += matrix->val[idx];
   }

   if(mpi_rank == 0)
   {
      MPI_Reduce(MPI_IN_PLACE, compl_rhs, global_dims->n_dofs, MPI_DOUBLE, MPI_SUM, 0, mpicomm);
      MPI_Reduce(MPI_IN_PLACE, compl_mat, global_dims->n_dofs * global_dims->n_dofs, MPI_DOUBLE, MPI_SUM, 0, mpicomm);
   }
   else
   {
      MPI_Reduce(compl_rhs, compl_rhs, global_dims->n_dofs, MPI_DOUBLE, MPI_SUM, 0, mpicomm);
      MPI_Reduce(compl_mat, compl_mat, global_dims->n_dofs * global_dims->n_dofs, MPI_DOUBLE, MPI_SUM, 0, mpicomm);
   }

   if(mpi_rank == 0)
   {
      FILE *f = fopen("discr.txt", "w");
      for(int dof_j = 0; dof_j < global_dims->n_dofs; dof_j++)
      {
         fprintf(f, "dof %d, rhs: %g\n", dof_j, compl_rhs[dof_j]);
         for(int dof_i = 0; dof_i < global_dims->n_dofs; dof_i++)
         {
            fprintf(f, "%g\n", compl_mat[dof_j*global_dims->n_dofs + dof_i]);
         }
      }
      fclose(f);
   }

   free(compl_rhs);
   free(compl_mat);
}


void ref_value_1D(int loc_id_1d, double x, double elem_len, double & value, double & der)
{
   if(loc_id_1d == 0)
   {
      value = (1-x)/2.;
      der = -1/elem_len;
   }
   else if(loc_id_1d == 1)
   {
      value = (1+x)/2.;
      der = 1/elem_len;
   }
   else
      assert(0);
}

void prepare_transformed_values(Quadrature q, double element_length,
                                vector<vector<double> > &values, vector<vector<vector<double> > > &gradients)
{
   values = vector<vector<double> >(P4EST_CHILDREN);
   gradients = vector<vector<vector<double> > >(P4EST_CHILDREN);

   for(int node = 0; node < P4EST_CHILDREN; node++)
   {
      int x_id_1D = node % 2;
      int y_id_1D = (node % 4) / 2;
      int z_id_1D = node / 4;

      for(unsigned int q_idx = 0; q_idx < q.weights.size(); q_idx++)
      {
         double value_x, der_x, value_y, der_y, value_z, der_z;
         ref_value_1D(x_id_1D, q.coords[q_idx][0], element_length, value_x, der_x);
         ref_value_1D(y_id_1D, q.coords[q_idx][1], element_length, value_y, der_y);

#ifdef P4_TO_P8
         ref_value_1D(z_id_1D, (q.coords[q_idx])[2], element_length, value_z, der_z);
#endif

         double value = value_x * value_y;
         double grad_1 = der_x * value_y;
         double grad_2 = value_x * der_y;

#ifdef P4_TO_P8
         value *= value_z;
         grad_1 *= value_z;
         grad_2 *= value_z;
         double grad_3 = value_x * value_y * der_z;
         vector<double> aux1;
         aux1.push_back(grad_1); 
         aux1.push_back(grad_2); 
         aux1.push_back(grad_3); 
         //gradients[node].push_back(vector<double>({grad_1, grad_2, grad_3}));
         gradients[node].push_back(aux1);
#else
         vector<double> aux1;
         aux1.push_back(grad_1); 
         aux1.push_back(grad_2); 
         //gradients[node].push_back(vector<double>({grad_1, grad_2}));
         gradients[node].push_back(aux1);
#endif
         values[node].push_back(value);
      }
   }
}

void zero_matrix_rhs(vector<vector<vector<vector<real> > > > &matrix, vector<vector<real> > &rhs, int num_comp)
{
   for(int i_node = 0; i_node < P4EST_CHILDREN; i_node++)
   {
      for(int i_comp = 0; i_comp < num_comp; i_comp++)
      {
         rhs[i_node][i_comp] = 0.0;
         for(int j_node = 0; j_node < P4EST_CHILDREN; j_node++)
         {
            for(int j_comp = 0; j_comp < num_comp; j_comp++)
            {
               matrix[i_node][i_comp][j_node][j_comp] = 0.0;
            }
         }
      }
   }
}

void print_matrix_rhs(vector<vector<vector<vector<real> > > > &matrix, vector<vector<real> > &rhs, int num_comp)
{
   cout << "LOCAL MATRIX:" << endl;
   for(int i_node = 0; i_node < P4EST_CHILDREN; i_node++)
   {
      for(int i_comp = 0; i_comp < num_comp; i_comp++)
      {
         for(int j_node = 0; j_node < P4EST_CHILDREN; j_node++)
         {
            for(int j_comp = 0; j_comp < num_comp; j_comp++)
            {
               cout << "nodes (" << i_node << ", " << j_node << "), comp (" << i_comp << ", " << j_comp << ") -> " <<
                       matrix[i_node][i_comp][j_node][j_comp] << endl;
            }
         }
      }
   }
   cout << "END LOCAL MATRIX:" << endl;
}

void assemble_local_laplace(Element* element, vector<vector<vector<vector<real> > > > &matrix,
                    vector<vector<real> > &rhs, RhsPtr rhs_ptr)
{
   GaussQuadrature q(P4EST_DIM, QUAD_ORDER);

   vector<vector<double> > values;
   vector<vector<vector<double> > > gradients;

   prepare_transformed_values(q, element->size, values, gradients);

   q.transform_to_physical(element);

   zero_matrix_rhs(matrix, rhs, 1);

   for(unsigned int q_idx = 0; q_idx < q.weights.size(); q_idx++)
   {
      for(int i_node = 0; i_node < P4EST_CHILDREN; i_node++)
      {

         rhs[i_node][0] += q.weights[q_idx] * values[i_node][q_idx] * rhs_ptr(q.coords[q_idx])[0];

         for(int j_node = 0; j_node < P4EST_CHILDREN; j_node++)
         {
            for(int idx_dim = 0; idx_dim < P4EST_DIM; idx_dim++)
            {
               matrix[i_node][0][j_node][0] += q.weights[q_idx] * gradients[i_node][q_idx][idx_dim] * gradients[j_node][q_idx][idx_dim];
            }
         }
      }
   }

}

void assemble_local_elasticity(Element* element, vector<vector<vector<vector<real> > > > &matrix,
                    vector<vector<real> > &rhs, RhsPtr rhs_ptr, Parameters p)
{
   GaussQuadrature q(P4EST_DIM, QUAD_ORDER);

   int num_comp = P4EST_DIM;
   vector<vector<double> > vals; //[node][quadrature_pt]
   vector<vector<vector<double> > > grads;//[node][quadrature_pt][component]

   prepare_transformed_values(q, element->size, vals, grads);

   q.transform_to_physical(element);

   zero_matrix_rhs(matrix, rhs, P4EST_DIM);

   for(unsigned int q_idx = 0; q_idx < q.weights.size(); q_idx++)
   {
      for(int i_node = 0; i_node < P4EST_CHILDREN; i_node++)
      {
         for(int i_comp = 0; i_comp < num_comp; i_comp++)
         {
            rhs[i_node][i_comp] += q.weights[q_idx] * vals[i_node][q_idx] * rhs_ptr(q.coords[q_idx])[i_comp];

            for(int j_node = 0; j_node < P4EST_CHILDREN; j_node++)
            {
               for(int j_comp = 0; j_comp < num_comp; j_comp++)
               {
                  double contrib = 0.0;

                  if(i_comp == j_comp)
                  {
                     contrib += (p.lambda + 2 * p.mu) * grads[i_node][q_idx][i_comp] * grads[j_node][q_idx][j_comp];
                     for(int other_component = 0; other_component < num_comp; other_component++)
                     {
                        if(other_component == i_comp)
                           continue;
                        contrib += p.mu * grads[i_node][q_idx][other_component] * grads[j_node][q_idx][other_component];

                     }
                  }
                  else
                  {
                     contrib += p.mu     * grads[i_node][q_idx][j_comp] * grads[j_node][q_idx][i_comp];
                     contrib += p.lambda * grads[i_node][q_idx][i_comp] * grads[j_node][q_idx][j_comp];
                  }

                  matrix[i_node][i_comp][j_node][j_comp] += q.weights[q_idx] * contrib;
               }
            }
         }
      }
   }

}


void assemble_matrix_rhs(p4est_lnodes_t *lnodes, BddcmlMesh *mesh, BddcmlFemSpace *femsp,
                         SparseMatrix *matrix, RealArray *rhss, RhsPtr rhs_ptr, Parameters params)
{
   Element element;
   real i_coeffs, j_coeffs;
   int n_components = mesh->subdomain_dims->n_node_dofs;
   vector<vector<vector<vector<real> > > > element_matrix(P4EST_CHILDREN, vector<vector<vector<real> > >(
                                           n_components, vector<vector<real> >(
                                           P4EST_CHILDREN, vector<real>(
                                           n_components,0.0))));
   vector<vector<real> > element_rhs(P4EST_CHILDREN, vector<real>(n_components, 0.0));
   p4est_locidx_t i_indep_nodes[4], j_indep_nodes[4];

   int element_offset = 0;
   for(int elem_idx = 0; elem_idx < mesh->subdomain_dims->n_elems; elem_idx++) {
      assert(mesh->num_nodes_of_elem.val[elem_idx] == P4EST_CHILDREN);

      mesh->get_element(elem_idx, &element);

      if(femsp->physicsType == PhysicsType::LAPLACE)
         assemble_local_laplace(&element, element_matrix, element_rhs, rhs_ptr);
      else if(femsp->physicsType == PhysicsType::ELASTICITY)
         assemble_local_elasticity(&element, element_matrix, element_rhs, rhs_ptr, params);
      else
         assert(0);

      //print_matrix_rhs(element_matrix, element_rhs, n_components);

      for(int i_node_loc = 0; i_node_loc < P4EST_CHILDREN; i_node_loc++) {
         int i_node = mesh->elem_node_indices.val[element_offset + i_node_loc];
         int i_nindep = independent_nodes(lnodes, elem_idx, i_node_loc, i_indep_nodes, &i_coeffs);
         assert((i_nindep != 1) || (i_node == i_indep_nodes[0]));

         for(int i_indep_nodes_idx = 0; i_indep_nodes_idx < i_nindep; i_indep_nodes_idx++)
         {
            int i_indep_node = i_indep_nodes[i_indep_nodes_idx];

            for(int i_comp = 0; i_comp < n_components; i_comp++)
            {
               int i_dof = femsp->node_num_dofs.val[i_indep_node] * i_indep_node + i_comp;
               for(int j_node_loc = 0; j_node_loc < P4EST_CHILDREN; j_node_loc++) {
                  //todo: dofs should be taken from femsp!
                  int j_node = mesh->elem_node_indices.val[element_offset + j_node_loc];
                  int j_nindep = independent_nodes(lnodes, elem_idx, j_node_loc, j_indep_nodes, &j_coeffs);
                  assert((j_nindep != 1) || (j_node == j_indep_nodes[0]));

                  for(int j_indep_nodes_idx = 0; j_indep_nodes_idx < j_nindep; j_indep_nodes_idx++)
                  {
                     int j_indep_node = j_indep_nodes[j_indep_nodes_idx];

                     for(int j_comp = 0; j_comp < n_components; j_comp++)
                     {
                        int j_dof = femsp->node_num_dofs.val[j_indep_node] * j_indep_node + j_comp;
                        double matrix_value = i_coeffs * j_coeffs * element_matrix[i_node_loc][i_comp][j_node_loc][j_comp];
                        add_matrix_entry(matrix, i_dof, j_dof, matrix_value);
//                        printf("adding entry loc (%d, %d), nodes orig (%d, %d), nodes indep (%d, %d), dofs (%d, %d), coefs (%3.2lf, %3.2lf), number indep (%d, %d), locstiff %lf, value %lf\n",
//                               i_node_loc, j_node_loc, i_node, j_node, i_indep_node, j_indep_node, i_dof, j_dof, i_coeffs, j_coeffs, i_nindep, j_nindep, element_matrix[i_node_loc][i_comp][j_node_loc][j_comp], matrix_value);
                     }
                  }
               }

               //double rhs_value = i_coeffs * 1./(real)P4EST_CHILDREN * elem_volume * 1;
               double rhs_value = i_coeffs * element_rhs[i_node_loc][i_comp];
               rhss->val[i_dof] += rhs_value;

            }
         }
      }
      element_offset += P4EST_CHILDREN;
   }

}

