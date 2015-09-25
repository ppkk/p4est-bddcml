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
#include "p4est_common.h"
#include "assemble.h"
#include "mesh.h"
#include "femspace.h"

using namespace std;

void prepare_dimmensions(p4est_t *p4est, p4est_lnodes_t *lnodes, PhysicsType physicsType,
                         BddcmlDimensions *subdomain_dims, BddcmlDimensions *global_dims,
                         sc_MPI_Comm mpicomm)
{
#ifndef P4_TO_P8
   init_dimmensions(subdomain_dims, 2, physicsType);
   init_dimmensions(global_dims, 2, physicsType);
#else
   init_dimmensions(subdomain_dims, 3, physicsType);
   init_dimmensions(global_dims, 3, physicsType);
#endif
   subdomain_dims->n_nodes = lnodes->num_local_nodes;
   subdomain_dims->n_dofs  = lnodes->num_local_nodes * subdomain_dims->n_node_dofs;
   subdomain_dims->n_elems = lnodes->num_local_elements;
   printf("proc %d, elems %d, nodes %d\n", mpi_rank, subdomain_dims->n_elems, subdomain_dims->n_nodes);

   int global_num_nodes;
   if(mpi_rank == mpi_size - 1)
   {
      global_num_nodes = lnodes->global_offset + lnodes->owned_count;
   }
   sc_MPI_Bcast(&global_num_nodes, 1, MPI_INT, mpi_size - 1, mpicomm);

   global_dims->n_nodes = global_num_nodes;
   global_dims->n_dofs = global_num_nodes * global_dims->n_node_dofs;
   global_dims->n_elems = p4est->global_num_quadrants;
}


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

vector<double> rhs(vector<double>)
{
   return {1};
}

struct Quadrature
{
   vector<double> weights;
   vector<vector<double> > coords;

   Quadrature(int dimmension, double element_length)
   {
      double scale = element_length / 2.;
      if(dimmension == 1)
      {
         weights = {5./9. * scale, 8./9. * scale, 5./9. * scale};
         coords = {vector<double>({-sqrt(3./5.)}), vector<double>({0.}), vector<double>({sqrt(3./5.)})};
      }
      else if(dimmension == 2)
      {
         Quadrature q1(1, element_length);
         product(q1, q1);
      }
      else if(dimmension == 3)
      {
         Quadrature q1(1, element_length), q2(2, element_length);
         product(q1, q2);
      }
      else
         assert(0);
   }

   void product(Quadrature quad1, Quadrature quad2)
   {
      for(unsigned int i1 = 0; i1 < quad1.weights.size(); i1++)
      {
         for(unsigned int i2 = 0; i2 < quad2.weights.size(); i2++)
         {
            weights.push_back(quad1.weights[i1] * quad2.weights[i2]);
            vector<double>product_coords;
            product_coords.insert(product_coords.end(), quad1.coords[i1].begin(), quad1.coords[i1].end());
            product_coords.insert(product_coords.end(), quad2.coords[i2].begin(), quad2.coords[i2].end());
            coords.push_back(product_coords);
         }
      }
   }

   void print()
   {
      assert(weights.size() == coords.size());
      double sum = 0.0;
      for(unsigned int i = 0; i < weights.size(); i++)
      {
         cout << "(";
         for(unsigned int j = 0; j < coords[i].size(); j++)
         {
            cout << coords[i][j] << ", ";
         }
         cout << "), " << weights[i] << endl;
         sum += weights[i];
      }
      cout << "sum of weights " << sum << endl;
   }

};

void ref_value_1D(int loc_id_1d, double x, double elem_len, double& value, double& der)
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
         ref_value_1D(z_id_1D, q.coords[q_idx][2], element_length, value_z, der_z);
#endif

         double value = value_x * value_y;
         double grad_1 = der_x * value_y;
         double grad_2 = value_x * der_y;

#ifdef P4_TO_P8
         value *= value_z;
         grad_1 *= value_z;
         grad_2 *= value_z;
         double grad_3 = value_x * value_y * der_z;
         gradients[node].push_back(vector<double>({grad_1, grad_2, grad_3}));
#else
         gradients[node].push_back(vector<double>({grad_1, grad_2}));
#endif
         values[node].push_back(value);
      }
   }
}

void assemble_local_laplace(double element_size, vector<vector<vector<vector<real> > > > &stiffness,
                    vector<vector<real> > &rhs, vector<double> (*rhs_ptr)(vector<double>))
{
   Quadrature q(P4EST_DIM, element_size);

   vector<vector<double> > values;
   vector<vector<vector<double> > > gradients;

   prepare_transformed_values(q, element_size, values, gradients);

   for(int i_node = 0; i_node < P4EST_CHILDREN; i_node++)
   {
      rhs[i_node][0] = 0.0;
      for(int j_node = 0; j_node < P4EST_CHILDREN; j_node++)
         stiffness[i_node][0][j_node][0] = 0.0;
   }

   for(unsigned int q_idx = 0; q_idx < q.weights.size(); q_idx++)
   {
      for(int i_node = 0; i_node < P4EST_CHILDREN; i_node++)
      {

         rhs[i_node][0] += q.weights[q_idx] * values[i_node][q_idx] * rhs_ptr(q.coords[q_idx])[0];

         for(int j_node = 0; j_node < P4EST_CHILDREN; j_node++)
         {
            for(int idx_dim = 0; idx_dim < P4EST_DIM; idx_dim++)
            {
               stiffness[i_node][0][j_node][0] += q.weights[q_idx] * gradients[i_node][q_idx][idx_dim] * gradients[j_node][q_idx][idx_dim];
            }
         }
      }
   }

}


void assemble_matrix_rhs(p4est_lnodes_t *lnodes, BddcmlMesh *mesh, BddcmlFemSpace *femsp,
                         SparseMatrix *matrix, RealArray *rhss)
{
   real i_coeffs, j_coeffs;
   int n_components = mesh->subdomain_dims->n_node_dofs;
   assert(n_components == 1);
   vector<vector<vector<vector<real> > > > element_matrix(P4EST_CHILDREN, vector<vector<vector<real> > >(
                                           n_components, vector<vector<real> >(
                                           P4EST_CHILDREN, vector<real>(
                                           n_components,0.0))));
   vector<vector<real> > element_rhs(P4EST_CHILDREN, vector<real>(n_components, 0.0));
   p4est_locidx_t i_indep_nodes[4], j_indep_nodes[4];

   int element_offset = 0;
   for(int elem_idx = 0; elem_idx < mesh->subdomain_dims->n_elems; elem_idx++) {
      assert(mesh->num_nodes_of_elem.val[elem_idx] == P4EST_CHILDREN);

      // TODO: elem_size and elem_volume is correct only when the mesh is obtained by refinements
      // from a UNIT SQUARE/CUBE
      // TODO: ONLY FROM ****UNIT****
      real element_size = mesh->element_lengths.val[elem_idx];

      assemble_local_laplace(element_size, element_matrix, element_rhs, &rhs);

      for(int i_node = 0; i_node < P4EST_CHILDREN; i_node++) {
         int idof = mesh->elem_node_indices.val[element_offset + i_node];
         int i_nindep = independent_nodes(lnodes, elem_idx, i_node, i_indep_nodes, &i_coeffs);
         assert((i_nindep != 1) || (idof == i_indep_nodes[0]));

         for(int i_indep_nodes_idx = 0; i_indep_nodes_idx < i_nindep; i_indep_nodes_idx++)
         {
            int i_indep_node = i_indep_nodes[i_indep_nodes_idx];

            for(int i_comp = 0; i_comp < n_components; i_comp++)
            {
               for(int j_node = 0; j_node < P4EST_CHILDREN; j_node++) {
                  //todo: dofs should be taken from femsp!
                  int jdof = mesh->elem_node_indices.val[element_offset + j_node];
                  assert(femsp->node_num_dofs.val[jdof] == 1);
                  int j_nindep = independent_nodes(lnodes, elem_idx, j_node, j_indep_nodes, &j_coeffs);
                  assert((j_nindep != 1) || (jdof == j_indep_nodes[0]));

                  for(int j_indep_nodes_idx = 0; j_indep_nodes_idx < j_nindep; j_indep_nodes_idx++)
                  {
                     int j_indep_node = j_indep_nodes[j_indep_nodes_idx];

                     for(int j_comp = 0; j_comp < n_components; j_comp++)
                     {
                        double matrix_value = i_coeffs * j_coeffs * element_matrix[i_node][i_comp][j_node][j_comp];
                        add_matrix_entry(matrix, i_indep_node, j_indep_node, matrix_value);
//                        printf("adding entry loc (%d, %d), nodes, (%d, %d), coefs (%3.2lf, %3.2lf), number indep (%d, %d), locstiff %lf, value %lf\n",
//                               j_node, i_node, j_indep_node, i_indep_node, j_coeffs, i_coeffs, j_nindep, i_nindep, element_matrix[i_node][i_comp][j_node][j_comp], matrix_value);
                     }
                  }
               }

               // TODO: integrate properly
               // TODO: elem_volume is correct only when the mesh is obtained by refinements
               // from a UNIT SQUARE/CUBE
               // TODO: ONLY FROM ****UNIT****


               //double rhs_value = i_coeffs * 1./(real)P4EST_CHILDREN * elem_volume * 1;
               double rhs_value = i_coeffs * element_rhs[i_node][i_comp];
               rhss->val[i_indep_node] += rhs_value;

            }
         }
      }
      element_offset += P4EST_CHILDREN;
   }

}

