#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <assert.h>

#include "assemble.h"
#include "bddcml/bddcml_mesh.h"
#include "bddcml/bddcml_femspace.h"
#include "quadrature.h"
#include "geometry_mesh.h"
#include "p4est/my_p4est_interface.h"


using namespace std;

void print_complete_matrix_rhs(const BddcmlFemSpace &femsp, const BddcmlDimensions &global_dims,
                               const SparseMatrix &matrix, const RealArray &rhss, MPI_Comm mpicomm) {
   real* compl_rhs = (real*) malloc(global_dims.n_dofs * sizeof(real));
   real* compl_mat = (real*) malloc(global_dims.n_dofs * global_dims.n_dofs * sizeof(real));

   if(!compl_rhs || !compl_mat) {
      printf("Unable to allocate dense matrix. Exiting print_complete_matrix_rhs....\n");
      return;
   }

   memset(compl_rhs, 0, global_dims.n_dofs * sizeof(real));
   memset(compl_mat, 0, global_dims.n_dofs * global_dims.n_dofs * sizeof(real));

   for(int loc_dof = 0; loc_dof < femsp.subdomain_dims->n_dofs; loc_dof++) {
      int glob_dof = femsp.dofs_global_map.val[loc_dof];
      compl_rhs[glob_dof] = rhss.val[loc_dof];
   }

   for(int idx = 0; idx < matrix.nnz; idx++) {
      int glob_dof_j = femsp.dofs_global_map.val[matrix.j[idx]];
      int glob_dof_i = femsp.dofs_global_map.val[matrix.i[idx]];
      compl_mat[global_dims.n_dofs * glob_dof_j + glob_dof_i] += matrix.val[idx];

      // adding also symmetric counterpart
      if(matrix.type != MatrixType::GENERAL)
         compl_mat[global_dims.n_dofs * glob_dof_i + glob_dof_j] += matrix.val[idx];
   }

   if(mpi_rank == 0) {
      MPI_Reduce(MPI_IN_PLACE, compl_rhs, global_dims.n_dofs, MPI_DOUBLE, MPI_SUM, 0, mpicomm);
      MPI_Reduce(MPI_IN_PLACE, compl_mat, global_dims.n_dofs * global_dims.n_dofs, MPI_DOUBLE, MPI_SUM, 0, mpicomm);
   }
   else {
      MPI_Reduce(compl_rhs, compl_rhs, global_dims.n_dofs, MPI_DOUBLE, MPI_SUM, 0, mpicomm);
      MPI_Reduce(compl_mat, compl_mat, global_dims.n_dofs * global_dims.n_dofs, MPI_DOUBLE, MPI_SUM, 0, mpicomm);
   }

   if(mpi_rank == 0) {
      FILE *f = fopen("discr.txt", "w");
      for(int dof_j = 0; dof_j < global_dims.n_dofs; dof_j++) {
         fprintf(f, "dof %d, rhs: %g\n", dof_j, compl_rhs[dof_j]);
         for(int dof_i = 0; dof_i < global_dims.n_dofs; dof_i++) {
            fprintf(f, "%g\n", compl_mat[dof_j*global_dims.n_dofs + dof_i]);
         }
      }
      fclose(f);
   }

   free(compl_rhs);
   free(compl_mat);
}


void ref_value_1D(int loc_id_1d, double x, double elem_len, double & value, double & der) {
   if(loc_id_1d == 0) {
      value = (1-x)/2.;
      der = -1/elem_len;
   }
   else if(loc_id_1d == 1) {
      value = (1+x)/2.;
      der = 1/elem_len;
   }
   else {
      assert(0);
   }
}

void prepare_transformed_values(const Quadrature &q, double element_length,
                                vector<vector<double> > &values, vector<vector<vector<double> > > &gradients) {
   values = vector<vector<double> >(P4estClass::children);
   gradients = vector<vector<vector<double> > >(P4estClass::children);

   for(int node = 0; node < P4estClass::children; node++) {
      int x_id_1D = node % 2;
      int y_id_1D = (node % 4) / 2;
      int z_id_1D = node / 4;

      for(unsigned int q_idx = 0; q_idx < q.np(); q_idx++) {
         double value_x, der_x, value_y, der_y, value_z, der_z;
         ref_value_1D(x_id_1D, q.coords[q_idx][0], element_length, value_x, der_x);
         ref_value_1D(y_id_1D, q.coords[q_idx][1], element_length, value_y, der_y);

         if(P4estClass::num_dim == 3) {
            ref_value_1D(z_id_1D, (q.coords[q_idx])[2], element_length, value_z, der_z);
         }

         double value = value_x * value_y;
         double grad_1 = der_x * value_y;
         double grad_2 = value_x * der_y;

         if(P4estClass::num_dim == 3) {
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
         }
         else {
            vector<double> aux1;
            aux1.push_back(grad_1);
            aux1.push_back(grad_2);
            //gradients[node].push_back(vector<double>({grad_1, grad_2}));
            gradients[node].push_back(aux1);
         }
         values[node].push_back(value);
      }
   }
}

void zero_matrix_rhs(LocalMatrix *matrix, vector<vector<real> > *rhs) {
   for(int i_node = 0; i_node < P4estClass::children; i_node++) {
      for(int i_comp = 0; i_comp < matrix->ncomponents; i_comp++) {
         (*rhs)[i_node][i_comp] = 0.0;
         for(int j_node = 0; j_node < P4estClass::children; j_node++) {
            for(int j_comp = 0; j_comp < matrix->ncomponents; j_comp++) {
               matrix->comps[i_comp][j_comp].mat[i_node][j_node] = 0.0;
            }
         }
      }
   }
}

void print_matrix_rhs(const vector<vector<vector<vector<real> > > > &matrix, const vector<vector<real> > &rhs, int num_comp) {
   cout << "LOCAL MATRIX:" << endl;
   for(int i_node = 0; i_node < P4estClass::children; i_node++) {
      for(int i_comp = 0; i_comp < num_comp; i_comp++) {
         for(int j_node = 0; j_node < P4estClass::children; j_node++) {
            for(int j_comp = 0; j_comp < num_comp; j_comp++) {
               cout << "nodes (" << i_node << ", " << j_node << "), comp (" << i_comp << ", " << j_comp << ") -> " <<
                       matrix[i_node][i_comp][j_node][j_comp] << endl;
            }
         }
      }
   }
   cout << "END LOCAL MATRIX:" << endl;
}

void assemble_local_laplace(const Element &element, RhsPtr rhs_ptr, LocalMatrix *matrix, vector<vector<real> > *rhs) {
   GaussQuadrature q(P4estClass::num_dim, QUAD_ORDER);

   vector<vector<double> > vals;
   vector<vector<vector<double> > > grads;

   prepare_transformed_values(q, element.size, vals, grads);

   Quadrature q_transformed(P4estClass::num_dim);
   q.transform_to_physical(element, &q_transformed);

   zero_matrix_rhs(matrix, rhs);

   for(unsigned int q_idx = 0; q_idx < q_transformed.np(); q_idx++) {
      for(int i_node = 0; i_node < P4estClass::children; i_node++) {

         (*rhs)[i_node][0] += q_transformed.weights[q_idx] * vals[i_node][q_idx] * rhs_ptr(q_transformed.coords[q_idx])[0];

         for(int j_node = 0; j_node < P4estClass::children; j_node++) {
            for(int idx_dim = 0; idx_dim < P4estClass::num_dim; idx_dim++) {
               matrix->comps[0][0].mat[i_node][j_node] += q_transformed.weights[q_idx] * grads[i_node][q_idx][idx_dim] * grads[j_node][q_idx][idx_dim];
            }
         }
      }
   }

}

void assemble_local_elasticity(const Element &element, RhsPtr rhs_ptr, Parameters p,
                               LocalMatrix *matrix, vector<vector<real> > *rhs) {
   GaussQuadrature q(P4estClass::num_dim, QUAD_ORDER);

   int num_comp = P4estClass::num_dim;
   vector<vector<double> > vals; //[node][quadrature_pt]
   vector<vector<vector<double> > > grads;//[node][quadrature_pt][component]

   prepare_transformed_values(q, element.size, vals, grads);

   Quadrature q_transformed(P4estClass::num_dim);
   q.transform_to_physical(element, &q_transformed);

   zero_matrix_rhs(matrix, rhs);

   for(unsigned int q_idx = 0; q_idx < q_transformed.np(); q_idx++) {
      for(int i_node = 0; i_node < P4estClass::children; i_node++) {
         for(int i_comp = 0; i_comp < num_comp; i_comp++) {
            (*rhs)[i_node][i_comp] += q_transformed.weights[q_idx] * vals[i_node][q_idx] * rhs_ptr(q_transformed.coords[q_idx])[i_comp];

            for(int j_node = 0; j_node < P4estClass::children; j_node++) {
               for(int j_comp = 0; j_comp < num_comp; j_comp++) {
                  double contrib = 0.0;

                  if(i_comp == j_comp) {
                     contrib += (p.lambda + 2 * p.mu) * grads[i_node][q_idx][i_comp] * grads[j_node][q_idx][j_comp];
                     for(int other_component = 0; other_component < num_comp; other_component++) {
                        if(other_component == i_comp)
                           continue;
                        contrib += p.mu * grads[i_node][q_idx][other_component] * grads[j_node][q_idx][other_component];

                     }
                  }
                  else {
                     contrib += p.mu     * grads[i_node][q_idx][j_comp] * grads[j_node][q_idx][i_comp];
                     contrib += p.lambda * grads[i_node][q_idx][i_comp] * grads[j_node][q_idx][j_comp];
                  }

                  matrix->comps[i_comp][j_comp].mat[i_node][j_node] += q_transformed.weights[q_idx] * contrib;
               }
            }
         }
      }
   }

}

LocalMatrixComponent::LocalMatrixComponent(int ndofs) : ndofs(ndofs){
   mat.resize(ndofs, vector<real>(ndofs, 0.0));
}

LocalMatrix::LocalMatrix(int ncomponents, int ndofs) : ncomponents(ncomponents), ndofs(ndofs) {
   comps.resize(ncomponents, vector<LocalMatrixComponent>(ncomponents, LocalMatrixComponent(ndofs)));
}


void assemble_matrix_rhs(const P4estClass &p4est, const GeometryMesh &geometry_mesh, const BddcmlMesh &bddcml_mesh, const BddcmlFemSpace &femsp,
                         SparseMatrix *matrix, RealArray *rhss, RhsPtr rhs_ptr, Parameters params) {

   assert(geometry_mesh.num_elements() == bddcml_mesh.subdomain_dims->n_elems);
   real i_coeffs, j_coeffs;
   int n_components = bddcml_mesh.subdomain_dims->n_node_dofs;

   LocalMatrix element_matrix(n_components, P4estClass::children);
   vector<vector<real> > element_rhs(P4estClass::children, vector<real>(n_components, 0.0));
   p4est_locidx_t i_indep_nodes[4], j_indep_nodes[4];

   int element_offset = 0;
   for(int elem_idx = 0; elem_idx < bddcml_mesh.subdomain_dims->n_elems; elem_idx++) {
      assert(bddcml_mesh.num_nodes_of_elem.val[elem_idx] == P4estClass::children);

      //mesh.get_element(elem_idx, &element);

      if(femsp.physicsType == PhysicsType::LAPLACE)
         assemble_local_laplace(geometry_mesh.elements[elem_idx], rhs_ptr, &element_matrix, &element_rhs);
      else if(femsp.physicsType == PhysicsType::ELASTICITY)
         assemble_local_elasticity(geometry_mesh.elements[elem_idx], rhs_ptr, params, &element_matrix, &element_rhs);
      else
         assert(0);

      //print_matrix_rhs(element_matrix, element_rhs, n_components);

      for(int i_node_loc = 0; i_node_loc < P4estClass::children; i_node_loc++) {
         int i_node = bddcml_mesh.elem_node_indices.val[element_offset + i_node_loc];
         int i_nindep = p4est.independent_nodes(elem_idx, i_node_loc, i_indep_nodes, &i_coeffs);
         assert((i_nindep != 1) || (i_node == i_indep_nodes[0]));

         for(int i_indep_nodes_idx = 0; i_indep_nodes_idx < i_nindep; i_indep_nodes_idx++) {
            int i_indep_node = i_indep_nodes[i_indep_nodes_idx];

            for(int i_comp = 0; i_comp < n_components; i_comp++) {
               int i_dof = femsp.node_num_dofs.val[i_indep_node] * i_indep_node + i_comp;
               for(int j_node_loc = 0; j_node_loc < P4estClass::children; j_node_loc++) {
                  //todo: dofs should be taken from femsp!
                  int j_node = bddcml_mesh.elem_node_indices.val[element_offset + j_node_loc];
                  int j_nindep = p4est.independent_nodes(elem_idx, j_node_loc, j_indep_nodes, &j_coeffs);
                  assert((j_nindep != 1) || (j_node == j_indep_nodes[0]));

                  for(int j_indep_nodes_idx = 0; j_indep_nodes_idx < j_nindep; j_indep_nodes_idx++) {
                     int j_indep_node = j_indep_nodes[j_indep_nodes_idx];

                     for(int j_comp = 0; j_comp < n_components; j_comp++) {
                        int j_dof = femsp.node_num_dofs.val[j_indep_node] * j_indep_node + j_comp;
                        double matrix_value = i_coeffs * j_coeffs * element_matrix.comps[i_comp][j_comp].mat[i_node_loc][j_node_loc];
                        add_matrix_entry(matrix, i_dof, j_dof, matrix_value);
//                        printf("adding entry loc (%d, %d), nodes orig (%d, %d), nodes indep (%d, %d), dofs (%d, %d), coefs (%3.2lf, %3.2lf), number indep (%d, %d), locstiff %lf, value %lf\n",
//                               i_node_loc, j_node_loc, i_node, j_node, i_indep_node, j_indep_node, i_dof, j_dof, i_coeffs, j_coeffs, i_nindep, j_nindep, element_matrix[i_node_loc][i_comp][j_node_loc][j_comp], matrix_value);
                     }
                  }
               }

               //double rhs_value = i_coeffs * 1./(real)P4estClass::children * elem_volume * 1;
               double rhs_value = i_coeffs * element_rhs[i_node_loc][i_comp];
               rhss->val[i_dof] += rhs_value;

            }
         }
      }
      element_offset += P4estClass::children;
   }

}


void shape_fun(int order, int idx) {
   int n_nodes = order + 1;
   assert((idx >= 0) && (idx < n_nodes));

   double distance = 2./order;
   vector<double> nodes;
   double node = -1.;
   for(int i = 0; i < n_nodes; i++) {
      nodes.push_back(node);
      node += distance;
   }

   double denominator = 1.;
   for(int i = 0; i < n_nodes; i++) {
      if(i != idx) {
         denominator *= (nodes[idx] - nodes[i]);
      }
   }

}
