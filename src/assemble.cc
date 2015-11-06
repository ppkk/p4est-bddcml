#include <string.h>
#include <stdlib.h>
#include <stdbool.h>

#include "assemble.h"
#include "element.h"
#include "quadrature.h"
#include "integration_cell.h"
#include "p4est/my_p4est_interface.h"
#include "shapefun.h"
#include "local_matrix.h"

using namespace std;

DiscreteSystem::DiscreteSystem(const ProblemDimensions &problem_dims, MatrixType matrix_type) {
   allocate_real_array(problem_dims.n_subdom_dofs, &rhss);
   zero_real_array(&rhss);

   int ndof_per_element = Def::d()->num_element_nodes * problem_dims.n_node_dofs;
   // how much space the upper triangle of the element matrix occupies
   int lelm = ndof_per_element * (ndof_per_element + 1) / 2;

   // todo: do it properly
   const int extra_space_for_hanging_nodes = 4 * (matrix_type == MatrixType::GENERAL ? 2 : 1);
   matrix.allocate(extra_space_for_hanging_nodes * problem_dims.n_subdom_elems * lelm, matrix_type);
   matrix.zero();
}

void DiscreteSystem::free() {
   free_real_array(&rhss);
   matrix.free_matrix();
}

//void print_complete_matrix_rhs(const BddcmlFemSpace &femsp, const ProblemDimensions &global_dims,
//                               const SparseMatrix &matrix, const RealArray &rhss, MPI_Comm mpicomm) {
//   real* compl_rhs = (real*) malloc(global_dims.n_dofs * sizeof(real));
//   real* compl_mat = (real*) malloc(global_dims.n_dofs * global_dims.n_dofs * sizeof(real));

//   if(!compl_rhs || !compl_mat) {
//      printf("Unable to allocate dense matrix. Exiting print_complete_matrix_rhs....\n");
//      return;
//   }

//   memset(compl_rhs, 0, global_dims.n_dofs * sizeof(real));
//   memset(compl_mat, 0, global_dims.n_dofs * global_dims.n_dofs * sizeof(real));

//   for(int loc_dof = 0; loc_dof < femsp.subdomain_dims->n_dofs; loc_dof++) {
//      int glob_dof = femsp.dofs_global_map.val[loc_dof];
//      compl_rhs[glob_dof] = rhss.val[loc_dof];
//   }

//   for(int idx = 0; idx < matrix.nnz; idx++) {
//      int glob_dof_j = femsp.dofs_global_map.val[matrix.j[idx]];
//      int glob_dof_i = femsp.dofs_global_map.val[matrix.i[idx]];
//      compl_mat[global_dims.n_dofs * glob_dof_j + glob_dof_i] += matrix.val[idx];

//      // adding also symmetric counterpart
//      if(matrix.type != MatrixType::GENERAL)
//         compl_mat[global_dims.n_dofs * glob_dof_i + glob_dof_j] += matrix.val[idx];
//   }

//   if(mpi_rank == 0) {
//      MPI_Reduce(MPI_IN_PLACE, compl_rhs, global_dims.n_dofs, MPI_DOUBLE, MPI_SUM, 0, mpicomm);
//      MPI_Reduce(MPI_IN_PLACE, compl_mat, global_dims.n_dofs * global_dims.n_dofs, MPI_DOUBLE, MPI_SUM, 0, mpicomm);
//   }
//   else {
//      MPI_Reduce(compl_rhs, compl_rhs, global_dims.n_dofs, MPI_DOUBLE, MPI_SUM, 0, mpicomm);
//      MPI_Reduce(compl_mat, compl_mat, global_dims.n_dofs * global_dims.n_dofs, MPI_DOUBLE, MPI_SUM, 0, mpicomm);
//   }

//   if(mpi_rank == 0) {
//      FILE *f = fopen("discr.txt", "w");
//      for(int dof_j = 0; dof_j < global_dims.n_dofs; dof_j++) {
//         fprintf(f, "dof %d, rhs: %g\n", dof_j, compl_rhs[dof_j]);
//         for(int dof_i = 0; dof_i < global_dims.n_dofs; dof_i++) {
//            fprintf(f, "%g\n", compl_mat[dof_j*global_dims.n_dofs + dof_i]);
//         }
//      }
//      fclose(f);
//   }

//   free(compl_rhs);
//   free(compl_mat);
//}


void assemble_local_laplace(const IntegrationCell &cell, const ReferenceElement &ref_elem, const GaussQuadrature &q, RhsPtr rhs_ptr,
                            LocalMatrix *matrix, LocalVector *rhs) {
   vector<vector<double> > vals;
   vector<vector<vector<double> > > grads;

   ref_elem.calc_transformed_values(q, cell.size, &vals, &grads);

   Quadrature q_transformed(Def::d()->num_dim);
   q.transform_to_physical(cell, &q_transformed);

   matrix->clear();
   rhs->clear();

   for(unsigned int q_idx = 0; q_idx < q_transformed.np(); q_idx++) {
      for(int i_node = 0; i_node < Def::d()->num_element_nodes; i_node++) {

         rhs->comps[0].vec[i_node] += q_transformed.weights[q_idx] * vals[i_node][q_idx] * rhs_ptr(q_transformed.coords[q_idx])[0];

         for(int j_node = 0; j_node < Def::d()->num_element_nodes; j_node++) {
            for(int idx_dim = 0; idx_dim < Def::d()->num_dim; idx_dim++) {
               matrix->comps[0][0].mat[i_node][j_node] += q_transformed.weights[q_idx] * grads[i_node][q_idx][idx_dim] * grads[j_node][q_idx][idx_dim];
            }
         }
      }
   }

}

void assemble_local_elasticity(const IntegrationCell &integ_cell, const ReferenceElement &ref_elem, const GaussQuadrature &q, RhsPtr rhs_ptr, Parameters p,
                               LocalMatrix *matrix, LocalVector *rhs) {
   int num_comp = Def::d()->num_dim;
   vector<vector<double> > vals; //[node][quadrature_pt]
   vector<vector<vector<double> > > grads;//[node][quadrature_pt][component]

   ref_elem.calc_transformed_values(q, integ_cell.size, &vals, &grads);

   Quadrature q_transformed(Def::d()->num_dim);
   q.transform_to_physical(integ_cell, &q_transformed);

   matrix->clear();
   rhs->clear();

   for(unsigned int q_idx = 0; q_idx < q_transformed.np(); q_idx++) {
      for(int i_node = 0; i_node < Def::d()->num_element_nodes; i_node++) {
         for(int i_comp = 0; i_comp < num_comp; i_comp++) {
            rhs->comps[i_comp].vec[i_node] += q_transformed.weights[q_idx] * vals[i_node][q_idx] * rhs_ptr(q_transformed.coords[q_idx])[i_comp];

            for(int j_node = 0; j_node < Def::d()->num_element_nodes; j_node++) {
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

void DiscreteSystem::assemble(const P4estClass &p4est, const IntegrationMesh &integration_mesh,
                         const NodalElementMesh &nodal_mesh, const ProblemDimensions &problem_dims,
                         RhsPtr rhs_ptr, Parameters params) {
   assert(integration_mesh.num_elements() == problem_dims.n_subdom_elems);

   GaussQuadrature q(Def::d()->num_dim, 2*Def::d()->order);
   ReferenceElement ref_elem(Def::d()->num_dim, Def::d()->order);

   HangingInfo hanging_info(p4est);
   int n_components = problem_dims.n_node_dofs;

   LocalMatrix element_matrix_nohang(n_components, Def::d()->num_element_nodes), element_matrix(n_components, Def::d()->num_element_nodes);
   LocalVector element_rhs_nohang(n_components, Def::d()->num_element_nodes), element_rhs(n_components, Def::d()->num_element_nodes);
   //vector<vector<real> > element_rhs(Def::d()->num_element_nodes, vector<real>(n_components, 0.0));

   int element_offset = 0;
   for(int elem_idx = 0; elem_idx < problem_dims.n_subdom_elems; elem_idx++) {
      const IntegrationCell& cell = integration_mesh.cells[elem_idx];
      const NodalElement& nodal_element = nodal_mesh.elements[elem_idx];

      //mesh.get_element(elem_idx, &element);

      if(nodal_mesh.physics_type == PhysicsType::LAPLACE)
         assemble_local_laplace(cell, ref_elem, q, rhs_ptr, &element_matrix_nohang, &element_rhs_nohang);
      else if(nodal_mesh.physics_type == PhysicsType::ELASTICITY)
         assemble_local_elasticity(cell, ref_elem, q, rhs_ptr, params, &element_matrix_nohang, &element_rhs_nohang);
      else
         assert(0);

      hanging_info.apply_constraints(elem_idx, cell, element_matrix_nohang, &element_matrix);
      hanging_info.apply_constraints(elem_idx, cell, element_rhs_nohang, &element_rhs);

//      element_matrix.print();
//      element_rhs.print();

      for(int i_node_loc = 0; i_node_loc < Def::d()->num_element_nodes; i_node_loc++) {
         for(int i_comp = 0; i_comp < n_components; i_comp++) {
            int i_dof = nodal_element.components[i_comp].dofs[i_node_loc];

            for(int j_node_loc = 0; j_node_loc < Def::d()->num_element_nodes; j_node_loc++) {
               for(int j_comp = 0; j_comp < n_components; j_comp++) {
                  int j_dof = nodal_element.components[j_comp].dofs[j_node_loc];

                  double matrix_value = element_matrix.comps[i_comp][j_comp].mat[i_node_loc][j_node_loc];
                  matrix.add_entry(i_dof, j_dof, matrix_value);
                  //                        printf("adding entry loc (%d, %d), nodes orig (%d, %d), nodes indep (%d, %d), dofs (%d, %d), coefs (%3.2lf, %3.2lf), number indep (%d, %d), locstiff %lf, value %lf\n",
                  //                               i_node_loc, j_node_loc, i_node, j_node, i_indep_node, j_indep_node, i_dof, j_dof, i_coeffs, j_coeffs, i_nindep, j_nindep, element_matrix[i_node_loc][i_comp][j_node_loc][j_comp], matrix_value);
               }
            }

            //double rhs_value = i_coeffs * 1./(real)Def::d()->num_children * elem_volume * 1;
            double rhs_value = element_rhs.comps[i_comp].vec[i_node_loc];
            rhss.val[i_dof] += rhs_value;
         }
      }
      element_offset += Def::d()->num_element_nodes;
   }

}

