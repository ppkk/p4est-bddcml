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

#include "helpers.h"
#include "bddcml_structs.h"
#include "p4est_common.h"
#include "p4est_bddcml_interaction.h"



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

void prepare_subdomain_mesh(p4est_t *p4est, p4est_lnodes_t *lnodes, BddcmlDimensions *subdomain_dims, BddcmlMesh *mesh)
{
   double              vxyz[3];  /* We embed the 2D vertices into 3D space. */
   p4est_quadrant_t   sp, node;

   init_mesh(subdomain_dims, mesh);
   for(int lnode = 0; lnode < lnodes->num_local_nodes; lnode++)
   {
      mesh->node_global_map.val[lnode] = node_loc_to_glob(lnodes, lnode);
   }
   //p4est->global_first_quadrant
   /* Loop over local quadrants to apply the element matrices. */
   p4est_locidx_t quad_idx = 0;
   p4est_quadrant_t *quad;

   for_all_quads(p4est, quad_idx, quad)
   {
      // element local to global mapping -- obtained by adding the offset from p4est
      mesh->elem_global_map.val[quad_idx] = (int)p4est->global_first_quadrant[p4est->mpirank] + quad_idx;
      mesh->num_nodes_of_elem.val[quad_idx] = P4EST_CHILDREN;

      for (int lnode = 0; lnode < P4EST_CHILDREN; ++lnode) {
         /* Cache some information on corner nodes. */
         p4est_locidx_t node_idx = lnodes->element_nodes[P4EST_CHILDREN * quad_idx + lnode];
         //      isboundary[i] = (bc == NULL ? 0 : bc[lni]);
         //       inloc[i] = !isboundary[i] ? in[lni] : 0.;
         mesh->elem_node_indices.val[P4EST_CHILDREN * quad_idx + lnode] = node_idx;
      }

      mesh->element_lengths.val[quad_idx] = pow(0.5, quad->level);

      /* Figure out the hanging corners on this element, if any. */
      int hanging_corner[P4EST_CHILDREN];
      int anyhang = lnodes_decode2 (lnodes->face_code[quad_idx], hanging_corner);

      p4est_quadrant_t* parent;
      if (!anyhang) {
         parent = NULL;          /* Defensive programming. */
      }
      else {
         /* At least one node is hanging.  We need the parent quadrant to
           * find the location of the corresponding non-hanging node. */
         parent = &sp;
         p4est_quadrant_parent (quad, parent);
      }


      for (int lnode = 0; lnode < P4EST_CHILDREN; ++lnode) {
         p4est_locidx_t node_idx = lnodes->element_nodes[P4EST_CHILDREN * quad_idx + lnode];
         P4EST_ASSERT (node_idx >= 0 && node_idx < subdomain_dims->n_nodes);
         if (anyhang && hanging_corner[lnode] >= 0) {
            /* This node is hanging; access the referenced node instead. */
            p4est_quadrant_corner_node (parent, lnode, &node);
         }
         else {
            p4est_quadrant_corner_node (quad, lnode, &node);
         }

         /* Transform per-tree reference coordinates into physical space. */
         p4est_qcoord_to_vertex (p4est->connectivity, tt, node.x, node.y,
#ifdef P4_TO_P8
                                 node.z,
#endif
                                 vxyz);

         mesh->coords.val[0][node_idx] = vxyz[0];
         mesh->coords.val[1][node_idx] = vxyz[1];
#ifdef P4_TO_P8
         mesh->coords.val[2][node_idx] = vxyz[2];
#endif



//         PPP printf("(%3.2lf, %3.2lf), ", vxyz[0], vxyz[1]);

      }
//      PPP printf("\n");

   }}}//for all quads

}


const double EPS = 1e-10;
bool real_equal(real a, real b)
{
   return fabs(a - b) < EPS;
}

void prepare_subdomain_fem_space(BddcmlMesh *mesh, BddcmlFemSpace *femsp)
{
   init_fem_space(mesh->subdomain_dims, femsp);
   int num_dofs_per_node = mesh->subdomain_dims->n_node_dofs;
   for(int node = 0; node < mesh->subdomain_dims->n_nodes; node++)
   {      
      bool is_on_boundary = ((real_equal(mesh->coords.val[0][node], 0.0)) || (real_equal(mesh->coords.val[0][node], 1.0))
            || (real_equal(mesh->coords.val[1][node], 0.0)) || (real_equal(mesh->coords.val[1][node], 1.0))
#ifdef P4_TO_P8
            || (real_equal(mesh->coords.val[2][node], 0.0)) || (real_equal(mesh->coords.val[2][node], 1.0))
#endif
            );

      femsp->node_num_dofs.val[node] = num_dofs_per_node;

      for(int local_dof_idx = 0; local_dof_idx < num_dofs_per_node; local_dof_idx++)
      {
         int subdomain_dof = num_dofs_per_node * node + local_dof_idx;
         int global_dof = num_dofs_per_node * mesh->node_global_map.val[node] + local_dof_idx;

         femsp->dofs_global_map.val[subdomain_dof] = global_dof;

         if(is_on_boundary)
         {
            femsp->fixs_code.val[subdomain_dof] = 1;
            femsp->fixs_values.val[subdomain_dof] = 0.0;
         }
         else
         {
            femsp->fixs_code.val[subdomain_dof] = 0;
            femsp->fixs_values.val[subdomain_dof] = 0.0;
         }

      }

   }
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

void assemble_matrix_rhs(p4est_lnodes_t *lnodes, BddcmlMesh *mesh, BddcmlFemSpace *femsp,
                         SparseMatrix *matrix, RealArray *rhss)
{
   real i_coeffs, j_coeffs;
//   real mass_ref[P4EST_CHILDREN][P4EST_CHILDREN];
//   real dudv_ref[P4EST_CHILDREN][P4EST_CHILDREN];
   real dudv_phys_elem[P4EST_CHILDREN][P4EST_CHILDREN];
   //generate_reference_matrices(dudv_ref, mass_ref);

   int element_offset = 0;
   for(int elem_idx = 0; elem_idx < mesh->subdomain_dims->n_elems; elem_idx++) {
      int num_nodes_of_elem = mesh->num_nodes_of_elem.val[elem_idx];
      int ndof_per_element = num_nodes_of_elem;
      assert(ndof_per_element == P4EST_CHILDREN);

      // TODO: elem_size and elem_volume is correct only when the mesh is obtained by refinements
      // from a UNIT SQUARE/CUBE
      // TODO: ONLY FROM ****UNIT****
      real elem_length = mesh->element_lengths.val[elem_idx];
      real elem_volume = pow(elem_length, mesh->subdomain_dims->n_problem_dims);

//      double reference_scaled =
//#ifndef P4_TO_P8
//            1.;
//#else
//            elem_length;
//#endif

      //scale_reference_matrix(dudv_ref, reference_scaled, dudv_phys_elem);
      generate_scaled_matrix_new(elem_length, dudv_phys_elem);


      for(int j = 0; j < ndof_per_element; j++) {

         //todo: dofs should be taken from femsp!
         int jdof = mesh->elem_node_indices.val[element_offset + j];
         assert(femsp->node_num_dofs.val[jdof] == 1);

         p4est_locidx_t j_nodes[4];
         int j_nindep = independent_nodes(lnodes, elem_idx, j, j_nodes, &j_coeffs);
         if(j_nindep == 1)
         {
            assert(jdof == j_nodes[0]);
         }
         for(int j_indep_nodes_idx = 0; j_indep_nodes_idx < j_nindep; j_indep_nodes_idx++)
         {
            int j_indep_node = j_nodes[j_indep_nodes_idx];

            for(int i = 0; i < ndof_per_element /*<= j*/; i++) {
               int idof = mesh->elem_node_indices.val[element_offset + i];

               p4est_locidx_t i_nodes[4];
               int i_nindep = independent_nodes(lnodes, elem_idx, i, i_nodes, &i_coeffs);

               if(i_nindep == 1)
               {
                  assert(idof == i_nodes[0]);
               }
               for(int i_indep_nodes_idx = 0; i_indep_nodes_idx < i_nindep; i_indep_nodes_idx++)
               {
                  int i_indep_node = i_nodes[i_indep_nodes_idx];

                  double matrix_value = i_coeffs * j_coeffs * dudv_phys_elem[j][i];
                  add_matrix_entry(matrix, i_indep_node, j_indep_node, matrix_value);
//                  printf("adding entry loc (%d, %d), nodes, (%d, %d), coefs (%3.2lf, %3.2lf), number indep (%d, %d), locstiff %lf, value %lf\n",
//                         j, i, j_indep_node, i_indep_node, j_coeffs, i_coeffs, j_nindep, i_nindep, stiffness_dd[j][i], matrix_value);
               }
            }

            // TODO: integrate properly
            // TODO: elem_volume is correct only when the mesh is obtained by refinements
            // from a UNIT SQUARE/CUBE
            // TODO: ONLY FROM ****UNIT****


            double rhs_value = j_coeffs * 1./(real)P4EST_CHILDREN * elem_volume * 1;
            rhss->val[j_indep_node] += rhs_value;

         }
      }
      element_offset += num_nodes_of_elem;
   }

}

