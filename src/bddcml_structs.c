#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include "bddcml_interface_c.h"
#include "bddcml_structs.h"

void set_implicit_general_params(BddcmlGeneralParams *params)
{
   params->numbase = 0;
   params->just_direct_solve_int = 0;
   params->verbose_level = 1;
   params->export_solution = 1;
   strcpy(params->output_file_prefix, "poisson_solution");
}

void set_implicit_krylov_params(BddcmlKrylovParams *params)
{
   params->krylov_method = 0;
   params->recycling_int = 1;
   params->max_number_of_stored_vectors = 50;
   params->maxit = 500;
   params->ndecrmax = 50;
   params->tol = 1.e-6;
}

void set_implicit_preconditioner_params(BddcmlPreconditionerParams *params)
{
   params->use_preconditioner_defaults = 0;
   params->use_corner_constraints = 1;
   params->use_arithmetic_constraints = 1;
   params->use_adaptive_constraints = 0;
   params->use_user_constraints = 0;
   params->weights_type = 0;
   params->parallel_division = 1;
   params->find_components_int = 1;
}

void init_levels(int n_subdomains_first_level, BddcmlLevelInfo *level_info)
{
   real coarsening;
   int i, ir;

   // initialize levels
   level_info->lnsublev = level_info->nlevels;
   level_info->nsublev = (int*) malloc(level_info->lnsublev * sizeof(int));
   if (level_info->nlevels == 2) {
      level_info->nsublev[0] = n_subdomains_first_level;
      level_info->nsublev[1] = 1;
   }
   else if (level_info->nlevels > 2) {
      // determine coarsening factor
      coarsening = pow(n_subdomains_first_level, 1./(level_info->nlevels-1));
      // prescribe number of subdomains on levels so that coarsening is fixed between levels
      level_info->nsublev[0] = n_subdomains_first_level;
      for( i = 1; i < level_info->nlevels - 1; i++) {
         ir = level_info->nlevels - i;
         level_info->nsublev[i] = (int)pow(coarsening, ir-1);
         if (level_info->nsublev[i] % 2 != 0) {
            level_info->nsublev[i] = level_info->nsublev[i] + 1;
         }
      }
      level_info->nsublev[level_info->nlevels-1] = 1;
   }
   else {
      printf("Unsupported number of levels: %d\n", level_info->nlevels);
      exit(0);
   }

   printf("XXX\nXXX\n nu sublev, %d, %d, %d\n", level_info->nsublev[0], level_info->nsublev[1], level_info->nsublev[2]);

}

void init_dimmensions(BddcmlDimensions* dimmensions, int mesh_dim)
{
   dimmensions->n_problem_dims = mesh_dim;
   dimmensions->n_mesh_dims = mesh_dim;
   dimmensions->n_dofs = 0;
   dimmensions->n_elems = 0;
   dimmensions->n_nodes = 0;
}

void init_mesh(BddcmlDimensions* dims, BddcmlMesh* mesh)
{
   mesh->subdomain_dims = dims;
   allocate_idx_array(dims->n_elems * 8, &mesh->elem_node_indices);
   allocate_idx_array(dims->n_elems, &mesh->num_nodes_of_elem);
   allocate_idx_array(dims->n_elems, &mesh->elem_global_map);
   allocate_idx_array(dims->n_nodes, &mesh->node_global_map);
   allocate_real_2D_array(dims->n_nodes, dims->n_problem_dims, &mesh->coords);
}

void free_mesh(BddcmlMesh* mesh)
{
   free_idx_array(&mesh->elem_node_indices);
   free_idx_array(&mesh->num_nodes_of_elem);
   free_idx_array(&mesh->elem_global_map);
   free_idx_array(&mesh->node_global_map);
   free_real_2D_array(&mesh->coords);
}

void init_fem_space(BddcmlDimensions* dims, BddcmlFemSpace* femsp)
{
   femsp->subdomain_dims = dims;
   allocate_idx_array(dims->n_nodes, &femsp->node_num_dofs);
   allocate_idx_array(dims->n_dofs, &femsp->dofs_global_map);
   allocate_idx_array(dims->n_nodes, &femsp->fixs_code);
   allocate_real_array(dims->n_nodes, &femsp->fixs_values);
}

void free_fem_space(BddcmlFemSpace* femsp)
{
   free_idx_array(&femsp->node_num_dofs);
   free_idx_array(&femsp->dofs_global_map);
   free_idx_array(&femsp->fixs_code);
   free_real_array(&femsp->fixs_values);
}

void bddcml_init(BddcmlGeneralParams *general_params, BddcmlLevelInfo *level_info, MPI_Comm communicator)
{
   int fortran_comm =  MPI_Comm_c2f(communicator);

   bddcml_init_c(&level_info->nlevels,
                 level_info->nsublev,
                 &level_info->lnsublev,
                 &level_info->nsub_loc_1,
                 &fortran_comm,
                 &general_params->verbose_level,
                 &general_params->numbase,
                 &general_params->just_direct_solve_int);
}

void bddcml_upload_subdomain_data(BddcmlDimensions *global_dims, BddcmlDimensions *subdomain_dims,
                                  int isub, BddcmlMesh *mesh, BddcmlFemSpace *femsp,
                                  RealArray *rhss, int is_rhs_complete, RealArray *sols, SparseMatrix *matrix,
                                  Real2DArray *user_constraints, Real2DArray *element_data,
                                  RealArray *dof_data, BddcmlPreconditionerParams* preconditioner_params)
{
   bddcml_upload_subdomain_data_c(&global_dims->n_elems,
                                  &global_dims->n_nodes,
                                  &global_dims->n_dofs,
                                  &global_dims->n_problem_dims,
                                  &global_dims->n_mesh_dims,
                                  &isub,
                                  &subdomain_dims->n_elems,
                                  &subdomain_dims->n_nodes,
                                  &subdomain_dims->n_dofs,
                                  mesh->elem_node_indices.val,
                                  &mesh->elem_node_indices.len,
                                  mesh->num_nodes_of_elem.val,
                                  &mesh->num_nodes_of_elem.len,
                                  femsp->node_num_dofs.val,
                                  &femsp->node_num_dofs.len,
                                  mesh->node_global_map.val,
                                  &mesh->node_global_map.len,
                                  femsp->dofs_global_map.val,
                                  &femsp->dofs_global_map.len,
                                  mesh->elem_global_map.val,
                                  &mesh->elem_global_map.len,
                                  mesh->coords.val,
                                  &mesh->coords.len1,
                                  &mesh->coords.len2,
                                  femsp->fixs_code.val,
                                  &femsp->fixs_code.len,
                                  femsp->fixs_values.val,
                                  &femsp->fixs_values.len,
                                  rhss->val,
                                  &rhss->len,
                                  &is_rhs_complete,
                                  sols->val,
                                  &sols->len,
                                  &matrix->type,
                                  matrix->i,
                                  matrix->j,
                                  matrix->val,
                                  &matrix->len,
                                  &matrix->is_assembled,
                                  user_constraints->val,
                                  &user_constraints->len1,
                                  &user_constraints->len2,
                                  element_data->val,
                                  &element_data->len1,
                                  &element_data->len2,
                                  dof_data->val,
                                  &dof_data->len,
                                  &preconditioner_params->find_components_int);

}

void bddcml_setup_preconditioner(int matrixtype, BddcmlPreconditionerParams *params)
{
   bddcml_setup_preconditioner_c(&matrixtype,
                                 &params->use_preconditioner_defaults,
                                 &params->parallel_division,
                                 &params->use_corner_constraints,
                                 &params->use_arithmetic_constraints,
                                 &params->use_adaptive_constraints,
                                 &params->use_user_constraints,
                                 &params->weights_type);
}

void bddcml_solve(BddcmlKrylovParams *krylov_params, BddcmlConvergenceInfo *convergence_info, MPI_Comm communicator)
{
   int fortran_comm =  MPI_Comm_c2f(communicator);

   bddcml_solve_c(&fortran_comm,
                  &krylov_params->krylov_method,
                  &krylov_params->tol,
                  &krylov_params->maxit,
                  &krylov_params->ndecrmax,
                  &krylov_params->recycling_int,
                  &krylov_params->max_number_of_stored_vectors,
                  &convergence_info->num_iter,
                  &convergence_info->converged_reason,
                  &convergence_info->condition_number);
}

