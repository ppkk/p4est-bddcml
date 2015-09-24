#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <mpi.h>

extern "C"{
#include "bddcml_interface_c.h"
}

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

//*******************************************************************************************

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

}

void init_dimmensions(BddcmlDimensions* dimmensions, int mesh_dim, PhysicsType physicsType)
{
   dimmensions->n_problem_dims = mesh_dim;
   dimmensions->n_mesh_dims = mesh_dim;
   dimmensions->n_dofs = 0;
   dimmensions->n_elems = 0;
   dimmensions->n_nodes = 0;

   if(mesh_dim == 2)
      dimmensions->n_elem_nodes = 4;
   else if(mesh_dim == 3)
      dimmensions->n_elem_nodes = 8;
   else
      assert(0);

   if(physicsType == PhysicsType::LAPLACE)
      dimmensions->n_node_dofs = 1;
   else if(physicsType == PhysicsType::LINEAR_ELASTICITY)
      dimmensions->n_node_dofs = 3;
   else
      assert(0);
}

//*******************************************************************************************

void init_mesh(BddcmlDimensions* subdomain_dims, BddcmlMesh* mesh)
{
   int nnodes_per_elem;
   if(subdomain_dims->n_mesh_dims == 2)
      nnodes_per_elem = 4;
   else if(subdomain_dims->n_mesh_dims == 3)
      nnodes_per_elem = 8;
   else
      assert(0);

   mesh->subdomain_dims = subdomain_dims;
   allocate_idx_array(subdomain_dims->n_elems * nnodes_per_elem, &mesh->elem_node_indices);
   allocate_idx_array(subdomain_dims->n_elems, &mesh->num_nodes_of_elem);
   allocate_idx_array(subdomain_dims->n_elems, &mesh->elem_global_map);
   allocate_idx_array(subdomain_dims->n_nodes, &mesh->node_global_map);
   allocate_real_array(subdomain_dims->n_elems, &mesh->element_lengths);
   allocate_real_2D_array(subdomain_dims->n_nodes, subdomain_dims->n_problem_dims, &mesh->coords);
}

void free_mesh(BddcmlMesh* mesh)
{
   free_idx_array(&mesh->elem_node_indices);
   free_idx_array(&mesh->num_nodes_of_elem);
   free_idx_array(&mesh->elem_global_map);
   free_idx_array(&mesh->node_global_map);
   free_real_array(&mesh->element_lengths);
   free_real_2D_array(&mesh->coords);
}

void print_bddcml_mesh(BddcmlMesh *mesh, int which_rank)
{
   print_rank = which_rank;
   PPP printf("\n*************** BEGIN BDDCML MESH ************************\n");
   PPP printf("elems: %d, nodes: %d\n", mesh->subdomain_dims->n_elems, mesh->subdomain_dims->n_nodes);
   PPP printf("linet: %d, lnnet: %d\n", mesh->elem_node_indices.len, mesh->num_nodes_of_elem.len);
   for(int elem = 0; elem < mesh->elem_global_map.len; elem++)
   {
      PPP printf("elem %d -> ", mesh->elem_global_map.val[elem]);
      for (int lnode = 0; lnode < mesh->num_nodes_of_elem.val[elem]; lnode++)
      {
         int node_local_idx = mesh->elem_node_indices.val[elem*mesh->subdomain_dims->n_elem_nodes + lnode];
         int node_global_idx = mesh->node_global_map.val[node_local_idx];
         PPP printf("%d[%d], ", node_global_idx, node_local_idx);
      }
      for (int lnode = 0; lnode < mesh->num_nodes_of_elem.val[elem]; lnode++)
      {
         int node_local_idx = mesh->elem_node_indices.val[elem*mesh->subdomain_dims->n_elem_nodes + lnode];
         PPP printf("(%3.2lf, %3.2lf), ", mesh->coords.val[0][node_local_idx], mesh->coords.val[1][node_local_idx]);
      }
      PPP printf("\n");
   }
   PPP printf("*************** END BDDCML MESH ************************\n\n");
}

//*******************************************************************************************


void init_fem_space(BddcmlDimensions* dims, BddcmlFemSpace* femsp)
{
   femsp->subdomain_dims = dims;
   allocate_idx_array(dims->n_nodes, &femsp->node_num_dofs);
   allocate_idx_array(dims->n_dofs, &femsp->dofs_global_map);
   allocate_idx_array(dims->n_dofs, &femsp->fixs_code);
   allocate_real_array(dims->n_dofs, &femsp->fixs_values);
}

void free_fem_space(BddcmlFemSpace* femsp)
{
   free_idx_array(&femsp->node_num_dofs);
   free_idx_array(&femsp->dofs_global_map);
   free_idx_array(&femsp->fixs_code);
   free_real_array(&femsp->fixs_values);
}

void print_bddcml_fem_space(BddcmlFemSpace *femsp, BddcmlMesh *mesh, int which_rank)
{
   print_rank = which_rank;
   PPP printf("\n*************** BEGIN BDDCML FEM SPACE ************************\n");
   PPP printf("dofs: %d\n", femsp->subdomain_dims->n_dofs);
   assert(femsp->subdomain_dims->n_dofs == femsp->subdomain_dims->n_nodes);
   for(int node = 0; node < femsp->subdomain_dims->n_dofs; node++)
   {
      assert(femsp->node_num_dofs.val[node] == 1);
      PPP printf("node %d -> num dofs %d, gdofs: (%d ), bc code %d, bc val %6.4lf, coords (%4.3lf, %4.3lf)\n", node, femsp->node_num_dofs.val[node],
                 femsp->dofs_global_map.val[node], femsp->fixs_code.val[node], femsp->fixs_values.val[node],
                 mesh->coords.val[0][node], mesh->coords.val[1][node]);
   }
   PPP printf("*************** END BDDCML FEM SPACE ************************\n\n");
}

//*******************************************************************************************


// Basic properties
void print_basic_properties(BddcmlDimensions *global_dims, int num_subdomains, BddcmlLevelInfo *level_info, BddcmlKrylovParams *krylov_params)
{
   if (mpi_rank == print_rank) {
      printf("Characteristics of the problem :\n");
      printf("  number of processors            nproc = %d\n" ,mpi_size);
      printf("  number of dimensions             ndim = %d\n", global_dims->n_problem_dims);
      printf("  mesh dimension                meshdim = %d\n", global_dims->n_mesh_dims);
      printf("  number of elements global       nelem = %d\n", global_dims->n_elems);
      printf("  number of subdomains             nsub = %d\n", num_subdomains);
      printf("  number of nodes global           nnod = %d\n", global_dims->n_nodes);
      printf("  number of DOF                    ndof = %d\n", global_dims->n_dofs);
      printf("  number of levels              nlevels = %d\n", level_info->nlevels);
      printf("  number of subdomains in levels        = ");
      for(int idx = 0; idx < level_info->nlevels; idx++) {
         printf("%d, ", level_info->nsublev[idx]);
      }
      printf("\n");
      printf("Characteristics of iterational process:\n");
      printf("  tolerance of error                tol = %lf\n", krylov_params->tol);
      printf("  maximum number of iterations    maxit = %d\n", krylov_params->maxit);
      printf("  number of incresing residual ndecrmax = %d\n", krylov_params->ndecrmax);
      printf("  using recycling of Krylov method ?      %d\n", krylov_params->recycling_int);
   }
}
//*******************************************************************************************


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
                                  mesh->coords.val_serialized,
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
                                  (int*) &matrix->type,
                                  matrix->i,
                                  matrix->j,
                                  matrix->val,
                                  &matrix->nnz,
                                  &matrix->is_assembled,
                                  user_constraints->val_serialized,
                                  &user_constraints->len1,
                                  &user_constraints->len2,
                                  element_data->val_serialized,
                                  &element_data->len1,
                                  &element_data->len2,
                                  dof_data->val,
                                  &dof_data->len,
                                  &preconditioner_params->find_components_int);

}

void bddcml_setup_preconditioner(MatrixType matrixtype, BddcmlPreconditionerParams *params)
{
   bddcml_setup_preconditioner_c((int*) &matrixtype,
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

void bddcml_download_local_solution(int isub, RealArray *sols)
{
   bddcml_download_local_solution_c(&isub,
                                    sols->val,
                                    &sols->len);
}

void bddcml_dotprod_subdomain(int isub, RealArray *sols1, RealArray *sols2, real *normRn2_sub)
{
   bddcml_dotprod_subdomain_c( &isub, sols1->val, &sols1->len, sols2->val, &sols2->len, normRn2_sub );
}

void bddcml_finalize()
{
   bddcml_finalize_c();
}
