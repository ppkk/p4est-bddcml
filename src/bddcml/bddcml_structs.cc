#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <mpi.h>

extern "C"{
#include "bddcml_interface_c.h"
}

#include "bddcml/bddcml_structs.h"
#include "bddcml/bddcml_mesh.h"
#include "bddcml/bddcml_femspace.h"

BddcmlGeneralParams::BddcmlGeneralParams()
{
   numbase = 0;
   just_direct_solve_int = 0;
   verbose_level = 1;
   export_solution = 1;
   strcpy(output_file_prefix, "poisson_solution");
}

BddcmlKrylovParams::BddcmlKrylovParams()
{
   krylov_method = 0;
   recycling_int = 1;
   max_number_of_stored_vectors = 50;
   maxit = 500;
   ndecrmax = 50;
   tol = 1.e-6;
}

BddcmlPreconditionerParams::BddcmlPreconditionerParams()
{
   use_preconditioner_defaults = 0;
   use_corner_constraints = 1;
   use_arithmetic_constraints = 1;
   use_adaptive_constraints = 0;
   use_user_constraints = 0;
   weights_type = 0;
   parallel_division = 1;
   find_components_int = 1;
}

//*******************************************************************************************

BddcmlLevelInfo::BddcmlLevelInfo(int n_levels, int n_subdomains_first_level)
{
   real coarsening;
   int i, ir;

   // initialize levels
   nlevels = n_levels;
   lnsublev = nlevels;
   nsublev = (int*) malloc(lnsublev * sizeof(int));
   if (nlevels == 2) {
      nsublev[0] = n_subdomains_first_level;
      nsublev[1] = 1;
   }
   else if (nlevels > 2) {
      // determine coarsening factor
      coarsening = pow(n_subdomains_first_level, 1./(nlevels-1));
      // prescribe number of subdomains on levels so that coarsening is fixed between levels
      nsublev[0] = n_subdomains_first_level;
      for( i = 1; i < nlevels - 1; i++) {
         ir = nlevels - i;
         nsublev[i] = (int)pow(coarsening, ir-1);
         if (nsublev[i] % 2 != 0) {
            nsublev[i] = nsublev[i] + 1;
         }
      }
      nsublev[nlevels-1] = 1;
   }
   else {
      printf("Unsupported number of levels: %d\n", nlevels);
      exit(0);
   }

}



//*******************************************************************************************


// Basic properties
void print_basic_properties(const BddcmlDimensions &global_dims, int num_subdomains,
                            const BddcmlLevelInfo &level_info, const BddcmlKrylovParams &krylov_params)
{
   if (mpi_rank == print_rank) {
      printf("Characteristics of the problem :\n");
      printf("  number of processors            nproc = %d\n" ,mpi_size);
      printf("  number of dimensions             ndim = %d\n", global_dims.n_problem_dims);
      printf("  mesh dimension                meshdim = %d\n", global_dims.n_mesh_dims);
      printf("  number of elements global       nelem = %d\n", global_dims.n_elems);
      printf("  number of subdomains             nsub = %d\n", num_subdomains);
      printf("  number of nodes global           nnod = %d\n", global_dims.n_nodes);
      printf("  number of DOF                    ndof = %d\n", global_dims.n_dofs);
      printf("  number of levels              nlevels = %d\n", level_info.nlevels);
      printf("  number of subdomains in levels        = ");
      for(int idx = 0; idx < level_info.nlevels; idx++) {
         printf("%d, ", level_info.nsublev[idx]);
      }
      printf("\n");
      printf("Characteristics of iterational process:\n");
      printf("  tolerance of error                tol = %lf\n", krylov_params.tol);
      printf("  maximum number of iterations    maxit = %d\n", krylov_params.maxit);
      printf("  number of incresing residual ndecrmax = %d\n", krylov_params.ndecrmax);
      printf("  using recycling of Krylov method ?      %d\n", krylov_params.recycling_int);
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
