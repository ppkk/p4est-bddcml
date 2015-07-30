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

