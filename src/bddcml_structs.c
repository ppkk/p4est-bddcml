#include <stdio.h>
#include <mpi.h>
#include "bddcml_structs.h"

void set_implicit(BddcmlGeneralParams *params)
{
   params->numbase = 0;
   params->just_direct_solve_int = 0;
   params->verbose_level = 1;
   params->export_solution = 1;
   strcpy(params->output_file_prefix, "poisson_solution");
}

void initialize_levels(int n_subdomains_first_level, BddcmlLevelInfo *level_info)
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
