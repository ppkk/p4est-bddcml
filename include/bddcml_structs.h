#ifndef BDDCML_INTERFACE_H
#define BDDCML_INTERFACE_H

#include <mpi.h>

// type used for floating point
typedef double real;

// **************************
// GENERAL BDDCML PARAMETERS:
// **************************
typedef struct BddcmlGeneralParams
{
   // beginning index of arrays ( 0 for C, 1 for Fortran )
   int numbase;

   // Just a direct solve by MUMPS?
   int just_direct_solve_int;

   // verbosity of BDDCML ( 0 - only fatal errors, 1 - mild output, 2 - detailed output )
   int verbose_level;

   // export solution to VTU files?
   int export_solution;

   // what is the name of that file (resp. collection of files)
   char output_file_prefix[255];

} BddcmlGeneralParams;

void set_implicit(BddcmlGeneralParams *params);

// **************************
// BDDCML LEVELS INFORMATION
// **************************
typedef struct BddcmlLevelInfo
{
   int nlevels; // number of levels

   // subdomains in levels
   int lnsublev;
   int *nsublev;
   int nsub_loc_1;
}BddcmlLevelInfo;

void initialize_levels(int n_subdomains_first_level, BddcmlLevelInfo *level_info);


void bddcml_init(BddcmlGeneralParams *general_params, BddcmlLevelInfo *level_info, MPI_Comm communicator);


#endif // BDDCML_INTERFACE_H
