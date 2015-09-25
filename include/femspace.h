#ifndef FEMSPACE_H
#define FEMSPACE_H

#include "arrays.h"
#include "bddcml_structs.h"

typedef struct BddcmlFemSpace
{
   BddcmlDimensions* subdomain_dims;

   IdxArray node_num_dofs;
   IdxArray fixs_code;
   RealArray fixs_values;
   IdxArray dofs_global_map;
}
BddcmlFemSpace;

void init_fem_space(BddcmlDimensions* dims, BddcmlFemSpace* femsp);
void free_fem_space(BddcmlFemSpace* femsp);
void print_bddcml_fem_space(BddcmlFemSpace *femsp, BddcmlMesh *mesh, int which_rank);

#endif // FEMSPACE_H
