#ifndef BDDCML_CUBE_EXAMPLE_H
#define BDDCML_CUBE_EXAMPLE_H

#include "bddcml_structs.h"

void prepare_subdomain_data(int isub, // global subdomain index
                            int num_sub_per_cube_edge, // number of subdomains in one edge of a cube
                            int num_el_per_sub_edge,  // number of elements in one edge of a subdomain
                            real hsize, // element size
                            BddcmlConnectivity *connectivity, IdxArray *node_num_dofs,
                            IdxArray *elem_global_map, IdxArray *node_global_map, IdxArray *dofs_global_map,
                            real* xyzs, int lxyzs1, int lxyzs2,
                            int* ifixs, int lifixs, real* fixvs, int lfixvs);

#endif // BDDCML_CUBE_EXAMPLE_H
