#ifndef BDDCML_CUBE_EXAMPLE_H
#define BDDCML_CUBE_EXAMPLE_H

#include "bddcml_structs.h"

void prepare_subdomain_data(int isub, // global subdomain index
                            int num_sub_per_cube_edge, // number of subdomains in one edge of a cube
                            int num_el_per_sub_edge,  // number of elements in one edge of a subdomain
                            real hsize, // element size
                            BddcmlMesh *mesh, BddcmlFemSpace *femsp);

#endif // BDDCML_CUBE_EXAMPLE_H
