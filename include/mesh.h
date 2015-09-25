#ifndef MESH_H
#define MESH_H

#include "arrays.h"
#include "bddcml_structs.h"
#include "p4est_common.h"

typedef struct BddcmlMesh
{
   BddcmlDimensions* subdomain_dims;

   IdxArray elem_node_indices;
   IdxArray num_nodes_of_elem;

   Real2DArray coords;

   IdxArray elem_global_map;
   IdxArray node_global_map;

   RealArray element_lengths;
}
BddcmlMesh;

void init_mesh(BddcmlDimensions* subdomain_dims, BddcmlMesh* mesh);
void free_mesh(BddcmlMesh* mesh);
void print_bddcml_mesh(BddcmlMesh* mesh, int which_rank);


// **************************
// BDDCML MESH DIMMENSIONS
// **************************
typedef struct BddcmlDimensions
{
   int n_elems;  // number of elements
   int n_nodes;  // number of nodes
   int n_dofs;   // number on degrees of freedom

   // spacial dimension
   int n_problem_dims;

   // topological dimension of elements elements, would be lower for shells or beams
   int n_mesh_dims;

   // number of nodes per element (4 or 8)
   int n_elem_nodes;

   // number of dofs per node (1 for Laplace, 3 for elasticity)
   int n_node_dofs;

} BddcmlDimensions;



void init_dimmensions(BddcmlDimensions* dimmensions, int mesh_dim, PhysicsType physicsType);
void prepare_dimmensions(p4est_t *p4est, p4est_lnodes_t *lnodes, PhysicsType physicsType,
                         BddcmlDimensions *subdomain_dims, BddcmlDimensions *global_dims,
                         sc_MPI_Comm mpicomm);


#endif // MESH_H
