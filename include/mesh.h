#ifndef MESH_H
#define MESH_H

#include "arrays.h"
#include "bddcml_structs.h"
#include "my_p4est_interface.h"

class Element;

class BddcmlMesh
{
public:
   BddcmlMesh(BddcmlDimensions* subdomain_dims) {init(subdomain_dims);}
   ~BddcmlMesh() {free(); }

   void print(int which_rank) const;

   // assuming that elements are squares/cubes aligned with cartesian grid in natural order of axes (x, y, z)
   // i.e. first node of the element has lowest x, y and z
   void get_element(int elem_idx, Element* element) const;

private:
   void init(BddcmlDimensions* subdomain_dims);
   void free();

public:
   BddcmlDimensions* subdomain_dims;

   IdxArray elem_node_indices;
   IdxArray num_nodes_of_elem;

   // remember, that in the case of hanging nodes, we store the idx (and thus coordinates) of the correspondent regular node
   Real2DArray coords;

   IdxArray elem_global_map;
   IdxArray node_global_map;

   RealArray element_lengths;
};



// **************************
// BDDCML MESH DIMMENSIONS
// **************************
class BddcmlDimensions
{   
public:
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

};



void init_dimmensions(BddcmlDimensions* dimmensions, int mesh_dim, PhysicsType physicsType);



#endif // MESH_H
