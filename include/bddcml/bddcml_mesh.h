#ifndef MESH_H
#define MESH_H

#include <stdio.h>
#include "arrays.h"

class ProblemDimensions;
class BddcmlMesh;
class BddcmlFemSpace;
class NodalElementMesh;
class P4estClass;

class BddcmlMesh
{
public:
   BddcmlMesh(const ProblemDimensions *subdomain_dims) {init(subdomain_dims);}
   ~BddcmlMesh() {free(); }
   void fill_nodes_info(const P4estClass &p4est, const NodalElementMesh &nodal_mesh);

   void print(int which_rank) const;

private:
   void init(const ProblemDimensions *subdomain_dims);
   void free();

public:
   const ProblemDimensions* subdomain_dims;

   IdxArray elem_node_indices;
   IdxArray num_nodes_of_elem;

   // remember, that in the case of hanging nodes, we store the idx (and thus coordinates) of the correspondent regular node
   Real2DArray coords;

   IdxArray elem_global_map;
   IdxArray node_global_map;
};

#endif // MESH_H
