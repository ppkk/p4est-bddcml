#ifndef FEMSPACE_H
#define FEMSPACE_H

#include "arrays.h"
#include "bddcml_structs.h"

// todo this is here because of exact solution, remove it!
#include "integral.h"

class BddcmlFemSpace
{
public:
   BddcmlFemSpace(const BddcmlMesh* mesh);
   ~BddcmlFemSpace();
   void prepare_subdomain_fem_space(PhysicsType physicsType, exact_fn dirichlet_bc_exact);
   void print(int which_rank) const;

public:
   const BddcmlDimensions* subdomain_dims;

   IdxArray node_num_dofs;
   IdxArray fixs_code;
   RealArray fixs_values;
   IdxArray dofs_global_map;
   PhysicsType physicsType;

private:
   const BddcmlMesh* mesh;
};



#endif // FEMSPACE_H
