#ifndef FEMSPACE_H
#define FEMSPACE_H

#include "arrays.h"
#include "bddcml_structs.h"

class BddcmlFemSpace
{
public:
   BddcmlFemSpace(const BddcmlMesh* mesh);
   ~BddcmlFemSpace();
   void prepare_subdomain_fem_space(PhysicsType physicsType);
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
