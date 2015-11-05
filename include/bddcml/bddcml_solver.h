#ifndef BDDCML_SOLVER_H
#define BDDCML_SOLVER_H

#include <memory>

#include "bddcml_structs.h"

class NodalElementMesh;
class DiscreteSystem;

class BddcmlSolver
{
public:
   BddcmlSolver(ProblemDimensions &subdomain_dims, ProblemDimensions &global_dims,
                BddcmlGeneralParams &general_params, BddcmlKrylovParams &krylov_params,
                BddcmlPreconditionerParams &preconditioner_params, const P4estClass &p4est_class, int num_levels);

   void solve(const NodalElementMesh &nodal_mesh, DiscreteSystem &system, std::vector<double> *sols);
   void clear();
private:
   ProblemDimensions &subdomain_dims;
   ProblemDimensions &global_dims;

   // todo make it const
   BddcmlGeneralParams &general_params;
   BddcmlKrylovParams &krylov_params;
   BddcmlPreconditionerParams &preconditioner_params;
   const P4estClass &p4est_class;
   int num_levels;

//   std::shared_ptr<BddcmlMesh> mesh;
//   std::shared_ptr<BddcmlFemSpace> femspace;

};

#endif // BDDCML_SOLVER_H
