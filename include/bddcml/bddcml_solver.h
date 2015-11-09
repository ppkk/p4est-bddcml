#ifndef BDDCML_SOLVER_H
#define BDDCML_SOLVER_H

#include <memory>

#include "bddcml_structs.h"
//#include "integral.h"

class NodalElementMesh;
class DiscreteSystem;

class BddcmlSolver
{
public:
   BddcmlSolver(ProblemDimensions &problem_dims, BddcmlGeneralParams &general_params, BddcmlKrylovParams &krylov_params,
                BddcmlPreconditionerParams &preconditioner_params, const P4estClass &p4est_class, int num_levels);

   void solve(const NodalElementMesh &nodal_mesh, DiscreteSystem &system, exact_fn dirichlet_bc_exact, std::vector<double> *sols);
   void clear();
private:
   ProblemDimensions &problem_dims;

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
