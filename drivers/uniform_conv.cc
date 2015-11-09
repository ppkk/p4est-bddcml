#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <assert.h>
#include <math.h>
#include <fstream>

#include "bddcml/bddcml_solver.h"
#include "p4est/my_p4est_interface.h"
#include "assemble.h"
#include "integration_cell.h"
#include "element.h"
#include "shapefun.h"
#include "vtk_output.h"
#include "integral.h"
#include "analytical_solution.h"

using namespace std;

const int num_dim = 3;
const int order = 1;
const int norm_order = 2 * order;
const PhysicsType physicsType = PhysicsType::LAPLACE;
const ExactSolID exact_sol_id = ExactSolID::InternalLayer;

vector<double> rhs_fn(vector<double> coords) {
   ExactSolution sol;
   sol.calculate(exact_sol_id, coords);
   return vector<double>(1, sol.rhs);
}

vector<double> exact_solution(const vector<double> &coords) {
   ExactSolution sol;
   sol.calculate(exact_sol_id, coords);
   return vector<double>(1, sol.u);
}


void read_command_line_params(int argc, char **argv, int *num_levels) {
   if(argc == 1 + 1) {
      *num_levels = atoi(argv[1]);
   }
   else {
      if ( mpi_rank == 0 ) {
         printf(" Usage: mpirun -np X ./p4est_bddcml NLEVELS\n");
      }
      exit(0);
   }
}

void run(const P4estClass &p4est_class, int num_levels, IntegralResults *integral_results)
{
   print_rank = 0;

   Def::d()->init(num_dim, order, physicsType, p4est_class);

   BddcmlGeneralParams general_params;
   //general_params.just_direct_solve_int = 1;
   BddcmlKrylovParams krylov_params;
   BddcmlPreconditionerParams preconditioner_params;

   ProblemDimensions problem_dims(Def::d()->num_dim, physicsType, p4est_class);
   BddcmlSolver bddcml_solver(problem_dims, general_params, krylov_params,
                              preconditioner_params, p4est_class, num_levels);

   ReferenceElement ref_elem(num_dim, order);
   IntegrationMesh integration_mesh(p4est_class);
   NodalElementMesh nodal_mesh(physicsType, Def::d()->num_components, integration_mesh, ref_elem, p4est_class);

   vector<double> sols(problem_dims.n_subdom_dofs, 0.0);
   DiscreteSystem discrete_system(problem_dims, MatrixType::SPD);
   discrete_system.assemble(p4est_class, integration_mesh, nodal_mesh, problem_dims, &rhs_fn);

   bddcml_solver.solve(nodal_mesh, discrete_system, exact_solution, &sols);

   VtkOutput vtk(p4est_class, nodal_mesh, sols);
   vtk.output_in_corners("out_corners");
   vtk.output_in_nodes("out_nodes");
   vtk.output_exact_sol_in_nodes("out_exact", exact_solution);

   Integrator integrator(p4est_class, nodal_mesh, ref_elem, sols);
   integral_results->l2_norm = integrator.l2_norm(norm_order);
   integral_results->l2_error = integrator.l2_error(norm_order, exact_solution);
   PPP cout << "**************************************" << endl;
   PPP cout << "L2 norm  " << integral_results->l2_norm << "L2 error " << integral_results->l2_error
            << ", L2 relative error " << integral_results->l2_rel_error() << endl;
   PPP cout << "**************************************" << endl;
}

int main (int argc, char **argv)
{
   int mpiret = sc_MPI_Init (&argc, &argv);
   SC_CHECK_MPI (mpiret);
   mpiret = sc_MPI_Comm_rank(mpicomm, &mpi_rank);
   mpiret = sc_MPI_Comm_size(mpicomm, &mpi_size);

   int num_levels;
   read_command_line_params(argc, argv, &num_levels);

   P4estClass *p4est_class = P4estClass::create(num_dim, order, mpicomm);

   p4est_class->refine_and_partition(1, RefineType::UNIFORM);

   ofstream conv_file("conv.txt");
   IntegralResults integral_results;
   for(int i = 0; i < 5; i++) {
      p4est_class->refine_and_partition(1, RefineType::UNIFORM);
      run(*p4est_class, num_levels, &integral_results);
      conv_file << integral_results.l2_norm << ", " << integral_results.l2_rel_error() << endl;
   }

   delete p4est_class;

   sc_finalize ();
   assert(get_num_allocations() == 0);
   mpiret = sc_MPI_Finalize ();
   SC_CHECK_MPI (mpiret);
   return 0;
}
