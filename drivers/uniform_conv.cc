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

const PhysicsType physicsType = PhysicsType::LAPLACE;
const ExactSolID exact_sol_id = ExactSolID::InternalLayer;
//const bool refine_adaptive = false;
const bool refine_adaptive = true;
const bool use_h1_seminorm_for_adapt = true;
//const int n_initial_uniform_refs = 4;
const int n_initial_uniform_refs = 3;
//const int max_ref_steps = 40;
const int max_ref_steps = 6;
const int max_glob_dofs = 1000000;
const double refine_fraction = 0.15;
const bool vtk_output = true;

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

vector<double> exact_gradient(const vector<double> &coords) {
   ExactSolution sol;
   sol.calculate(exact_sol_id, coords);
   return sol.grad;
}


void read_command_line_params(int argc, char **argv, int *num_dim, int *order, int *num_levels) {
   if(argc == 1 + 3) {
      *num_dim = atoi(argv[1]);
      *order = atoi(argv[2]);
      *num_levels = atoi(argv[3]);
   }
   else {
      if ( mpi_rank == 0 ) {
         printf(" Usage: mpirun -np X ./p4est_bddcml num_dims order, n_levels\n");
      }
      exit(0);
   }
}

void run(const P4estClass &p4est_class, int num_dim, int order, int num_levels, IntegralResults *integral_results)
{
   const int norm_order = 2 * order;
   print_rank = 0;

   clock_t assemble_begin = clock();

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

   PPP cout << "assembling... " << endl;
   vector<double> sols(problem_dims.n_subdom_dofs, 0.0);
   DiscreteSystem discrete_system(problem_dims, MatrixType::SPD);
   discrete_system.assemble(p4est_class, integration_mesh, nodal_mesh, problem_dims, &rhs_fn);

   clock_t assemble_end = clock();
   double assemble_time = double(assemble_end - assemble_begin) / CLOCKS_PER_SEC;
   PPP printf("Time of FEM assembly: %lf s\n", assemble_time);

   PPP cout << "done. Solving... " << endl;
   bddcml_solver.solve(nodal_mesh, discrete_system, exact_solution, &sols);

   if(vtk_output) {
      PPP cout << "done. vtk output..." << endl;
      VtkOutput vtk(p4est_class, nodal_mesh, sols);
      vtk.output_in_corners("out_corners");
      vtk.output_in_nodes("out_nodes");
      vtk.output_exact_sol_in_nodes("out_exact", exact_solution);
   }

   PPP cout << "done. Calculating integrals..." << endl;
   Integrator integrator(p4est_class, nodal_mesh, ref_elem, sols);
   integral_results->init(problem_dims);
   integral_results->l2_norm = integrator.l2_norm(norm_order);
   integral_results->h1_seminorm = integrator.h1_seminorm(norm_order);
   integral_results->l2_error = integrator.l2_error(norm_order, exact_solution, &integral_results->element_l2_error);
   integral_results->h1_semierror = integrator.h1_semierror(norm_order, exact_gradient, &integral_results->element_h1_semierror);
   PPP cout << "**************************************" << endl;
   PPP cout << "L2 norm  " << integral_results->l2_norm << "L2 error " << integral_results->l2_error
            << ", L2 relative error " << integral_results->l2_rel_error() << endl;
   PPP cout << "H1 seminorm  " << integral_results->h1_seminorm << "H1 semierror " << integral_results->h1_semierror
            << ", H1 relative semierror " << integral_results->h1_rel_semierror() << endl;
   PPP cout << "**************************************" << endl;
}

int main (int argc, char **argv)
{
   int mpiret = sc_MPI_Init (&argc, &argv);
   SC_CHECK_MPI (mpiret);
   mpiret = sc_MPI_Comm_rank(mpicomm, &mpi_rank);
   mpiret = sc_MPI_Comm_size(mpicomm, &mpi_size);

   int num_dim, order, num_levels;
   read_command_line_params(argc, argv, &num_dim, &order, &num_levels);

   P4estClass *p4est_class = P4estClass::create(num_dim, order, mpicomm);

   p4est_class->refine_and_partition(n_initial_uniform_refs, RefineType::UNIFORM);

   ofstream conv_file("conv.txt");
   IntegralResults integral_results;
   for(int i = 0; i < max_ref_steps ; i++) {
      run(*p4est_class, num_dim, order, num_levels, &integral_results);
      conv_file << integral_results.n_glob_elems << ", " << integral_results.n_glob_dofs << ", " <<
                   integral_results.l2_norm << ", " << integral_results.l2_rel_error() << ", " <<
                   integral_results.h1_seminorm << ", " << integral_results.h1_rel_semierror() << endl;

      if(integral_results.n_glob_dofs > max_glob_dofs)
         break;

      if(refine_adaptive) {
         if(use_h1_seminorm_for_adapt)
            p4est_class->refine_and_partition(integral_results.element_h1_semierror,
                                              refine_fraction * integral_results.n_glob_elems);
         else
            p4est_class->refine_and_partition(integral_results.element_l2_error,
                                              refine_fraction * integral_results.n_glob_elems);
      }
      else {
         p4est_class->refine_and_partition(1, RefineType::UNIFORM);
      }
   }

   delete p4est_class;

   sc_finalize ();
   assert(get_num_allocations() == 0);
   mpiret = sc_MPI_Finalize ();
   SC_CHECK_MPI (mpiret);
   return 0;
}
