#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <assert.h>
#include <math.h>

#include "bddcml/bddcml_solver.h"
#include "p4est/my_p4est_interface.h"
#include "assemble.h"
#include "integration_cell.h"
#include "element.h"
#include "shapefun.h"
#include "vtk_output.h"
#include "integral.h"

using namespace std;

const int num_dim = 2;
const int order = 1;
const int norm_order = 2 * order;
const PhysicsType physicsType = PhysicsType::LAPLACE;

const double slope = 60;
const double center[3] = {1.25, -0.25, -0.25 };

vector<double> rhs_fn(vector<double> coords) {

   double rr = 0.0;
   for(int i = 0; i < Def::d()->num_dim; i++) {
      coords[i] -= center[i];
      rr += coords[i] * coords[i];
   }
   double r = sqrt(rr);

//   double u = atan(slope * (r - M_PI / 3.));

//   double t_u_1 = r * (slope*slope * pow(r - M_PI/3, 2) + 1);
//   double dudx = slope * coords[0] / t_u_1;
//   double dudy = slope * coords[1] / t_u_1;

   double t_u_2 = (pow(M_PI - 3.0 * r, 2) * slope * slope + 9.0);
   double laplace = 27.0 * 2.0 * rr * (M_PI - 3.0*r) * pow(slope, 3.0) / (pow(t_u_2,2) * rr) - 9.0 * rr * slope / (t_u_2 * pow(r,3.0)) + 9.0 * slope / (t_u_2 * r);

   return {-laplace};
}

void exact_solution(const vector<double> &coords, vector<double> *result) {
   double rr = 0.0;
   vector<double> my_coords(coords);
   for(int i = 0; i < Def::d()->num_dim; i++) {
      my_coords[i] -= center[i];
      rr += my_coords[i] * my_coords[i];
   }
   double r = sqrt(rr);

   double u = atan(slope * (r - M_PI / 3.));

   *result = {u};
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

void run(const P4estClass &p4est_class, int num_levels)
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

   bddcml_solver.solve(nodal_mesh, discrete_system, &sols);

   VtkOutput vtk(p4est_class, nodal_mesh, sols);
   vtk.output_in_corners("out_corners");
   vtk.output_in_nodes("out_nodes");

   vector<double> l2errors;
   Integrator integrator(p4est_class, nodal_mesh, ref_elem, sols);
   double l2_norm = integrator.l2_norm(norm_order);
   double l2_error = integrator.l2_error(norm_order, exact_solution, &l2errors);
   PPP cout << "**************************************" << endl;
   PPP cout << "L2 norm  " << l2_norm << endl;
   PPP cout << "L2 error " << l2_error << endl;
   PPP cout << PrintVec<double>(l2errors) << endl;
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

   p4est_class->refine_and_partition(2, RefineType::UNIFORM);

   for(int i = 0; i < 5; i++) {
      p4est_class->refine_and_partition(1, RefineType::UNIFORM);
      run(*p4est_class, num_levels);
   }

   delete p4est_class;

   sc_finalize ();
   assert(get_num_allocations() == 0);
   mpiret = sc_MPI_Finalize ();
   SC_CHECK_MPI (mpiret);
   return 0;
}