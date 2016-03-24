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

// refinements: uniform, circle, square
// 2D elasticity on 4 procs gives 11 PCG iterations and condition number 0.433169186E+01, refs: 4, 3, 3
// 3D elasticity on 4 procs gives 19 iterations and condition number 0.190076921E+02, refs: 2, 3, 3
// 2D - anselm, refs: 11, 10, 6
// 3D - anselm, refs: 6, 5, 5

using namespace std;

// ugly global settings
int num_dim = -1;
int order = -1;
int norm_order = -1;
PhysicsType physicsType;

// 0 for poisson, 1 for elasticity
int use_corner_constraints;

const bool vtk_output = false;

vector<double> rhs_fn(vector<double>)
{
   if(physicsType == PhysicsType::LAPLACE)
      return {1};
   else if (physicsType == PhysicsType::ELASTICITY)
      if(num_dim == 2)
         return {0,-1e5};
      else
         return {0,0,-1e5};
   else
      assert(0);
}

const double slope = 60;
const double center[3] = {1.25, -0.25, -0.25 };

//vector<double> rhs_fn(vector<double> coords) {
//   double value = 0;
//   for(int dim = 0; dim < Def::d()->num_dim; dim++) {
//      double add = -2;
//      for(int dim_add = 0; dim_add < Def::d()->num_dim; dim_add++) {
//         if(dim != dim_add)
//         add *= coords[dim_add] * (1-coords[dim_add]);
//      }
//      value += add;
//   }

//   //1
//   double laplace = value;

////   double rr = 0.0;
////   for(int i = 0; i < Def::d()->num_dim; i++) {
////      coords[i] -= center[i];
////      rr += coords[i] * coords[i];
////   }
////   double r = sqrt(rr);

////   double u = atan(slope * (r - M_PI / 3.));

////   double t_u_1 = r * (slope*slope * pow(r - M_PI/3, 2) + 1);
////   double dudx = slope * coords[0] / t_u_1;
////   double dudy = slope * coords[1] / t_u_1;

////   double t_u_2 = (pow(M_PI - 3.0 * r, 2) * slope * slope + 9.0);

////   double laplace = 27.0 * 2.0 * rr * (M_PI - 3.0*r) * pow(slope, 3.0) / (pow(t_u_2,2) * rr) - 9.0 * rr * slope / (t_u_2 * pow(r,3.0)) + 9.0 * slope / (t_u_2 * r);


//   return {-laplace};
//}

vector<double> exact_solution(const vector<double> &coords) {

   //   //1
   double value = 1.;
   for(int dim = 0; dim < Def::d()->num_dim; dim++)
      value *= coords[dim] * (1-coords[dim]);
   return vector<double>(1, value);
   //   //2
   //   //*result = {(1-x)*x};



//   double rr = 0.0;
//   vector<double> my_coords(coords);
//   for(int i = 0; i < Def::d()->num_dim; i++) {
//      my_coords[i] -= center[i];
//      rr += my_coords[i] * my_coords[i];
//   }
//   double r = sqrt(rr);

//   double u = atan(slope * (r - M_PI / 3.));

//   *result = {u};
}



Parameters params(1e10, 0.33);

void read_command_line_params(int argc, char **argv, int *num_levels, int *unif, int* circle, int* square, int* extra_unif) {
   if(argc == 1 + 8) {
      num_dim = atoi(argv[1]);
      order = atoi(argv[2]);
      *num_levels = atoi(argv[3]);
      norm_order = 2 * order;
      if(strcmp(argv[4], "laplace") == 0) {
         physicsType = PhysicsType::LAPLACE;
         use_corner_constraints = 0;
      }
      else if(strcmp(argv[4], "elasticity") == 0) {
         physicsType = PhysicsType::ELASTICITY;
         use_corner_constraints = 1;
      }
      else
         assert(0);
      *unif = atoi(argv[5]);
      *circle = atoi(argv[6]);
      *square = atoi(argv[7]);
      *extra_unif = atoi(argv[8]);

      if ( mpi_rank == 0 ) {
         printf("NUM DIMENSIONS %d\n", num_dim);
         printf("ELEMENT ORDER %d\n", order);
         printf("NUM LEVELS %d\n", *num_levels);
         if(physicsType == PhysicsType::LAPLACE)
            printf("PHYSICS TYPE LAPLACE\n");
         else if(physicsType == PhysicsType::ELASTICITY)
            printf("PHYSICS TYPE ELASTICITY\n");
         else
            assert(0);
         printf("USE CORNER CONSTRAINTS %d\n", use_corner_constraints);
         printf("REFINEMENTS %d UNIFORM, %d CIRCLE, %d SQUARE, %d EXTRA UNIFORM\n", *unif, *circle, *square, *extra_unif);
      }
      MPI_Barrier(mpicomm);
   }
   else {
      if ( mpi_rank == 0 ) {
         printf(" Usage: mpirun -np X ./p4est_bddcml DIM ORDER NLEVELS PHYSICS UNIF_REF CIRCLE_REF SQUARE_REF EXTRA_UNIF_REF\n");
      }
      exit(0);
   }
}


void run(const P4estClass &p4est_class, int num_levels)
{
   clock_t assemble_begin = clock();

   print_rank = 0;

   Def::d()->init(num_dim, order, physicsType, p4est_class);

   BddcmlGeneralParams general_params;
   //general_params.just_direct_solve_int = 1;
   BddcmlKrylovParams krylov_params;
   BddcmlPreconditionerParams preconditioner_params;
   preconditioner_params.use_corner_constraints = use_corner_constraints;

   ProblemDimensions problem_dims(Def::d()->num_dim, physicsType, p4est_class);
   BddcmlSolver bddcml_solver(problem_dims, general_params, krylov_params,
                              preconditioner_params, p4est_class, num_levels);

   ReferenceElement ref_elem(num_dim, order);
   IntegrationMesh integration_mesh(p4est_class);
   NodalElementMesh nodal_mesh(physicsType, Def::d()->num_components, integration_mesh, ref_elem, p4est_class);

   vector<double> sols(problem_dims.n_subdom_dofs, 0.0);
   DiscreteSystem discrete_system(problem_dims, MatrixType::SPD);
   discrete_system.assemble(p4est_class, integration_mesh, nodal_mesh, problem_dims, &rhs_fn, params);

   clock_t assemble_end = clock();
   double assemble_time = double(assemble_end - assemble_begin) / CLOCKS_PER_SEC;
   PPP printf("Time of FEM assembly: %lf s\n", assemble_time);

   if(physicsType == PhysicsType::ELASTICITY)
      bddcml_solver.solve(nodal_mesh, discrete_system, nullptr, &sols);
   else
      bddcml_solver.solve(nodal_mesh, discrete_system, exact_solution, &sols);

   if(vtk_output)
   {
      VtkOutput vtk(p4est_class, nodal_mesh, sols);
      vtk.output_in_corners("out_corners");
      vtk.output_in_nodes("out_nodes");
   }

   Integrator integrator(p4est_class, nodal_mesh, ref_elem, sols);
   double l2_norm = integrator.l2_norm(norm_order);
   double l2_error = integrator.l2_error(norm_order, exact_solution);
   PPP cout << "**************************************" << endl;
   PPP cout << "L2 norm  " << l2_norm << endl;
   PPP cout << "L2 error " << l2_error << endl;
   PPP cout << "**************************************" << endl;
}

int main (int argc, char **argv)
{
   int num_levels, unif_ref, circle_ref, square_ref, extra_unif_ref;
   int mpiret = sc_MPI_Init (&argc, &argv);
   SC_CHECK_MPI (mpiret);
   mpiret = sc_MPI_Comm_rank(mpicomm, &mpi_rank);
   mpiret = sc_MPI_Comm_size(mpicomm, &mpi_size);

   if(((physicsType == PhysicsType::LAPLACE) && (use_corner_constraints == 1)) ||
      ((physicsType == PhysicsType::ELASTICITY) && (use_corner_constraints == 0)))
      printf("Warning: strange setting of use_corner_constraints!\n");

   read_command_line_params(argc, argv, &num_levels, &unif_ref, &circle_ref, &square_ref, &extra_unif_ref);

   P4estClass* p4est_class = P4estClass::create(num_dim, order, mpicomm);

   p4est_class->refine_and_partition(unif_ref, RefineType::UNIFORM);
   p4est_class->refine_and_partition(circle_ref, RefineType::CIRCLE);
   p4est_class->refine_and_partition(square_ref, RefineType::SQUARE);
   p4est_class->refine_and_partition(extra_unif_ref, RefineType::UNIFORM);

   run(*p4est_class, num_levels);

   delete p4est_class;

   sc_finalize ();
   assert(get_num_allocations() == 0);

   mpiret = sc_MPI_Finalize ();
   SC_CHECK_MPI (mpiret);
   return 0;
}
