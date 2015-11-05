#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <assert.h>
#include <math.h>

//todo: remove
#include "bddcml/bddcml_mesh.h"
#include "bddcml/bddcml_femspace.h"

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
const int order = 3;
const int norm_order = 2 * order;
const PhysicsType physicsType = PhysicsType::ELASTICITY;

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

void exact_solution(const vector<double> &coords, vector<double> *result) {

   //   //1
   double value = 1.;
   for(int dim = 0; dim < Def::d()->num_dim; dim++)
      value *= coords[dim] * (1-coords[dim]);
   *result = {value};
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


void run(int argc, char **argv)
{
   int mpiret = 0, num_levels;
   int print_rank_l = 0;
   print_rank = print_rank_l;

   P4estClass* p4est_class = P4estClass::create(num_dim, order, mpicomm);
   Def::d()->init(num_dim, order, physicsType, p4est_class);

   p4est_class->refine_and_partition(2, RefineType::UNIFORM);
   p4est_class->refine_and_partition(2, RefineType::CIRCLE);


   // 2D
//   p4est_class->refine_and_partition(4, RefineType::UNIFORM);
//   p4est_class->refine_and_partition(5, RefineType::CIRCLE);
//   p4est_class->refine_and_partition(6, RefineType::SQUARE);

   // 2D elasticity on 4 procs gives 11 PCG iterations and condition number 0.433169186E+01
//   p4est_class->refine_and_partition(4, RefineType::UNIFORM);
//   p4est_class->refine_and_partition(3, RefineType::CIRCLE);
//   p4est_class->refine_and_partition(3, RefineType::SQUARE);

   // 3D elasticity on 4 procs gives 19 iterations and condition number 0.190076921E+02
//   p4est_class->refine_and_partition(2, RefineType::UNIFORM);
//   p4est_class->refine_and_partition(3, RefineType::CIRCLE);
//   p4est_class->refine_and_partition(3, RefineType::SQUARE);

   // Number of elements in an edge of a subdomain and number of subdomains in an edge of the unit cube
   if(argc == 1 + 1) {
      num_levels = atoi(argv[1]);
   }
   else {
      if ( mpi_rank == 0 ) {
         printf(" Usage: mpirun -np X ./p4est_bddcml NLEVELS\n");
      }
      exit(0);
   }

   BddcmlGeneralParams general_params;
   //general_params.just_direct_solve_int = 1;
   BddcmlKrylovParams krylov_params;
   BddcmlPreconditionerParams preconditioner_params;

   ProblemDimensions subdomain_dims(Def::d()->num_dim, physicsType);
   ProblemDimensions global_dims(Def::d()->num_dim, physicsType);
   // todo: using MPI Bcast in the following, should be possible to do without
   p4est_class->prepare_dimmensions(&subdomain_dims, &global_dims);
   //print_p4est_mesh(p4est, lnodes, print_rank_l);

   BddcmlSolver bddcml_solver(subdomain_dims, global_dims, general_params, krylov_params,
                              preconditioner_params, *p4est_class, num_levels);

   ReferenceElement ref_elem(num_dim, order);

   IntegrationMesh integration_mesh;
   p4est_class->prepare_integration_mesh(&integration_mesh);

   NodalElementMesh nodal_mesh(physicsType);
   p4est_class->prepare_nodal_mesh(Def::d()->num_components, integration_mesh, ref_elem, &nodal_mesh);


   vector<double> sols(subdomain_dims.n_dofs, 0.0);
   DiscreteSystem discrete_system(subdomain_dims, MatrixType::SPD);
   discrete_system.assemble(*p4est_class, integration_mesh, nodal_mesh, subdomain_dims, &rhs_fn, params);

   bddcml_solver.solve(nodal_mesh, discrete_system, &sols);

   VtkOutput vtk(*p4est_class, nodal_mesh, sols);
   vtk.output_in_corners("out_corners");
   vtk.output_in_nodes("out_nodes");

   Integrator integrator(*p4est_class, nodal_mesh, ref_elem, sols);
   double l2_norm = integrator.l2_norm(norm_order);
   double l2_error = integrator.l2_error(norm_order, exact_solution);
   PPP cout << "**************************************" << endl;
   PPP cout << "L2 norm  " << l2_norm << endl;
   PPP cout << "L2 error " << l2_error << endl;
   PPP cout << "**************************************" << endl;

   delete p4est_class;

   SC_CHECK_MPI (mpiret);
}

int main (int argc, char **argv)
{
   int mpiret = sc_MPI_Init (&argc, &argv);
   SC_CHECK_MPI (mpiret);
   mpiret = sc_MPI_Comm_rank(mpicomm, &mpi_rank);
   mpiret = sc_MPI_Comm_size(mpicomm, &mpi_size);

//   /* These functions are optional.  If called they store the MPI rank as a
//   * static variable so subsequent global p4est log messages are only issued
//   * from processor zero.  Here we turn off most of the logging; see sc.h. */
//   sc_init (mpicomm, 1, 1, NULL, SC_LP_ESSENTIAL);
//   p4est_init (NULL, SC_LP_PRODUCTION);  /* SC_LP_ERROR for silence. */

   run(argc, argv);

   /* Verify that allocations internal to p4est and sc do not leak memory.
   * This should be called if sc_init () has been called earlier. */
   sc_finalize ();

   assert(get_num_allocations() == 0);

   /* This is standard MPI programs.  Without --enable-mpi, this is a dummy. */
   mpiret = sc_MPI_Finalize ();
   SC_CHECK_MPI (mpiret);
   return 0;
}
