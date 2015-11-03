#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <assert.h>

#include "bddcml/bddcml_mesh.h"
#include "bddcml/bddcml_femspace.h"
#include "p4est/my_p4est_interface.h"
#include "assemble.h"
#include "integration_cell.h"
#include "element.h"
#include "shapefun.h"
#include "vtk_output.h"

using namespace std;

const int num_dim = 2;
const int order = 3;
const PhysicsType physicsType = PhysicsType::LAPLACE;

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

Parameters params(1e10, 0.33);

sc_MPI_Comm mpicomm = sc_MPI_COMM_WORLD;

void run(int argc, char **argv)
{
   int mpiret = 0, num_levels;

   P4estClass* p4est_class = P4estClass::create(num_dim, order, mpicomm);
   Def::d()->init(num_dim, order, physicsType, p4est_class);

//   p4est_class->refine_and_partition(2, RefineType::UNIFORM);
//   p4est_class->refine_and_partition(2, RefineType::SQUARE);


   // 2D
//   p4est_class->refine_and_partition(4, RefineType::UNIFORM);
//   p4est_class->refine_and_partition(5, RefineType::CIRCLE);
//   p4est_class->refine_and_partition(6, RefineType::SQUARE);

   // 2D elasticity on 4 procs gives 11 PCG iterations and condition number 0.433169186E+01
   p4est_class->refine_and_partition(4, RefineType::UNIFORM);
   p4est_class->refine_and_partition(3, RefineType::CIRCLE);
   p4est_class->refine_and_partition(3, RefineType::SQUARE);

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

   // number of subdomains == mpi_size
   BddcmlLevelInfo level_info(num_levels, mpi_size);
   BddcmlGeneralParams general_params;
   //general_params.just_direct_solve_int = 1;
   BddcmlKrylovParams krylov_params;
   BddcmlPreconditionerParams preconditioner_params;

   BddcmlDimensions subdomain_dims(Def::d()->num_dim, physicsType);
   BddcmlDimensions global_dims(Def::d()->num_dim, physicsType);

   int print_rank_l = 3;

   // todo: using MPI Bcast in the following, should be possible to do without
   p4est_class->prepare_dimmensions(&subdomain_dims, &global_dims);
   //print_p4est_mesh(p4est, lnodes, print_rank_l);

   ReferenceElement ref_elem(num_dim, order);

   IntegrationMesh integration_mesh;
   p4est_class->prepare_integration_mesh(&integration_mesh);

   NodalElementMesh nodal_mesh;
   p4est_class->prepare_nodal_mesh(Def::d()->num_components, integration_mesh, ref_elem, &nodal_mesh);

   BddcmlMesh bddcml_mesh(&subdomain_dims);
   p4est_class->prepare_bddcml_mesh_global_mappings(&bddcml_mesh);
//   p4est_class->prepare_bddcml_mesh_nodes_old(&bddcml_mesh);
   bddcml_mesh.fill_nodes_info(*p4est_class, nodal_mesh);

   //print_bddcml_mesh(&mesh, print_rank_l);

   BddcmlFemSpace femsp(&bddcml_mesh);
   femsp.prepare_subdomain_fem_space(physicsType);
   //print_bddcml_fem_space(&femsp, &mesh, print_rank_l);

   print_rank = print_rank_l;
   print_basic_properties(global_dims, mpi_size, level_info, krylov_params);
   PPP printf("Initializing BDDCML ...");
   // tell me how much subdomains should I load
   level_info.nsub_loc_1 = -1;

   bddcml_init(&general_params, &level_info, mpicomm);
   // should be 1 subdomain per processor
   assert(level_info.nsub_loc_1 == 1);

   mpiret = MPI_Barrier(mpicomm);
   PPP printf("Initializing BDDCML done.\n");

   RealArray rhss;
   allocate_real_array(subdomain_dims.n_dofs, &rhss);
   zero_real_array(&rhss);

   int is_rhs_complete = 0;

   RealArray sols;
   allocate_real_array(subdomain_dims.n_dofs, &sols);
   zero_real_array(&sols);

   SparseMatrix matrix;
   int ndof_per_element = Def::d()->num_element_nodes * femsp.subdomain_dims->n_node_dofs;
   // how much space the upper triangle of the element matrix occupies
   int lelm = ndof_per_element * (ndof_per_element + 1) / 2;

   MatrixType matrix_type = MatrixType::SPD;
   // todo: do it properly
   const int extra_space_for_hanging_nodes = 4 * (matrix_type == MatrixType::GENERAL ? 2 : 1);
   allocate_sparse_matrix(extra_space_for_hanging_nodes * subdomain_dims.n_elems * lelm, matrix_type, &matrix);
   zero_matrix(&matrix);

   assemble_matrix_rhs(*p4est_class, integration_mesh, bddcml_mesh, femsp, &matrix, &rhss, &rhs_fn, params);
   //print_complete_matrix_rhs(&femsp, &global_dims, &matrix, &rhss, mpicomm);

   // user constraints - not really used here
   Real2DArray user_constraints;
   allocate_real_2D_array(0, 0, &user_constraints);

   // data for elements - not really used here
   Real2DArray element_data;
   allocate_real_2D_array(0, 0, &element_data);

   // data for dofs - not really used here
   RealArray dof_data;
   allocate_real_array(0, &dof_data);


   PPP printf("Loading data ...\n");

   int subdomain_idx = mpi_rank;
   bddcml_upload_subdomain_data(&global_dims, &subdomain_dims,
                                     subdomain_idx, &bddcml_mesh, &femsp,
                                     &rhss, is_rhs_complete, &sols, &matrix,
                                     &user_constraints, &element_data,
                                     &dof_data, &preconditioner_params);

   PPP printf("Loading data done.\n");


   mpiret = MPI_Barrier(mpicomm);

   PPP printf("Preconditioner set-up ...\n");

   // PRECONDITIONER SETUP
   mpiret = MPI_Barrier(mpicomm);
   // TODO: call time_start
   bddcml_setup_preconditioner(matrix.type, &preconditioner_params);

   mpiret = MPI_Barrier(mpicomm);
   // TODO: call time_end(t_pc_setup)

   PPP printf("Preconditioner set-up done.\n");

   PPP printf("Calling Krylov method ...\n");

   mpiret = MPI_Barrier(mpicomm);
   // TODO: call time_start
   // call with setting of iterative properties

   BddcmlConvergenceInfo convergence_info;

//   real normRn_sol, normRn2, normRn2_loc, normRn2_sub;
//   real normL2_sol, normL2_loc, normL2_sub;
//   real normLinf_sol, normLinf_loc;

   bddcml_solve(&krylov_params, &convergence_info, mpicomm);
   mpiret = MPI_Barrier(mpicomm);

   // TODO: call time_end(t_krylov)

   PPP printf("Krylov method done.\n");

   PPP printf(" Output of PCG: ==============\n");
   PPP printf(" Number of iterations: %d\n", convergence_info.num_iter);
   PPP printf(" Convergence reason:   %d\n", convergence_info.converged_reason);
   if ( convergence_info.condition_number >= 0. ) {
      PPP printf(" Condition number: %lf\n", convergence_info.condition_number);
   }
   PPP printf(" =============================\n");


   bddcml_download_local_solution(subdomain_idx, &sols);

   VtkOutput vtk(*p4est_class, nodal_mesh, sols);
   vtk.output("out");

   free_real_array(&rhss);
   free_real_array(&sols);
   free_sparse_matrix(&matrix);

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
