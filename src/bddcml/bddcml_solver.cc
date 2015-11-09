#include "definitions.h"
#include "p4est/my_p4est_interface.h"
#include "element.h"
#include "assemble.h"
#include "bddcml/bddcml_mesh.h"
#include "bddcml/bddcml_femspace.h"
#include "bddcml/bddcml_solver.h"

class BddcmlKrylovParams;

using namespace std;

BddcmlSolver::BddcmlSolver(ProblemDimensions &problem_dims, BddcmlGeneralParams &general_params,
                           BddcmlKrylovParams &krylov_params, BddcmlPreconditionerParams &preconditioner_params,
                           const P4estClass &p4est_class, int num_levels):
      problem_dims(problem_dims), general_params(general_params), krylov_params(krylov_params),
      preconditioner_params(preconditioner_params), p4est_class(p4est_class), num_levels(num_levels) {

}


void BddcmlSolver::solve(const NodalElementMesh &nodal_mesh, DiscreteSystem &system, exact_fn dirichlet_bc_exact, vector<double> *sols) {
   int mpiret = 0;
   // number of subdomains == mpi_size
   BddcmlLevelInfo level_info(num_levels, mpi_size);

   BddcmlMesh bddcml_mesh(problem_dims);
   p4est_class.prepare_bddcml_mesh_global_mappings(&bddcml_mesh);
//   p4est_class->prepare_bddcml_mesh_nodes_old(&bddcml_mesh);
   bddcml_mesh.fill_nodes_info(p4est_class, nodal_mesh);

   //print_bddcml_mesh(&mesh, print_rank_l);

   BddcmlFemSpace femsp(bddcml_mesh);
   femsp.prepare_subdomain_fem_space(nodal_mesh.physics_type, dirichlet_bc_exact);
   //print_bddcml_fem_space(&femsp, &mesh, print_rank_l);

   print_basic_properties(problem_dims, mpi_size, level_info, krylov_params);
   PPP printf("Initializing BDDCML ...");
   // tell me how much subdomains should I load
   level_info.nsub_loc_1 = -1;

   bddcml_init(&general_params, &level_info, mpicomm);
   // should be 1 subdomain per processor
   assert(level_info.nsub_loc_1 == 1);

   mpiret = MPI_Barrier(mpicomm);
   PPP printf("Initializing BDDCML done.\n");

   int is_rhs_complete = 0;

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
   bddcml_upload_subdomain_data(&problem_dims, subdomain_idx, &bddcml_mesh, &femsp,
                                     &system.rhss, is_rhs_complete, sols, &system.matrix,
                                     &user_constraints, &element_data,
                                     &dof_data, &preconditioner_params);

   PPP printf("Loading data done.\n");

   mpiret = MPI_Barrier(mpicomm);

   PPP printf("Preconditioner set-up ...\n");

   // PRECONDITIONER SETUP
   mpiret = MPI_Barrier(mpicomm);
   // TODO: call time_start
   bddcml_setup_preconditioner(system.matrix.type, &preconditioner_params);

   mpiret = MPI_Barrier(mpicomm);
   // TODO: call time_end(t_pc_setup)

   PPP printf("Preconditioner set-up done.\n");

   PPP printf("Calling Krylov method ...\n");

   mpiret = MPI_Barrier(mpicomm);
   // TODO: call time_start
   // call with setting of iterative properties

   BddcmlConvergenceInfo convergence_info;

   bddcml_solve(/*(BddcmlKrylovParams*)*/&krylov_params, &convergence_info, mpicomm);
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

   bddcml_download_local_solution(subdomain_idx, sols);
   bddcml_finalize();

   SC_CHECK_MPI (mpiret);
}
