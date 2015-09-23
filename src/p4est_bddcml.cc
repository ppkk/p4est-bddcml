#ifndef P4_TO_P8
#include <p4est_bits.h>
#include <p4est_ghost.h>
#include <p4est_lnodes.h>
#include <p4est_vtk.h>
#else
#include <p8est_bits.h>
#include <p8est_ghost.h>
#include <p8est_lnodes.h>
#include <p8est_vtk.h>
#endif

#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <assert.h>

#include "helpers.h"
#include "bddcml_structs.h"
#include "p4est_common.h"
#include "p4est_bddcml_interaction.h"

const int degree = 1;
const PhysicsType physicsType = LAPLACE;

int main (int argc, char **argv)
{
   init_corner_to_hanging();

   int mpiret = sc_MPI_Init (&argc, &argv);
   SC_CHECK_MPI (mpiret);
   sc_MPI_Comm mpicomm = sc_MPI_COMM_WORLD;
   mpiret = sc_MPI_Comm_rank(mpicomm, &mpi_rank);
   mpiret = sc_MPI_Comm_size(mpicomm, &mpi_size);

   /* These functions are optional.  If called they store the MPI rank as a
   * static variable so subsequent global p4est log messages are only issued
   * from processor zero.  Here we turn off most of the logging; see sc.h. */
   sc_init (mpicomm, 1, 1, NULL, SC_LP_ESSENTIAL);
   p4est_init (NULL, SC_LP_PRODUCTION);  /* SC_LP_ERROR for silence. */
   P4EST_GLOBAL_PRODUCTIONF
         ("This is the p4est %dD demo example/steps/%s_step4\n",
          P4EST_DIM, P4EST_STRING);

   // WARNING: integration is based on the fact, that the domain is unit square or cube!
   // WARNING: if the domain is changed, so has to be the integration!
   // TODO: do it properly
#ifndef P4_TO_P8
   p4est_connectivity_t *conn = p4est_connectivity_new_unitsquare ();
#else
   p4est_connectivity_t *conn = p8est_connectivity_new_unitcube ();
#endif

   BddcmlLevelInfo level_info;
   // Number of elements in an edge of a subdomain and number of subdomains in an edge of the unit cube
   if(argc == 1 + 1) {
      level_info.nlevels = atoi(argv[1]);
   }
   else {
      if ( mpi_rank == 0 ) {
         printf(" Usage: mpirun -np X ./p4est_bddcml NLEVELS\n");
      }
      exit(0);
   }

   // number of subdomains == mpi_size
   init_levels(mpi_size, &level_info);

   BddcmlGeneralParams general_params;
   set_implicit_general_params(&general_params);
   //general_params.just_direct_solve_int = 1;

   BddcmlKrylovParams krylov_params;
   set_implicit_krylov_params(&krylov_params);

   BddcmlPreconditionerParams preconditioner_params;
   set_implicit_preconditioner_params(&preconditioner_params);

   p4est_t *p4est = p4est_new (mpicomm, conn, 0, NULL, NULL);

   refine_and_partition(p4est, 1, refine_uniform);
   refine_and_partition(p4est, 3, refine_circle);
   refine_and_partition(p4est, 4, refine_square);
   refine_and_partition(p4est, 0, refine_point);
   refine_and_partition(p4est, 0, refine_diagonal);

   /* Create the ghost layer to learn about parallel neighbors. */
   p4est_ghost_t *ghost = p4est_ghost_new (p4est, P4EST_CONNECT_FULL);

   /* Create a node numbering for continuous linear finite elements. */
   p4est_lnodes_t *lnodes = p4est_lnodes_new (p4est, ghost, degree);

   /* Destroy the ghost structure -- no longer needed after node creation. */
   p4est_ghost_destroy (ghost);
   ghost = NULL;

   BddcmlDimensions subdomain_dims, global_dims;
   BddcmlMesh mesh;

   int print_rank_l = 3;

   // todo: using MPI Bcast in the following, should be possible to do without
   prepare_dimmensions(p4est, lnodes, physicsType, &subdomain_dims, &global_dims, mpicomm);

   //print_p4est_mesh(p4est, lnodes, print_rank_l);

   // TODO: elem_volume is correct only when the mesh is obtained by refinements
   // from a UNIT SQUARE/CUBE

   prepare_subdomain_mesh(p4est, lnodes, &subdomain_dims, &mesh);
   //print_bddcml_mesh(&mesh, print_rank_l);

   BddcmlFemSpace femsp;
   prepare_subdomain_fem_space(&mesh, &femsp);
   //print_bddcml_fem_space(&femsp, &mesh, print_rank_l);

   plot_solution(p4est, lnodes, NULL, NULL, NULL);

   print_rank = print_rank_l;
   print_basic_properties(&global_dims, mpi_size, &level_info, &krylov_params);
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
   int ndof_per_element = P4EST_CHILDREN;
   // how much space the upper triangle of the element matrix occupies
   int lelm = ndof_per_element * (ndof_per_element + 1) / 2;

   MatrixType matrix_type = SPD;
   // todo: do it properly
   const int extra_space_for_hanging_nodes = 4 * (matrix_type == GENERAL ? 2 : 1);
   allocate_sparse_matrix(extra_space_for_hanging_nodes * subdomain_dims.n_elems*lelm, matrix_type, &matrix);
   zero_matrix(&matrix);

   assemble_matrix_rhs(lnodes, &mesh, &femsp, &matrix, &rhss);
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
                                     subdomain_idx, &mesh, &femsp,
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


   plot_solution(p4est, lnodes, sols.val, NULL, NULL); //uexact_eval, NULL);


   free_mesh(&mesh);
   free_fem_space(&femsp);

   free_real_array(&rhss);
   free_real_array(&sols);
   free_sparse_matrix(&matrix);

   assert(get_num_allocations() == 0);

   /* Destroy the p4est and the connectivity structure. */
   p4est_lnodes_destroy (lnodes);
   p4est_destroy(p4est);
   p4est_connectivity_destroy(conn);

   /* Verify that allocations internal to p4est and sc do not leak memory.
   * This should be called if sc_init () has been called earlier. */
   sc_finalize ();


   /* This is standard MPI programs.  Without --enable-mpi, this is a dummy. */
   mpiret = sc_MPI_Finalize ();
   SC_CHECK_MPI (mpiret);
   return 0;
}
