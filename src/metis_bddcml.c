#ifndef P4_TO_P8
#include <p4est_bits.h>
#include <p4est_ghost.h>
#include <p4est_lnodes.h>
#include <p4est_vtk.h>
#include <p4est_iterate.h>
#else
#include <p8est_bits.h>
#include <p8est_ghost.h>
#include <p8est_lnodes.h>
#include <p8est_vtk.h>
#include <p8est_iterate.h>
#endif

#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <assert.h>

#include <metis.h>
//#include <GKlib/gk_proto.h>

#include "helpers.h"
#include "bddcml_structs.h"
#include "p4est_common.h"
#include "p4est_bddcml_interaction.h"

const int degree = 1;

#ifndef P4_TO_P8
   const int MAX_GRAPH_NEIGHBOURS = 16-4;
#else
   const int MAX_GRAPH_NEIGHBOURS = 64-8;
#endif

idx_t *graph_xadj, *graph_adjncy;

void add_graph_edge(idx_t first, idx_t second)
{
   //printf("adding graph edge %d -> %d\n", first, second);

   // we do not store a->a
   if(first == second)
      return;

   // check whether has been allready added (when considering connections through vertices, we try to add them multiple times
   for(int i = 0; i < graph_xadj[first]; i++)
   {
      if(graph_adjncy[first*MAX_GRAPH_NEIGHBOURS + i] == second)
         return;
   }

   assert(graph_xadj[first] < MAX_GRAPH_NEIGHBOURS);
   graph_adjncy[first*MAX_GRAPH_NEIGHBOURS + graph_xadj[first]] = second;
   graph_xadj[first]++;
}

void corners_action (p4est_iter_corner_info_t * info, void *user_data)
{
   sc_array_t *sides = &(info->sides);

#ifndef P4_TO_P8
   const int MAX_CORNER_ADJ = 4;
#else
   const int MAX_CORNER_ADJ = 8;
#endif

   int elem_num = sides->elem_count;
   int elem_idx[MAX_CORNER_ADJ];

   assert(sides->elem_count <= MAX_CORNER_ADJ);
   if(sides->elem_count > 1)
   {
      for(unsigned int side_idx = 0; side_idx < sides->elem_count; side_idx++)
      {
         p4est_iter_corner_side_t *side = p4est_iter_cside_array_index_int (sides, side_idx);
         assert(!side->is_ghost);
         elem_idx[side_idx] = side->quadid;
      }
   }

   for(int i = 0; i < elem_num; i++)
   {
      for(int j = 0; j < elem_num; j++)
      {
         add_graph_edge(elem_idx[i], elem_idx[j]);
      }
   }
}

void face_action (p4est_iter_face_info_t * info, void *user_data)
{
   //printf("face action \n");
   sc_array_t *sides = &(info->sides);

   int elem_nums[2]; // number of elements on each side (1-2, for reg/hanging)
   int elem_idx[2][2]; // element(s) on each side

   if(sides->elem_count > 1)
   {
      assert(sides->elem_count == 2);
      // inner face
      for(int side_idx = 0; side_idx < 2; side_idx++)
      {
         p4est_iter_face_side_t *side = p4est_iter_fside_array_index_int (sides, side_idx);
         // for more initial elements (trees) we need to convert to global indices
         assert(side->treeid == 0);

         if(side->is_hanging)
         {
            // not ready for distributed yet
            assert(!side->is.hanging.is_ghost[0]);
            assert(!side->is.hanging.is_ghost[1]);
            elem_nums[side_idx] = 2;
            elem_idx[side_idx][0] = side->is.hanging.quadid[0];
            elem_idx[side_idx][1] = side->is.hanging.quadid[1];
         }
         else
         {
            // not ready for distributed yet
            assert(!side->is.full.is_ghost);
            elem_nums[side_idx] = 1;
            elem_idx[side_idx][0] = side->is.full.quadid;
         }
      }
      for(int i = 0; i < elem_nums[0]; i++)
      {
         for(int j = 0; j < elem_nums[1]; j++)
         {
            add_graph_edge(elem_idx[0][i], elem_idx[1][j]);
            add_graph_edge(elem_idx[1][j], elem_idx[0][i]);
         }
      }
   }


}


void print_graph(int num_vert)
{
   printf("\n************ BEGIN GRAPH *********************\n");
   printf("xadj (%d) : ", num_vert+1);
   for(int i = 0; i <= num_vert; i++)
      printf("%d, ", graph_xadj[i]);
   printf("\n");
   printf("adjncy (%d) : ", MAX_GRAPH_NEIGHBOURS * num_vert);
   for(int i = 0; i < MAX_GRAPH_NEIGHBOURS * num_vert; i++)
      printf("%d, ", graph_adjncy[i]);
   printf("\n\n");
   for(int vert = 0; vert < num_vert; vert++)
   {
      printf("%d : (%d) -> ", vert, graph_xadj[vert]);
      for(int neighbour = graph_xadj[vert]; neighbour < graph_xadj[vert+1]; neighbour++)
      {
         printf("%d, ", graph_adjncy[neighbour]);
      }
      printf("\n");
   }
   printf("************   END GRAPH *********************\n\n");

}

enum ConnectType{
   CONNECT_CORNERS,
   CONNECT_FACES
};

void create_graph(p4est_t *p4est, p4est_lnodes_t *lnodes, idx_t nparts, mptype_et metis_algorithm, ConnectType connect, idx_t* part)
{
   // vertices of the graph are mesh elements
   idx_t num_vert = lnodes->num_local_elements;
   graph_xadj = (idx_t*) malloc((num_vert + 1) * sizeof(idx_t));
   memset(graph_xadj, 0, (num_vert + 1) * sizeof(idx_t));
   graph_adjncy = (idx_t*) malloc(MAX_GRAPH_NEIGHBOURS * num_vert * sizeof(idx_t));
   memset(graph_adjncy, 0, MAX_GRAPH_NEIGHBOURS * num_vert * sizeof(idx_t));


   if(connect == CONNECT_CORNERS)
   {
      p4est_iterate (p4est,                 /* the forest */
                     NULL,                 /* the ghost layer */
                     NULL,   /* the synchronized ghost data */
                     NULL,
                     NULL,
#ifdef P4_TO_P8
                     NULL,           // edges
#endif
                     corners_action);          // corners

   }
   else if(connect == CONNECT_FACES)
   {
      p4est_iterate (p4est,                 /* the forest */
                     NULL,                 /* the ghost layer */
                     NULL,   /* the synchronized ghost data */
                     NULL,
                     face_action,
#ifdef P4_TO_P8
                     NULL,           // edges
#endif
                     NULL);          // corners
   }
   else
   {
      assert(0);
   }

   // the following estimates might not work for more complex geometries (more than 4/8 elements meeting in 1 point)
   assert(p4est->trees->elem_count == 1);

   //print_graph(num_vert);

   // compress graph_adjncy
   int new_position = 0;
   for(int vert = 0; vert < num_vert; vert++)
   {
      for(int vert_edge = 0; vert_edge < graph_xadj[vert]; vert_edge++)
      {
         int old_position = MAX_GRAPH_NEIGHBOURS * vert + vert_edge;
         //printf("moving from %d to %d\n", old_position, new_position);
         graph_adjncy[new_position] = graph_adjncy[old_position];
         new_position++;
      }
   }

   //now graph_xadj has to be changed to running sum
   for(int i = num_vert; i > 0; i--)
      graph_xadj[i] = graph_xadj[i-1];
   graph_xadj[0] = 0;
   for(int i = 0; i < num_vert; i++)
      graph_xadj[i+1] += graph_xadj[i];

   //print_graph(num_vert);

   idx_t options[METIS_NOPTIONS];
   idx_t ncon = 1;
   idx_t objval;
   METIS_SetDefaultOptions(options);

//   gk_malloc_init();
//   real_t part_timer;
//   gk_startcputimer(part_timer);

   int status;
   switch (metis_algorithm) {
     case METIS_PTYPE_RB:
       status = METIS_PartGraphRecursive(&num_vert, &ncon, graph_xadj, graph_adjncy, NULL, NULL, NULL,
                    &nparts, NULL, NULL, options, &objval, part);
       break;

     case METIS_PTYPE_KWAY:
       status = METIS_PartGraphKway(&num_vert, &ncon, graph_xadj, graph_adjncy, NULL, NULL, NULL,
                                    &nparts, NULL, NULL, options, &objval, part);
       break;

   }

//   gk_stopcputimer(part_timer);

//   if (gk_GetCurMemoryUsed() != 0)
//     printf("***It seems that Metis did not free all of its memory! Report this.\n");
//   size_t maxmemory = gk_GetMaxMemoryUsed();
//   gk_malloc_cleanup(0);


   if (status == METIS_OK)
     printf("\n***Metis OK.\n");
   else if (status == METIS_ERROR_INPUT)
      printf("\n***Metis input error.\n");
   else if (status == METIS_ERROR_MEMORY)
      printf("\n***Metis memory error.\n");
   else if (status == METIS_ERROR)
      printf("\n***Metis some other error.\n");
   else
      assert(0);

//   for(int i = 0; i < num_vert; i++)
//      printf("%d, ", part[i]);
//   printf("\n");

   free(graph_xadj);
   free(graph_adjncy);
}

int main (int argc, char **argv)
{
   init_corner_to_hanging();

   int mpiret = sc_MPI_Init (&argc, &argv);
   SC_CHECK_MPI (mpiret);
   sc_MPI_Comm mpicomm_bddcml = sc_MPI_COMM_WORLD;
   sc_MPI_Comm mpicomm_p4est = sc_MPI_COMM_SELF;
   mpiret = sc_MPI_Comm_rank(mpicomm_bddcml, &mpi_rank);
   mpiret = sc_MPI_Comm_size(mpicomm_bddcml, &mpi_size);

   /* These functions are optional.  If called they store the MPI rank as a
   * static variable so subsequent global p4est log messages are only issued
   * from processor zero.  Here we turn off most of the logging; see sc.h. */
   sc_init (mpicomm_p4est, 1, 1, NULL, SC_LP_ESSENTIAL);
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

   p4est_t *p4est = p4est_new (mpicomm_p4est, conn, 0, NULL, NULL);

   refine_and_partition(p4est, 3, refine_uniform);
   refine_and_partition(p4est, 6, refine_circle);
   refine_and_partition(p4est, 8, refine_square);
   refine_and_partition(p4est, 0, refine_point);
   refine_and_partition(p4est, 0, refine_diagonal);

   /* Create the ghost layer to learn about parallel neighbors. */
   p4est_ghost_t *ghost = p4est_ghost_new (p4est, P4EST_CONNECT_FULL);

   /* Create a node numbering for continuous linear finite elements. */
   p4est_lnodes_t *lnodes = p4est_lnodes_new (p4est, ghost, degree);

   /* Destroy the ghost structure -- no longer needed after node creation. */
   p4est_ghost_destroy (ghost);
   ghost = NULL;

   idx_t *metis_part = NULL;
   if(mpi_rank == 0)
   {
      metis_part = (idx_t*) malloc(lnodes->num_local_elements * sizeof(idx_t));
      create_graph(p4est, lnodes, 10, METIS_PTYPE_KWAY, CONNECT_FACES, metis_part);
      plot_solution(p4est, lnodes, NULL, NULL, metis_part);
      free(metis_part);
      metis_part = NULL;
   }


   BddcmlDimensions subdomain_dims, global_dims;
   BddcmlMesh mesh;

   int print_rank_l = 2;

   // todo: using MPI Bcast in the following, should be possible to do without
   prepare_dimmensions(p4est, lnodes, &subdomain_dims, &global_dims, mpicomm_bddcml);

   print_p4est_mesh(p4est, lnodes, print_rank_l);

   // TODO: elem_volume is correct only when the mesh is obtained by refinements
   // from a UNIT SQUARE/CUBE
   real *element_volumes = (real*) malloc(subdomain_dims.n_elems * sizeof(real));

   prepare_subdomain_mesh(p4est, lnodes, &subdomain_dims, &mesh, element_volumes);
   //print_bddcml_mesh(&mesh, print_rank_l);

   BddcmlFemSpace femsp;
   prepare_subdomain_fem_space(&mesh, &femsp);
   //print_bddcml_fem_space(&femsp, &mesh, print_rank_l);

   //plot_solution(p4est, lnodes, NULL, NULL, metis_part);

   print_rank = print_rank_l;
   print_basic_properties(&global_dims, mpi_size, &level_info, &krylov_params);
   PPP printf("Initializing BDDCML ...");
   // tell me how much subdomains should I load
   level_info.nsub_loc_1 = -1;

   bddcml_init(&general_params, &level_info, mpicomm_bddcml);
   // should be 1 subdomain per processor
   assert(level_info.nsub_loc_1 == 1);

   mpiret = MPI_Barrier(mpicomm_bddcml);
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

   assemble_matrix_rhs(lnodes, &mesh, element_volumes, &femsp, &matrix, &rhss);
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


   mpiret = MPI_Barrier(mpicomm_bddcml);

   PPP printf("Preconditioner set-up ...\n");

   // PRECONDITIONER SETUP
   mpiret = MPI_Barrier(mpicomm_bddcml);
   // TODO: call time_start
   bddcml_setup_preconditioner(matrix.type, &preconditioner_params);

   mpiret = MPI_Barrier(mpicomm_bddcml);
   // TODO: call time_end(t_pc_setup)

   PPP printf("Preconditioner set-up done.\n");



   PPP printf("Calling Krylov method ...\n");

   mpiret = MPI_Barrier(mpicomm_bddcml);
   // TODO: call time_start
   // call with setting of iterative properties

   BddcmlConvergenceInfo convergence_info;

//   real normRn_sol, normRn2, normRn2_loc, normRn2_sub;
//   real normL2_sol, normL2_loc, normL2_sub;
//   real normLinf_sol, normLinf_loc;

   bddcml_solve(&krylov_params, &convergence_info, mpicomm_bddcml);
   mpiret = MPI_Barrier(mpicomm_bddcml);

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


   //plot_solution(p4est, lnodes, sols.val, NULL, NULL); //uexact_eval, NULL);


   free_mesh(&mesh);
   free_fem_space(&femsp);

   free_real_array(&rhss);
   free_real_array(&sols);
   free_sparse_matrix(&matrix);

   free(element_volumes);

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

