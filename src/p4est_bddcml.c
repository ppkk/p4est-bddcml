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

#include <stdbool.h>
#include <assert.h>

#include "helpers.h"
#include "bddcml_structs.h"
#include "p4est_common.h"

const int degree = 1;

void prepare_dimmensions(p4est_t *p4est, p4est_lnodes_t *lnodes,
                         BddcmlDimensions *subdomain_dims, BddcmlDimensions *global_dims,
                         sc_MPI_Comm mpicomm)
{
#ifndef P4_TO_P8
   init_dimmensions(subdomain_dims, 2);
   init_dimmensions(global_dims, 2);
#else
   init_dimmensions(subdomain_dims, 3);
   init_dimmensions(global_dims, 3);
#endif
   subdomain_dims->n_nodes = lnodes->num_local_nodes;
   subdomain_dims->n_dofs  = lnodes->num_local_nodes;
   subdomain_dims->n_elems = lnodes->num_local_elements;

   int global_num_nodes;
   if(mpi_rank == mpi_size - 1)
   {
      global_num_nodes = lnodes->global_offset + lnodes->owned_count;
   }
   sc_MPI_Bcast(&global_num_nodes, 1, MPI_INT, mpi_size - 1, mpicomm);

   global_dims->n_nodes = global_num_nodes;
   global_dims->n_dofs = global_num_nodes;
   global_dims->n_elems = p4est->global_num_quadrants;
}

void prepare_subdomain_mesh(p4est_t *p4est, p4est_lnodes_t *lnodes, BddcmlDimensions *subdomain_dims, BddcmlMesh *mesh)
{
   double              vxyz[3];  /* We embed the 2D vertices into 3D space. */
   p4est_quadrant_t   sp, node;

   init_mesh(subdomain_dims, mesh);
   for(int lnode = 0; lnode < lnodes->num_local_nodes; lnode++)
   {
      mesh->node_global_map.val[lnode] = node_loc_to_glob(lnodes, lnode);
   }
   //p4est->global_first_quadrant
   /* Loop over local quadrants to apply the element matrices. */
   p4est_locidx_t quad_idx = 0;
   p4est_quadrant_t *quad;

   for_all_quads(p4est, quad_idx, quad)
   {
      mesh->elem_global_map.val[quad_idx] = (int)p4est->global_first_quadrant[p4est->mpirank] + quad_idx;
      mesh->num_nodes_of_elem.val[quad_idx] = P4EST_CHILDREN;

      for (int lnode = 0; lnode < P4EST_CHILDREN; ++lnode) {
         /* Cache some information on corner nodes. */
         p4est_locidx_t node_idx = lnodes->element_nodes[P4EST_CHILDREN * quad_idx + lnode];
         //      isboundary[i] = (bc == NULL ? 0 : bc[lni]);
         //       inloc[i] = !isboundary[i] ? in[lni] : 0.;
         mesh->elem_node_indices.val[P4EST_CHILDREN * quad_idx + lnode] = node_idx;
      }


      /* Figure out the hanging corners on this element, if any. */
      int hanging_corner[P4EST_CHILDREN];
      int anyhang = lnodes_decode2 (lnodes->face_code[quad_idx], hanging_corner);

      p4est_quadrant_t* parent;
      if (!anyhang) {
         parent = NULL;          /* Defensive programming. */
      }
      else {
         /* At least one node is hanging.  We need the parent quadrant to
           * find the location of the corresponding non-hanging node. */
         parent = &sp;
         p4est_quadrant_parent (quad, parent);
      }
      for (int lnode = 0; lnode < P4EST_CHILDREN; ++lnode) {
         p4est_locidx_t node_idx = lnodes->element_nodes[P4EST_CHILDREN * quad_idx + lnode];
         P4EST_ASSERT (node_idx >= 0 && node_idx < subdomain_dims->n_nodes);
         if (anyhang && hanging_corner[lnode] >= 0) {
            /* This node is hanging; access the referenced node instead. */
            p4est_quadrant_corner_node (parent, lnode, &node);
         }
         else {
            p4est_quadrant_corner_node (quad, lnode, &node);
         }

         /* Transform per-tree reference coordinates into physical space. */
         p4est_qcoord_to_vertex (p4est->connectivity, tt, node.x, node.y,         
#ifdef P4_TO_P8
                                 node.z,
#endif
                                 vxyz);

         mesh->coords.val[0][node_idx] = vxyz[0];
         mesh->coords.val[1][node_idx] = vxyz[1];
#ifdef P4_TO_P8
         mesh->coords.val[2][node_idx] = vxyz[2];
#endif



//         PPP printf("(%3.2lf, %3.2lf), ", vxyz[0], vxyz[1]);

      }
//      PPP printf("\n");

   }}}//for all quads

}


const double EPS = 1e-10;
bool real_equal(real a, real b)
{
   return fabs(a - b) < EPS;
}

void prepare_subdomain_fem_space(BddcmlMesh *mesh, BddcmlFemSpace *femsp)
{
   init_fem_space(mesh->subdomain_dims, femsp);
   mesh->subdomain_dims->n_dofs = mesh->subdomain_dims->n_nodes;
   for(int node = 0; node < mesh->subdomain_dims->n_dofs; node++)
   {
      femsp->node_num_dofs.val[node] = 1;
      femsp->dofs_global_map.val[node] = mesh->node_global_map.val[node];
      if((real_equal(mesh->coords.val[0][node], 0.0)) || (real_equal(mesh->coords.val[0][node], 1.0))
            || (real_equal(mesh->coords.val[0][node], 0.0)) || (real_equal(mesh->coords.val[0][node], 1.0))
            || (real_equal(mesh->coords.val[0][node], 0.0)) || (real_equal(mesh->coords.val[0][node], 1.0))
            )
      {
         femsp->fixs_code.val[node] = 1;
         femsp->fixs_values.val[node] = 0.0;
      }
      else
      {
         femsp->fixs_code.val[node] = 0;
      }

   }
}

void assemble_matrix(p4est_t * p4est, SparseMatrix *matrix)
{

}


static int refine_uniform (p4est_t * p4est, p4est_topidx_t which_tree,
                p4est_quadrant_t * quadrant)
{
   return 1;
}


int main (int argc, char **argv)
{
   int                 mpiret;
   sc_MPI_Comm         mpicomm;
   p4est_connectivity_t *conn;

   mpiret = sc_MPI_Init (&argc, &argv);
   SC_CHECK_MPI (mpiret);
   mpicomm = sc_MPI_COMM_WORLD;
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

#ifndef P4_TO_P8
   conn = p4est_connectivity_new_unitsquare ();
#else
   conn = p8est_connectivity_new_unitcube ();
#endif

   /* Create a forest that is not refined; it consists of the root octant.
   * The p4est_new_ext function can take a startlevel for a load-balanced
   * initial uniform refinement.  Here we refine adaptively instead. */
   int startlevel = 0;
   p4est_t *p4est = p4est_new (mpicomm, conn, 0, NULL, NULL);

   int endlevel = 3;
   for (int level = startlevel; level < endlevel; ++level) {
      p4est_refine (p4est, 0, refine_uniform, NULL);
      p4est_partition (p4est, 0, NULL);
   }
   if (startlevel < endlevel) {
      p4est_balance (p4est, P4EST_CONNECT_FULL, NULL);
      p4est_partition (p4est, 0, NULL);
   }


//   SparseMatrix matrix;
//   assemble_matrix(p4est, &matrix);


   /* Create the ghost layer to learn about parallel neighbors. */
   p4est_ghost_t *ghost = p4est_ghost_new (p4est, P4EST_CONNECT_FULL);

   /* Create a node numbering for continuous linear finite elements. */
   p4est_lnodes_t *lnodes = p4est_lnodes_new (p4est, ghost, degree);

   /* Destroy the ghost structure -- no longer needed after node creation. */
   p4est_ghost_destroy (ghost);
   ghost = NULL;

   BddcmlDimensions subdomain_dims, global_dims;
   BddcmlMesh mesh;

   int print_rank_l = 2;

   // todo: using MPI Bcast in the following, should be possible to do without
   prepare_dimmensions(p4est, lnodes, &subdomain_dims, &global_dims, mpicomm);

   print_p4est_mesh(p4est, lnodes, print_rank_l);
   prepare_subdomain_mesh(p4est, lnodes, &subdomain_dims, &mesh);
   print_bddcml_mesh(&mesh, print_rank_l);

   BddcmlFemSpace femsp;
   prepare_subdomain_fem_space(&mesh, &femsp);
   print_bddcml_fem_space(&femsp, &mesh, print_rank_l);

   plot_solution(p4est, lnodes, NULL, NULL);

   /* Destroy the p4est and the connectivity structure. */
   p4est_lnodes_destroy (lnodes);
   p4est_destroy(p4est);
   p4est_connectivity_destroy(conn);

   /* Verify that allocations internal to p4est and sc do not leak memory.
   * This should be called if sc_init () has been called earlier. */
   sc_finalize ();

   //*************************************************************************************
   //*************************************************************************************
   //    BDDCML part of the example
   //*************************************************************************************
   //*************************************************************************************

   BddcmlLevelInfo level_info;
   // Number of elements in an edge of a subdomain and number of subdomains in an edge of the unit cube
   if(argc == 1 + 1) {
      level_info.nlevels = atoi(argv[1]);
   }
   else {
      if ( mpi_rank == 0 ) {
         printf(" Usage: mpirun -np X ./p4est_bddcml NLEVELS");
      }
      exit(0);
   }

   // number of subdomains == mpi_size
   init_levels(mpi_size, &level_info);

   BddcmlGeneralParams general_params;
   set_implicit_general_params(&general_params);

   BddcmlKrylovParams krylov_params;
   set_implicit_krylov_params(&krylov_params);

   BddcmlPreconditionerParams preconditioner_params;
   set_implicit_preconditioner_params(&preconditioner_params);

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

   PPP printf("Loading data ...\n");

   bddcml_upload_subdomain_data(&global_dims, &subdomain_dims,
                                     mpi_rank, &mesh, &femsp,
                                     &rhss, is_rhs_complete, &sols, &matrix,
                                     &user_constraints, &element_data,
                                     &dof_data, &preconditioner_params);

   mpiret = MPI_Barrier(mpicomm);
   PPP printf("Loading data done.\n");


   /* This is standard MPI programs.  Without --enable-mpi, this is a dummy. */
   mpiret = sc_MPI_Finalize ();
   SC_CHECK_MPI (mpiret);
   return 0;
}
