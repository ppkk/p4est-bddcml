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

#include "helpers.h"
#include "bddcml_structs.h"
#include "p4est_common.h"

int mpi_rank;
int print_rank = 0;

void prepare_subdomain_mesh(p4est_t *p4est, p4est_lnodes_t *lnodes, BddcmlDimensions *subdomain_dims, BddcmlMesh *mesh)
{
   int                 anyhang, hanging_corner[P4EST_CHILDREN];
   int                 i;
   p4est_locidx_t      all_lni[P4EST_CHILDREN];

   double              vxyz[3];  /* We embed the 2D vertices into 3D space. */
   int8_t             *bc;
   p4est_topidx_t      tt;       /* Connectivity variables have this type. */
   p4est_locidx_t      k, q, Q;  /* Process-local counters have this type. */
   p4est_locidx_t      lni;      /* Node index relative to this processor. */
   p4est_tree_t       *tree;     /* Pointer to one octree */
   p4est_quadrant_t   *quad, *parent, sp, node;
   sc_array_t         *tquadrants;       /* Quadrant array for one tree */

   subdomain_dims->n_nodes = lnodes->num_local_nodes;
   subdomain_dims->n_dofs  = lnodes->num_local_nodes;
   subdomain_dims->n_elems = lnodes->num_local_elements;
#ifndef P4_TO_P8
   subdomain_dims->n_problem_dims = 2;
   subdomain_dims->n_mesh_dims = 2;
#else
   subdomain_dims->n_problem_dims = 3;
   subdomain_dims->n_mesh_dims = 3;
#endif

   init_mesh(subdomain_dims, mesh);

   /* Loop over local quadrants to apply the element matrices. */
   for (tt = p4est->first_local_tree, k = 0;
        tt <= p4est->last_local_tree; ++tt) {
      tree = p4est_tree_array_index (p4est->trees, tt);
      tquadrants = &tree->quadrants;
      Q = (p4est_locidx_t) tquadrants->elem_count;

      print_rank = 3;
      PPP printf("rank %d, elems %d, nodes %d (owned %d), global offset = %d\n", print_rank, lnodes->num_local_elements, lnodes->num_local_nodes, lnodes->owned_count, (int)lnodes->global_offset);
      for(i = 0; i < lnodes->num_local_nodes; i++)
      {
         PPP printf("loc %d -> glob %ld\n", i, node_loc_to_glob(lnodes, i));
      }
      for (q = 0; q < Q; ++q, ++k) {
         quad = p4est_quadrant_array_index (tquadrants, q);

         for (i = 0; i < P4EST_CHILDREN; ++i) {
            /* Cache some information on corner nodes. */
            lni = lnodes->element_nodes[P4EST_CHILDREN * k + i];
            //      isboundary[i] = (bc == NULL ? 0 : bc[lni]);
            //       inloc[i] = !isboundary[i] ? in[lni] : 0.;
            all_lni[i] = lni;
         }

         PPP printf("%d, %d, %d, %d", all_lni[0], all_lni[1], all_lni[2], all_lni[3]);
         //PPP printf("%d, %d, %d, %d", glob_all_lni[0], glob_all_lni[1], glob_all_lni[2], glob_all_lni[3]);


         /* Figure out the hanging corners on this element, if any. */
         anyhang = lnodes_decode2 (lnodes->face_code[k], hanging_corner);


         if (!anyhang) {
            parent = NULL;          /* Defensive programming. */
         }
         else {
            /* At least one node is hanging.  We need the parent quadrant to
           * find the location of the corresponding non-hanging node. */
            parent = &sp;
            p4est_quadrant_parent (quad, parent);
         }
         for (i = 0; i < P4EST_CHILDREN; ++i) {
            lni = lnodes->element_nodes[P4EST_CHILDREN * k + i];
            P4EST_ASSERT (lni >= 0 && lni < subdomain_dims->n_nodes);
            if (anyhang && hanging_corner[i] >= 0) {
               /* This node is hanging; access the referenced node instead. */
               p4est_quadrant_corner_node (parent, i, &node);
            }
            else {
               p4est_quadrant_corner_node (quad, i, &node);
            }

            /* Transform per-tree reference coordinates into physical space. */
            p4est_qcoord_to_vertex (p4est->connectivity, tt, node.x, node.y,
                        #ifdef P4_TO_P8
                                    node.z,
                        #endif
                                    vxyz);
            PPP printf("(%3.2lf, %3.2lf), ", vxyz[0], vxyz[1]);

         }
         PPP printf("\n");
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
   p4est_lnodes_t *lnodes = p4est_lnodes_new (p4est, ghost, 1);

   /* Destroy the ghost structure -- no longer needed after node creation. */
   p4est_ghost_destroy (ghost);
   ghost = NULL;

   BddcmlDimensions subdomain_dims;
   BddcmlMesh mesh;

   prepare_subdomain_mesh(p4est, lnodes, &subdomain_dims, &mesh);
   plot_solution(p4est, lnodes, NULL, NULL);

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
