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

#include <vector>

#include "arrays.h"
#include "p4est_common.h"

using namespace std;

int main (int argc, char **argv)
{
   int                 mpiret;
   sc_MPI_Comm         mpicomm;
   p4est_t            *p4est;
   p4est_connectivity_t *conn;

   mpiret = sc_MPI_Init (&argc, &argv);
   SC_CHECK_MPI (mpiret);
   mpicomm = sc_MPI_COMM_WORLD;
   mpiret = sc_MPI_Comm_rank(mpicomm, &mpi_rank);
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
   p4est = p4est_new (mpicomm, conn, 0, NULL, NULL);

   refine_and_partition(p4est, 2, refine_uniform);
   refine_and_partition(p4est, 2, refine_circle);
   refine_and_partition(p4est, 0, refine_square);
   refine_and_partition(p4est, 0, refine_point);
   refine_and_partition(p4est, 0, refine_diagonal);


   /* Destroy the p4est and the connectivity structure. */
   p4est_destroy (p4est);
   p4est_connectivity_destroy (conn);

   /* Verify that allocations internal to p4est and sc do not leak memory.
   * This should be called if sc_init () has been called earlier. */
   sc_finalize ();

   /* This is standard MPI programs.  Without --enable-mpi, this is a dummy. */
   mpiret = sc_MPI_Finalize ();
   SC_CHECK_MPI (mpiret);
   return 0;
}
