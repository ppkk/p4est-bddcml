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

#include <assert.h>
#include <vector>
#include <mpi.h>

#include "my_p4est_interface.h"

sc_MPI_Comm         mpicomm;

#define TYPE_FIRST 0
#define TYPE_MIDDLE 1
#define TYPE_LAST 2

struct UserData
{
   int source_rank;
   int type;
};

void print_status(p4est_t *p4est)
{
   MPI_Barrier(mpicomm);
   printf("before, rank: %d, first: %ld, last: %ld, num: %d\n", mpi_rank, p4est->global_first_quadrant[mpi_rank],
          p4est->global_first_quadrant[mpi_rank]+p4est->local_num_quadrants - 1, p4est->local_num_quadrants);
   MPI_Barrier(mpicomm);
}


//
// equips each element with user data : source mpi_rank and type (FIRST/MIDDLE/LAST)
//
void fill_user_data(p4est_t *p4est)
{
   UserData* user_data_ptr;
   p4est_locidx_t quad_idx;
   p4est_quadrant_t *quad;
   bool is_first = true;

   for_all_quads(p4est, quad_idx, quad)
   {
      user_data_ptr = (UserData*)quad->p.user_data;
      user_data_ptr->source_rank = mpi_rank;
      if(is_first)
      {
         user_data_ptr->type = TYPE_FIRST;
         is_first = false;
      }
      else
      {
         user_data_ptr->type = TYPE_MIDDLE;
      }
   }
   end_for_all_quads
   user_data_ptr->type = TYPE_LAST;
}

#define MAX_RANKS_INVOLVED 20

//
// finds source processor for each element
//
void find_sources(p4est_t *p4est)
{
   UserData* user_data_ptr;
   p4est_locidx_t quad_idx;
   p4est_quadrant_t *quad;
   int target_displacement;
   MPI_Request req;

   int bounds[2];
   MPI_Win win;

   // bounds will contain mpi ranks of my first and last element after p4est migration
   MPI_Win_create(bounds, 2*sizeof(int), sizeof(int), MPI_INFO_NULL, mpicomm, &win);
   MPI_Win_fence(0, win);

   for_all_quads(p4est, quad_idx, quad)
   {
      user_data_ptr = (UserData*)quad->p.user_data;
      //printf("rank: %d, quad: \n", mpi_rank);
      if(user_data_ptr->type != TYPE_MIDDLE)
      {
         if(user_data_ptr->type == TYPE_FIRST)
            target_displacement = 0;
         if(user_data_ptr->type == TYPE_LAST)
            target_displacement = 1;

//         printf("rank: %d, writing to rank %d, type %d\n", mpi_rank, user_data_ptr->source_rank, target_displacement);
         //I know, that I have FIRST or LAST element of process user_data_ptr->source_rank
         //I can write it directly to its memory
         MPI_Put(&mpi_rank, 1, MPI_INT, user_data_ptr->source_rank, target_displacement, 1, MPI_INT, win);
      }
   }
   end_for_all_quads


   MPI_Win_fence(0, win);
   MPI_Win_free(&win);

   printf("rank %d: sending to %d ... %d\n", mpi_rank, bounds[0], bounds[1]);

   int gid_ranges_p4est_sent[MAX_RANKS_INVOLVED][2];
   int num_sent_to = 0;
   for(int rank_p4est_sent_to = bounds[0]; rank_p4est_sent_to <= bounds[1]; rank_p4est_sent_to++, num_sent_to++)
   {
      // Using the previous information, I know to which processes my elements went, so listen to them
      MPI_Irecv(gid_ranges_p4est_sent[num_sent_to], 2, MPI_INT, rank_p4est_sent_to, 0, mpicomm, &req);
   }
   assert(num_sent_to < MAX_RANKS_INVOLVED);


   int gid_ranges_p4est_obtained[MAX_RANKS_INVOLVED][2];
   int obtained_from_ranks[MAX_RANKS_INVOLVED];
   int num_obtained_from = 0;
   int global_offset = (int)p4est->global_first_quadrant[mpi_rank];
   int previous_rank = -1;

   // go through all elements and find ranges obtained from different ranks
   for_all_quads(p4est, quad_idx, quad)
   {
      user_data_ptr = (UserData*)quad->p.user_data;
      if(quad_idx == 0)
      {
         // first quadrant on this processor is allways a start of block from 1 sending processor
         gid_ranges_p4est_obtained[num_obtained_from][0] = global_offset + quad_idx;
         obtained_from_ranks[num_obtained_from] = user_data_ptr->source_rank;
      }
      else if(user_data_ptr->source_rank != previous_rank)
      {
         // this is a point, where block obtained from 1 processor ends and another begins
         gid_ranges_p4est_obtained[num_obtained_from][1] = global_offset + quad_idx - 1;
         num_obtained_from++;
         gid_ranges_p4est_obtained[num_obtained_from][0] = global_offset + quad_idx;
         obtained_from_ranks[num_obtained_from] = user_data_ptr->source_rank;
      }
      previous_rank = user_data_ptr->source_rank;
   }
   end_for_all_quads

   // close the last sequence
   gid_ranges_p4est_obtained[num_obtained_from][1] = global_offset + quad_idx - 1;
   num_obtained_from++;

   // send the colected data (ranges of elements originating on different processors)
   for(int i = 0; i < num_obtained_from; i++)
   {
      MPI_Isend(gid_ranges_p4est_obtained[i], 2, MPI_INT, obtained_from_ranks[i], 0, mpicomm, &req);
   }

   MPI_Barrier(mpicomm);

   // output the data
   char msg[1000];
   sprintf(msg, "rank %d: sent to ", mpi_rank);
   for(int i = 0; i < num_sent_to; i++)
      sprintf(msg, "%s p %d -> range (%d, %d), ", msg, bounds[0]+i, gid_ranges_p4est_sent[i][0], gid_ranges_p4est_sent[i][1]);
   printf("%s\n", msg);
   MPI_Barrier(mpicomm);

   sprintf(msg, "rank %d: obtained from ", mpi_rank);
   for(int i = 0; i < num_obtained_from; i++)
      sprintf(msg, "%s p %d -> range (%d, %d), ", msg, obtained_from_ranks[i], gid_ranges_p4est_obtained[i][0], gid_ranges_p4est_obtained[i][1]);
   printf("%s\n", msg);
   MPI_Barrier(mpicomm);

}



using namespace std;

int degree = 1;

int main (int argc, char **argv)
{
   int                 mpiret;
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
   p4est = p4est_new (mpicomm, conn, sizeof(UserData), NULL, NULL);

   // first few uniform refinements and partition
   p4est_refine (p4est, 0, refine_uniform, NULL);
   p4est_refine (p4est, 0, refine_uniform, NULL);
   p4est_partition (p4est, 0, NULL);

   // now the elements are distributed "uniformly" among processors

   // do some nonuniform refinements to cause disbalance
   p4est_refine (p4est, 0, refine_circle, NULL);
   p4est_refine (p4est, 0, refine_circle, NULL);
   // enforce 2:1 rule
   p4est_balance (p4est, P4EST_CONNECT_FULL, NULL);

   // before partitioning, fill "fake" element data (used for retrieveing migration information only
   fill_user_data(p4est);
   print_status(p4est);

   // now let p4est do the partitioning
   p4est_partition (p4est, 0, NULL);
   print_status(p4est);

   // reconstruct which elements were sent to which processor, using the element data which were provided
   find_sources(p4est);

   // finalize
   p4est_destroy (p4est);
   p4est_connectivity_destroy (conn);
   sc_finalize ();

   /* This is standard MPI programs.  Without --enable-mpi, this is a dummy. */
   mpiret = sc_MPI_Finalize ();
   SC_CHECK_MPI (mpiret);
   return 0;
}

