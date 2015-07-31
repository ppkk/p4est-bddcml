#ifndef P4EST_COMMON_H
#define P4EST_COMMON_H

#define PPP if(mpi_rank == print_rank)

/** List number of possible independent nodes for each hanging node. */
static const int    corner_num_hanging[P4EST_CHILDREN] =
#ifndef P4_TO_P8
   { 1, 2, 2, 1 }
#else
   { 1, 2, 2, 4, 2, 4, 4, 1 }
#endif
;

static const int    zero = 0;           /**< Constant zero. */
static const int    ones = P4EST_CHILDREN - 1;  /**< One bit per dimension. */

/** For each node i of the reference quadrant, corner_num_hanging[i] many. */
static const int   *corner_to_hanging[P4EST_CHILDREN];


void
interpolate_hanging_nodes (p4est_lnodes_code_t face_code,
                           double inplace[P4EST_CHILDREN]);

#endif // P4EST_COMMON_H
