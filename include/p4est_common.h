#ifndef P4EST_COMMON_H
#define P4EST_COMMON_H

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

#define for_all_quads(p4est, quad_idx, quad) \
   quad_idx = 0;\
   for (p4est_topidx_t tt = p4est->first_local_tree;\
        tt <= p4est->last_local_tree; ++tt) {\
      p4est_tree_t *tree = p4est_tree_array_index (p4est->trees, tt);\
      sc_array_t *tquadrants = &tree->quadrants;\
      p4est_locidx_t num_tree_quads = (p4est_locidx_t) tquadrants->elem_count;\
      for (p4est_locidx_t tree_quad_idx = 0; tree_quad_idx < num_tree_quads; ++tree_quad_idx, ++quad_idx) {\
         quad = p4est_quadrant_array_index (tquadrants, tree_quad_idx);\

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


void interpolate_hanging_nodes (p4est_lnodes_code_t face_code,
                           double inplace[P4EST_CHILDREN]);

int lnodes_decode2 (p4est_lnodes_code_t face_code,
                int hanging_corner[P4EST_CHILDREN]);

void plot_solution(p4est_t * p4est, p4est_lnodes_t * lnodes, double* u_sol, double* u_exact);

void print_p4est_mesh (p4est_t * p4est, p4est_lnodes_t * lnodes, int which_rank);

p4est_gloidx_t node_loc_to_glob(p4est_lnodes_t * lnodes, p4est_locidx_t loc_idx);

#endif // P4EST_COMMON_H
