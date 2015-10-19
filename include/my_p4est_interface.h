#ifndef P4EST_COMMON_H
#define P4EST_COMMON_H

#include <stdbool.h>
#include <vector>

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

#include "definitions.h"

#define for_all_quads(p4est, quad_idx, quad) \
   quad_idx = 0;\
   for (p4est_topidx_t tt = p4est->first_local_tree;\
        tt <= p4est->last_local_tree; ++tt) {\
      p4est_tree_t *tree = p4est_tree_array_index (p4est->trees, tt);\
      sc_array_t *tquadrants = &tree->quadrants;\
      p4est_locidx_t num_tree_quads = (p4est_locidx_t) tquadrants->elem_count;\
      for (p4est_locidx_t tree_quad_idx = 0; tree_quad_idx < num_tree_quads; ++tree_quad_idx, ++quad_idx) {\
         quad = p4est_quadrant_array_index (tquadrants, tree_quad_idx);\

#define end_for_all_quads }}

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
extern const int   *corner_to_hanging[P4EST_CHILDREN];


void interpolate_hanging_nodes (p4est_lnodes_code_t face_code,
                           double inplace[P4EST_CHILDREN]);

int lnodes_decode2 (p4est_lnodes_code_t face_code,
                int hanging_corner[P4EST_CHILDREN]);

void plot_solution(p4est_t * p4est, p4est_lnodes_t * lnodes, int num_components, double* u_sol, double* u_exact, int *partition);

void print_p4est_mesh (p4est_t * p4est, p4est_lnodes_t * lnodes, int which_rank);

p4est_gloidx_t node_loc_to_glob(p4est_lnodes_t * lnodes, p4est_locidx_t loc_idx);

/* 1D mass matrix on the reference element [0, 1]. */
static const double m_1d[2][2] = {
   {1 / 3., 1 / 6.},
   {1 / 6., 1 / 3.},
};
/* 1D stiffness matrix on the reference element [0, 1]. */
static const double s_1d[2][2] = {
   {1., -1.},
   {-1., 1.},
};

void generate_reference_matrices(real stiffness_dd[P4EST_CHILDREN][P4EST_CHILDREN], real mass_dd[P4EST_CHILDREN][P4EST_CHILDREN]);
void scale_reference_matrix(real ref_mat[P4EST_CHILDREN][P4EST_CHILDREN], real coefficient, real phys_elem_mat[P4EST_CHILDREN][P4EST_CHILDREN]);
//void prepare_transformed_values(int dimmension, std::vector<double> element_lengths);

void init_corner_to_hanging();

int independent_nodes(p4est_lnodes_t *lnodes, p4est_locidx_t quadrant, int lnode, p4est_locidx_t *nodes, real *coeffs);

int refine_uniform (p4est_t * p4est, p4est_topidx_t which_tree, p4est_quadrant_t * quadrant);

/** Callback function to decide on refinement.
 *
 * This function is called for every processor-local quadrant in order; its
 * return value is understood as a boolean refinement flag.  We refine around a
 * h = 1/8 block with left front lower corner (5/8, 2/8, 6/8).
 */
int refine_square (p4est_t * p4est, p4est_topidx_t which_tree, p4est_quadrant_t * quadrant);
int refine_center (p4est_t * p4est, p4est_topidx_t which_tree, p4est_quadrant_t * quadrant);
int refine_diagonal (p4est_t * p4est, p4est_topidx_t which_tree, p4est_quadrant_t * quadrant);
int refine_point (p4est_t * p4est, p4est_topidx_t which_tree, p4est_quadrant_t * quadrant);
int refine_circle (p4est_t * p4est, p4est_topidx_t which_tree, p4est_quadrant_t * quadrant);

void refine_and_partition(p4est_t* p4est, int num, p4est_refine_t fn);

class P4estClass
{

};

#endif // P4EST_COMMON_H
