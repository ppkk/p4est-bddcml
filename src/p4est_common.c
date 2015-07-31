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

#include "p4est_common.h"

/** Compute values at hanging nodes by interpolation.
 * A face hanging node in 3D depends on the four corner nodes of the face,
 * edge hanging nodes or face hanging nodes in 2D depend on two nodes.
 * This function works in place, we have to be careful about the ordering.
 * Face hanging node values are not reused, so they are overwritten first.
 * \param [in] face_code    This number encodes the child id of the quadrant
 *                          and the hanging status of faces and edges.
 * \param [in,out] inplace  On input, the values at the independent nodes.
 *                          On output, interpolated to hanging node locations.
 */
void
interpolate_hanging_nodes (p4est_lnodes_code_t face_code,
                           double inplace[P4EST_CHILDREN])
{
   const int           c = (int) (face_code & ones);
   int                 i, j;
   int                 ef;
   int                 work = (int) (face_code >> P4EST_DIM);
   double              sum;
   const double        factor = 1. / P4EST_HALF;

   /* Compute face hanging nodes first (this is all there is in 2D). */
   for (i = 0; i < P4EST_DIM; ++i) {
      if (work & 1) {
         ef = p4est_corner_faces[c][i];
         sum = 0.;
         for (j = 0; j < P4EST_HALF; ++j) {
            sum += inplace[p4est_face_corners[ef][j]];
         }
         inplace[c ^ ones ^ (1 << i)] = factor * sum;
      }
      work >>= 1;
   }

#ifdef P4_TO_P8
   /* Compute edge hanging nodes afterwards */
   for (i = 0; i < P4EST_DIM; ++i) {
      if (work & 1) {
         ef = p8est_corner_edges[c][i];
         inplace[c ^ (1 << i)] = .5 * (inplace[p8est_edge_corners[ef][0]] +
               inplace[p8est_edge_corners[ef][1]]);
      }
      work >>= 1;
   }
#endif
}
