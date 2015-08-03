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
void interpolate_hanging_nodes (p4est_lnodes_code_t face_code,
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

/** Decode the information from p{4,8}est_lnodes_t for a given element.
 *
 * \see p4est_lnodes.h for an in-depth discussion of the encoding.
 * \param [in] face_code         Bit code as defined in p{4,8}est_lnodes.h.
 * \param [out] hanging_corner   Undefined if no node is hanging.
 *                               If any node is hanging, this contains
 *                               one integer per corner, which is -1
 *                               for corners that are not hanging,
 *                               and the number of the non-hanging
 *                               corner on the hanging face/edge otherwise.
 *                               For faces in 3D, it is diagonally opposite.
 * \return true if any node is hanging, false otherwise.
 */
int lnodes_decode2 (p4est_lnodes_code_t face_code,
                int hanging_corner[P4EST_CHILDREN])
{
   if (face_code) {
      const int           c = (int) (face_code & ones);
      int                 i, h;
      int                 work = (int) (face_code >> P4EST_DIM);

      /* These two corners are never hanging by construction. */
      hanging_corner[c] = hanging_corner[c ^ ones] = -1;
      for (i = 0; i < P4EST_DIM; ++i) {
         /* Process face hanging corners. */
         h = c ^ (1 << i);
         hanging_corner[h ^ ones] = (work & 1) ? c : -1;
#ifdef P4_TO_P8
         /* Process edge hanging corners. */
         hanging_corner[h] = (work & P4EST_CHILDREN) ? c : -1;
#endif
         work >>= 1;
      }
      return 1;
   }
   return 0;
}


void plot_solution(p4est_t * p4est, p4est_lnodes_t * lnodes, double* u_sol, double* u_exact)
{
   p4est_topidx_t      tt;       /* Connectivity variables have this type. */
   p4est_locidx_t      k, q, Q, node_total;  /* Process-local counters have this type. */
   p4est_locidx_t      lni;      /* Node index relative to this processor. */
   p4est_tree_t       *tree;     /* Pointer to one octree */
   p4est_quadrant_t   *quad, *parent, sp, node;
   sc_array_t         *tquadrants;       /* Quadrant array for one tree */
   int                 i;
   double              loc_vertex_values_sol[P4EST_CHILDREN], loc_vertex_values_exact[P4EST_CHILDREN];


   /* Write the forest to disk for visualization, one file per processor. */
   double *u_interp_sol = NULL;
   double *u_interp_exact = NULL;
   if(u_sol)
      u_interp_sol = P4EST_ALLOC (double, p4est->local_num_quadrants * P4EST_CHILDREN);
   if(u_exact)
      u_interp_exact = P4EST_ALLOC (double, p4est->local_num_quadrants * P4EST_CHILDREN);

   /* Loop over local quadrants to apply the element matrices. */
   for (tt = p4est->first_local_tree, k = 0, node_total = 0;
        tt <= p4est->last_local_tree; ++tt) {
      tree = p4est_tree_array_index (p4est->trees, tt);
      tquadrants = &tree->quadrants;
      Q = (p4est_locidx_t) tquadrants->elem_count;

      for (q = 0; q < Q; ++q, ++k) {
         quad = p4est_quadrant_array_index (tquadrants, q);

         for (i = 0; i < P4EST_CHILDREN; ++i) {
            lni = lnodes->element_nodes[P4EST_CHILDREN * k + i];
            if(u_sol)
               loc_vertex_values_sol[i] = u_sol[lni];
            if(u_exact)
               loc_vertex_values_exact[i] = u_exact[lni];
         }

         if(u_sol)
            interpolate_hanging_nodes (lnodes->face_code[k], loc_vertex_values_sol);
         if(u_exact)
            interpolate_hanging_nodes (lnodes->face_code[k], loc_vertex_values_exact);

         for (i = 0; i < P4EST_CHILDREN; ++i) {
            if(u_sol)
               u_interp_sol[node_total]   = loc_vertex_values_sol[i];
            if(u_exact)
               u_interp_exact[node_total] = loc_vertex_values_exact[i];
            ++node_total;
         }
      }
   }

   if(u_sol)
   {
      if(u_exact)
      {
         p4est_vtk_write_all (p4est, NULL, 0.99999, 1, 1, 1, 0, 2, 0, "output",
                              "solution", u_interp_sol, "exact", u_interp_exact);
      }
      else
      {
         p4est_vtk_write_all (p4est, NULL, 0.99999, 1, 1, 1, 0, 1, 0, "output",
                              "solution", u_interp_sol);
      }
   }
   else
   {
      p4est_vtk_write_all (p4est, NULL, 0.99999, 1, 1, 1, 0, 0, 0, "output");
   }

   if(u_sol)
      P4EST_FREE (u_interp_sol);
   if(u_exact)
      P4EST_FREE (u_interp_exact);
}


/** info.
 * \param [in] p4est          The forest is not changed.
 * \param [in] lnodes         The node numbering is not changed.
 */
void print_p4est_mesh (p4est_t * p4est, p4est_lnodes_t * lnodes, int which_rank)
{
   const int           nloc = lnodes->num_local_nodes;
   int                 anyhang, hanging_corner[P4EST_CHILDREN];

   double              vxyz[3];  /* We embed the 2D vertices into 3D space. */
   int8_t             *bc;
   p4est_topidx_t      tt;       /* Connectivity variables have this type. */
   p4est_locidx_t      node_idx;      /* Node index relative to this processor. */
   p4est_quadrant_t   *quad, *parent, sp, node;

   /* Loop over local quadrants to apply the element matrices. */
   print_rank = which_rank;
   PPP printf("\n*************** BEGIN P4EST MESH ************************\n");
   PPP printf("rank %d, nodes owned count = %d, nodes global offset = %d\n", print_rank, lnodes->owned_count, (int)lnodes->global_offset);
   PPP printf("local num quadrants %d\n", p4est->local_num_quadrants);
   PPP printf("global nodes: ");
   for(int i = 0; i < lnodes->num_local_nodes; i++)
   {
      PPP printf("%ld, ", node_loc_to_glob(lnodes, i));
   }
   PPP printf("\nglobal[local]\n");

  p4est_locidx_t quad_idx;
  for_all_quads(p4est, quad_idx, quad)
  {

      PPP printf("elem %d: ", (int)p4est->global_first_quadrant[p4est->mpirank] + quad_idx);
      for(int lnode = 0; lnode < P4EST_CHILDREN; lnode++)
      {
         node_idx = lnodes->element_nodes[P4EST_CHILDREN * quad_idx + lnode];
         PPP printf("%d[%d], ", (int)node_loc_to_glob(lnodes, node_idx), node_idx);
      }
      //PPP printf("%d, %d, %d, %d", glob_all_lni[0], glob_all_lni[1], glob_all_lni[2], glob_all_lni[3]);


      /* Figure out the hanging corners on this element, if any. */
      anyhang = lnodes_decode2 (lnodes->face_code[quad_idx], hanging_corner);


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
         node_idx = lnodes->element_nodes[P4EST_CHILDREN * quad_idx + lnode];
         P4EST_ASSERT (node_idx >= 0 && node_idx < nloc);
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
         PPP printf("(%3.2lf, %3.2lf), ", vxyz[0], vxyz[1]);

      }
      PPP printf("\n");
   }}} // for all quads
   PPP printf("*************** END P4EST MESH ************************\n\n");


}



p4est_gloidx_t node_loc_to_glob(p4est_lnodes_t * lnodes, p4est_locidx_t loc_idx)
{
   if(loc_idx < lnodes->owned_count)
      return loc_idx + lnodes->global_offset;
   else
      return lnodes->nonlocal_nodes[loc_idx - lnodes->owned_count];
}

