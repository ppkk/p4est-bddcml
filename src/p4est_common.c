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

#include "definitions.h"
#include "helpers.h"
#include "p4est_common.h"

const int   *corner_to_hanging[P4EST_CHILDREN];

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

void generate_reference_matrices(real stiffness_dd[P4EST_CHILDREN][P4EST_CHILDREN], real mass_dd[P4EST_CHILDREN][P4EST_CHILDREN])
{
   /* Compute entries of reference mass and stiffness matrices in 2D.
   * In this example we can proceed without numerical integration. */
   for (int l = 0; l < 2; ++l) {
      for (int k = 0; k < 2; ++k) {
         for (int j = 0; j < 2; ++j) {
            for (int i = 0; i < 2; ++i) {
#ifndef P4_TO_P8
               mass_dd[2 * j + i][2 * l + k] = m_1d[i][k] * m_1d[j][l];
               stiffness_dd[2 * j + i][2 * l + k] =
                     s_1d[i][k] * m_1d[j][l] + m_1d[i][k] * s_1d[j][l];
#else
               for (int n = 0; n < 2; ++n) {
                  for (int m = 0; m < 2; ++m) {
                     mass_dd[4 * i + 2 * n + m][4 * l + 2 * k + j] =
                           m_1d[m][j] * m_1d[n][k] * m_1d[i][l];
                     stiffness_dd[4 * i + 2 * n + m][4 * l + 2 * k + j] =
                           s_1d[m][j] * m_1d[n][k] * m_1d[i][l] +
                           m_1d[m][j] * s_1d[n][k] * m_1d[i][l] +
                           m_1d[m][j] * m_1d[n][k] * s_1d[i][l];
                  }
               }
#endif
            }
         }
      }
   }

//   for(int row = 0; row < P4EST_CHILDREN; row++)
//   {
//      for(int col = 0; col < P4EST_CHILDREN; col++)
//      {
//         PPP printf("%6.4lf, ", stiffness_dd[row][col]);
//      }
//      PPP printf("\n");
//   }
}


void init_corner_to_hanging()
{
   corner_to_hanging[0] = &zero;
#ifdef P4_TO_P8
   corner_to_hanging[1] = p8est_edge_corners[0];
   corner_to_hanging[2] = p8est_edge_corners[4];
   corner_to_hanging[3] = p8est_face_corners[4];
   corner_to_hanging[4] = p8est_edge_corners[8];
#endif
   corner_to_hanging[ones - 2] = p4est_face_corners[2];
   corner_to_hanging[ones - 1] = p4est_face_corners[0];
   corner_to_hanging[ones] = &ones;
}

static const real coeffs_regular[1] = {1};
static const real coeffs_hang_edge[2] = {0.5, 0.5};
static const real coeffs_hang_face[4] = {0.25, 0.25, 0.25, 0.25};

int independent_nodes(p4est_lnodes_t *lnodes, p4est_locidx_t quadrant, int lnode, p4est_locidx_t *nodes, real* coeffs)
{
   int anyhang, hanging_corner[P4EST_CHILDREN];

   /* Figure out the hanging corners on this element, if any. */
   anyhang = lnodes_decode2 (lnodes->face_code[quadrant], hanging_corner);

   if ((!anyhang) ||  (hanging_corner[lnode] == -1))
   {
      *coeffs = 1.;
      nodes[0] = lnodes->element_nodes[P4EST_CHILDREN * quadrant + lnode];
      return 1;
   }
   else
   {
      int c = hanging_corner[lnode];      /* Child id of quadrant. */
      int ncontrib = corner_num_hanging[lnode ^ c];
      const int *contrib_corner = corner_to_hanging[lnode ^ c];

      for (int j = 0; j < ncontrib; ++j) {
         int h = contrib_corner[j] ^ c;  /* Inverse transform of node number. */
         nodes[j] = lnodes->element_nodes[P4EST_CHILDREN * quadrant + h];
      }
      if(ncontrib == 2)
         *coeffs = 0.5;
      else if(ncontrib == 4)
         *coeffs = 0.25;
      else
         assert(0);

      return ncontrib;
   }
}

void get_quad_coords(p4est_quadrant_t * quadrant, double* coords, double *length)
{
   int quad_coords[3];
   quad_coords[0] = quadrant->x;
   quad_coords[1] = quadrant->y;
   int ndims = 2;
#ifdef P4_TO_P8
   quad_coords[2] = quadrant->z;
   ndims = 3;
#endif

   for(int dim = 0; dim < ndims; dim++)
   {
      int mask = 1 << P4EST_QMAXLEVEL;
      double add = 0.5;
      coords[dim] = 0.;
      while(mask > 0)
      {
         if(quad_coords[dim] & mask)
         {
            coords[dim] += add;
         }
         add /= 2.;
         mask >>= 1;
      }
   }

   *length = pow(0.5, quadrant->level);
}

int refine_uniform (p4est_t * p4est, p4est_topidx_t which_tree, p4est_quadrant_t * quadrant)
{
   return 1;
}

int refine_square (p4est_t * p4est, p4est_topidx_t which_tree, p4est_quadrant_t * quadrant)
{
   /* Compute the integer coordinate extent of a quadrant of length 2^(-3). */
   const p4est_qcoord_t one_over_64 = P4EST_QUADRANT_LEN (6);

   /* Compute the length of the current quadrant in integer coordinates. */
   const p4est_qcoord_t length = P4EST_QUADRANT_LEN (quadrant->level);

   /* Refine if the quadrant intersects the block in question. */
   return ((quadrant->x + length > 17 * one_over_64 && quadrant->x < 18 * one_over_64) &&
           (quadrant->y + length > 17 * one_over_64 && quadrant->y < 18 * one_over_64) &&
 #ifdef P4_TO_P8
           (quadrant->z + length > 17 * one_over_64 && quadrant->z < 18 * one_over_64) &&
 #endif
           1);
}

int refine_center (p4est_t * p4est, p4est_topidx_t which_tree, p4est_quadrant_t * quadrant)
{
//   printf("morton %d, %d, level %d; num %d\n", quadrant->x, quadrant->y, quadrant->level, 1<<29);

   return ((quadrant->x == 1 << P4EST_QMAXLEVEL) &&
           (quadrant->y == 1 << P4EST_QMAXLEVEL) &&
        #ifdef P4_TO_P8
           (quadrant->z == 0) &&
        #endif
           1);
}

int refine_points (p4est_t * p4est, p4est_topidx_t which_tree, p4est_quadrant_t * quadrant)
{
   //printf("morton %d, %d, level %d; num %d, %d\n", quadrant->x, quadrant->y, quadrant->level, 1<<29, (1<<27) + (1<<28));
   return ((quadrant->x == 1 << (P4EST_QMAXLEVEL-1)) &&
           (quadrant->y == 1 << (P4EST_QMAXLEVEL-1)) &&
          // (quadrant->y == (1 << P4EST_QMAXLEVEL) + (1 << (P4EST_QMAXLEVEL-1))) &&
        #ifdef P4_TO_P8
           (quadrant->z == 1 << (P4EST_QMAXLEVEL-1)) &&
        #endif
           1);
}

int refine_diagonal (p4est_t * p4est, p4est_topidx_t which_tree, p4est_quadrant_t * quadrant)
{
   return ((quadrant->x == quadrant->y) &&
        #ifdef P4_TO_P8
           (quadrant->x == quadrant->z) &&
        #endif
           1);
}

int refine_circle (p4est_t * p4est, p4est_topidx_t which_tree, p4est_quadrant_t * quadrant)
{
   double coords[3], work_coords[3], length;
   get_quad_coords(quadrant, coords, &length);

   int ndims = 2;
#ifdef P4_TO_P8
   ndims = 3;
#endif

   int num_inside = 0;
   int dim[3];
   work_coords[2] = 0.0;

   for(dim[0] = 0; dim[0] < 2; dim[0]++)
   {
      for(dim[1] = 0; dim[1] < 2; dim[1]++)
      {
#ifdef P4_TO_P8
         for(dim[2] = 0; dim[2] < 2; dim[2]++)
         {
#endif
            for(int i = 0; i < ndims; i++)
            {
               work_coords[i] = coords[i] + dim[i] * length;
            }

            double dist = sqrt(work_coords[0]*work_coords[0] + work_coords[1]*work_coords[1] + work_coords[2]*work_coords[2]);

            if(dist < 0.85)
               num_inside++;
#ifdef P4_TO_P8
         }
#endif
      }
   }

   if((num_inside == 0) || (num_inside == P4EST_CHILDREN))
      return 0;
   else
      return 1;

}

void refine_and_partition(p4est_t* p4est, int num, p4est_refine_t fn)
{
   for (int level = 0; level < num; ++level) {
      p4est_refine (p4est, 0, fn, NULL);
      p4est_partition (p4est, 0, NULL);
   }
}

