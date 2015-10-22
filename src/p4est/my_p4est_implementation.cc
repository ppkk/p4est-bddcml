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

#include <vector>

#include "p4est/my_p4est_implementation.h"
#include "bddcml/bddcml_mesh.h"
#include "geometry_mesh.h"

using namespace std;

// this file defines both P4estClass2D and P4estClass3D using p4est macros
#ifndef P4_TO_P8
   #define P4estClassDim P4estClass2D
   #define P4estNamespaceDim P4estNamespace2D
#else
   #define P4estClassDim P4estClass3D
   #define P4estNamespaceDim P4estNamespace3D
#endif

#define for_all_quads(p4est, quad_idx, quad) \
   quad_idx = 0;\
   for (p4est_topidx_t tt = p4est->first_local_tree;\
        tt <= p4est->last_local_tree; ++tt) {\
      p4est_tree_t *tree_ptr = p4est_tree_array_index (p4est->trees, tt);\
      sc_array_t *tquadrants = &tree_ptr->quadrants;\
      p4est_locidx_t num_tree_quads = (p4est_locidx_t) tquadrants->elem_count;\
      for (p4est_locidx_t tree_quad_idx = 0; tree_quad_idx < num_tree_quads; ++tree_quad_idx, ++quad_idx) {\
         quad = p4est_quadrant_array_index (tquadrants, tree_quad_idx);\
         assert(quad);

#define for_all_quads_with_tree(p4est, tt, quad_idx, quad) \
   quad_idx = 0;\
   for (tt = p4est->first_local_tree;\
        tt <= p4est->last_local_tree; ++tt) {\
      p4est_tree_t *tree_ptr = p4est_tree_array_index (p4est->trees, tt);\
      sc_array_t *tquadrants = &tree_ptr->quadrants;\
      p4est_locidx_t num_tree_quads = (p4est_locidx_t) tquadrants->elem_count;\
      for (p4est_locidx_t tree_quad_idx = 0; tree_quad_idx < num_tree_quads; ++tree_quad_idx, ++quad_idx) {\
         quad = p4est_quadrant_array_index (tquadrants, tree_quad_idx);\
         assert(quad);

#define end_for_all_quads }}

namespace P4estNamespaceDim
{

const int    zero = 0;           /**< Constant zero. */
const int    ones = P4EST_CHILDREN - 1;  /**< One bit per dimension. */

/** For each node i of the reference quadrant, corner_num_hanging[i] many. */
const int   *corner_to_hanging[P4EST_CHILDREN];


void init_corner_to_hanging() {
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
                           double inplace[P4EST_CHILDREN]) {
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
                int hanging_corner[P4EST_CHILDREN]) {
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


/** List number of possible independent nodes for each hanging node. */
static const int    corner_num_hanging[P4EST_CHILDREN] =
#ifndef P4_TO_P8
   { 1, 2, 2, 1}
#else
   { 1, 2, 2, 4, 2, 4, 4, 1 }
#endif
;

} //P4estNamespaceDim

//****************************************************************************************

P4estClassDim::P4estClassDim(int degree, sc_MPI_Comm mpicomm) : P4estClass(degree, mpicomm) {
   P4estNamespaceDim::init_corner_to_hanging();

#ifndef P4_TO_P8
   num_dim = 2;
   children = 4;
   conn = p4est_connectivity_new_unitsquare ();
#else
   num_dim = 3;
   children = 8;
   conn = p8est_connectivity_new_unitcube ();
#endif

   p4est = p4est_new (mpicomm, conn, 0, NULL, NULL);
   __lnodes_possibly_not_consistent = nullptr;
}

//****************************************************************************************

P4estClassDim::~P4estClassDim() {
   /* Destroy the p4est and the connectivity structure. */
   destroy_lnodes();
   p4est_destroy(p4est);
   p4est_connectivity_destroy(conn);
}

void P4estClassDim::init() {

}

//****************************************************************************************

void P4estClassDim::update_lnodes() const {
   if(__lnodes_possibly_not_consistent != nullptr)
      return;

   /* Create the ghost layer to learn about parallel neighbors. */
   p4est_ghost_t *ghost = p4est_ghost_new (p4est, P4EST_CONNECT_FULL);

   /* Create a node numbering for continuous linear finite elements. */
   __lnodes_possibly_not_consistent = p4est_lnodes_new (p4est, ghost, degree);

   /* Destroy the ghost structure -- no longer needed after node creation. */
   p4est_ghost_destroy (ghost);

}

//****************************************************************************************

void P4estClassDim::destroy_lnodes() {
   if(__lnodes_possibly_not_consistent == nullptr)
      return;

   p4est_lnodes_destroy (__lnodes_possibly_not_consistent);
   __lnodes_possibly_not_consistent = nullptr;
}

//****************************************************************************************

p4est_gloidx_t P4estClassDim::node_loc_to_glob(p4est_locidx_t loc_idx) const {
   if(loc_idx < lnodes()->owned_count)
      return loc_idx + lnodes()->global_offset;
   else
      return lnodes()->nonlocal_nodes[loc_idx - lnodes()->owned_count];
}

//****************************************************************************************

void P4estClassDim::plot_solution(int num_comp, double* u_sol, double* u_exact, int* partition) const {
   p4est_locidx_t      quad_idx, node_total = 0;  /* Process-local counters have this type. */
   double              loc_vertex_values_sol[num_comp][P4EST_CHILDREN], loc_vertex_values_exact[P4EST_CHILDREN];
   p4est_quadrant_t *quad;

   /* Write the forest to disk for visualization, one file per processor. */
   double *u_interp_sol[num_comp];
   double *u_interp_exact = NULL;
   double *interp_partition = NULL;
   if(u_sol)
      for(int comp = 0; comp < num_comp; comp++)
         u_interp_sol[comp] = P4EST_ALLOC (double, p4est->local_num_quadrants * P4EST_CHILDREN);
   if(u_exact)
      u_interp_exact = P4EST_ALLOC (double, p4est->local_num_quadrants * P4EST_CHILDREN);
   if(partition)
      interp_partition = P4EST_ALLOC (double, p4est->local_num_quadrants * P4EST_CHILDREN);

   for_all_quads(p4est, quad_idx, quad){
      for (int i = 0; i < P4EST_CHILDREN; ++i) {
         int node = lnodes()->element_nodes[P4EST_CHILDREN * quad_idx + i];
         if(u_sol) {
            for(int comp = 0; comp < num_comp; comp++) {
               int dof = num_comp * node + comp;
               loc_vertex_values_sol[comp][i] = u_sol[dof];
            }
         }
         if(u_exact)
            loc_vertex_values_exact[i] = u_exact[node];
      }

      if(u_sol) {
         for(int comp = 0; comp < num_comp; comp++) {
            P4estNamespaceDim::interpolate_hanging_nodes (lnodes()->face_code[quad_idx], loc_vertex_values_sol[comp]);
         }
      }
      if(u_exact)
         P4estNamespaceDim::interpolate_hanging_nodes (lnodes()->face_code[quad_idx], loc_vertex_values_exact);

      for (int i = 0; i < P4EST_CHILDREN; ++i) {
         if(u_sol) {
            for(int comp = 0; comp < num_comp; comp++) {
               u_interp_sol[comp][node_total]   = loc_vertex_values_sol[comp][i];
            }
         }
         if(u_exact)
            u_interp_exact[node_total] = loc_vertex_values_exact[i];
         if(partition)
            interp_partition[node_total] = partition[quad_idx];
         ++node_total;
      }
   }
   end_for_all_quads

   if(u_sol) {
      if(u_exact) {
         p4est_vtk_write_all (p4est, NULL, 0.99999, 1, 1, 1, 0, 2, 0, "output",
                              "solution", u_interp_sol, "exact", u_interp_exact);
      }
      else {
         if(num_comp == 1)
            p4est_vtk_write_all (p4est, NULL, 0.99999, 1, 1, 1, 0, 1, 0, "output",
                                 "solution", u_interp_sol[0]);
         else if(num_comp == 2)
            p4est_vtk_write_all (p4est, NULL, 0.99999, 1, 1, 1, 0, 2, 0, "output",
                                 "solution_1", u_interp_sol[0], "solution_2", u_interp_sol[1]);
         else if(num_comp == 3)
            p4est_vtk_write_all (p4est, NULL, 0.99999, 1, 1, 1, 0, 3, 0, "output",
                                 "solution_1", u_interp_sol[0], "solution_2", u_interp_sol[1], "solution_3", u_interp_sol[2]);
         else
            assert(0);
      }
   }
   else {
      if(partition) {
         p4est_vtk_write_all (p4est, NULL, 0.99999, 1, 1, 1, 0, 1, 0, "output",
                              "metis_part", interp_partition);
      }
      else {
         p4est_vtk_write_all (p4est, NULL, 0.99999, 1, 1, 1, 0, 0, 0, "output");
      }
   }

   if(u_sol)
      for(int comp = 0; comp < num_comp; comp++)
         P4EST_FREE (u_interp_sol[comp]);
   if(u_exact)
      P4EST_FREE (u_interp_exact);
   if(partition)
      P4EST_FREE (interp_partition);
}

//****************************************************************************************

/** info.
 * \param [in] p4est          The forest is not changed.
 * \param [in] lnodes         The node numbering is not changed.
 */
void P4estClassDim::print_p4est_mesh(int which_rank) const {
   const int           nloc = lnodes()->num_local_nodes;
   int                 anyhang, hanging_corner[P4EST_CHILDREN];

   double              vxyz[3];  /* We embed the 2D vertices into 3D space. */
   p4est_locidx_t      node_idx;      /* Node index relative to this processor. */
   p4est_quadrant_t   *quad, *parent, sp, node;

   /* Loop over local quadrants to apply the element matrices. */
   print_rank = which_rank;
   PPP printf("\n*************** BEGIN P4EST MESH ************************\n");
   PPP printf("rank %d, nodes owned count = %d, nodes global offset = %d\n", print_rank, lnodes()->owned_count, (int)lnodes()->global_offset);
   PPP printf("local num quadrants %d\n", p4est->local_num_quadrants);
   PPP printf("global nodes: ");
   for(int i = 0; i < lnodes()->num_local_nodes; i++) {
      PPP printf("%ld, ", node_loc_to_glob(i));
   }
   PPP printf("\nglobal[local]\n");

  p4est_locidx_t quad_idx;
  for_all_quads(p4est, quad_idx, quad) {
      PPP printf("elem %d: ", (int)p4est->global_first_quadrant[p4est->mpirank] + quad_idx);
      for(int lnode = 0; lnode < P4EST_CHILDREN; lnode++) {
         node_idx = lnodes()->element_nodes[P4EST_CHILDREN * quad_idx + lnode];
         PPP printf("%d[%d], ", (int)node_loc_to_glob(node_idx), node_idx);
      }
      //PPP printf("%d, %d, %d, %d", glob_all_lni[0], glob_all_lni[1], glob_all_lni[2], glob_all_lni[3]);


      /* Figure out the hanging corners on this element, if any. */
      anyhang = P4estNamespaceDim::lnodes_decode2 (lnodes()->face_code[quad_idx], hanging_corner);


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
         node_idx = lnodes()->element_nodes[P4EST_CHILDREN * quad_idx + lnode];
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
   }
  end_for_all_quads

   PPP printf("*************** END P4EST MESH ************************\n\n");


}

//****************************************************************************************

int P4estClassDim::independent_nodes(p4est_locidx_t quadrant, int lnode, p4est_locidx_t *nodes, real* coeffs) const {
   int anyhang, hanging_corner[P4EST_CHILDREN];

   /* Figure out the hanging corners on this element, if any. */
   anyhang = P4estNamespaceDim::lnodes_decode2 (lnodes()->face_code[quadrant], hanging_corner);

   if ((!anyhang) ||  (hanging_corner[lnode] == -1)) {
      *coeffs = 1.;
      nodes[0] = lnodes()->element_nodes[P4EST_CHILDREN * quadrant + lnode];
      return 1;
   }
   else {
      int c = hanging_corner[lnode];      /* Child id of quadrant. */
      int ncontrib = P4estNamespaceDim::corner_num_hanging[lnode ^ c];
      const int *contrib_corner = P4estNamespaceDim::corner_to_hanging[lnode ^ c];

      for (int j = 0; j < ncontrib; ++j) {
         int h = contrib_corner[j] ^ c;  /* Inverse transform of node number. */
         nodes[j] = lnodes()->element_nodes[P4EST_CHILDREN * quadrant + h];
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

//****************************************************************************************

void get_quad_coords(p4est_quadrant_t * quadrant, double* coords, double *length) {
   int quad_coords[3];
   quad_coords[0] = quadrant->x;
   quad_coords[1] = quadrant->y;
   int ndims = 2;
#ifdef P4_TO_P8
   quad_coords[2] = quadrant->z;
   ndims = 3;
#endif

   for(int dim = 0; dim < ndims; dim++) {
      int mask = 1 << P4EST_QMAXLEVEL;
      double add = 0.5;
      coords[dim] = 0.;
      while(mask > 0) {
         if(quad_coords[dim] & mask) {
            coords[dim] += add;
         }
         add /= 2.;
         mask >>= 1;
      }
   }

   *length = pow(0.5, quadrant->level);
}

int refine_uniform (p4est_t * /*p4est*/, p4est_topidx_t /*which_tree*/, p4est_quadrant_t * /*quadrant*/) {
   return 1;
}

int refine_square (p4est_t * p4est, p4est_topidx_t which_tree, p4est_quadrant_t * quadrant) {
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

int refine_center (p4est_t * p4est, p4est_topidx_t which_tree, p4est_quadrant_t * quadrant) {
//   printf("morton %d, %d, level %d; num %d\n", quadrant->x, quadrant->y, quadrant->level, 1<<29);

   return ((quadrant->x == 1 << P4EST_QMAXLEVEL) &&
           (quadrant->y == 1 << P4EST_QMAXLEVEL) &&
        #ifdef P4_TO_P8
           (quadrant->z == 0) &&
        #endif
           1);
}

int refine_point (p4est_t * p4est, p4est_topidx_t which_tree, p4est_quadrant_t * quadrant) {
   //printf("morton %d, %d, level %d; num %d, %d\n", quadrant->x, quadrant->y, quadrant->level, 1<<29, (1<<27) + (1<<28));
   return ((quadrant->x == 1 << (P4EST_QMAXLEVEL-1)) &&
           (quadrant->y == 1 << (P4EST_QMAXLEVEL-1)) &&
          // (quadrant->y == (1 << P4EST_QMAXLEVEL) + (1 << (P4EST_QMAXLEVEL-1))) &&
        #ifdef P4_TO_P8
           (quadrant->z == 1 << (P4EST_QMAXLEVEL-1)) &&
        #endif
           1);
}

int refine_diagonal (p4est_t * p4est, p4est_topidx_t which_tree, p4est_quadrant_t * quadrant) {
   return ((quadrant->x == quadrant->y) &&
        #ifdef P4_TO_P8
           (quadrant->x == quadrant->z) &&
        #endif
           1);
}

int refine_circle (p4est_t * p4est, p4est_topidx_t which_tree, p4est_quadrant_t * quadrant) {
   double coords[3], work_coords[3], length;
   get_quad_coords(quadrant, coords, &length);

   int ndims = 2;
#ifdef P4_TO_P8
   ndims = 3;
#endif

   int num_inside = 0;
   int dim[3];
   work_coords[2] = 0.0;

   for(dim[0] = 0; dim[0] < 2; dim[0]++) {
      for(dim[1] = 0; dim[1] < 2; dim[1]++) {
#ifdef P4_TO_P8
         for(dim[2] = 0; dim[2] < 2; dim[2]++) {
#endif
            for(int i = 0; i < ndims; i++) {
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

//****************************************************************************************

void P4estClassDim::refine_and_partition(int num, RefineType type) {
   if(num == 0)
      return;

   destroy_lnodes();

   const char *name;
   p4est_refine_t fn;

   switch(type) {
   case RefineType::UNIFORM :
      name = "UNIFORM";
      fn = refine_uniform;
      break;
   case RefineType::CIRCLE :
      name = "CIRCLE";
      fn = refine_circle;
      break;
   case RefineType::SQUARE :
      name = "SQUARE";
      fn = refine_square;
      break;
   case RefineType::POINT :
      name = "POINT";
      fn = refine_point;
      break;
   case RefineType::DIAGONAL :
      name = "DIAGONAL";
      fn = refine_diagonal;
      break;

   default:
      assert(0);
   }


   int nelems_before = p4est->global_num_quadrants;
   for (int level = 0; level < num; ++level) {
      p4est_refine (p4est, 0, fn, NULL);
      p4est_partition (p4est, 0, NULL);
   }

   p4est_balance (p4est, P4EST_CONNECT_FULL, NULL);
   p4est_partition (p4est, 0, NULL);


   if(mpi_rank == 0) {
      int added = p4est->global_num_quadrants - nelems_before;
      printf("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
      printf("Mesh refinement %s, %d levels, added %d elements (%3.2lf %%), now %ld elements\n",
             name, num, added, 100*(double)added/p4est->global_num_quadrants, p4est->global_num_quadrants);
      printf("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
  }
}

//****************************************************************************************

void P4estClassDim::prepare_dimmensions(BddcmlDimensions *subdomain_dims, BddcmlDimensions *global_dims) const {
#ifndef P4_TO_P8
   assert(num_dim == 2);
#else
   assert(num_dim == 3);
#endif
   subdomain_dims->n_nodes = lnodes()->num_local_nodes;
   subdomain_dims->n_dofs  = lnodes()->num_local_nodes * subdomain_dims->n_node_dofs;
   subdomain_dims->n_elems = lnodes()->num_local_elements;
   printf("proc %d, elems %d, nodes %d\n", mpi_rank, subdomain_dims->n_elems, subdomain_dims->n_nodes);

   int global_num_nodes;
   if(mpi_rank == mpi_size - 1) {
      global_num_nodes = lnodes()->global_offset + lnodes()->owned_count;
   }
   sc_MPI_Bcast(&global_num_nodes, 1, MPI_INT, mpi_size - 1, mpicomm);

   global_dims->n_nodes = global_num_nodes;
   global_dims->n_dofs = global_num_nodes * global_dims->n_node_dofs;
   global_dims->n_elems = p4est->global_num_quadrants;
}

//****************************************************************************************

void P4estClassDim::get_node_coords(p4est_topidx_t tree, const p4est_quadrant_t &node, vector<double> *coords) const {
   double vxyz[3];  /* We embed the 2D vertices into 3D space. */

   coords->clear();
   /* Transform per-tree reference coordinates into physical space. */
   p4est_qcoord_to_vertex (p4est->connectivity, tree, node.x, node.y,
#ifdef P4_TO_P8
                           node.z,
#endif
                           vxyz);
   for(int i = 0; i < num_dim; i++) {
      coords->push_back(vxyz[i]);
   }

}

//****************************************************************************************

void P4estClassDim::prepare_subdomain_bddcml_mesh(BddcmlMesh* mesh) const {
   vector<double> coords;
   p4est_quadrant_t   sp, node;

   for(int lnode = 0; lnode < lnodes()->num_local_nodes; lnode++) {
      mesh->node_global_map.val[lnode] = node_loc_to_glob(lnode);
   }
   //p4est->global_first_quadrant
   /* Loop over local quadrants to apply the element matrices. */
   p4est_locidx_t quad_idx = 0;
   p4est_quadrant_t *quad;
   p4est_topidx_t tree;

   for_all_quads_with_tree(p4est, tree, quad_idx, quad) {
      // element local to global mapping -- obtained by adding the offset from p4est
      mesh->elem_global_map.val[quad_idx] = (int)p4est->global_first_quadrant[p4est->mpirank] + quad_idx;
      mesh->num_nodes_of_elem.val[quad_idx] = P4EST_CHILDREN;

      for (int lnode = 0; lnode < P4EST_CHILDREN; ++lnode) {
         /* Cache some information on corner nodes. */
         p4est_locidx_t node_idx = lnodes()->element_nodes[P4EST_CHILDREN * quad_idx + lnode];
         //      isboundary[i] = (bc == NULL ? 0 : bc[lni]);
         //       inloc[i] = !isboundary[i] ? in[lni] : 0.;
         mesh->elem_node_indices.val[P4EST_CHILDREN * quad_idx + lnode] = node_idx;
      }

      mesh->element_lengths.val[quad_idx] = pow(0.5, quad->level);

      /* Figure out the hanging corners on this element, if any. */
      int hanging_corner[P4EST_CHILDREN];
      int anyhang = P4estNamespaceDim::lnodes_decode2 (lnodes()->face_code[quad_idx], hanging_corner);

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
         p4est_locidx_t node_idx = lnodes()->element_nodes[P4EST_CHILDREN * quad_idx + lnode];
         P4EST_ASSERT (node_idx >= 0 && node_idx < mesh->subdomain_dims->n_nodes);
         if (anyhang && hanging_corner[lnode] >= 0) {
            /* This node is hanging; access the referenced node instead. */
            p4est_quadrant_corner_node (parent, lnode, &node);
         }
         else {
            p4est_quadrant_corner_node (quad, lnode, &node);
         }

         get_node_coords(tree, node, &coords);
         for(int i = 0; i < num_dim; i++)
            mesh->coords.val[i][node_idx] = coords[i];
      }
   }
   end_for_all_quads
}

//****************************************************************************************

const int nodes_in_axes_dirs[3] = {1,2,4};
const double tolerance = 1e-10;

void P4estClassDim::prepare_subdomain_geometry_mesh(GeometryMesh *mesh) const {
   p4est_quadrant_t node;
   p4est_locidx_t quad_idx = 0;
   p4est_quadrant_t *quad;
   p4est_topidx_t tree;
   Element element;
   vector<double> other_point;

   mesh->clear();

   for_all_quads_with_tree(p4est, tree, quad_idx, quad) {
      element.clear();
      p4est_quadrant_corner_node (quad, 0, &node);
      get_node_coords(tree, node, &element.position);

      p4est_quadrant_corner_node (quad, nodes_in_axes_dirs[0], &node);
      get_node_coords(tree, node, &other_point);
      element.size = other_point[0] - element.position[0];

      for(int i = 1; i < num_dim; i++) {
         p4est_quadrant_corner_node (quad, nodes_in_axes_dirs[i], &node);
         get_node_coords(tree, node, &other_point);
         assert(element.size - (other_point[i] - element.position[i]) < tolerance * element.size);
      }

      mesh->elements.push_back(element);
   }
   end_for_all_quads
}

//****************************************************************************************