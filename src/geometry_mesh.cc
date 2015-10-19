#include "geometry_mesh.h"
#include "my_p4est_interface.h"


void GeometryMesh::prepare_subdomain_mesh(p4est_t *p4est, p4est_lnodes_t *lnodes)
{
   double              vxyz[3];  /* We embed the 2D vertices into 3D space. */
   p4est_quadrant_t   sp, node;

   p4est_locidx_t quad_idx = 0;
   p4est_quadrant_t *quad;

   for_all_quads(p4est, quad_idx, quad)
   {

      for (int lnode = 0; lnode < P4EST_CHILDREN; ++lnode) {
            p4est_quadrant_corner_node (quad, lnode, &node);

         /* Transform per-tree reference coordinates into physical space. */
         p4est_qcoord_to_vertex (p4est->connectivity, tt, node.x, node.y,
#ifdef P4_TO_P8
                                 node.z,
#endif
                                 vxyz);

         printf("(%lf, %lf), ", vxyz[0], vxyz[1]);
      }
      printf("\n");
   }
   end_for_all_quads
}

