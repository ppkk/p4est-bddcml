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

#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <assert.h>

#include "arrays.h"
#include "bddcml_structs.h"
#include "my_p4est_interface.h"
#include "mesh.h"
#include "femspace.h"
#include "assemble.h"

void BddcmlMesh::init(BddcmlDimensions* subdomain_dims)
{
   int nnodes_per_elem;
   if(subdomain_dims->n_mesh_dims == 2)
      nnodes_per_elem = 4;
   else if(subdomain_dims->n_mesh_dims == 3)
      nnodes_per_elem = 8;
   else
      assert(0);

   this->subdomain_dims = subdomain_dims;
   allocate_idx_array(subdomain_dims->n_elems * nnodes_per_elem, &elem_node_indices);
   allocate_idx_array(subdomain_dims->n_elems, &num_nodes_of_elem);
   allocate_idx_array(subdomain_dims->n_elems, &elem_global_map);
   allocate_idx_array(subdomain_dims->n_nodes, &node_global_map);
   allocate_real_array(subdomain_dims->n_elems, &element_lengths);
   allocate_real_2D_array(subdomain_dims->n_nodes, subdomain_dims->n_problem_dims, &coords);
}

void BddcmlMesh::free()
{
   free_idx_array(&elem_node_indices);
   free_idx_array(&num_nodes_of_elem);
   free_idx_array(&elem_global_map);
   free_idx_array(&node_global_map);
   free_real_array(&element_lengths);
   free_real_2D_array(&coords);
}

void BddcmlMesh::print(int which_rank)
{
   print_rank = which_rank;
   PPP printf("\n*************** BEGIN BDDCML MESH ************************\n");
   PPP printf("elems: %d, nodes: %d\n", subdomain_dims->n_elems, subdomain_dims->n_nodes);
   PPP printf("linet: %d, lnnet: %d\n", elem_node_indices.len, num_nodes_of_elem.len);
   for(int elem = 0; elem < elem_global_map.len; elem++)
   {
      PPP printf("elem %d -> ", elem_global_map.val[elem]);
      for (int lnode = 0; lnode < num_nodes_of_elem.val[elem]; lnode++)
      {
         int node_local_idx = elem_node_indices.val[elem*subdomain_dims->n_elem_nodes + lnode];
         int node_global_idx = node_global_map.val[node_local_idx];
         PPP printf("%d[%d], ", node_global_idx, node_local_idx);
      }
      for (int lnode = 0; lnode < num_nodes_of_elem.val[elem]; lnode++)
      {
         int node_local_idx = elem_node_indices.val[elem*subdomain_dims->n_elem_nodes + lnode];
         PPP printf("(%3.2lf, %3.2lf), ", coords.val[0][node_local_idx], coords.val[1][node_local_idx]);
      }
      PPP printf("\n");
   }
   PPP printf("*************** END BDDCML MESH ************************\n\n");
}


void BddcmlMesh::prepare_subdomain_mesh(p4est_t *p4est, p4est_lnodes_t *lnodes)
{
   double              vxyz[3];  /* We embed the 2D vertices into 3D space. */
   p4est_quadrant_t   sp, node;

   for(int lnode = 0; lnode < lnodes->num_local_nodes; lnode++)
   {
      node_global_map.val[lnode] = node_loc_to_glob(lnodes, lnode);
   }
   //p4est->global_first_quadrant
   /* Loop over local quadrants to apply the element matrices. */
   p4est_locidx_t quad_idx = 0;
   p4est_quadrant_t *quad;

   for_all_quads(p4est, quad_idx, quad)
   {
      // element local to global mapping -- obtained by adding the offset from p4est
      elem_global_map.val[quad_idx] = (int)p4est->global_first_quadrant[p4est->mpirank] + quad_idx;
      num_nodes_of_elem.val[quad_idx] = P4EST_CHILDREN;

      for (int lnode = 0; lnode < P4EST_CHILDREN; ++lnode) {
         /* Cache some information on corner nodes. */
         p4est_locidx_t node_idx = lnodes->element_nodes[P4EST_CHILDREN * quad_idx + lnode];
         //      isboundary[i] = (bc == NULL ? 0 : bc[lni]);
         //       inloc[i] = !isboundary[i] ? in[lni] : 0.;
         elem_node_indices.val[P4EST_CHILDREN * quad_idx + lnode] = node_idx;
      }

      element_lengths.val[quad_idx] = pow(0.5, quad->level);

      /* Figure out the hanging corners on this element, if any. */
      int hanging_corner[P4EST_CHILDREN];
      int anyhang = lnodes_decode2 (lnodes->face_code[quad_idx], hanging_corner);

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
         p4est_locidx_t node_idx = lnodes->element_nodes[P4EST_CHILDREN * quad_idx + lnode];
         P4EST_ASSERT (node_idx >= 0 && node_idx < subdomain_dims->n_nodes);
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

         coords.val[0][node_idx] = vxyz[0];
         coords.val[1][node_idx] = vxyz[1];
#ifdef P4_TO_P8
         coords.val[2][node_idx] = vxyz[2];
#endif



//         PPP printf("(%3.2lf, %3.2lf), ", vxyz[0], vxyz[1]);

      }
//      PPP printf("\n");

   }
   end_for_all_quads

}

void BddcmlMesh::get_element(int elem_idx, Element *element)
{
   assert((elem_idx >= 0) && (elem_idx < subdomain_dims->n_elems));
   element->position.clear();
   for(int dim_idx = 0; dim_idx < subdomain_dims->n_problem_dims; dim_idx++)
      element->position.push_back(coords.val[dim_idx][elem_idx]);
   element->size = element_lengths.val[elem_idx];
}


void init_dimmensions(BddcmlDimensions* dimmensions, int mesh_dim, PhysicsType physicsType)
{
   dimmensions->n_problem_dims = mesh_dim;
   dimmensions->n_mesh_dims = mesh_dim;
   dimmensions->n_dofs = 0;
   dimmensions->n_elems = 0;
   dimmensions->n_nodes = 0;

   if(mesh_dim == 2)
      dimmensions->n_elem_nodes = 4;
   else if(mesh_dim == 3)
      dimmensions->n_elem_nodes = 8;
   else
      assert(0);

   if(physicsType == PhysicsType::LAPLACE)
      dimmensions->n_node_dofs = 1;
   else if(physicsType == PhysicsType::ELASTICITY)
      dimmensions->n_node_dofs = mesh_dim;
   else
      assert(0);
}


void prepare_dimmensions(p4est_t *p4est, p4est_lnodes_t *lnodes, PhysicsType physicsType,
                         BddcmlDimensions *subdomain_dims, BddcmlDimensions *global_dims,
                         sc_MPI_Comm mpicomm)
{
#ifndef P4_TO_P8
   init_dimmensions(subdomain_dims, 2, physicsType);
   init_dimmensions(global_dims, 2, physicsType);
#else
   init_dimmensions(subdomain_dims, 3, physicsType);
   init_dimmensions(global_dims, 3, physicsType);
#endif
   subdomain_dims->n_nodes = lnodes->num_local_nodes;
   subdomain_dims->n_dofs  = lnodes->num_local_nodes * subdomain_dims->n_node_dofs;
   subdomain_dims->n_elems = lnodes->num_local_elements;
   printf("proc %d, elems %d, nodes %d\n", mpi_rank, subdomain_dims->n_elems, subdomain_dims->n_nodes);

   int global_num_nodes;
   if(mpi_rank == mpi_size - 1)
   {
      global_num_nodes = lnodes->global_offset + lnodes->owned_count;
   }
   sc_MPI_Bcast(&global_num_nodes, 1, MPI_INT, mpi_size - 1, mpicomm);

   global_dims->n_nodes = global_num_nodes;
   global_dims->n_dofs = global_num_nodes * global_dims->n_node_dofs;
   global_dims->n_elems = p4est->global_num_quadrants;
}

