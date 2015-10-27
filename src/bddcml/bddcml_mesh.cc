#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <assert.h>

#include "bddcml/bddcml_mesh.h"

using namespace std;

void BddcmlMesh::init(const BddcmlDimensions *subdomain_dims) {
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
   allocate_real_2D_array(subdomain_dims->n_nodes, subdomain_dims->n_problem_dims, &coords);
}

/**********************************************************************************************************/

void BddcmlMesh::free() {
   free_idx_array(&elem_node_indices);
   free_idx_array(&num_nodes_of_elem);
   free_idx_array(&elem_global_map);
   free_idx_array(&node_global_map);
   free_real_2D_array(&coords);
}

/**********************************************************************************************************/

void BddcmlMesh::print(int which_rank) const {
   print_rank = which_rank;
   PPP printf("\n*************** BEGIN BDDCML MESH ************************\n");
   PPP printf("elems: %d, nodes: %d\n", subdomain_dims->n_elems, subdomain_dims->n_nodes);
   PPP printf("linet: %d, lnnet: %d\n", elem_node_indices.len, num_nodes_of_elem.len);
   for(int elem = 0; elem < elem_global_map.len; elem++) {
      PPP printf("elem %d -> ", elem_global_map.val[elem]);
      for (int lnode = 0; lnode < num_nodes_of_elem.val[elem]; lnode++) {
         int node_local_idx = elem_node_indices.val[elem*subdomain_dims->n_elem_nodes + lnode];
         int node_global_idx = node_global_map.val[node_local_idx];
         PPP printf("%d[%d], ", node_global_idx, node_local_idx);
      }
      for (int lnode = 0; lnode < num_nodes_of_elem.val[elem]; lnode++) {
         int node_local_idx = elem_node_indices.val[elem*subdomain_dims->n_elem_nodes + lnode];
         PPP printf("(%3.2lf, %3.2lf), ", coords.val[0][node_local_idx], coords.val[1][node_local_idx]);
      }
      PPP printf("\n");
   }
   PPP printf("*************** END BDDCML MESH ************************\n\n");
}

/**********************************************************************************************************/

void BddcmlMesh::fill_nodes_info(const NodalElementMesh &nodal_mesh) {
   
}

/**********************************************************************************************************/
/**********************************************************************************************************/

BddcmlDimensions::BddcmlDimensions(int mesh_dim, PhysicsType physicsType) {
   n_problem_dims = mesh_dim;
   n_mesh_dims = mesh_dim;
   n_dofs = 0;
   n_elems = 0;
   n_nodes = 0;

   if(mesh_dim == 2)
      n_elem_nodes = 4;
   else if(mesh_dim == 3)
      n_elem_nodes = 8;
   else
      assert(0);

   if(physicsType == PhysicsType::LAPLACE)
      n_node_dofs = 1;
   else if(physicsType == PhysicsType::ELASTICITY)
      n_node_dofs = mesh_dim;
   else
      assert(0);
}



