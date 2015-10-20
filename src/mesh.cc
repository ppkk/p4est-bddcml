#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <assert.h>

#include "arrays.h"
#include "bddcml_structs.h"
#include "my_p4est_interface.h"
#include "mesh.h"
#include "femspace.h"
#include "element.h"

using namespace std;

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



void BddcmlMesh::get_element(int elem_idx, Element *element)
{
   assert((elem_idx >= 0) && (elem_idx < subdomain_dims->n_elems));
   element->position.clear();

   // assuming cartesian grid...
   int first_node_idx = elem_node_indices.val[elem_idx * subdomain_dims->n_elem_nodes];
   for(int dim_idx = 0; dim_idx < subdomain_dims->n_problem_dims; dim_idx++)
   {
      element->position.push_back(coords.val[dim_idx][first_node_idx]);
   }
   element->size = element_lengths.val[elem_idx];

   // this cannot be used because of hanging nodes - we do not get the actual hanging node coordinates,
   // but the corresponding regular node coords.
   //element->size = coords.val[0][first_node_idx + 1] - coords.val[0][first_node_idx];

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


