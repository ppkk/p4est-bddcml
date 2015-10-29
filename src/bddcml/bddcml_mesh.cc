#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <assert.h>

#include "bddcml/bddcml_mesh.h"
#include "p4est/my_p4est_interface.h"
#include "local_matrix.h"
#include "shapefun.h"
#include "element.h"
#include "integration_cell.h"

using namespace std;

void BddcmlMesh::init(const BddcmlDimensions *subdomain_dims) {
   this->subdomain_dims = subdomain_dims;
   allocate_idx_array(subdomain_dims->n_elems * Def::d()->d()->num_element_nodes, &elem_node_indices);
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

// this is not nice, but is here from testing purposes, can be removed later
const double FORBIDDEN = -1.2345e40;
bool reals_equal(real x, real y) {
   if(fabs(x) + fabs(y) < 1e-10)
         return true;
   return fabs(x-y) / (fabs(x) + fabs(y)) < 1e-10;
}

void BddcmlMesh::fill_nodes_info(const P4estClass &p4est, const NodalElementMesh &nodal_mesh) {
   // First fill in "forbidden" value to coordinates
   for(int node = 0; node < subdomain_dims->n_nodes; node++) {
      for(int dim = 0; dim < subdomain_dims->n_mesh_dims; dim++) {
         coords.val[dim][node] = FORBIDDEN;
      }
   }

   ReferenceElement ref_elem(Def::d()->d()->num_dim, Def::d()->d()->order);
   HangingInfo hanging_info(p4est);
   vector<vector<real> > nodes_coords(Def::d()->num_element_nodes, vector<real>(Def::d()->num_dim, 0.0));
   vector<vector<real> > parent_nodes_coords(Def::d()->num_element_nodes, vector<real>(Def::d()->num_dim, 0.0));

   for(unsigned elem_idx = 0; elem_idx < nodal_mesh.elements.size(); elem_idx++) {
      // first nodes number and indices
      num_nodes_of_elem.val[elem_idx] = Def::d()->num_element_nodes;
      const NodalElement& nodal_element(nodal_mesh.elements[elem_idx]);
//      cout << "elem start " << nodal_element.cell.position[0] << ", " << nodal_element.cell.position[1] << endl;
      for (int lnode = 0; lnode < Def::d()->num_element_nodes; ++lnode) {
         int node_idx = nodal_element.nodes[lnode];
         elem_node_indices.val[Def::d()->num_element_nodes * elem_idx + lnode] = node_idx;
      }

      // now find information of nodes coordinates.
      p4est.get_hanging_info(elem_idx, &hanging_info);
      ref_elem.find_nodes_coords(nodal_element.cell, &nodes_coords);

      IntegrationCell parent_cell;
      nodal_element.cell.fill_parent_cell(&parent_cell);
      ref_elem.find_nodes_coords(parent_cell, &parent_nodes_coords);

//      cout << "nodes coords ";
//      for(auto x : nodes_coords) {
//         cout << "(";
//         for(double c : x) {
//            cout << c << ", ";
//         }
//         cout << "),";
//      }
//      cout << endl;


      for(int lnode = 0; lnode < Def::d()->num_element_nodes; lnode++)
      {
         int node_idx = nodal_element.nodes[lnode];
//         cout << "node " << node_idx << endl;
         // first check, if this node is on hanging face or edge
         bool is_hanging = false;
         for(int face = 0; face < Def::d()->num_faces; face++) {
            if(hanging_info.is_face_hanging(face) && ref_elem.face_contains(face, lnode)) {
//               cout << "is hanging face " << face << endl;
               is_hanging = true;
               break;
            }
         }
         if(! is_hanging) {
            for(int edge = 0; edge < Def::d()->num_edges; edge++) {
               if(hanging_info.is_edge_hanging(edge) && ref_elem.edge_contains(edge, lnode)) {
//                  cout << "is hanging edge " << edge << endl;
                  is_hanging = true;
                  break;
               }
            }
         }
         // this node is not hanging, use its coordinates
         for(int dim = 0; dim < Def::d()->num_dim; dim++) {
            double coord;
            if (is_hanging)
               coord = parent_nodes_coords[lnode][dim];
            else
               coord = nodes_coords[lnode][dim];


            //   cout << coords.val[dim][node_idx] << " <-> " << nodes_coords[lnode][dim] << endl;
            assert((coords.val[dim][node_idx] == FORBIDDEN) || reals_equal(coords.val[dim][node_idx], coord));
            coords.val[dim][node_idx] = coord;
         }
      }
   }

   // check that all values have been assigned
   for(int node = 0; node < subdomain_dims->n_nodes; node++) {
      for(int dim = 0; dim < subdomain_dims->n_mesh_dims; dim++) {
         assert(coords.val[dim][node] != FORBIDDEN);
      }
   }
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



