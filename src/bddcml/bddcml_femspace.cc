#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <assert.h>
#include <math.h>

#include "bddcml/bddcml_mesh.h"
#include "bddcml/bddcml_femspace.h"

using namespace std;

BddcmlFemSpace::BddcmlFemSpace(const BddcmlMesh *mesh) : mesh(mesh) {
   subdomain_dims = mesh->subdomain_dims;
   allocate_idx_array(subdomain_dims->n_nodes, &node_num_dofs);
   allocate_idx_array(subdomain_dims->n_dofs, &dofs_global_map);
   allocate_idx_array(subdomain_dims->n_dofs, &fixs_code);
   allocate_real_array(subdomain_dims->n_dofs, &fixs_values);
}

BddcmlFemSpace::~BddcmlFemSpace() {
   free_idx_array(&node_num_dofs);
   free_idx_array(&dofs_global_map);
   free_idx_array(&fixs_code);
   free_real_array(&fixs_values);
}

void BddcmlFemSpace::print(int which_rank) const {
   print_rank = which_rank;
   PPP printf("\n*************** BEGIN BDDCML FEM SPACE ************************\n");
   PPP printf("dofs: %d\n", subdomain_dims->n_dofs);
   assert(subdomain_dims->n_dofs == subdomain_dims->n_nodes);
   for(int node = 0; node < subdomain_dims->n_dofs; node++) {
      assert(node_num_dofs.val[node] == 1);
      PPP printf("node %d -> num dofs %d, gdofs: (%d ), bc code %d, bc val %6.4lf, coords (%4.3lf, %4.3lf)\n", node, node_num_dofs.val[node],
                 dofs_global_map.val[node], fixs_code.val[node], fixs_values.val[node],
                 mesh->coords.val[0][node], mesh->coords.val[1][node]);
   }
   PPP printf("*************** END BDDCML FEM SPACE ************************\n\n");
}

const double EPS = 1e-10;
bool real_equal(real a, real b) {
   return fabs(a - b) < EPS;
}

void BddcmlFemSpace::prepare_subdomain_fem_space(PhysicsType physicsType, exact_fn dirichlet_bc_exact) {
   this->physicsType = physicsType;
   int num_dofs_per_node = mesh->subdomain_dims->n_node_dofs;
   for(int node = 0; node < mesh->subdomain_dims->n_nodes; node++) {
      bool is_on_boundary = ((real_equal(mesh->coords.val[0][node], 0.0)) || (real_equal(mesh->coords.val[0][node], 1.0))
            || (real_equal(mesh->coords.val[1][node], 0.0)) || (real_equal(mesh->coords.val[1][node], 1.0)));

      if(subdomain_dims->n_mesh_dims == 3) {
         is_on_boundary = is_on_boundary || ((real_equal(mesh->coords.val[2][node], 0.0)) || (real_equal(mesh->coords.val[2][node], 1.0)));
      }

      node_num_dofs.val[node] = num_dofs_per_node;

      // todo: this should be alligned with NodalElementMesh creation...
      for(int local_dof_idx = 0; local_dof_idx < num_dofs_per_node; local_dof_idx++) {
         int subdomain_dof = num_dofs_per_node * node + local_dof_idx;
         int global_dof = num_dofs_per_node * mesh->node_global_map.val[node] + local_dof_idx;

         dofs_global_map.val[subdomain_dof] = global_dof;

         if(is_on_boundary) {
            double dirichlet_value = 0.0;
            if(dirichlet_bc_exact != nullptr) {
               // this is hack... improve...
               assert(physicsType == PhysicsType::LAPLACE);
               vector<double> coords(Def::d()->num_dim, 0.0);
               for(int dim = 0; dim < Def::d()->num_dim; dim++) {
                  coords[dim] = mesh->coords.val[dim][node];
               }
               vector<double> result(1);
               dirichlet_bc_exact(coords, &result);
               dirichlet_value = result[0];
            }
            fixs_code.val[subdomain_dof] = 1;
            fixs_values.val[subdomain_dof] = dirichlet_value;
         }
         else {
            fixs_code.val[subdomain_dof] = 0;
            fixs_values.val[subdomain_dof] = 0.0;
         }

      }

   }
}
