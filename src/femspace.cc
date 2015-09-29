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
#include "p4est_common.h"
#include "mesh.h"
#include "femspace.h"


void init_fem_space(BddcmlDimensions* dims, BddcmlFemSpace* femsp)
{
   femsp->subdomain_dims = dims;
   allocate_idx_array(dims->n_nodes, &femsp->node_num_dofs);
   allocate_idx_array(dims->n_dofs, &femsp->dofs_global_map);
   allocate_idx_array(dims->n_dofs, &femsp->fixs_code);
   allocate_real_array(dims->n_dofs, &femsp->fixs_values);
}

void free_fem_space(BddcmlFemSpace* femsp)
{
   free_idx_array(&femsp->node_num_dofs);
   free_idx_array(&femsp->dofs_global_map);
   free_idx_array(&femsp->fixs_code);
   free_real_array(&femsp->fixs_values);
}

void print_bddcml_fem_space(BddcmlFemSpace *femsp, BddcmlMesh *mesh, int which_rank)
{
   print_rank = which_rank;
   PPP printf("\n*************** BEGIN BDDCML FEM SPACE ************************\n");
   PPP printf("dofs: %d\n", femsp->subdomain_dims->n_dofs);
   assert(femsp->subdomain_dims->n_dofs == femsp->subdomain_dims->n_nodes);
   for(int node = 0; node < femsp->subdomain_dims->n_dofs; node++)
   {
      assert(femsp->node_num_dofs.val[node] == 1);
      PPP printf("node %d -> num dofs %d, gdofs: (%d ), bc code %d, bc val %6.4lf, coords (%4.3lf, %4.3lf)\n", node, femsp->node_num_dofs.val[node],
                 femsp->dofs_global_map.val[node], femsp->fixs_code.val[node], femsp->fixs_values.val[node],
                 mesh->coords.val[0][node], mesh->coords.val[1][node]);
   }
   PPP printf("*************** END BDDCML FEM SPACE ************************\n\n");
}

const double EPS = 1e-10;
bool real_equal(real a, real b)
{
   return fabs(a - b) < EPS;
}

void prepare_subdomain_fem_space(BddcmlMesh *mesh, BddcmlFemSpace *femsp, PhysicsType physicsType)
{
   init_fem_space(mesh->subdomain_dims, femsp);
   femsp->physicsType = physicsType;
   int num_dofs_per_node = mesh->subdomain_dims->n_node_dofs;
   for(int node = 0; node < mesh->subdomain_dims->n_nodes; node++)
   {
      bool is_on_boundary = ((real_equal(mesh->coords.val[0][node], 0.0)) || (real_equal(mesh->coords.val[0][node], 1.0))
            || (real_equal(mesh->coords.val[1][node], 0.0)) || (real_equal(mesh->coords.val[1][node], 1.0))
#ifdef P4_TO_P8
            || (real_equal(mesh->coords.val[2][node], 0.0)) || (real_equal(mesh->coords.val[2][node], 1.0))
#endif
            );

      femsp->node_num_dofs.val[node] = num_dofs_per_node;

      for(int local_dof_idx = 0; local_dof_idx < num_dofs_per_node; local_dof_idx++)
      {
         int subdomain_dof = num_dofs_per_node * node + local_dof_idx;
         int global_dof = num_dofs_per_node * mesh->node_global_map.val[node] + local_dof_idx;

         femsp->dofs_global_map.val[subdomain_dof] = global_dof;

         if(is_on_boundary)
         {
            femsp->fixs_code.val[subdomain_dof] = 1;
            femsp->fixs_values.val[subdomain_dof] = 0.0;
         }
         else
         {
            femsp->fixs_code.val[subdomain_dof] = 0;
            femsp->fixs_values.val[subdomain_dof] = 0.0;
         }

      }

   }
}
