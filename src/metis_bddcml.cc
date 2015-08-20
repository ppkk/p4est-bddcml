#ifndef P4_TO_P8
#include <p4est_bits.h>
#include <p4est_ghost.h>
#include <p4est_lnodes.h>
#include <p4est_vtk.h>
#include <p4est_iterate.h>
#else
#include <p8est_bits.h>
#include <p8est_ghost.h>
#include <p8est_lnodes.h>
#include <p8est_vtk.h>
#include <p8est_iterate.h>
#endif

#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <assert.h>

#include <metis.h>

#include "helpers.h"
#include "bddcml_structs.h"
#include "p4est_common.h"
#include "p4est_bddcml_interaction.h"

const int degree = 1;

// this is the main difference of the metis example
// p4est is global - not distributed.
// thus all indices (node and element) in p4est are global
int *elems_glob_to_loc;
int *nodes_glob_to_loc;



#ifndef P4_TO_P8
   const int MAX_GRAPH_NEIGHBOURS = 16-4;
#else
   const int MAX_GRAPH_NEIGHBOURS = 64-8;
#endif

idx_t *graph_xadj, *graph_adjncy;

void add_graph_edge(idx_t first, idx_t second)
{
   //printf("adding graph edge %d -> %d\n", first, second);

   // we do not store a->a
   if(first == second)
      return;

   // check whether has been allready added (when considering connections through vertices, we try to add them multiple times
   for(int i = 0; i < graph_xadj[first]; i++)
   {
      if(graph_adjncy[first*MAX_GRAPH_NEIGHBOURS + i] == second)
         return;
   }

   assert(graph_xadj[first] < MAX_GRAPH_NEIGHBOURS);
   graph_adjncy[first*MAX_GRAPH_NEIGHBOURS + graph_xadj[first]] = second;
   graph_xadj[first]++;
}

void corners_action (p4est_iter_corner_info_t * info, void *user_data)
{
   sc_array_t *sides = &(info->sides);

#ifndef P4_TO_P8
   const int MAX_CORNER_ADJ = 4;
#else
   const int MAX_CORNER_ADJ = 8;
#endif

   int elem_num = sides->elem_count;
   int elem_idx[MAX_CORNER_ADJ];

   assert(sides->elem_count <= MAX_CORNER_ADJ);
   if(sides->elem_count > 1)
   {
      for(unsigned int side_idx = 0; side_idx < sides->elem_count; side_idx++)
      {
         p4est_iter_corner_side_t *side = p4est_iter_cside_array_index_int (sides, side_idx);
         assert(!side->is_ghost);
         elem_idx[side_idx] = side->quadid;
      }
   }

   for(int i = 0; i < elem_num; i++)
   {
      for(int j = 0; j < elem_num; j++)
      {
         add_graph_edge(elem_idx[i], elem_idx[j]);
      }
   }
}

void face_action (p4est_iter_face_info_t * info, void *user_data)
{
   //printf("face action \n");
   sc_array_t *sides = &(info->sides);

   int elem_nums[2]; // number of elements on each side (1-2, for reg/hanging)
   int elem_idx[2][2]; // element(s) on each side

   if(sides->elem_count > 1)
   {
      assert(sides->elem_count == 2);
      // inner face
      for(int side_idx = 0; side_idx < 2; side_idx++)
      {
         p4est_iter_face_side_t *side = p4est_iter_fside_array_index_int (sides, side_idx);
         // for more initial elements (trees) we need to convert to global indices
         assert(side->treeid == 0);

         if(side->is_hanging)
         {
            // not ready for distributed yet
            assert(!side->is.hanging.is_ghost[0]);
            assert(!side->is.hanging.is_ghost[1]);
            elem_nums[side_idx] = 2;
            elem_idx[side_idx][0] = side->is.hanging.quadid[0];
            elem_idx[side_idx][1] = side->is.hanging.quadid[1];
         }
         else
         {
            // not ready for distributed yet
            assert(!side->is.full.is_ghost);
            elem_nums[side_idx] = 1;
            elem_idx[side_idx][0] = side->is.full.quadid;
         }
      }
      for(int i = 0; i < elem_nums[0]; i++)
      {
         for(int j = 0; j < elem_nums[1]; j++)
         {
            add_graph_edge(elem_idx[0][i], elem_idx[1][j]);
            add_graph_edge(elem_idx[1][j], elem_idx[0][i]);
         }
      }
   }


}


void print_graph(int num_vert)
{
   printf("\n************ BEGIN GRAPH *********************\n");
   printf("xadj (%d) : ", num_vert+1);
   for(int i = 0; i <= num_vert; i++)
      printf("%d, ", graph_xadj[i]);
   printf("\n");
   printf("adjncy (%d) : ", MAX_GRAPH_NEIGHBOURS * num_vert);
   for(int i = 0; i < MAX_GRAPH_NEIGHBOURS * num_vert; i++)
      printf("%d, ", graph_adjncy[i]);
   printf("\n\n");
   for(int vert = 0; vert < num_vert; vert++)
   {
      printf("%d : (%d) -> ", vert, graph_xadj[vert]);
      for(int neighbour = graph_xadj[vert]; neighbour < graph_xadj[vert+1]; neighbour++)
      {
         printf("%d, ", graph_adjncy[neighbour]);
      }
      printf("\n");
   }
   printf("************   END GRAPH *********************\n\n");

}

enum ConnectType{
   CONNECT_CORNERS,
   CONNECT_FACES
};

void create_graph(p4est_t *p4est, p4est_lnodes_t *lnodes, idx_t nparts, mptype_et metis_algorithm, ConnectType connect_type, idx_t* part)
{
   // vertices of the graph are mesh elements
   idx_t num_vert = lnodes->num_local_elements;
   graph_xadj = (idx_t*) malloc((num_vert + 1) * sizeof(idx_t));
   memset(graph_xadj, 0, (num_vert + 1) * sizeof(idx_t));
   graph_adjncy = (idx_t*) malloc(MAX_GRAPH_NEIGHBOURS * num_vert * sizeof(idx_t));
   memset(graph_adjncy, 0, MAX_GRAPH_NEIGHBOURS * num_vert * sizeof(idx_t));


   if(connect_type == CONNECT_CORNERS)
   {
      p4est_iterate (p4est,                 /* the forest */
                     NULL,                 /* the ghost layer */
                     NULL,   /* the synchronized ghost data */
                     NULL,
                     NULL,
#ifdef P4_TO_P8
                     NULL,           // edges
#endif
                     corners_action);          // corners

   }
   else if(connect_type == CONNECT_FACES)
   {
      p4est_iterate (p4est,                 /* the forest */
                     NULL,                 /* the ghost layer */
                     NULL,   /* the synchronized ghost data */
                     NULL,
                     face_action,
#ifdef P4_TO_P8
                     NULL,           // edges
#endif
                     NULL);          // corners
   }
   else
   {
      assert(0);
   }

   // the following estimates might not work for more complex geometries (more than 4/8 elements meeting in 1 point)
   assert(p4est->trees->elem_count == 1);

   //print_graph(num_vert);

   // compress graph_adjncy
   int new_position = 0;
   for(int vert = 0; vert < num_vert; vert++)
   {
      for(int vert_edge = 0; vert_edge < graph_xadj[vert]; vert_edge++)
      {
         int old_position = MAX_GRAPH_NEIGHBOURS * vert + vert_edge;
         //printf("moving from %d to %d\n", old_position, new_position);
         graph_adjncy[new_position] = graph_adjncy[old_position];
         new_position++;
      }
   }

   //now graph_xadj has to be changed to running sum
   for(int i = num_vert; i > 0; i--)
      graph_xadj[i] = graph_xadj[i-1];
   graph_xadj[0] = 0;
   for(int i = 0; i < num_vert; i++)
      graph_xadj[i+1] += graph_xadj[i];

   //print_graph(num_vert);

   idx_t options[METIS_NOPTIONS];
   idx_t ncon = 1;
   idx_t objval;
   METIS_SetDefaultOptions(options);

//   gk_malloc_init();
//   real_t part_timer;
//   gk_startcputimer(part_timer);

   int status;
   switch (metis_algorithm) {
     case METIS_PTYPE_RB:
       status = METIS_PartGraphRecursive(&num_vert, &ncon, graph_xadj, graph_adjncy, NULL, NULL, NULL,
                    &nparts, NULL, NULL, options, &objval, part);
       break;

     case METIS_PTYPE_KWAY:
       status = METIS_PartGraphKway(&num_vert, &ncon, graph_xadj, graph_adjncy, NULL, NULL, NULL,
                    &nparts, NULL, NULL, options, &objval, part);
       break;

   }

//   gk_stopcputimer(part_timer);

//   if (gk_GetCurMemoryUsed() != 0)
//     printf("***It seems that Metis did not free all of its memory! Report this.\n");
//   size_t maxmemory = gk_GetMaxMemoryUsed();
//   gk_malloc_cleanup(0);


   if (status == METIS_OK)
     printf("\n***Metis OK.\n");
   else if (status == METIS_ERROR_INPUT)
      printf("\n***Metis input error.\n");
   else if (status == METIS_ERROR_MEMORY)
      printf("\n***Metis memory error.\n");
   else if (status == METIS_ERROR)
      printf("\n***Metis some other error.\n");
   else
      assert(0);

//   for(int i = 0; i < num_vert; i++)
//      printf("%d, ", part[i]);
//   printf("\n");

   free(graph_xadj);
   free(graph_adjncy);
}

void prepare_subdomain_data(p4est_t *p4est, p4est_lnodes_t *lnodes, int *metis_part, BddcmlDimensions *subdomain_dims,
                       BddcmlDimensions *global_dims, BddcmlMesh* mesh, real** element_volumes, sc_MPI_Comm mpicomm)
{
   double              vxyz[3];  /* We embed the 2D vertices into 3D space. */
   p4est_quadrant_t   sp, node;

   // here the local p4est is the global mesh...
#ifndef P4_TO_P8
   init_dimmensions(subdomain_dims, 2);
   init_dimmensions(global_dims, 2);
#else
   init_dimmensions(subdomain_dims, 3);
   init_dimmensions(global_dims, 3);
#endif
   global_dims->n_nodes = lnodes->num_local_nodes;
   global_dims->n_dofs = lnodes->num_local_nodes;
   global_dims->n_elems = lnodes->num_local_elements;


   subdomain_dims->n_nodes = 0;
   subdomain_dims->n_dofs  = 0;
   subdomain_dims->n_elems = 0;

   elems_glob_to_loc = (int*) malloc(global_dims->n_elems * sizeof(int));
   for(int gelem = 0; gelem < global_dims->n_elems; gelem++)
   {
      if(metis_part[gelem] == mpi_rank)
      {
         elems_glob_to_loc[gelem] = subdomain_dims->n_elems;
         subdomain_dims->n_elems++;
      }
      else
      {
         elems_glob_to_loc[gelem] = -1;
      }
   }

   // upper bound for initialization
   // todo: we waste memory, but this is just a simple test with metis...
   subdomain_dims->n_nodes = P4EST_CHILDREN * subdomain_dims->n_elems;
   init_mesh(subdomain_dims, mesh);
   subdomain_dims->n_nodes = 0;

   // I do not have the right dimmensions yet
   free_real_2D_array(&mesh->coords);

    *element_volumes = (real*) malloc(subdomain_dims->n_elems * sizeof(real));

   nodes_glob_to_loc = (int*) malloc(global_dims->n_nodes * sizeof(int));
   for(int i = 0; i < global_dims->n_nodes; i++)
      nodes_glob_to_loc[i] = -1;

   p4est_locidx_t glob_quad_idx = 0;
   p4est_quadrant_t *quad;

   for_all_quads(p4est, glob_quad_idx, quad)
   {
      int loc_quad_idx = elems_glob_to_loc[glob_quad_idx];
      if(loc_quad_idx == -1)
         continue;

      mesh->elem_global_map.val[loc_quad_idx] = glob_quad_idx;
      mesh->num_nodes_of_elem.val[loc_quad_idx] = P4EST_CHILDREN;

      for (int lnode = 0; lnode < P4EST_CHILDREN; ++lnode) {
         p4est_locidx_t glob_node_idx = lnodes->element_nodes[P4EST_CHILDREN * glob_quad_idx + lnode];
         if(nodes_glob_to_loc[glob_node_idx] == -1)
         {
            nodes_glob_to_loc[glob_node_idx] = subdomain_dims->n_nodes;
            mesh->node_global_map.val[subdomain_dims->n_nodes] = glob_node_idx;
            subdomain_dims->n_nodes++;
         }
         int loc_node_idx = nodes_glob_to_loc[glob_node_idx];
         mesh->elem_node_indices.val[P4EST_CHILDREN * loc_quad_idx + lnode] = loc_node_idx;
      }

      (*element_volumes)[loc_quad_idx] = pow(0.5, quad->level * mesh->subdomain_dims->n_problem_dims);

   }}} //for all quads


   // fix the dimmensions and mesh allocation
   subdomain_dims->n_dofs  = subdomain_dims->n_nodes;
   mesh->node_global_map.len = subdomain_dims->n_nodes;
   allocate_real_2D_array(subdomain_dims->n_nodes, subdomain_dims->n_problem_dims, &mesh->coords);


   for_all_quads(p4est, glob_quad_idx, quad)
   {
      int loc_quad_idx = elems_glob_to_loc[glob_quad_idx];
      if(loc_quad_idx == -1)
         continue;
      /* Figure out the hanging corners on this element, if any. */
      int hanging_corner[P4EST_CHILDREN];
      int anyhang = lnodes_decode2 (lnodes->face_code[glob_quad_idx], hanging_corner);

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
         p4est_locidx_t glob_node_idx = lnodes->element_nodes[P4EST_CHILDREN * glob_quad_idx + lnode];
         int loc_node_idx = nodes_glob_to_loc[glob_node_idx];
         P4EST_ASSERT (glob_node_idx >= 0 && glob_node_idx < global_dims->n_nodes);
         P4EST_ASSERT (loc_node_idx  >= 0 && loc_node_idx  < subdomain_dims->n_nodes);
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

         mesh->coords.val[0][loc_node_idx] = vxyz[0];
         mesh->coords.val[1][loc_node_idx] = vxyz[1];
#ifdef P4_TO_P8
         mesh->coords.val[2][loc_node_idx] = vxyz[2];
#endif



//         PPP printf("(%3.2lf, %3.2lf), ", vxyz[0], vxyz[1]);

      }
//      PPP printf("\n");

   }}}//for all quads

}

// a variant of the same function from general use
// the difference is, that the indices in p4est are GLOBAL here and have to be tramsformed to local
void assemble_matrix_rhs_metis(p4est_lnodes_t *lnodes, BddcmlMesh *mesh, double *element_volumes, BddcmlFemSpace *femsp,
                            SparseMatrix *matrix, RealArray *rhss)
{
   real i_coeffs, j_coeffs;
   real mass_dd[P4EST_CHILDREN][P4EST_CHILDREN];
   real stiffness_dd[P4EST_CHILDREN][P4EST_CHILDREN];
   generate_reference_matrices(stiffness_dd, mass_dd);

   int element_offset = 0;
   for(int loc_elem_idx = 0; loc_elem_idx < mesh->subdomain_dims->n_elems; loc_elem_idx++) {
      int glob_elem_idx = mesh->elem_global_map.val[loc_elem_idx];
      int num_nodes_of_elem = mesh->num_nodes_of_elem.val[loc_elem_idx];
      int ndof_per_element = num_nodes_of_elem;
      assert(ndof_per_element == P4EST_CHILDREN);

      // TODO: elem_size and elem_volume is correct only when the mesh is obtained by refinements
      // from a UNIT SQUARE/CUBE
      // TODO: ONLY FROM ****UNIT****
      real elem_volume = element_volumes[loc_elem_idx];

      double reference_scaled =
      #ifndef P4_TO_P8
            1.;
#else
            pow(elem_volume, 1./(double)mesh->subdomain_dims->n_problem_dims);
#endif
      for(int j = 0; j < ndof_per_element; j++) {

         //todo: dofs should be taken from femsp!
         int jdof = mesh->elem_node_indices.val[element_offset + j];
         assert(femsp->node_num_dofs.val[jdof] == 1);

         p4est_locidx_t j_nodes[4];
         int j_nindep = independent_nodes(lnodes, glob_elem_idx, j, j_nodes, &j_coeffs);
         if(j_nindep == 1)
         {
            assert(jdof == nodes_glob_to_loc[j_nodes[0]]);
            assert(mesh->node_global_map.val[jdof] == j_nodes[0]);
         }
         for(int j_indep_nodes_idx = 0; j_indep_nodes_idx < j_nindep; j_indep_nodes_idx++)
         {
            int j_indep_node_glob = j_nodes[j_indep_nodes_idx];
            int j_indep_node_loc = nodes_glob_to_loc[j_indep_node_glob];

            for(int i = 0; i < ndof_per_element /*<= j*/; i++) {
               int idof = mesh->elem_node_indices.val[element_offset + i];

               p4est_locidx_t i_nodes[4];
               int i_nindep = independent_nodes(lnodes, glob_elem_idx, i, i_nodes, &i_coeffs);

               if(i_nindep == 1)
               {
                  assert(idof == nodes_glob_to_loc[i_nodes[0]]);
                  assert(mesh->node_global_map.val[idof] == i_nodes[0]);
               }
               for(int i_indep_nodes_idx = 0; i_indep_nodes_idx < i_nindep; i_indep_nodes_idx++)
               {
                  int i_indep_node_glob = i_nodes[i_indep_nodes_idx];
                  int i_indep_node_loc = nodes_glob_to_loc[i_indep_node_glob];

                  double matrix_value = i_coeffs * j_coeffs * reference_scaled * stiffness_dd[j][i];
                  add_matrix_entry(matrix, i_indep_node_loc, j_indep_node_loc, matrix_value);
                  //                  printf("adding entry loc (%d, %d), nodes, (%d, %d), coefs (%3.2lf, %3.2lf), number indep (%d, %d), locstiff %lf, value %lf\n",
                  //                         j, i, j_indep_node, i_indep_node, j_coeffs, i_coeffs, j_nindep, i_nindep, stiffness_dd[j][i], matrix_value);
               }
            }

            // TODO: integrate properly
            // TODO: elem_volume is correct only when the mesh is obtained by refinements
            // from a UNIT SQUARE/CUBE
            // TODO: ONLY FROM ****UNIT****


            double rhs_value = j_coeffs * 1./(real)P4EST_CHILDREN * elem_volume * 1;
            rhss->val[j_indep_node_loc] += rhs_value;

         }
      }
      element_offset += num_nodes_of_elem;
   }

}

void plot_solution_metis(p4est_t * p4est, p4est_lnodes_t * lnodes, double* u_sol, double* u_exact, int* partition)
{
   p4est_topidx_t      tt;       /* Connectivity variables have this type. */
   p4est_locidx_t      glob_quad_idx, q, Q, node_total;  /* Process-local counters have this type. */
   p4est_locidx_t      lni;      /* Node index relative to this processor. */
   p4est_tree_t       *tree;     /* Pointer to one octree */
   sc_array_t         *tquadrants;       /* Quadrant array for one tree */
   int                 i;
   double              loc_vertex_values_sol[P4EST_CHILDREN], loc_vertex_values_exact[P4EST_CHILDREN];


   /* Write the forest to disk for visualization, one file per processor. */
   double *u_interp_sol = NULL;
   double *u_interp_exact = NULL;
   double *interp_partition = NULL;
   if(u_sol)
      u_interp_sol = P4EST_ALLOC (double, p4est->local_num_quadrants * P4EST_CHILDREN);
   if(u_exact)
      u_interp_exact = P4EST_ALLOC (double, p4est->local_num_quadrants * P4EST_CHILDREN);
   if(partition)
      interp_partition = P4EST_ALLOC (double, p4est->local_num_quadrants * P4EST_CHILDREN);

   /* Loop over local quadrants to apply the element matrices. */
   for (tt = p4est->first_local_tree, glob_quad_idx = 0, node_total = 0;
        tt <= p4est->last_local_tree; ++tt) {
      tree = p4est_tree_array_index (p4est->trees, tt);
      tquadrants = &tree->quadrants;
      Q = (p4est_locidx_t) tquadrants->elem_count;

      for (q = 0; q < Q; ++q, ++glob_quad_idx) {
//         int loc_quad_idx = elems_glob_to_loc[glob_quad_idx];

         for (i = 0; i < P4EST_CHILDREN; ++i) {
            lni = lnodes->element_nodes[P4EST_CHILDREN * glob_quad_idx + i];
            if(u_sol)
               loc_vertex_values_sol[i] = u_sol[lni];
            if(u_exact)
               loc_vertex_values_exact[i] = u_exact[lni];
         }

         if(u_sol)
            interpolate_hanging_nodes (lnodes->face_code[glob_quad_idx], loc_vertex_values_sol);
         if(u_exact)
            interpolate_hanging_nodes (lnodes->face_code[glob_quad_idx], loc_vertex_values_exact);

         for (i = 0; i < P4EST_CHILDREN; ++i) {
            if(u_sol)
               u_interp_sol[node_total]   = loc_vertex_values_sol[i];
            if(u_exact)
               u_interp_exact[node_total] = loc_vertex_values_exact[i];
            if(partition)
               interp_partition[node_total] = partition[glob_quad_idx];
            ++node_total;
         }
      }
   }

   if(u_sol)
   {
      if(partition)
      {
         p4est_vtk_write_all (p4est, NULL, 0.99999, 1, 1, 1, 0, 2, 0, "output",
                              "solution", u_interp_sol, "metis_part", interp_partition);
      }
      else
      {
         p4est_vtk_write_all (p4est, NULL, 0.99999, 1, 1, 1, 0, 1, 0, "output",
                              "solution", u_interp_sol);
      }
   }
   else
   {
      if(partition)
      {
         p4est_vtk_write_all (p4est, NULL, 0.99999, 1, 1, 1, 0, 1, 0, "output",
                              "metis_part", interp_partition);
      }
      else
      {
         p4est_vtk_write_all (p4est, NULL, 0.99999, 1, 1, 1, 0, 0, 0, "output");
      }
   }

   if(u_sol)
      P4EST_FREE (u_interp_sol);
   if(u_exact)
      P4EST_FREE (u_interp_exact);
   if(partition)
      P4EST_FREE (interp_partition);
}


int main (int argc, char **argv)
{
   init_corner_to_hanging();

   int mpiret = sc_MPI_Init (&argc, &argv);
   SC_CHECK_MPI (mpiret);
   sc_MPI_Comm mpicomm_bddcml = sc_MPI_COMM_WORLD;
   sc_MPI_Comm mpicomm_p4est = sc_MPI_COMM_SELF;
   mpiret = sc_MPI_Comm_rank(mpicomm_bddcml, &mpi_rank);
   mpiret = sc_MPI_Comm_size(mpicomm_bddcml, &mpi_size);

   /* These functions are optional.  If called they store the MPI rank as a
   * static variable so subsequent global p4est log messages are only issued
   * from processor zero.  Here we turn off most of the logging; see sc.h. */
   sc_init (mpicomm_p4est, 1, 1, NULL, SC_LP_ESSENTIAL);
   p4est_init (NULL, SC_LP_PRODUCTION);  /* SC_LP_ERROR for silence. */
   P4EST_GLOBAL_PRODUCTIONF
         ("This is the p4est %dD demo example/steps/%s_step4\n",
          P4EST_DIM, P4EST_STRING);

   // WARNING: integration is based on the fact, that the domain is unit square or cube!
   // WARNING: if the domain is changed, so has to be the integration!
   // TODO: do it properly
#ifndef P4_TO_P8
   p4est_connectivity_t *conn = p4est_connectivity_new_unitsquare ();
#else
   p4est_connectivity_t *conn = p8est_connectivity_new_unitcube ();
#endif

   BddcmlLevelInfo level_info;
   // Number of elements in an edge of a subdomain and number of subdomains in an edge of the unit cube
   if(argc == 1 + 1) {
      level_info.nlevels = atoi(argv[1]);
   }
   else {
      if ( mpi_rank == 0 ) {
         printf(" Usage: mpirun -np X ./p4est_bddcml NLEVELS\n");
      }
      exit(0);
   }

   // number of subdomains == mpi_size
   init_levels(mpi_size, &level_info);

   BddcmlGeneralParams general_params;
   set_implicit_general_params(&general_params);
   //general_params.just_direct_solve_int = 1;

   BddcmlKrylovParams krylov_params;
   set_implicit_krylov_params(&krylov_params);

   BddcmlPreconditionerParams preconditioner_params;
   set_implicit_preconditioner_params(&preconditioner_params);

   p4est_t *p4est = p4est_new (mpicomm_p4est, conn, 0, NULL, NULL);

   // pro metis: 3, 6, 8, np=10
   refine_and_partition(p4est, 3, refine_uniform);
   refine_and_partition(p4est, 6, refine_circle);
   refine_and_partition(p4est, 8, refine_square);
   refine_and_partition(p4est, 0, refine_point);
   refine_and_partition(p4est, 0, refine_diagonal);

   /* Create the ghost layer to learn about parallel neighbors. */
   p4est_ghost_t *ghost = p4est_ghost_new (p4est, P4EST_CONNECT_FULL);

   /* Create a node numbering for continuous linear finite elements. */
   p4est_lnodes_t *lnodes = p4est_lnodes_new (p4est, ghost, degree);

   /* Destroy the ghost structure -- no longer needed after node creation. */
   p4est_ghost_destroy (ghost);
   ghost = NULL;

   idx_t *metis_part = NULL;
   metis_part = (idx_t*) malloc(lnodes->num_local_elements * sizeof(idx_t));
   create_graph(p4est, lnodes, mpi_size, METIS_PTYPE_KWAY, CONNECT_FACES, metis_part);
   if(mpi_rank == 0)
   {
      plot_solution(p4est, lnodes, NULL, NULL, metis_part);
   }


   BddcmlDimensions subdomain_dims, global_dims;
   BddcmlMesh mesh;

   int print_rank_l = 0;

   // TODO: elem_volume is correct only when the mesh is obtained by refinements
   // from a UNIT SQUARE/CUBE
   real *element_volumes = NULL;

   // todo: using MPI Bcast in the following, should be possible to do without
   prepare_subdomain_data(p4est, lnodes, metis_part, &subdomain_dims, &global_dims, &mesh, &element_volumes, mpicomm_bddcml);

   //print_p4est_mesh(p4est, lnodes, print_rank_l);
   //print_bddcml_mesh(&mesh, print_rank_l);


   BddcmlFemSpace femsp;
   prepare_subdomain_fem_space(&mesh, &femsp);
   //print_bddcml_fem_space(&femsp, &mesh, print_rank_l);

   //plot_solution(p4est, lnodes, NULL, NULL, metis_part);

   print_rank = print_rank_l;
   print_basic_properties(&global_dims, mpi_size, &level_info, &krylov_params);
   PPP printf("Initializing BDDCML ...");
   // tell me how much subdomains should I load
   level_info.nsub_loc_1 = -1;

   bddcml_init(&general_params, &level_info, mpicomm_bddcml);
   // should be 1 subdomain per processor
   assert(level_info.nsub_loc_1 == 1);

   mpiret = MPI_Barrier(mpicomm_bddcml);
   PPP printf("Initializing BDDCML done.\n");

   RealArray rhss;
   allocate_real_array(subdomain_dims.n_dofs, &rhss);
   zero_real_array(&rhss);

   int is_rhs_complete = 0;

   RealArray sols;
   allocate_real_array(subdomain_dims.n_dofs, &sols);
   zero_real_array(&sols);

   SparseMatrix matrix;
   int ndof_per_element = P4EST_CHILDREN;
   // how much space the upper triangle of the element matrix occupies
   int lelm = ndof_per_element * (ndof_per_element + 1) / 2;

   MatrixType matrix_type = SPD;
   // todo: do it properly
   const int extra_space_for_hanging_nodes = 4 * (matrix_type == GENERAL ? 2 : 1);
   allocate_sparse_matrix(extra_space_for_hanging_nodes * subdomain_dims.n_elems*lelm, matrix_type, &matrix);
   zero_matrix(&matrix);

   assemble_matrix_rhs_metis(lnodes, &mesh, element_volumes, &femsp, &matrix, &rhss);
   //print_complete_matrix_rhs(&femsp, &global_dims, &matrix, &rhss, mpicomm);

   // user constraints - not really used here
   Real2DArray user_constraints;
   allocate_real_2D_array(0, 0, &user_constraints);

   // data for elements - not really used here
   Real2DArray element_data;
   allocate_real_2D_array(0, 0, &element_data);

   // data for dofs - not really used here
   RealArray dof_data;
   allocate_real_array(0, &dof_data);


   PPP printf("Loading data ...\n");

   int subdomain_idx = mpi_rank;
   bddcml_upload_subdomain_data(&global_dims, &subdomain_dims,
                                     subdomain_idx, &mesh, &femsp,
                                     &rhss, is_rhs_complete, &sols, &matrix,
                                     &user_constraints, &element_data,
                                     &dof_data, &preconditioner_params);

   PPP printf("Loading data done.\n");


   mpiret = MPI_Barrier(mpicomm_bddcml);

   PPP printf("Preconditioner set-up ...\n");

   // PRECONDITIONER SETUP
   mpiret = MPI_Barrier(mpicomm_bddcml);
   // TODO: call time_start
   bddcml_setup_preconditioner(matrix.type, &preconditioner_params);

   mpiret = MPI_Barrier(mpicomm_bddcml);
   // TODO: call time_end(t_pc_setup)

   PPP printf("Preconditioner set-up done.\n");



   PPP printf("Calling Krylov method ...\n");

   mpiret = MPI_Barrier(mpicomm_bddcml);
   // TODO: call time_start
   // call with setting of iterative properties

   BddcmlConvergenceInfo convergence_info;

//   real normRn_sol, normRn2, normRn2_loc, normRn2_sub;
//   real normL2_sol, normL2_loc, normL2_sub;
//   real normLinf_sol, normLinf_loc;

   bddcml_solve(&krylov_params, &convergence_info, mpicomm_bddcml);
   mpiret = MPI_Barrier(mpicomm_bddcml);

   // TODO: call time_end(t_krylov)

   PPP printf("Krylov method done.\n");

   PPP printf(" Output of PCG: ==============\n");
   PPP printf(" Number of iterations: %d\n", convergence_info.num_iter);
   PPP printf(" Convergence reason:   %d\n", convergence_info.converged_reason);
   if ( convergence_info.condition_number >= 0. ) {
      PPP printf(" Condition number: %lf\n", convergence_info.condition_number);
   }
   PPP printf(" =============================\n");


   bddcml_download_local_solution(subdomain_idx, &sols);
   //printf("mpirank %d, sols len %d\n", mpi_rank, sols.len);


   // communicate the solution for plotting
   int *recvbuff_lens = (int*) malloc(mpi_size * sizeof(int));
   int *displs = (int*) malloc(mpi_size * sizeof(int));
   assert(sols.len == subdomain_dims.n_dofs);
   MPI_Gather(&sols.len, 1, MPI_INT, recvbuff_lens, 1, MPI_INT, 0, MPI_COMM_WORLD);
   int recv_len = 0;
   if(mpi_rank == 0)
   {
      for(int i = 0; i < mpi_size; i++)
      {
         recv_len += recvbuff_lens[i];
      }

      displs[0] = 0;
      for(int i = 1; i < mpi_size; i++)
      {
         displs[i] = displs[i-1] + recvbuff_lens[i-1];
      }
   }
   double *recvbuff_sol = (double*) malloc(recv_len * sizeof(double));
   int *recvbuff_node_map = (int*) malloc(recv_len * sizeof(int));
   MPI_Gatherv(sols.val, sols.len, MPI_DOUBLE, recvbuff_sol, recvbuff_lens, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
   MPI_Gatherv(mesh.node_global_map.val, sols.len, MPI_INT, recvbuff_node_map, recvbuff_lens, displs, MPI_INT, 0, MPI_COMM_WORLD);

   double *complete_sol = (double*) malloc(global_dims.n_dofs * sizeof(double));

   const double NOTHING = 1e100;
   for(int i = 0; i < global_dims.n_dofs; i++)
      complete_sol[i] = NOTHING;

   if(mpi_rank == 0)
   {
      int recbuff_position = 0;
      for(int subd_idx = 0; subd_idx < mpi_size; subd_idx++)
      {
         for(int loc_dof = 0; loc_dof < recvbuff_lens[subd_idx]; loc_dof++)
         {
            int glob_dof = recvbuff_node_map[recbuff_position];
            if(complete_sol[glob_dof] != NOTHING)
            {
               assert(fabs(complete_sol[glob_dof] - recvbuff_sol[recbuff_position]) < 1e-8);
            }
            complete_sol[glob_dof] = recvbuff_sol[recbuff_position];
            recbuff_position++;
         }
      }
   }

   free(recvbuff_sol);
   free(recvbuff_node_map);
   free(recvbuff_lens);
   free(displs);

   if(mpi_rank == 0)
   {
//      for(int i = 0; i < global_dims.n_dofs; i++)
//         printf("globsol %d : %lf\n", i, complete_sol[i]);
      plot_solution_metis(p4est, lnodes, complete_sol, NULL, metis_part); //uexact_eval, NULL);
   }

   free(complete_sol);

   free_mesh(&mesh);
   free_fem_space(&femsp);

   free_real_array(&rhss);
   free_real_array(&sols);
   free_sparse_matrix(&matrix);

   free(element_volumes);

   free(metis_part);
   metis_part = NULL;
   free(nodes_glob_to_loc);
   free(elems_glob_to_loc);

   assert(get_num_allocations() == 0);

   /* Destroy the p4est and the connectivity structure. */
   p4est_lnodes_destroy (lnodes);
   p4est_destroy(p4est);
   p4est_connectivity_destroy(conn);

   /* Verify that allocations internal to p4est and sc do not leak memory.
   * This should be called if sc_init () has been called earlier. */
   sc_finalize ();


   /* This is standard MPI programs.  Without --enable-mpi, this is a dummy. */
   mpiret = sc_MPI_Finalize ();
   SC_CHECK_MPI (mpiret);
   return 0;
}

