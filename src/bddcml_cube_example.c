#include <stdio.h>
#include <stdlib.h>
#include "bddcml_structs.h"
#include "bddcml_cube_example.h"

//************************************************************************************************
// subroutine creating data for one subdomain
void prepare_subdomain_data(int isub, // global subdomain index
                            int num_sub_per_cube_edge, // number of subdomains in one edge of a cube
                            int num_el_per_sub_edge,  // number of elements in one edge of a subdomain
                            real hsize, // element size
                            BddcmlConnectivity *connectivity, int* nndfs, int lnndfs,
                            int* isegns, int lisegns, int* isngns, int lisngns, int* isvgvns, int lisvgvns,
                            real* xyzs, int lxyzs1, int lxyzs2,
                            int* ifixs, int lifixs, real* fixvs, int lfixvs)
{
   const char routine_name[] = "PREPARE_SUBDOMAIN_DATA";
   int num_sub_xy;
   int ind_sub_x, ind_sub_y, ind_sub_z;
   int num_el_per_cube_edge;
   int num_nodes_per_sub_edge, num_nodes_per_cube_edge;
   int i, j, k;
   int ig, jg, kg;
   int indng, indns;
   int indelg, indels;
   int indinets;
   int nne;
   int n1, n2, n3, n4, n5, n6, n7, n8;
   int idx;

   // number of elements on one edge of the cube
   num_el_per_cube_edge = num_el_per_sub_edge * num_sub_per_cube_edge;

   // determine subdomain indices along coordinate axes
   // intentional integer divisons
   num_sub_xy = num_sub_per_cube_edge * num_sub_per_cube_edge;

   // TODO: zkontrolovat +- 1 v nasledujicich 3
   ind_sub_z = isub / num_sub_xy;
   ind_sub_y = (isub - ind_sub_z * num_sub_xy)/num_sub_per_cube_edge;
   ind_sub_x =  isub - (ind_sub_z * num_sub_xy) - (ind_sub_y * num_sub_per_cube_edge);

   // debug

   printf("subdomain index and coord indices: %d, %d, %d, %d\n", isub, ind_sub_x, ind_sub_y, ind_sub_z);

   // initialize boundary conditions
   for(idx = 0; idx < lifixs; idx++) {
      ifixs[idx] = 0;
   }
   for(idx = 0; idx < lfixvs; idx++) {
      fixvs[idx] = 0.;
   }

   // number nodes and degrees of freedom
   num_nodes_per_sub_edge  = num_el_per_sub_edge + 1;
   num_nodes_per_cube_edge = num_el_per_cube_edge + 1;

   // local number of nodes
   int lxyz1 = num_nodes_per_sub_edge * num_nodes_per_sub_edge * num_nodes_per_sub_edge;

   indns = 0;
   for(k = 0; k < num_nodes_per_sub_edge; k++) {
      kg = ind_sub_z * (num_nodes_per_sub_edge - 1) + k;
      for(j = 0; j < num_nodes_per_sub_edge; j++) {
         jg = ind_sub_y * (num_nodes_per_sub_edge - 1) + j;
         for(i = 0; i < num_nodes_per_sub_edge; i++) {
            ig = ind_sub_x * (num_nodes_per_sub_edge - 1) + i;

            // compute global node index
            indng = ig + jg * num_nodes_per_cube_edge + kg * num_nodes_per_cube_edge * num_nodes_per_cube_edge;

            isngns[indns] = indng;

            // compute coordinates. In C interface, first all x, then all y, then all z.
            xyzs[indns          ] = ig * hsize;
            xyzs[indns +   lxyz1] = jg * hsize;
            xyzs[indns + 2*lxyz1] = kg * hsize;

            // for Poisson problem, there is only one dof per node,
            nndfs[indns] = 1;
            //and thus the numbering of nodes and dofs is the same,
            isvgvns[indns] = indng;

            // if node is on the boundary, fix boundary conditions
            if ( (ig == 0) || (ig == num_nodes_per_cube_edge-1) ||
                 (jg == 0) || (jg == num_nodes_per_cube_edge-1) ||
                 (kg == 0) || (kg == num_nodes_per_cube_edge-1))  {

               ifixs[indns] = 1;
               fixvs[indns] = 0.;
            }

            // increase counter of local nodes
            indns = indns + 1;
         }
      }
   }
   if (indns != lisngns) {
      printf("%s : Some bug in node index computing for sub %d\n", routine_name, isub);
      exit(0);
   }
   // debug
   // !write(*,*) 'isub',isub,'isngns',isngns


   // create element connectivity
   indels = 0;
   indinets = 0;
   nne = 8;
   for(k = 0; k < num_el_per_sub_edge; k++) {
      kg = ind_sub_z * num_el_per_sub_edge + k;
      for(j = 0; j < num_el_per_sub_edge; j++) {
         jg = ind_sub_y * num_el_per_sub_edge + j;
         for(i = 0; i < num_el_per_sub_edge; i++) {
            ig = ind_sub_x * num_el_per_sub_edge + i;

            // compute global element index
            indelg = ig + jg * num_el_per_cube_edge + kg * num_el_per_cube_edge * num_el_per_cube_edge;

            // compute local node index of the first node
            indns  = i + j * num_nodes_per_sub_edge + k * num_nodes_per_sub_edge * num_nodes_per_sub_edge;

            // compute indices of the eight nodes of each element
            n1 = indns;
            n2 = indns + 1;
            n3 = n2 + num_nodes_per_sub_edge;
            n4 = n3 - 1;
            n5 = n1 + num_nodes_per_sub_edge * num_nodes_per_sub_edge;
            n6 = n2 + num_nodes_per_sub_edge * num_nodes_per_sub_edge;
            n7 = n3 + num_nodes_per_sub_edge * num_nodes_per_sub_edge;
            n8 = n4 + num_nodes_per_sub_edge * num_nodes_per_sub_edge;

            connectivity->elem_node_indices[indinets + 0] = n1;
            connectivity->elem_node_indices[indinets + 1] = n2;
            connectivity->elem_node_indices[indinets + 2] = n3;
            connectivity->elem_node_indices[indinets + 3] = n4;
            connectivity->elem_node_indices[indinets + 4] = n5;
            connectivity->elem_node_indices[indinets + 5] = n6;
            connectivity->elem_node_indices[indinets + 6] = n7;
            connectivity->elem_node_indices[indinets + 7] = n8;

            indinets = indinets + nne;

            // number of nodes on element is constant for all elements
            connectivity->n_nodes_of_elem[indels] = nne;

            // embedding of local elements into global numbers
            isegns[indels] = indelg;

            // increase counter of local elements
            indels = indels + 1;
         }
      }
   }
   // debug
   //write(*,*) 'isub',isub,'isegns',isegns
   //write(*,*) 'isub',isub,'inets',inets
   //write(*,*) 'isub',isub,'xyzs',xyzs
   //write(*,*) 'isub',isub,'ifixs',ifixs
   //write(*,*) 'isub',isub,'fixvs',fixvs

}


