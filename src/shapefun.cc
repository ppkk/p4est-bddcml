#include <assert.h>
#include <iostream>

#include "shapefun.h"
#include "p4est/my_p4est_interface.h"
#include "quadrature.h"

using namespace std;

//void ref_value_1D(int loc_id_1d, double x, double elem_len, double *value, double *der) {
//   if(loc_id_1d == 0) {
//      *value = (1-x)/2.;
//      *der = -1/elem_len;
//   }
//   else if(loc_id_1d == 1) {
//      *value = (1+x)/2.;
//      *der = 1/elem_len;
//   }
//   else {
//      assert(0);
//   }
//}

void ReferenceElement::find_node_coords(const std::vector<double> &start, double element_len,
                                   std::vector<std::vector<double> > *coords) const {
   int num_points_1d = order + 1;
   double node_distance = element_len / order;
   vector<double> node(num_dim, 0.0);
   int difs[num_dim];
   for(difs[2] = 0; difs[2] < ((num_dim == 3) ? num_points_1d : 1); difs[2]++) {
      for(difs[1] = 0; difs[1] < num_points_1d; difs[1]++) {
         for(difs[0] = 0; difs[0] < num_points_1d; difs[0]++) {
            // now we know which of the dim^(order+1) points we want, construct its coordinates
            cout << difs[0] << ", " << difs[1] << ", " << difs[2] << endl;
            for(int dim = 0; dim < num_dim; dim++) {
               node[dim] = start[dim] + difs[dim] * node_distance;
            }
            coords->push_back(node);
         }
      }
   }

}

ReferenceElement::ReferenceElement(int num_dim, int order) : num_dim(num_dim), order(order)
{
   vector<double> starting_corner(num_dim, -1.0);
   // use for reference element [-1,1]^dim
   find_node_coords(starting_corner, 2., &node_coords);
}

void ReferenceElement::shape_fun_1d(int idx, double x, double *value, double* der) const {
   int n_nodes = order + 1;
   assert((idx >= 0) && (idx < n_nodes));

   double distance = 2./order;
   vector<double> nodes;
   double node = -1.;
   for(int i = 0; i < n_nodes; i++) {
      nodes.push_back(node);
      node += distance;
   }

   double denominator = 1.;
   double nominator = 1.;
   vector<double> derivative_terms(n_nodes, 1.0);

   for(int i = 0; i < n_nodes; i++) {
      if(i != idx) {
         denominator *= (nodes[idx] - nodes[i]);
         nominator *= (x - nodes[i]);
         for(int j = 0; j < n_nodes; j++) {
            if((j != idx) && (j != i)) {
               derivative_terms[j] *= (x - nodes[i]);
            }
         }
      }
   }

   *value = nominator / denominator;
   *der = 0.0;
   for(int i = 0; i < n_nodes; i++) {
      if(i != idx) {
         *der += derivative_terms[i];
      }
   }
   *der /= denominator;
}


void ReferenceElement::prepare_transformed_values(const Quadrature &q, double element_length,
                                vector<vector<double> > *values, vector<vector<vector<double> > > *gradients) const {
   *values = vector<vector<double> >(Def::num_children);
   *gradients = vector<vector<vector<double> > >(Def::num_children);

   for(int node = 0; node < Def::num_children; node++) {
      int x_id_1D = node % 2;
      int y_id_1D = (node % 4) / 2;
      int z_id_1D = node / 4;

      for(unsigned int q_idx = 0; q_idx < q.np(); q_idx++) {
         double value_x, der_x, value_y, der_y, value_z, der_z;
         shape_fun_1d(x_id_1D, q.coords[q_idx][0], &value_x, &der_x);
         der_x /= element_length;
         shape_fun_1d(y_id_1D, q.coords[q_idx][1], &value_y, &der_y);
         der_y /= element_length;

         if(Def::num_dim == 3) {
            shape_fun_1d(z_id_1D, (q.coords[q_idx])[2], &value_z, &der_z);
            der_z /= element_length;
         }

         double value = value_x * value_y;
         double grad_1 = der_x * value_y;
         double grad_2 = value_x * der_y;

         if(Def::num_dim == 3) {
            value *= value_z;
            grad_1 *= value_z;
            grad_2 *= value_z;
            double grad_3 = value_x * value_y * der_z;
            vector<double> aux1;
            aux1.push_back(grad_1);
            aux1.push_back(grad_2);
            aux1.push_back(grad_3);
            //gradients[node].push_back(vector<double>({grad_1, grad_2, grad_3}));
            (*gradients)[node].push_back(aux1);
         }
         else {
            vector<double> aux1;
            aux1.push_back(grad_1);
            aux1.push_back(grad_2);
            //gradients[node].push_back(vector<double>({grad_1, grad_2}));
            (*gradients)[node].push_back(aux1);
         }
         (*values)[node].push_back(value);
      }
   }
}

