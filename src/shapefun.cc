#include <assert.h>

#include "shapefun.h"
#include "p4est/my_p4est_interface.h"
#include "quadrature.h"

using namespace std;

void ref_value_1D(int loc_id_1d, double x, double elem_len, double *value, double *der) {
   if(loc_id_1d == 0) {
      *value = (1-x)/2.;
      *der = -1/elem_len;
   }
   else if(loc_id_1d == 1) {
      *value = (1+x)/2.;
      *der = 1/elem_len;
   }
   else {
      assert(0);
   }
}

void shape_fun(int order, int idx, double x, double elem_len, double *value, double* der) {
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

   *der /= elem_len; // todo ne takhle
}


void prepare_transformed_values(const Quadrature &q, double element_length,
                                vector<vector<double> > *values, vector<vector<vector<double> > > *gradients) {
   *values = vector<vector<double> >(Def::num_children);
   *gradients = vector<vector<vector<double> > >(Def::num_children);

   for(int node = 0; node < Def::num_children; node++) {
      int x_id_1D = node % 2;
      int y_id_1D = (node % 4) / 2;
      int z_id_1D = node / 4;

      for(unsigned int q_idx = 0; q_idx < q.np(); q_idx++) {
         double value_x, der_x, value_y, der_y, value_z, der_z;
//         ref_value_1D(x_id_1D, q.coords[q_idx][0], element_length, &value_x, &der_x);
//         ref_value_1D(y_id_1D, q.coords[q_idx][1], element_length, &value_y, &der_y);
         shape_fun(1, x_id_1D, q.coords[q_idx][0], element_length, &value_x, &der_x);
         shape_fun(1, y_id_1D, q.coords[q_idx][1], element_length, &value_y, &der_y);


         if(Def::num_dim == 3) {
            //ref_value_1D(z_id_1D, (q.coords[q_idx])[2], element_length, &value_z, &der_z);
            shape_fun(1, z_id_1D, (q.coords[q_idx])[2], element_length, &value_z, &der_z);
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

