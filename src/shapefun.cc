#include <assert.h>
#include <iostream>
#include <algorithm>

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
   int num_points_1d = order + 1;
   vector<double> starting_corner(num_dim, -1.0);
   // use for reference element [-1,1]^dim
   find_node_coords(starting_corner, 2., &node_coords);

   // order of faces in my difs[] ordering...
   int face_ids[3][2] = {{0,1}, {2,3}, {4,5}};

   // find face nodes
   face_nodes.resize(Def::num_faces);
   int difs[num_dim];
   int node_idx = 0;
   for(difs[2] = 0; difs[2] < ((num_dim == 3) ? num_points_1d : 1); difs[2]++) {
      for(difs[1] = 0; difs[1] < num_points_1d; difs[1]++) {
         for(difs[0] = 0; difs[0] < num_points_1d; difs[0]++) {
            for(int dim = 0; dim < num_dim; dim++) {
               if(difs[dim] == 0)
                  face_nodes[face_ids[dim][0]].push_back(node_idx);
               if(difs[dim] == num_points_1d - 1)
                  face_nodes[face_ids[dim][1]].push_back(node_idx);
            }

            node_idx++;
         }
      }
   }
   assert(node_idx == Def::num_element_nodes);

   // find corner nodes
   corner_nodes.resize(Def::num_corners);
   std::vector<int> intersection;
   for(int corner = 0; corner < Def::num_corners; corner++) {
      intersection.clear();
      for(int node = 0; node < Def::num_element_nodes; node++) {
         bool contained = true;
         for(int face : Def::corner_faces[corner]) {
            if(std::find(face_nodes[face].begin(), face_nodes[face].end(), node) == face_nodes[face].end()) {
               contained = false;
            }
         }
         if(contained) {
            intersection.push_back(node);
         }
      }
      assert((int)intersection.size() == Def::num_corner_nodes);
      corner_nodes[corner] = intersection[0];
   }

   if(num_dim == 2) {
      // find face interior nodes (faces minus corners)
      face_interior_nodes.resize(Def::num_faces);
      for(int face = 0; face < Def::num_faces; face++) {
         for(int node : face_nodes[face]) {
            bool contained_inside = true;
            for(int corner : Def::face_corners[face]) {
               if(corner_nodes[corner] == node) {
                  contained_inside = false;
               }
            }
            if(contained_inside) {
               face_interior_nodes[face].push_back(node);
            }
         }
         assert((int)face_interior_nodes[face].size() == Def::num_face_interior_nodes);
      }
   }
   else if(num_dim == 3) {
      // find edge nodes
      edge_nodes.resize(Def::num_edges);
      std::vector<int> intersection;
      for(int edge = 0; edge < Def::num_edges; edge++) {
         intersection.clear();
         for(int node = 0; node < Def::num_element_nodes; node++) {
            bool contained = true;
            for(int face : Def::edge_faces[edge]) {
               if(std::find(face_nodes[face].begin(), face_nodes[face].end(), node) == face_nodes[face].end()) {
                  contained = false;
               }
            }
            if(contained) {
               intersection.push_back(node);
            }
         }
         assert((int)intersection.size() == Def::num_edge_nodes);
         edge_nodes[edge] = intersection;
      }

      // find face interior nodes (faces minus its edges)
      face_interior_nodes.resize(Def::num_faces);
      for(int face = 0; face < Def::num_faces; face++) {
         for(int node : face_nodes[face]) {
            bool contained_inside = true;
            for(int edge : Def::face_edges[face]) {
               if(find(edge_nodes[edge].begin(), edge_nodes[edge].end(), node) != edge_nodes[edge].end()) {
                  contained_inside = false;
               }
            }
            if(contained_inside) {
               face_interior_nodes[face].push_back(node);
            }
         }
         assert((int)face_interior_nodes[face].size() == Def::num_face_interior_nodes);
      }

      // find edge interior nodes (edge minus corners)
      edge_interior_nodes.resize(Def::num_edges);
      for(int edge = 0; edge < Def::num_edges; edge++) {
         for(int node : edge_nodes[edge]) {
            bool contained_inside = true;
            for(int corner : Def::edge_corners[edge]) {
               if(corner_nodes[corner] == node) {
                  contained_inside = false;
               }
            }
            if(contained_inside) {
               edge_interior_nodes[edge].push_back(node);
            }
         }
         assert((int)edge_interior_nodes[edge].size() == Def::num_edge_interior_nodes);
      }
   }
   else {
      assert(0);
   }
}

void ReferenceElement::print_node_types() const
{
   for(int face = 0; face < Def::num_faces; face++) {
      cout << "Face " << face << ": ";
      for (int node : face_nodes[face]) {
         cout << node << ", ";
      }
      cout << endl;
   }

   for(int face = 0; face < Def::num_faces; face++) {
      cout << "Face interior " << face << ": ";
      for (int node : face_interior_nodes[face]) {
         cout << node << ", ";
      }
      cout << endl;
   }

   for(int edge = 0; edge < Def::num_edges; edge++) {
      cout << "Edge " << edge << ": ";
      for (int node : edge_nodes[edge]) {
         cout << node << ", ";
      }
      cout << endl;
   }

   for(int edge = 0; edge < Def::num_edges; edge++) {
      cout << "Edge interior " << edge << ": ";
      for (int node : edge_interior_nodes[edge]) {
         cout << node << ", ";
      }
      cout << endl;
   }

   for(int corner = 0; corner < Def::num_corners; corner++) {
      cout << "Corner " << corner << ": " << corner_nodes[corner] << endl;
   }
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

