#include <algorithm>

#include "shapefun.h"
#include "p4est/my_p4est_interface.h"
#include "quadrature.h"
#include "integration_cell.h"

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

void ReferenceElement::find_nodes_coords(const std::vector<double> &start, double element_len,
                                   std::vector<std::vector<double> > *coords) const {
   coords->clear();
   int num_points_1d = order + 1;
   double node_distance = element_len / order;
   vector<double> node(num_dim, 0.0);
   int difs[num_dim];
   for(difs[2] = 0; difs[2] < ((num_dim == 3) ? num_points_1d : 1); difs[2]++) {
      for(difs[1] = 0; difs[1] < num_points_1d; difs[1]++) {
         for(difs[0] = 0; difs[0] < num_points_1d; difs[0]++) {
            // now we know which of the dim^(order+1) points we want, construct its coordinates
            for(int dim = 0; dim < num_dim; dim++) {
               node[dim] = start[dim] + difs[dim] * node_distance;
            }
            coords->push_back(node);
         }
      }
   }
}

void ReferenceElement::find_nodes_coords(const IntegrationCell &cell, std::vector<std::vector<double> > *coords) const {
   find_nodes_coords(cell.position, cell.size, coords);
}


ReferenceElement::ReferenceElement(int num_dim, int order) : num_dim(num_dim), order(order)
{   
   prepare_node_categories();
   prepare_children_nodes_values();
}

void ReferenceElement::prepare_node_categories()
{
   int num_points_1d = order + 1;
   vector<double> starting_corner(num_dim, -1.0);
   // use for reference element [-1,1]^dim
   find_nodes_coords(starting_corner, 2., &node_coords);

   // order of faces in my difs[] ordering...
   int face_ids[3][2] = {{0,1}, {2,3}, {4,5}};

   // find face nodes
   face_nodes.resize(Def::d()->num_faces);
   int node_idx = 0;

   for(auto difs : Def::d()->cartesian_ids_nodes) {
      for(int dim = 0; dim < num_dim; dim++) {
         if(difs[dim] == 0)
            face_nodes[face_ids[dim][0]].push_back(node_idx);
         if(difs[dim] == num_points_1d - 1)
            face_nodes[face_ids[dim][1]].push_back(node_idx);
      }

      node_idx++;
   }
   assert(node_idx == Def::d()->num_element_nodes);

   // find corner nodes
   corner_nodes.resize(Def::d()->num_corners);
   std::vector<int> intersection;
   for(int corner = 0; corner < Def::d()->num_corners; corner++) {
      intersection.clear();
      for(int node = 0; node < Def::d()->num_element_nodes; node++) {
         bool contained = true;
         for(int face : Def::d()->corner_faces[corner]) {
            if(std::find(face_nodes[face].begin(), face_nodes[face].end(), node) == face_nodes[face].end()) {
               contained = false;
            }
         }
         if(contained) {
            intersection.push_back(node);
         }
      }
      assert((int)intersection.size() == Def::d()->num_corner_nodes);
      corner_nodes[corner] = intersection[0];
   }

   if(num_dim == 2) {
      // find face interior nodes (faces minus corners)
      face_interior_nodes.resize(Def::d()->num_faces);
      for(int face = 0; face < Def::d()->num_faces; face++) {
         for(int node : face_nodes[face]) {
            bool contained_inside = true;
            for(int corner : Def::d()->face_corners[face]) {
               if(corner_nodes[corner] == node) {
                  contained_inside = false;
               }
            }
            if(contained_inside) {
               face_interior_nodes[face].push_back(node);
            }
         }
         assert((int)face_interior_nodes[face].size() == Def::d()->num_face_interior_nodes);
      }
   }
   else if(num_dim == 3) {
      // find edge nodes
      edge_nodes.resize(Def::d()->num_edges);
      std::vector<int> intersection;
      for(int edge = 0; edge < Def::d()->num_edges; edge++) {
         intersection.clear();
         for(int node = 0; node < Def::d()->num_element_nodes; node++) {
            bool contained = true;
            for(int face : Def::d()->edge_faces[edge]) {
               if(std::find(face_nodes[face].begin(), face_nodes[face].end(), node) == face_nodes[face].end()) {
                  contained = false;
               }
            }
            if(contained) {
               intersection.push_back(node);
            }
         }
         assert((int)intersection.size() == Def::d()->num_edge_nodes);
         edge_nodes[edge] = intersection;
      }

      // find face interior nodes (faces minus its edges)
      face_interior_nodes.resize(Def::d()->num_faces);
      for(int face = 0; face < Def::d()->num_faces; face++) {
         for(int node : face_nodes[face]) {
            bool contained_inside = true;
            for(int edge : Def::d()->face_edges[face]) {
               if(find(edge_nodes[edge].begin(), edge_nodes[edge].end(), node) != edge_nodes[edge].end()) {
                  contained_inside = false;
               }
            }
            if(contained_inside) {
               face_interior_nodes[face].push_back(node);
            }
         }
         assert((int)face_interior_nodes[face].size() == Def::d()->num_face_interior_nodes);
      }

      // find edge interior nodes (edge minus corners)
      edge_interior_nodes.resize(Def::d()->num_edges);
      for(int edge = 0; edge < Def::d()->num_edges; edge++) {
         for(int node : edge_nodes[edge]) {
            bool contained_inside = true;
            for(int corner : Def::d()->edge_corners[edge]) {
               if(corner_nodes[corner] == node) {
                  contained_inside = false;
               }
            }
            if(contained_inside) {
               edge_interior_nodes[edge].push_back(node);
            }
         }
         assert((int)edge_interior_nodes[edge].size() == Def::d()->num_edge_interior_nodes);
      }
   }
   else {
      assert(0);
   }

   // find element interior nodes (all minus faces)
   for(int node = 0; node < Def::d()->num_element_nodes; node++) {
      bool contained_inside = true;
      for(int face = 0; face < Def::d()->num_faces; face++) {
         if(find(face_nodes[face].begin(), face_nodes[face].end(), node) != face_nodes[face].end()) {
            contained_inside = false;
         }
      }
      if(contained_inside) {
         element_interior_nodes.push_back(node);
      }
   }
   assert((int)element_interior_nodes.size() == Def::d()->num_element_interior_nodes);
}

void ReferenceElement::prepare_children_nodes_values()
{
   const int num_nodes_1D = order + 1;
   const double child_elem_size = 1.;
   for(auto difs : Def::d()->cartesian_ids_corners) {
      vector<double> child_begining;
      assert(child_begining.empty());
      vector<vector<double> > table(Def::d()->num_element_nodes, vector<double>(Def::d()->num_element_nodes));
      for(int dim = 0; dim < Def::d()->num_dim; dim++) {
         child_begining.push_back(-1 + difs[dim]);
      }
      IntegrationCell cell(child_begining, child_elem_size);
      vector<vector<double> > nodes_coords = cell.nodes_coords(num_nodes_1D);

      for(int child_node = 0; child_node < Def::d()->num_element_nodes; child_node++) {
         for(int parent_node = 0; parent_node < Def::d()->num_element_nodes; parent_node++) {
            table[child_node][parent_node] = shape_value(parent_node, nodes_coords[child_node]);
         }
      }

      children_nodes_parent_basis_values.push_back(table);
   }
}

void ReferenceElement::print_node_types() const
{
   for(int face = 0; face < Def::d()->num_faces; face++) {
      cout << "Face " << face << ": ";
      for (int node : face_nodes[face]) {
         cout << node << ", ";
      }
      cout << endl;
   }

   for(int face = 0; face < Def::d()->num_faces; face++) {
      cout << "Face interior " << face << ": ";
      for (int node : face_interior_nodes[face]) {
         cout << node << ", ";
      }
      cout << endl;
   }

   for(int edge = 0; edge < Def::d()->num_edges; edge++) {
      cout << "Edge " << edge << ": ";
      for (int node : edge_nodes[edge]) {
         cout << node << ", ";
      }
      cout << endl;
   }

   for(int edge = 0; edge < Def::d()->num_edges; edge++) {
      cout << "Edge interior " << edge << ": ";
      for (int node : edge_interior_nodes[edge]) {
         cout << node << ", ";
      }
      cout << endl;
   }

   for(int corner = 0; corner < Def::d()->num_corners; corner++) {
      cout << "Corner " << corner << ": " << corner_nodes[corner] << endl;
   }
   cout << "element interior: ";
   for(int node : element_interior_nodes) {
      cout << node << ", ";
   }
   cout << endl;
   assert(0);
}

bool ReferenceElement::face_contains(int face, int node) const {
   return std::find(face_nodes[face].begin(), face_nodes[face].end(), node) != face_nodes[face].end();
}

bool ReferenceElement::edge_contains(int edge, int node) const {
   return std::find(edge_nodes[edge].begin(), edge_nodes[edge].end(), node) != edge_nodes[edge].end();
}

void ReferenceElement::shape_fun_1d(int idx_1d, double x, double *value, double* der) const {
   int n_nodes = order + 1;
   assert((idx_1d >= 0) && (idx_1d < n_nodes));

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
      if(i != idx_1d) {
         denominator *= (nodes[idx_1d] - nodes[i]);
         nominator *= (x - nodes[i]);
         for(int j = 0; j < n_nodes; j++) {
            if((j != idx_1d) && (j != i)) {
               derivative_terms[j] *= (x - nodes[i]);
            }
         }
      }
   }

   *value = nominator / denominator;
   *der = 0.0;
   for(int i = 0; i < n_nodes; i++) {
      if(i != idx_1d) {
         *der += derivative_terms[i];
      }
   }
   *der /= denominator;
}

double ReferenceElement::shape_fun_1d(int idx_1d, double x) const {
   double value, derivative;
   shape_fun_1d(idx_1d, x, &value, &derivative);
   return value;
}

double ReferenceElement::shape_value(int node_idx, vector<double> coords) const {
   assert((node_idx >= 0) && (node_idx < Def::d()->num_element_nodes));
   assert((int)coords.size() == num_dim);
   double value = 1.;
   for (int dim = 0; dim < num_dim; dim++) {
      value *= shape_fun_1d(node_idx % (order + 1), coords[dim]);
      node_idx /= (order + 1);
   }
   return value;
}

void ReferenceElement::fill_transformed_values(const Quadrature &q, double element_length,
                                vector<vector<double> > *values, vector<vector<vector<double> > > *gradients) const {
   *values = vector<vector<double> >(Def::d()->num_element_nodes);
   *gradients = vector<vector<vector<double> > >(Def::d()->num_element_nodes);

   // todo pass order and do not use Def::d()->order

   for(int node = 0; node < Def::d()->num_element_nodes; node++) {
      int node_tmp = node;
      int ids_1D[Def::d()->num_dim];
      for(int dim = 0; dim < Def::d()->num_dim; dim++) {
         ids_1D[dim] = node_tmp % (Def::d()->order + 1);
         node_tmp /= (Def::d()->order + 1);
      }

      for(unsigned int q_idx = 0; q_idx < q.np(); q_idx++) {
         double values_1D[Def::d()->num_dim], derivatives_1D[Def::d()->num_dim];

         for(int dim = 0; dim < Def::d()->num_dim; dim++) {
            shape_fun_1d(ids_1D[dim], q.coords[q_idx][dim], &values_1D[dim], &derivatives_1D[dim]);
            // partial transformation
            // todo: should it be here?
            derivatives_1D[dim] /= element_length;
         }

         double value = 1.;
         vector<double> grad(Def::d()->num_dim, 1.);

         for(int dim_out = 0; dim_out < Def::d()->num_dim; dim_out++) {
            value *= values_1D[dim_out];
            for(int dim_in = 0; dim_in < Def::d()->num_dim; dim_in++) {
               if(dim_out == dim_in)
                  grad[dim_out] *= derivatives_1D[dim_in];
               else
                  grad[dim_out] *= values_1D[dim_in];
            }
         }
         (*values)[node].push_back(value);
         (*gradients)[node].push_back(grad);
      }
   }
}

