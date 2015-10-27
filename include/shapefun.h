#ifndef SHAPEFUN_H
#define SHAPEFUN_H

#include <vector>

class Quadrature;

class ReferenceElement
{
public:
   ReferenceElement(int num_dim, int order);

   void shape_fun_1d(int idx, double x, double *value, double* der) const;
   void prepare_transformed_values(const Quadrature &q, double element_length,
                                   std::vector<std::vector<double> > *values,
                                   std::vector<std::vector<std::vector<double> > > *gradients) const;
   void find_node_coords(const std::vector<double> &start, double element_len,
                         std::vector<std::vector<double> > *coords) const;

   void print_node_types() const;

public:      
   int num_dim;
   int order;

   int num_nodes;
   int num_corner_nodes;
   int num_edge_nodes;
   int num_face_nodes;
   int num_interior_nodes;

   std::vector<std::vector<double> > node_coords;

   std::vector<int> corner_nodes;
   std::vector<std::vector<int> > edge_nodes;
   std::vector<std::vector<int> > face_nodes;
   std::vector<std::vector<int> > edge_interior_nodes;
   std::vector<std::vector<int> > face_interior_nodes;
   std::vector<int> element_interior_nodes;
};

//void ref_value_1D(int loc_id_1d, double x, double elem_len, double *value, double *der);

#endif // SHAPEFUN_H
