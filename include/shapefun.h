#ifndef SHAPEFUN_H
#define SHAPEFUN_H

#include "definitions.h"

class Quadrature;
class IntegrationCell;

class ReferenceElement
{
public:
   ReferenceElement(int num_dim, int order);

   // 1D lagrange basis function, 0 <= idx_1d < order + 1
   void shape_fun_1d(int idx_1d, double x, double *value, double* der) const;
   double shape_fun_1d(int idx_1d, double x) const;

   // value of 2D/3D basis function corresponding to node_idx on reference domain [-1,1]^dim
   // 0 <= node_idx < (order+1)^dim
   double shape_value(int node_idx, std::vector<double> coords) const;

   void calc_transformed_values(const Quadrature &q, double element_length,
                                   std::vector<std::vector<double> > *values,
                                   std::vector<std::vector<std::vector<double> > > *gradients) const;
   void find_nodes_coords(const std::vector<double> &start, double element_len,
                         std::vector<std::vector<double> > *coords) const;
   void find_nodes_coords(const IntegrationCell &cell, std::vector<std::vector<double> > *coords) const;

   void print_node_types() const;

   bool face_contains(int face, int node) const;
   bool edge_contains(int face, int node) const;

private:
   void prepare_node_categories();
   void prepare_children_nodes_values();

public:      
   int num_dim;
   int order;

   std::vector<std::vector<double> > node_coords;

   // nodes in those arrays are sorted by construction
   // this is important : allowes to use std::find
   std::vector<int> corner_nodes;
   std::vector<std::vector<int> > edge_nodes;
   std::vector<std::vector<int> > face_nodes;
   std::vector<std::vector<int> > edge_interior_nodes;
   std::vector<std::vector<int> > face_interior_nodes;
   std::vector<int> element_interior_nodes;

   // [child_type][child_node_id][parent_node(shapefn)_id]
   // child_type is 0..3 or 0..7 for 2D or 3D
   std::vector<std::vector<std::vector<double> > > children_nodes_parent_basis_values;
};

//void ref_value_1D(int loc_id_1d, double x, double elem_len, double *value, double *der);

#endif // SHAPEFUN_H
