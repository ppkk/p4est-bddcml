#include <math.h>
#include <assert.h>
#include <iostream>

#include "quadrature.h"
#include "element.h"

using namespace std;

void Quadrature::tensor_product(const Quadrature &quad1, const Quadrature &quad2)
{
   for(unsigned int i1 = 0; i1 < quad1.weights.size(); i1++)
   {
      for(unsigned int i2 = 0; i2 < quad2.weights.size(); i2++)
      {
         weights.push_back(quad1.weights[i1] * quad2.weights[i2]);
         vector<double>product_coords;
         product_coords.insert(product_coords.end(), quad1.coords[i1].begin(), quad1.coords[i1].end());
         product_coords.insert(product_coords.end(), quad2.coords[i2].begin(), quad2.coords[i2].end());
         coords.push_back(product_coords);
      }
   }
}

void Quadrature::print()
{
   assert(weights.size() == coords.size());
   double sum = 0.0;
   for(unsigned int i = 0; i < weights.size(); i++)
   {
      cout << "(";
      for(unsigned int j = 0; j < coords[i].size(); j++)
      {
         cout << coords[i][j] << ", ";
      }
      cout << "), " << weights[i] << endl;
      sum += weights[i];
   }
   cout << "sum of weights " << sum << endl;
}

void Quadrature::transform_to_physical(const Element &element, Quadrature *transformed) const
{
   transformed->clear();
   transformed->weights.resize(weights.size(), 0.0);
   transformed->coords.resize(weights.size(), vector<double>(dimension, 0.0));
   double scale_coeff = pow(element.size * 0.5, dimension);
   for(unsigned qi = 0; qi < weights.size(); qi++)
   {
      transformed->weights[qi] = weights[qi] * scale_coeff;
      for(int dim_idx = 0; dim_idx < dimension; dim_idx++)
      {
         transformed->coords[qi][dim_idx] = element.position[dim_idx] + 0.5 * (coords[qi][dim_idx] - 1) * element.size;
      }
   }
}

void Quadrature::clear()
{
   weights.clear();
   coords.clear();
}

GaussQuadrature::GaussQuadrature(int dimension, int order)
   : Quadrature(dimension)
{
   if(dimension == 1)
   {
      GaussTables1D gauss_quad;

      for(int i = 0; i < gauss_quad.get_num_points(order); i++)
      {
         coords.push_back(vector<double>(1, gauss_quad.get_points(order)[i][0]));
         weights.push_back(gauss_quad.get_points(order)[i][1]);
      }
   }
   else if(dimension == 2)
   {
      GaussQuadrature q1(1, order);
      tensor_product(q1, q1);
   }
   else if(dimension == 3)
   {
      GaussQuadrature q1(1, order), q2(2, order);
      tensor_product(q1, q2);
   }
   else
      assert(0);
}

EquidistantQuadrature::EquidistantQuadrature(int dimension, int order)
   : Quadrature(dimension)
{
   if(dimension == 1)
   {
      // equidistant points. The number of points is the same as for Gauss quadrature
      GaussTables1D gauss_quad;

      int num_points =  gauss_quad.get_num_points(order);
      double weight = 2. / num_points;
      double point_distance = 2. / (num_points + 1);
      double x = point_distance;
      for(int i = 0; i < num_points; i++)
      {
         coords.push_back(vector<double>(1, x));
         x += point_distance;
         weights.push_back(weight);
      }

   }
   else if(dimension == 2)
   {
      EquidistantQuadrature q1(1, order);
      tensor_product(q1, q1);
   }
   else if(dimension == 3)
   {
      EquidistantQuadrature q1(1, order), q2(2, order);
      tensor_product(q1, q2);
   }
   else
      assert(0);
}
