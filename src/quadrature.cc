#include <math.h>
#include <assert.h>
#include <iostream>

#include "quadrature.h"

using namespace std;

Quadrature::Quadrature(int dimension, int order, double element_length)
{
   double scale = element_length / 2.;
   if(dimension == 1)
   {
      GaussQuad1D gauss_quad;

      for(int i = 0; i < gauss_quad.get_num_points(order); i++)
      {
         coords.push_back(vector<double>(1, gauss_quad.get_points(order)[i][0]));
         weights.push_back(gauss_quad.get_points(order)[i][1] * scale);
      }

      //weights = {5./9. * scale, 8./9. * scale, 5./9. * scale};
      //coords = {vector<double>({-sqrt(3./5.)}), vector<double>({0.}), vector<double>({sqrt(3./5.)})};
   }
   else if(dimension == 2)
   {
      Quadrature q1(1, order, element_length);
      product(q1, q1);
   }
   else if(dimension == 3)
   {
      Quadrature q1(1, order, element_length), q2(2, order, element_length);
      product(q1, q2);
   }
   else
      assert(0);
}

void Quadrature::product(const Quadrature &quad1, const Quadrature &quad2)
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

