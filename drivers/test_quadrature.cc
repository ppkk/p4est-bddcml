#include <iostream>
#include "math.h"

#include "element.h"
#include "quadrature.h"
#include "level_set.h"

using namespace std;

const int dimension = 2;
const int order = 15;
const int elems_in_dir = 40;

double level_set_fn_circle(vector<double> point)
{
   double rr = 0.0;
   for(double x : point)
      rr += x*x;

   return 1 - rr;
}

double integrate_element(const Quadrature &quadrature, const LevelSet &level_set, const Element &element)
{
   double result = 0.0;
   Quadrature quad_trans(dimension);
   quadrature.transform_to_physical(element, &quad_trans);
   for(unsigned qi = 0; qi < quad_trans.np(); qi++)
   {
      if(level_set.apply(quad_trans.coords[qi]) == LevelSetValue::Inside)
      {
         result += quad_trans.weights[qi];

//         double rr = quad_trans.coords[qi][0]*quad_trans.coords[qi][0] + quad_trans.coords[qi][1]*quad_trans.coords[qi][1];
//         cout << "including point (" << quad_trans.coords[qi][0] << ", " << quad_trans.coords[qi][1] << "), rr is " << rr << ", contrib " << quad_trans.weights[qi] << endl;
      }
   }

   return result;
}

double integrate(const Quadrature &quadrature, const LevelSet &level_set, const vector<Element> &elements)
{
   double result = 0.0;
   for(const Element &element : elements)
   {
      result += integrate_element(quadrature, level_set, element);
   }

   return result;
}


void generate_elements(int num_in_dir, vector<Element> *elements)
{
   double size = 2./num_in_dir;
   vector<double> coords;
   int offsets[dimension];
#ifdef P4_TO_P8
   for(offsets[2] = 0; offsets[2] < num_in_dir; offsets[2]++)
   {
#endif
      for(offsets[1] = 0; offsets[1] < num_in_dir; offsets[1]++)
      {
         for(offsets[0] = 0; offsets[0] < num_in_dir; offsets[0]++)
         {
            coords.clear();
            for(int i = 0; i < dimension; i++)
               coords.push_back(-1. + offsets[i] * size);
            elements->push_back(Element(coords, size));
         }
      }
#ifdef P4_TO_P8
   }
#endif
}

int main()
{
   GaussQuadrature gauss_quad(dimension, order);
   EquidistantQuadrature equidist_quad(dimension, order);
//   gauss_quad.print();
//   equidist_quad.print();

   LevelSet level_set(level_set_fn_circle);

   vector<Element> elements;

   generate_elements(elems_in_dir, &elements);

   double gauss_result = integrate(gauss_quad, level_set, elements);
   double equidist_result = integrate(equidist_quad, level_set, elements);

   cout << "gauss error: " << fabs(gauss_result - M_PI) << ", equidist: " << fabs(equidist_result - M_PI) << endl;
}
