#include "element.h"
#include "quadrature.h"
#include "level_set.h"

using namespace std;

const int dimension = 2;
const int order = 3;

double level_set_fn_circle(vector<double> point)
{
   double rr = 0.0;
   for(double x : point)
      rr += x*x;

   return 1 - rr;
}

double integrate(const Quadrature &quadrature, const LevelSet &level_set, const Element &element)
{
   Quadrature quad_trans(dimension);
   quadrature.transform_to_physical(element, &quad_trans);
}


void generate_elements(int num_in_dir, vector<Element>& elements)
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
               coords.push_back(offsets[i] * size);
            elements.push_back(Element(coords, size));
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

   LevelSet level_set(level_set_fn_circle);
}
