#include "element.h"

using namespace std;

Element::Element(double x, double y, double size)
{
   position.push_back(x);
   position.push_back(y);
   this->size = size;
}

Element::Element(double x, double y, double z, double size) : Element(x, y, size)
{
   position.push_back(z);
}

void Element::clear()
{
   position.clear();
   size = 0.0;
}

vector<vector<double> > Element::nodes() const
{
   vector<vector<double> > nodes_tmp;
   vector<double> node;
   int difs[n_dimensions()];
   for(difs[2] = 0; difs[2] < (n_dimensions() == 3) ? 2 : 1; difs[2]++)
   {
      for(difs[1] = 0; difs[1] < 2; difs[1]++)
      {
         for(difs[0] = 0; difs[0] < 2; difs[0]++)
         {
            // now we know which of the 4/8 points we want, construct its coordinates
            for(int dim = 0; dim < n_dimensions(); dim++)
            {
               node.push_back(node[dim] + difs[dim] * size);
            }
            nodes_tmp.push_back(node);
         }
      }
   }

   return nodes_tmp;
}

