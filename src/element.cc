#include "element.h"

using namespace std;


vector<vector<double> > Element::nodes() const
{
   vector<vector<double> > nodes_tmp;
   vector<double> node;
   int difs[n_dimensions()];
#ifdef P4_TO_P8
   for(difs[2] = 0; difs[2] < 2; difs[2]++)
   {
#endif
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
#ifdef P4_TO_P8
   }
#endif

   return nodes_tmp;
}

