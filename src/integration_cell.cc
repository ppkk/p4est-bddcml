#include "integration_cell.h"

using namespace std;

IntegrationCell::IntegrationCell(double x, double y, double size) {
   position.push_back(x);
   position.push_back(y);
   this->size = size;
}

IntegrationCell::IntegrationCell(double x, double y, double z, double size) : IntegrationCell(x, y, size) {
   position.push_back(z);
}

void IntegrationCell::clear() {
   position.clear();
   size = 0.0;
}

vector<vector<double> > IntegrationCell::corners() const {
   vector<vector<double> > corners_tmp;
   vector<double> corner(n_dimensions(), 0.0);
   int difs[n_dimensions()];
   for(difs[2] = 0; difs[2] < ((n_dimensions() == 3) ? 2 : 1); difs[2]++) {
      for(difs[1] = 0; difs[1] < 2; difs[1]++) {
         for(difs[0] = 0; difs[0] < 2; difs[0]++) {
            // now we know which of the 4/8 points we want, construct its coordinates
            for(int dim = 0; dim < n_dimensions(); dim++) {
               corner[dim] = corner[dim] + difs[dim] * size;
            }
            corners_tmp.push_back(corner);
         }
      }
   }

   return corners_tmp;
}


