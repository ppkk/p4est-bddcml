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
   child_position = -1;
}

vector<vector<double> > IntegrationCell::nodes_coords(int num_nodes_1d) const {
   assert(num_nodes_1d >= 2);
   vector<vector<double> > corners_tmp;
   vector<double> corner(n_dimensions(), 0.0);
   double nodes_distance = size / (num_nodes_1d - 1);
   int difs[n_dimensions()];
   for(difs[2] = 0; difs[2] < ((n_dimensions() == 3) ? num_nodes_1d : 1); difs[2]++) {
      for(difs[1] = 0; difs[1] < num_nodes_1d; difs[1]++) {
         for(difs[0] = 0; difs[0] < num_nodes_1d; difs[0]++) {
            // now we know which of the 4/8 points we want, construct its coordinates
            for(int dim = 0; dim < n_dimensions(); dim++) {
               corner[dim] = corner[dim] + difs[dim] * nodes_distance;
            }
            corners_tmp.push_back(corner);
         }
      }
   }

   return corners_tmp;
}

vector<vector<double> > IntegrationCell::corners_coords() const {
   return nodes_coords(2);
}

void IntegrationCell::fill_parent_cell(IntegrationCell *parent) const {
   parent->clear();
   int child_position_tmp = child_position;

   for(int dim = 0; dim < n_dimensions(); dim++) {
      if(child_position_tmp % 2 == 0) {
         parent->position.push_back(position[dim]);
      }
      else {
         parent->position.push_back(position[dim] - size);
      }
      child_position_tmp /= 2;
   }
   parent->size = 2 * size;
}
