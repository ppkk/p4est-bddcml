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

void IntegrationCell::nodes_coords(int num_nodes_1d, vector<vector<double> > *coords) const {
   assert(num_nodes_1d >= 2);
   coords->clear();
   vector<double> node_coord(n_dimensions(), 0.0);
   double nodes_distance = size / (num_nodes_1d - 1);

   // if I want anything else, I have to extend cartesian_ids in Defs
   assert((num_nodes_1d == 2) || (num_nodes_1d == Def::d()->order + 1));
   const std::vector<std::vector<int> >& cartesian_ids =
         ((num_nodes_1d == 2) ? Def::d()->cartesian_ids_corners : Def::d()->cartesian_ids_nodes);

   for(auto difs : cartesian_ids) {
      for(int dim = 0; dim < n_dimensions(); dim++) {
         node_coord[dim] = position[dim] + difs[dim] * nodes_distance;
      }
      coords->push_back(node_coord);
   }
}

void IntegrationCell::corners_coords(vector<vector<double> > *coords) const {
   nodes_coords(2, coords);
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
