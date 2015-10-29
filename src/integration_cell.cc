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
   vector<vector<double> > nodes_coords_ret;
   vector<double> node_coord(n_dimensions(), 0.0);
   double nodes_distance = size / (num_nodes_1d - 1);
   for(auto difs : Def::d()->cartesian_ids_nodes) {
      for(int dim = 0; dim < n_dimensions(); dim++) {
         node_coord[dim] = node_coord[dim] + difs[dim] * nodes_distance;
      }
      nodes_coords_ret.push_back(node_coord);
   }
   return nodes_coords_ret;
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
