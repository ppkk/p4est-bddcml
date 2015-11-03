#include "level_set.h"
#include "integration_cell.h"

using namespace std;

LevelSetValue LevelSet::apply(const std::vector<double> point) const {
   if(fn(point) > 0.)
      return LevelSetValue::Inside;
   if(fn(point) == 0.)
      return LevelSetValue::Border;
   if(fn(point) < 0.)
      return LevelSetValue::Outside;
   else
      assert(0);
}

LevelSetValue LevelSet::apply(const IntegrationCell* cell) const {
   LevelSetValue result = LevelSetValue::Border;

   vector<vector<double> > corners_coords;
   cell->corners_coords(&corners_coords);
   for (auto node : corners_coords) {
      LevelSetValue node_result = apply(node);
      if(result == LevelSetValue::Border) {
         // first node or the previous were also on the Border
         result = node_result;
      }
      // one of the previous has allready been set to Interior/Exterior. Now if this node does not correspond
      else if(result != node_result) {
         // and if it is not on the Border (which we tolerate - do not count as intersected element)
         if(node_result != LevelSetValue::Border) {
            // this element is than intersected
            return LevelSetValue::Border;
         }
      }
   }

   // all nodes were of the same type (or on the border)
   return result;
}
