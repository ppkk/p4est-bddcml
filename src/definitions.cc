#include <assert.h>

#include "definitions.h"

void Def::init(int num_dim, int order)
{
   Def::num_dim = num_dim;
   Def::order = order;
   if(num_dim == 2) {
      num_children = 4;
      num_corners = 4;
      num_edges = 0;
      num_faces = 4;

      num_loc_dofs = (order + 1) * (order + 1);
   }
   else if(num_dim == 3) {
      num_children = 8;
      num_corners = 8;
      num_edges = 12;
      num_faces = 6;

      num_loc_dofs = (order + 1) * (order + 1) * (order + 1);
   }
   else {
      assert(0);
   }
}

int Def::num_dim;
int Def::order;
int Def::num_children;
int Def::num_corners;
int Def::num_edges;
int Def::num_faces;
int Def::num_loc_dofs;
