#include <assert.h>

#include "definitions.h"
#include "p4est/my_p4est_interface.h"

void Def::init(int num_dim, int order)
{
   Def::num_dim = num_dim;
   Def::order = order;
   if(num_dim == 2) {
      num_children = 4;
      num_corners = 4;
      num_edges = 0;
      num_faces = 4;

      num_face_corners = 2;
      num_corner_faces = 2;

      num_edge_corners = -1;
      num_corner_edges = -1;

      num_element_nodes = (order + 1) * (order + 1);
   }
   else if(num_dim == 3) {
      num_children = 8;
      num_corners = 8;
      num_edges = 12;
      num_faces = 6;

      num_face_corners = 4;
      num_corner_faces = 3;

      num_edge_corners = 2;
      num_corner_edges = 3;

      num_element_nodes = (order + 1) * (order + 1) * (order + 1);
   }
   else {
      assert(0);
   }

   // copy important connectivity information froma p4est
   // has to be done in my_p4est_implementation to distinguish 2D and 3D
   P4estClass::init_definitions();
}

int Def::num_dim;
int Def::order;
int Def::num_children;
int Def::num_corners;
int Def::num_edges;
int Def::num_faces;
int Def::num_element_nodes;

int Def::num_face_corners;
int Def::num_corner_faces;

int Def::num_edge_corners;
int Def::num_corner_edges;

std::vector<std::vector<int> > Def::edge_corners;
std::vector<std::vector<int> > Def::face_corners;
