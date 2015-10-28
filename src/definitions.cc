#include "definitions.h"
#include "p4est/my_p4est_interface.h"

using namespace std;

void Def::init(int num_dim, int order, PhysicsType physicsType)
{
   Def::num_dim = num_dim;
   Def::order = order;
   if(num_dim == 2) {
      num_children = num_children_2D;
      num_corners = num_corners_2D;
      num_edges = num_edges_2D;
      num_faces = num_faces_2D;

      num_face_corners = num_face_corners_2D;
      num_corner_faces = num_corner_faces_2D;

      num_edge_corners = num_edge_corners_2D;
      num_corner_edges = num_corner_edges_2D;

      num_face_edges = num_face_edges_2D;
      num_edge_faces = num_edge_faces_2D;

      num_element_nodes = (order + 1) * (order + 1);
      num_face_nodes = order + 1;
      num_edge_nodes = 0;
      num_corner_nodes = 1;
      num_element_interior_nodes = (order - 1) * (order - 1);
      num_face_interior_nodes = order - 1;
      num_edge_interior_nodes = 0;
   }
   else if(num_dim == 3) {
      num_children = num_children_3D;
      num_corners = num_corners_3D;
      num_edges = num_edges_3D;
      num_faces = num_faces_3D;

      num_face_corners = num_face_corners_3D;
      num_corner_faces = num_corner_faces_3D;

      num_edge_corners = num_edge_corners_3D;
      num_corner_edges = num_corner_edges_3D;

      num_face_edges = num_face_edges_3D;
      num_edge_faces = num_edge_faces_3D;

      num_element_nodes = (order + 1) * (order + 1) * (order + 1);
      num_face_nodes = (order + 1) * (order + 1);
      num_edge_nodes = order + 1;
      num_corner_nodes = 1;
      num_element_interior_nodes = (order - 1) * (order - 1) * (order - 1);
      num_face_interior_nodes = (order - 1) * (order - 1);
      num_edge_interior_nodes = order - 1;
   }
   else {
      assert(0);
   }

   if(physicsType == PhysicsType::LAPLACE)
      num_components = 1;
   else if (physicsType == PhysicsType::ELASTICITY)
      num_components = num_dim;
   else
      assert(0);

   // copy important connectivity information froma p4est
   // has to be done in my_p4est_implementation to distinguish 2D and 3D
   P4estClass::init_definitions();
}

int Def::num_dim;
int Def::order;
int Def::num_components;
int Def::num_children;
int Def::num_corners;
int Def::num_edges;
int Def::num_faces;

int Def::num_element_nodes;
int Def::num_face_nodes;
int Def::num_edge_nodes;
int Def::num_corner_nodes;
int Def::num_element_interior_nodes;
int Def::num_face_interior_nodes;
int Def::num_edge_interior_nodes;

int Def::num_face_corners;
int Def::num_corner_faces;

int Def::num_edge_corners;
int Def::num_corner_edges;

int Def::num_edge_faces;
int Def::num_face_edges;

vector<vector<int> > Def::edge_corners;
vector<vector<int> > Def::face_corners;

vector<vector<int> > Def::corner_edges;
vector<vector<int> > Def::corner_faces;

vector<vector<int> > Def::face_edges;
vector<vector<int> > Def::edge_faces;
