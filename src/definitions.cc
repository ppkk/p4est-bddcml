#include "definitions.h"
#include "p4est/my_p4est_interface.h"

using namespace std;

void Def::init(int num_dim, int order, PhysicsType physicsType, const P4estClass *p4est)
{
   singleton = new Def;
   singleton->num_dim = num_dim;
   singleton->order = order;

   singleton->prepare_cartesian_ids(2, &singleton->cartesian_ids_corners);
   singleton->prepare_cartesian_ids(order + 1, &singleton->cartesian_ids_nodes);
   singleton->prepare_cartesian_ids(order, &singleton->cartesian_ids_plot_subelements);

   if(num_dim == 2) {
      singleton->num_children = num_children_2D;
      singleton->num_corners = num_corners_2D;
      singleton->num_edges = num_edges_2D;
      singleton->num_faces = num_faces_2D;

      singleton->num_face_corners = num_face_corners_2D;
      singleton->num_corner_faces = num_corner_faces_2D;

      singleton->num_edge_corners = num_edge_corners_2D;
      singleton->num_corner_edges = num_corner_edges_2D;

      singleton->num_face_edges = num_face_edges_2D;
      singleton->num_edge_faces = num_edge_faces_2D;

      singleton->num_element_nodes = (order + 1) * (order + 1);
      singleton->num_face_nodes = order + 1;
      singleton->num_edge_nodes = 0;
      singleton->num_corner_nodes = 1;
      singleton->num_element_interior_nodes = (order - 1) * (order - 1);
      singleton->num_face_interior_nodes = order - 1;
      singleton->num_edge_interior_nodes = 0;
   }
   else if(num_dim == 3) {
      singleton->num_children = num_children_3D;
      singleton->num_corners = num_corners_3D;
      singleton->num_edges = num_edges_3D;
      singleton->num_faces = num_faces_3D;

      singleton->num_face_corners = num_face_corners_3D;
      singleton->num_corner_faces = num_corner_faces_3D;

      singleton->num_edge_corners = num_edge_corners_3D;
      singleton->num_corner_edges = num_corner_edges_3D;

      singleton->num_face_edges = num_face_edges_3D;
      singleton->num_edge_faces = num_edge_faces_3D;

      singleton->num_element_nodes = (order + 1) * (order + 1) * (order + 1);
      singleton->num_face_nodes = (order + 1) * (order + 1);
      singleton->num_edge_nodes = order + 1;
      singleton->num_corner_nodes = 1;
      singleton->num_element_interior_nodes = (order - 1) * (order - 1) * (order - 1);
      singleton->num_face_interior_nodes = (order - 1) * (order - 1);
      singleton->num_edge_interior_nodes = order - 1;
   }
   else {
      assert(0);
   }

   if(physicsType == PhysicsType::LAPLACE)
      singleton->num_components = 1;
   else if (physicsType == PhysicsType::ELASTICITY)
      singleton->num_components = num_dim;
   else
      assert(0);

   // copy important connectivity information froma p4est
   // has to be done in my_p4est_implementation to distinguish 2D and 3D
   if(p4est != nullptr) {
      p4est->init_definitions(singleton);
   }
}

void Def::prepare_cartesian_ids(int num_points_1d, std::vector<std::vector<int> > *ids) const {
   ids->clear();
   int difs[3];
   for(difs[2] = 0; difs[2] < ((num_dim == 3) ? num_points_1d : 1); difs[2]++) {
      for(difs[1] = 0; difs[1] < num_points_1d; difs[1]++) {
         for(difs[0] = 0; difs[0] < num_points_1d; difs[0]++) {
            vector<int> cart;
            cart.clear();
            for(int dim = 0; dim < num_dim; dim++) {
               cart.push_back(difs[dim]);
            }
            ids->push_back(cart);
         }
      }
   }
}

Def* Def::singleton;

/**********************************************************************************************************/
/**********************************************************************************************************/

ProblemDimensions::ProblemDimensions(int mesh_dim, PhysicsType physicsType) {
   n_problem_dims = mesh_dim;
   n_mesh_dims = mesh_dim;
   n_dofs = 0;
   n_elems = 0;
   n_nodes = 0;

   if(mesh_dim == 2)
      n_elem_nodes = 4;
   else if(mesh_dim == 3)
      n_elem_nodes = 8;
   else
      assert(0);

   if(physicsType == PhysicsType::LAPLACE)
      n_node_dofs = 1;
   else if(physicsType == PhysicsType::ELASTICITY)
      n_node_dofs = mesh_dim;
   else
      assert(0);
}

/**********************************************************************************************************/
/**********************************************************************************************************/

// todo: move somewhere else
int print_rank = 0;
int mpi_rank;
int mpi_size;
MPI_Comm mpicomm = MPI_COMM_WORLD;

