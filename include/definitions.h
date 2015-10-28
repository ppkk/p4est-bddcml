#ifndef DEFINITIONS_H
#define DEFINITIONS_H

#include <assert.h>
#include <vector>
#include <iostream>

// type used for floating point
typedef double real;

extern int mpi_rank;
extern int mpi_size;
extern int print_rank;

#define PPP if(mpi_rank == print_rank)


enum class MatrixType{
   GENERAL,      // general (full storage)
   SPD,          //symmetric positive definite (only triangle stored)
   SYM_GENERAL  // symmetric general (only triangle stored)
};

enum class PhysicsType{
   LAPLACE = 0,
   ELASTICITY = 1
};

struct Parameters
{
   real mu;
   real lambda;

   Parameters(real young_mod, real poisson_num)
   {
      this->lambda = (young_mod * poisson_num) / ((1 + poisson_num) * (1 - 2 * poisson_num));
      this->mu = young_mod / (2 * (1 + poisson_num));
   }
};


class Def {
public:
   static void init(int num_dim, int order, PhysicsType physicsType);
   static int num_dim;
   static int num_children;
   static int num_corners;
   static int num_edges;
   static int num_faces;

   static int num_face_corners;
   static int num_corner_faces;

   static int num_edge_corners;
   static int num_corner_edges;

   static int num_edge_faces;
   static int num_face_edges;

   static int num_components;

   static int order;   // this actually should not be here...
   static int num_element_nodes;
   static int num_face_nodes;
   static int num_edge_nodes;
   static int num_corner_nodes;
   static int num_element_interior_nodes;
   static int num_face_interior_nodes;
   static int num_edge_interior_nodes;

   static std::vector<std::vector<int> > edge_corners;
   static std::vector<std::vector<int> > face_corners;

   static std::vector<std::vector<int> > corner_edges;
   static std::vector<std::vector<int> > corner_faces;

   static std::vector<std::vector<int> > face_edges;
   static std::vector<std::vector<int> > edge_faces;

private:

   static const int num_children_2D = 4;
   static const int num_corners_2D = 4;
   static const int num_edges_2D = 0;
   static const int num_faces_2D = 4;

   static const int num_face_corners_2D = 2;
   static const int num_corner_faces_2D = 2;

   static const int num_edge_corners_2D = -1;
   static const int num_corner_edges_2D = -1;

   static const int num_edge_faces_2D = -1;
   static const int num_face_edges_2D = -1;


   static const int num_children_3D = 8;
   static const int num_corners_3D = 8;
   static const int num_edges_3D = 12;
   static const int num_faces_3D = 6;

   static const int num_face_corners_3D = 4;
   static const int num_corner_faces_3D = 3;

   static const int num_edge_corners_3D = 2;
   static const int num_corner_edges_3D = 3;

   static const int num_edge_faces_3D = 2;
   static const int num_face_edges_3D = 4;

   friend class P4estClass2D;
   friend class P4estClass3D;
};


#endif // DEFINITIONS_H
