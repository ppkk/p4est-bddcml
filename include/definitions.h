#ifndef DEFINITIONS_H
#define DEFINITIONS_H

// todo: remove mpi from here
#include <mpi.h>

#include <assert.h>
#include <vector>
#include <iostream>

class P4estClass;

// type used for floating point
typedef double real;

extern int mpi_rank;
extern int mpi_size;
extern int print_rank;
extern MPI_Comm mpicomm;

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

   Parameters() {
      lambda = 0.0;
      mu = 0.0;
   }

   Parameters(real young_mod, real poisson_num){
      this->lambda = (young_mod * poisson_num) / ((1 + poisson_num) * (1 - 2 * poisson_num));
      this->mu = young_mod / (2 * (1 + poisson_num));
   }
};

class Def {
public:
   static void init(int num_dim, int order, PhysicsType physicsType, const P4estClass &p4est);
   static const Def *d() {return singleton; }

   int num_dim;
   int num_children;
   int num_corners;
   int num_edges;
   int num_faces;

   int num_face_corners;
   int num_corner_faces;

   int num_edge_corners;
   int num_corner_edges;

   int num_edge_faces;
   int num_face_edges;

   int num_components;

   int order;   // this actually should not be here...
   int num_element_nodes;
   int num_face_nodes;
   int num_edge_nodes;
   int num_corner_nodes;
   int num_element_interior_nodes;
   int num_face_interior_nodes;
   int num_edge_interior_nodes;

   std::vector<std::vector<int> > edge_corners;
   std::vector<std::vector<int> > face_corners;

   std::vector<std::vector<int> > corner_edges;
   std::vector<std::vector<int> > corner_faces;

   std::vector<std::vector<int> > face_edges;
   std::vector<std::vector<int> > edge_faces;

   // when going through corners of square/cube in order first x, than y, than z (2^dim)
   // this gives the vector [[0,0], [1,0], [0,1], .. ] or its 3D version
   std::vector<std::vector<int> > cartesian_ids_corners;

   // like cartesian_ids_corners, but for nodes (order^dim)
   // gives [[0,0], [1,0], ... [order,0], [0,1], ... ] or its 3D version
   std::vector<std::vector<int> > cartesian_ids_nodes;

   // the same, but ((order-1)^dim)
   std::vector<std::vector<int> > cartesian_ids_plot_subelements;

private:
   static Def *singleton;

   void prepare_cartesian_ids(int num_points_1d, std::vector<std::vector<int> > *ids) const;

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


class ProblemDimensions
{
public:
   ProblemDimensions(int mesh_dim, PhysicsType physicsType, const P4estClass &p4est);

   int n_glob_elems;  // number of global elements
   int n_glob_nodes;  // number of global nodes
   int n_glob_dofs;   // number on global degrees of freedom

   int n_subdom_elems;  // number of subdomain elements
   int n_subdom_nodes;  // number of subdomain nodes
   int n_subdom_dofs;   // number on subdomain degrees of freedom

   // number of dofs per node (1 for Laplace, 3 for elasticity)
   int n_node_dofs;

};

typedef std::vector<double> (*exact_fn)(const std::vector<double> &coords);

template <typename number>
class PrintVec
{
public:
   PrintVec<number>(const std::vector<number> &vec, std::string separator = ", ", std::string prefix = "")
      : vec(vec), separator(separator), prefix(prefix) {}
   const std::vector<number> &vec;
   std::string separator, prefix;
};

template <typename number>
std::ostream& operator<<(std::ostream& os, const PrintVec<number>& pv)
{
    for(auto x : pv.vec)
       os << pv.prefix << x << pv.separator;
//    os << std::endl;
    return os;
}

template <typename number>
class PrintVec2D
{
public:
   PrintVec2D<number>(const std::vector<std::vector<number> > &vec, std::string separator = ", ", std::string prefix = "")
                      : vec(vec), separator(separator), prefix(prefix) {}
   const std::vector<std::vector<number> > &vec;
   std::string separator, prefix;
};

template <typename number>
std::ostream& operator<<(std::ostream& os, const PrintVec2D<number>& pv)
{
   for(auto row : pv.vec) {
      os << pv.prefix;
      for(auto x : row) {
         os << x << pv.separator;
      }
      os << std::endl;
   }
   //os << std::endl;
   return os;
}

#endif // DEFINITIONS_H
