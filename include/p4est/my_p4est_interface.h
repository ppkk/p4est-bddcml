#ifndef MY_P4EST_INTERFACE_H
#define MY_P4EST_INTERFACE_H

#include "definitions.h"
#include "p4est_base.h"

class BddcmlMesh;
class BddcmlDimensions;
class GeometryMesh;
class HangingInfo;

enum class RefineType
{
   UNIFORM,
   CIRCLE,
   SQUARE,
   CENTER,
   DIAGONAL,
   POINT
};

class P4estClass
{
public:
   P4estClass(int degree, sc_MPI_Comm mpicomm) : degree(degree), mpicomm(mpicomm) {}
   virtual ~P4estClass() {}

   void init();

   virtual void plot_solution(int num_components, double* u_sol, double* u_exact, int *partition) const = 0;
   virtual void print_p4est_mesh (int which_rank) const = 0;

   virtual void prepare_dimmensions(BddcmlDimensions *subdomain_dims, BddcmlDimensions *global_dims) const = 0;
   virtual void prepare_subdomain_bddcml_mesh(BddcmlMesh* mesh) const = 0;
   virtual void prepare_subdomain_geometry_mesh(GeometryMesh *mesh) const = 0;

   virtual bool get_hanging_info(int quad_idx, HangingInfo *hanging_info) const = 0;

   virtual void refine_and_partition(int num, RefineType type) = 0;

   static P4estClass *create(int num_dim, int degree, MPI_Comm mpicomm);

   static void init_definitions();

   inline int get_num_dim() const {return num_dim; }

protected:
   int degree;
   int num_dim;
   sc_MPI_Comm mpicomm;

   // I cannot store p4est or lnodes here, since it is only a #define and the actuall types are
   // different for 2D and 3D. They are stored in 2 descendants of this class, P4estClass2D and P4estClass3D,
   // which are defined in a single file using p4est tricks
};

#endif // MY_P4EST_INTERFACE_H
