#ifndef MY_P4EST_INTERFACE_H
#define MY_P4EST_INTERFACE_H

#include "definitions.h"
#include "p4est_base.h"

class BddcmlMesh;
class BddcmlDimensions;

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

   virtual void prepare_bddcml_subdomain_mesh(BddcmlMesh* mesh) const = 0;
   virtual void prepare_dimmensions(PhysicsType physicsType,
                                    BddcmlDimensions *subdomain_dims, BddcmlDimensions *global_dims) const = 0;
   virtual int independent_nodes(p4est_locidx_t quadrant, int lnode, p4est_locidx_t *nodes, real* coeffs) const = 0;

   virtual void refine_and_partition(int num, RefineType type) = 0;

   static P4estClass *create(int num_dim, int degree, MPI_Comm mpicomm);
   static int num_dim;
   static int children;

protected:
   virtual p4est_gloidx_t node_loc_to_glob(p4est_locidx_t loc_idx) const = 0;

protected:
   int degree;
   sc_MPI_Comm mpicomm;
};

#endif // MY_P4EST_INTERFACE_H
