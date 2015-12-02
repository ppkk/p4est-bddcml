#ifndef MY_P4EST_INTERFACE_H
#define MY_P4EST_INTERFACE_H

#include "definitions.h"
#include "p4est_base.h"

class BddcmlMesh;
class ProblemDimensions;
class IntegrationMesh;
class NodalElementMesh;
class HangingInfo;
class ReferenceElement;
class IntegrationCell;

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

   static P4estClass *create(int num_dim, int degree, MPI_Comm mpicomm);

   virtual void plot_solution(int num_components, double* u_sol, double* u_exact, int *partition) const = 0;
   virtual void print_p4est_mesh (int which_rank) const = 0;

   virtual void prepare_dimmensions(ProblemDimensions *problem_dims) const = 0;
   virtual void prepare_bddcml_mesh_global_mappings(BddcmlMesh *mesh) const = 0;
   virtual void prepare_bddcml_mesh_nodes_old(BddcmlMesh *mesh) const = 0;
   virtual void prepare_integration_mesh(IntegrationMesh *mesh) const = 0;

   // todo: rather then ncomponents, it should take vector of ReferenceElements to allow different (velocity + pressure)
   // todo: consider const ref -> (smart) pointers...
   virtual void prepare_nodal_mesh(int ncomponents, const IntegrationMesh &integration_mesh,
                                   const ReferenceElement &reference_element, NodalElementMesh *nodal_mesh) const = 0;

   virtual bool get_hanging_info(int quad_idx, HangingInfo *hanging_info) const = 0;

   virtual void refine_and_partition(int num, RefineType type) = 0;
   virtual void refine_and_partition(const std::vector<double> &element_errors, double refine_approx_elems) = 0;

   virtual void init_definitions(Def *def) const = 0;

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
