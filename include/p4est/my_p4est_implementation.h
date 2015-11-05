
#ifndef MY_P4EST_IMPLEMENTATION_H
#define MY_P4EST_IMPLEMENTATION_H

#include <vector>

#ifndef P4_TO_P8
#include <p4est_bits.h>
#include <p4est_ghost.h>
#include <p4est_lnodes.h>
#include <p4est_vtk.h>
#else
#include <p8est_bits.h>
#include <p8est_ghost.h>
#include <p8est_lnodes.h>
#include <p8est_vtk.h>
#endif

#include "my_p4est_interface.h"

// this file defines both P4estClass2D and P4estClass3D using p4est macros
#ifndef P4_TO_P8
   #define P4estClassDim P4estClass2D
#else
   #define P4estClassDim P4estClass3D
#endif

class P4estClassDim : public P4estClass
{

public:
   P4estClassDim(int degree, sc_MPI_Comm mpicomm);
   virtual ~P4estClassDim();
   void init();


   virtual void prepare_dimmensions(ProblemDimensions *subdomain_dims, ProblemDimensions *global_dims) const;
   virtual void prepare_bddcml_mesh_global_mappings(BddcmlMesh *mesh) const;


   virtual void prepare_integration_mesh(IntegrationMesh *mesh) const;
   virtual void prepare_nodal_mesh(int ncomponents, const IntegrationMesh &integration_mesh,
                                   const ReferenceElement &reference_element, NodalElementMesh *nodal_mesh) const;

   virtual bool get_hanging_info(int quad_idx, HangingInfo *hanging_info) const;

   virtual void refine_and_partition(int num, RefineType type);

   virtual void init_definitions(Def* def) const;

   virtual void print_p4est_mesh (int which_rank) const;

   /*@unused@*/
   virtual void plot_solution(int num_components, double* u_sol, double* u_exact, int *partition) const;

   /*@unused@*/
   virtual void prepare_bddcml_mesh_nodes_old(BddcmlMesh *mesh) const;

protected:
   p4est_gloidx_t node_loc_to_glob(p4est_locidx_t loc_idx) const;
   void get_node_coords(p4est_topidx_t tree, const p4est_quadrant_t &node, std::vector<double> *coords) const;

   void fill_integration_cell(const p4est_quadrant_t *quad, p4est_topidx_t tree, IntegrationCell *cell) const;

   // should be called at the beginning of each function, which uses lnodes
   void update_lnodes() const ;

   // it is done whenewer the mesh changes (lnodes are not up-to-date any more)
   // actually every non-const method should call it
   void destroy_lnodes();

private:
   p4est_connectivity_t *conn;
   p4est_t *p4est;

   // access through the function allways
   inline p4est_lnodes_t *lnodes() const {update_lnodes(); return __lnodes_possibly_not_consistent; }
   mutable p4est_lnodes_t *__lnodes_possibly_not_consistent;
};

#undef P4estClassDim

#endif // MY_P4EST_IMPLEMENTATION_H
