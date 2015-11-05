#include "local_solution.h"
#include "local_matrix.h"
#include "element.h"
#include "shapefun.h"
#include "quadrature.h"
#include "integration_cell.h"

using namespace std;

LocalSolution::LocalSolution(const P4estClass &p4est, const NodalElement &nodal_elem,
                             const ReferenceElement &ref_elem, const double * const sol) :
               loc_vec(Def::d()->num_components, Def::d()->num_element_nodes),
               ref_elem(ref_elem), nodal_elem(nodal_elem){
   LocalVector vec_before(Def::d()->num_components, Def::d()->num_element_nodes);
   assert(nodal_elem.ncomponents() == Def::d()->num_components);
   int nnodes = nodal_elem.components[0].ndofs();
   for(int comp = 0; comp < nodal_elem.ncomponents(); comp++) {
      assert(vec_before.comps[comp].ndofs == nnodes);
      assert(nodal_elem.components[comp].ndofs() == nnodes);
      for(int node = 0; node < nnodes; node++) {
         vec_before.comps[comp].vec[node] = sol[nodal_elem.components[comp].dofs[node]];
      }
   }

   HangingInfo hanging_info(p4est);
   hanging_info.apply_constraints_inverse(nodal_elem.elem_idx, nodal_elem.cell, vec_before, &loc_vec);
}

void LocalSolution::get_value_reference_coords(const std::vector<double> &ref_coords,
                                               std::vector<double> *values) const {
   //ref_elem.fill_transformed_values();
}

void LocalSolution::get_values_in_quadpoints(const Quadrature &q, vector<vector<double> > *values,
                                             vector<vector<vector<double> > > *grads) const {
   vector<vector<double> > shape_values;
   vector<vector<vector<double> > > shape_grads;
   ref_elem.calc_transformed_values(q, nodal_elem.cell.size, &shape_values, &shape_grads);

   *values = vector<vector<double> >(q.np(), vector<double>(nodal_elem.ncomponents(), 0.0));
   *grads = vector<vector<vector<double > > >(q.np(), vector<vector<double> >
                                              (nodal_elem.ncomponents(), vector<double>(Def::d()->num_dim, 0.0)));

   for(unsigned q_idx = 0; q_idx < q.np(); q_idx++) {
      for(unsigned node_idx = 0; node_idx < nodal_elem.nodes.size(); node_idx++) {
         for(int comp = 0; comp < nodal_elem.ncomponents(); comp++) {
            (*values)[q_idx][comp] += loc_vec.comps[comp].vec[node_idx] * shape_values[node_idx][q_idx];
            for(int dim = 0; dim < Def::d()->num_dim; dim++) {
               (*grads)[q_idx][comp][dim] += loc_vec.comps[comp].vec[node_idx] * shape_grads[node_idx][q_idx][dim];
            }
         }
      }
   }
}
