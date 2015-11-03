#include "local_solution.h"
#include "local_matrix.h"
#include "element.h"
#include "shapefun.h"

using namespace std;

LocalSolution::LocalSolution(const P4estClass &p4est, const NodalElement &nodal_elem,
                             const ReferenceElement &ref_elem, const double * const sol) :
               loc_vec(Def::d()->num_components, Def::d()->num_element_nodes), ref_elem(ref_elem) {
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
//   hanging_info.apply_constraints(nodal_elem.elem_idx, nodal_elem.cell, vec_before, &loc_vec);
}

void LocalSolution::get_value_reference_coords(const std::vector<double> &ref_coords,
                                               std::vector<double> *values) const {
   //ref_elem.fill_transformed_values();
}
