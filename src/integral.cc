#include <math.h>
#include <mpi.h>

#include "integral.h"
#include "element.h"
#include "shapefun.h"
#include "quadrature.h"
#include "local_solution.h"

using namespace std;

double l2_norm_callback(const std::vector<double> &coords,
                        const vector<double> &pt_val,
                        const vector<vector<double> > &pt_grad,
                        exact_fn function) {
   assert((int)pt_val.size() == Def::d()->num_components);
   double result = 0;
   for(double val : pt_val) {
      result += val * val;
   }
   return result;
}

double l2_error_callback(const std::vector<double> &coords,
                         const vector<double> &pt_val,
                         const vector<vector<double> > &pt_grad,
                         exact_fn function) {
   assert((int)pt_val.size() == Def::d()->num_components);
   double result = 0;
   vector<double> fun_val(pt_val.size(), 0.0);
   function(coords, &fun_val);
   for(unsigned comp = 0; comp < pt_val.size(); comp++) {
      result += sqr(pt_val[comp] - fun_val[comp]);
   }
   return result;
}


Integrator::Integrator(const P4estClass &p4est, const NodalElementMesh &mesh,
                       const ReferenceElement &ref_elem, const double * const sol)
               : p4est(p4est), mesh(mesh), ref_elem(ref_elem), sol(sol) {

}


double Integrator::calculate(int quad_order, integral_callback callback, exact_fn function_callback) const {
   double result = 0.0;
   GaussQuadrature quad(Def::d()->num_dim, quad_order);
   vector<vector<double> > values;
   vector<vector<vector<double> > > grads;
   for(const NodalElement& elem : mesh.elements) {
      Quadrature q_transformed(Def::d()->num_dim);
      quad.transform_to_physical(elem.cell, &q_transformed);

      LocalSolution loc_sol(p4est, elem, ref_elem, sol);
      loc_sol.get_values_in_quadpoints(quad, &values, &grads);
      assert(values.size() == quad.np());
      assert(grads.size() == quad.np());
      for(unsigned q_id = 0; q_id < quad.np(); q_id++) {
         result += q_transformed.weights[q_id] * callback(q_transformed.coords[q_id], values[q_id], grads[q_id], function_callback);
      }
   }

   MPI_Allreduce(MPI_IN_PLACE, &result, 1, MPI_DOUBLE, MPI_SUM, mpicomm);

   return result;
}

double Integrator::l2_norm(int quad_order) const {
   return sqrt(calculate(quad_order, l2_norm_callback, nullptr));
}

double Integrator::l2_error(int quad_order, exact_fn exact_solution) const {
   return sqrt(calculate(quad_order, l2_error_callback, exact_solution));
}
