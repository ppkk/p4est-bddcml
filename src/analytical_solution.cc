#include <math.h>

#include "analytical_solution.h"

using namespace std;

void ExactSolution::calculate(ExactSolID sol_id, const std::vector<double> &coords) {
   switch(sol_id) {
   case ExactSolID::LapQuadratic : {
      u = 1.;
      for(int dim = 0; dim < Def::d()->num_dim; dim++)
         u *= coords[dim] * (1-coords[dim]);

      double laplace = 0;
      for(int dim = 0; dim < Def::d()->num_dim; dim++) {
         double add = -2;
         for(int dim_add = 0; dim_add < Def::d()->num_dim; dim_add++) {
            if(dim != dim_add)
               add *= coords[dim_add] * (1-coords[dim_add]);
         }
         laplace += add;
      }
      rhs = -laplace;
      break;
   }
   case ExactSolID::InternalLayer : {
      //const double slope = 90;
      const double slope = 60;
      const double center[3] = {1.25, -0.25, -0.25 };

      double rr = 0.0;
      double coords_sum = 0.0;
      vector<double> my_coords = coords;
      for(int i = 0; i < Def::d()->num_dim; i++) {
         my_coords[i] -= center[i];
         rr += my_coords[i] * my_coords[i];
         coords_sum += my_coords[i];
      }
      double r = sqrt(rr);

      u = atan(slope * (r - M_PI / 3.));

      double term = 1 / ( 1 + sqr(slope * (r - M_PI/3)));
      grad.resize(Def::d()->num_dim);
      for(int i= 0; i < Def::d()->num_dim; i++) {
         grad[i] = term * slope / r * my_coords[i];
      }

      double laplace = Def::d()->num_dim * term * slope / r +
                       term *  (- slope / r) +
                       (-2 * pow(slope, 3) * (r-M_PI/3) * term * term);

      rhs = -laplace;

      break;
   }
   default : {
      assert(0);
      break;
   }
   }
}
