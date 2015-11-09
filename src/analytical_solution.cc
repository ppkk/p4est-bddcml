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
      const double slope = 60;
      const double center[3] = {1.25, -0.25, -0.25 };

      double rr = 0.0;
      vector<double> my_coords = coords;
      for(int i = 0; i < Def::d()->num_dim; i++) {
         my_coords[i] -= center[i];
         rr += my_coords[i] * my_coords[i];
      }
      double r = sqrt(rr);

      u = atan(slope * (r - M_PI / 3.));

      double t_u_1 = r * (slope*slope * pow(r - M_PI/3, 2) + 1);
      dudx = slope * coords[0] / t_u_1;
      dudy = slope * coords[1] / t_u_1;

      double t_u_2 = (pow(M_PI - 3.0 * r, 2) * slope * slope + 9.0);

      double laplace = 27.0 * 2.0 * rr * (M_PI - 3.0*r) * pow(slope, 3.0) / (pow(t_u_2,2) * rr) -
            9.0 * rr * slope / (t_u_2 * pow(r,3.0)) + 9.0 * 2.0 * slope / (t_u_2 * r);

      rhs = -laplace;

      break;
   }
   default : {
      assert(0);
      break;
   }
   }
}
