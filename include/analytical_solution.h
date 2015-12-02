#ifndef ANALYTICAL_SOLUTION_H
#define ANALYTICAL_SOLUTION_H

#include "definitions.h"

enum class ExactSolID
{
   // bi or tri quadratic function with zero BC on ref. square/cube
   LapQuadratic,

   // Internal layer with slope defined in .cc file
   InternalLayer,
};

class ExactSolution
{
public:
   void calculate(ExactSolID sol_id, const std::vector<double> &coords);
public:
   double u, rhs;
   std::vector<double> grad;
};

#endif // ANALYTICAL_SOLUTION_H
