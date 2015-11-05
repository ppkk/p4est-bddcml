#ifndef LOCAL_SOLUTION_H
#define LOCAL_SOLUTION_H

#include <memory>

#include "definitions.h"
#include "local_matrix.h"

class NodalElement;
class ReferenceElement;
class Quadrature;

class LocalSolution
{
public:
   LocalSolution(const P4estClass &p4est, const NodalElement &nodal_elem, const ReferenceElement &ref_elem, const double * const sol);

   void get_value_reference_coords(const std::vector<double> &ref_coords, std::vector<double> *values) const;

   // values[quad_pt][component]; grads[quad_pt][component][grad_comp]
   void get_values_in_quadpoints(const Quadrature &q, std::vector<std::vector<double> > *values,
                                 std::vector<std::vector<std::vector<double> > > *grads) const;

public:
   LocalVector loc_vec;
   const ReferenceElement &ref_elem;
   const NodalElement &nodal_elem;
};

#endif // LOCAL_SOLUTION_H
