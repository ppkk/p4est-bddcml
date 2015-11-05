#ifndef INTEGRAL_H
#define INTEGRAL_H

#include "definitions.h"

#define sqr(x) ((x)*(x))

class NodalElementMesh;
class ReferenceElement;



typedef void (*exact_fn)(const std::vector<double> &coords,
                         std::vector<double> *result);

typedef double (*integral_callback)(const std::vector<double> &coords,
                                    const std::vector<double> &pt_val,
                                    const std::vector<std::vector<double> > &pt_grad,
                                    exact_fn function);

class Integrator
{
public:
   Integrator(const P4estClass &p4est, const NodalElementMesh &mesh, const ReferenceElement & ref_elem, const double * const sol);

   double calculate(int quad_order, integral_callback callback, exact_fn function_callback) const;
   double l2_norm(int quad_order) const;
   double l2_error(int quad_order, exact_fn exact_solution) const;

public:
   const P4estClass &p4est;
   const NodalElementMesh &mesh;
   const ReferenceElement & ref_elem;
   const double * const sol;
};


#endif // INTEGRAL_H
