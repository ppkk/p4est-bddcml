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
   Integrator(const P4estClass &p4est, const NodalElementMesh &mesh, const ReferenceElement & ref_elem, const std::vector<double> &sol);

   double calculate(int quad_order, integral_callback callback, exact_fn function_callback, std::vector<double> *element_values) const;
   double l2_norm(int quad_order, std::vector<double> *element_values = nullptr) const;
   double l2_error(int quad_order, exact_fn exact_solution, std::vector<double> *element_values = nullptr) const;

public:
   const P4estClass &p4est;
   const NodalElementMesh &mesh;
   const ReferenceElement & ref_elem;
   const std::vector<double> &sol;
};


#endif // INTEGRAL_H
