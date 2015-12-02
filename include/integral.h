#ifndef INTEGRAL_H
#define INTEGRAL_H

#include "definitions.h"

class NodalElementMesh;
class ReferenceElement;

typedef double (*integral_callback)(const std::vector<double> &coords,
                                    const std::vector<double> &pt_val,
                                    const std::vector<std::vector<double> > &pt_grad,
                                    exact_fn function);

class Integrator
{
public:
   Integrator(const P4estClass &p4est, const NodalElementMesh &mesh, const ReferenceElement & ref_elem, const std::vector<double> &sol);

   double l2_norm(int quad_order, std::vector<double> *element_values = nullptr) const;
   double l2_error(int quad_order, exact_fn exact_solution, std::vector<double> *element_values = nullptr) const;

   double h1_seminorm(int quad_order, std::vector<double> *element_values = nullptr) const;
   double h1_semierror(int quad_order, exact_fn exact_gradient, std::vector<double> *element_values = nullptr) const;

private:
   double calculate(int quad_order, integral_callback callback, exact_fn function_callback, std::vector<double> *element_values) const;

private:
   const P4estClass &p4est;
   const NodalElementMesh &mesh;
   const ReferenceElement & ref_elem;
   const std::vector<double> &sol;
};

class IntegralResults
{
public:
   void init(const ProblemDimensions &problem_dims);
   double h1_norm() {return sqrt(sqr(l2_norm) + sqr(h1_seminorm)); }
   double h1_error() {return sqrt(sqr(l2_error) + sqr(h1_semierror)); }
   double l2_rel_error() {return l2_error / l2_norm; }
   double h1_rel_error() {return h1_error() / h1_norm(); }
   double h1_rel_semierror() {return h1_semierror / h1_seminorm; }

public:
   int n_glob_dofs;
   int n_glob_elems;

   double l2_norm, l2_error;
   double h1_seminorm, h1_semierror;
   std::vector<double> element_l2_error;
   std::vector<double> element_h1_semierror;
};

#endif // INTEGRAL_H
