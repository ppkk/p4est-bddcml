#ifndef SHAPEFUN_H
#define SHAPEFUN_H

#include <vector>

class Quadrature;

void ref_value_1D(int loc_id_1d, double x, double elem_len, double *value, double *der);
void shape_fun(int order, int idx, double x, double elem_len, double *value, double* der);
void prepare_transformed_values(const Quadrature &q, double element_length,
                                std::vector<std::vector<double> > &values,
                                std::vector<std::vector<std::vector<double> > > &gradients);

#endif // SHAPEFUN_H
