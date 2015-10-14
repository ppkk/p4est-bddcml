#ifndef QUADRATURE_H
#define QUADRATURE_H

#include <vector>

/*
 *  1D Quadrature rules of orders (capable of exact integration of polynomial of orders) 0-99.
 */
typedef double double2[2];
class GaussQuad1D
{
public:
   GaussQuad1D();
   inline double2* get_points(int order) const { return tables[order]; }
   inline unsigned char get_num_points(int order) const { return np[order]; }

   inline unsigned short get_max_order() const { return max_order; }
   inline double get_ref_vertex(int n) const { return ref_vert[n]; }

protected:

   double2** tables;
   unsigned char* np;

   double ref_vert[2];
   unsigned short max_order;

};

struct Quadrature
{
   std::vector<double> weights;
   std::vector<std::vector<double> > coords;

   Quadrature(int dimension, int order, double element_length);
   void product(const Quadrature &quad1, const Quadrature &quad2);
   void print();

};


#endif // QUADRATURE_H
