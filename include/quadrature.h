#ifndef QUADRATURE_H
#define QUADRATURE_H

#include <vector>

class Element;
class LevelSet;

/*
 *  1D Quadrature rules of orders (capable of exact integration of polynomial of orders) 0-99.
 */
typedef double double2[2];
class GaussTables1D
{
public:
   GaussTables1D();
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

class Quadrature
{
public:
   Quadrature(int dimmension) : dimension(dimmension) {}
   void tensor_product(const Quadrature &quad1, const Quadrature &quad2);
   void print();

   void transform_to_physical(Element* element);

public:
   std::vector<double> weights;
   std::vector<std::vector<double> > coords;

   int dimension;
};

class GaussQuadrature : public Quadrature
{
public:
   GaussQuadrature(int dimension, int order);

};

class EquidistantQuadrature : public Quadrature
{
public:
   EquidistantQuadrature(int dimension, int order);

};

class AdaptiveQuadrature : public Quadrature
{
public:
   AdaptiveQuadrature(int dimension, int order);
   void refine_by_level_set(LevelSet* level_set);
};

#endif // QUADRATURE_H
