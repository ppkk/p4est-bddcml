#ifndef LOCAL_MATRIX_H
#define LOCAL_MATRIX_H

#include "definitions.h"

class IntegrationCell;

// todo: in the future, it could be rectangular (pressure has less dofs than velocity)
// todo: -> ndofs1, ndofs2
class LocalMatrixComponent
{
public:
   LocalMatrixComponent(int ndofs);
   void copy(const LocalMatrixComponent &copy_from);
   void clear();

public:
   int ndofs;
   std::vector<std::vector<real> > mat;
};

class LocalVectorComponent
{
public:
   LocalVectorComponent(int ndofs);
   void copy(const LocalVectorComponent &copy_from);
   void clear();

public:
   int ndofs;
   std::vector<real> vec;
};

class LocalMatrix
{
public:
   LocalMatrix(int ncomponents, int ndofs);
   void clear();

public:
   int ncomponents;
   int ndofs;
   std::vector<std::vector<LocalMatrixComponent> > comps;
};

class LocalVector
{
public:
   LocalVector(int ncomponents, int ndofs);
   void clear();

public:
   int ncomponents;
   int ndofs;
   std::vector<LocalVectorComponent> comps;
};

class HangingInfo
{
public:
   HangingInfo(const P4estClass &p4est);

   void apply_constraints(int elem_idx, const IntegrationCell &cell,
                          const LocalMatrixComponent &in, LocalMatrixComponent *out);
   void apply_constraints(int elem_idx, const IntegrationCell &cell,
                          const LocalMatrix &in, LocalMatrix *out);

   void apply_constraints(int elem_idx, const IntegrationCell &cell,
                          const LocalVectorComponent &in, LocalVectorComponent *out);
   void apply_constraints(int elem_idx, const IntegrationCell &cell,
                          const LocalVector &in, LocalVector *out);

   bool inline is_face_hanging(int face) const {return faces[face] != -1; }
   bool inline is_edge_hanging(int edge) const {return edges[edge] != -1; }

private:
   void init_coefs(int elem_idx, const IntegrationCell &cell);

private:
   int active_elem_idx;
   bool active_elem_anyhang;

   std::vector<int> faces;
   std::vector<int> edges;
   const P4estClass &p4est;

   std::vector<std::vector<double> > coefs;

   friend class P4estClass2D;
   friend class P4estClass3D;
};


#endif // LOCAL_MATRIX_H
