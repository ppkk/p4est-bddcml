#ifndef LOCAL_MATRIX_H
#define LOCAL_MATRIX_H

#include <vector>

#include "definitions.h"

class P4estClass;

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
   HangingInfo(P4estClass &p4est);

   void apply_constraints(int elem_idx, const LocalMatrixComponent &in, LocalMatrixComponent *out);

private:
   void init_coefs(int elem_idx);

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
