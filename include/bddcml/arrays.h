#ifndef HELPERS_H
#define HELPERS_H

#include "definitions.h"


// matrix in coordinate format - triplets (i,j,a_ij)
class SparseMatrix
{
public:
   ~SparseMatrix() {free_matrix(); }
   void allocate(int length, MatrixType type);
   void free_matrix();
   void zero();

   // if matrix type is symmetric, it will add the entry to the upper triangle
   void add_entry(int i, int j, real value );

public:
   // numerical properties of the matrix (MUMPS-like notation)
   MatrixType type;

   int len;
   int *ii;
   int *jj;
   real *val;
   int nnz;

   int is_assembled;
};

typedef struct IdxArray
{
   int len;
   int *val;
}
IdxArray;

typedef struct RealArray
{
   int len;
   real *val;
}
RealArray;

// todo: helping 2D array
typedef struct Real2DArray
{
   int len1, len2;
   real **val;
   real *val_serialized;
}
Real2DArray;

void allocate_idx_array(int length, IdxArray* array);
void free_idx_array(IdxArray* array);
void zero_idx_array(IdxArray *array);
void print_idx_array(IdxArray* array, char name[]);

void allocate_real_array(int length, RealArray* array);
void free_real_array(RealArray* array);
void zero_real_array(RealArray *array);
void print_real_array(RealArray* array, char name[]);

void allocate_real_2D_array(int length1, int length2, Real2DArray *array);
void free_real_2D_array(Real2DArray* array);

int get_num_allocations();


#endif // HELPERS_H
