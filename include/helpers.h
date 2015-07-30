#ifndef HELPERS_H
#define HELPERS_H

#include "definitions.h"

// general (full storage)
#define GENERAL 0
//symmetric positive definite (only triangle stored)
#define SPD 1
// symmetric general (only triangle stored)
#define SYM_GENERAL 2

// matrix in coordinate format - triplets (i,j,a_ij)
typedef struct SparseMatrix
{
   // numerical properties of the matrix (MUMPS-like notation)
   int type;

   int len;
   int *i;
   int *j;
   real *val;

   int is_assembled;
}
SparseMatrix;

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

typedef struct Real2DArray
{
   int len1, len2;
   real *val;
}
Real2DArray;

void allocate_sparse_matrix(int length, int type, SparseMatrix* matrix);
void free_sparse_matrix(SparseMatrix* matrix);

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
