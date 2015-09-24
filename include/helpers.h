#ifndef HELPERS_H
#define HELPERS_H

#include "definitions.h"

extern int mpi_rank;
extern int mpi_size;
extern int print_rank;

#define PPP if(mpi_rank == print_rank)


enum class MatrixType{
   GENERAL,      // general (full storage)
   SPD,          //symmetric positive definite (only triangle stored)
   SYM_GENERAL  // symmetric general (only triangle stored)
};

enum class PhysicsType{
   LAPLACE = 0,
   LINEAR_ELASTICITY = 1
};

// matrix in coordinate format - triplets (i,j,a_ij)
typedef struct SparseMatrix
{
   // numerical properties of the matrix (MUMPS-like notation)
   MatrixType type;

   int len;
   int *i;
   int *j;
   real *val;
   int nnz;

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

// todo: helping 2D array
typedef struct Real2DArray
{
   int len1, len2;
   real **val;
   real *val_serialized;
}
Real2DArray;

void allocate_sparse_matrix(int length, MatrixType type, SparseMatrix* matrix);
void free_sparse_matrix(SparseMatrix* matrix);
void zero_matrix(SparseMatrix* matrix);

// if matrix type is symmetric, it will add the entry to the upper triangle
void add_matrix_entry(SparseMatrix *matrix, int i, int j, real value );

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
