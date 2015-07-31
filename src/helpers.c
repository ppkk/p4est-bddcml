#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <stdio.h>
#include "helpers.h"

int print_rank = 0;
int mpi_rank;

static int num_allocations = 0;

void allocate_sparse_matrix(int length, int type, SparseMatrix *matrix)
{
   matrix->len = length;
   matrix->type = type;
   matrix->is_assembled = 0;
   matrix->i = (int*) malloc(length * sizeof(int));
   matrix->j = (int*) malloc(length * sizeof(int));
   matrix->val = (real*) malloc(length * sizeof(real));

   num_allocations++;
}

void zero_matrix(SparseMatrix *matrix)
{
   matrix->nnz = 0;
   memset(matrix->i, 0, matrix->len * sizeof(int));
   memset(matrix->j, 0, matrix->len * sizeof(int));
   memset(matrix->val, 0, matrix->len * sizeof(double));
}

void add_matrix_entry(SparseMatrix *matrix, int i, int j, real value)
{
   assert(matrix->nnz < matrix->len);
   if((matrix->type == SPD) || (matrix->type == SYM_GENERAL))
   {
      if(i > j)
      {
         int help = i;
         i = j;
         j = help;
      }
   }

   matrix->i[matrix->nnz] = i;
   matrix->j[matrix->nnz] = j;
   matrix->val[matrix->nnz] = value;
   matrix->nnz++;
}

void allocate_idx_array(int length, IdxArray *array)
{
   array->len = length;
   if(length > 0)
   {
      array->val = (int*) malloc(length * sizeof(int));
      num_allocations++;
   }
   else
   {
      array->val = NULL;
   }
}

void allocate_real_array(int length, RealArray *array)
{
   array->len = length;
   if(length > 0)
   {
      array->val = (real*) malloc(length * sizeof(real));
      num_allocations++;
   }
   else
   {
      array->val = NULL;
   }
}

void allocate_real_2D_array(int length1, int length2, Real2DArray *array)
{
   array->len1 = length1;
   array->len2 = length2;
   if(length1 * length2 > 0)
   {
      array->val = (real*) malloc(length1 * length2 * sizeof(real));
      num_allocations++;
   }
   else
   {
      array->val = NULL;
   }
}

void free_sparse_matrix(SparseMatrix *matrix)
{
   free(matrix->val);
   matrix->val = NULL;
   free(matrix->i);
   matrix->i = NULL;
   free(matrix->j);
   matrix->j = NULL;
   matrix->len = 0;

   num_allocations--;
}

void free_idx_array(IdxArray *array)
{
   free(array->val);
   array->val = NULL;
   array->len = 0;
   num_allocations--;
}

void free_real_array(RealArray *array)
{
   free(array->val);
   array->val = NULL;
   array->len = 0;
   num_allocations--;
}

void free_real_2D_array(Real2DArray *array)
{
   free(array->val);
   array->val = NULL;
   array->len1 = 0;
   array->len2 = 0;
   num_allocations--;
}

void print_idx_array(IdxArray *array, char name[])
{
   printf("%s = [", name);
   int i;
   for(i = 0; i < array->len; i++) {
      printf("%d, ", array->val[i]);
   }
   printf("]\n");
}

void print_real_array(RealArray *array, char name[])
{
   printf("%s = [", name);
   int i;
   for(i = 0; i < array->len; i++) {
      printf("%lf, ", array->val[i]);
   }
   printf("]\n");
}

void zero_idx_array(IdxArray *array)
{
   for(int i = 0; i < array->len; i++) {
      array->val[i] = 0;
   }
}

void zero_real_array(RealArray *array)
{
   for(int i = 0; i < array->len; i++) {
      array->val[i] = 0;
   }
}



 int get_num_allocations()
{
   return num_allocations;
}

