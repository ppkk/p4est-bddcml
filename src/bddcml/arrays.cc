#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <stdio.h>
#include "bddcml/arrays.h"

static int num_allocations = 0;

void SparseMatrix::allocate(int length, MatrixType type)
{
   this->len = length;
   this->type = type;
   this->is_assembled = 0;
   this->ii = (int*) malloc(length * sizeof(int));
   this->jj = (int*) malloc(length * sizeof(int));
   this->val = (real*) malloc(length * sizeof(real));

   num_allocations++;
}

void SparseMatrix::zero()
{
   nnz = 0;
   memset(ii, 0, len * sizeof(int));
   memset(jj, 0, len * sizeof(int));
   memset(val, 0, len * sizeof(double));
}

void SparseMatrix::add_entry(int i, int j, real value)
{
   assert(nnz < len);
   if((type == MatrixType::SPD) || (type == MatrixType::SYM_GENERAL))
   {
      if(i > j)
      {
//         int help = i;
//         i = j;
//         j = help;
         return;
      }
   }

   ii[nnz] = i;
   jj[nnz] = j;
   val[nnz] = value;
   nnz++;
}


void SparseMatrix::free_matrix()
{
   if((val == nullptr) && (ii == nullptr) && (jj == nullptr))
      return;

   free(val);
   val = nullptr;
   free(ii);
   ii = nullptr;
   free(jj);
   jj = nullptr;
   len = 0;

   num_allocations--;
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
      array->val = nullptr;
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
      array->val = nullptr;
   }
}

void allocate_real_2D_array(int length1, int length2, Real2DArray *array)
{
   array->len1 = length1;
   array->len2 = length2;
   if(length1 * length2 > 0)
   {
      array->val_serialized = (real*) malloc(length1 * length2 * sizeof(real));
      array->val = (real**) malloc(length2 * sizeof(real*));
      for(int i = 0; i < length2; i++)
      {
         array->val[i] = array->val_serialized + length1 * i;
      }
      num_allocations++;
   }
   else
   {
      array->val = nullptr;
      array->val_serialized = nullptr;
   }
}


void free_idx_array(IdxArray *array)
{
   free(array->val);
   array->val = nullptr;
   array->len = 0;
   num_allocations--;
}

void free_real_array(RealArray *array)
{
   free(array->val);
   array->val = nullptr;
   array->len = 0;
   num_allocations--;
}

void free_real_2D_array(Real2DArray *array)
{
   free(array->val);
   array->val = nullptr;
   free(array->val_serialized);
   array->val = nullptr;
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

