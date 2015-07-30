#include <stdlib.h>
#include <stdio.h>
#include "helpers.h"

static int num_allocations = 0;

void allocate_idx_array(int length, IdxArray *array)
{
   array->len = length;
   array->val = (int*) malloc(length * sizeof(int));
   num_allocations++;
}

void allocate_real_array(int length, RealArray *array)
{
   array->len = length;
   array->val = (real*) malloc(length * sizeof(real));
   num_allocations++;
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

 int get_num_allocations()
{
   return num_allocations;
}

