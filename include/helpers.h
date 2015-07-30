#ifndef HELPERS_H
#define HELPERS_H

#include "definitions.h"

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

void allocate_idx_array(int length, IdxArray* array);
void allocate_real_array(int length, RealArray* array);
void free_idx_array(IdxArray* array);
void free_real_array(RealArray* array);

void print_idx_array(IdxArray* array, char name[]);
void print_real_array(RealArray* array, char name[]);

int get_num_allocations();


#endif // HELPERS_H
