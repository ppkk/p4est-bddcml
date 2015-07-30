#include "helpers.h"

void print_array(int size, int* array, char name[])
{
   printf("%s = [", name);
   int i;
   for(i = 0; i < size; i++) {
      printf("%d, ", array[i]);
   }
   printf("]\n");
}

void print_f_array(int size, real* array, char name[])
{
   printf("%s = [", name);
   int i;
   for(i = 0; i < size; i++) {
      printf("%lf, ", array[i]);
   }
   printf("]\n");
}

