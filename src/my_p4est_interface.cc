#include <assert.h>

#include "my_p4est_interface.h"

// this is not too nice, but it is hard to deal with p4est without something like that
#undef P4_TO_P8
// first include p4estClass2D
#include "my_p4est_implementation.h"
#undef MY_P4EST_IMPLEMENTATION_H
#define P4_TO_P8
// second include p4estClass3D
#include "my_p4est_implementation.h"

P4estClass *P4estClass::create(int num_dim, int degree, sc_MPI_Comm mpicomm)
{
   if(num_dim == 2)
      return new P4estClass2D(degree, mpicomm);
   else if(num_dim == 3)
      return new P4estClass3D(degree, mpicomm);
   else
      assert(0);
}

int P4estClass::num_dim;
int P4estClass::children;
