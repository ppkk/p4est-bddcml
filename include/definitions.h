#ifndef DEFINITIONS_H
#define DEFINITIONS_H

// type used for floating point
typedef double real;

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


#endif // DEFINITIONS_H
