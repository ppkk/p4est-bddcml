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
   ELASTICITY = 1
};

struct Parameters
{
   real mu;
   real lambda;

//   set_params(real mu, real lambda)
//   {
//      this->mu = mu;
//      this->lambda = lambda;
//   }


   Parameters(real young_mod, real poisson_num)
   {
      this->lambda = (young_mod * poisson_num) / ((1 + poisson_num) * (1 - 2 * poisson_num));
      this->mu = young_mod / (2 * (1 + poisson_num));
   }
};


#endif // DEFINITIONS_H
