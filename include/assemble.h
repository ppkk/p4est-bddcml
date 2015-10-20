#ifndef P4EST_BDDCML_INTERACTION_H
#define P4EST_BDDCML_INTERACTION_H

#include <vector>

#include "definitions.h"
#include "bddcml_structs.h"
#include "my_p4est_interface.h"

#define QUAD_ORDER 2

typedef std::vector<double> (*RhsPtr)(std::vector<double>);

void print_complete_matrix_rhs(const BddcmlFemSpace &femsp, const BddcmlDimensions &global_dims, const SparseMatrix &matrix, const RealArray &rhss, MPI_Comm mpicomm);

void assemble_matrix_rhs(const P4estClass &p4est, const BddcmlMesh &mesh, const BddcmlFemSpace &femsp,
                         SparseMatrix *matrix, RealArray *rhss, RhsPtr rhs_ptr, Parameters params);






#endif // P4EST_BDDCML_INTERACTION_H
