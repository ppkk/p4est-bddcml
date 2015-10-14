#ifndef P4EST_BDDCML_INTERACTION_H
#define P4EST_BDDCML_INTERACTION_H

#include "definitions.h"
#include "bddcml_structs.h"
#include "my_p4est_interface.h"

#define QUAD_ORDER 2

typedef std::vector<double> (*RhsPtr)(std::vector<double>);

void print_complete_matrix_rhs(BddcmlFemSpace *femsp, BddcmlDimensions *global_dims, SparseMatrix *matrix, RealArray *rhss, MPI_Comm mpicomm);

void assemble_matrix_rhs(p4est_lnodes_t *lnodes, BddcmlMesh *mesh, BddcmlFemSpace *femsp,
                         SparseMatrix *matrix, RealArray *rhss, RhsPtr rhs_ptr, Parameters params);






#endif // P4EST_BDDCML_INTERACTION_H
