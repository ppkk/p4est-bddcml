#ifndef P4EST_BDDCML_INTERACTION_H
#define P4EST_BDDCML_INTERACTION_H

#include "definitions.h"
#include "bddcml_structs.h"
#include "p4est_common.h"

typedef std::vector<double> (*RhsPtr)(std::vector<double>);

void prepare_subdomain_mesh(p4est_t *p4est, p4est_lnodes_t *lnodes, BddcmlDimensions *subdomain_dims, BddcmlMesh *mesh);

void prepare_subdomain_fem_space(BddcmlMesh *mesh, BddcmlFemSpace *femsp);

void print_complete_matrix_rhs(BddcmlFemSpace *femsp, BddcmlDimensions *global_dims, SparseMatrix *matrix, RealArray *rhss, MPI_Comm mpicomm);

void assemble_matrix_rhs(p4est_lnodes_t *lnodes, BddcmlMesh *mesh, BddcmlFemSpace *femsp,
                         SparseMatrix *matrix, RealArray *rhss, RhsPtr rhs_ptr);






#endif // P4EST_BDDCML_INTERACTION_H
