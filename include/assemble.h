#ifndef P4EST_BDDCML_INTERACTION_H
#define P4EST_BDDCML_INTERACTION_H

#include "definitions.h"
#include "bddcml/bddcml_structs.h"

class P4estClass;
class IntegrationMesh;

typedef std::vector<double> (*RhsPtr)(std::vector<double>);

void print_complete_matrix_rhs(const BddcmlFemSpace &femsp, const BddcmlDimensions &global_dims, const SparseMatrix &matrix, const RealArray &rhss, MPI_Comm mpicomm);

void assemble_matrix_rhs(const P4estClass &p4est, const IntegrationMesh &integration_mesh, const BddcmlMesh &bddcml_mesh, const BddcmlFemSpace &femsp,
                         SparseMatrix *matrix, RealArray *rhss, RhsPtr rhs_ptr, Parameters params);


#endif // P4EST_BDDCML_INTERACTION_H
