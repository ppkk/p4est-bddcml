#ifndef P4EST_BDDCML_INTERACTION_H
#define P4EST_BDDCML_INTERACTION_H

void prepare_dimmensions(p4est_t *p4est, p4est_lnodes_t *lnodes,
                         BddcmlDimensions *subdomain_dims, BddcmlDimensions *global_dims,
                         sc_MPI_Comm mpicomm);

void prepare_subdomain_mesh(p4est_t *p4est, p4est_lnodes_t *lnodes, BddcmlDimensions *subdomain_dims, BddcmlMesh *mesh, real* volumes);

void prepare_subdomain_fem_space(BddcmlMesh *mesh, BddcmlFemSpace *femsp);

void print_complete_matrix_rhs(BddcmlFemSpace *femsp, BddcmlDimensions *global_dims, SparseMatrix *matrix, RealArray *rhss, MPI_Comm mpicomm);

void assemble_matrix_rhs(p4est_lnodes_t *lnodes, BddcmlMesh *mesh, double *element_volumes, BddcmlFemSpace *femsp,
                         SparseMatrix *matrix, RealArray *rhss);






#endif // P4EST_BDDCML_INTERACTION_H
