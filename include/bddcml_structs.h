#ifndef BDDCML_INTERFACE_H
#define BDDCML_INTERFACE_H

#include <mpi.h>
#include "definitions.h"
#include "helpers.h"

// **************************
// GENERAL BDDCML PARAMETERS:
// **************************
typedef struct BddcmlGeneralParams
{
   // beginning index of arrays ( 0 for C, 1 for Fortran )
   int numbase;

   // Just a direct solve by MUMPS?
   int just_direct_solve_int;

   // verbosity of BDDCML ( 0 - only fatal errors, 1 - mild output, 2 - detailed output )
   int verbose_level;

   // export solution to VTU files?
   int export_solution;

   // what is the name of that file (resp. collection of files)
   char output_file_prefix[255];
}
BddcmlGeneralParams;


// *************************
// KRYLOV METHOD PARAMETERS:
// *************************
typedef struct BddcmlKrylovParams
{
   // Krylov subspace iterative method to be used
   //   -1 - use solver defaults
   //   0 - PCG
   //   1 - BICGSTAB (choose for general symmetric and general matrices)
   //   2 - steepest descent method
   //   5 - direct solve by MUMPS
   int krylov_method;

   // use recycling of Krylov subspace
   //   0 - no recycling used
   //   1 - basis of the Krylov subspace will be orthogonalized and also used for new right hand sides
   int recycling_int;

   // size of the Krylov subspace basis to store
   int max_number_of_stored_vectors;

   // maximum number of iterations of a Krylov subspace method
   int maxit;

   // maximum number of iterations of a Krylov subspace method with non-decreasing residual
   int ndecrmax;

   // relative precision of the Krylov subspace method ||residual||/||right-hand side||
   real tol;
}
BddcmlKrylovParams;


// *******************************
// BDDC PRECONDITIONER PARAMETERS:
// *******************************
typedef struct BddcmlPreconditionerParams
{
   // use default values in preconditioner? In such case, all other parameters are ignored
   int use_preconditioner_defaults;

   // use continuity at corners as constraints?
   int use_corner_constraints;

   // use arithmetic constraints on edges and faces?
   int use_arithmetic_constraints;

   // use adaptive constraints on faces?
   int use_adaptive_constraints;

   // use user constraints? - not used in this example
   int use_user_constraints;

   // what type of weights use on interface?
   //   0 - weights by cardinality
   //   1 - weights by diagonal stiffness
   //   2 - weights based on first row of element data
   //   3 - weights based on dof data
   //   4 - weights by Marta Certikova - unit load
   //   5 - weights by Marta Certikova - unit jump
   //   6 - weights by Schur row sums for whole subdomain
   //   7 - weights by Schur row sums computed face by face
   int weights_type;

   // should parallel division be used (ParMETIS instead of METIS) on the first level?
   int parallel_division;

   // find components of the mesh and handle them as independent subdomains when selecting coarse dofs
   // recommended for unstructured meshes, but could be switched off for these simple cubes
   int find_components_int;
}
BddcmlPreconditionerParams;


// **************************
// BDDCML LEVELS INFORMATION
// **************************
typedef struct BddcmlLevelInfo
{
   int nlevels; // number of levels

   // subdomains in levels
   int lnsublev;
   int *nsublev;
   int nsub_loc_1;
}
BddcmlLevelInfo;


// **************************
// BDDCML CONVERGENCE INFO
// **************************
typedef struct BddcmlConvergenceInfo
{
   int num_iter;
   int converged_reason;
   real condition_number;
}
BddcmlConvergenceInfo;


// **************************
// BDDCML MESH DIMMENSIONS
// **************************
typedef struct BddcmlDimensions
{
   int n_elems;  // number of elements
   int n_nodes;  // number of nodes
   int n_dofs;   // number on degrees of freedom

   // spacial dimension
   int n_problem_dims;

   // topological dimension of elements elements, would be lower for shells or beams
   int n_mesh_dims;

   int n_elem_nodes;

} BddcmlDimensions;



// **************************
// BDDCML MESH CONNECTIVITY
// **************************
typedef struct BddcmlMesh
{
   BddcmlDimensions* subdomain_dims;

   IdxArray elem_node_indices;
   IdxArray num_nodes_of_elem;

   Real2DArray coords;

   IdxArray elem_global_map;
   IdxArray node_global_map;
}
BddcmlMesh;


// **************************
// BDDCML FEM SPACE
// **************************
typedef struct BddcmlFemSpace
{
   BddcmlDimensions* subdomain_dims;

   IdxArray node_num_dofs;
   IdxArray fixs_code;
   RealArray fixs_values;
   IdxArray dofs_global_map;
}
BddcmlFemSpace;




void set_implicit_general_params(BddcmlGeneralParams *params);
void set_implicit_krylov_params(BddcmlKrylovParams *params);
void set_implicit_preconditioner_params(BddcmlPreconditionerParams *params);
void init_levels(int n_subdomains_first_level, BddcmlLevelInfo *level_info);
void init_dimmensions(BddcmlDimensions* dimmensions, int mesh_dim);

void init_mesh(BddcmlDimensions* subdomain_dims, BddcmlMesh* mesh);
void free_mesh(BddcmlMesh* mesh);
void print_bddcml_mesh(BddcmlMesh* mesh, int which_rank);

void init_fem_space(BddcmlDimensions* dims, BddcmlFemSpace* femsp);
void free_fem_space(BddcmlFemSpace* femsp);


// **************************************************************
// INTERFACE TO BDDCML FUNCTIONS USING DEFINED STRUCTURES
// **************************************************************
void bddcml_init(BddcmlGeneralParams *general_params, BddcmlLevelInfo *level_info, MPI_Comm communicator);
void bddcml_upload_subdomain_data(BddcmlDimensions *global_dims, BddcmlDimensions *subdomain_dims,
                                  int isub, BddcmlMesh *mesh, BddcmlFemSpace *femsp,
                                  RealArray *rhss, int is_rhs_complete, RealArray *sols, SparseMatrix *matrix,
                                  Real2DArray *user_constraints, Real2DArray *element_data,
                                  RealArray *dof_data, BddcmlPreconditionerParams* preconditioner_params);
void bddcml_setup_preconditioner(MatrixType matrixtype, BddcmlPreconditionerParams *params);
void bddcml_solve(BddcmlKrylovParams *krylov_params, BddcmlConvergenceInfo *convergence_info, MPI_Comm communicator);
void bddcml_download_local_solution(int isub, RealArray *sols);
void bddcml_dotprod_subdomain(int isub, RealArray *sols1, RealArray *sols2, real *normRn2_sub);
void bddcml_finalize();

#endif // BDDCML_INTERFACE_H
