#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include "bddcml_interface_c.h"

#include "bddcml_structs.h"
#include "helpers.h"

// *************************
// KRYLOV METHOD PARAMETERS:
// *************************

// Krylov subspace iterative method to be used
//   -1 - use solver defaults
//   0 - PCG
//   1 - BICGSTAB (choose for general symmetric and general matrices)
//   2 - steepest descent method
//   5 - direct solve by MUMPS
   /*const*/ int krylov_method = 0; 

// use recycling of Krylov subspace
//   0 - no recycling used
//   1 - basis of the Krylov subspace will be orthogonalized and also used for new right hand sides
   /*const*/ int recycling_int = 1;
   
// size of the Krylov subspace basis to store
   /*const*/ int max_number_of_stored_vectors = 50;

// maximum number of iterations of a Krylov subspace method
   /*const*/ int maxit = 500;

// maximum number of iterations of a Krylov subspace method with non-decreasing residual
   /*const*/ int ndecrmax = 50;

// relative precision of the Krylov subspace method ||residual||/||right-hand side||
   /*const*/ real tol = 1.e-6;

// *******************************
// BDDC PRECONDITIONER PARAMETERS:
// *******************************

// use default values in preconditioner? In such case, all other parameters are ignored
   /*const*/ int use_preconditioner_defaults = 0;

// use continuity at corners as constraints?
   /*const*/ int use_corner_constraints = 1;

// use arithmetic constraints on edges and faces?
   /*const*/ int use_arithmetic_constraints = 1;

// use adaptive constraints on faces?
   /*const*/ int use_adaptive_constraints = 0;

// use user constraints? - not used in this example
   /*const*/ int use_user_constraints = 0;

// what type of weights use on interface?
//   0 - weights by cardinality
//   1 - weights by diagonal stiffness
//   2 - weights based on first row of element data
//   3 - weights based on dof data
//   4 - weights by Marta Certikova - unit load
//   5 - weights by Marta Certikova - unit jump
//   6 - weights by Schur row sums for whole subdomain
//   7 - weights by Schur row sums computed face by face
   /*const*/ int weights_type = 0;

// should parallel division be used (ParMETIS instead of METIS) on the first level?
   /*const*/ int parallel_division = 1;

// find components of the mesh and handle them as independent subdomains when selecting coarse dofs 
// recommended for unstructured meshes, but could be switched off for these simple cubes
   /*const*/ int find_components_int = 1;

// *******************
// PROBLEM PARAMETERS:
// *******************
// numerical properties of the matrix (MUMPS-like notation)
//   0 - general (full storage)
//   1 - symmetric positive definite (only triangle stored)
//   2 - symmetric general (only triangle stored)
   /*const*/ int matrixtype = 1;

// assuming tri-linear hexahedral finite elements
//     z
//   ^
//   |                                                                        
//   |                                                                        
//   5----------8            
//   |\         |\                 |          
//   | \        | \                |          
//   |  \       |  \               |       
//   |   6----------7        
//   |   |      |   |        
//   1---|------4---|--> y   
//    \  |       \  |        
//     \ |        \ |        
//      \|         \|        
//       2----------3        
//        \                       |
//         \                      |
//         `'
//           x

// number of degrees of freedom on element
#define NDOF_PER_ELEMENT 8
   
// element matrix on reference element [0,1]^3
    const  real element_matrix_ref[NDOF_PER_ELEMENT * NDOF_PER_ELEMENT] = 
      { 1./3. , 0.       ,-1./12., 0.       , 0.       ,-1./12.,-1./12.,-1./12.,
         0.       , 1./3. , 0.       ,-1./12.,-1./12., 0.       ,-1./12.,-1./12.,
        -1./12., 0.       , 1./3. , 0.       ,-1./12.,-1./12., 0.       ,-1./12.,
         0.       ,-1./12., 0.       , 1./3. ,-1./12.,-1./12.,-1./12., 0.       ,
         0.       ,-1./12.,-1./12.,-1./12., 1./3. , 0.       ,-1./12., 0.       ,
        -1./12., 0.       ,-1./12.,-1./12., 0.       , 1./3. , 0.       ,-1./12.,
        -1./12.,-1./12., 0.       ,-1./12.,-1./12., 0.       , 1./3. , 0.       ,
        -1./12.,-1./12.,-1./12., 0.       , 0.       ,-1./12., 0.       , 1./3. };

// spacial dimension
   /*const*/ int ndim = 3;

// topological dimension of elements elements, would be lower for shells or beams
   /*const*/ int meshdim = 3;

   const char routine_name[] = "POISSON_ON_CUBE";

// input read from command line
   int num_el_per_sub_edge, num_sub_per_cube_edge; // basic properties of the cubes-in-cubes problem

   int num_el_per_cube_edge;
   real hsize, el_vol;                              // elements size and volume

// parallel variables
   MPI_Comm comm_all;
   int myid, nproc, ierr;

   int nsub;  // number of subdomains on the first level 
   int nelem; // number of elements 
   int ndof;  // number of degrees of freedom 
   int nnod;  // number of nodes
   
// scaled element matrix
   real element_matrix[NDOF_PER_ELEMENT * NDOF_PER_ELEMENT];

// local subdomain data
   int nelems;  // subdomain number of elements
   int ndofs;   // subdomain number on degrees of freedom
   int nnods;   // subdomain number of nodes
   int linets,   lnnets,   lnndfs;
   int *inets, *nnets, *nndfs;
   int lxyzs1,   lxyzs2;
   real *xyzs;
   int lifixs;
   int *ifixs;
   int lfixvs;
   real *fixvs;
   int lrhss;
   real *rhss;
   int lsols;
   real *sols;
   int lisegns, lisngns, lisvgvns;
   int *isegns, *isngns, *isvgvns;

// matrix in coordinate format - triplets (i,j,a_ij)
   int la;
   int *i_sparse;
   int *j_sparse;
   real *a_sparse;

// user constraints - not really used here
   int luser_constraints1;
   int luser_constraints2;
   real *user_constraints;

// data for elements - not really used here
   int lelement_data1;
   int lelement_data2;
   real *element_data;

// data for dofs - not really used here
   int ldof_data;
   real *dof_data;

// data about resulting convergence
   int num_iter, converged_reason;
   real condition_number;
   real normRn_sol, normRn2, normRn2_loc, normRn2_sub;
   real normL2_sol, normL2_loc, normL2_sub;
   real normLinf_sol, normLinf_loc;

// time variables
   real t_total, t_load, t_pc_setup, t_krylov;

// small variables - indices, etc.
   int ia, indinets, nne, idof, jdof, lelm;
   int ie, i, isub, j;
   int is_rhs_complete_int;
   int is_assembled_int;
   char aux[32];
   char command[256];
   int idx;


int main(int argc, char **argv)
{
   BddcmlGeneralParams general_params;
   BddcmlLevelInfo level_info;

   set_implicit(&general_params);

   // MPI initialization
//***************************************************************PARALLEL
   ierr = MPI_Init(NULL, NULL);
   // Communicator
   comm_all = MPI_COMM_WORLD;
   ierr = MPI_Comm_rank(comm_all, &myid);
   ierr = MPI_Comm_size(comm_all, &nproc);
//***************************************************************PARALLEL

// Initial screen
   if (myid == 0) {
      printf("===========Possion on cube solver===========\n");
      printf("| Solves problem                             |\n");
      printf("|        -laplace u = 1 on D = [0,1]^3,            |\n");
      printf("|           u = 0 on dD,                     |\n");
      printf("| using FEM and the BDDCML solver.           |\n");
      printf(" ============================================ \n");
   }

   // Number of elements in an edge of a subdomain and number of subdomains in an edge of the unit cube
   if ( myid == 0 ) {
      if(argc == 3 + 1) {
         num_el_per_sub_edge = atoi(argv[1]);
         num_sub_per_cube_edge = atoi(argv[2]);
         level_info.nlevels = atoi(argv[3]);
      }
      else {
         printf(" Usage: mpirun -np X ./poisson_on_cube NUM_EL_PER_SUB_EDGE NUM_SUB_PER_CUBE_EDGE NLEVELS");
         exit(0);
         //call error(routine_name,'trouble getting problem sizes')
      }
   }
      
   // Broadcast of name of the problem      
//***************************************************************PARALLEL
   ierr = MPI_Bcast(&num_el_per_sub_edge,   1, MPI_INT,   0, comm_all);
   ierr = MPI_Bcast(&num_sub_per_cube_edge, 1, MPI_INT,   0, comm_all);
   ierr = MPI_Bcast(&level_info.nlevels,    1, MPI_INT,   0, comm_all);
//***************************************************************PARALLEL

   // measuring time
   ierr = MPI_Barrier(comm_all);
     
   //TODO: timing
   //call time_start

   // number of elements on an edge of the unit cube
   num_el_per_cube_edge = num_el_per_sub_edge * num_sub_per_cube_edge;
   // element size
   hsize = 1./ num_el_per_cube_edge;
   // element volume
   el_vol = pow(hsize, ndim);
   // total number of elements
   nelem = pow(num_el_per_cube_edge, ndim);
   // total number of subdomains
   nsub  = pow(num_sub_per_cube_edge, ndim);
   // total number of nodes
   nnod  = pow(num_el_per_cube_edge + 1, ndim);
   // total number of degrees of freedom - equal to number of nodes for scalar problem
   ndof  = nnod;

   initialize_levels(nsub, &level_info);

   printf("YYY\nYYY\n nu sublev, %d, %d, %d\n", level_info.nsublev[0], level_info.nsublev[1], level_info.nsublev[2]);

// Basic properties 
   if (myid == 0) {
      printf("Characteristics of the problem :\n");
      printf("  number of processors            nproc = %d\n" ,nproc);
      printf("  number of dimensions             ndim = %d\n", ndim);
      printf("  mesh dimension                meshdim = %d\n", meshdim);
      printf("  number of elements global       nelem = %d\n", nelem);
      printf("  number of subdomains             nsub = %d\n", nsub);
      printf("  number of nodes global           nnod = %d\n", nnod);
      printf("  number of DOF                    ndof = %d\n", ndof);
      printf("  number of levels              nlevels = %d\n", level_info.nlevels);
      printf("  number of subdomains in levels        = ");
      for(idx = 0; idx < level_info.nlevels; idx++) {
         printf("%d, ", level_info.nsublev[idx]);
      }
      printf("\n");
      printf("Characteristics of iterational process:\n");
      printf("  tolerance of error                tol = %lf\n", tol);
      printf("  maximum number of iterations    maxit = %d\n", maxit);
      printf("  number of incresing residual ndecrmax = %d\n", ndecrmax);
      printf("  using recycling of Krylov method ?      %d\n", recycling_int);
   }
   if (myid == 0) {
         printf("Initializing BDDCML ...");
   }
   // tell me how much subdomains should I load
   level_info.nsub_loc_1 = -1;
   
   bddcml_init(&general_params, &level_info, comm_all);
   if (myid == 0) {
      printf("Initializing BDDCML done.\n");
   }

   // processor i takes care of subdomains sub2proc[i], ... , sub2proc[i+1] - 1
   int lsub2proc = nproc + 1;
   int sub2proc[lsub2proc];
        
   // first find out starting subdomain for each processor, it was given by bddcml_init
   // since there is less processors than subdomains, we can store it to sub2proc
//***************************************************************PARALLEL
   ierr = MPI_Allgather( &level_info.nsub_loc_1, 1, MPI_INT, sub2proc, 1, MPI_INT, comm_all);
//***************************************************************PARALLEL
   //the array now contains counts, change it to starts
   for(i = 1; i <= nproc; i++) {
      sub2proc[i] = sub2proc[i-1] + sub2proc[i];
   }

   // shift it one back and add one 
   for(i = nproc; i > 0; i--) {
      sub2proc[i] = sub2proc[i-1];
   }
   sub2proc[0] = 0;
   
   
   // create and load subdomains
   if (myid == 0) {
      printf("Loading data ...\n");
   }
   //TODO: call time_start
   // Loop over subdomains and load them to BDDCML
   for(isub = sub2proc[myid]; isub <  sub2proc[myid+1]; isub++) {
      // create mesh for subdomain
      nelems   = pow(num_el_per_sub_edge, ndim);
      nnods    = pow(num_el_per_sub_edge+1, ndim);
      ndofs    = nnods;
      linets   = nelems * 8;
      lnnets   = nelems;
      lnndfs   = nnods;
      lisegns  = nelems;
      lisngns  = nnods;
      lisvgvns = ndofs;
      inets = (int*) malloc(linets * sizeof(int));
      nnets = (int*) malloc(lnnets * sizeof(int));
      nndfs = (int*) malloc(lnndfs * sizeof(int));
      isegns = (int*) malloc(lisegns * sizeof(int));
      isngns = (int*) malloc(lisngns * sizeof(int));
      isvgvns = (int*) malloc(lisvgvns * sizeof(int));

      lxyzs1   = nnods;
      lxyzs2   = ndim;
      xyzs = (real*) malloc(lxyzs1 * lxyzs2 * sizeof(real));
      lifixs   = nnods;
      lfixvs   = nnods;
      ifixs = (int*) malloc(lifixs * sizeof(int));
      fixvs = (real*) malloc(lfixvs * sizeof(real));

      // create subdomain mesh and boundary conditions
      prepare_subdomain_data(isub, num_sub_per_cube_edge, num_el_per_sub_edge, hsize,
                                    inets,linets, nnets,lnnets, nndfs,lnndfs,
                                    isegns,lisegns, isngns,lisngns, isvgvns,lisvgvns,
                                    xyzs,lxyzs1,lxyzs2,
                                    ifixs,lifixs, fixvs,lfixvs);

      // create local right hand side
      lrhss = ndofs;
      rhss = (real*) malloc(lrhss * sizeof(real));
      for(idx = 0; idx < lrhss; idx++) {
         if(ifixs[idx] > 0) {
            // rhs is zero on boundary
            rhss[idx] = 0.;
         }
         else {
            rhss[idx] = 1. * el_vol;
         }
      }
      is_rhs_complete_int = 1;

      // create local initial solution
      lsols = ndofs;
      sols = (real*) malloc(lsols * sizeof(real));
      for(idx = 0; idx < lsols; idx++) {
         sols[idx] = 0.;
      }

      // create local subdomain matrix for each subdomain
      // full element matrix scaled based on element size
      for(idx = 0; idx < NDOF_PER_ELEMENT * NDOF_PER_ELEMENT; idx++) {
         element_matrix[idx] = hsize * element_matrix_ref[idx];
      }

      // how much space the upper triangle of the element matrix occupies
      lelm = NDOF_PER_ELEMENT * (NDOF_PER_ELEMENT + 1) / 2;
      // space for all upper triangles of element matrics
      la = nelems*lelm;
      i_sparse = (int*) malloc(la * sizeof(int));
      j_sparse = (int*) malloc(la * sizeof(int));
      a_sparse = (real*) malloc(la * sizeof(real));

      // copy the upper triangle of the element matrix to the sparse triplet
      ia = 0;
      indinets = 0;
      for(ie = 0; ie < nelems; ie++) {
         nne = nnets[ie];
         for(j = 0; j < NDOF_PER_ELEMENT; j++) {
            jdof = inets[indinets + j];
            for(i = 0; i <= j; i++) {
               idof = inets[indinets + i];

               if (idof <= jdof) {
                  i_sparse[ia] = idof;
                  j_sparse[ia] = jdof;
               }
               else {
                  // transpose the entry
                  i_sparse[ia] = jdof;
                  j_sparse[ia] = idof;
               }
               a_sparse[ia] = element_matrix[j * NDOF_PER_ELEMENT + i];

               ia = ia + 1;
            }
         }
         indinets = indinets + nne;
      }
      is_assembled_int = 0;

      // prepare user constraints - not really used here
      luser_constraints1 = 0;
      luser_constraints2 = 0;
      user_constraints = NULL;

      // prepare element data - not really used
      lelement_data1 = 0;
      lelement_data2 = 0;
      element_data = NULL;

      // prepare dof data - not really used
      ldof_data = 0 ;
      dof_data = NULL;

//       if(myid == 0)
//       {
//          print_array(linets, inets, "inets");
//          print_array(lnnets, nnets, "nnets");
//          print_array(la, i_sparse, "i");
//          print_array(la, j_sparse, "j");
//          print_f_array(la, a_sparse, "a");
//       }

      bddcml_upload_subdomain_data_c(&nelem, &nnod, &ndof, &ndim, &meshdim,
                                          &isub, &nelems, &nnods, &ndofs,
                                          inets, &linets, nnets, &lnnets, nndfs, &lnndfs,
                                          isngns, &lisngns, isvgvns, &lisvgvns, isegns, &lisegns,
                                          xyzs, &lxyzs1, &lxyzs2,
                                          ifixs, &lifixs, fixvs, &lfixvs,
                                          rhss, &lrhss, &is_rhs_complete_int,
                                          sols, &lsols,
                                          &matrixtype, i_sparse, j_sparse, a_sparse, &la, &is_assembled_int,
                                          user_constraints, &luser_constraints1, &luser_constraints2,
                                          element_data, &lelement_data1, &lelement_data2,
                                          dof_data, &ldof_data, &find_components_int);

      free(inets);
      free(nnets);
      free(nndfs);
      free(isegns);
      free(isngns);
      free(isvgvns);
      free(xyzs);
      free(ifixs);
      free(fixvs);
      free(rhss);
      free(sols);
      free(i_sparse);
      free(j_sparse);
      free(a_sparse);
   } // end loop over subdomains

   ierr = MPI_Barrier(comm_all);
   //TODO: call time_end(t_load)
   if (myid == 0) {
      printf("Loading data done.\n");
   }

   if (myid == 0) {
      printf("Preconditioner set-up ...\n");
   }

   // PRECONDITIONER SETUP
   ierr = MPI_Barrier(comm_all);
   // TODO: call time_start
   bddcml_setup_preconditioner_c(&matrixtype,
                                 &use_preconditioner_defaults,
                                 &parallel_division,
                                 &use_corner_constraints,
                                 &use_arithmetic_constraints,
                                 &use_adaptive_constraints,
                                 &use_user_constraints,
                                 &weights_type);

   ierr = MPI_Barrier(comm_all);
   // TODO: call time_end(t_pc_setup)

   if (myid == 0) {
      printf("Preconditioner set-up done.\n");
   }



   // call Krylov method
   if (myid == 0) {
      printf("Calling Krylov method ...\n");
   }
   ierr = MPI_Barrier(comm_all);
   // TODO: call time_start
   // call with setting of iterative properties

   int fortran_comm =  MPI_Comm_c2f(comm_all);

   bddcml_solve_c(&fortran_comm, &krylov_method, &tol,&maxit,&ndecrmax, &recycling_int, &max_number_of_stored_vectors,
                     &num_iter, &converged_reason, &condition_number);
   ierr = MPI_Barrier(comm_all);
   // TODO: call time_end(t_krylov)
   if (myid == 0) {
      printf("Krylov method done.\n");
   }
   
   if (myid == 0) {
      printf(" Output of PCG: ==============\n");
      printf(" Number of iterations: %d\n", num_iter);
      printf(" Convergence reason:   %d\n", converged_reason);
      if ( condition_number >= 0. ) {
         printf(" Condition number: %lf\n", condition_number);
      }
      printf(" =============================\n");
   }

   normRn2_loc  = 0.;
   normL2_loc   = 0.;
   normLinf_loc = 0.;
   if (general_params.export_solution) {
      // make sure the directory exists
      if (myid == 0) {
         sprintf(command, "mkdir -p %s", general_params.output_file_prefix);
         system(command);
      }
      ierr = MPI_Barrier(comm_all);
   }
   for(isub = sub2proc[myid]; isub < sub2proc[myid+1]; isub++) {
      // download local solution
      nnods    = pow(num_el_per_sub_edge+1, ndim);
      ndofs    = nnods;
      lsols = ndofs;
      sols = (real*) malloc(lsols * sizeof(real));

      bddcml_download_local_solution_c(&isub, sols, &lsols);

      // compute norm of local solution
      if (nsub > 1) {
         if (general_params.just_direct_solve_int == 0) {
            bddcml_dotprod_subdomain_c( &isub, sols, &lsols, sols, &lsols, &normRn2_sub );
         }
         else {
            // cannot determine norm for solution by a direct solver
            normRn2_sub = 0.;
         }
      }
      else {
         //TODO: call blas... normRn2_sub = dot_product(sols,sols);
         normRn2_sub = 0.;
         for(i = 0; i < lsols; i++){
            normRn2_sub += sols[i] * sols[i];
         }
      }
         
      normRn2_loc = normRn2_loc + normRn2_sub;

      // re-create mesh for subdomain
      nelems   = pow(num_el_per_sub_edge, ndim);
      nnods    = pow(num_el_per_sub_edge+1, ndim);
      ndofs    = nnods;
      linets   = nelems * 8;
      lnnets   = nelems;
      lnndfs   = nnods;
      lisegns  = nelems;
      lisngns  = nnods;
      lisvgvns = ndofs;
      
      inets = (int*) malloc(linets * sizeof(int));
      nnets = (int*) malloc(lnnets * sizeof(int));
      nndfs = (int*) malloc(lnndfs * sizeof(int));
      isegns = (int*) malloc(lisegns * sizeof(int));
      isngns = (int*) malloc(lisngns * sizeof(int));
      isvgvns = (int*) malloc(lisvgvns * sizeof(int));

      lxyzs1   = nnods;
      lxyzs2   = ndim;
      xyzs = (real*) malloc(lxyzs1 * lxyzs2 * sizeof(real));
      lifixs   = nnods;
      lfixvs   = nnods;
      ifixs = (int*) malloc(lifixs * sizeof(int));
      fixvs = (real*) malloc(lfixvs * sizeof(real));

      prepare_subdomain_data(isub, num_sub_per_cube_edge, num_el_per_sub_edge, hsize,
                                    inets, linets, nnets, lnnets, nndfs, lnndfs, 
                                    isegns, lisegns, isngns, lisngns, isvgvns, lisvgvns,
                                    xyzs,lxyzs1,lxyzs2, 
                                    ifixs,lifixs, fixvs,lfixvs);

      // compute L_2 norm of the solution
      normL2_sub = 0.;
      indinets = 0;
      for (ie = 0; ie < nelems; ie++) {
         // number of nodes on element
         nne = nnets[ie];
         
         // sum(sols(inets(indinets+1:indinets+nne)))
         real sum_help = 0.;
         int sum_idx;
         for (sum_idx = indinets; sum_idx < indinets + nne; sum_idx++) { // TODO: nebo +1???
            sum_help += sols[inets[sum_idx]];
         }
         
         normL2_sub = normL2_sub + el_vol * sum_help / nne;
         indinets = indinets + nne;
      }
      normL2_loc = normL2_loc + normL2_sub;

      // find maximum of the solution
      int max_idx;
      for(max_idx = 0; max_idx < lsols; max_idx++) {
         if(sols[max_idx] > normLinf_loc) {
            normLinf_loc = sols[max_idx];
         }
      }
//      normLinf_loc = max(normLinf_loc, maxval(sols));

      if (general_params.export_solution) {
         //TODO: 
//          export_vtu_file_with_solution(output_file_prefix, 
//                                              isub, nelems, nnods, 
//                                              inets,linets, nnets,lnnets, 
//                                              xyzs,lxyzs1,lxyzs2, 
//                                              sols,lsols);
      }
      free(inets);
      free(nnets);
      free(nndfs);
      free(isegns);
      free(isngns);
      free(isvgvns);
      free(xyzs);
      free(ifixs);
      free(fixvs);
      free(sols);
   }
   if (general_params.export_solution) {
      // export the umbrella PVD file
      if (myid == 0) {
         //TODO: paraview_export_pvd_file("poisson_solution", nsub);
      }
   }

   if (myid == 0) {
      printf("Finalizing BDDCML ...\n");
   }
   
   bddcml_finalize_c();      
   if (myid == 0) {
      printf("Finalizing BDDCML done.\n");
   }

   // find global R^n norm of solution
   ierr = MPI_Allreduce(&normRn2_loc, &normRn2, 1, MPI_REAL, MPI_SUM, comm_all);
   normRn_sol = sqrt( normRn2 );
   if (myid == 0) {
   }

   // find global L_2 norm of solution
   ierr = MPI_Allreduce(&normL2_loc, &normL2_sol, 1, MPI_REAL, MPI_SUM, comm_all);
   if (myid == 0) {
   }

   // find global L_inf norm of solution
   ierr = MPI_Allreduce(&normLinf_loc, &normLinf_sol, 1, MPI_REAL, MPI_MAX, comm_all);
   if (myid == 0) {
   }
   if (myid == 0) {
      printf(" Solution properties========\n");
      printf(" L_2 norm:   %g\n", normL2_sol);
      printf(" L_inf norm: %g\n", normLinf_sol);
      if (general_params.just_direct_solve_int == 0) {
         printf(" R^n norm:   %g\n", normRn_sol);
      }
      printf(" ===========================\n");
   }

   ierr = MPI_Barrier(comm_all);
   // TODO: call time_end(t_total);

   // TODO:
//       ! Information about times
//       if (myid.eq.0) then
//          write(*,'(a)')         ' Profiling infirmation==========='
//          write(*,'(a)')         ' TIMES OF RUN OF ROUTINES:'
//          write(*,'(a,f11.3,a)') '  loading           ',t_load, ' s'
//          write(*,'(a,f11.3,a)') '  pc_setup          ',t_pc_setup, ' s'
//          write(*,'(a,f11.3,a)') '  Krylov method     ',t_krylov, ' s'
//          write(*,'(a)')         '  _______________________________'
//          write(*,'(a,f11.3,a)') '  total             ',t_total,    ' s'
//          write(*,'(a)')         ' ================================'
//       end if

   // MPI finalization
//***************************************************************PARALLEL
   ierr = MPI_Finalize();
//***************************************************************PARALLEL

//     free(nsublev);
//     free(sub2proc);
   
   return(0);
}
