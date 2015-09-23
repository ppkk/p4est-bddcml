#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <assert.h>
#include "bddcml_structs.h"
#include "helpers.h"
#include "bddcml_cube_example_subdomain.h"


// *******************
// PROBLEM PARAMETERS:
// *******************

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


   const char routine_name[] = "POISSON_ON_CUBE";

// input read from command line
   int num_el_per_sub_edge, num_sub_per_cube_edge; // basic properties of the cubes-in-cubes problem

   int num_el_per_cube_edge;
   real hsize, el_vol;                              // elements size and volume

// parallel variables
   MPI_Comm comm_all;
   int myid, nproc, ierr;

   int nsub;  // number of subdomains on the first level 
   
// scaled element matrix
   real element_matrix[NDOF_PER_ELEMENT * NDOF_PER_ELEMENT];

   SparseMatrix matrix;
   RealArray rhss;
   RealArray sols;


// user constraints - not really used here
   Real2DArray user_constraints;

// data for elements - not really used here
   Real2DArray element_data;

// data for dofs - not really used here
   RealArray dof_data;

// data about resulting convergence
   BddcmlConvergenceInfo convergence_info;

   real normRn_sol, normRn2, normRn2_loc, normRn2_sub;
   real normL2_sol, normL2_loc, normL2_sub;
   real normLinf_sol, normLinf_loc;

// time variables
   real t_total, t_load, t_pc_setup, t_krylov;

// small variables - indices, etc.
   int element_offset, num_nodes_of_elem, idof, jdof, lelm;
   int ie, i, isub, j;
   char aux[32];
   char command[256];
   int idx;


int main(int argc, char **argv)
{
   BddcmlGeneralParams general_params;
   set_implicit_general_params(&general_params);

   BddcmlKrylovParams krylov_params;
   set_implicit_krylov_params(&krylov_params);

   BddcmlPreconditionerParams preconditioner_params;
   set_implicit_preconditioner_params(&preconditioner_params);

   BddcmlLevelInfo level_info;

   BddcmlDimensions global_dims, subdomain_dims;
   init_dimmensions(&global_dims, 3, LAPLACE);
   init_dimmensions(&subdomain_dims, 3, LAPLACE);

   BddcmlMesh mesh;
   BddcmlFemSpace femsp;

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
   el_vol = pow(hsize, global_dims.n_problem_dims);
   // total number of elements
   global_dims.n_elems = pow(num_el_per_cube_edge, global_dims.n_problem_dims);
   // total number of subdomains
   nsub  = pow(num_sub_per_cube_edge, global_dims.n_problem_dims);
   // total number of nodes
   global_dims.n_nodes  = pow(num_el_per_cube_edge + 1, global_dims.n_problem_dims);
   // total number of degrees of freedom - equal to number of nodes for scalar problem
   global_dims.n_dofs  = global_dims.n_nodes * global_dims.n_node_dofs;

   init_levels(nsub, &level_info);

   print_basic_properties(&global_dims, nsub, &level_info, &krylov_params);

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
      subdomain_dims.n_elems   = pow(num_el_per_sub_edge, global_dims.n_problem_dims);
      subdomain_dims.n_nodes    = pow(num_el_per_sub_edge+1, global_dims.n_problem_dims);
      subdomain_dims.n_dofs    = subdomain_dims.n_nodes * subdomain_dims.n_node_dofs;

      init_mesh(&subdomain_dims, &mesh);
      init_fem_space(&subdomain_dims, &femsp);

      // create subdomain mesh and boundary conditions
      prepare_subdomain_data(isub, num_sub_per_cube_edge, num_el_per_sub_edge, hsize,
                                    &mesh, &femsp);

      // create local right hand side
      allocate_real_array(subdomain_dims.n_dofs, &rhss);
      for(idx = 0; idx < rhss.len; idx++) {
         if(femsp.fixs_code.val[idx] > 0) {
            // rhs is zero on boundary
            rhss.val[idx] = 0.;
         }
         else {
            rhss.val[idx] = 1. * el_vol;
         }
      }
      int is_rhs_complete = 1;

      // create local initial solution
      allocate_real_array(subdomain_dims.n_dofs, &sols);
      zero_real_array(&sols);

      // create local subdomain matrix for each subdomain
      // full element matrix scaled based on element size
      for(idx = 0; idx < NDOF_PER_ELEMENT * NDOF_PER_ELEMENT; idx++) {
         element_matrix[idx] = hsize * element_matrix_ref[idx];
      }

      // how much space the upper triangle of the element matrix occupies
      lelm = NDOF_PER_ELEMENT * (NDOF_PER_ELEMENT + 1) / 2;
      // space for all upper triangles of element matrics

      allocate_sparse_matrix(subdomain_dims.n_elems*lelm, SPD, &matrix);
      zero_matrix(&matrix);

      // copy the upper triangle of the element matrix to the sparse triplet
      element_offset = 0;
      for(ie = 0; ie < subdomain_dims.n_elems; ie++) {
         num_nodes_of_elem = mesh.num_nodes_of_elem.val[ie];
         for(j = 0; j < NDOF_PER_ELEMENT; j++) {
            jdof = mesh.elem_node_indices.val[element_offset + j];
            // this has been changed, we send all entries to the function, and it discards those where idof > jdof
            // it is done because of the hanging nodes
            for(i = 0; i < NDOF_PER_ELEMENT /*<= j*/; i++) {
               idof = mesh.elem_node_indices.val[element_offset + i];

               add_matrix_entry(&matrix, idof, jdof, element_matrix[j * NDOF_PER_ELEMENT + i]);
            }
         }
         element_offset += num_nodes_of_elem;
      }

      // prepare user constraints - not really used here
      allocate_real_2D_array(0, 0, &user_constraints);

      // prepare element data - not really used
      allocate_real_2D_array(0, 0, &element_data);

      // prepare dof data - not really used
      allocate_real_array(0, &dof_data);

//       if(myid == 0)
//       {
//          print_array(linets, inets, "inets");
//          print_array(lnnets, nnets, "nnets");
//          print_array(la, i_sparse, "i");
//          print_array(la, j_sparse, "j");
//          print_f_array(la, a_sparse, "a");
//       }

      bddcml_upload_subdomain_data(&global_dims, &subdomain_dims,
                                        isub, &mesh, &femsp,
                                        &rhss, is_rhs_complete, &sols, &matrix,
                                        &user_constraints, &element_data,
                                        &dof_data, &preconditioner_params);

      free_mesh(&mesh);
      free_fem_space(&femsp);

      free_real_array(&rhss);
      free_real_array(&sols);
      free_sparse_matrix(&matrix);
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
   bddcml_setup_preconditioner(matrix.type, &preconditioner_params);

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

   bddcml_solve(&krylov_params, &convergence_info, comm_all);
   ierr = MPI_Barrier(comm_all);
   // TODO: call time_end(t_krylov)
   if (myid == 0) {
      printf("Krylov method done.\n");
   }
   
   if (myid == 0) {
      printf(" Output of PCG: ==============\n");
      printf(" Number of iterations: %d\n", convergence_info.num_iter);
      printf(" Convergence reason:   %d\n", convergence_info.converged_reason);
      if ( convergence_info.condition_number >= 0. ) {
         printf(" Condition number: %lf\n", convergence_info.condition_number);
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
      subdomain_dims.n_nodes    = pow(num_el_per_sub_edge+1, global_dims.n_problem_dims);
      subdomain_dims.n_dofs    = subdomain_dims.n_nodes * subdomain_dims.n_node_dofs;

      allocate_real_array(subdomain_dims.n_dofs, &sols);

      bddcml_download_local_solution(isub, &sols);

      // compute norm of local solution
      if (nsub > 1) {
         if (general_params.just_direct_solve_int == 0) {
            bddcml_dotprod_subdomain( isub, &sols, &sols, &normRn2_sub );
         }
         else {
            // cannot determine norm for solution by a direct solver
            normRn2_sub = 0.;
         }
      }
      else {
         //TODO: call blas... normRn2_sub = dot_product(sols,sols);
         normRn2_sub = 0.;
         for(i = 0; i < sols.len; i++){
            normRn2_sub += sols.val[i] * sols.val[i];
         }
      }
         
      normRn2_loc = normRn2_loc + normRn2_sub;

      // re-create mesh for subdomain
      subdomain_dims.n_elems   = pow(num_el_per_sub_edge, global_dims.n_problem_dims);
      subdomain_dims.n_nodes    = pow(num_el_per_sub_edge+1, global_dims.n_problem_dims);
      subdomain_dims.n_dofs    = subdomain_dims.n_nodes * subdomain_dims.n_node_dofs;
      
      init_mesh(&subdomain_dims, &mesh);
      init_fem_space(&subdomain_dims, &femsp);

      prepare_subdomain_data(isub, num_sub_per_cube_edge, num_el_per_sub_edge, hsize,
                                    &mesh, &femsp);

      // compute L_2 norm of the solution
      normL2_sub = 0.;
      element_offset = 0;
      for (ie = 0; ie < subdomain_dims.n_elems; ie++) {
         // number of nodes on element
         num_nodes_of_elem = mesh.num_nodes_of_elem.val[ie];
         
         // sum(sols(inets(indinets+1:indinets+nne)))
         real sum_help = 0.;
         int sum_idx;
         for (sum_idx = element_offset; sum_idx < element_offset + num_nodes_of_elem; sum_idx++) { // TODO: nebo +1???
            sum_help += sols.val[mesh.elem_node_indices.val[sum_idx]];
         }
         
         normL2_sub = normL2_sub + el_vol * sum_help / num_nodes_of_elem;
         element_offset = element_offset + num_nodes_of_elem;
      }
      normL2_loc = normL2_loc + normL2_sub;

      // find maximum of the solution
      int max_idx;
      for(max_idx = 0; max_idx < sols.len; max_idx++) {
         if(sols.val[max_idx] > normLinf_loc) {
            normLinf_loc = sols.val[max_idx];
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

      free_mesh(&mesh);
      free_fem_space(&femsp);
      free_real_array(&sols);
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
   
   bddcml_finalize();
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
   
   assert(get_num_allocations() == 0);
   return(0);
}
