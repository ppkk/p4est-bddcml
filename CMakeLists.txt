cmake_minimum_required (VERSION 2.6)
project(p4est-bddcml)
include(CMake.vars)

#debug, remove...
set (CMAKE_C_FLAGS " -g --std=gnu11 -pedantic-errors -Wall ${CMAKE_C_FLAGS}")
set (CMAKE_CXX_FLAGS " -g -Wall ${CMAKE_CXX_FLAGS}")
set (CMAKE_CXX_COMPILER mpic++)

set_source_files_properties(src/metis_bddcml.c  src/p4est_common.c src/p8est_common.c src/bddcml_structs.c  src/helpers.c src/p4est_bddcml_interaction.c src/p8est_bddcml_interaction.c src/p8est_bddcml.c src/p4est_bddcml.c src/bddcml_cube_example.c src/bddcml_cube_example_subdomain.c src/p4est_example.c src/p8est_example.c PROPERTIES LANGUAGE CXX )

include_directories("${BDDCML_ROOT}/src")
link_directories("${BDDCML_ROOT}/lib")
set(LIBS ${LIBS} bddcml)

include_directories("${MUMPS_ROOT}/include")
link_directories("${MUMPS_ROOT}/lib")
set(LIBS ${LIBS}  dmumps mumps_common pord)

include_directories("${PARMETIS_ROOT}/include")
link_directories("${PARMETIS_ROOT}/build/Linux-x86_64/libparmetis")
set(LIBS ${LIBS}  parmetis)

include_directories("${PARMETIS_ROOT}/metis/include")
link_directories("${PARMETIS_ROOT}/build/Linux-x86_64/libmetis")
set(LIBS ${LIBS}  metis)

include_directories("${P4EST_ROOT}/include")
link_directories("${P4EST_ROOT}/local/lib")
set(LIBS ${LIBS} sc p4est )


# BLOPEX objects
include_directories("${BLOPEX_ROOT}/blopex_serial_double/multivector")
include_directories("${BLOPEX_ROOT}/blopex_serial_double/pcg_multi")
include_directories("${BLOPEX_ROOT}/blopex_serial_double/matmultivec")
include_directories("${BLOPEX_ROOT}//blopex_abstract/include")
link_directories("${BLOPEX_ROOT}/blopex_abstract/lib")
set(LIBS ${LIBS}  BLOPEX)

set(BLOPEX_OBJECTS ${BLOPEX_OBJECTS} "${BLOPEX_ROOT}/blopex_serial_double/multivector/multi_vector.o")
set(BLOPEX_OBJECTS ${BLOPEX_OBJECTS} "${BLOPEX_ROOT}/blopex_serial_double/pcg_multi/pcg_multi.o")
set(BLOPEX_OBJECTS ${BLOPEX_OBJECTS} "${BLOPEX_ROOT}/blopex_serial_double/matmultivec/matmultivec.o")

set(LIBS ${LIBS}  scalapack-openmpi m gfortran blas lapack blacs-openmpi blacsCinit-openmpi blacsF77init-openmpi blacs-openmpi)

include_directories("${MPI_ROOT}/include")
set(LIBS ${LIBS}  mpi pthread mpi_f77)


include_directories("include")


add_executable(p4est_example src/p4est_example.c  src/p4est_common.c src/helpers.c)
add_executable(p8est_example src/p8est_example.c src/p8est_common.c src/helpers.c)
add_executable(bddcml_example ${BLOPEX_OBJECTS} src/bddcml_cube_example.c src/bddcml_structs.c src/bddcml_cube_example_subdomain.c src/helpers.c)
add_executable(p4est_bddcml ${BLOPEX_OBJECTS} src/p4est_bddcml.c  src/p4est_common.c src/bddcml_structs.c  src/helpers.c src/p4est_bddcml_interaction.c)
add_executable(p8est_bddcml ${BLOPEX_OBJECTS} src/p8est_bddcml.c  src/p8est_common.c src/bddcml_structs.c  src/helpers.c src/p8est_bddcml_interaction.c)
add_executable(metis_bddcml ${BLOPEX_OBJECTS} src/metis_bddcml.c  src/p4est_common.c src/bddcml_structs.c  src/helpers.c src/p4est_bddcml_interaction.c)
target_link_libraries(p4est_example ${LIBS})
target_link_libraries(p8est_example ${LIBS})
target_link_libraries(bddcml_example ${LIBS})
target_link_libraries(p4est_bddcml ${LIBS})
target_link_libraries(p8est_bddcml ${LIBS})
target_link_libraries(metis_bddcml ${LIBS})
