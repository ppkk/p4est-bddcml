cmake_minimum_required (VERSION 2.6)
project(p4est-bddcml)
include(CMake.vars)

#debug, remove...
set (CMAKE_C_FLAGS " -g --std=gnu11 -pedantic-errors -Wall ${CMAKE_C_FLAGS}")
set (CMAKE_CXX_FLAGS " -g --std=c++11 -Wall -Wno-literal-suffix ${CMAKE_CXX_FLAGS}")
set (CMAKE_CXX_COMPILER mpic++)

include_directories("${BDDCML_ROOT}/src")
link_directories("${BDDCML_ROOT}/lib")
set(LIBS ${LIBS} bddcml)

include_directories("${MUMPS_ROOT}/include")
link_directories("${MUMPS_ROOT}/lib")
set(LIBS ${LIBS}  dmumps mumps_common pord)

include_directories("${PARMETIS_INCLUDE}")
link_directories("${PARMETIS_LIB}")
set(LIBS ${LIBS}  parmetis)

include_directories("${METIS_INCLUDE}")
link_directories("${METIS_LIB}")
set(LIBS ${LIBS}  metis)

include_directories("${P4EST_ROOT}/local/include")
link_directories("${P4EST_ROOT}/local/lib")
set(LIBS ${LIBS} sc p4est )

set(SOURCES
src/bddcml/arrays.cc
src/bddcml/bddcml_structs.cc
src/bddcml/bddcml_mesh.cc
src/bddcml/bddcml_femspace.cc
src/p4est/my_p4est_interface.cc
src/p4est/my_p4est_implementation.cc
src/p4est/my_p8est_implementation.cc
src/definitions.cc
src/quadrature.cc
src/gauss_quadrature.cc
src/level_set.cc
src/element.cc
src/integration_cell.cc
src/assemble.cc
src/geometry_mesh.cc
src/shapefun.cc
src/local_matrix.cc
src/local_solution.cc
src/vtk_output.cc
)

set(HEADERS
include/bddcml/arrays.h
include/bddcml/bddcml_structs.h
include/bddcml/bddcml_mesh.h
include/bddcml/bddcml_femspace.h
include/p4est/my_p4est_interface.h
include/p4est/my_p4est_implementation.h
include/definitions.h
include/quadrature.h
include/level_set.h
include/element.h
include/integration_cell.h
include/assemble.h
include/geometry_mesh.h
include/shapefun.h
include/local_matrix.h
include/local_solution.h
include/vtk_output.h
)

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

set(LIBS ${LIBS}
scalapack-openmpi
m
gfortran
blas
lapack
blacs-openmpi
blacsCinit-openmpi
blacsF77init-openmpi
blacs-openmpi
)

include_directories("${MPI_ROOT}/include")
set(LIBS ${LIBS}
mpi
pthread mpi_f77
)


include_directories("include")

add_library(mylib ${HEADERS} ${SOURCES})

add_executable(bddcml_example drivers/bddcml_cube_example.cc)
add_executable(p4est_bddcml drivers/p4est_bddcml)
#add_executable(metis_bddcml drivers/metis_bddcml.cc)
#add_executable(p4est_experiment drivers/p4est_experiment.cc)
#add_executable(p4est_migration drivers/p4est_migration.cc)
add_executable(test_quadrature drivers/test_quadrature.cc)

target_link_libraries(bddcml_example ${BLOPEX_OBJECTS} mylib ${LIBS})
target_link_libraries(p4est_bddcml ${BLOPEX_OBJECTS} mylib ${LIBS})
#target_link_libraries(metis_bddcml ${BLOPEX_OBJECTS} mylib ${LIBS})
#target_link_libraries(p4est_experiment ${BLOPEX_OBJECTS} mylib ${LIBS})
#target_link_libraries(p4est_migration ${BLOPEX_OBJECTS} mylib ${LIBS})
target_link_libraries(test_quadrature mylib)
