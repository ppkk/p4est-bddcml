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

include_directories("${PARMETIS_ROOT}/include")
link_directories("${PARMETIS_ROOT}/build/Linux-x86_64/libparmetis")
set(LIBS ${LIBS}  parmetis)

include_directories("${PARMETIS_ROOT}/metis/include")
link_directories("${PARMETIS_ROOT}/build/Linux-x86_64/libmetis")
set(LIBS ${LIBS}  metis)

include_directories("${P4EST_ROOT}/include")
link_directories("${P4EST_ROOT}/local/lib")
set(LIBS ${LIBS} sc p4est )

set(SOURCES_DIM_INDEP
src/arrays.cc
src/bddcml_structs.cc
)

set(SOURCES_2D
src/p4est_common.cc
src/mesh.cc
src/femspace.cc
src/assemble.cc
)

set(SOURCES_3D
3d_wrappers/p8est_common.cc
3d_wrappers/p8est_mesh.cc
3d_wrappers/p8est_femspace.cc
3d_wrappers/p8est_assemble.cc
)

set(HEADERS
include/definitions.h
include/arrays.h
include/assemble.h
include/p4est_common.h
include/bddcml_structs.h
include/mesh.h
include/femspace.h
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

add_library(mylib_2D ${HEADERS} ${SOURCES_DIM_INDEP} ${SOURCES_2D})
add_library(mylib_3D ${HEADERS} ${SOURCES_DIM_INDEP} ${SOURCES_3D})

add_executable(p4est_example drivers/p4est_example.cc)
add_executable(p8est_example 3d_wrappers/p8est_example.cc)
add_executable(bddcml_example drivers/bddcml_cube_example.cc)
add_executable(p4est_bddcml drivers/p4est_bddcml)
add_executable(p8est_bddcml 3d_wrappers/p8est_bddcml.cc)
add_executable(metis_bddcml drivers/metis_bddcml.cc)
add_executable(metis8_bddcml 3d_wrappers/metis8_bddcml.cc)
add_executable(p4est_experiment drivers/p4est_experiment.cc)
target_link_libraries(p4est_example ${BLOPEX_OBJECTS} mylib_2D ${LIBS} )
target_link_libraries(p8est_example ${BLOPEX_OBJECTS} mylib_3D ${LIBS} )
target_link_libraries(bddcml_example ${BLOPEX_OBJECTS} mylib_2D ${LIBS})
target_link_libraries(p4est_bddcml ${BLOPEX_OBJECTS} mylib_2D ${LIBS})
target_link_libraries(p8est_bddcml ${BLOPEX_OBJECTS} mylib_3D ${LIBS})
target_link_libraries(metis_bddcml ${BLOPEX_OBJECTS} mylib_2D ${LIBS})
target_link_libraries(metis8_bddcml ${BLOPEX_OBJECTS} mylib_3D ${LIBS})
target_link_libraries(p4est_experiment ${BLOPEX_OBJECTS} mylib_2D ${LIBS})
