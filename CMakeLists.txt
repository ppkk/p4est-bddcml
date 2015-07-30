cmake_minimum_required (VERSION 2.6)
include(CMake.vars)

set (CMAKE_C_FLAGS "--std=gnu11 -pedantic-errors ${CMAKE_C_FLAGS}")

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


add_executable(p4est-bddcml src/p4est-bddcml.c)
add_executable(p8est-bddcml src/p8est-bddcml.c)
add_executable(bddcml-example ${BLOPEX_OBJECTS} src/bddcml.c src/bddcml_structs.c src/bddcml_cube_example.c src/helpers.c)
target_link_libraries(p4est-bddcml ${LIBS})
target_link_libraries(p8est-bddcml ${LIBS})
target_link_libraries(bddcml-example ${LIBS})
