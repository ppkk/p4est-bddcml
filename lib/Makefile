include ../make.inc

INC = -I../include \
  $(P4B_INCP4EST) \
  $(P4B_INCBDDCML)

OBJ = \
  analytical_solution.o \
  assemble.o \
  definitions.o \
  element.o \
  gauss_quadrature.o \
  geometry_mesh.o \
  integral.o \
  integration_cell.o \
  level_set.o \
  local_matrix.o \
  local_solution.o \
  quadrature.o \
  shapefun.o \
  vtk_output.o \
  arrays.o \
  bddcml_mesh.o \
  bddcml_femspace.o \
  bddcml_structs.o \
  bddcml_solver.o \
  my_p4est_implementation.o \
  my_p4est_interface.o \
  my_p8est_implementation.o

all: modules p4blib

modules: $(OBJ)

p4blib: $(LIBP4B)

%.o : ../src/%.cc
	$(P4B_MPCC) $(P4B_MPCFLAGS) $(INC) -c $< -o $@

%.o : ../src/bddcml/%.cc
	$(P4B_MPCC) $(P4B_MPCFLAGS) $(INC) -c $< -o $@

%.o : ../src/p4est/%.cc
	$(P4B_MPCC) $(P4B_MPCFLAGS) $(INC) -c $< -o $@

$(LIBP4B): $(OBJ)
	$(P4B_AR) $@ $^
	$(P4B_RANLIB) $@

clean:
	-rm *.o $(LIBP4B)

