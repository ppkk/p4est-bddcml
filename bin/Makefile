include ../make.inc

INC = -I../include \
  $(P4B_INCP4EST) \
  $(P4B_INCBDDCML)
  
LIBS = \
  $(P4B_LIBBDDCML) $(P4B_LIBP4EST) $(P4B_LIBMUMPS) $(P4B_LIBPARMETIS) \
  $(P4B_LIBMETIS) $(P4B_LIBSCALAPACK) $(P4B_LIBBLACS) \
  $(P4B_LIBPORD) $(P4B_LIBLAPACK) $(P4B_LIBBLAS) \
  $(P4B_BLOPEXLINK) $(P4B_LIBOTHERS) 
  
DRIVERS = p4est_bddcml adapt

all: $(DRIVERS)
#all: p4est_bddcml.o $(DRIVERS)

%.o : ../drivers/%.cc
	$(P4B_MPCC) $(P4B_MPCFLAGS) $(INC) -c $< -o $@

p4est_bddcml : p4est_bddcml.o $(LIBP4B)
	$(P4B_MPCC) $(P4B_MPCFLAGS) -o $@ $^  -L../lib/ -lp4b $(LIBS)

adapt : uniform_conv.o $(LIBP4B)
	$(P4B_MPCC) $(P4B_MPCFLAGS) -o $@ $^  -L../lib/ -lp4b $(LIBS)

# %.x : %.o
# 	$(P4B_MPCC) $(P4B_MPCFLAGS) -o $@ $^  -L../lib/ -lp4b $(LIBS)

clean: 
	-rm -f *.o $(DRIVERS)
