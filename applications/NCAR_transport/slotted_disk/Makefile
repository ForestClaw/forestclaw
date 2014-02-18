# Makefile for slotted_disk example

TARGET = slotted_disk
GEO = sphere

include $(FORESTCLAW)/Make.defs.local

# -----------------------------------------------------------
# FC = gfortran  # Now defined in Make.defs.local
FFLAGS = $(LOCAL_FFLAGS)

# -----------------------------------------------------------
# Define local objects here

BENCH = ../all

exampleobjs = \
	slotted_disk_user.o \
	slotted_disk_output.o \
	tag4refinement.o \
	${BENCH}/setprob_transport.o \
	${BENCH}/qinit_transport.o \
	${BENCH}/slotted_disk.o \
	${BENCH}/gaussian.o \
	${BENCH}/cosine_bell.o \
	${BENCH}/velocity_wind.o \
	${BENCH}/setaux_transport.o \
	${BENCH}/mapc2m.o \
	${BENCH}/rpn2noncons.o \
	${BENCH}/rpt2noncons.o \
	${BENCH}/exact_solution.o \
	${BENCH}/b4step2_transport.o \
	${BENCH}/ghost_twopatch.o

# Solver routines (several solvers can be added here)
include Makefile.waveprop.inc
solverobjs = $(waveprop_objs)
solverincludes = $(waveprop_includes)

# Mapping routines
geoobjs = $(FORESTCLAW)/src/mappings/sphere/set_maptype_${GEO}.o

OBJECTS = $(TARGET).o $(exampleobjs) $(solverobjs) $(geoobjs)

# -----------------------------------------------------------
# The rest of this is more or less default.
# -----------------------------------------------------------

include $(P4EST_DIR)/etc/Makefile.p4est.mk

# make sure that the first make target is the executable
default: $(TARGET)

.cpp.o:
	$(P4EST_CXX) $(P4EST_CFLAGS) $(P4EST_CPPFLAGS) \
	-I$(FORESTCLAW)/src -I. ${solverincludes} -c $< -o $@

.c.o:
	$(P4EST_CXX) $(P4EST_CFLAGS) $(P4EST_CPPFLAGS) \
	-I$(FORESTCLAW)/src -I. ${solverincludes} -c $< -o $@

.f.o:
	$(FC) -Wall $(FFLAGS) -c $< ${solverincludes} -o $@

.PHONY: clean realclean default $(TARGET)

$(TARGET): $(OBJECTS) $(FORESTCLAW)/src/libforestclaw.a
	$(P4EST_CXX) $(P4EST_CFLAGS) $^ $(P4EST_LDFLAGS) $(FORT_LIB) \
	$(P4EST_LIBS) -o $@

clean:
	rm -f *.o $(TARGET) $(exampleobjs)

realclean:
	rm -f *.o $(TARGET) $(OBJECTS) fort.*
