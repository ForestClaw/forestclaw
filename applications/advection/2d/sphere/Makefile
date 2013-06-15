# Makefile for sphere example

TARGET = sphere
GEO = sphere

# -----------------------------------------------------------
FC = gfortran
FFLAGS = -O0 -g -Wno-unused-dummy-argument -Wno-unused-variable

# -----------------------------------------------------------
# Define local objects here
exampleobjs = \
	ghost_twopatch.o \
	velocity_wind.o \
	tag4refinement.o \
	mapc2m.o \
	rpn2noncons.o \
	rpt2noncons.o \
	qinit_manifold.o \
	b4step2_manifold.o \
	setaux_manifold.o \
	setprob.o

# Solver routines (several solvers can be added here)
include Makefile.waveprop.inc
solverobjs = $(waveprop_objs)
solverincludes = $(waveprop_includes)

# Mapping routines
geoobjs = $(FORESTCLAW)/src/mappings/sphere/set_maptype_sphere.o

OBJECTS = $(TARGET).o $(exampleobjs) $(solverobjs) $(geoobjs)

# Specify location of fortran library
FORT_LIB = -L/usr/local/gfortran/lib -lgfortran

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
