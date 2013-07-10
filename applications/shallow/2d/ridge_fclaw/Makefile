# Makefile for sphere example

TARGET = ridge
GEO = sphere

include $(FORESTCLAW)/Make.defs.local

# -----------------------------------------------------------
FC = gfortran
FFLAGS = -O0 -g -Wno-unused-dummy-argument -Wno-unused-variable

# -----------------------------------------------------------
# Define local objects here
exampleobjs = \
	ridge_user.o \
	ghost_twopatch.o \
	tag4refinement.o \
	mapc2m.o \
	rpn2swsphere_fw.o \
	rpt2swq.o \
	qinit.o \
	interpolate_geo.o \
	b4step2.o \
	setgauges.o \
	settsunami_geo.o \
	setaux_ridge.o \
	src2.o \
	setprob.o \
	output_geo.o \
	bmount.o

# Solver routines (several solvers can be added here)
include Makefile.waveprop.inc
solverobjs = $(waveprop_objs)
solverincludes = $(waveprop_includes)

# Mapping routines
geoobjs = $(FORESTCLAW)/src/mappings/sphere/set_maptype_sphere.o

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
