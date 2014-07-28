# Makefile for sphere example

TARGET = sphere
GEO = sphere

include $(FORESTCLAW)/Make.defs.local

# -----------------------------------------------------------
FFLAGS = $(LOCAL_FFLAGS)

# -----------------------------------------------------------
# Define local objects here
exampleobjs = \
	sphere_user.o \
	setprob.o \
	ghost_twopatch.o \
	output_metric.o \
	tag4refinement.o \
	mapc2m.o

# Solver routines (several solvers can be added here)
solverobjs =
solverincludes =

# Mapping routines
geoobjs = $(FORESTCLAW)/src/mappings/set_maptype.o

OBJECTS = $(TARGET).o $(exampleobjs) $(geoobjs)

# Specify location of fortran library
FORT_INC = $(FORESTCLAW)/src/fortran_source2d

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
	$(FC) $(FFLAGS) -I${FORT_INC} -c $<  -o $@

.PHONY: clean realclean default $(TARGET)

$(TARGET): $(OBJECTS) $(FORESTCLAW)/src/libforestclaw.a
	$(P4EST_CXX) $(P4EST_CFLAGS) $^ $(P4EST_LDFLAGS) $(FORT_LIB) \
	$(P4EST_LIBS) -o $@

clean:
	rm -f *.o $(TARGET) $(exampleobjs)

realclean:
	rm -f *.o $(TARGET) $(OBJECTS) fort.*
