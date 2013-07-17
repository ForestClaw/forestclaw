# Makefile for sphere example

TARGET = swirl
GEO = cart

# -----------------------------------------------------------
FC = gfortran
FFLAGS = -O0 -g -Wno-unused-dummy-argument -Wno-unused-variable

# -----------------------------------------------------------
# Define local objects here
exampleobjs = \
	ghost_twopatch.o \
	tag4refinement.o \
	simple_user.o \
	qinit.o \
	b4step2.o \
	setaux.o \
	setprob.o

# Solver routines (several solvers can be added here)
include Makefile.manyclaw.inc
solverobjs = $(manyclaw_objs)
solverincludes = $(manyclaw_includes)

# Mapping routines
geoobjs = $(FORESTCLAW)/src/mappings/sphere/set_maptype_sphere.o

OBJECTS = $(TARGET).o $(exampleobjs) $(solverobjs) $(geoobjs)

# Specify location of gfortran library
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
