# Makefile for sphere example

TARGET = test_parms
GEO = cart

# -----------------------------------------------------------
FC = gfortran
FFLAGS = -O0 -g -Wno-unused-dummy-argument -Wno-unused-variable

# -----------------------------------------------------------
# Define local objects here
exampleobjs = \
	test_parms_user.o \
	ghost_twopatch.o

# Solver routines (several solvers can be added here)
include Makefile.waveprop.inc

solverobjs =
solverincludes = $(waveprop_includes)

# Mapping routines
geoobjs = $(FORESTCLAW)/src/mappings/cart/set_maptype_cart.o

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
