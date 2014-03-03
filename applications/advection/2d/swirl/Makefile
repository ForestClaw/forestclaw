# Makefile for sphere example

TARGET = swirl
GEO = cart

include $(FORESTCLAW)/Make.defs.local

# -----------------------------------------------------------
FFLAGS = $(LOCAL_FFLAGS)

# -----------------------------------------------------------
# Define local objects here
exampleobjs = \
	swirl_user.o \
	ghost_twopatch.o \
	tag4refinement.o \
	psi.o \
	rpn2ad1.o \
	rpt2ad1.o \
	qinit.o \
	b4step2.o \
	setaux.o \
	setprob.o \
	swirl_output.o

# Solver routines (several solvers can be added here)
include Makefile.waveprop.inc
solverobjs = $(waveprop_objs)
solverincludes = $(waveprop_includes)


OBJECTS = $(TARGET).o $(exampleobjs) $(solverobjs) $(localobjs)

# -----------------------------------------------------------
# The rest of this is more or less default.
# -----------------------------------------------------------

include $(P4EST_DIR)/etc/Makefile.p4est.mk

# make sure that the first make target is the executable
default: $(TARGET)

.cpp.o:
	$(P4EST_CXX) $(P4EST_CFLAGS) $(P4EST_CPPFLAGS) $(CPPFLAGS) \
	-I$(FORESTCLAW)/src -I. ${solverincludes} -c $< -o $@

.c.o:
	$(CC) $(P4EST_CFLAGS) $(P4EST_CPPFLAGS) \
	-I$(FORESTCLAW)/src -I. ${solverincludes} -c $< -o $@
.f.o:
	$(FC) $(FFLAGS) -c $< ${solverincludes} -o $@

.PHONY: clean realclean default $(TARGET)

$(TARGET): $(OBJECTS) $(FORESTCLAW)/src/libforestclaw.a
	$(P4EST_CXX) $(P4EST_CFLAGS) $^ $(P4EST_LDFLAGS) $(FORT_LIB) \
	$(P4EST_LIBS) -o $@

clean:
	rm -f *.o $(TARGET) $(exampleobjs)

realclean:
	rm -f *.o $(TARGET) $(OBJECTS) fort.*
