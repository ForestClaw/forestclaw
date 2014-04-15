# Makefile for swirl example

TARGET = swirl
GEO = cart

include $(FORESTCLAW)/Make.defs.local

# -----------------------------------------------------------
FFLAGS = $(LOCAL_FFLAGS)

# -----------------------------------------------------------
# Define local objects here

exampleobjs = \
	swirl_user.o \
	psi.o \
	user_${version}/ghost_twopatch.o \
	user_${version}/tag4refinement.o \
	user_${version}/rpn2ad1.o \
	user_${version}/rpt2ad1.o \
	user_${version}/qinit.o \
	user_${version}/b4step2.o \
	user_${version}/setaux.o \
	user_${version}/setprob.o \
	user_${version}/swirl_output.o

solverobjs = $(clawpack_mods) $(clawpack_objs)
solverincludes = $(clawpack_includes)


OBJECTS = $(TARGET).o $(exampleobjs) $(solverobjs) $(localobjs)

# -----------------------------------------------------------
# The rest of this is more or less default.
# -----------------------------------------------------------

include $(P4EST_DIR)/etc/Makefile.p4est.mk

.SUFFIXES : .o .c .cpp .f .f90 .mod

# make sure that the first make target is the executable
default: $(TARGET)

.cpp.o:
	$(P4EST_CXX)  $(CPPFLAGS) $(P4EST_CFLAGS) $(P4EST_CPPFLAGS) \
	-I$(FORESTCLAW)/src -I. ${solverincludes} -c $< -o $@

.c.o:
	$(CC) $(P4EST_CFLAGS) $(P4EST_CPPFLAGS) \
	-I$(FORESTCLAW)/src -I. ${solverincludes} -c $< -o $@
.f90.mod:
	$(FC) -Wall $(FFLAGS) ${solverincludes} -c $<  -o $@

.f90.o:
	$(FC) -Wall $(FFLAGS) ${solverincludes} -c $<  -o $@

.f.o:
	$(FC) $(FFLAGS) -c $< ${solverincludes} -o $@


.PHONY: clean realclean default $(TARGET)

$(TARGET): $(OBJECTS) $(FORESTCLAW)/src/libforestclaw.a
	$(P4EST_CXX) $(P4EST_CFLAGS) $^ $(P4EST_LDFLAGS) $(FORT_LIB) \
	$(P4EST_LIBS) -o $@

clean:
	rm -f *.o $(TARGET) $(exampleobjs) $(solverobjs)

realclean:
	rm -f *.o $(TARGET) $(OBJECTS) fort.*
