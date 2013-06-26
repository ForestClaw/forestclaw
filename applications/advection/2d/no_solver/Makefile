# Makefile in src

TARGET = no_solver
GEO = cart



FC = gfortran
FFLAGS = -O0 -g -Wno-unused-dummy-argument -Wno-unused-variable
FORT_LIB = -L/usr/local/gfortran/lib -lgfortran

FORTRAN_SRC = $(FORESTCLAW)/src/fortran_source2d
FORTRAN_INC = $(FORTRAN_SRC)

include $(P4EST_DIR)/etc/Makefile.p4est.mk

solver_objs :=
solver_includes :=

# Include solver makefiles here....
#....

# Add addtional objects here
solver_objs +=
solver_includes +=

exampleobjs = \
	no_solver_user.o \
	initialize_ac.o \
	initialize.o \
	ghost_twopatch.o \
	tag4refinement.o


# make sure that the first make target is the executable
default: $(TARGET)

.cpp.o:
	$(P4EST_CXX) $(P4EST_CFLAGS) $(P4EST_CPPFLAGS) \
	-I. -I$(FORESTCLAW)/src $(solver_objs) -c $< -o $@

.c.o:
	$(P4EST_CXX) $(P4EST_CFLAGS) $(P4EST_CPPFLAGS) \
	-I$(FORESTCLAW)/AMRClaw/src -I. -I$(FORESTCLAW)/src $(solver_objs) -c $< -o $@

.f.o:
	$(FC) -Wall $(FFLAGS) -I$(FORTRAN_INC) -c $< $(solver_objs) -o $@


.PHONY: clean realclean default $(TARGET)

GeoFiles := \
      $(FORESTCLAW)/src/mappings/cart/set_maptype_cart.o

OBJECTS = $(TARGET).o $(exampleobjs) $(amrclawobjs) $(GeoFiles)

$(TARGET): $(OBJECTS) $(FORESTCLAW)/src/libforestclaw.a
	$(P4EST_CXX) $(P4EST_CFLAGS) $^ $(P4EST_LDFLAGS) $(FORT_LIB) \
	$(P4EST_LIBS) -o $@

clean:
	rm -f *.o $(TARGET) $(OBJECTS)

realclean: clean
	rm -f $(amrclawobjs) fort.*
