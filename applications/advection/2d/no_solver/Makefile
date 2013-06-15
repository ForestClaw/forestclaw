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

# make sure that the first make target is the executable
default: $(TARGET)

.cpp.o:
	$(P4EST_CXX) $(P4EST_CFLAGS) $(P4EST_CPPFLAGS) \
	-I$(FORESTCLAW)/AMRClaw/src -I. -I$(FORESTCLAW)/src $(solver_objs) -c $< -o $@

.c.o:
	$(P4EST_CXX) $(P4EST_CFLAGS) $(P4EST_CPPFLAGS) \
	-I$(FORESTCLAW)/AMRClaw/src -I. -I$(FORESTCLAW)/src $(solver_objs) -c $< -o $@

.f.o:
	$(FC) -Wall $(FFLAGS) -c $< $(solver_objs) -o $@


.PHONY: clean realclean default $(TARGET)

exampleobjs = \
	no_solver_user.o \
	initialize.o \
	ghost_twopatch.o \
	tag4refinement.o


# Lots of these routines are also only used for Clawpack, but I haven't
# and should all really be listed under 'solver_routines'.  But for now
# I've left them here.
amrclawobjs = \
	$(FORTRAN_SRC)/philim.o \
	$(FORTRAN_SRC)/claw_out2.o \
	$(FORTRAN_SRC)/mapc2p.o \
	$(FORTRAN_SRC)/mapc2m.o \
	$(FORTRAN_SRC)/setup_mesh.o \
	$(FORTRAN_SRC)/claw2d_utils_mapped.o \
	$(FORTRAN_SRC)/claw2d_utils.o \
	$(FORTRAN_SRC)/debug.o \
	$(FORTRAN_SRC)/cellave2.o \
	$(FORESTCLAW)/AMRClaw/src2d/set_maptype.o \
	$(FORESTCLAW)/AMRClaw/src2d/scale_rot_map.o \
	$(FORESTCLAW)/AMRClaw/src2d/mappings/block.o \
	$(FORESTCLAW)/AMRClaw/src2d/mappings/mapc2m_cart.o \
	$(FORESTCLAW)/AMRClaw/src2d/mappings/mapc2m_disk.o \
	$(FORESTCLAW)/AMRClaw/src2d/mappings/mapc2m_mb_sphere.o

GeoFiles := \
      $(FORESTCLAW)/AMRClaw/src2d/mappings/set_maptype_cart.o

OBJECTS = $(TARGET).o $(exampleobjs) $(amrclawobjs) $(GeoFiles)

$(TARGET): $(OBJECTS) $(FORESTCLAW)/src/libforestclaw.a
	$(P4EST_CXX) $(P4EST_CFLAGS) $^ $(P4EST_LDFLAGS) $(FORT_LIB) \
	$(P4EST_LIBS) -o $@

clean:
	rm -f *.o $(TARGET) $(OBJECTS)

realclean: clean
	rm -f $(amrclawobjs) fort.*
