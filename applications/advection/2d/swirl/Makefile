# Makefile in src

TARGET = swirl

GEO = cart

FC = gfortran
FFLAGS = -O0 -g -Wno-unused-dummy-argument -Wno-unused-variable
FORT_LIB = -L/usr/local/gfortran/lib -lgfortran

FORTRAN_SRC = $(FORESTCLAW)/src/fortran_source2d

include $(P4EST_DIR)/etc/Makefile.p4est.mk

solver_objs :=
solver_includes :=

include Makefile.waveprop.inc

solver_objs += $(waveprop_objs)
solver_includes += $(waveprop_includes)

# make sure that the first make target is the executable
default: $(TARGET)

.cpp.o:
	$(P4EST_CXX) $(P4EST_CFLAGS) $(P4EST_CPPFLAGS) \
	-I$(FORESTCLAW)/AMRClaw/src -I$(FORESTCLAW)/src/ $(solver_includes) -c $< -o $@

.c.o:
	$(P4EST_CXX) $(P4EST_CFLAGS) $(P4EST_CPPFLAGS) \
	-I$(FORESTCLAW)/AMRClaw/src -I. -I$(FORESTCLAW)/src/ -c $< -o $@

.f.o:
	$(FC) -Wall $(FFLAGS) -c $< $(solver_includes) -o $@


.PHONY: clean realclean default $(TARGET)

exampleobjs = \
	ghost_twopatch.o \
	tag4refinement.o \
	psi.o \
	rpn2ad1.o \
	rpt2ad1.o \
	swirl_user.o \
	$(solver_objs)


amrclawobjs = \
	$(FORTRAN_SRC)/philim.o \
	$(FORTRAN_SRC)/claw_out2.o \
	$(FORTRAN_SRC)/mapc2p.o \
	$(FORTRAN_SRC)/mapc2m.o \
	$(FORTRAN_SRC)/setup_mesh.o \
	$(FORTRAN_SRC)/claw2d_utils_mapped.o \
	$(FORTRAN_SRC)/claw2d_utils.o \
	$(FORTRAN_SRC)/debug.o \
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
	rm -f $(TARGET) $(example_objs) $(GeoFiles) fort.*
