# Makefile in src

TARGET = hemisphere

GEO = sphere

FC = gfortran
FFLAGS = -O0 -g -Wno-unused-dummy-argument -Wno-unused-variable
FORT_LIB = -L/usr/local/gfortran/lib -lgfortran

include $(P4EST_DIR)/etc/Makefile.p4est.mk

solver_objs :=
solver_includes :=

# ---------------------------------
include Makefile.waveprop.inc
solver_objs += $(waveprop_objs)
solver_includes += $(waveprop_includes)

# ---------------------------------

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
	velocity_wind.o \
	tag4refinement.o \
	rpn2noncons.o \
	rpt2noncons.o \
	mapc2m.o \
	$(solver_objs)

amrclawobjs = \
	$(FORESTCLAW)/AMRClaw/src/philim.o \
	$(FORESTCLAW)/AMRClaw/src2d/mapc2p.o \
	$(FORESTCLAW)/AMRClaw/src2d/setup_mesh.o \
	$(FORESTCLAW)/AMRClaw/src2d/claw2d_utils.o \
	$(FORESTCLAW)/AMRClaw/src2d/claw_out2.o \
	$(FORESTCLAW)/AMRClaw/src2d/set_maptype.o \
	$(FORESTCLAW)/AMRClaw/src2d/scale_rot_map.o \
	$(FORESTCLAW)/AMRClaw/src2d/claw2d_utils_mapped.o \
	$(FORESTCLAW)/AMRClaw/src2d/mappings/block.o \
	$(FORESTCLAW)/AMRClaw/src2d/mappings/mapc2m_cart.o \
	$(FORESTCLAW)/AMRClaw/src2d/mappings/mapc2m_disk.o \
	$(FORESTCLAW)/AMRClaw/src2d/mappings/mapc2m_mb_sphere.o

ifeq ($(GEO),sphere)
   GeoFiles := \
      $(FORESTCLAW)/AMRClaw/src2d/mappings/set_maptype_sphere.o
endif

ifeq ($(GEO),disk)
   GeoFiles := \
      $(FORESTCLAW)/AMRClaw/src2d/mappings/set_maptype_disk.o
endif

ifeq ($(GEO),cart)
   GeoFiles := \
      $(FORESTCLAW)/AMRClaw/src2d/mappings/set_maptype_cart.o
endif

OBJECTS = $(TARGET).o $(exampleobjs) $(amrclawobjs) $(GeoFiles)

$(TARGET): $(OBJECTS) $(FORESTCLAW)/src/libforestclaw.a
	$(P4EST_CXX) $(P4EST_CFLAGS) $^ $(P4EST_LDFLAGS) $(FORT_LIB) \
	$(P4EST_LIBS) -o $@

clean:
	rm -f *.o $(TARGET) $(OBJECTS)

realclean:
	rm -f *.o $(OBJECTS) fort.* $(TARGET)
