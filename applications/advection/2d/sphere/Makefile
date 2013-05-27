# Makefile in src

TARGET = sphere

GEO = sphere

FC = gfortran
FFLAGS = -O0 -g -Wno-unused-dummy-argument -Wno-unused-variable
FORT_LIB = -L/usr/local/gfortran/lib -lgfortran

FORTRAN_SRC = $(FORESTCLAW)/src/fortran_source2d
FORTRAN_INC = $(FORTRAN_SRC)

include $(P4EST_DIR)/etc/Makefile.p4est.mk

solver_objs :=
solver_includes :=
debug_files :=

include Makefile.waveprop.inc

solver_objs += $(waveprop_objs)
solver_includes += $(waveprop_includes)
debug_files += $(FORTRAN_SRC)/debug.o

# make sure that the first make target is the executable
default: $(TARGET)

.cpp.o:
	$(P4EST_CXX) $(P4EST_CFLAGS) $(P4EST_CPPFLAGS) \
	-I$(FORESTCLAW)/src -I$(FORESTCLAW)/AMRClaw/src -I. ${solver_includes} -c $< -o $@

.c.o:
	$(P4EST_CXX) $(P4EST_CFLAGS) $(P4EST_CPPFLAGS) \
	-I$(FORESTCLAW)/src -I$(FORESTCLAW)/AMRClaw/src -I. ${solver_includes} -c $< -o $@

.f.o:
	$(FC) -Wall $(FFLAGS) -I$(FORTRAN_INC) -c $< ${solver_includes} -o $@

.PHONY: clean realclean default $(TARGET)

exampleobjs = \
	ghost_twopatch.o \
	velocity_wind.o \
	tag4refinement.o \
	mapc2m.o \
	rpn2noncons.o \
	rpt2noncons.o \
	$(solver_objs)

amrclawobjs = \
	$(FORTRAN_SRC)/philim.o \
	$(FORTRAN_SRC)/setup_mesh.o \
	$(FORTRAN_SRC)/claw2d_utils.o \
	$(FORTRAN_SRC)/claw_out2.o \
	$(FORTRAN_SRC)/claw2d_utils_mapped.o \
	$(FORESTCLAW)/AMRClaw/src2d/set_maptype.o \
	$(FORESTCLAW)/AMRClaw/src2d/scale_rot_map.o \
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

OBJECTS = $(TARGET).o $(exampleobjs) $(amrclawobjs) $(GeoFiles) $(debug_files)

$(TARGET): $(OBJECTS) $(FORESTCLAW)/src/libforestclaw.a
	$(P4EST_CXX) $(P4EST_CFLAGS) $^ $(P4EST_LDFLAGS) $(FORT_LIB) \
	$(P4EST_LIBS) -o $@

clean:
	rm -f *.o $(TARGET) $(OBJECTS)

realclean:
	rm -f *.o $(amrclawobjs) fort.* $(TARGET)
