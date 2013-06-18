# Makefile in src

TARGET = sphere

GEO = sphere

FC = gfortran
FFLAGS = -O0 -g -Wno-unused-dummy-argument -Wno-unused-variable
FORT_LIB = -L/usr/local/gfortran/lib -lgfortran

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
	-I$(FORESTCLAW)/src -I$(FORESTCLAW)/AMRClaw/src -I. ${solver_includes} -c $< -o $@

.c.o:
	$(P4EST_CXX) $(P4EST_CFLAGS) $(P4EST_CPPFLAGS) \
	-I$(FORESTCLAW)/src -I$(FORESTCLAW)/AMRClaw/src -I. ${solver_includes} -c $< -o $@

.f.o:
	$(FC) -Wall $(FFLAGS) -c $< ${solver_includes} -o $@

.PHONY: clean realclean default $(TARGET)

exampleobjs = \
	sphere_user.o \
	ghost_twopatch.o \
	velocity_wind.o \
	tag4refinement.o \
	mapc2m.o \
	rpn2noncons.o \
	rpt2noncons.o \
	$(solver_objs)

GeoFiles := $(FORESTCLAW)/AMRClaw/src2d/mappings/set_maptype_sphere.o

OBJECTS = $(TARGET).o $(exampleobjs) $(GeoFiles) $(debug_files)

$(TARGET): $(OBJECTS) $(FORESTCLAW)/src/libforestclaw.a
	$(P4EST_CXX) $(P4EST_CFLAGS) $^ $(P4EST_LDFLAGS) $(FORT_LIB) \
	$(P4EST_LIBS) -o $@

clean:
	rm -f *.o $(TARGET) $(OBJECTS)

realclean:
	rm -f *.o $(exampleobjs) fort.* $(TARGET)
