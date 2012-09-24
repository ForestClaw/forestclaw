# Makefile in src

TARGET = swirl

GEO = cart

FC = gfortran
FFLAGS = -O0 -g -Wno-unused-dummy-argument -Wno-unused-variable
FORT_LIB = -L/usr/local/gfortran/lib -lgfortran


.cpp.o:
	$(P4EST_CXX) $(P4EST_CFLAGS) $(P4EST_CPPFLAGS) \
	-I$(FORESTCLAW)/AMRClaw/src -I. -I$(FORESTCLAW)/src/ -c $< -o $@

.f.o:
	$(FC) -Wall $(FFLAGS) -c $< -I$(FORESTCLAW)/AMRClaw/src2d -o $@


include $(P4EST_DIR)/etc/Makefile.p4est.mk

.PHONY: clean $(TARGET)

exampleobjs = \
	qinit.o \
	b4step2.o \
	setaux.o \
	setprob.o \
        psi.o \
	ghost_twopatch.o \
	rpn2ad1.o \
	rpt2ad1.o \
	init_routines.o \
	tag4refinement.o

amrclawobjs = \
	$(FORESTCLAW)/AMRClaw/src/philim.o \
	$(FORESTCLAW)/AMRClaw/src/limiter.o \
	$(FORESTCLAW)/AMRClaw/src2d/mapc2p.o \
	$(FORESTCLAW)/AMRClaw/src2d/b4step2_mapped.o \
	$(FORESTCLAW)/AMRClaw/src2d/setaux_mapped.o \
	$(FORESTCLAW)/AMRClaw/src2d/qinit_mapped.o \
	$(FORESTCLAW)/AMRClaw/src2d/input_parms.o \
	$(FORESTCLAW)/AMRClaw/src2d/src2.o \
	$(FORESTCLAW)/AMRClaw/src2d/bc2.o \
	$(FORESTCLAW)/AMRClaw/src2d/fdisc.o \
	$(FORESTCLAW)/AMRClaw/src2d/setup_mesh.o \
	$(FORESTCLAW)/AMRClaw/src2d/claw2d_utils.o \
	$(FORESTCLAW)/AMRClaw/src2d/ClawPatch2.o \
	$(FORESTCLAW)/AMRClaw/src2d/cellave2.o \
	$(FORESTCLAW)/AMRClaw/src2d/flux2.o \
	$(FORESTCLAW)/AMRClaw/src2d/qad2.o \
	$(FORESTCLAW)/AMRClaw/src2d/step2.o \
	$(FORESTCLAW)/AMRClaw/src2d/claw_out2.o \
	$(FORESTCLAW)/AMRClaw/src2d/set_maptype.o \
	$(FORESTCLAW)/AMRClaw/src2d/scale_rot_map.o \
	$(FORESTCLAW)/AMRClaw/src2d/claw2d_utils_mapped.o \
	$(FORESTCLAW)/AMRClaw/src2d/mappings/block.o \
	$(FORESTCLAW)/AMRClaw/src2d/mapc2m.o \
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

realclean:
	rm -f *.o $(amrclawobjs) fort.* $(TARGET)
