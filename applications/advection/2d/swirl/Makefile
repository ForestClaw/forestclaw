# Makefile in src

TARGET = swirl

include $(P4EST_DIR)/etc/Makefile.p4est.mk

.PHONY: clean $(TARGET)

exampleobjs = \
	qinit.o \
	rpn2.o \
	rpt2.o \
	setaux.o \
	setprob.o \
	b4step2.o \
	psi.o

amrclawobjs = \
	$(FORESTCLAW)/AMRClaw/src/philim.o \
	$(FORESTCLAW)/AMRClaw/src/limiter.o \
	$(FORESTCLAW)/AMRClaw/src2d/mapc2p.o \
	$(FORESTCLAW)/AMRClaw/src2d/src2.o \
	$(FORESTCLAW)/AMRClaw/src2d/bc2.o \
	$(FORESTCLAW)/AMRClaw/src2d/fdisc.o \
	$(FORESTCLAW)/AMRClaw/src2d/estimate_error2.o \
	$(FORESTCLAW)/AMRClaw/src2d/claw2d_utils.o \
	$(FORESTCLAW)/AMRClaw/src2d/ClawPatch2.o \
	$(FORESTCLAW)/AMRClaw/src2d/average2.o \
	$(FORESTCLAW)/AMRClaw/src2d/cellave2.o \
	$(FORESTCLAW)/AMRClaw/src2d/fixcapaq2.o \
	$(FORESTCLAW)/AMRClaw/src2d/flux2.o \
	$(FORESTCLAW)/AMRClaw/src2d/qad2.o \
	$(FORESTCLAW)/AMRClaw/src2d/step2.o \
	$(FORESTCLAW)/AMRClaw/src2d/claw_out2.o

FORT_LIB = -L/usr/local/gfortran/lib -lgfortran

# .c.o:
# 	$(P4EST_CC) $(P4EST_CFLAGS) $(P4EST_CPPFLAGS) -c $< -o $@

.cpp.o:
	$(P4EST_CXX) $(P4EST_CFLAGS) $(P4EST_CPPFLAGS) \
	-I$(FORESTCLAW)/AMRClaw/src -I. -I$(FORESTCLAW)/src/ -c $< -o $@

.f.o:
	$(FC) -Wall -Wno-unused-dummy-argument -c $< -I$(FORESTCLAW)/AMRClaw/src2d -o $@

$(TARGET): $(TARGET).o $(FORESTCLAW)/src/libforestclaw.a $(exampleobjs) $(amrclawobjs)
	$(P4EST_CXX) $(P4EST_CFLAGS) $^ $(P4EST_LDFLAGS) $(FORT_LIB) \
	$(P4EST_LIBS) -o $@

clean:
	rm -f *.o

realclean:
	rm -f *.o $(amrclawobjs) fort.*
