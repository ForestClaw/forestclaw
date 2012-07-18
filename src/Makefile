# Makefile in src

include $(P4EST_DIR)/etc/Makefile.p4est.mk

.PHONY: clean examples

examples: example1

amrclawobjs = \
	qinit.o \
	rpn2.o \
	rpt2.o \
	$(FORESTCLAW)/forestclaw/AMRClaw/src/philim.o \
	$(FORESTCLAW)/forestclaw/AMRClaw/src/limiter.o \
	$(FORESTCLAW)/forestclaw/AMRClaw/src/ClawPatch.o \
	$(FORESTCLAW)/forestclaw/AMRClaw/src2d/claw2d_utils.o \
	$(FORESTCLAW)/forestclaw/AMRClaw/src2d/bc2.o \
	$(FORESTCLAW)/forestclaw/AMRClaw/src2d/ClawPatch2.o \
	$(FORESTCLAW)/forestclaw/AMRClaw/src2d/setaux.o \
	$(FORESTCLAW)/forestclaw/AMRClaw/src2d/average2.o \
	$(FORESTCLAW)/forestclaw/AMRClaw/src2d/b4step2.o \
	$(FORESTCLAW)/forestclaw/AMRClaw/src2d/cellave2.o \
	$(FORESTCLAW)/forestclaw/AMRClaw/src2d/estimate_error2.o \
	$(FORESTCLAW)/forestclaw/AMRClaw/src2d/fdisc.o \
	$(FORESTCLAW)/forestclaw/AMRClaw/src2d/fixcapaq2.o \
	$(FORESTCLAW)/forestclaw/AMRClaw/src2d/flux2.o \
	$(FORESTCLAW)/forestclaw/AMRClaw/src2d/flux2fw.o \
	$(FORESTCLAW)/forestclaw/AMRClaw/src2d/mapc2p.o \
	$(FORESTCLAW)/forestclaw/AMRClaw/src2d/qad2.o \
	$(FORESTCLAW)/forestclaw/AMRClaw/src2d/setprob.o \
	$(FORESTCLAW)/forestclaw/AMRClaw/src2d/step2.o \
	$(FORESTCLAW)/forestclaw/AMRClaw/src2d/src2.o


FORT_LIB = -L/usr/local/gfortran/lib  -lgfortran

libforestclaw.a: fclaw2d_convenience.o forestclaw2d.o amr_forestclaw.o $(amrclawobjs)
	ar cru $@ $^

.c.o:
	$(P4EST_CC) $(P4EST_CFLAGS) $(P4EST_CPPFLAGS) -c $< -o $@

.cpp.o:
	$(P4EST_CXX) $(P4EST_CFLAGS) $(P4EST_CPPFLAGS) -I$(FORESTCLAW)/forestclaw/AMRClaw/src -I. -c $< -o $@

.f.o:
	$(FC) -c $< -o $@

example1: example1.o libforestclaw.a
	$(P4EST_CXX) $(P4EST_CFLAGS) $^ $(P4EST_LDFLAGS) $(FORT_LIB) $(P4EST_LIBS) -o $@

clean:
	rm -f *.o *.a
