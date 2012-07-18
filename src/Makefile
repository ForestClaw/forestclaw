# Makefile in src

include $(P4EST_DIR)/etc/Makefile.p4est.mk

.PHONY: clean examples

examples: example1

amrclawobjs = \
	$(FORESTCLAW)/forestclaw/AMRClaw/src/philim.o \
	$(FORESTCLAW)/forestclaw/AMRClaw/src/limiter.o \
	$(FORESTCLAW)/forestclaw/AMRClaw/src2d/claw2d_utils.o \
	$(FORESTCLAW)/forestclaw/AMRClaw/src/ClawPatch.o


libforestclaw.a: fclaw2d_convenience.o forestclaw2d.o amr_forestclaw.o $(amrclawobjs)
	ar cru $@ $^

.c.o:
	$(P4EST_CC) $(P4EST_CFLAGS) $(P4EST_CPPFLAGS) -c $< -o $@

.cpp.o:
	$(P4EST_CXX) $(P4EST_CFLAGS) $(P4EST_CPPFLAGS) -I$(FORESTCLAW)/forestclaw/AMRClaw/src -I. -c $< -o $@

.f.o:
	$(FC) -c $< -o $@

example1: example1.o libforestclaw.a
	$(P4EST_CXX) $(P4EST_CFLAGS) $^ $(P4EST_LDFLAGS) $(P4EST_LIBS) -o $@

clean:
	rm -f *.o *.a
