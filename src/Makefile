# Makefile in src

include $(P4EST_DIR)/etc/Makefile.p4est.mk

amrclawobjs = \
	$(FORESTCLAW)/AMRClaw/src/ClawPatch.o

libforestclaw.a: fclaw2d_convenience.o forestclaw2d.o \
        amr_forestclaw.o amr_utils.o amr_options.o $(amrclawobjs)
	ar cru $@ $^

.c.o:
	$(P4EST_CC) $(P4EST_CFLAGS) $(P4EST_CPPFLAGS) -c $< -o $@

.cpp.o:
	$(P4EST_CXX) $(P4EST_CFLAGS) $(P4EST_CPPFLAGS) \
	-I$(FORESTCLAW)/AMRClaw/src -I. -c $< -o $@

clean:
	rm -f *.o *.a $(amrclawobjs)
