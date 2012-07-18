# Makefile in src

include $(P4EST_DIR)/etc/Makefile.p4est.mk

.PHONY: clean examples

examples: example1

libforestclaw.a: amr_forestclaw.o fclaw2d_convenience.o forestclaw2d.o
	ar cru $@ $^

.c.o:
	$(P4EST_CC) $(P4EST_CFLAGS) $(P4EST_CPPFLAGS) -c $<

.cpp.o:
	$(P4EST_CXX) $(P4EST_CFLAGS) $(P4EST_CPPFLAGS) -I$(FORESTCLAW)/AMRClaw/src -I. -c $<

example1: example1.o libforestclaw.a
	$(P4EST_CXX) $(P4EST_CFLAGS) $^ $(P4EST_LDFLAGS) $(P4EST_LIBS) -o $@

clean:
	rm -f *.o *.a
