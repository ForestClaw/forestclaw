# Makefile in src

include $(P4EST_DIR)/etc/Makefile.p4est.mk

.PHONY: clean

libforestclaw.a: amr_forestclaw.o fclaw2d_convenience.o forestclaw2d.o
	ar cru $@ $^

.c.o:
	$(P4EST_CC) $(P4EST_CFLAGS) $(P4EST_CPPFLAGS) -c $<

.cpp.o:
	$(P4EST_CXX) $(P4EST_CFLAGS) $(P4EST_CPPFLAGS) -I$(FORESTCLAW)/AMRClaw/src -I. -c $<

clean:
	rm -f *.o *.a
