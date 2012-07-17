# Makefile in src

include $(P4EST_DIR)/etc/Makefile.p4est.mk

.PHONY: clean

libforestclaw.a: forestclaw.o
	ar cru $@ $<

.c.o:
	$(P4EST_CC) $(P4EST_CFLAGS) $(P4EST_CPPFLAGS) -c $^

clean:
	rm -f *.o
