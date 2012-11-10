# Makefile in src

include $(P4EST_DIR)/etc/Makefile.p4est.mk

amrclawobjs = \
	$(FORESTCLAW)/AMRClaw/src/ClawPatch.o \
	$(FORESTCLAW)/AMRClaw/src/clawpack_fort.o \
	$(FORESTCLAW)/AMRClaw/src/parser.o \
	$(FORESTCLAW)/AMRClaw/src/parse_inputs.o

fclawobjs = \
	amrinit.o \
	amrrun.o \
	amrreset.o \
	subcycle_manager.o \
	amr_advance.o \
	amr_regrid.o \
	amr_output.o \
	amr_diagnostics.o \
	amr_exchangecoarse.o \
	amr_levelexchange.o \
	amr_neighbors.o \
	amr_physbc.o \
	amr_utils.o \
	amr_options.o \
	fclaw2d_convenience.o \
	forestclaw2d.o

libforestclaw.a: $(amrclawobjs) $(fclawobjs)
	ar cru $@ $^

.c.o:
	$(P4EST_CC) $(P4EST_CFLAGS) $(P4EST_CPPFLAGS) -c $< -o $@

.cpp.o:
	$(P4EST_CXX) $(P4EST_CFLAGS) $(P4EST_CPPFLAGS) \
	-I$(FORESTCLAW)/AMRClaw/src -I. -c $< -o $@

clean:
	rm -f *.o *.a $(amrclawobjs) $(fclawobjs)
