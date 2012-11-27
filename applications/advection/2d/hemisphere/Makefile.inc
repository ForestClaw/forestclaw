# -*- Makefile -*-


waveprop := \
	$(FORESTCLAW)/src/amr_waveprop.o \
	$(FORESTCLAW)/AMRClaw/src/limiter.o \
	$(FORESTCLAW)/AMRClaw/src2d/src2.o \
	$(FORESTCLAW)/AMRClaw/src2d/fdisc.o \
	$(FORESTCLAW)/AMRClaw/src2d/flux2.o \
	$(FORESTCLAW)/AMRClaw/src2d/qad2.o \
	$(FORESTCLAW)/AMRClaw/src2d/step2.o \
	$(FORESTCLAW)/AMRClaw/src2d/ClawPatch2.o



solver_routines := \
	b4step2_mapped.o \
	setaux_mapped.o \
	setprob_advection.o \
	rpn2noncons.o \
	rpt2noncons.o \
	$(waveprop)
