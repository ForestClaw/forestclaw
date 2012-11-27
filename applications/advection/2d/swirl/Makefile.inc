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
	b4step2.o \
	setaux.o \
	setprob.o \
	psi.o \
	rpn2ad1.o \
	rpt2ad1.o \
	$(waveprop)
