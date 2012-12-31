# -*- Makefile -*-

AMR_WAVEPROP = $(FORESTCLAW)/src/solvers/amr_waveprop
AMR_WAVEPROP_FORTRAN = $(AMR_WAVEPROP)/fortran_source

solver_includes = -I$(AMR_WAVEPROP) -I$(AMR_WAVEPROP_FORTRAN)

waveprop_objs := \
	$(AMR_WAVEPROP)/amr_waveprop.o \
	$(AMR_WAVEPROP_FORTRAN)/limiter.o \
	$(AMR_WAVEPROP_FORTRAN)/src2.o \
	$(AMR_WAVEPROP_FORTRAN)/src2_manifold.o \
	$(AMR_WAVEPROP_FORTRAN)/b4step2.o \
	$(AMR_WAVEPROP_FORTRAN)/setaux.o \
	$(AMR_WAVEPROP_FORTRAN)/qinit.o \
	$(AMR_WAVEPROP_FORTRAN)/bc2.o \
	$(AMR_WAVEPROP_FORTRAN)/bc2_manifold.o \
	$(AMR_WAVEPROP_FORTRAN)/fdisc.o \
	$(AMR_WAVEPROP_FORTRAN)/flux2.o \
	$(AMR_WAVEPROP_FORTRAN)/qad2.o \
	$(AMR_WAVEPROP_FORTRAN)/step2.o \
	$(AMR_WAVEPROP_FORTRAN)/ClawPatch2.o




solver_objs := \
	qinit_manifold.o \
	b4step2_manifold.o \
	setaux_manifold.o \
	setprob_advection.o \
	rpn2noncons.o \
	rpt2noncons.o \
	$(waveprop_objs)
