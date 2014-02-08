# -*- Makefile -*-

AMR_WAVEPROP = $(FORESTCLAW)/src/solvers/amr_waveprop
AMR_WAVEPROP_FORTRAN = $(AMR_WAVEPROP)/fortran_source

waveprop_includes := -I$(AMR_WAVEPROP) -I$(AMR_WAVEPROP_FORTRAN)

waveprop_objs = \
	$(AMR_WAVEPROP)/amr_waveprop.o \
	$(AMR_WAVEPROP_FORTRAN)/limiter.o \
	$(AMR_WAVEPROP_FORTRAN)/src2.o \
	$(AMR_WAVEPROP_FORTRAN)/bc2.o \
	$(AMR_WAVEPROP_FORTRAN)/flux2.o \
	$(AMR_WAVEPROP_FORTRAN)/step2.o \
	$(AMR_WAVEPROP_FORTRAN)/ClawPatch2.o
