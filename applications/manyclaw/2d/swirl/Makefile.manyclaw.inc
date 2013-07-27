# -*- Makefile -*-

MANYCLAW = $(FORESTCLAW)/../ManyClaw

AMR_MANYCLAW = $(FORESTCLAW)/src/solvers/amr_manyclaw
AMR_MANYCLAW_FORTRAN = $(AMR_MANYCLAW)/fortran_source

manyclaw_linkers := -lmanyclaw -L$(MANYCLAW)/lib -L/opt/local/lib -ltbb
manyclaw_includes := -I$(AMR_MANYCLAW) -I$(AMR_MANYCLAW_FORTRAN) -I$(MANYCLAW) -I/opt/local/include

manyclaw_objs = \
	$(AMR_MANYCLAW)/amr_manyclaw.o \
	$(AMR_MANYCLAW)/fortran_source/manyclaw_reorder.o \
	$(AMR_MANYCLAW_FORTRAN)/limiter.o \
	$(AMR_MANYCLAW_FORTRAN)/src2.o \
	$(AMR_MANYCLAW_FORTRAN)/bc2.o \
	$(AMR_MANYCLAW_FORTRAN)/flux2.o \
	$(AMR_MANYCLAW_FORTRAN)/qad2.o \
	$(AMR_MANYCLAW_FORTRAN)/step2.o \
	$(AMR_MANYCLAW_FORTRAN)/ClawPatch2.o
