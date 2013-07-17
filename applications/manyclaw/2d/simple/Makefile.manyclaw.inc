# -*- Makefile -*-

AMR_MANYCLAW = $(FORESTCLAW)/src/solvers/amr_manyclaw
AMR_MANYCLAW_FORTRAN = $(AMR_MANYCLAW)/fortran_source

manyclaw_linkers := -lmanyclaw -L$(MANYCLAW)/lib
manyclaw_includes := -I$(AMR_MANYCLAW) -I$(AMR_MANYCLAW_FORTRAN)

manyclaw_objs = \
	$(AMR_MANYCLAW)/amr_manyclaw.o \
