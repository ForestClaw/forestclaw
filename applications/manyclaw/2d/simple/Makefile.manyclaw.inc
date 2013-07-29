# -*- Makefile -*-

MANYCLAW = $(MANYCLAW_DIR)

AMR_MANYCLAW = $(FORESTCLAW)/src/solvers/amr_manyclaw
AMR_MANYCLAW_FORTRAN = $(AMR_MANYCLAW)/fortran_source

manyclaw_linkers := -lmanyclaw -L$(MANYCLAW)/lib -L$(TBBROOT)/lib/intel64/gcc4.4 -ltbb
manyclaw_includes := -I$(AMR_MANYCLAW) -I$(AMR_MANYCLAW_FORTRAN) -I$(MANYCLAW) -I$(TBBROOT)/include

manyclaw_objs = \
	$(AMR_MANYCLAW)/amr_manyclaw.o \
	$(AMR_MANYCLAW)/fortran_source/manyclaw_reorder.o \
	$(AMR_MANYCLAW)/fortran_source/ClawPatch2.o \
	$(AMR_MANYCLAW)/fortran_source/src2.o \
	$(AMR_MANYCLAW)/fortran_source/step2.o \
	$(AMR_MANYCLAW)/fortran_source/bc2.o \
	$(AMR_MANYCLAW)/fortran_source/flux2.o \
	$(AMR_MANYCLAW)/fortran_source/limiter.o
