# -*- Makefile -*-

version := 4.6

FCLAW2D_CLAWPACK = $(FORESTCLAW)/src/solvers/fclaw2d_clawpack${version}
FCLAW2D_CLAWPACK_FORTRAN = $(FCLAW2D_CLAWPACK)/fortran_source

clawpack_includes := -I$(FCLAW2D_CLAWPACK) -I$(FCLAW2D_CLAWPACK_FORTRAN)

clawpack_objs = \
	$(FCLAW2D_CLAWPACK)/fclaw2d_clawpack.o \
	$(FCLAW2D_CLAWPACK_FORTRAN)/limiter.o \
	$(FCLAW2D_CLAWPACK_FORTRAN)/src2.o \
	$(FCLAW2D_CLAWPACK_FORTRAN)/b4step2.o \
	$(FCLAW2D_CLAWPACK_FORTRAN)/setaux.o \
	$(FCLAW2D_CLAWPACK_FORTRAN)/setprob.o \
	$(FCLAW2D_CLAWPACK_FORTRAN)/qinit.o \
	$(FCLAW2D_CLAWPACK_FORTRAN)/bc2.o \
	$(FCLAW2D_CLAWPACK_FORTRAN)/flux2.o \
	$(FCLAW2D_CLAWPACK_FORTRAN)/step2.o \
	$(FCLAW2D_CLAWPACK_FORTRAN)/ClawPatch2.o
