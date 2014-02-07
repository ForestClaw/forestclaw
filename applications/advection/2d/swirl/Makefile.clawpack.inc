# -*- Makefile -*-

# Use
#     gmake "version=4.6"
# to specify the version
# Default : version 5.0

version := 4.6

FCLAW2D_CLAWPACK = $(FORESTCLAW)/src/solvers/fclaw2d_clawpack${version}
FCLAW2D_CLAWPACK_FORTRAN = $(FCLAW2D_CLAWPACK)/fortran_source

clawpack_includes :=  -I$(FCLAW2D_CLAWPACK_FORTRAN) -I$(FCLAW2D_CLAWPACK)


ifeq ($(version),4.6)
clawpack_mods =
clawpack_objs = \
	$(FCLAW2D_CLAWPACK)/fclaw2d_clawpack.o \
	$(FCLAW2D_CLAWPACK_FORTRAN)/claw2d_utils.o \
	$(FCLAW2D_CLAWPACK_FORTRAN)/inlinelimiter.o \
	$(FCLAW2D_CLAWPACK_FORTRAN)/src2.o \
	$(FCLAW2D_CLAWPACK_FORTRAN)/bc2.o \
	$(FCLAW2D_CLAWPACK_FORTRAN)/flux2.o \
	$(FCLAW2D_CLAWPACK_FORTRAN)/step2.o \
	$(FCLAW2D_CLAWPACK_FORTRAN)/ClawPatch2.o
else
clawpack_mods = \
	$(FCLAW2D_CLAWPACK_FORTRAN)/amr_module.mod
clawpack_objs = \
	$(FCLAW2D_CLAWPACK)/fclaw2d_clawpack.o \
	$(FCLAW2D_CLAWPACK_FORTRAN)/claw2d_utils.o \
	$(FCLAW2D_CLAWPACK_FORTRAN)/set_amr_module.o \
	$(FCLAW2D_CLAWPACK_FORTRAN)/inlinelimiter.o \
	$(FCLAW2D_CLAWPACK_FORTRAN)/src2.o \
	$(FCLAW2D_CLAWPACK_FORTRAN)/bc2.o \
	$(FCLAW2D_CLAWPACK_FORTRAN)/flux2.o \
	$(FCLAW2D_CLAWPACK_FORTRAN)/step2.o \
	$(FCLAW2D_CLAWPACK_FORTRAN)/ClawPatch2.o
endif
