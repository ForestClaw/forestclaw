
# This file is part of ForestClaw.
# Makefile.am in src
# included non-recursively from toplevel directory

libforestclaw_generated_headers = src/fclaw_config.h
libforestclaw_installed_headers = \
	src/fclaw_base.h \
	src/amr_options.h \
	src/fclaw2d_base.h \
	src/forestclaw2d.h \
	src/fclaw2d_convenience.h \
	src/fclaw2d_map.h \
	src/clawpack_fort.H

libforestclaw_internal_headers =

# -------------------------------------------------------------------
# Main source files for ForestClaw are here - these need to be compiled
# -------------------------------------------------------------------

libforestclaw_compiled_sources = \
	src/fclaw_base.c \
	src/amr_options.c \
	src/fclaw2d_base.c \
	src/forestclaw2d.c \
	src/fclaw2d_convenience.c \
	src/fclaw2d_map.c \
	src/clawpack_fort.cpp \
	src/mappings/block.f \
	src/mappings/set_maptype.f \
	src/mappings/scale_rot_map.f \
	src/mappings/cart/mapc2m_cart.f \
	src/mappings/disk/mapc2m_disk.f \
	src/mappings/sphere/mapc2m_sphere.f \
	src/ClawPatch.cpp \
	src/amrinit.cpp \
	src/amrrun.cpp \
	src/amrreset.cpp \
	src/subcycle_manager.cpp \
	src/amr_advance.cpp \
	src/amr_regrid.cpp \
	src/amr_output.cpp \
	src/amr_diagnostics.cpp \
	src/amr_repartition_domain.cpp \
	src/amr_exchangecoarse.cpp \
	src/amr_levelexchange.cpp \
	src/amr_mol.cpp \
	src/amr_single_step.cpp \
	src/amr_neighbors.cpp \
	src/amr_timeinterp.cpp \
	src/amr_physbc.cpp \
	src/amr_utils.cpp \
	src/regrid_functions.cpp \
	src/output_functions.cpp \
	src/fclaw2d_solvers.cpp \
	src/fortran_source2d/philim.f \
	src/fortran_source2d/claw_out2.f \
	src/fortran_source2d/mapc2p.f \
	src/fortran_source2d/mapc2m.f \
	src/fortran_source2d/setup_mesh.f \
	src/fortran_source2d/claw2d_utils_mapped.f \
	src/fortran_source2d/claw2d_utils.f \
	src/fortran_source2d/debug.f \
	src/fortran_source2d/cellave2.f


# this variable is used for headers that are not publicly installed
FCLAW_CPPFLAGS =

lib_LTLIBRARIES += src/libforestclaw.la
src_libforestclaw_la_SOURCES = \
 $(libforestclaw_internal_headers) \
 $(libforestclaw_compiled_sources)
src_libforestclaw_la_CPPFLAGS = $(AM_CPPFLAGS) $(FCLAW_CPPFLAGS)
LDADD += @top_builddir@/src/libforestclaw.la @FCLAW_P4EST_LDADD@ @FCLAW_SC_LDADD@
EXTRA_src_libforestclaw_la_DEPENDENCIES = @FCLAW_P4EST_LDADD@ @FCLAW_SC_LDADD@

nodist_include_HEADERS += $(libforestclaw_generated_headers)
include_HEADERS += $(libforestclaw_installed_headers)

AM_CPPFLAGS += -I@top_srcdir@/src @FCLAW_SC_CPPFLAGS@ @FCLAW_P4EST_CPPFLAGS@
