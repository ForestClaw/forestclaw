
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
	src/clawpack_fort.H \
	src/fclaw2d_transform.h \
	src/fortran_source2d/metric_terms.i \
	src/fclaw2d_map_query.h \
	src/fclaw2d_map_query_defs.h

libforestclaw_internal_headers = \
	src/fortran_source2d/debug.i

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
	src/fclaw2d_vtk.c \
	src/fclaw2d_transform.c \
	src/fclaw2d_map_query.c \
	src/mappings/block.f \
	src/mappings/context.f \
	src/mappings/scale_rot_map.f \
	src/mappings/cart/mapc2m_cart.f \
	src/mappings/pillowdisk/mapc2m_pillowdisk.f \
	src/mappings/squareddisk/mapc2m_squareddisk.f \
	src/mappings/map_query.f \
	src/mappings/cubedsphere/mapc2m_cubedsphere.f \
	src/mappings/pillowsphere/mapc2m_pillowsphere.f \
	src/mappings/pillowsphere/pillowsphere_corner_ghost.f \
	src/clawpack_fort.cpp \
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
	src/fclaw2d_corner_neighbors.cpp \
	src/fclaw2d_face_neighbors.cpp \
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
	src/fortran_source2d/average_fv1.f \
	src/fortran_source2d/interpolate_fv1.f \
	src/fortran_source2d/exchange_fv1.f \
	src/fortran_source2d/diagnostics.f \
	src/fortran_source2d/phys_corner.f \
	src/fortran_source2d/setup_mesh.f \
	src/fortran_source2d/claw_out2.f \
	src/fortran_source2d/debug.f \
	src/fortran_source2d/cellave2.f \
	src/fortran_source2d/fdisc.f

# this variable is used for headers that are not publicly installed
FCLAW_CPPFLAGS = -I@top_srcdir@/src -I@top_srcdir@/src/fortran_source2d
FCLAW_FFLAGS = $(FCLAW_CPPFLAGS)

lib_LTLIBRARIES += src/libforestclaw.la
src_libforestclaw_la_SOURCES = \
 $(libforestclaw_internal_headers) \
 $(libforestclaw_compiled_sources)
# src_libforestclaw_la_CPPFLAGS = $(AM_CPPFLAGS)
LDADD += @top_builddir@/src/libforestclaw.la @FCLAW_P4EST_LDADD@ @FCLAW_SC_LDADD@
EXTRA_src_libforestclaw_la_DEPENDENCIES = @FCLAW_P4EST_LDADD@ @FCLAW_SC_LDADD@

nodist_include_HEADERS += $(libforestclaw_generated_headers)
include_HEADERS += $(libforestclaw_installed_headers)

AM_CPPFLAGS += $(FCLAW_CPPFLAGS) @FCLAW_SC_CPPFLAGS@ @FCLAW_P4EST_CPPFLAGS@
AM_FFLAGS += $(FCLAW_FFLAGS)
