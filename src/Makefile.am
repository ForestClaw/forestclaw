# This file is part of ForestClaw.
# Makefile.am in src
# included non-recursively from toplevel directory

libforestclaw_generated_headers = src/fclaw_config.h
libforestclaw_installed_headers = \
	src/fclaw_base.h \
	src/fclaw2d_global.h \
	src/fclaw2d_base.h \
	src/fclaw2d_domain.h \
	src/fclaw2d_defs.h \
	src/fclaw2d_block.h \
	src/fclaw2d_vtable.h \
	src/fclaw2d_output.h \
	src/ClawPatch.hpp \
	src/fclaw2d_output_ascii.h \
	src/fclaw2d_output_ascii_fort.h \
	src/fclaw2d_regrid_default.h \
	src/fclaw2d_regrid_default_fort.h \
	src/fclaw2d_update_single_step.h \
	src/fclaw2d_advance.hpp \
	src/subcycle_manager.hpp \
	src/fclaw2d_physical_bc.h \
	src/fclaw2d_ghost_fill.h \
	src/fclaw2d_partition.h \
	src/fclaw2d_exchange.h \
	src/fclaw2d_map.h \
	src/forestclaw2d.h \
	src/fclaw2d_convenience.h \
	src/fclaw2d_clawpatch.h \
	src/fclaw2d_clawpatch.hpp \
	src/fclaw2d_patch.hpp \
	src/fclaw2d_patch.h \
	src/fclaw2d_neighbors_fort.h \
	src/fclaw2d_transform.h \
	src/fclaw_package.h \
	src/fclaw_register.h \
	src/fortran_source2d/metric_terms.i \
	src/fclaw2d_map_query.h \
	src/fclaw2d_map_query_defs.h \
	src/fclaw2d_diagnostics.h \
	src/fp_exception_glibc_extension.h

libforestclaw_internal_headers =

# -------------------------------------------------------------------
# Main source files for ForestClaw are here - these need to be compiled
# -------------------------------------------------------------------

libforestclaw_compiled_sources = \
	src/fclaw_base.c \
	src/fclaw_options.c \
	src/fclaw_package.c \
	src/fclaw_register.c \
	src/fclaw_math.c \
	src/fclaw_timer.c \
	src/fclaw_mpi.c \
	src/fclaw2d_global.c \
	src/fclaw2d_base.c \
	src/fclaw2d_domain.cpp \
	src/fclaw2d_block.c \
	src/fclaw2d_patch.cpp \
	src/fclaw2d_convenience.c \
	src/forestclaw2d.c \
	src/fclaw2d_vtable.c \
	src/fclaw2d_map.c \
	src/fclaw2d_vtk.c \
	src/fclaw2d_transform.c \
	src/fclaw2d_map_query.c \
	src/fclaw2d_map_brick.c \
	src/fp_exception_glibc_extension.c \
	src/fclaw2d_output.cpp \
	src/fclaw2d_output_ascii.cpp \
	src/fclaw2d_regrid_default.cpp \
	src/fclaw2d_regrid.cpp \
	src/fclaw2d_manifold_default.cpp \
	src/fclaw2d_partition.cpp \
	src/fclaw2d_exchange.cpp \
	src/fclaw2d_diagnostics.cpp \
	src/fclaw2d_corner_neighbors.cpp \
	src/fclaw2d_face_neighbors.cpp \
	src/fclaw2d_ghost_fill.cpp \
	src/fclaw2d_physical_bc.cpp \
	src/fclaw2d_timeinterp.cpp \
	src/fclaw2d_clawpatch.cpp \
	src/fclaw2d_farraybox.cpp \
	src/mappings/context.f \
	src/mappings/scale_rot_map.f \
	src/mappings/map_query.f \
	src/mappings/cart/mapc2m_cart.f \
	src/mappings/annulus/mapc2m_annulus.f \
	src/mappings/brick/write_brick_data.f \
	src/mappings/identity/mapc2m_identity.f \
	src/mappings/pillowdisk/mapc2m_pillowdisk.f \
	src/mappings/squareddisk/mapc2m_squareddisk.f \
	src/mappings/pillowdisk5/mapc2m_pillowdisk5.f \
	src/mappings/cubedsphere/mapc2m_cubedsphere.f \
	src/mappings/pillowsphere/mapc2m_pillowsphere.f \
	src/mappings/fivepatch/mapc2m_fivepatch.f \
	src/mappings/torus/mapc2m_torus.f \
	src/mappings/latlong/mapc2m_latlong.f \
	src/mappings/pillowsphere/pillowsphere_corner_ghost.f \
	src/ClawPatch.cpp \
	src/fclaw2d_initialize.cpp \
	src/fclaw2d_run.cpp \
	src/fclaw2d_finalize.cpp \
	src/subcycle_manager.cpp \
	src/fclaw2d_advance.cpp \
	src/fclaw2d_update_single_step.cpp \
	src/fortran_source2d/fclaw2d_ghost_pack.f \
	src/fortran_source2d/fclaw2d_output_ascii_fort.f \
	src/fortran_source2d/fclaw2d_regrid_fort.f \
	src/fortran_source2d/philim.f \
	src/fortran_source2d/average_fv1.f \
	src/fortran_source2d/interpolate_fv1.f \
	src/fortran_source2d/exchange_fv1.f \
	src/fortran_source2d/diagnostics.f \
	src/fortran_source2d/phys_corner.f \
	src/fortran_source2d/setup_mesh.f \
	src/fortran_source2d/cellave2.f \
	src/fortran_source2d/timeinterp_interior.f \
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
