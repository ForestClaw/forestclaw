# This file is part of ForestClaw.
# Makefile.am in src
# included non-recursively from toplevel directory

libforestclaw_generated_headers = src/fclaw_config.h

libforestclaw_installed_headers = \
	src/fclaw_base.h \
	src/fclaw_timer.h \
	src/fclaw_package.h \
	src/fclaw_packing.h \
	src/fclaw_filesystem.h \
	src/fclaw_pointer_map.h \
	src/fclaw_options.h \
	src/fclaw_gauges.h \
	src/fclaw_mpi.h \
	src/fclaw_math.h \
	src/fclaw_wrap.h \
	src/fclaw_convenience.h \
	src/forestclaw.h \
	src/forestclaw2d.h \
  	src/fclaw_diagnostics.h \
  	src/fclaw_domain.h \
  	src/fclaw_patch.h \
  	src/fclaw_exchange.h \
  	src/fclaw_global.h \
	src/fclaw_ghost_fill.h \
    src/fclaw_corner_neighbors.h \
    src/fclaw_edge_neighbors.h \
    src/fclaw_face_neighbors.h \
    src/fclaw_block.h \
    src/fclaw_physical_bc.h \
    src/fclaw_partition.h \
    src/fclaw_regrid.h \
  	src/fclaw_vtable.h \
  	src/fclaw_output.h \
	src/fclaw_advance.h \
  	src/fclaw_timeinterp.h \
  	src/fclaw_update_single_step.h \
  	src/fclaw_time_sync.h \
	src/fclaw_forestclaw.h \
	src/fclaw_include_all.h \
  	src/fclaw_rays.h \
  	src/fclaw_elliptic_solver.h \
	src/fclaw_farraybox.hpp \
	src/fp_exception_glibc_extension.h \
	src/fclaw2d_include_all.h \
	src/fclaw2d_defs.h \
	src/fclaw2d_options.h \
	src/fclaw2d_global.h \
	src/fclaw2d_domain.h \
	src/fclaw2d_patch.h \
	src/fclaw2d_vtable.h \
	src/fclaw2d_output.h \
	src/fclaw2d_physical_bc.h \
	src/fclaw2d_ghost_fill.h \
	src/fclaw2d_corner_neighbors.h \
	src/fclaw2d_face_neighbors.h \
	src/fclaw2d_map.h \
	src/fclaw2d_regrid.h \
	src/fclaw2d_map_brick.h \
	src/fclaw2d_partition.h \
	src/fclaw2d_exchange.h \
	src/fclaw2d_convenience.h \
	src/fclaw2d_map_query.h \
	src/fclaw2d_map_query_defs.h \
	src/fclaw2d_diagnostics.h \
	src/fclaw2d_file.h \
	src/fclaw2d_wrap.h \
    src/fclaw2d_to_3d.h \
	src/forestclaw3d.h \
	src/fclaw3d_include_all.h \
	src/fclaw3d_defs.h \
	src/fclaw3d_exchange.h \
	src/fclaw3d_options.h \
	src/fclaw3d_patch.h \
	src/fclaw3d_domain.h \
	src/fclaw3d_global.h \
	src/fclaw3d_convenience.h \
	src/fclaw3d_map.h \
	src/fclaw3d_file.h \
	src/fclaw3d_wrap.h 

libforestclaw_internal_headers =

# -------------------------------------------------------------------
# Main source files for ForestClaw are here - these need to be compiled
# -------------------------------------------------------------------

libforestclaw_compiled_sources = \
	src/fclaw_base.c \
	src/fclaw_options.c \
	src/fclaw_gauges.c \
	src/fclaw_package.c \
	src/fclaw_packing.c \
	src/fclaw_filesystem.cpp \
	src/fclaw_pointer_map.c \
	src/fclaw_math.c \
	src/fclaw_timer.c \
	src/fclaw_mpi.c \
    src/fclaw_domain.c \
    src/fclaw_patch.c \
    src/fclaw_exchange.c \
    src/fclaw_diagnostics.c \
	src/fclaw_convenience.c \
	src/fclaw_ghost_fill.c \
	src/fclaw_global.c \
    src/fclaw_corner_neighbors.c \
    src/fclaw_edge_neighbors.c \
    src/fclaw_face_neighbors.c \
    src/fclaw_block.c \
    src/fclaw_physical_bc.c \
    src/fclaw_partition.c \
    src/fclaw_regrid.c \
  	src/fclaw_vtable.c \
  	src/fclaw_output.c \
	src/fclaw_advance.c \
  	src/fclaw_timeinterp.c \
  	src/fclaw_update_single_step.c \
  	src/fclaw_time_sync.c \
	src/fclaw_finalize.c \
	src/fclaw_initialize.c \
  	src/fclaw_rays.c \
    src/fclaw_forestclaw.c \
  	src/fclaw_elliptic_solver.c \
	src/fclaw_run.c \
	src/fclaw_farraybox.cpp \
	src/forestclaw.c \
	src/fclaw2d_options.c \
	src/fclaw2d_global.c \
	src/fclaw2d_vtable.c \
	src/fclaw2d_map.c \
	src/fclaw2d_map_query.c \
	src/fclaw2d_map_brick.c \
	src/fclaw2d_convenience.c \
	src/fclaw2d_exchange.c \
	src/fclaw2d_physical_bc.c \
	src/fclaw2d_ghost_fill.c \
	src/fclaw2d_output.c \
	src/fclaw2d_diagnostics.c \
	src/fclaw2d_domain.c \
	src/fclaw2d_regrid.c \
	src/fclaw2d_patch.c \
	src/fclaw2d_partition.c \
	src/fclaw2d_corner_neighbors.c \
	src/fclaw2d_face_neighbors.c \
	src/fclaw2d_output_tikz.c \
	src/fclaw2d_file.c \
	src/fclaw2d_wrap.c \
	src/forestclaw2d.c \
	src/forestclaw3d.c \
	src/fclaw3d_convenience.c \
	src/fclaw3d_domain.c \
	src/fclaw3d_exchange.c \
	src/fclaw3d_options.c \
	src/fclaw3d_patch.c \
	src/fclaw3d_global.c \
	src/fclaw3d_map.c \
	src/fclaw3d_file.c \
	src/fclaw3d_wrap.c \
	src/fp_exception_glibc_extension.c \
	src/mappings/context.f \
	src/mappings/scale_rot_map.f \
	src/mappings/map_query.f \
	src/mappings/cart/mapc2m_cart.f \
	src/mappings/bilinear/mapc2m_bilinear.f \
	src/mappings/annulus/mapc2m_annulus.f \
	src/mappings/identity/mapc2m_identity.f \
	src/mappings/pillowdisk/mapc2m_pillowdisk.f \
	src/mappings/squareddisk/mapc2m_squareddisk.f \
	src/mappings/pillowdisk5/mapc2m_pillowdisk5.f \
	src/mappings/cubedsphere/mapc2m_cubedsphere.f \
	src/mappings/pillowsphere/mapc2m_pillowsphere.f \
	src/mappings/fivepatch/mapc2m_fivepatch.f \
	src/mappings/torus/mapc2m_torus.f \
	src/mappings/torus/mapc2m_twisted_torus.f \
	src/mappings/latlong/mapc2m_latlong.f \
	src/mappings/fclaw2d_map_nomap.c \
	src/mappings/fclaw2d_map_nomap_brick.c \
	src/fortran_source2d/fclaw_dopri5.f \
	src/fortran_source2d/cellave2.f \
	src/fortran_source2d/cellave3.f \
	src/fortran_source2d/fdisc.f \
	src/fortran_source2d/philim.f

# this variable is used for headers that are not publicly installed
FCLAW_CPPFLAGS = \
    -I@top_srcdir@/src \
	-I@top_srcdir@/src/fortran_source2d

FCLAW_FFLAGS = $(FCLAW_CPPFLAGS)   # (orignally used for metric.i;  not needed now)

lib_LTLIBRARIES += src/libforestclaw.la
src_libforestclaw_la_SOURCES = \
	$(libforestclaw_internal_headers) \
 	$(libforestclaw_compiled_sources)

src_libforestclaw_la_CPPFLAGS = $(AM_CPPFLAGS) $(FCLAW_CPPFLAGS)

FCLAW_LDADD = @top_builddir@/src/libforestclaw.la \
	          @FCLAW_P4EST_LDADD@ \
              @FCLAW_SC_LDADD@


LDADD += @top_builddir@/src/libforestclaw.la \
	     @FCLAW_P4EST_LDADD@ \
	     @FCLAW_SC_LDADD@

EXTRA_src_libforestclaw_la_DEPENDENCIES = @FCLAW_P4EST_LDADD@ @FCLAW_SC_LDADD@

nodist_include_HEADERS += $(libforestclaw_generated_headers)
include_HEADERS += $(libforestclaw_installed_headers)

AM_CPPFLAGS += $(FCLAW_CPPFLAGS) @FCLAW_SC_CPPFLAGS@ @FCLAW_P4EST_CPPFLAGS@
AM_FFLAGS += $(FCLAW_FFLAGS)

## UNIT TESTS
check_PROGRAMS += src/forestclaw.TEST
TESTS += src/forestclaw.TEST

src_forestclaw_TEST_SOURCES = \
    src/fclaw_gauges.h.TEST.cpp \
    src/fclaw_pointer_map.h.TEST.cpp \
    src/fclaw_global.h.TEST.cpp \
    src/fclaw_options.h.TEST.cpp \
    src/fclaw_patch.h.TEST.cpp \
	src/fclaw_elliptic_solver.h.TEST.cpp \
	src/fclaw2d_diagnostics.h.TEST.cpp \
	src/fclaw2d_global.h.TEST.cpp \
	src/fclaw2d_options.h.TEST.cpp \
	src/fclaw2d_patch.h.TEST.cpp \
	src/fclaw2d_vtable.h.TEST.cpp 

src_forestclaw_TEST_CPPFLAGS = \
	$(test_libtestutils_la_CPPFLAGS) \
    $(src_libforestclaw_la_CPPFLAGS)

src_forestclaw_TEST_LDADD = \
	test/libtestutils.la \
	$(test_libtestutils_la_LDADD) \
    $(LDADD) \
    $(FCLAW_LDADD)
