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
	src/forestclaw2d.h \
	src/fp_exception_glibc_extension.h \
	src/fclaw2d_include_all.h \
	src/fclaw2d_defs.h \
	src/fclaw2d_options.h \
	src/fclaw2d_advance.h \
	src/fclaw2d_elliptic_solver.h \
	src/fclaw2d_farraybox.hpp \
	src/fclaw2d_global.h \
	src/fclaw2d_forestclaw.h \
	src/fclaw2d_domain.h \
	src/fclaw2d_block.h \
	src/fclaw2d_patch.h \
	src/fclaw2d_vtable.h \
	src/fclaw2d_output.h \
	src/fclaw2d_time_sync.h \
	src/fclaw2d_update_single_step.h \
	src/fclaw2d_physical_bc.h \
	src/fclaw2d_rays.h \
	src/fclaw2d_ghost_fill.h \
	src/fclaw2d_corner_neighbors.h \
	src/fclaw2d_face_neighbors.h \
	src/fclaw2d_map.h \
	src/fclaw2d_timeinterp.h \
	src/fclaw2d_regrid.h \
	src/fclaw2d_map_brick.h \
	src/fclaw2d_partition.h \
	src/fclaw2d_exchange.h \
	src/fclaw2d_convenience.h \
	src/fclaw2d_map_query.h \
	src/fclaw2d_map_query_defs.h \
	src/fclaw2d_diagnostics.h \
	src/fclaw2d_file.h \
        src/fclaw2d_to_3d.h \
	src/forestclaw3d.h \
	src/fclaw3d_defs.h \
	src/fclaw3d_convenience.h \
	src/fclaw3d_file.h

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
	src/fclaw2d_block.c \
	src/fclaw2d_options.c \
	src/fclaw2d_global.c \
	src/fclaw2d_forestclaw.c \
	src/fclaw2d_vtable.c \
	src/fclaw2d_map.c \
	src/fclaw2d_map_query.c \
	src/fclaw2d_map_brick.c \
	src/fclaw2d_convenience.c \
	src/fclaw2d_exchange.c \
	src/fclaw2d_advance.c \
	src/fclaw2d_elliptic_solver.c \
	src/fclaw2d_physical_bc.c \
	src/fclaw2d_rays.c \
	src/fclaw2d_ghost_fill.c \
	src/fclaw2d_output.c \
	src/fclaw2d_run.c \
	src/fclaw2d_diagnostics.c \
	src/fclaw2d_update_single_step.c \
	src/fclaw2d_domain.c \
	src/fclaw2d_regrid.c \
	src/fclaw2d_patch.c \
	src/fclaw2d_partition.c \
	src/fclaw2d_timeinterp.c \
	src/fclaw2d_initialize.c \
	src/fclaw2d_finalize.c \
	src/fclaw2d_time_sync.c \
	src/fclaw2d_corner_neighbors.c \
	src/fclaw2d_face_neighbors.c \
	src/fclaw2d_farraybox.cpp \
	src/fclaw2d_output_tikz.c \
	src/fclaw2d_file.c \
	src/forestclaw2d.c \
	src/forestclaw3d.c \
	src/fclaw3d_convenience.c \
	src/fclaw3d_file.c \
	src/fp_exception_glibc_extension.c \
	src/mappings/context.f \
	src/mappings/scale_rot_map.f \
	src/mappings/map_query.f \
	src/mappings/cart/mapc2m_cart.f \
	src/mappings/bilinear/mapc2m_bilinear.f e
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

## Extra flags that are not provided by general AM_CPPFLAGS
FCLAW_CPPFLAGS =

## In case we like to build only the library and no programs
justlibs-local: src/libforestclaw.la

## Using libtool to build the ForestClaw library
lib_LTLIBRARIES += src/libforestclaw.la
src_libforestclaw_la_SOURCES = \
	$(libforestclaw_internal_headers) \
 	$(libforestclaw_compiled_sources)
src_libforestclaw_la_CPPFLAGS = $(AM_CPPFLAGS) $(FCLAW_CPPFLAGS)

## This is the official API versioning scheme of libtool.  Please see:
## Read https://www.gnu.org/software/libtool/manual/libtool.html#Versioning
src_libforestclaw_la_LDFLAGS = -version-info 2:0:0
src_libforestclaw_la_LIBADD = @FCLAW_P4EST_LIBADD@ @FCLAW_SC_LIBADD@

## Explicitly depend on other libraries
EXTRA_src_libforestclaw_la_DEPENDENCIES = @FCLAW_P4EST_EDEPS@ @FCLAW_SC_EDEPS@

## Code p4est and libsc installation path into ForestClaw executables
AM_LDFLAGS += @FCLAW_P4EST_RPATH@ @FCLAW_SC_RPATH@

## Flags to support linking against ForestClaw and support libraries
LDADD += src/libforestclaw.la \
	@FCLAW_P4EST_LDADD@ \
	@FCLAW_SC_LDADD@
FCLAW_LDADD = src/libforestclaw.la \
	@FCLAW_P4EST_LDADD@ \
	@FCLAW_SC_LDADD@

nodist_include_HEADERS += $(libforestclaw_generated_headers)
include_HEADERS += $(libforestclaw_installed_headers)

## These flags are passed to populate the include path
AM_CPPFLAGS += \
	-I@top_srcdir@/src \
	-I@top_srcdir@/src/fortran_source2d \
	@FCLAW_SC_CPPFLAGS@ @FCLAW_P4EST_CPPFLAGS@

## UNIT TESTS
check_PROGRAMS += src/forestclaw.TEST
TESTS += src/forestclaw.TEST

src_forestclaw_TEST_SOURCES = \
    src/fclaw_gauges.h.TEST.cpp \
    src/fclaw_pointer_map.h.TEST.cpp \
	src/fclaw2d_elliptic_solver.h.TEST.cpp \
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
    $(LDADD)
