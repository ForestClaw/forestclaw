# -*- Makefile -*-
libgeoclaw_generated_headers =
libgeoclaw_installed_headers = \
	src/solvers/fc2d_geoclaw/fc2d_geoclaw.h \
	src/solvers/fc2d_geoclaw/fc2d_geoclaw_options.h

libgeoclaw_compiled_sources = \
	src/solvers/fc2d_geoclaw/fc2d_geoclaw.cpp \
	src/solvers/fc2d_geoclaw/fc2d_geoclaw_options.c \
	src/solvers/fc2d_geoclaw/amrlib_source/amr_module.f90 \
	src/solvers/fc2d_geoclaw/geolib_source/geoclaw_module.f90 \
	src/solvers/fc2d_geoclaw/geolib_source/utility_module.f90 \
	src/solvers/fc2d_geoclaw/geolib_source/topo_module.f90 \
	src/solvers/fc2d_geoclaw/amrlib_source/regions_module.f90 \
	src/solvers/fc2d_geoclaw/geolib_source/gauges_module.f90 \
	src/solvers/fc2d_geoclaw/geolib_source/qinit_module.f90 \
	src/solvers/fc2d_geoclaw/geolib_source/surge/holland_storm_module.f90 \
	src/solvers/fc2d_geoclaw/geolib_source/surge/constant_storm_module.f90 \
	src/solvers/fc2d_geoclaw/geolib_source/surge/stommel_storm_module.f90 \
	src/solvers/fc2d_geoclaw/geolib_source/surge/storm_module.f90 \
	src/solvers/fc2d_geoclaw/geolib_source/friction_module.f90 \
	src/solvers/fc2d_geoclaw/fortran_source/geoclaw_set_modules.f90 \
	src/solvers/fc2d_geoclaw/fortran_source/geoclaw_bc2.f \
	src/solvers/fc2d_geoclaw/amrlib_source/opendatafile.f \
	src/solvers/fc2d_geoclaw/geolib_source/topointegral.f \
	src/solvers/fc2d_geoclaw/geolib_source/bilinearintegral.f \
	src/solvers/fc2d_geoclaw/geolib_source/setaux.f90 \
	src/solvers/fc2d_geoclaw/geolib_source/src2.f90 \
	src/solvers/fc2d_geoclaw/geolib_source/b4step2.f90 \
	src/solvers/fc2d_geoclaw/geolib_source/qinit.f90 \
	src/solvers/fc2d_geoclaw/geolib_source/cellgridintegrate2.f \
	src/solvers/fc2d_geoclaw/fortran_source/geoclaw_block.f

lib_LTLIBRARIES += src/solvers/fc2d_geoclaw/libgeoclaw.la

src_solvers_fc2d_geoclaw_libgeoclaw_la_SOURCES = \
  $(libgeoclaw_internal_headers) \
  $(libgeoclaw_compiled_sources)

src_solvers_fc2d_geoclaw_libgeoclaw_la_CPPFLAGS = \
  $(AM_CPPFLAGS) $(FCLAW_CPPFLAGS)

FCLAW_GEOCLAW_CPPFLAGS = -I@top_srcdir@/src/solvers/fc2d_geoclaw
FCLAW_GEOCLAW_LDADD = @top_builddir@/src/solvers/fc2d_geoclaw/libgeoclaw.la

include_HEADERS += $(libgeoclaw_installed_headers)
