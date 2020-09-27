# -*- Makefile -*-
libactiveflux_generated_headers =

libactiveflux_installed_headers = \
	src/solvers/fc2d_activeflux/fc2d_activeflux.h \
	src/solvers/fc2d_activeflux/fc2d_activeflux_options.h \
	src/solvers/fc2d_activeflux/activeflux_user_fort.h \
	src/solvers/fc2d_activeflux/fc2d_activeflux_fort.h 

libactiveflux_compiled_sources = \
	src/solvers/fc2d_activeflux/fc2d_activeflux.cpp \
	src/solvers/fc2d_activeflux/fc2d_activeflux_options.c \
    src/solvers/fc2d_activeflux/fortran_source/activeflux_time_sync.f \
	src/solvers/fc2d_activeflux/fortran_source/activeflux_set_capacity.f \
	src/solvers/fc2d_activeflux/fortran_source/activeflux_block.f \
	src/solvers/fc2d_activeflux/fortran_source/activeflux_bc2_default.f \
	src/solvers/fc2d_activeflux/fortran_source/activeflux_flux2.f \
	src/solvers/fc2d_activeflux/fortran_source/activeflux_flux2fw.f \
	src/solvers/fc2d_activeflux/fortran_source/activeflux_inlinelimiter.f \
	src/solvers/fc2d_activeflux/fortran_source/activeflux_step2.f \
	src/solvers/fc2d_activeflux/fortran_source/activeflux_step2_wrap.f


## Name of library to build
lib_LTLIBRARIES += src/solvers/fc2d_activeflux/libactiveflux.la

## Named variables that can be referenced from other libraries/apps
FCLAW_activeflux_CPPFLAGS = -I@top_srcdir@/src/solvers/fc2d_activeflux 
FCLAW_activeflux_LDADD    = @top_builddir@/src/solvers/fc2d_activeflux/libactiveflux.la

## Sources needed to build this library
src_solvers_fc2d_activeflux_libactiveflux_la_SOURCES = \
          $(libactiveflux_internal_headers) \
          $(libactiveflux_compiled_sources)

## Headers needed to build this library.  Only AM_CPPFLAGS required; others are 
## included in AM_CPPFLAGS, but are included here for clarity (although makes 
## for lengthy compile lines)
src_solvers_fc2d_activeflux_libactiveflux_la_CPPFLAGS = \
          $(AM_CPPFLAGS) \
          $(FCLAW_CPPFLAGS) \
          $(FCLAW_activeflux_CPPFLAGS) \
          $(FCLAW_CLAWPATCH_CPPFLAGS)

## Add headers from this library to general header pool, to be used by other apps
AM_CPPFLAGS  += $(FCLAW_activeflux_CPPFLAGS)

LDADD        += $(FCLAW_activeflux_LDADD)

include_HEADERS += $(libactiveflux_installed_headers)
