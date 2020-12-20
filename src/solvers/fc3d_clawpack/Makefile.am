# -*- Makefile -*-
libclawpack_generated_headers =

libclawpack_installed_headers = \
	src/solvers/fc3d_clawpack/fc3d_clawpack.h \
	src/solvers/fc3d_clawpack/fc3d_clawpack_options.h \
	src/solvers/fc3d_clawpack/clawpack_user_fort.h \
	src/solvers/fc3d_clawpack/fc3d_clawpack_fort.h 

libclawpack_compiled_sources = \
	src/solvers/fc3d_clawpack/fc3d_clawpack.cpp \
	src/solvers/fc3d_clawpack/fc3d_clawpack_options.c \
    src/solvers/fc3d_clawpack/fortran_source/clawpack_time_sync.f \
	src/solvers/fc3d_clawpack/fortran_source/clawpack_inlinelimiter.f \
	src/solvers/fc3d_clawpack/fortran_source/clawpack_set_capacity.f \
	src/solvers/fc3d_clawpack/fortran_source/clawpack_block.f \
	src/solvers/fc3d_clawpack/fortran_source/clawpack_bc3_default.f \
	src/solvers/fc3d_clawpack/fortran_source/clawpack_flux3.f \
	src/solvers/fc3d_clawpack/fortran_source/clawpack_flux3fw.f \
	src/solvers/fc3d_clawpack/fortran_source/clawpack_step3.f \
	src/solvers/fc3d_clawpack/fortran_source/clawpack_step3_wrap.f


## Name of library to build
lib_LTLIBRARIES += src/solvers/fc3d_clawpack/libclawpack.la

## Named variables that can be referenced from other libraries/apps
FCLAW_CLAWPACK_CPPFLAGS = -I@top_srcdir@/src/solvers/fc3d_clawpack 
FCLAW_CLAWPACK_LDADD    = @top_builddir@/src/solvers/fc3d_clawpack/libclawpack.la

## Sources needed to build this library
src_solvers_fc3d_clawpack_libclawpack_la_SOURCES = \
          $(libclawpack_internal_headers) \
          $(libclawpack_compiled_sources)

## Headers needed to build this library.  Only AM_CPPFLAGS required; others are 
## included in AM_CPPFLAGS, but are included here for clarity (although makes 
## for lengthy compile lines)
src_solvers_fc3d_clawpack_libclawpack_la_CPPFLAGS = \
          $(AM_CPPFLAGS) \
          $(FCLAW_CPPFLAGS) \
          $(FCLAW_CLAWPACK_CPPFLAGS) \
          $(FCLAW_CLAWPATCH_CPPFLAGS)

## Add headers from this library to general header pool, to be used by other apps
AM_CPPFLAGS  += $(FCLAW_CLAWPACK_CPPFLAGS)

LDADD        += $(FCLAW_CLAWPACK_LDADD)

include_HEADERS += $(libclawpack_installed_headers)
