# -*- Makefile -*-
libclawpack25d_46_generated_headers =

libclawpack25d_46_installed_headers = \
	src/solvers/fc3d_clawpack46/fc3d_clawpack46.h \
	src/solvers/fc3d_clawpack46/fc3d_clawpack46_options.h \
	src/solvers/fc3d_clawpack46/fc3d_clawpack46_user_fort.h \
	src/solvers/fc3d_clawpack46/fc3d_clawpack46_fort.h 

libclawpack25d_46_compiled_sources = \
	src/solvers/fc3d_clawpack46/fc3d_clawpack46.cpp \
	src/solvers/fc3d_clawpack46/fc3d_clawpack46_options.c \
    src/solvers/fc3d_clawpack46/fortran_source/fc3d_clawpack46_time_sync.f \
	src/solvers/fc3d_clawpack46/fortran_source/fc3d_clawpack46_inlinelimiter.f \
	src/solvers/fc3d_clawpack46/fortran_source/fc3d_clawpack46_set_capacity.f \
	src/solvers/fc3d_clawpack46/fortran_source/fc3d_clawpack46_block.f \
	src/solvers/fc3d_clawpack46/fortran_source/fc3d_clawpack46_bc2_default.f \
	src/solvers/fc3d_clawpack46/fortran_source/fc3d_clawpack46_flux2.f \
	src/solvers/fc3d_clawpack46/fortran_source/fc3d_clawpack46_flux2fw.f \
	src/solvers/fc3d_clawpack46/fortran_source/fc3d_clawpack46_step2.f \
	src/solvers/fc3d_clawpack46/fortran_source/fc3d_clawpack46_step2_wrap.f


## Name of library to build
lib_LTLIBRARIES += src/solvers/fc3d_clawpack46/libclawpack25d_46.la

## Named variables that can be referenced from other libraries/apps
FCLAW_CLAWPACK25d_46_CPPFLAGS = -I@top_srcdir@/src/solvers/fc3d_clawpack46 
FCLAW_CLAWPACK25d_46_LDADD    = @top_builddir@/src/solvers/fc3d_clawpack46/libclawpack25d_46.la

## Sources needed to build this library
src_solvers_fc3d_clawpack46_libclawpack25d_46_la_SOURCES = \
          $(libclawpack25d_46_internal_headers) \
          $(libclawpack25d_46_compiled_sources)

## Headers needed to build this library.  Only AM_CPPFLAGS required; others are 
## included in AM_CPPFLAGS, but are included here for clarity (although makes 
## for lengthy compile lines)
src_solvers_fc3d_clawpack46_libclawpack25d_46_la_CPPFLAGS = \
          $(AM_CPPFLAGS) \
          $(FCLAW25d_46_CPPFLAGS) \
          $(FCLAW25d_46_CLAWPACK_CPPFLAGS) \
          $(FCLAW_CLAWPATCH_CPPFLAGS)

## Add headers from this library to general header pool, to be used by other apps
AM_CPPFLAGS  += $(FCLAW_CLAWPACK25d_46_CPPFLAGS)

LDADD        += $(FCLAW_CLAWPACK25d_46_LDADD)

include_HEADERS += $(libclawpack25d_46_installed_headers)
