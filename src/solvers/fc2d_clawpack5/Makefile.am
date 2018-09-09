# -*- Makefile -*-
libclawpack5_generated_headers =
libclawpack5_installed_headers = \
	src/solvers/fc2d_clawpack5/fc2d_clawpack5.h \
	src/solvers/fc2d_clawpack5/fc2d_clawpack5_options.h \
	src/solvers/fc2d_clawpack5/clawpack5_user_fort.h \
	src/solvers/fc2d_clawpack5/fc2d_clawpack5_fort.h

libclawpack5_compiled_sources = \
	src/solvers/fc2d_clawpack5/fc2d_clawpack5.cpp \
	src/solvers/fc2d_clawpack5/fc2d_clawpack5_options.c \
	src/solvers/fc2d_clawpack5/fortran_source/clawpack5_amr_module.f90 \
	src/solvers/fc2d_clawpack5/fortran_source/clawpack5_set_amr_module.f90 \
	src/solvers/fc2d_clawpack5/fortran_source/clawpack5_inlinelimiter.f \
	src/solvers/fc2d_clawpack5/fortran_source/clawpack5_limiter.f \
	src/solvers/fc2d_clawpack5/fortran_source/clawpack5_set_capacity.f \
	src/solvers/fc2d_clawpack5/fortran_source/clawpack5_block.f \
	src/solvers/fc2d_clawpack5/fortran_source/clawpack5_bc2_default.f \
	src/solvers/fc2d_clawpack5/fortran_source/clawpack5_flux2.f \
	src/solvers/fc2d_clawpack5/fortran_source/clawpack5_step2.f90 \
	src/solvers/fc2d_clawpack5/fortran_source/clawpack5_step2_wrap.f

## Name of library to build
lib_LTLIBRARIES += src/solvers/fc2d_clawpack5/libclawpack5.la

## Named variables that can be reference from other libraries/apps
FCLAW_CLAWPACK5_CPPFLAGS = -I@top_srcdir@/src/solvers/fc2d_clawpack5
FCLAW_CLAWPACK5_LDADD    = @top_builddir@/src/solvers/fc2d_clawpack5/libclawpack5.la

## Sources needed to build this library
src_solvers_fc2d_clawpack5_libclawpack5_la_SOURCES = \
  $(libclawpack5_internal_headers) \
  $(libclawpack5_compiled_sources)

## Headers needed to build this library.  Only AM_CPPFLAGS required; others are 
## included in AM_CPPFLAGS, but are included here for clarity (although makes 
## for lengthy compile lines)
src_solvers_fc2d_clawpack5_libclawpack5_la_CPPFLAGS = \
              $(AM_CPPFLAGS) \
              $(FCLAW_CPPFLAGS) \
              $(FCLAW_CLAWPACK5_CPPFLAGS) \
              $(FCLAW_CLAWPATCH_CPPFLAGS)


## Add headers/libraries from this library to general header pool for other libraries/apps
LDADD        += $(FCLAW_CLAWPACK5_LDADD)
AM_CPPFLAGS  += $(FCLAW_CLAWPACK5_CPPFLAGS)

include_HEADERS += $(libclawpack5_installed_headers)



