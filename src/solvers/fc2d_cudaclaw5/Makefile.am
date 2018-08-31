# -*- Makefile -*-
libcudaclaw5_generated_headers =
libcudaclaw5_installed_headers = \
	src/solvers/fc2d_cudaclaw5/fc2d_cudaclaw5.h \
	src/solvers/fc2d_cudaclaw5/fc2d_cudaclaw5_options.h \
	src/solvers/fc2d_cudaclaw5/cudaclaw5_user_fort.h \
	src/solvers/fc2d_cudaclaw5/fc2d_cudaclaw5_fort.h

libcudaclaw5_compiled_sources = \
	src/solvers/fc2d_cudaclaw5/fc2d_cudaclaw5.cpp \
	src/solvers/fc2d_cudaclaw5/fc2d_cudaclaw5_options.c \
	src/solvers/fc2d_cudaclaw5/fortran_source/amr_module.f90 \
	src/solvers/fc2d_cudaclaw5/fortran_source/set_amr_module.f90 \
	src/solvers/fc2d_cudaclaw5/fortran_source/cudaclaw5_inlinelimiter.f \
	src/solvers/fc2d_cudaclaw5/fortran_source/cudaclaw5_limiter.f \
	src/solvers/fc2d_cudaclaw5/fortran_source/cudaclaw5_set_capacity.f \
	src/solvers/fc2d_cudaclaw5/fortran_source/cudaclaw5_block.f \
	src/solvers/fc2d_cudaclaw5/fortran_source/cudaclaw5_bc2_default.f \
	src/solvers/fc2d_cudaclaw5/fortran_source/cudaclaw5_flux2.f \
	src/solvers/fc2d_cudaclaw5/fortran_source/cudaclaw5_step2.f90 \
	src/solvers/fc2d_cudaclaw5/fortran_source/cudaclaw5_step2_wrap.f \
	src/solvers/fc2d_cudaclaw5/fortran_source/philim.f

## Name of library to build
lib_LTLIBRARIES += src/solvers/fc2d_cudaclaw5/libcudaclaw5.la

## Named variables that can be reference from other libraries/apps
FCLAW_CUDACLAW5_CPPFLAGS = -I@top_srcdir@/src/solvers/fc2d_cudaclaw5
FCLAW_CUDACLAW5_LDADD    = @top_builddir@/src/solvers/fc2d_cudaclaw5/libcudaclaw5.la

## Sources needed to build this library
src_solvers_fc2d_cudaclaw5_libcudaclaw5_la_SOURCES = \
  $(libcudaclaw5_internal_headers) \
  $(libcudaclaw5_compiled_sources)

## Headers needed to build this library.  Only AM_CPPFLAGS required; others are 
## included in AM_CPPFLAGS, but are included here for clarity (although makes 
## for lengthy compile lines)
src_solvers_fc2d_cudaclaw5_libcudaclaw5_la_CPPFLAGS = \
              $(AM_CPPFLAGS) \
              $(FCLAW_CPPFLAGS) \
              $(FCLAW_CUDACLAW5_CPPFLAGS) \
              $(FCLAW_CLAWPATCH_CPPFLAGS)


## Add headers/libraries from this library to general header pool for other libraries/apps
LDADD        += $(FCLAW_CUDACLAW5_LDADD)
AM_CPPFLAGS  += $(FCLAW_CUDACLAW5_CPPFLAGS)

include_HEADERS += $(libcudaclaw5_installed_headers)



