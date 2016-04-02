# -*- Makefile -*-
libclawpack5_0_generated_headers =
libclawpack5_0_installed_headers = \
	src/solvers/fc2d_clawpack5.0/fc2d_clawpack50.h \
	src/solvers/fc2d_clawpack5.0/fc2d_clawpack50_options.h

libclawpack5_0_compiled_sources = \
	src/solvers/fc2d_clawpack5.0/fc2d_clawpack50.cpp \
	src/solvers/fc2d_clawpack5.0/fc2d_clawpack50_options.c \
	src/solvers/fc2d_clawpack5.0/fortran_source/clawpack50_inlinelimiter.f \
	src/solvers/fc2d_clawpack5.0/fortran_source/clawpack50_limiter.f \
	src/solvers/fc2d_clawpack5.0/fortran_source/clawpack50_set_capacity.f \
	src/solvers/fc2d_clawpack5.0/fortran_source/clawpack50_block.f \
	src/solvers/fc2d_clawpack5.0/fortran_source/clawpack50_bc2.f \
	src/solvers/fc2d_clawpack5.0/fortran_source/clawpack50_flux2.f \
	src/solvers/fc2d_clawpack5.0/fortran_source/clawpack50_flux2fw.f \
	src/solvers/fc2d_clawpack5.0/fortran_source/clawpack50_step2.f \
	src/solvers/fc2d_clawpack5.0/fortran_source/clawpack50_step2_wrap.f

lib_LTLIBRARIES += src/solvers/fc2d_clawpack5.0/libclawpack5.0.la

src_solvers_fc2d_clawpack5_0_libclawpack5_0_la_SOURCES = \
  $(libclawpack5_0_internal_headers) \
  $(libclawpack5_0_compiled_sources)

src_solvers_fc2d_clawpack5.0_libclawpack5.0_la_CPPFLAGS = \
  $(AM_CPPFLAGS) $(FCLAW_CPPFLAGS)

FCLAW_CLAWPACK50_CPPFLAGS = -I@top_srcdir@/src/solvers/fc2d_clawpack5.0
FCLAW_CLAWPACK50_LDADD = @top_builddir@/src/solvers/fc2d_clawpack5.0/libclawpack5.0.la

include_HEADERS += $(libclawpack5_0_installed_headers)
