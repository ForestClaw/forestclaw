# -*- Makefile -*-
libclawpack4_6_generated_headers =
libclawpack4_6_installed_headers = \
	src/solvers/fclaw2d_clawpack4.6/fclaw2d_clawpack.H

libclawpack4_6_internal_headers = src/solvers/fclaw2d_clawpack4.6/fortran_source/claw.i

libclawpack4_6_compiled_sources = \
	src/solvers/fclaw2d_clawpack4.6/fc2d_clawpack46.cpp \
	src/solvers/fclaw2d_clawpack4.6/fc2d_clawpack46_options.c \
	src/solvers/fclaw2d_clawpack4.6/fortran_source/limiter.f \
	src/solvers/fclaw2d_clawpack4.6/fortran_source/block.f \
	src/solvers/fclaw2d_clawpack4.6/fortran_source/corner_info.f \
	src/solvers/fclaw2d_clawpack4.6/fortran_source/clawpack46_bc2.f \
	src/solvers/fclaw2d_clawpack4.6/fortran_source/flux2.f \
	src/solvers/fclaw2d_clawpack4.6/fortran_source/step2.f \
	src/solvers/fclaw2d_clawpack4.6/fortran_source/ClawPatch2.f

lib_LTLIBRARIES += src/solvers/fclaw2d_clawpack4.6/libclawpack4.6.la

src_solvers_fclaw2d_clawpack4_6_libclawpack4_6_la_SOURCES = \
  $(libclawpack4_6_internal_headers) \
  $(libclawpack4_6_compiled_sources)

src_solvers_fclaw2d_clawpack4.6_libclawpack4.6_la_CPPFLAGS = \
  $(AM_CPPFLAGS) $(FCLAW_CPPFLAGS)

FCLAW_CLAWPACK46_CPPFLAGS = -I@top_srcdir@/src/solvers/fclaw2d_clawpack4.6
FCLAW_CLAWPACK46_LDADD = @top_builddir@/src/solvers/fclaw2d_clawpack4.6/libclawpack4.6.la

include_HEADERS += $(libclawpack4_6_installed_headers)
