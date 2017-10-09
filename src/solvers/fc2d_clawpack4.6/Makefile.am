# -*- Makefile -*-
libclawpack4_6_generated_headers =
libclawpack4_6_installed_headers = \
	src/solvers/fc2d_clawpack4.6/fc2d_clawpack46.h \
	src/solvers/fc2d_clawpack4.6/fc2d_clawpack46_options.h \
	src/solvers/fc2d_clawpack4.6/clawpack46_user_fort.h \
	src/solvers/fc2d_clawpack4.6/fc2d_clawpack46_fort.h 

libclawpack4_6_compiled_sources = \
	src/solvers/fc2d_clawpack4.6/fc2d_clawpack46.cpp \
	src/solvers/fc2d_clawpack4.6/fc2d_clawpack46_options.c \
	src/solvers/fc2d_clawpack4.6/fortran_source/clawpack46_inlinelimiter.f \
	src/solvers/fc2d_clawpack4.6/fortran_source/clawpack46_set_capacity.f \
	src/solvers/fc2d_clawpack4.6/fortran_source/clawpack46_block.f \
	src/solvers/fc2d_clawpack4.6/fortran_source/clawpack46_bc2_default.f \
	src/solvers/fc2d_clawpack4.6/fortran_source/clawpack46_flux2.f \
	src/solvers/fc2d_clawpack4.6/fortran_source/clawpack46_flux2fw.f \
	src/solvers/fc2d_clawpack4.6/fortran_source/clawpack46_step2.f \
	src/solvers/fc2d_clawpack4.6/fortran_source/clawpack46_step2_wrap.f


lib_LTLIBRARIES += src/solvers/fc2d_clawpack4.6/libclawpack4.6.la

src_solvers_fc2d_clawpack4_6_libclawpack4_6_la_SOURCES = \
  $(libclawpack4_6_internal_headers) \
  $(libclawpack4_6_compiled_sources)

src_solvers_fc2d_clawpack4.6_libclawpack4.6_la_CPPFLAGS = \
  $(AM_CPPFLAGS) $(FCLAW_CPPFLAGS) $(FCLAW_CLAWPATCH_CPPFLAGS)

FCLAW_CLAWPACK46_CPPFLAGS = -I@top_srcdir@/src/solvers/fc2d_clawpack4.6 

LDADD += @top_builddir@/src/solvers/fc2d_clawpack4.6/libclawpack4.6.la

AM_CPPFLAGS += -I@top_srcdir@/src/solvers/fcd_clawpack4.6 $(FCLAW_CLAWPATCH_CPPFLAGS)

include_HEADERS += $(libclawpack4_6_installed_headers)
