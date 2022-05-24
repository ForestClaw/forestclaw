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
    src/solvers/fc2d_clawpack4.6/fortran_source/clawpack46_time_sync.f \
	src/solvers/fc2d_clawpack4.6/fortran_source/clawpack46_inlinelimiter.f \
	src/solvers/fc2d_clawpack4.6/fortran_source/clawpack46_set_capacity.f \
	src/solvers/fc2d_clawpack4.6/fortran_source/clawpack46_block.f \
	src/solvers/fc2d_clawpack4.6/fortran_source/clawpack46_bc2_default.f \
	src/solvers/fc2d_clawpack4.6/fortran_source/clawpack46_flux2.f \
	src/solvers/fc2d_clawpack4.6/fortran_source/clawpack46_flux2fw.f \
	src/solvers/fc2d_clawpack4.6/fortran_source/clawpack46_step2.f \
	src/solvers/fc2d_clawpack4.6/fortran_source/clawpack46_step2_wrap.f


## Name of library to build
lib_LTLIBRARIES += src/solvers/fc2d_clawpack4.6/libclawpack4.6.la

## Named variables that can be referenced from other libraries/apps
FCLAW_CLAWPACK46_CPPFLAGS = -I@top_srcdir@/src/solvers/fc2d_clawpack4.6 
FCLAW_CLAWPACK46_LDADD    = @top_builddir@/src/solvers/fc2d_clawpack4.6/libclawpack4.6.la

## Sources needed to build this library
src_solvers_fc2d_clawpack4_6_libclawpack4_6_la_SOURCES = \
          $(libclawpack4_6_internal_headers) \
          $(libclawpack4_6_compiled_sources)

## Headers needed to build this library.  Only AM_CPPFLAGS required; others are 
## included in AM_CPPFLAGS, but are included here for clarity (although makes 
## for lengthy compile lines)
src_solvers_fc2d_clawpack4.6_libclawpack4.6_la_CPPFLAGS = \
          $(AM_CPPFLAGS) \
          $(FCLAW_CPPFLAGS) \
          $(FCLAW_CLAWPACK46_CPPFLAGS) \
          $(FCLAW_CLAWPATCH_CPPFLAGS)

## Add headers from this library to general header pool, to be used by other apps
AM_CPPFLAGS  += $(FCLAW_CLAWPACK46_CPPFLAGS)

LDADD        += $(FCLAW_CLAWPACK46_LDADD)

include_HEADERS += $(libclawpack4_6_installed_headers)

## UNIT TESTS
check_PROGRAMS += src/solvers/fc2d_clawpack4.6/libclawpack4.6_TEST.MPI1
TESTS += src/solvers/fc2d_clawpack4.6/libclawpack4.6_TEST.MPI1

src_solvers_fc2d_clawpack4_6_libclawpack4_6_TEST_MPI1_SOURCES = \
    src/solvers/fc2d_clawpack4.6/fc2d_clawpack46.h.TEST.cpp

src_solvers_fc2d_clawpack4_6_libclawpack4_6_TEST_MPI1_CPPFLAGS = \
	$(test_libtestutils_la_CPPFLAGS) \
    $(FCLAW_CLAWPACK46_CPPFLAGS)

src_solvers_fc2d_clawpack4_6_libclawpack4_6_TEST_MPI1_LDADD = \
	test/libtestutils.la \
	$(test_libtestutils_la_LDADD) \
    $(LDADD) \
    $(FCLAW_CLAWPACK46_LDADD) \
    $(FCLAW_CLAWPATCH_LDADD) \
    $(FCLAW_LDADD)