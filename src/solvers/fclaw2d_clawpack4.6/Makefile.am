# -*- Makefile -*-

## TODO: are these comments still current?
##
## This variable should list the .h and .H files from the source tree that
## should be installed in the "include" directory.
##
## The first part of the variable, "solver," must refer to a library or a
## program that is to be built.  That is, a target file that itself has
## been mentioned in bin_PROGRAMS for example.  Any variable used with +=
## must also first be initialized with = (usually to empty).

libclawpack46_generated_headers =
libclawpack46_installed_headers = fclaw2d_clawpack.H

## TODO: there could be a name conflict with multiple claw.i files present
libclawpack46_internal_headers = src/solvers/fclaw2d_clawpack4.6/fortran_source/claw.i

libclawpack46_compiled_sources = \
	src/solvers/fclaw2d_clawpack4.6/fclaw2d_clawpack.cpp \
	src/solvers/fclaw2d_clawpack4.6/fortran_source/limiter.f \
	src/solvers/fclaw2d_clawpack4.6/fortran_source/qinit.f \
	src/solvers/fclaw2d_clawpack4.6/fortran_source/b4step2.f \
	src/solvers/fclaw2d_clawpack4.6/fortran_source/claw2d_utils.f \
	src/solvers/fclaw2d_clawpack4.6/fortran_source/setaux.f \
	src/solvers/fclaw2d_clawpack4.6/fortran_source/src2.f \
	src/solvers/fclaw2d_clawpack4.6/fortran_source/bc2.f \
	src/solvers/fclaw2d_clawpack4.6/fortran_source/flux2.f \
	src/solvers/fclaw2d_clawpack4.6/fortran_source/step2.f \
	src/solvers/fclaw2d_clawpack4.6/fortran_source/ClawPatch2.f

lib_LTLIBRARIES += @top_builddir@/src/solvers/fclaw2d_clawpack4.6/libclawpack46.la

src_solvers_fclaw2d_clawpack4.6_libclawpack46_la_SOURCES = \
  $(libclawpack46_internal_headers) \
  $(libclawpack46_compiled_sources)

src_solvers_fclaw2d_clawpack4.6_libclawpack46_la_CPPFLAGS = $(AM_CPPFLAGS) $(FCLAW_CPPFLAGS)

LDADD += @top_builddir@/src/solvers/fclaw2d_clawpack4.6/libclawpack46.la

include_HEADERS += $(libclawpack4.6_installed_headers)

AM_CPPFLAGS += -I@top_srcdir@/src/solvers/fclaw2d_clawpack4.6
