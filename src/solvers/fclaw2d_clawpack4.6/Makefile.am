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

## TODO: Preserve some subdirectory structure for the make install
libclawpack4_6_generated_headers =
libclawpack4_6_installed_headers = \
	src/solvers/fclaw2d_clawpack4.6/fclaw2d_clawpack.H

## TODO: there could be a name conflict with multiple claw.i files present
##       but that's ok for internal headers (and ONLY for these)
libclawpack4_6_internal_headers = src/solvers/fclaw2d_clawpack4.6/fortran_source/claw.i

## TODO: remove dependencies of this library on any sources under applications
libclawpack4_6_compiled_sources = \
	src/solvers/fclaw2d_clawpack4.6/fclaw2d_clawpack.cpp \
	src/solvers/fclaw2d_clawpack4.6/fortran_source/limiter.f \
	src/solvers/fclaw2d_clawpack4.6/fortran_source/qinit.f \
	src/solvers/fclaw2d_clawpack4.6/fortran_source/b4step2.f \
	src/solvers/fclaw2d_clawpack4.6/fortran_source/setaux.f \
	src/solvers/fclaw2d_clawpack4.6/fortran_source/src2.f \
	src/solvers/fclaw2d_clawpack4.6/fortran_source/corner_info.f \
	src/solvers/fclaw2d_clawpack4.6/fortran_source/bc2.f \
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
