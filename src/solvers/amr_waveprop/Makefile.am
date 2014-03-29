# -*- Makefile -*-

## This variable should list the .h and .H files from the source tree that
## should be installed in the "include" directory.

## The first part of the variable, "solver," must refer to a library or a
## program that is to be built.  That is, a target file that itself has
## been mentioned in bin_PROGRAMS for example.  Any variable used with +=
## must also first be initialized with = (usually to empty).

libwaveprop_generated_headers =
libwaveprop_installed_headers =

## TODO: there could be a name conflict with multiple claw.i files present
libwaveprop_internal_headers = src/solvers/amr_waveprop/fortran_source/claw.i

libwaveprop_compiled_sources = \
	src/solvers/amr_waveprop/amr_waveprop.cpp \
	src/solvers/amr_waveprop/fortran_source/limiter.f \
	src/solvers/amr_waveprop/fortran_source/src2.f \
	src/solvers/amr_waveprop/fortran_source/bc2.f \
	src/solvers/amr_waveprop/fortran_source/flux2.f \
	src/solvers/amr_waveprop/fortran_source/step2.f \
	src/solvers/amr_waveprop/fortran_source/ClawPatch2.f

lib_LTLIBRARIES += src/solvers/amr_waveprop/libwaveprop.la
src_solvers_amr_waveprop_libwaveprop_la_SOURCES = \
 $(libwaveprop_internal_headers) \
 $(libwaveprop_compiled_sources)
src_solvers_amr_waveprop_libwaveprop_la_CPPFLAGS = \
	$(AM_CPPFLAGS) -I@top_srcdir@/src/solvers/amr_waveprop

nodist_include_HEADERS += $(libwaveprop_generated_headers)
include_HEADERS += $(libwaveprop_installed_headers)

## We don't need to change global LDADD or AM_CPPFLAGS variables here.
## This is only done in src/Makefile.am, which is special in this regard.
