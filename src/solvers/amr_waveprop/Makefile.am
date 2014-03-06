# -*- Makefile -*-

## This variable should list the .h and .H files from the source tree that
## should be installed in the "include" directory.

## The first part of the variable, "solver," must refer to a library or a
## program that is to be built.  That is, a target file that itself has
## been mentioned in bin_PROGRAMS for example.  Any variable used with +=
## must also first be initialized with = (usually to empty).

waveprop_compiled_sources = \
	src/solvers/amr_waveprop/amr_waveprop.cpp \
	src/solvers/amr_waveprop/fortran_source/limiter.f \
	src/solvers/amr_waveprop/fortran_source/src2.f \
	src/solvers/amr_waveprop/fortran_source/bc2.f \
	src/solvers/amr_waveprop/fortran_source/flux2.f \
	src/solvers/amr_waveprop/fortran_source/step2.f \
	src/solvers/amr_waveprop/fortran_source/ClawPatch2.f

include_HEADERS += src/solvers/amr_waveprop/amr_waveprop.H src/solvers/amr_waveprop/fortran_source/claw.i

lib_LTLIBRARIES += src/solvers/amr_waveprop/libwaveprop.la
src_solvers_amr_waveprop_libwaveprop_la_SOURCES = \
 $(libwaveprop_internal_headers) \
 $(waveprop_compiled_sources)

# src_libforestclaw_la_CPPFLAGS = $(AM_CPPFLAGS) $(FCLAW_CPPFLAGS)
