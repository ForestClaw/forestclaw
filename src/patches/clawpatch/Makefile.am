# -*- Makefile -*-
libclawpatch_generated_headers =

libclawpatch_installed_headers = \
	src/patches/clawpatch/fclaw2d_clawpatch.h \
	src/patches/clawpatch/fclaw2d_clawpatch.hpp \
	src/patches/clawpatch/fclaw2d_clawpatch_options.h \
	src/patches/clawpatch/fclaw2d_clawpatch_transform.h \
	src/patches/clawpatch/fclaw2d_clawpatch_diagnostics.h \
	src/patches/clawpatch/fclaw2d_clawpatch_pillow.h \
	src/patches/clawpatch/fclaw2d_clawpatch_fort.h \
	src/patches/clawpatch/fclaw2d_clawpatch46_fort.h \
	src/patches/clawpatch/fclaw2d_clawpatch5_fort.h \
	src/patches/clawpatch/fclaw2d_clawpatch_output_ascii.h \
	src/patches/clawpatch/fclaw2d_clawpatch_output_vtk.h \
	src/patches/metric/fclaw2d_metric.h \
	src/patches/metric/fclaw2d_metric.hpp \
	src/patches/metric/fclaw2d_metric_default.h \
	src/patches/metric/fclaw2d_metric_default_fort.h \
	src/patches/metric/fortran_source2d/metric_terms.i 

libclawpatch_compiled_sources = \
	src/patches/clawpatch/fclaw2d_clawpatch.cpp \
	src/patches/clawpatch/fclaw2d_clawpatch_options.c \
	src/patches/clawpatch/fclaw2d_clawpatch_diagnostics.c \
	src/patches/clawpatch/fclaw2d_clawpatch_diagnostics_default.c \
	src/patches/clawpatch/fclaw2d_clawpatch_pillow.c \
	src/patches/clawpatch/fclaw2d_clawpatch_transform.c \
	src/patches/clawpatch/fclaw2d_clawpatch_output_ascii.c \
	src/patches/clawpatch/fclaw2d_clawpatch_output_vtk.c \
	src/patches/clawpatch/fclaw2d_clawpatch_utils.f \
	src/patches/clawpatch/fort_4.6/fclaw2d_clawpatch46_average.f \
	src/patches/clawpatch/fort_4.6/fclaw2d_clawpatch46_copy.f \
	src/patches/clawpatch/fort_4.6/fclaw2d_clawpatch46_interpolate.f \
	src/patches/clawpatch/fort_4.6/fclaw2d_clawpatch46_local_ghost_pack.f \
	src/patches/clawpatch/fort_4.6/fclaw2d_clawpatch46_diagnostics.f \
	src/patches/clawpatch/fort_4.6/fclaw2d_clawpatch46_pillowsphere.f \
	src/patches/clawpatch/fort_4.6/fclaw2d_clawpatch46_output_ascii.f \
	src/patches/clawpatch/fort_4.6/fclaw2d_clawpatch46_tag4refinement.f \
	src/patches/clawpatch/fort_4.6/fclaw2d_clawpatch46_tag4coarsening.f \
	src/patches/clawpatch/fort_4.6/fclaw2d_clawpatch46_timeinterp.f  \
	src/patches/clawpatch/fort_5.0/fclaw2d_clawpatch5_average.f \
	src/patches/clawpatch/fort_5.0/fclaw2d_clawpatch5_copy.f \
	src/patches/clawpatch/fort_5.0/fclaw2d_clawpatch5_interpolate.f \
	src/patches/clawpatch/fort_5.0/fclaw2d_clawpatch5_local_ghost_pack.f \
	src/patches/clawpatch/fort_5.0/fclaw2d_clawpatch5_diagnostics.f \
	src/patches/clawpatch/fort_5.0/fclaw2d_clawpatch5_pillowsphere.f \
	src/patches/clawpatch/fort_5.0/fclaw2d_clawpatch5_output_ascii.f \
	src/patches/clawpatch/fort_5.0/fclaw2d_clawpatch5_tag4refinement.f \
	src/patches/clawpatch/fort_5.0/fclaw2d_clawpatch5_tag4coarsening.f \
	src/patches/clawpatch/fort_5.0/fclaw2d_clawpatch5_timeinterp.f  \
	src/patches/clawpatch/fclaw2d_clawpatch_conservation.c \
	src/patches/metric/fclaw2d_metric.cpp \
	src/patches/metric/fclaw2d_metric_default.c \
	src/patches/metric/fortran_source2d/fclaw2d_metric_default_fort.f


## Name of library to build
lib_LTLIBRARIES += src/patches/clawpatch/libclawpatch.la

## Named variables that can be referenced from other libraries/apps
FCLAW_CLAWPATCH_CPPFLAGS = -I@top_srcdir@/src/patches/clawpatch \
                           -I@top_srcdir@/src/patches/metric

## Needed for metric_terms.i
FCLAW_CLAWPATCH_FFLAGS = -I@top_srcdir@/src/patches/metric/fortran_source2d

FCLAW_CLAWPATCH_LDADD    = @top_builddir@/src/patches/clawpatch/libclawpatch.la

## Sources needed to build this library
src_patches_clawpatch_libclawpatch_la_SOURCES = \
        $(libclawpatch_internal_headers) \
        $(libclawpatch_compiled_sources)

## Headers needed to build this library.  
src_patches_clawpatch_libclawpatch_la_CPPFLAGS = \
                $(AM_CPPFLAGS)

## Add headers from this library to general header "pool", to be used by other apps
AM_CPPFLAGS += $(FCLAW_CLAWPATCH_CPPFLAGS)
LDADD       += $(FCLAW_CLAWPATCH_LDADD)
AM_FFLAGS   += $(FCLAW_CLAWPATCH_FFLAGS)


include_HEADERS += $(libclawpatch_installed_headers)
