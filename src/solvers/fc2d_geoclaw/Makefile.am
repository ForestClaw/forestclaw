# -*- Makefile -*-
libgeoclaw_generated_headers =
libgeoclaw_installed_headers = \
	src/solvers/fc2d_clawpack5/fc2d_clawpack5.h \
	src/solvers/fc2d_geoclaw/fc2d_geoclaw.h \
	src/solvers/fc2d_geoclaw/types.h \
	src/solvers/fc2d_geoclaw/fc2d_geoclaw_options.h \
	src/solvers/fc2d_geoclaw/fc2d_geoclaw_gauges_default.h

libgeoclaw_compiled_sources = \
	src/solvers/fc2d_geoclaw/fc2d_geoclaw.cpp \
	src/solvers/fc2d_geoclaw/fc2d_geoclaw_options.c \
	src/solvers/fc2d_geoclaw/fc2d_geoclaw_gauges_default.c \
	src/solvers/fc2d_geoclaw/fc2d_geoclaw_run.c \
	src/solvers/fc2d_geoclaw/fc2d_geoclaw_output_ascii.c \
	src/solvers/fc2d_geoclaw/amrlib_source/amr_module.f90 \
	src/solvers/fc2d_geoclaw/geolib_source/utility_module.f90 \
	src/solvers/fc2d_geoclaw/geolib_source/geoclaw_module.f90 \
	src/solvers/fc2d_geoclaw/amrlib_source/regions_module.f90 \
	src/solvers/fc2d_geoclaw/geolib_source/topo_module.f90 \
	src/solvers/fc2d_geoclaw/geolib_source/qinit_module.f90 \
	src/solvers/fc2d_geoclaw/geolib_source/fixedgrids_module.f90 \
	src/solvers/fc2d_geoclaw/geolib_source/refinement_module.f90 \
	src/solvers/fc2d_geoclaw/geolib_source/surge/model_storm_module.f90 \
	src/solvers/fc2d_geoclaw/geolib_source/surge/data_storm_module.f90 \
	src/solvers/fc2d_geoclaw/geolib_source/surge/storm_module.f90 \
	src/solvers/fc2d_geoclaw/geolib_source/friction_module.f90 \
	src/solvers/fc2d_geoclaw/riemann_source/rpn2_geoclaw.f \
	src/solvers/fc2d_geoclaw/riemann_source/rpt2_geoclaw.f \
	src/solvers/fc2d_geoclaw/riemann_source/geoclaw_riemann_utils.f \
	src/solvers/fc2d_geoclaw/amrlib_source/opendatafile.f \
	src/solvers/fc2d_geoclaw/amrlib_source/check4nans.f90 \
	src/solvers/fc2d_geoclaw/amrlib_source/quick_sort1.f \
	src/solvers/fc2d_geoclaw/geolib_source/topointegral.f \
	src/solvers/fc2d_geoclaw/geolib_source/topo_update.f90 \
	src/solvers/fc2d_geoclaw/geolib_source/cellgridintegrate2.f \
	src/solvers/fc2d_geoclaw/geolib_source/allowflag.f \
	src/solvers/fc2d_geoclaw/geolib_source/bilinearintegral.f \
	src/solvers/fc2d_geoclaw/fortran_source/philim.f \
	src/solvers/fc2d_geoclaw/fortran_source/geoclaw_set_modules_fort.f90 \
	src/solvers/fc2d_geoclaw/fortran_source/geoclaw_step2_wrap_fort.f \
	src/solvers/fc2d_geoclaw/fortran_source/geoclaw_flux2fw_fort.f \
	src/solvers/fc2d_geoclaw/fortran_source/geoclaw_limiter_fort.f \
	src/solvers/fc2d_geoclaw/fortran_source/geoclaw_block_fort.f \
	src/solvers/fc2d_geoclaw/fortran_source/geoclaw_gauges_fort.f90  \
	src/solvers/fc2d_geoclaw/fortran_source/geoclaw_topo_update_fort.f \
	src/solvers/fc2d_geoclaw/fortran_source/geoclaw_bc2_fort.f \
	src/solvers/fc2d_geoclaw/fortran_source/geoclaw_step2_fort.f90 \
	src/solvers/fc2d_geoclaw/fortran_source/geoclaw_setaux_fort.f90 \
	src/solvers/fc2d_geoclaw/fortran_source/geoclaw_src2_fort.f90 \
	src/solvers/fc2d_geoclaw/fortran_source/geoclaw_b4step2_fort.f90 \
	src/solvers/fc2d_geoclaw/fortran_source/geoclaw_qinit_fort.f90 \
	src/solvers/fc2d_geoclaw/fclaw2d_source/fc2d_geoclaw_copy_fort.f \
	src/solvers/fc2d_geoclaw/fclaw2d_source/fc2d_geoclaw_average_fort.f \
	src/solvers/fc2d_geoclaw/fclaw2d_source/fc2d_geoclaw_interpolate_fort.f \
	src/solvers/fc2d_geoclaw/fclaw2d_source/fc2d_geoclaw_output_ascii_fort.f \
	src/solvers/fc2d_geoclaw/fclaw2d_source/fc2d_geoclaw_tag4refinement_fort.f90 \
	src/solvers/fc2d_geoclaw/fclaw2d_source/fc2d_geoclaw_tag4coarsening_fort.f90 \
	src/solvers/fc2d_geoclaw/fclaw2d_source/fc2d_geoclaw_flag2refine_fort.f90 \
	src/solvers/fc2d_geoclaw/fclaw2d_source/fc2d_geoclaw_test_regions.f90 \
	src/solvers/fc2d_geoclaw/fclaw2d_source/fc2d_geoclaw_local_ghost_pack_fort.f \
	src/solvers/fc2d_geoclaw/fclaw2d_source/fc2d_geoclaw_local_ghost_pack_aux_fort.f \
	src/solvers/fc2d_geoclaw/fclaw2d_source/fc2d_geoclaw_diagnostics_fort.f \
	src/solvers/fc2d_geoclaw/fclaw2d_source/fc2d_geoclaw_timeinterp_fort.f

lib_LTLIBRARIES += src/solvers/fc2d_geoclaw/libgeoclaw.la

## Named variables that can be reference from other libraries/apps
FCLAW_GEOCLAW_CPPFLAGS = -I@top_srcdir@/src/solvers/fc2d_geoclaw
FCLAW_GEOCLAW_LDADD    = @top_builddir@/src/solvers/fc2d_geoclaw/libgeoclaw.la

## Sources needed to build this library
src_solvers_fc2d_geoclaw_libgeoclaw_la_SOURCES = \
  $(libgeoclaw_internal_headers) \
  $(libgeoclaw_compiled_sources)

## Headers needed to build this library.  
src_solvers_fc2d_geoclaw_libgeoclaw_la_CPPFLAGS = \
  	$(AM_CPPFLAGS) \
  	$(FCLAW_CPPFLAGS) \
	$(FCLAW_GEOCLAW_CPPFLAGS) \
    $(FCLAW_CLAWPATCH_CPPFLAGS)


AM_CPPFLAGS += $(FCLAW_GEOCLAW_CPPFLAGS)
LDADD += $(FCLAW_GEOCLAW_LDADD)

include_HEADERS += $(libgeoclaw_installed_headers)
