# -*- Makefile -*-
libclawpack5_generated_headers =
libclawpack5_installed_headers = \
	src/solvers/fc2d_clawpack5/fc2d_clawpack5.h \
	src/solvers/fc2d_clawpack5/fc2d_clawpack5_options.h \
	src/solvers/fc2d_clawpack5/clawpack5_user_fort.h \
	src/solvers/fc2d_clawpack5/fc2d_clawpack5_fort.h

libclawpack5_compiled_sources = \
	src/solvers/fc2d_clawpack5/fc2d_clawpack5.cpp \
	src/solvers/fc2d_clawpack5/fc2d_clawpack5_options.c \
	src/solvers/fc2d_clawpack5/fortran_source/amr_module.f90 \
	src/solvers/fc2d_clawpack5/fortran_source/set_amr_module.f90 \
	src/solvers/fc2d_clawpack5/fortran_source/clawpack5_inlinelimiter.f \
	src/solvers/fc2d_clawpack5/fortran_source/clawpack5_limiter.f \
	src/solvers/fc2d_clawpack5/fortran_source/clawpack5_set_capacity.f \
	src/solvers/fc2d_clawpack5/fortran_source/clawpack5_block.f \
	src/solvers/fc2d_clawpack5/fortran_source/clawpack5_bc2_default.f \
	src/solvers/fc2d_clawpack5/fortran_source/clawpack5_flux2.f \
	src/solvers/fc2d_clawpack5/fortran_source/clawpack5_step2.f90 \
	src/solvers/fc2d_clawpack5/fortran_source/clawpack5_step2_wrap.f \
	src/solvers/fc2d_clawpack5/fortran_source/philim.f \
	src/solvers/fc2d_clawpack5/fclaw2d_source/fc2d_clawpack5_average.f \
	src/solvers/fc2d_clawpack5/fclaw2d_source/fc2d_clawpack5_copy.f \
	src/solvers/fc2d_clawpack5/fclaw2d_source/fc2d_clawpack5_interpolate.f \
	src/solvers/fc2d_clawpack5/fclaw2d_source/fc2d_clawpack5_local_ghost_pack.f \
	src/solvers/fc2d_clawpack5/fclaw2d_source/fc2d_clawpack5_diagnostics.f \
	src/solvers/fc2d_clawpack5/fclaw2d_source/fc2d_clawpack5_tag4refinement.f \
	src/solvers/fc2d_clawpack5/fclaw2d_source/fc2d_clawpack5_tag4coarsening.f \
	src/solvers/fc2d_clawpack5/fclaw2d_source/fc2d_clawpack5_output_ascii.f \
	src/solvers/fc2d_clawpack5/fclaw2d_source/fc2d_clawpack5_timeinterp.f

lib_LTLIBRARIES += src/solvers/fc2d_clawpack5/libclawpack5.la

src_solvers_fc2d_clawpack5_libclawpack5_la_SOURCES = \
  $(libclawpack5_internal_headers) \
  $(libclawpack5_compiled_sources)

src_solvers_fc2d_clawpack5_libclawpack5_la_CPPFLAGS = \
  $(AM_CPPFLAGS) $(FCLAW_CPPFLAGS)

FCLAW_CLAWPACK5_CPPFLAGS = -I@top_srcdir@/src/solvers/fc2d_clawpack5
FCLAW_CLAWPACK5_LDADD = @top_builddir@/src/solvers/fc2d_clawpack5/libclawpack5.la

include_HEADERS += $(libclawpack5_installed_headers)
