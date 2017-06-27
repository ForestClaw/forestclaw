# -*- Makefile -*-
libclawpack5_3d_generated_headers =
libclawpack5_3d_installed_headers = \
	src/solvers/fc3d_clawpack5/fc3d_clawpack5.h \
	src/solvers/fc3d_clawpack5/clawpack5_user_fort.h \
	src/solvers/fc3d_clawpack5/fc3d_clawpack5_options.h

libclawpack5_3d_compiled_sources = \
    src/solvers/fc3d_clawpack5/fc3d_clawpack5.cpp \
    src/solvers/fc3d_clawpack5/fc3d_clawpack5_options.c \
	src/solvers/fc3d_clawpack5/fclaw2d_source/fc3d_clawpack5_average.f \
	src/solvers/fc3d_clawpack5/fclaw2d_source/fc3d_clawpack5_copy.f \
	src/solvers/fc3d_clawpack5/fclaw2d_source/fc3d_clawpack5_interpolate.f \
	src/solvers/fc3d_clawpack5/fclaw2d_source/fc3d_clawpack5_ghostpack_qarea.f \
	src/solvers/fc3d_clawpack5/fclaw2d_source/fc3d_clawpack5_diagnostics.f \
	src/solvers/fc3d_clawpack5/fclaw2d_source/fc3d_clawpack5_tag4refinement.f \
	src/solvers/fc3d_clawpack5/fclaw2d_source/fc3d_clawpack5_tag4coarsening.f \
	src/solvers/fc3d_clawpack5/fclaw2d_source/fc3d_clawpack5_output_ascii.f \
	src/solvers/fc3d_clawpack5/fclaw2d_source/fc3d_clawpack5_timeinterp.f

lib_LTLIBRARIES += src/solvers/fc3d_clawpack5/libclawpack5_3d.la

src_solvers_fc3d_clawpack5_libclawpack5_3d_la_SOURCES = \
  $(libclawpack5_3d_internal_headers) \
  $(libclawpack5_3d_compiled_sources)

src_solvers_fc3d_clawpack5_libclawpack5_3d_la_CPPFLAGS = \
  $(AM_CPPFLAGS) $(FCLAW_CPPFLAGS)

FCLAW_CLAWPACK5_3D_CPPFLAGS = -I@top_srcdir@/src/solvers/fc3d_clawpack5
FCLAW_CLAWPACK5_3D_LDADD = @top_builddir@/src/solvers/fc3d_clawpack5/libclawpack5_3d.la

include_HEADERS += $(libclawpack5_installed_headers)
