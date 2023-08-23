## Makefile.am in applications/clawpack/acoustics_2d_radial

bin_PROGRAMS += applications/clawpack/acoustics/2d/radial/radial

TESTS += applications/clawpack/acoustics/2d/radial/regressions.txt

applications_clawpack_acoustics_2d_radial_radial_SOURCES = \
	applications/clawpack/acoustics/2d/radial/radial.cpp \
	applications/clawpack/acoustics/2d/radial/radial_user.cpp \
	applications/clawpack/acoustics/2d/radial/radial_options.c \
	applications/clawpack/acoustics/2d/radial/radial_user.h \
	applications/clawpack/acoustics/2d/radial/fclaw2d_map_pillowdisk5.c \
	applications/clawpack/acoustics/2d/radial/fclaw2d_map_pillowdisk.c \
	applications/clawpack/acoustics/2d/radial/setprob.f \
	applications/clawpack/acoustics/2d/radial/user_4.6/qinit.f \
	applications/clawpack/acoustics/2d/radial/user_4.6/setaux.f90 \
	applications/clawpack/acoustics/2d/radial/user_5.0/qinit.f90 \
	applications/clawpack/acoustics/2d/radial/user_5.0/setaux.f90 \
	applications/clawpack/acoustics/2d/rp/clawpack46_rpn2_acoustics.f \
	applications/clawpack/acoustics/2d/rp/clawpack46_rpt2_acoustics.f \
	applications/clawpack/acoustics/2d/rp/clawpack46_rpn2_acoustics_manifold.f90 \
	applications/clawpack/acoustics/2d/rp/clawpack46_rpt2_acoustics_manifold.f90 \
	applications/clawpack/acoustics/2d/rp/clawpack5_rpn2_acoustics.f90 \
	applications/clawpack/acoustics/2d/rp/clawpack5_rpt2_acoustics.f90 \
	applications/clawpack/acoustics/2d/rp/clawpack5_rpn2_acoustics_manifold.f90 \
	applications/clawpack/acoustics/2d/rp/clawpack5_rpt2_acoustics_manifold.f90

EXTRA_DIST += \
	applications/clawpack/acoustics/2d/radial/fclaw_options.ini \
	applications/clawpack/acoustics/2d/radial/regression.ini \
	applications/clawpack/acoustics/2d/radial/afterframe.m \
	applications/clawpack/acoustics/2d/radial/setplot2.m \
	applications/clawpack/acoustics/2d/radial/mapc2m.m \
	applications/clawpack/acoustics/2d/radial/map1d.m \
	applications/clawpack/acoustics/2d/radial/make_data.py \
	applications/clawpack/acoustics/2d/radial/make_plots.py \
	applications/clawpack/acoustics/2d/radial/1drad/setplot.py \
	applications/clawpack/acoustics/2d/radial/1drad/setrun.py

radial_headers = \
	applications/clawpack/acoustics/2d/radial/radial_user.h

radial_inifiles = \
	applications/clawpack/acoustics/2d/radial/fclaw_options.ini \
	applications/clawpack/acoustics/2d/radial/regression.ini

applications_clawpack_acoustics_2d_radial_radial_SOURCES += \
	$(radial_headers) \
	$(radial_inifiles)

applications_clawpack_acoustics_2d_radial_radial_CPPFLAGS = \
        $(AM_CPPFLAGS) \
        $(FCLAW_CLAWPACK46_CPPFLAGS) \
        $(FCLAW_CLAWPACK5_CPPFLAGS) \
        $(FCLAW_CLAWPATCH_CPPFLAGS)

applications_clawpack_acoustics_2d_radial_radial_LDADD = \
        $(LDADD) \
        $(FCLAW_CLAWPACK46_LDADD) \
        $(FCLAW_CLAWPACK5_LDADD) \
        $(FCLAW_CLAWPATCH_LDADD) \
        $(FCLAW_LDADD)
