bin_PROGRAMS += applications/clawpack/advection/2d/latlong/latlong

TESTS += applications/clawpack/advection/2d/latlong/regressions.sh

applications_clawpack_advection_2d_latlong_latlong_SOURCES = \
	applications/clawpack/advection/2d/latlong/latlong.cpp \
	applications/clawpack/advection/2d/latlong/latlong_options.c \
	applications/clawpack/advection/2d/latlong/latlong_user.h \
	applications/clawpack/advection/2d/latlong/latlong_user.cpp \
	applications/clawpack/advection/2d/latlong/fdisc.f \
	applications/clawpack/advection/2d/latlong/setprob.f90 \
	applications/clawpack/advection/2d/latlong/psi.f \
	applications/clawpack/advection/2d/latlong/user_4.6/qinit.f \
	applications/clawpack/advection/2d/latlong/user_5.0/qinit.f90 \
	applications/clawpack/advection/2d/all/advection_patch_setup_manifold.c \
	applications/clawpack/advection/2d/all/fclaw2d_map_latlong.c \
	applications/clawpack/advection/2d/all/user46_setaux_manifold.f \
	applications/clawpack/advection/2d/all/user5_setaux_manifold.f90 \
	applications/clawpack/advection/2d/all/clawpack46_rpn2adv_manifold.f \
	applications/clawpack/advection/2d/all/clawpack46_rpt2adv_manifold.f \
	applications/clawpack/advection/2d/all/clawpack5_rpn2adv_manifold.f90 \
	applications/clawpack/advection/2d/all/clawpack5_rpt2adv_manifold.f90

EXTRA_DIST += \
	applications/clawpack/advection/2d/latlong/fclaw_options.ini \
	applications/clawpack/advection/2d/latlong/regression.ini \
	applications/clawpack/advection/2d/latlong/afterframe.m \
	applications/clawpack/advection/2d/latlong/setplot2.m \
	applications/clawpack/advection/2d/latlong/mapc2m.m

applications_clawpack_advection_2d_latlong_latlong_CPPFLAGS = \
	$(AM_CPPFLAGS) \
	$(FCLAW_CLAWPACK46_CPPFLAGS) \
	$(FCLAW_CLAWPACK5_CPPFLAGS) \
	$(FCLAW_CLAWPATCH_CPPFLAGS)


applications_clawpack_advection_2d_latlong_latlong_LDADD = \
        $(LDADD) \
        $(FCLAW_CLAWPACK46_LDADD) \
        $(FCLAW_CLAWPACK5_LDADD)  \
        $(FCLAW_CLAWPATCH_LDADD)  \
        $(FCLAW_LDADD)
