bin_PROGRAMS += applications/clawpack/advection/3d/latlong/latlong

TESTS += applications/clawpack/advection/3d/latlong/regressions.sh

applications_clawpack_advection_3d_latlong_latlong_SOURCES = \
	applications/clawpack/advection/3d/latlong/latlong.cpp \
	applications/clawpack/advection/3d/latlong/latlong_options.c \
	applications/clawpack/advection/3d/latlong/latlong_user.h \
	applications/clawpack/advection/3d/latlong/latlong_user.cpp \
	applications/clawpack/advection/3d/latlong/latlong_map_extrude.c \
	applications/clawpack/advection/3d/latlong/fdisc.f \
	applications/clawpack/advection/3d/latlong/setprob.f90 \
	applications/clawpack/advection/3d/latlong/psi.f90 \
	applications/clawpack/advection/3d/latlong/user_4.6/qinit.f90 \
	applications/clawpack/advection/3d/all/advection_patch_setup_manifold.c \
	applications/clawpack/advection/2d/all/fclaw2d_map_latlong.c \
	applications/clawpack/advection/3d/all/claw3_setaux_manifold.f90 \
	applications/clawpack/advection/3d/all/rpn3_vc_advection.f90 \
	applications/clawpack/advection/3d/all/rpt3_vc_advection.f90 \
	applications/clawpack/advection/3d/all/rptt3_vc_advection.f90

EXTRA_DIST += \
	applications/clawpack/advection/3d/latlong/fclaw_options.ini \
	applications/clawpack/advection/3d/latlong/regression.ini \
	applications/clawpack/advection/3d/latlong/afterframe.m \
	applications/clawpack/advection/3d/latlong/setplot2.m \
	applications/clawpack/advection/3d/latlong/mapc2m.m

applications_clawpack_advection_3d_latlong_latlong_CPPFLAGS = \
	$(AM_CPPFLAGS) \
	$(FCLAW_CLAWPACK46_CPPFLAGS) \
	$(FCLAW_CLAWPATCH_CPPFLAGS)


applications_clawpack_advection_3d_latlong_latlong_LDADD = \
        $(LDADD) \
        $(FCLAW_CLAWPACK46_LDADD) \
        $(FCLAW_CLAWPATCH_LDADD)  \
        $(FCLAW_LDADD)
