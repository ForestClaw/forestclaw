bin_PROGRAMS += applications/clawpack/advection/2d/sphere/sphere

TESTS += applications/clawpack/advection/2d/sphere/regressions.sh

applications_clawpack_advection_2d_sphere_sphere_SOURCES = \
	applications/clawpack/advection/2d/sphere/sphere.cpp \
	applications/clawpack/advection/2d/sphere/sphere_user.cpp \
	applications/clawpack/advection/2d/sphere/sphere_user.h \
	applications/clawpack/advection/2d/sphere/sphere_options.c \
	applications/clawpack/advection/2d/sphere/psi.f \
	applications/clawpack/advection/2d/sphere/setprob.f \
	applications/clawpack/advection/2d/sphere/fdisc.f \
	applications/clawpack/advection/2d/sphere/user_4.6/qinit.f \
	applications/clawpack/advection/2d/sphere/user_5.0/qinit.f90 \
	applications/clawpack/advection/2d/all/advection_patch_setup_manifold.c \
	applications/clawpack/advection/2d/all/fclaw2d_map_cubedsphere.c \
	applications/clawpack/advection/2d/all/fclaw2d_map_pillowsphere.c \
	applications/clawpack/advection/2d/all/user46_setaux_manifold.f \
	applications/clawpack/advection/2d/all/user5_setaux_manifold.f90 \
	applications/clawpack/advection/2d/all/clawpack46_rpn2adv_manifold.f \
	applications/clawpack/advection/2d/all/clawpack46_rpt2adv_manifold.f \
	applications/clawpack/advection/2d/all/clawpack5_rpn2adv_manifold.f90 \
	applications/clawpack/advection/2d/all/clawpack5_rpt2adv_manifold.f90

EXTRA_DIST += \
	applications/clawpack/advection/2d/sphere/fclaw_options.ini \
	applications/clawpack/advection/2d/sphere/regression.ini \
	applications/clawpack/advection/2d/sphere/afterframe.m \
	applications/clawpack/advection/2d/sphere/setplot2.m \
	applications/clawpack/advection/2d/sphere/mapc2m.m

applications_clawpack_advection_2d_sphere_sphere_CPPFLAGS = \
        $(AM_CPPFLAGS) \
        $(FCLAW_CLAWPACK46_CPPFLAGS) \
        $(FCLAW_CLAWPACK5_CPPFLAGS) \
        $(FCLAW_CLAWPATCH_CPPFLAGS)

applications_clawpack_advection_2d_sphere_sphere_LDADD = \
        $(LDADD) \
        $(FCLAW_CLAWPACK46_LDADD) \
        $(FCLAW_CLAWPACK5_LDADD) \
        $(FCLAW_CLAWPATCH_LDADD) \
        $(FCLAW_LDADD)
