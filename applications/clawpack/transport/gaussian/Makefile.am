## Makefile.am in applications/clawpack/transport/gaussian

bin_PROGRAMS += applications/clawpack/transport/gaussian/gaussian

applications_clawpack_transport_gaussian_gaussian_SOURCES = \
	applications/clawpack/transport/gaussian/gaussian.cpp \
	applications/clawpack/transport/gaussian/gaussian_user.h \
	applications/clawpack/transport/gaussian/gradient_exceeds_th.f90 \
	applications/clawpack/transport/gaussian/user_4.6/qinit.f \
	applications/clawpack/transport/gaussian/user_5.0/qinit.f90 \
	applications/clawpack/transport/gaussian/user_4.6/tag4refinement.f \
	applications/clawpack/transport/gaussian/user_5.0/tag4refinement.f \
	applications/clawpack/transport/gaussian/user_4.6/tag4coarsening.f \
	applications/clawpack/transport/gaussian/user_5.0/tag4coarsening.f \
	applications/clawpack/transport/all/transport_setprob.f \
	applications/clawpack/transport/all/velocity_wind.f \
	applications/clawpack/transport/all/gaussian.f \
	applications/clawpack/transport/all/transport_user.h \
	applications/clawpack/transport/all/transport_options.h \
	applications/clawpack/transport/all/transport_user.cpp \
	applications/clawpack/transport/all/transport_domain.c \
	applications/clawpack/transport/all/transport_options.c \
	applications/clawpack/transport/all/user46_b4step2_manifold.f \
	applications/clawpack/transport/all/user5_b4step2_manifold.f90 \
	applications/clawpack/transport/all/fclaw2d_map_cubedsphere.c \
	applications/clawpack/transport/all/fclaw2d_map_pillowsphere.c \
	applications/clawpack/advection/2d/all/advection_user_fort.h \
	applications/clawpack/advection/2d/all/user46_setaux_manifold.f \
	applications/clawpack/advection/2d/all/user5_setaux_manifold.f90 \
	applications/clawpack/advection/2d/all/clawpack46_rpn2adv_manifold.f \
	applications/clawpack/advection/2d/all/clawpack46_rpt2adv_manifold.f \
	applications/clawpack/advection/2d/all/clawpack5_rpn2adv_manifold.f90 \
	applications/clawpack/advection/2d/all/clawpack5_rpt2adv_manifold.f90

EXTRA_DIST += \
	applications/clawpack/transport/gaussian/fclaw_options.ini \
	applications/clawpack/transport/gaussian/regression.ini \
	applications/clawpack/transport/gaussian/afterframe.m \
	applications/clawpack/transport/gaussian/setplot2.m \
	applications/clawpack/transport/gaussian/mapc2m.m

applications_clawpack_transport_gaussian_gaussian_CPPFLAGS = \
        $(AM_CPPFLAGS) \
        $(FCLAW_CLAWPACK46_CPPFLAGS) \
        $(FCLAW_CLAWPACK5_CPPFLAGS) \
        $(FCLAW_CLAWPATCH_CPPFLAGS)

applications_clawpack_transport_gaussian_gaussian_LDADD = \
	    $(LDADD) \
	    $(FCLAW_CLAWPACK46_LDADD) \
        $(FCLAW_CLAWPACK5_LDADD) \
        $(FCLAW_CLAWPATCH_LDADD) \
	    $(FCLAW_LDADD)